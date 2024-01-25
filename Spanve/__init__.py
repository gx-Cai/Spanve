import time
import warnings
from collections import defaultdict
from itertools import combinations, combinations_with_replacement

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from joblib import Parallel, delayed
from scipy.sparse import issparse
try:
    from scipy.sparse import csr_array
except ImportError:
    from scipy.sparse import csr_matrix as csr_array
    warnings.warn("scipy.sparse.csr_array is not available. Recommend to install scipy >= 1.8")
import pickle
import os

from scipy.special import iv
from scipy.stats import chi2, entropy
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.metrics import calinski_harabasz_score as chs
from sklearn.metrics import davies_bouldin_score as dbs
from sklearn.metrics import silhouette_score as sil
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm, trange

__all__ = ['Spanve','adata_preprocess','adata_preprocess_int','Spanve_gpu','AutoCluster']

## ----- utils -----
class Distribution(defaultdict):
    def __init__(self, value=None, prob=None) -> None:
        super().__init__(float)
        if (value is not None) and (prob is not None):
            if not isinstance(value, (list, np.ndarray)) or \
                not isinstance(prob, (list, np.ndarray)
            ):
                value = list(value)
                prob = list(prob)

            for i, j in zip(value, prob):
                self[i] = j

    @property
    def counts(self):
        return np.array(list(super().keys()))

    @property
    def prob(self):
        return np.array(list(super().values()))

    def prob_check(self):
        assert self.prob.sum() == 1.0, "Probability does not sum to 1.0"

    def from_obs(self, obs):
        ind, count = np.unique(obs, return_counts=True)
        self.__init__(value=ind, prob=count / count.sum())
        return self

    def dist_abs_subtract(self):
        new_dist = Distribution()
        counts = self.counts
        prob = self.prob
        for i, j in combinations_with_replacement(range(len(self)), 2):
            new_dist[np.abs(counts[i] - counts[j])] += (
                prob[i] * prob[j] * (2 - (i == j))
            )
        return new_dist

    def __getitem__(self, __k):
        # if __k is a iterable
        if isinstance(__k, (list, tuple, np.ndarray)):
            return np.array([self[i] for i in __k])
        # if __k is a scalar
        else:
            return super().__getitem__(__k)

    def isf(self, n, alpha=0.05):
        try:
            return np.where(np.cumsum(self[np.arange(n)]) > (1 - alpha))[0][0] + 1
        except:
            return n+1
def poisson_dispersion_test(x):
    n = x.shape[0]
    D = np.std(x, axis=0) ** 2 / np.mean(x, axis=0) * n
    p = 2 * np.minimum(chi2.sf(D, n - 1), chi2.cdf(D, n - 1))
    return p

def ASP_pdf(x, lam):
    # Abosulte Substracted Possibility Distribution
    return (2 - (x == 0)) * np.exp(-2 * lam) * iv(x, 2 * lam)

def outer_subtract(x):
    return np.abs(x[:, np.newaxis] - x[np.newaxis, :])

def adata_preprocess(anndata,copy=True):
    if copy:
        anndata = anndata.copy()

    sc.pp.normalize_per_cell(anndata, counts_per_cell_after=10000)
    sc.pp.log1p(anndata)
    sc.pp.scale(anndata)
    # anndata.X = (anndata.X - anndata.X.mean(0)) / anndata.X.std(0)
    return anndata

def adata_preprocess_int(anndata,eps = 1e-7,exclude_highly_expressed=True,copy=True):
    if copy:
        adata = anndata.copy()

    sc.pp.normalize_total(
        adata,
        exclude_highly_expressed=exclude_highly_expressed,
        )
    sc.pp.log1p(adata)
    expr_median = eps+np.median(adata.X,axis=0)
    adata.X = ((adata.X / expr_median) * expr_median).astype(int)

    return adata

def elbow(X: np.ndarray) -> int:
    max_idx = np.argmax(X)

    X = X[max_idx:]  # truncate data from max (not origin) to endpoint.
    b = np.array([len(X), X[-1] - X[0]])  # Vector from origin to end.
    norm_vec = [0]  # Initial point ignored.

    for i in range(1, len(X)):
        p = np.array([i, X[i] - X[0]])  # Vector from origin to current point on curve.
        d = np.linalg.norm(p - (np.dot(p, b) / np.dot(b, b)) * b)  # Distance from point to b.

        norm_vec.append(d)

    # Pick the longest connecting line - note max_idx added to slice back into original data.
    return max_idx + np.argmax(norm_vec) 

def persent_select(arr,prop=0.8):
    arr_sort = np.sort(arr)
    i = 1
    all_sum = arr_sort.sum()
    if all_sum==0:
        ret = np.zeros_like(arr,dtype=bool)
        return ret
    while True:
        prop_sum = arr_sort[-i:].sum()
        if prop_sum/all_sum > prop:
            break
        i += 1 
    ret = np.zeros_like(arr,dtype=bool)
    ret[np.argsort(arr)[-i:]] = True
    return ret

## ----- Main class -----

class Spanve(object):
    def __init__(
        self, adata,
        spatial_info=None,
        neighbor_finder=None, # 'knn' or 'Delaunay'
        K: int = None, 
        hypoth_type: str = "nodist",        
        n_jobs: int = -1, 
        verbose:bool=False,
        **kwargs
    ) -> None:
        """spanve model.

        :param adata: AnnData object
        :type adata: AnnData
        :param spatial_info: spatial_infomation,default to adata.obsm['spatial'], defaults to None
        :type spatial_info: str;np.ndarray;pd.DataFrame, optional
        :param neighbor_finder: neighbor finder,could be 'knn' or 'Delaunay', defaults to None
        :type neighbor_finder: str, optional
        :param hypoth_type: distribution hypoth,could be 'nodist' or 'possion', defaults to "nodist"
        :type hypoth_type: str, optional
        :param n_jobs: number of paralle workers, defaults to -1
        :type n_jobs: int, optional
        :param verbose: verbose, defaults to False
        :type verbose: bool, optional
        """
        super().__init__()
        n_genes = adata.shape[1]
        sc.pp.filter_genes(adata,min_counts=1)
        if adata.shape[1] < n_genes:
            print(f'Filter genes with min_counts=1, {n_genes-adata.shape[1]} genes removed.')
        self.adata = adata
        self.K = max(K if K is not None else self.adata.shape[0]//100,5)
        self.n_jobs = n_jobs
        self.hypoth_type = hypoth_type
        self.verbose = verbose
        
        if neighbor_finder is None:
            if self.adata.shape[0] < 10000:
                self.neighbor_finder = "knn"
            else:
                self.neighbor_finder = "Delaunay"
        else:
            self.neighbor_finder = neighbor_finder

        X = adata.X.astype(int)
        if issparse(X):
            X = X.toarray()
        if spatial_info is None:
            assert 'spatial' in adata.obsm.keys(), "'spatial' is not in obsm keys, try set param `spatial_info`" 
            L = adata.obsm["spatial"]
        elif type(spatial_info) == str:
            L = adata.obsm[spatial_info]
        elif type(spatial_info) == np.ndarray:
            L = spatial_info
        elif type(spatial_info) == pd.DataFrame:
            L = spatial_info.loc[adata.obs_names,:].values
        else:
            raise TypeError(f'spatial_info is not valid. Now get type {type(spatial_info)}; spatial_info can be str[key of obsm], numpy.ndarry and pd.DataFrame.')

        self.X = X
        self.L = L

        for k, v in kwargs.items():
            setattr(self, k, v)

        self.__input_check(verbose=verbose)

    def __input_check(self,verbose=False):
        if not hasattr(self.adata, "X"):
            raise ValueError("adata.X is not defined")
        
        assert self.neighbor_finder in ["knn","Delaunay"], f"neighbor_finder should be 'knn' or 'Delaunay', now get {self.neighbor_finder}"
        assert self.hypoth_type in ["nodist","possion"], f"hypoth_type should be 'nodist' or 'possion', now get {self.hypoth_type}"
        assert self.X.shape[0] == self.L.shape[0], f"expression data is not consistent with spatial data, now get {self.X.shape[0]} != {self.L.shape[0]}"

        if self.adata.X.dtype not in [np.int64, np.int32, np.int16, np.int8, np.int0]:
            warnings.warn("""
            WARNNING: X must be an int matrix; 
            ----------------------------------------
            Will NOT automatically convert to int. Inputs can be Raw Counts or use `adata_preprocess_int` to get a normalized data with int dtype. """
            )

    def spatial_coexp(self,search_space,groupby=None,verbose=False):
        
        def spatial_coexp_single(x,y):
            sample_corr = (x - x.mean()) * (y - y.mean()) / (x.std() * y.std()+1e-7)
            return sample_corr.astype(int)
        
        def spatial_coexp_group(adata):
            newdf = pd.DataFrame(
                spatial_coexp_single(adata[:,list(var1)].X.toarray(),adata[:,list(var2)].X.toarray()),
                index = adata.obs_names,
                columns = [f"{i}~{j}" for i,j in search_space]
            )

            newad = sc.AnnData(newdf,obsm={'spatial':adata.obsm['spatial']},dtype=int)
            return newad
        
        adata = self.adata.copy()
        var1,var2 = zip(*search_space)
        
        if groupby is not None:
            assert groupby in adata.obs.columns, "groupby should be obs columns and should be categories." 
            groups = adata.obs_vector(groupby)
            n_groups = np.unique(groups).size
            newads = []
            
            bar = tqdm(total=n_groups,disable=not verbose, desc = f'There are {n_groups} groups. Will cal coexp strength separately.') 
            for g in np.unique(groups):
                adata_ = adata[groups==g,:]
                newad_ = spatial_coexp_group(adata_)
                newads.append(newad_)
                bar.update(1)
            newad = sc.concat(newads)
            bar.close()
        else:
            newad = spatial_coexp_group(adata)
        
        sc.pp.filter_genes(newad,min_counts=1)
        # show the number of genes pair filtered
        print("Delete gene pair with pseduo counts < 1:",newad.shape[1]-self.adata.shape[1],"pairs", "\nRemaining:",self.adata.shape[1],"pairs")
        self.adata = newad
        self.X = newad.X

    def possion_hypoth(self, X, verbose=False):
        overall_max = X.max(axis=0)
        n_features = X.shape[1]
        lams = np.std(X, axis=0) ** 2
        overall_dist = [
            Distribution(
                value=np.arange(0, overall_max[i] + 1),
                prob=ASP_pdf(np.arange(0, overall_max[i] + 1), lams[i]),
            )
            for i in trange(n_features, desc="#1 Expected Dist within Possion Hypoth", disable=not verbose)
        ]

        return overall_dist, overall_max

    def nodist_hypoth(self, X, verbose=False):
        n_features = X.shape[1]
        overall_dist = [
            Distribution().from_obs(obs=X[:, i]).dist_abs_subtract()
            for i in trange(n_features, desc="#1 Expected Dist within Nodist Hypoth", disable=not verbose)
        ]
        overall_max = X.max(axis=0)
        return overall_dist, overall_max

    def ent2gtest(self, Ents, ddof=0):

        df = self.overall_max
        pvals = chi2.sf(2 * len(self.adata) * Ents, df - ddof)
        pvals[np.isnan(pvals)] = 1
        rejects, fdrs, _1, _2 = multipletests(pvals, method="fdr_bh")
        return {"pvals": pvals, "rejects": rejects, "fdrs": fdrs}

    def finding_spatial_neibors(self, K,finder=None):
        finder = self.neighbor_finder if finder is None else finder
        nbr = NearestNeighbors(n_neighbors=K)
        nbr.fit(self.L)
        if finder =='knn':
            graph = nbr.kneighbors_graph()
            diag_mask = ~np.eye(*graph.shape).astype(bool)
            nbr_indices = np.where((graph == 1).todense() & diag_mask)
            return nbr_indices
        elif finder =='Delaunay':
            tri = scipy.spatial.Delaunay(self.L)
            nbr_idx1 = np.zeros((0), dtype=int)
            nbr_idx2 = np.zeros((0), dtype=int)
            for i,j in combinations(range(tri.simplices.shape[1]),2):
                nbr_idx1 = np.append(nbr_idx1,tri.simplices[:,i])
                nbr_idx2 = np.append(nbr_idx2,tri.simplices[:,j])
            return (nbr_idx1,nbr_idx2)

    def _AbsSubstract(self,X,indices,verbose):
        def computed_r(i):
            r = np.abs(X[indices[0], i] - X[indices[1], i])
            return np.unique(r, return_counts=True)

        Rs = Parallel(n_jobs=self.n_jobs)(
            delayed(computed_r)(i) for i in trange(X.shape[1], desc="#3 Computing Absolute Substract Value",disable=not verbose)
        )

        return Rs

    def fit(self, verbose=None, force_reject = False):
        # count time
        start = time.time()
        if verbose is None:
            verbose = self.verbose
        X = self.X
        n_features = X.shape[1]

        if self.hypoth_type == "possion":
            overall_dist, overall_max = self.possion_hypoth(X, verbose=verbose)
        elif self.hypoth_type == "nodist":
            overall_dist, overall_max = self.nodist_hypoth(X, verbose=verbose)
        else:
            raise ValueError("Unknown hypothesis type")

        self.overall_dist = overall_dist
        self.overall_max = overall_max
        # finding nearest k neighbors of each sample
        # from graph to get indices: recoder the index where the graph is 1

        indices = self.finding_spatial_neibors(K=self.K)
        self.nbr_indices = indices

        if verbose:
            print("#2 Nearest Neighbors Found")

        Rs = self._AbsSubstract(X,indices,verbose)

        def computed_G(i):
            ind, counts = Rs[i]
            obs_dist = Distribution(value=ind, prob=counts / counts.sum())
            inds = np.arange(overall_max[i] + 1)
            x = obs_dist[inds]
            y = overall_dist[i][inds]
            ent = entropy(x, y)
            return ent

        Ents = np.array([computed_G(i) for i in range(n_features)])
        Ents[np.isnan(Ents)] = np.inf

        self.ent = Ents
        if verbose:
            print("#4 Entropy Calculated")

        gtest_result = self.ent2gtest(Ents)

        for k, v in gtest_result.items():
            setattr(self, k, v)
        
        if self.rejects.sum() < 1 and force_reject:
            warnings.warn(
            """
            WARNNING: little significant features found.
            -----------------------------------------------
            Adjusted `rejetcs` by 20-80 rule. You can still see fdrs in attribute `self.fdrs`.
            Or try to change params `neighbor_finder`; Or try to set a proper K (recommend to `int(0.1*n_cells/n_clusters)`). 
            Or number of observation are too small.
            """)
            self.rejects = persent_select(self.ent)
        
        if verbose:
            print("#5 G-test Performed")
        result_df = pd.DataFrame(
            {
                "ent": self.ent,
                "pvals": self.pvals,
                "rejects": self.rejects,
                "fdrs": self.fdrs,
                "max_expr": self.overall_max,
            },
            index= self.adata.var_names,
        )
        self.result_df = result_df
        self.adata.var['spanve_spatial_features'] = result_df['rejects']
        self.adata.uns['spanve_running_parmas'] = {
            key : getattr(self,key) for key in ['K','hypoth_type','neighbor_finder','n_jobs','verbose']
            }
        if verbose:
            print("Write results to adata.var['spanve_spatial_features']")
            print(f"#--- Done, using time {time.time()-start:.2f} sec ---#")
        
        return self
    
    def genelayer_graph(self,alpha=0.05,select=None,K=None,verbose=None):
        if verbose is None:
            verbose = self.verbose
        select = select if select is not None else self.rejects
        n_samples = self.X.shape[0]

        # --- neighborhood---
        if K is None or K == self.K:
            K = self.K
            nbr_indices = self.nbr_indices
        else:
            nbr_indices = self.finding_spatial_neibors(K)

        # --- computing ---

        graph = csr_array((n_samples,n_samples),dtype=int)
        for i in tqdm(np.where(select)[0],disable=not verbose,desc='generate graph from spatial genes'):
            thres = self.overall_dist[i].isf(n=self.overall_max[i],alpha=alpha)
            x = self.X[:,i]
            osx = np.abs(x[nbr_indices[0]] - x[nbr_indices[1]])
            idx = np.where(osx >= thres)[0]
            cell_id0 = nbr_indices[0][idx]
            cell_id1 = nbr_indices[1][idx]
            graph += csr_array((np.ones(len(cell_id0)),(cell_id0,cell_id1)),shape=(n_samples,n_samples),dtype=int)
        graph = graph + graph.T
        return graph

    def impute_from_graph(
        self,X,
        n_circle=2,
        graph=None,
        verbose=None,
        ):
        assert n_circle >= 0 and isinstance(n_circle,int), 'n_circle must be a positive integer'
        if n_circle == 0:
            return X
        X = self.X[:,self.rejects] if X is None else X
        n_samples = X.shape[0]
        verbose = self.verbose if verbose is None else verbose
        if verbose:
            print('Impute data, there are',n_circle,'circles')
        graph = self.genelayer_graph(alpha=0.05,select=self.rejects,K=self.K,verbose=verbose) if graph is None else graph

        assert graph.shape == (n_samples,n_samples)

        graph = csr_array(np.eye(n_samples) + graph / (1+graph.sum(axis=0)))
        graph_imputed_X =( X.T @ graph ).T
        np.nan_to_num(graph_imputed_X,copy=False)
        
        if n_circle==1:
            return graph_imputed_X
        else:
            return self.impute_from_graph(X=graph_imputed_X,graph=graph,verbose=verbose,n_circle=n_circle-1)

    def plot_spatial(self,value,anndata=None,ax=None,):
        if anndata is None:
            anndata = self.adata
        if ax is None:
            fig,ax = plt.subplots()
        if value.dtype=='object':
            value = LabelEncoder().fit_transform(value)
        spatial_info = anndata.obsm['spatial']
        assert spatial_info.shape[0]==value.shape[0]
        ax.scatter(
            spatial_info[:,0],
            spatial_info[:,1],
            c=value,
            cmap='viridis',
            s=5,
            )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')
        # set color bar
        plt.colorbar(
            mappable= ax.collections[0],
            ax=ax)
        return ax

    def save(self,path,format='df'):
        warnings.warn('This function will not save input data. If you want to save input data, please use `anndata.AnnData.write_h5ad` function.')
        if format == 'df':
            self.result_df.to_csv(path,index=True)
        elif format == 'pickle':
            attr = self.__dict__.copy()
            # delete adata and X
            del attr['adata']
            del attr['X']
            del attr['L']

            with open(path,'wb') as f:
                pickle.dump(attr,f)
        else:
            raise ValueError('format must be df or pickle')

    def load(self,obj,verbose=True):
        if type(obj) == str:
            assert os.path.exists(obj), f'file {obj} not found'
            if obj.endswith('.pkl'):
                with open(obj,'rb') as f:
                    attr = pickle.load(f)
                self.load_dict(attr,verbose=verbose)
            elif obj.endswith('.csv'):
                df = pd.read_csv(obj,index_col=0)
                self.load_df(df, verbose)
        elif type(obj) == pd.DataFrame:
            self.load_df(obj,verbose)
        elif type(obj) == dict:
            self.load_dict(obj,verbose)

    def load_df(self, df, verbose):
        self.result_df = df
        self.rejects = self.result_df['rejects'].values
        self.fdrs = self.result_df['fdrs'].values
        self.pvals = self.result_df['pvals'].values
        self.overall_max = self.result_df['max_expr'].values
        self.ent = self.result_df['ent'].values

        self.overall_dist,self.overall_max = self.nodist_hypoth(self.X, verbose=verbose)
        self.nbr_indices = self.finding_spatial_neibors(self.K)

    def load_dict(self,attr,verbose=True):

        def print_attr(name,a,b=None):
            if not verbose:
                return
            if type(a) in [int, str, float, bool]:
                if b is None:
                    print(f'Load {name}: {a}')
                else:
                    print(f'Load {name}: {a} -> {b}')
            elif type(a) in [list, tuple]:
                print(f'Load {name}: {type(a)} with length {len(a)}')
            
            elif type(a) in [np.ndarray,csr_array]:
                print(f'Load {name}: {type(a)} with shape {a.shape}')
            else:
                print(f'Load {name}: {type(a)}')

        # verbose of the changed attributes
        for k in attr.keys():
            if k in self.__dict__:
                print_attr(k,self.__dict__[k],attr[k])
                self.__dict__[k] = attr[k]
            else:
                self.__dict__[k] = attr[k]
                print_attr(k,attr[k])

class Spanve_gpu(Spanve):
    def __init__(
        self,
        adata,
        K: int = None,
        device: int = 1,
        batch_size: int = 1024,
        hypoth_type: str = "nodist",
        neighbor_finder:str="knn", # or 'Delaunay'
        verbose=True):
        super().__init__(adata=adata,K=K,hypoth_type=hypoth_type,neighbor_finder=neighbor_finder,verbose=verbose)
        self.device = device
        self.batch_size = batch_size
    
    def _AbsSubstract(self,X,indices,verbose):
        try:
            import cupy as cp
            print(f'using cupy {cp.__version__} with {cp.cuda.Device(self.device).use()}')
        except:
            print(f'gpu is supportted by cupy package, follow the instruction (https://docs.cupy.dev/en/stable/install.html) to install cupy and set correct device id(now get {self.device}).')
            raise
        batch_size = self.batch_size
        n_features = X.shape[1]
        n_batches = int(np.ceil(n_features / batch_size))
        X = cp.array(X)
        indices = cp.array(indices)
        Rs = []

        def cpunique(x):
            y = cp.unique(x,return_counts=True)
            return y[0].get(),y[1].get()

        for i in tqdm(range(n_batches),disable=not verbose,desc=f'#3 Computing Absolute Substract Value(batch={batch_size})'):
            start = i * batch_size
            end = min((i+1) * batch_size,n_features)
            X_batch = X[:,start:end]
            Rs_batch = cp.abs(X_batch[indices[0],:] - X_batch[indices[1],:])
            Rs.extend([cpunique(Rs_batch[:,ii]) for ii in range(Rs_batch.shape[1])])
        return Rs


## ----- Cluster -----
# Automatic determine the number of clusters

class AutoCluster():
    def __init__(self,criteria='inertia',init_k = 3,max_k=10) -> None:
        self.criteria = criteria
        self.max_k = max_k
        self.init_k = init_k

        criteria_funcs = {
            'bic': self.bic,
            'inertia': self.inertia,
            'sh': sil,
            'ch': chs,
            'db': dbs,
        }
        self.criteria_func = criteria_funcs[criteria]

    def bic(self,model,X,labels):
        n_params = len(np.unique(labels)) * (X.shape[1] + 1)
        return -2 * model.score(X) * X.shape[0] + n_params * np.log(X.shape[0])

    def inertia(self,model,X=None,labels=None):
        return model.inertia_

    def fit(self,X,model=None,verbose=False,**kwargs):
        if model is None and X.shape[0] > 10000:
            self.model = MiniBatchKMeans
        elif model is None and X.shape[0] <= 10000:
            self.model = KMeans
        else:
            self.model = model
        if verbose:
            print(f'Sample size: {X.shape[0]}, using model: {self.model.__name__}')

        self.scores = []
        for k in trange(self.init_k,self.max_k+1,disable=not verbose, desc='finding best cluster number'):
            cluster = self.model(n_clusters=k,**kwargs)
            cluster.fit(X)
            if self.criteria in ['bic','inertia']:
                self.scores.append(
                    self.criteria_func(
                        model=cluster,labels=cluster.labels_,X=X))
            else:
                self.scores.append(
                    self.criteria_func(
                        X=X,labels=cluster.labels_))
        self.scores = np.array(self.scores)
        if self.criteria in ['bic','inertia']:
            self.best_k = elbow(self.scores) + self.init_k
        else:
            self.best_k = np.argmax(self.scores) + self.init_k

    def predict(self,X,**kwargs):
        cluster = self.model(n_clusters=self.best_k,**kwargs)
        cluster.fit(X)
        return cluster.labels_
    
    def fit_predict(self,X,model=None,verbose=False,**kwargs):
        self.fit(X,model=model,verbose=verbose,**kwargs)
        if verbose:
            print(f'Best k: {self.best_k}, Now predicting')
        return self.predict(X,**kwargs)

    def plot_elbow(self,ax=None):
        if ax is None:
            fig,ax = plt.subplots()
        ax.plot(range(self.init_k,self.max_k+1),self.scores)
        # show best k
        ax.plot(self.best_k,self.scores[self.best_k-self.init_k],'o',color='red')
        ax.text(self.best_k,self.scores[self.best_k-self.init_k],f'Best k: {self.best_k}')
        ax.set_xlabel('k')
        ax.set_ylabel(self.criteria)
        ax.set_title(f'Elbow for {self.criteria}')
        return ax