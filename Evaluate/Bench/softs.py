import sys
import os
sys.path.append('../Softs')
import pandas as pd

def load_soft(soft):
    return SpatialSoft(soft)

class SpatialSoft():
    def __init__(self,soft) -> None:
        self.name = soft.lower()
        self.func = self.load_func()

    def __call__(self,data,saving_path,**kwargs):
        self.func(data=data,saving_path=saving_path,**kwargs)

    def load_func(self):
        if self.name == 'spanve-k':
            return lambda data,saving_path: run_spanve(data,saving_path,neighbor_finder='knn')
        elif self.name == 'spanve-d':
            return lambda data,saving_path: run_spanve(data,saving_path,neighbor_finder='Delaunay')
        elif self.name == 'p.spanve-k':
            return lambda data,saving_path: run_spanve(data,saving_path,preprocessed=True,neighbor_finder='knn')
        elif self.name == 'p.spanve-d':
            return lambda data,saving_path: run_spanve(data,saving_path,preprocessed=True,neighbor_finder='Delaunay')
        elif self.name == 'scgco':
            return run_scgco
        elif self.name == 'sepal':
            return run_sepal
        elif self.name == 'somde':
            return run_somde
        elif self.name == 'spatialde':
            return run_spatialde
        elif self.name == 'moran':
            return run_moran_test
        elif self.name == 'geary':
            return run_geary_test
        elif self.name == 'r-gitto-rank':
            return lambda data,saving_path: run_gitto(data,saving_path,method='rank')
        elif self.name == 'r-gitto-sirank':
            return lambda data,saving_path: run_gitto(data,saving_path,method='sirank')
        elif self.name == 'r-sparkx':
            return run_sparkx
        elif self.name == 'r-meringue':
            return run_meringue
        else:
            raise NotImplementedError

def run_spanve(data,saving_path,preprocessed=False,**kwargs):
    from Spanve import Spanve,adata_preprocess_int
    if preprocessed:
        data = adata_preprocess_int(data)
    model = Spanve(adata=data,verbose=False,**kwargs)
    model.fit()
    model.save(os.path.join(saving_path,'spanve.model.csv'))

def run_scgco(data,saving_path,**kwargs):
    # raise NotImplementedError('Bugs Meets When running scgco. Waiting for resonse in Github.')
    from scGCO import normalize_count_cellranger,create_graph_with_weight,multiGMM,identify_spatial_genes,write_result_to_csv
    locs = data.obsm['spatial']
    data_norm = normalize_count_cellranger(pd.DataFrame(data.X))
    exp= data_norm.iloc[:,0]
    cellGraph= create_graph_with_weight(locs, exp)
    
    gmmDict=multiGMM(data_norm)
    result_df= identify_spatial_genes(locs, data_norm, cellGraph ,gmmDict)
    write_result_to_csv(result_df,os.path.join(saving_path,'scgco.out.csv'))
    
def run_sepal(data,saving_path,data_dir,**kwargs):
    import sepal
    import sepal.datasets as d
    import sepal.models as m
    import sepal.utils as ut
    raw_data = d.RawData(data_dir,)
    data = m.UnstructuredData(raw_data)
    times = m.propagate(data,normalize = True,scale =True)
    times.sort_values(by='average').to_csv(os.path.join(saving_path,'sepal.score.csv'))
    
def run_somde(data,saving_path,**kwargs):
    from somde import SomNode
    som = SomNode(data.obsm['spatial'],20)
    ndf,ninfo = som.mtx(data.to_df().T)
    nres = som.norm()
    result, SVnum =som.run()
    result.to_csv(os.path.join(saving_path,'somde.out.csv'))
    
def run_spatialde(data,saving_path,**kwargs):
    import NaiveDE
    import SpatialDE
    sample_info = pd.DataFrame(index = data.obs_names)
    sample_info[['x','y']] = data.obsm['spatial']
    sample_info['total_counts'] = data.X.sum(axis=1)
    norm_expr = NaiveDE.stabilize(data.to_df().T).T
    resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T

    sample_resid_expr = resid_expr.sample(n=min(1000,data.shape[1]), axis=1, random_state=1)
    X = sample_info[['x', 'y']]
    results = SpatialDE.run(X, sample_resid_expr)

    results.sort_values(by='qval').to_csv(os.path.join(saving_path,'spatialde.out.csv'))
    
def run_moran_test(data,saving_path,**kwargs):
    from squidpy import gr
    gr.spatial_neighbors(data)
    gr.spatial_autocorr(data,mode='moran')
    
    data.uns['moranI'].to_csv(os.path.join(saving_path,'moranI.test.csv'))

def run_geary_test(data,saving_path,**kwargs):
    from squidpy import gr
    gr.spatial_neighbors(data)
    gr.spatial_autocorr(data,mode='geary')
    
    data.uns['gearyC'].to_csv(os.path.join(saving_path,'gearyC.test.csv'))
    
def run_gitto(data,saving_path,method='rank',**kwargs):
    import os
    command = f'Rscript ../Softs/run_gitto.R -d {data} -w {saving_path} -m "{method}"'
    print(command)
    os.system(command)
    
def run_sparkx(data,saving_path,**kwargs):
    import os
    command = f'Rscript ../Softs/run_sparkx.R -d {data} -w {saving_path}'
    print(command)
    os.system(command)

def run_meringue(data,saving_path,**kwargs):
    import os
    command = f'Rscript ../Softs/run_meringue.R -d {data} -w {saving_path}'
    print(command)
    os.system(command)