import numpy as np
import pandas as pd
import os
import scanpy as sc

class GeneSelector():
    def __init__(self,stdata, scdata, n_genes, cell_type_key = 'celltype'):
        genes_st = stdata.var_names.values
        genes_sc = scdata.var_names.values
        st_sc_intersect = list(set(genes_sc).intersection(set(genes_st)))
        self.valid_genes = set(st_sc_intersect)
        
        self.stdata = stdata
        self.scdata = self.process_scdata(scdata)
        self.n_genes = n_genes
        self.cell_type_key = cell_type_key
        
        
        self.FUNSTR = {
            'spanve': self.run_spanve,
            'seurat': lambda **kwargs: self.run_hvg(flavor = 'seurat',**kwargs),
            'seurat_v3': lambda **kwargs: self.run_hvg(flavor = 'seurat_v3',**kwargs),
            'cell_ranger': lambda **kwargs: self.run_hvg(flavor = 'cell_ranger',**kwargs),
            'rank_genes_group': self.run_rnk_from_sc,
            'moranI':self.run_moran_test,
            'spatialDE': self.run_spatialde
        }
        
    def process_scdata(self, scdata):
        scdata.layers['counts'] = scdata.X.copy()
        sc.pp.normalize_total(scdata)
        sc.pp.log1p(scdata)
        scdata.layers['log'] = scdata.X.copy()
        scdata.X = scdata.layers['counts'].copy()
        return scdata[:, list(self.valid_genes)].copy()
    
    def run_moran_test(self):
        from squidpy import gr
        data = self.stdata.copy()
        gr.spatial_neighbors(data)
        gr.spatial_autocorr(data,mode='moran')
        return data.uns['moranI'].sort_values(by='I', ascending=False).index[0:self.n_genes]

    def run_spatialde(self):
        import NaiveDE
        import SpatialDE
        
        data = self.stdata.copy()
        sample_info = pd.DataFrame(index = data.obs_names)
        sample_info[['x','y']] = data.obsm['spatial']
        sample_info['total_counts'] = data.X.sum(axis=1)
        norm_expr = NaiveDE.stabilize(data.to_df().T).T
        resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
        sample_resid_expr = resid_expr.sample(n=min(1000,data.shape[1]), axis=1, random_state=1)
        X = sample_info[['x', 'y']].values
        results = SpatialDE.run(X, sample_resid_expr)
        return results.sort_values(by='qval')['g'][0:self.n_genes]
    
    def run_spanve(self):
        from Spanve import adata_preprocess_int, Spanve
        stdata = adata_preprocess_int(self.stdata,copy=True)
        spanve = Spanve(stdata)
        spanve.fit(verbose=False)
        spanve_selected = spanve.result_df.sort_values(by='ent',ascending = False)[0:self.n_genes].index.tolist()
        return spanve_selected

    def run_hvg(self, flavor):
        if flavor in ['seurat', 'cell_ranger']:
            sc.pp.highly_variable_genes(self.scdata, flavor=flavor, layer = 'log', n_top_genes = self.n_genes)
        elif flavor == 'seurat_v3': 
            sc.pp.highly_variable_genes(self.scdata, flavor=flavor, layer = 'counts', n_top_genes = self.n_genes)
        else:
            raise
        selected = self.scdata.var[self.scdata.var['highly_variable'] == True].index.tolist()
        return selected
    
    def run_rnk_from_sc(self):
        sc.tl.rank_genes_groups(self.scdata, groupby=self.cell_type_key, use_raw=False, layer = 'log')
        scores = sc.get.rank_genes_groups_df(self.scdata, group=None).groupby('names').apply(lambda x: np.abs(x['scores']).sum()).sort_values(ascending=False)
        return scores.index[0:self.n_genes]
    
    def fit(self, method):
        if isinstance(method, list):
            selected = set()
            for m in method:
                select_sep = self.fit(m)
                selected |= set(select_sep)
        elif isinstance(method, str):
            func = self.FUNSTR[method]
            selected = set(func())
        selected = list(selected & self.valid_genes)
        self.selected = selected
        return selected
    
    def fit_transform(self, method):
        selected_genes = self.fit(method)
        return self.scdata[:, selected_genes].copy()
