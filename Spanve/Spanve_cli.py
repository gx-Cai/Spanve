## ----- Cli usage -----
import click
from . import Spanve,adata_preprocess,AutoCluster,adata_preprocess_int
import anndata
import time
import os
from sklearn.decomposition import PCA

__all__ = ['main']

@click.command()
@click.option('--input_file','-i',type=click.Path(exists=True),help='input anndata file(h5ad file.)')
@click.option('--running_mode','-r',type=click.STRING,default='f',help='running mode, default is f(c:cluster; i:imputation; f:fitting)')
@click.option('--save_path','-s',type=click.Path(),default=None,help='save path')
@click.option('--verbose','-v',type=click.BOOL,default=False,help='verbose')
@click.option('--njobs','-n',type=click.INT,default=-1,help='number of parallel jobs,default is -1(use all cpus)')
@click.option('--preprocessed','-p',type=click.BOOL,default=0,help='int preprocessed or not.')
def main(
    input_file,
    running_mode='f', # 'c:cluster; i:imputation; f:fitting'
    preprocessed=False,
    save_path=None,
    verbose=True,
    njobs=-1
):  
    if save_path is None:
        # file names
        save_path = 'spanve_outs_' + os.path.basename(input_file).split('.')[0]

    if not os.path.exists(save_path):
        os.mkdirs(save_path)
        print('save path created: ', save_path)

    adata = anndata.read_h5ad(input_file)
    if preprocessed:
        adata = adata_preprocess_int(adata)
    model = Spanve(adata,n_jobs=njobs)
    time_log = {}
    st = time.time()
    model.fit(verbose=verbose)
    et = time.time()
    time_log['fitting_time'] = et - st
    model.save(
        os.path.join(save_path,'spanve_model.csv'),
    )
    anndata.var[model.rejects].index.to_csv(os.path.join(save_path,'spatial_expressed_genes.csv'),index=False)
    if 'i' in running_mode:
        pre_adata = adata_preprocess(adata)
        st = time.time()
        X_ = model.impute_from_graph(X=pre_adata.X,verbose=verbose)
        et = time.time()
        time_log['imputation_time'] = et - st
        X_.tofile(os.path.join(save_path,'imputed_data.npy'))

    if 'c' in running_mode:
        assert 'i' in running_mode, "imputation must be performed before clustering"
        pca = PCA(n_components=50)
        X_ = pca.fit_transform(X_)
        st = time.time()
        cluster = AutoCluster(init_k=2)
        labels = cluster.fit_predict(X_[:,model.rejects],verbose=verbose)
        et = time.time()

        time_log['clustering_time'] = et - st

        labels.tofile(os.path.join(save_path,'cluster_labels.npy'))
        cluster.plot_elbow().get_figure().savefig(os.path.join(save_path,'cluster_elbow.png'))
    
    with open(os.path.join(save_path,'spanve_time_log.txt'),'w') as f:
        f.write(str(time_log))

if __name__ == '__main__':
    
    main()
