import click
import logging
import anndata
import os
import time
from softs import *
from softs import load_soft

class Evaluator():
    def __init__(self,dataset_dir,soft,saving_path,) -> None:
        self.dataset_dir = dataset_dir
        assert dataset_dir.endswith('h5ad')
        self.soft_name = soft
        self.saving_path = saving_path
        self.logger = logging.getLogger(f'{dataset_dir}_{soft}')
        self.logger.setLevel(logging.INFO)

        logger_path = os.path.join(saving_path,f'run.log')
        fh = logging.FileHandler(logger_path)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.soft = load_soft(soft)
        

    def load_data(self,data_path):
        return anndata.read_h5ad(data_path)
    
    def update_kwargs(self,kwargs):
        if self.soft_name=='sepal':
            kwargs.update(data_dir = self.dataset_dir)

    def evaluate_py(self,**kwargs):
        self.logger.log(msg=f'running {self.dataset_dir}',level=logging.INFO)
        data = self.load_data(self.dataset_dir)
        self.update_kwargs(kwargs)
        st = time.time()
        self.soft(
            data=data,saving_path=self.saving_path,
            **kwargs
        )
        et = time.time()
        self.logger.log(msg='time cost:{}'.format(et-st),level=logging.INFO)
        
        return self
    
    def evaluate_r(self,**kwargs):
        self.logger.log(msg=f'running {self.dataset_dir}',level=logging.INFO)
        self.update_kwargs(kwargs)
        
        self.soft(
            data=self.dataset_dir,saving_path=self.saving_path,
            **kwargs
        )
        
        return self
    
    
    def evaluate(self,**kwargs):
        if self.soft_name.startswith('r'):
            self.evaluate_r(**kwargs)
        else:
            self.evaluate_py(**kwargs)

@click.command()
@click.option('--dataset','-d',default=None)
@click.option('--soft','-s',default=None)
@click.option('--saving_path','-sp',default=None)
def main(dataset=None,soft=None,saving_path=None):
    evaluator = Evaluator(dataset,soft,saving_path).evaluate()

if __name__ == '__main__':
    main()