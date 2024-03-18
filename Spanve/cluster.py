import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.metrics import calinski_harabasz_score as chs
from sklearn.metrics import davies_bouldin_score as dbs
from sklearn.metrics import silhouette_score as sil
from tqdm import trange
from .utils import elbow

class AutoCluster():
    def __init__(self,criteria='inertia',init_k = 3,max_k=10) -> None:
        
        self.max_k = max_k
        self.init_k = init_k
        criteria_funcs = {
            'bic': self.bic,
            'inertia': self.inertia,
            'sh': sil,
            'ch': chs,
            'db': dbs,
        }
        self.criteria = criteria

        if isinstance(criteria,str):
            self.criteria_func = criteria_funcs[criteria]
        else:
            self.criteria_func = criteria

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
            ax = plt.gca()
        ax.plot(range(self.init_k,self.max_k+1),self.scores)
        # show best k
        ax.plot(self.best_k,self.scores[self.best_k-self.init_k],'o',color='red')
        ax.text(self.best_k,self.scores[self.best_k-self.init_k],f'Best k: {self.best_k}')
        ax.set_xlabel('k')
        ax.set_ylabel(self.criteria)
        ax.set_title(f'Elbow for {self.criteria}')
        return ax