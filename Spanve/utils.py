from collections import defaultdict
from itertools import combinations_with_replacement

import numpy as np
import scanpy as sc


class Distribution(defaultdict):
    """Distribution used to model the gene expression.
    """
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

def adata_preprocess(anndata,copy=True):
    """Preprocess the AnnData object.

    Args:
        anndata (AnnData): AnnData object.
        copy (bool, optional): Defaults to True.

    Returns:
        AnnData: preprocessed AnnData object.
    """
    if copy:
        anndata = anndata.copy()

    sc.pp.normalize_per_cell(anndata, counts_per_cell_after=10000)
    sc.pp.log1p(anndata)
    sc.pp.scale(anndata)
    # anndata.X = (anndata.X - anndata.X.mean(0)) / anndata.X.std(0)
    return anndata

def adata_preprocess_int(anndata,eps = 1e-7,exclude_highly_expressed=True,copy=True):
    """Special preprocess the AnnData object to convert the data to int. This process will keep the median of gene expression after the normalization.

    Args:
        anndata (AnnData): AnnData object.
        eps (float, optional): Used to escape the devide 0 error. Defaults to 1e-7.
        exclude_highly_expressed (bool, optional): Exclude the HVG or not in `sc.pp.normalize_total`. Defaults to True.
        copy (bool, optional): Defaults to True.

    Returns:
        AnnData: preprocessed AnnData object.
    """
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
    """Elbow method to find the optimal number of clusters.

    Args:
        X (np.ndarray): Metrics array.

    Returns:
        int: The optimal number of clusters.
    """
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

def modularity(X, labels):
    """Modularity metrics.

    Args:
        X (np.ndarray): Affinity matrix. (n_sample, n_sample)
        labels (list like): Labels of the nodes. (n_sample, )

    Returns:
        float: Modularity value.
    """
    # convert the affinity matrix to an adjacency matrix
    adj = (X > 0).astype(int)
    m = adj.sum() / 2
    k = adj.sum(axis=1)
    # comms is a list of arrays, where each array contains the indices of the nodes in the community
    comms = [np.where(labels == i)[0] for i in np.unique(labels)]
    Q = 0
    for comm in comms:
        # get the submatrix of the adjacency matrix corresponding to the community
        sub_adj = adj[np.ix_(comm, comm)]
        # get the subvector of the degree vector corresponding to the community
        sub_k = k[comm]
        # compute the modularity contribution of the community
        Q += (sub_adj.sum() / (2 * m)) - ((sub_k.sum() / (2 * m)) ** 2)
    # return the modularity value
    return Q

