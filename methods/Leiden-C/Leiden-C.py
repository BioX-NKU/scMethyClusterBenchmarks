import episcanpy as epi
import numpy as np
import scanpy as sc
import sklearn
import pandas as pd
import sys

n_clusters = int(sys.argv[1])

def get_N_clusters(adata, n_cluster, cluster_method='louvain', range_min=0, range_max=3, max_steps=30, tolerance=0):
    """
    Tune the resolution parameter in clustering to make the number of clusters and the specified number as close as possible.
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    n_cluster
        Specified number of clusters.
    cluster_method
        Method (`louvain` or `leiden`) used for clustering. By default, cluster_method='louvain'.
    range_min
        Minimum clustering resolution for the binary search.
    range_max
        Maximum clustering resolution for the binary search.
    max_steps
        Maximum number of steps for the binary search.
    tolerance
        Tolerance of the difference between the number of clusters and the specified number.
    Returns
    -------
    adata
        AnnData object with clustering assignments in `adata.obs`:
        - `adata.obs['louvain']` - Louvain clustering assignments if `cluster_method='louvain'`.
        - `adata.obs['leiden']` - Leiden clustering assignments if `cluster_method='leiden'`.
    """ 
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        if cluster_method=='leiden':
            sc.tl.leiden(adata, resolution=this_resolution)
        if cluster_method=='louvain':
            sc.tl.louvain(adata, resolution=this_resolution)
        this_clusters = adata.obs[cluster_method].nunique()
        if this_clusters > n_cluster+tolerance:
            this_max = this_resolution
        elif this_clusters < n_cluster-tolerance:
            this_min = this_resolution
        else:
            print(this_step)
            print("Succeed to find %d clusters at resolution %.3f."%(n_cluster, this_resolution))
            return adata
        this_step += 1
    print('Cannot find the number of clusters.')
    print(this_step)
    return adata

adata = sc.read_h5ad('constructed.h5ad')
print(adata)

adata = adata[:, ~np.all(np.isnan(adata.X), axis=0)].copy()
# Step 2: Calculate the median for each column, ignoring NaNs
medians = np.nanmedian(adata.X, axis=0)
# Step 3: Replace NaNs with the column medians
for i in range(adata.X.shape[1]):
    adata.X[:, i] = np.where(np.isnan(adata.X[:, i]), medians[i], adata.X[:, i])
epi.pp.lazy(adata, nb_pcs=min(adata.n_obs - 1, 50), perplexity=min(adata.n_obs - 1, 30))

adata = get_N_clusters(adata, n_clusters, 'leiden')

cluster_key='cell_type'
label_key='leiden'
print('Leiden')
print(adata.obs[cluster_key])
print(adata.obs[label_key])
ARI = epi.tl.ARI(adata, cluster_key, label_key)
AMI = epi.tl.AMI(adata, cluster_key, label_key)
NMI = sklearn.metrics.normalized_mutual_info_score(adata.obs[cluster_key], adata.obs[label_key])
FMI = sklearn.metrics.fowlkes_mallows_score(adata.obs[cluster_key], adata.obs[label_key])
Comp = sklearn.metrics.completeness_score(adata.obs[cluster_key], adata.obs[label_key])
Homo = epi.tl.homogeneity(adata, cluster_key, label_key)
print('Leiden')
print('ARI:%.6f\tAMI:%.6f\tNMI:%.6f\tFMI:%.6f\tComp:%.6f\tHomoL%.6f'%(ARI, AMI, NMI, FMI, Comp, Homo))
df = pd.DataFrame({'Leiden-C': adata.obs[label_key]})
df.to_csv('./result/Leiden-C.csv', header=True, index=False)
df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ARI, AMI, NMI, FMI, Comp, Homo]})
df.to_csv('./result/Evaluation_Leiden-C.csv', sep='\t', header=True, index=False)
