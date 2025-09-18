import episcanpy as epi
import numpy as np
import scanpy as sc
import sklearn
import pandas as pd

adata = sc.read_h5ad('constructed.h5ad')

adata = adata[:, ~np.all(np.isnan(adata.X), axis=0)].copy()
# Step 2: Calculate the median for each column, ignoring NaNs
medians = np.nanmedian(adata.X, axis=0)
# Step 3: Replace NaNs with the column medians
for i in range(adata.X.shape[1]):
    adata.X[:, i] = np.where(np.isnan(adata.X[:, i]), medians[i], adata.X[:, i])
epi.pp.lazy(adata, nb_pcs=min(adata.n_obs - 1, 50), perplexity=min(adata.n_obs - 1, 30))
epi.tl.leiden(adata)

cluster_key='cell_type'
label_key='leiden'
print('Leiden')
print(adata.obs[cluster_key])
print(adata.obs[label_key])
df_clustering = pd.DataFrame()
df_clustering['Leiden'] = adata.obs[label_key]
df_clustering.to_csv('./result/Leiden.csv', header=True, index=False)
ARI = epi.tl.ARI(adata, cluster_key, label_key)
AMI = epi.tl.AMI(adata, cluster_key, label_key)
NMI = sklearn.metrics.normalized_mutual_info_score(adata.obs[cluster_key], adata.obs[label_key])
FMI = sklearn.metrics.fowlkes_mallows_score(adata.obs[cluster_key], adata.obs[label_key])
Comp = sklearn.metrics.completeness_score(adata.obs[cluster_key], adata.obs[label_key])
Homo = epi.tl.homogeneity(adata, cluster_key, label_key)
print('Leiden')
print('ARI:%.6f\tAMI:%.6f\tNMI:%.6f\tFMI:%.6f\tComp:%.6f\tHomoL%.6f'%(ARI, AMI, NMI, FMI, Comp, Homo))
df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ARI, AMI, NMI, FMI, Comp, Homo]})
df.to_csv('./result/Evaluation_Leiden.csv', header=True, index=False, sep='\t')

