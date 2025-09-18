import pandas as pd
import numpy as np

true_label = pd.read_csv('../input/true_clone_membership.txt.gz', sep='\t', header=0, index_col=None, compression='gzip').iloc[:, 1].values
print(true_label)
num_cells = len(true_label)

df = pd.read_csv("./result/Melissa.tsv", sep='\t', header=0)
df = df.iloc[0:num_cells, :].values
pred = [None] * num_cells
for i in range(0, num_cells):
    Max = 0
    Maxj = 0
    for j in range(0, df.shape[1]):
        if df[i][j] > Max:
            Max = df[i][j]
            Maxj = j
    pred[i] = Maxj
pred = np.array(pred)
print(pred)

from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import completeness_score

def cluster_metrics(target, pred):
    target = np.array(target)
    pred = np.array(pred)

    ari = adjusted_rand_score(target, pred)
    ami = adjusted_mutual_info_score(target, pred)
    nmi = normalized_mutual_info_score(target, pred)
    fmi = fowlkes_mallows_score(target, pred)
    comp = completeness_score(target, pred)
    homo = homogeneity_score(target, pred)
    print(f"ari=%.6f" % (ari))
    print(f"ami=%.6f" % (ami))
    print(f"nmi=%.6f" % (nmi))
    print(f"fmi=%.6f" % (fmi))
    print(f"comp=%.6f" % (comp))
    print(f"homo=%.6f" % (homo))
    df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ari, ami, nmi, fmi, comp, homo]})
    df.to_csv('./result/Evaluation_Melissa.csv', sep='\t', header=True, index=False)

cluster_metrics(true_label, pred)
