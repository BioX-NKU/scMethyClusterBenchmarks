import pandas as pd
import numpy as np
import glob

outdir = "./config/"
df_clustering = pd.DataFrame()
cell_files = glob.glob("../input/cell_files/*")
num_cells = len(cell_files)

#####Basic part
df = pd.read_csv(outdir + "runs_epiclomal/0_0.95_10000/result_basic/DIC_LINE_ELBOW_gainthr0.02_0.02/all_results_bestrun_basic.tsv", sep='\t', header=0)
file_path = df.loc[0, 'best_cluster']
file_path = file_path.replace("../", "./config/", 1)
print(file_path)
clusters_data = pd.read_csv(file_path, sep='\t', header=0, compression='gzip')
clusters_data = clusters_data.iloc[0:num_cells, 1:].values
pred = [None] * num_cells
for i in range(0, num_cells):
    Max = 0
    Maxj = 0
    for j in range(0, clusters_data.shape[1]):
        if clusters_data[i][j] > Max:
            Max = clusters_data[i][j]
            Maxj = j
    pred[i] = Maxj
pred = np.array(pred)
df_clustering['EpiclomalBasic'] = pred
print(pred)
true_label = pd.read_csv("../input/true_clone_membership.txt.gz", sep='\t', header=0, compression='gzip')
true_label = true_label.iloc[0:num_cells, 1].values
print(true_label)


from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import completeness_score

def cluster_metrics_basic(target, pred):
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
    df.to_csv('./result/Evaluation_EpiclomalBasic.csv', sep='\t', header=True, index=False)

def cluster_metrics_region(target, pred):
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
    df.to_csv('./result/Evaluation_EpiclomalRegion.csv', sep='\t', header=True, index=False)

print("Epiclomal basic:")
cluster_metrics_basic(true_label, pred)

#####Region part

df = pd.read_csv(outdir + "runs_epiclomal/0_0.95_10000/result_region/DIC_LINE_ELBOW_gainthr0.02_0.02/all_results_bestrun_region.tsv", sep='\t', header=0)
file_path = df.loc[0, 'best_cluster']
file_path = file_path.replace("../", "./config/", 1)
print(file_path)
clusters_data = pd.read_csv(file_path, sep='\t', header=0, compression='gzip')
clusters_data = clusters_data.iloc[0:num_cells, 1:].values
pred = [None] * num_cells
for i in range(0, num_cells):
    Max = 0
    Maxj = 0
    for j in range(0, clusters_data.shape[1]):
        if clusters_data[i][j] > Max:
            Max = clusters_data[i][j]
            Maxj = j
    pred[i] = Maxj
pred = np.array(pred)
df_clustering['EpiclomalRegion'] = pred
print(pred)
true_label = pd.read_csv("../input/true_clone_membership.txt.gz", sep='\t', header=0, compression='gzip')
true_label = true_label.iloc[0:num_cells, 1].values
print(true_label)
print("Epiclomal region:")
cluster_metrics_region(true_label, pred)

df_clustering.to_csv('./result/Epiclomal.csv', header=True, index=False)

