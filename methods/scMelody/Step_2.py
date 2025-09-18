from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
import numpy as np

###Calculate the basic cell partitions using spectral clustering for each basic similarity measure
###Input : the basic similarity matrices; the number of clusters
###Output : the basic cell partitions
def sc_pre(sm,C):
    sc = SpectralClustering(n_clusters=C, n_neighbors=2*C, affinity='precomputed',n_init=100,assign_labels='kmeans')
    out = sc.fit_predict(sm)
    return out

###Calculate the basic cell partitions using hierarchical clustering for each basic similarity measure
###Input : the basic similarity matrices; the number of clusters
###Output : the basic cell partitions
def hc_pre(sm,C):
    hc = AgglomerativeClustering(n_clusters=C, affinity='precomputed',linkage='complete')
    out = hc.fit_predict(1-sm)
    return out

###Calculate the weight of basic cell partitions based on the clustering separability
###Input : the basic similarity matrices and the respective inferred cell partitions
###Output : the weights for basic cell partitions
def sil(s1,s2,s3,c1,c2,c3):
    score1 = np.exp(metrics.silhouette_score(1-s1, c1, metric='precomputed'))
    score2 = np.exp(metrics.silhouette_score(1-s2, c2, metric='precomputed'))
    score3 = np.exp(metrics.silhouette_score(1-s3, c3, metric='precomputed'))
    sil_sum = score1+score2+score3
    w1 = score1/sil_sum
    w2 = score2/sil_sum
    w3 = score3/sil_sum
    return np.array([w1,w2,w3])

###Calculate the weight of basic cell partitions based on the clustering diversity
###Input : the basic cell partitions
###Output : the weights for basic cell partitions
def nmi(c1,c2,c3):
    p12 = metrics.normalized_mutual_info_score(c1,c2)
    p13 = metrics.normalized_mutual_info_score(c1,c3)
    p23 = metrics.normalized_mutual_info_score(c2,c3)
    p1 = np.exp(-1/2 *(p12+p13))
    p2 = np.exp(-1/2 *(p12+p23))
    p3 = np.exp(-1/2 *(p23+p13))
    pw_sum = p1+p2+p3
    w1 = p1/pw_sum
    w2 = p2/pw_sum
    w3 = p3/pw_sum
    return np.array([w1,w2,w3])

###Calculate the binary co-occurrence matrix to convert the basic cell cluster to cell-to-cell similarity measure
###Input : inferred cell clusters
###Output : the binary co-occurrence matrix
def COM(y):
    S = np.zeros((len(y), len(y)))
    for i in range(len(y)):
        for j in range(i, len(y)):
            if(y[i]==y[j]):
                S[i][j] = 1
            else:
                S[i][j] = 0
            S[j][i]=S[i][j]
    return S

###Calculate the reconstructed similarity matrix to generate final cell partitions
###Input : the basic similarity matrices; number of clusters
###Output : the reconstructed similarity matrix
def scM(s1,s2,s3,C):
    c1 = sc_pre(s1,C)
    c2 = sc_pre(s2,C)
    c3 = sc_pre(s3,C)
    co1= np.multiply(s1,COM(c1))
    co2= np.multiply(s2,COM(c2))
    co3= np.multiply(s3,COM(c3))
    cn1 = 0.5*(co1*sil(s1,s2,s3,c1,c2,c3)[0]+
               co2*sil(s1,s2,s3,c1,c2,c3)[1]+
               co3*sil(s1,s2,s3,c1,c2,c3)[2])
    cn2 = 0.5*(co1*nmi(c1,c2,c3)[0]+
               co2*nmi(c1,c2,c3)[1]+
               co3*nmi(c1,c2,c3)[2])
    return cn1+cn2

###Calculate the optimal number of clusters for the spectral clustering
###Input : possible maximum number of clusters; the basic similarity matrices
###Output : the number of clusters and the respective SI score
def find_kspcluster(k_max,sm1,sm2,sm3):
    sil_list = []
    for k in range(2,k_max+1):
        out1 = sc_pre(sm1,k)
        out2 = sc_pre(sm2,k)
        out3 = sc_pre(sm3,k)
        sil1 = metrics.silhouette_score(1-sm1, out1, metric='precomputed')
        sil2 = metrics.silhouette_score(1-sm2, out2, metric='precomputed')
        sil3 = metrics.silhouette_score(1-sm3, out3, metric='precomputed')
        sil_sum = sil1+sil2+sil3
        sil_list.append([k,sil_sum])    
    return(sil_list)

###Calculate the final optimal number of clusters for the HC algorithm
###Input : possible maximum number of clusters; the reconstructed similarity matrix
###Output : the number of clusters and the respective SI score
def find_kopcluster(k_max,co):
    sil_list = []
    for k in range(2,k_max+1):
        out = hc_pre(co,k)
        sil = metrics.silhouette_score(1-co, out, metric='precomputed')
        sil_list.append([k,sil])    
    return(sil_list)
