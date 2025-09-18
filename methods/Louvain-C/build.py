import pandas as pd
import numpy as np
import episcanpy as epi
import scanpy as sc
import anndata as ad
import os
import sys
import gzip
import scipy.io as sio

MOUSE = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19','X', 'Y', 'M']

mm9_size = [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296, 15902555, 16299]

mm10_size = [195471971, 182113224, 160039680, 156508116,
              151834684, 149736546, 145441459, 129401213, 124595110,
              130694993, 122082543, 120129022, 120421639, 124902244,
              104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299]

HUMAN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','X', 'Y', 'M']

hg38_size = [248956422, 242193529, 198295559, 190214555,
              181538259, 170805979, 159345973, 145138636, 138394717,
              133797422, 135086622, 133275309, 114364328, 107043718,
              101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
              46709983, 50818468, 156040895, 57227415, 16569]

hg19_size = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 63025520,  59128983, 51304566, 48129895,155270560,59373566, 16571]


def read_meth_fileCG(sample_name, path, chromosome,col_dict ):
    # met = 5 refer to methylation level in bed.gz file
    """
    Read file from which you want to extract the methylation level and
    (assuming it is like the Ecker/Methylpy format) extract the number of
    methylated read and the total number of read for the cytosines covered and
    in the right genomic context (CG or CH)
    Parameters
    ----------
    sample_name:
        name of the file to read to extract key information.
    meth_type:
        CG, CH or not specified
    head: if there is header that you don't want to read. An annotation in the
        file you are reading. The default value is the Methylpy/Ecker header
    path: path of the to access the file of the sample you want to read. 
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes)
    """
    chrom=col_dict['chrom']; pos=col_dict['pos'];met=col_dict['met']; tot=col_dict['tot']
    reduced_cyt = {key: [] for key in chromosome} # to store cyt for each chrom (intermed output)
    with gzip.open(path+sample_name) as sample:
        for line in sample:
            line = str(line, encoding="utf-8").split('\t')
            print(line)
            if(line[chrom][3:] in chromosome):
                reduced_cyt[line[chrom][3:]].append((int(line[pos]), float(line[met]), int(line[tot])))
    return(reduced_cyt)
def methylation_level(reduced_cyt, feature, chromosome, threshold=1):
    """
    Measure the methylation for the feature you give as input using the reduce
    representation of the sample cytosine (output of read_methylation_file)
    Parameters
    ----------
    reduced_cyt: datatype that contained processed sample file. It only contains
        the cytosines that were in the genomic context you wanted to filter for.
        (output of read_methylation_file function).
    feature: the feature in the right datatype for which you want to determine
        the methylation level.
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes).
    """
    ## actually, to write sparse matrix I need a list, not a dictionary
    #meth_levels_bins = {key:[] for key in chromosome}
    meth_levels_bins = []

    for c in chromosome:
        meth_reads = np.zeros(len(feature[c]))
        tot_reads = np.zeros(len(feature[c]))
        nb_cyt = np.zeros(len(feature[c]))
        cytosines = reduced_cyt[c] 
        i = 0
        for j in range(len(feature[c])): # for every bins in a given chrom
            meth_reads = 0.0
            tot_reads = 0.0
            nb_cyt = 0
            # I am skipping cytosine that are before the beginning 
            # of my current bin. 
            while (i < len(cytosines)) and cytosines[i][0] < feature[c][j][0]:
                i += 1
            # Once I got one cytosine that is after the beginning of my feature 
            # I need to check if this feature is within the enhancer limits
            # so if the position of the cytosine is not > to the end of the feature
            if i<len(cytosines) and cytosines[i][0] <= feature[c][j][1]:
                meth_reads += cytosines[i][-2] # meth cyt read
                tot_reads += cytosines[i][-1] # tot cyt read
                nb_cyt += 1 # nb of cyt

            # check if the next cytosine fall into the current feature is important
            # to this end, I have another pointer/iterator k. 
            # at the next feature I will have to check from i but for the current 
            # feature I need to check the next cytosine and I use the variable k for 
            # this.
            k = i+1
            while k < len(cytosines) and cytosines[k][0] <= feature[c][j][1]:
                meth_reads += cytosines[k][-2] # meth cyt read
                tot_reads += cytosines[k][-1] # tot cyt read
                nb_cyt += 1  # nb of cyt
                k += 1
            ## actually, to write sparse matrix I need a list, not a dictionary
            if nb_cyt >= threshold:
                #meth_levels_bins[c].append(format(meth_reads/tot_reads, '.3f'))
                meth_levels_bins.append(format(meth_reads/tot_reads, '.3f'))
            else:
                #meth_levels_bins[c].append(np.nan)
                meth_levels_bins.append(np.nan)


    return(meth_levels_bins)
def build_count_mtx(cells, annotation,col_dict, path="", output_file=None, writing_option="a",
                    meth_context="CG", chromosome=HUMAN, feature_names=None,
                   threshold=1, ct_mtx=None, sparse=False):
    """
    Build methylation count matrix for a given annotation.
    It either write the count matrix (if given an output file) or return it as a variable (numpy matrix). 
    If you want to add cells to an already existing matrix (with the same annotations),  you put the initial matrix as ct_mtx or you specify the matrix to write +
    writing option = a
        
    if you want to write down the matrix as a sparse matrix you have to specify it (not implented yet)
        
    I need to pay attention to where I am writing the output file.
        
    Also, verbosity..
        
    Pay attention, it does not average variables. If you want to process many small features such as
    tfbs, we advise to use the dedicated function.
        
    Parameters
    ----------
    cells:
        list of the file names to read to build the count matrix.
    annotation:
        loaded annotation to use to build the count matrix
        'str' or 'list' depending of the number of matrices to build
    path:
        path to the input data. 
    output_file:
        name files to write. 'str' or 'list' depending of the number of matrices to build
    writing_option:
        either 'w' if you want to erase potentialy already existing file or 'a' to append.
        'str' or 'list' if you have a list of matrices and the writing options are differents
    meth_context:
        read methylation in 'CG' of 'CH' context
    chromosome:
        'MOUSE' and 'HUMAN' (without mitochondrial genome) or list with chromosomes.
    feature_names:
        If you want to write down the name of the annotation features. 
        'Int' (or 'list' if you have multiple annotations)
    thereshold:
        the minimum of cytosines covered per annotation to calculate a methylation level.
        default=1 'Int' (or 'list' if you have multiple annotations with different thresholds)
    ct_mtx:
        numpy matrix containing the same set of annotations and for which you want to append.
        default: None
    sparse:
        Boolean, writing option as a normal or sparse matrix.
        default: False
    
    """
    #verbosity
    i = 0
    
    #################################
    if type(annotation) != list:
        annotation = [annotation]
        output_file = [output_file]
        ct_mtx = [ct_mtx]
        feature_names = [feature_names]
        
    nb_annotation = len(annotation)
    
    if type(writing_option) != list:
        writing_option = [writing_option for x in range(nb_annotation)]
    if type(threshold) != list:
        threshold = [threshold for x in range(nb_annotation)]
    if (output_file != None):
        if (type(output_file) != list):
            output_file = [output_file]
        
    #################################
    for cell in cells:
        #verbosity
        print(i, cell)
        i += 1
        
        # read the file to extract cytosines in the right context
        if meth_context == 'CG':
            tmp_file = read_meth_fileCG(cell, path, chromosome,col_dict)
        elif meth_context == 'CH':
            tmp_file = read_meth_fileCH(cell, path, chromosome)
        else:
            break
        
        # build the cell vector for the count matrix at every set of annotations
        for index_annot in range(nb_annotation):
            meth_level_annot = methylation_level(tmp_file, annotation[index_annot], chromosome, threshold[index_annot])
            if type(output_file) == list:
                write_methlevel(meth_level_annot, output_file[index_annot], cell, writing_option[index_annot], feature_names[index_annot])
            else:
                if ct_mtx == None:
                    ct_mtx = [np.matrix(meth_level_annot)]
                    #ct_mtx[index_annot] = np.matrix(meth_level_annot)
                elif index_annot>=len(ct_mtx):
                    ct_mtx.append(np.matrix(meth_level_annot))
                else:
                    ct_mtx[index_annot] = np.vstack([ct_mtx[index_annot], meth_level_annot])
                
    if ct_mtx != None:
        return(ct_mtx)
    else:
        return()



def replace_bedgz(name):
    #replace_name = name.replace('.single.CpG.txt.gz','')
    replace_name = name.replace('.txt.gz', '')
    return replace_name

def count_matrix(bed_path,output_path,col_dict,species,windows_width=10000,mm_refgenome=None,hg_refgenome=None):
    
    #定义参考基因组size
    if mm_refgenome =='mm9':
        refgenome_size = mm9_size
    elif mm_refgenome =='mm10' :
        refgenome_size = mm10_size
        
    if hg_refgenome =='hg19':
        refgenome_size = hg19_size
    elif hg_refgenome =='hg38':
        refgenome_size = hg38_size
        
    
    # 读取metadata
    #metadata = pd.read_csv('/home/stclassify/data/lisiyu/methylation/metadata/%s.csv'%dataset)

    # 获取bed.gz列表        
    bed_path = bed_path
    bed_dir = os.listdir(bed_path)
    bed_dir = np.sort(bed_dir)
    print(bed_dir)
    
    if '.ipynb_checkpoints' in bed_dir:
        bed_dir.remove('.ipynb_checkpoints')
    # 分类并构建count_matrix
    
    
    windows = epi.ct.make_windows(windows_width, chromosomes=species, chromosome_sizes=refgenome_size)
    w_names = epi.ct.name_features(windows) # extract name of the features
    w_mtx = build_count_mtx(cells=bed_dir,
                               annotation=[windows],
                               col_dict=col_dict,
                               path=bed_path,
                               chromosome=species,
                               output_file=None,
                               meth_context='CG',
                               threshold=[1])# minimum number of cytosine/reads to have at any given feature to 
                                # not consider the feature to have a missing methylation level


    cell_name = list(map(replace_bedgz,bed_dir))
#     metadata.index = cell_name
#     sio.savemat('/home/stclassify/data/lisiyu/methylation/23_science_data/h5ad/%s.mat', w_mtx[0])
    w_mtx[0]=np.array(w_mtx[0]).astype('float64')
    adata_w = ad.AnnData(w_mtx[0],  var=pd.DataFrame(index=w_names),obs=pd.DataFrame(index=cell_name))
    cell_type = pd.read_csv('../input/true_clone_membership.txt.gz', index_col=None, header=0, sep='\t')
    adata_w.obs['cell_type'] = cell_type.iloc[:, 1].values
    
    
    adata_w.write_h5ad(output_path)

    del w_mtx

bed_path = '../input/cell_files/'
output_path = 'constructed.h5ad'
col_dict ={'chrom':0, 'pos':1, 'met':2, 'tot':3}
species = HUMAN
mm_refgenome
hg_refgenome = 'hg38'
    
count_matrix(bed_path = bed_path,
         output_path = output_path,
         col_dict = col_dict,
         species = species,
         windows_width=10000,
         mm_refgenome=mm_refgenome,
         hg_refgenome=hg_refgenome)

