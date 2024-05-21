import scanpy as sc
import infercnvpy # type: ignore
import celltypist
from  celltypist import models
import pandas as pd 
import numpy as np
import re

#Samples name list:
short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']

#Path:
to_data = '../raw/'
save_data = '../data/'
save_graph = '../fig/'

print("Open files...")
#Open files names:
file_names = open (to_data+'file_list.txt','r')
file_names_list = []
for name in file_names:
    file_names_list.append(name.replace('\n',''))
print("Files Opened")

#Quality control parameters:
# Quality cells:
n_genes_threshold = 500
n_barcodes_threshold = 1500
pct_mito_threshold = 15
# Cells with >1% of transcripts representing erythroid genes (HBA1, HBA2, HBB, HBM, and ALAS2) were removed forrm the analissis.
e_genes = ['HBA1','HBA1','HBB','HBM','ALAS1']
e_genes_threshold = 0.01

print("Functions definition...")
#Calcul QC parameters funcction:
def QC(selected):

    #Anndata oppening:
    adata = sc.read_10x_mtx(to_data , var_names='gene_symbols',cache=True,prefix=file_names_list[selected])


    #Mitocondrial genes:    
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    #Metrics calculation
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    return(adata)
#Filtering cells:
def filtering (adata):
    #Filter cells:
    sc.pp.filter_cells(adata, min_counts = n_barcodes_threshold, min_genes = None, inplace=True) #Select cells by molecular identiers
    sc.pp.filter_cells(adata, min_counts = None, min_genes = n_genes_threshold, inplace=True) #select cells by detectable gens


    #Cells selection:
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt<pct_mito_threshold,:]

    #Erythroid genes:
    adata.obs['erythroid_gene_counts'] = adata[:, e_genes].X.sum(axis=1)
    adata.obs['erythroid_gene_percentage'] = (
        adata.obs['erythroid_gene_counts'] / adata.obs['total_counts'] * 100)
    num_cells_excluded = sum(adata.obs['erythroid_gene_percentage'] > 1)
    #print(f"Number of cells excluded: {num_cells_excluded}")
    adata = adata[adata.obs['erythroid_gene_percentage'] < 1]
    adata.obs = adata.obs.drop(['erythroid_gene_counts', 'erythroid_gene_percentage'], axis=1)


    adata.layers['counts'] = adata.X

    return (adata)

#Annotation function:

#Celltypist model downolad:
models.download_models(
    force_update=True, model= "Immune_All_High.pkl")

#Annotation function:
def annotation(adata):
    #Annotation:
    adata_celltypist = adata.copy()  # make a copy of our adata
    adata_celltypist.X = adata.layers['counts'] # set adata.X to raw counts
    sc.pp.normalize_per_cell(
        adata_celltypist, counts_per_cell_after=10**4
    )  # normalize to 10,000 counts per cell
    sc.pp.log1p(adata_celltypist)  # log-transform
    # make .X dense instead of sparse, for compatibility with celltypist:
    adata_celltypist.X = adata_celltypist.X.toarray()
    
    model_high = models.Model.load(model="Immune_All_High.pkl")
    predictions_high = celltypist.annotate(
        adata_celltypist, model=model_high, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata()
    adata.obs["celltypist_cell_label"] = predictions_high_adata.obs.loc[
        adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score"] = predictions_high_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]

    return (adata)


def gtf(gtf):
    df = pd.read_csv(to_data+'hg38_gencode_v27.txt',sep='\t',names=['gene_name','chromosome','start','end'] )
    print(df)
    return (df)
print("Functions defined")
print("GTF generation...")
# DF for gene localization:
try:
    gencode_genes=pd.read_csv(to_data+"genes_dataframe.csv")
except IOError:
    gencode_genes=gtf('genes.gtf')
gencode_genes.set_index('gene_name', inplace=True)
print("GTF generated")
print(gencode_genes)

#Sample analyssis loop:
gene_count_per_sample = []
sample_list = []
for i,sample in enumerate(short_names):
    print("sample: ",short_names[i])
    adata = QC(i)

    adata_QC = filtering(adata)

    sc.pp.normalize_total(adata_QC, inplace=True)
    sc.pp.log1p(adata_QC)
    sc.pp.highly_variable_genes(adata_QC, flavor="seurat", n_top_genes=2000, inplace=True)
    
    adata_PP = annotation(adata_QC)

    #Gene location
    chr_list = []
    start_list = []
    stop_list  = []

    count = 0
    gene_count = 0
    print("localitzation")
    adata_gene_names = adata_PP.var_names.tolist()
    for gene in adata_gene_names:
        count += 1
        if gene in gencode_genes.index.values:
            print('Gene: ',gene)
            line = gencode_genes.loc[gene]
            print('Line: ',line)
            chr_list.append( line['chromosome'])
            start_list.append( line['start'])
            stop_list.append( line['end']  )  
            gene_count += 1  
        else:
            adata_PP = adata_PP[:, adata_PP.var_names != gene]

    adata_PP.var['chromosome'] = chr_list
    adata_PP.var['start'] = start_list
    adata_PP.var['end'] = stop_list

    final_count = str(short_names[i])+': '+str(gene_count)+'/'+str(count)
    gene_count_per_sample.append(final_count)
    print("Filtering done")
    print(final_count)
    sc.tl.pca(adata_PP,random_state=1)
    sc.pp.neighbors(adata_PP,random_state=1)
    sc.tl.umap(adata_PP,random_state=1)
    sc.tl.leiden(adata_PP, key_added="clusters", directed=False,resolution = 0.1, random_state=1)
    

    print("infercnv")
    infercnvpy.tl.infercnv(adata_PP, reference_key='celltypist_cell_label', reference_cat='T cells', reference=None
                           , lfc_clip=3, window_size=100, step=10, dynamic_threshold=1.5
                           , exclude_chromosomes=('chrX', 'chrY'), chunksize=5000, n_jobs=None
                           , inplace=True, layer=None, key_added='cnv')
    
    infercnvpy.tl.pca(adata_PP,random_state=1)
    infercnvpy.pp.neighbors(adata_PP,random_state=1)
    infercnvpy.tl.leiden(adata_PP,random_state=1)
    infercnvpy.tl.umap(adata_PP,random_state=1)
    infercnvpy.tl.cnv_score(adata_PP)
    print("save")
    adata_PP.write(save_data + sample + "_postinfercnv.h5ad")
    
    
    sample_list.append(adata_QC)
    infercnvpy.pl.umap(
    adata_PP, color=["cnv_leiden", "celltypist_cell_label"], palette=sc.pl.palettes.default_20
    , save = 'infercnv_Inmuno_High_'+sample+'_infercnv.png')
    infercnvpy.pl.umap(
    adata_PP, color=[ "celltypist_cell_label","cnv_score"], palette=sc.pl.palettes.default_20
    , save ='infercnv_tumoral_'+sample+'_infercnv.png')
    sc.pl.umap(
    adata_PP, color=["cnv_leiden", "celltypist_cell_label"], palette=sc.pl.palettes.default_20
    , save = 'scanpy_Inmuno_High_'+sample+'_infercnv.png')
    sc.pl.umap(
    adata_PP, color=[ "celltypist_cell_label","cnv_score"], palette=sc.pl.palettes.default_20
    , save ='scanpy_tumoral_'+sample+'_infercnv.png')
try:
    with open(save_data + 'localization_resume.txt', 'w') as file:
        for item in gene_count_per_sample:
            file.write("%s\n" % item)
    print("Resumen de localización guardado correctamente.")
except IOError as e:
    print("Error al escribir el archivo de resumen de localización:", e)
