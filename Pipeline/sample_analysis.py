
#Libraries:
import scanpy as sc
import celltypist
from celltypist import models
import scanorama
import pandas as pd
import matplotlib.pyplot as plt

#Samples name list:
short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']
genes = ['KRT19','CD3E','VWF','CD68','MKI67','CD79A','EPCAM']

#Path:
to_data = '/raw/'
save_graph = '/fig/'


#Open files names:
file_names = open (to_data+'file_list.txt','r')
file_names_list = []
for name in file_names:
    file_names_list.append(name.replace('\n',''))

#Quality control parameters:
# Quality cells:
n_genes_threshold = 500
n_barcodes_threshold = 1500
pct_mito_threshold = 15
# Cells with >1% of transcripts representing erythroid genes (HBA1, HBA2, HBB, HBM, and ALAS2) were removed forrm the analissis.
e_genes = ['HBA1','HBA1','HBB','HBM','ALAS1']
e_genes_threshold = 0.01


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

#Sample analyssis loop:

sample_list = []
for i,sample in enumerate(short_names):
    adata = QC(i)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save =save_graph+ 'VIO_'+sample+'_raw.png')

    adata_QC = filtering(adata)

    sc.pl.violin(adata_QC, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True, save =save_graph+ 'VIO_' +sample+ '_qc.png')

    sample_list.append(adata_QC)


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

#Preprocessing sample loop:
for i,adata in enumerate(sample_list):
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)
    
    
    adata.X = adata.X.tocsr()
    sample_list[i] = annotation(adata)

#Scanorama bach correction:
adatas_cor = scanorama.correct_scanpy(sample_list, return_dimred=True)

#Integration of the datasets:
adata_spatial = sc.concat(
    adatas_cor,
    label="library_id",
    uns_merge="unique",
    keys=short_names,
    index_unique="-",
)

#Sample clustering:
sc.pp.neighbors(adata_spatial, use_rep="X_scanorama", random_state= 1)
sc.tl.umap(adata_spatial, random_state=1)
sc.tl.leiden(
    adata_spatial, key_added="clusters", directed=False,resolution = 0.1, random_state=1
)
sc.pp.scale(adata_spatial, max_value=10,zero_center=False)

#Umap generation:
sc.pl.umap(
    adata_spatial, color=["clusters", "library_id"], palette=sc.pl.palettes.default_20
    ,save= save_graph+'scanorma_umap.png')

sc.pl.umap(
    adata_spatial, color=["clusters", "celltypist_cell_label"], palette=sc.pl.palettes.default_20
    , save =save_graph+ 'umap_Inmuno_High.png')

for e in genes:
    sc.pl.umap(
    adata_spatial, color=["celltypist_cell_label", e], palette=sc.pl.palettes.default_20
        , save =save_graph+ e+'.png')
    
