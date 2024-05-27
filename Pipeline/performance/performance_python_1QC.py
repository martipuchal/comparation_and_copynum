import scanpy as sc

#Path
to_data = '../raw/'
save_data = '../data/performance_test/'

short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']

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


for i,sample in enumerate(short_names):
    adata = QC(i)

    adata = filtering(adata)

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)
    
    adata.write(save_data + sample+".h5ad")
