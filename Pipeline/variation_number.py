import scanpy as sc
import numpy as np
import pandas as pd

short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']
to_data = '../data/'
save_graph = '../fig/'

cnv_score_threshold = 0.9 #Top 10% are cancer cells

def discriminator (cnv_score,num):
    if cnv_score<num:
        return 'normal'
    else:
        return 'tumoral'

list_resume =[]
list_resume.append(('sample','total_epithelial','epithelial_and_tumoral','total_tumoral'))

sample_list = []

for sample_name in short_names:
    adata = sc.read_h5ad(to_data+sample_name + "_postinfercnv.h5ad")

    cnv_total = np.max(adata.obs['cnv_score'])
    threshold = cnv_total*cnv_score_threshold
    print(sample_name)
    print(threshold)
    adata.obs['diagnosis'] = adata.obs['cnv_score'].map(lambda x : discriminator(x,threshold) )
    
    total_cell=adata.n_obs
    total_cancer= (adata.obs['diagnosis'] == 'tumoral').sum()
    print(total_cancer)
    total_epithlial_cells = (adata.obs['celltypist_cell_label'] == 'Epithelial cells').sum()
    total_cancer_epithelial = ((adata.obs['diagnosis'] == 'tumoral') & (adata.obs['celltypist_cell_label'] == 'Epithelial cells')).sum()    
    
    list_resume.append((sample_name,total_epithlial_cells,total_cancer_epithelial,total_cancer))

    sample_list.append(adata)

    adata.write(to_data + sample_name + "_postinfercnv.h5ad")

    sc.pl.umap(
    adata, color='diagnosis', palette=sc.pl.palettes.default_20, size = 1
    , save =sample_name+ 'tumor_normal.png')

print('file')
try:
    with open(to_data + 'copynum_per_sample_resum.txt', 'w') as file:
        for item in list_resume:
            file.write("%s\n" % str(item))
    print("Resumen de localización guardado correctamente.")
except IOError as e:
    print("Error al escribir el archivo de resumen de localización:", e)


integred = sc.read_h5ad(to_data+"integrated_object.h5ad")

print(integred.obs_names)
integred.obs['diagnosis'] = 'normal'

for i, sample in enumerate(sample_list):
    barcodes = sample.obs_names[sample.obs['diagnosis'] == 'tumoral'] + '-' + short_names[i]
    integred.obs.loc[barcodes, 'diagnosis'] = 'tumoral'

integred.write(to_data +  "integrated_object_cn.h5ad")

sc.pl.umap(
    integred, color="celltypist_cell_label", palette=sc.pl.palettes.default_20, size = 1
    , save = 'final_Inmuno_High.png')

sc.pl.umap(
    integred, color='diagnosis', palette=sc.pl.palettes.default_20, size = 1
    , save ='final_normal_tumor.png')



