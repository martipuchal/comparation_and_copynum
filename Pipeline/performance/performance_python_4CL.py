import scanpy as sc


short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']


#Path
save_data = '../data/performance_test/'


adata_spatial = sc.read_h5ad(save_data + "allsamples.h5ad")

#Sample clustering:
sc.pp.neighbors(adata_spatial, use_rep="X_scanorama", random_state= 1)
sc.tl.umap(adata_spatial, random_state=1)
sc.tl.leiden(
    adata_spatial, key_added="clusters", directed=False,resolution = 0.1, random_state=1
)
sc.pp.scale(adata_spatial, max_value=10,zero_center=False)
adata_spatial.write(save_data + "allsamples.h5ad")