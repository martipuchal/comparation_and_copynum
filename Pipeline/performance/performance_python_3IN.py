import scanpy as sc
import scanorama
import pickle

short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']


#Path
save_data = '../data/performance_test/'

'''
sample_list=[]
for sample in short_names:
    adata = sc.read_h5ad(save_data+sample+".h5ad")
    adata.X = adata.X.tocsr()
    sample_list.append(adata)
    print(sample)
'''
with open("GSE205013_python.pkl", "rb") as file:
    sample_list = pickle.load(file)
print('start integration')

adatas_cor = scanorama.correct_scanpy(sample_list, return_dimred=True)


#Integration of the datasets:
adata_spatial = sc.concat(
    adatas_cor,
    label="library_id",
    uns_merge="unique",
    keys=short_names,
    index_unique="-",
)

adata_spatial.write(save_data + "allsamples.h5ad")