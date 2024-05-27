import scanpy as sc
import celltypist
from celltypist import models
import pickle

short_names = ['P03','P04','P05','P06','P07','P08','P09','P10','P12','P13','P14','P15','P19','P20','P22','P23','P26']

#Path
save_data = '../data/performance_test/'


#Celltypist model downolad:
models.download_models(
    force_update=True, model= "Immune_All_High.pkl")

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

data_to_save =[]
for i,sample in enumerate(short_names):
   
    adata = sc.read_h5ad(save_data+sample+".h5ad")
    

    anotated = annotation(adata)
    anotated.X = anotated.X.tocsr()

    #anotated.write(save_data + sample+".h5ad")
    #data_to_save[sample]=anotated # Correct way to save data
    data_to_save.append(anotated)

with open("GSE205013_python.pkl", "wb") as file:
    pickle.dump(data_to_save, file)