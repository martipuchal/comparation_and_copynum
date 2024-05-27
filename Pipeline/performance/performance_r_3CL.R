#Libraries:
library(Seurat)
library(SingleR)
library(celldex)
library(stringr)
wd <- file.path("/home", "data", "performance_test")
resultsDir<- wd
setwd(wd) 


load(file = file.path(wd, "dwIntegrated.RData"))
dwIntegrated <- ScaleData(dwIntegrated)
dwIntegrated <- RunPCA(dwIntegrated, npcs = 30)

dwIntegrated <- FindNeighbors(dwIntegrated, reduction="pca", k.param=30) #minimum distance of 0.3, I do not see the param??
dwIntegrated <- FindClusters(dwIntegrated, resolution=0.1) #checked 0.3 and 0.6 but many clusters were created
dwIntegrated <- RunUMAP(dwIntegrated, dims = 1:16, n.neighbors = 30L)

saveRDS(dwIntegrated,file.path(wd,"dwIntegrated.processed.rds"))