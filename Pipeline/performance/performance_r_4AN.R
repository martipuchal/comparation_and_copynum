#Libraries:
library(Seurat)
library(SingleR)
library(celldex)
library(stringr)

wd <- file.path("/home", "data", "performance_test")
resultsDir <- wd
setwd(wd)
print('1')
dwIntegrated <- readRDS(file.path(wd, "dwIntegrated.processed.rds"))


Idents(dwIntegrated) <- "orig.ident"
Idents(dwIntegrated) <- "seurat.clusters"
DefaultAssay(dwIntegrated) <- "RNA"
print(dwIntegrated)


ref <- celldex::BlueprintEncodeData()

DefaultAssay(dwIntegrated) <- "integrated"
print('2')

te <- as.SingleCellExperiment(dwIntegrated)
print('te')
clu <- dwIntegrated$seurat_clusters
print('clu')
la <- ref$label.main
print('la')

predClust <- SingleR(test = as.SingleCellExperiment(dwIntegrated), ref = ref, clusters = dwIntegrated$seurat_clusters, labels = ref$label.main)
print(predClust)
datatable(as.data.frame(table(predClust$labels)), colnames = "", rownames = "")
print('3')
dwIntegrated$bpClust <- dwIntegrated$seurat_clusters
levels(dwIntegrated$bpClust) <- predClust$labels
print('4')

saveRDS(dwIntegrated, file.path(wd, "dwIntegrated.annotated.rds"))