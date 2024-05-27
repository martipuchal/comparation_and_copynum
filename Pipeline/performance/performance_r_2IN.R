#Libraries:
library(Seurat)
library(SingleR)
library(celldex)
library(stringr)

wd <- file.path("/home", "data", "performance_test")
resultsDir <- wd
setwd(wd) 

load(file="GSE205013.norm.seurat.RData")


dwAnchors <- FindIntegrationAnchors(object.list=c(P03,P04,P05,P06,P07,P08,P09,P10,P12,P13,P14,P15,P19,P20,P22,P23,P26)) #default= dims = 1:30
dwIntegrated <- IntegrateData(anchorset=dwAnchors)

save(dwIntegrated, file=file.path(wd,"dwIntegrated.RData"))
