#Libraries:
library(Seurat)
library(SingleR)
library(celldex)
library(stringr)

#Directories:
# path to raw: ../raw/
# path to working dir.: ../data/performance_test
dataDir <- file.path("/home","raw")
wd <- file.path("/home", "data", "performance_test")
resultsDir<- wd
setwd(wd) 

print(getwd())

#Open data:
P03 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204111_P03_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204111_P03_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204111_P03_barcodes.tsv.gz")) 
P04 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204112_P04_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204112_P04_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204112_P04_barcodes.tsv.gz")) 
P05 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204113_P05_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204113_P05_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204113_P05_barcodes.tsv.gz")) 
P06 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204114_P06_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204114_P06_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204114_P06_barcodes.tsv.gz")) 
P07 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204115_P07_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204115_P07_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204115_P07_barcodes.tsv.gz")) 
P08 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204116_P08_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204116_P08_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204116_P08_barcodes.tsv.gz")) 
P09 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204117_P09_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204117_P09_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204117_P09_barcodes.tsv.gz")) 
P10 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204118_P10_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204118_P10_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204118_P10_barcodes.tsv.gz")) 
P12 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204120_P12_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204120_P12_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204120_P12_barcodes.tsv.gz")) 
P13 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204121_P13_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204121_P13_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204121_P13_barcodes.tsv.gz")) 
P14 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204122_P14_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204122_P14_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204122_P14_barcodes.tsv.gz")) 
P15 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204123_P15_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204123_P15_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204123_P15_barcodes.tsv.gz")) 
P19 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204127_P19_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204127_P19_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204127_P19_barcodes.tsv.gz")) 
P20 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204128_P20_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204128_P20_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204128_P20_barcodes.tsv.gz")) 
P22 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204130_P22_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204130_P22_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204130_P22_barcodes.tsv.gz")) 
P23 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204131_P23_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204131_P23_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204131_P23_barcodes.tsv.gz")) 
P26 <- ReadMtx(mtx=paste0(dataDir,"/GSM6204134_P26_matrix.mtx.gz"), features=paste0(dataDir,"/GSM6204134_P26_features.tsv.gz"), cells=paste0(dataDir,"/GSM6204134_P26_barcodes.tsv.gz")) 

print("Build sample objects:")
#Build sample objects:
P03 <- CreateSeuratObject(counts=P03, project = 'P03', min.cells=3, min.features=200)
P03[["percent.mt"]] <- PercentageFeatureSet(P03, pattern = "^MT-")
P04 <- CreateSeuratObject(counts=P04, project = 'P04', min.cells=3, min.features=200)
P04[["percent.mt"]] <- PercentageFeatureSet(P04, pattern = "^MT-")
P05 <- CreateSeuratObject(counts=P05, project = 'P05', min.cells=3, min.features=200)
P05[["percent.mt"]] <- PercentageFeatureSet(P05, pattern = "^MT-")
P06 <- CreateSeuratObject(counts=P06, project = 'P06', min.cells=3, min.features=200)
P06[["percent.mt"]] <- PercentageFeatureSet(P06, pattern = "^MT-")
P07 <- CreateSeuratObject(counts=P07, project = 'P07', min.cells=3, min.features=200)
P07[["percent.mt"]] <- PercentageFeatureSet(P07, pattern = "^MT-")
P08 <- CreateSeuratObject(counts=P08, project = 'P08', min.cells=3, min.features=200)
P08[["percent.mt"]] <- PercentageFeatureSet(P08, pattern = "^MT-")
P09 <- CreateSeuratObject(counts=P09, project = 'P09', min.cells=3, min.features=200)
P09[["percent.mt"]] <- PercentageFeatureSet(P09, pattern = "^MT-")
P10 <- CreateSeuratObject(counts=P10, project = 'P10', min.cells=3, min.features=200)
P10[["percent.mt"]] <- PercentageFeatureSet(P10, pattern = "^MT-")
P12 <- CreateSeuratObject(counts=P12, project = 'P12', min.cells=3, min.features=200)
P12[["percent.mt"]] <- PercentageFeatureSet(P12, pattern = "^MT-")
P13 <- CreateSeuratObject(counts=P13, project = 'P13', min.cells=3, min.features=200)
P13[["percent.mt"]] <- PercentageFeatureSet(P13, pattern = "^MT-")
P14 <- CreateSeuratObject(counts=P14, project = 'P14', min.cells=3, min.features=200)
P14[["percent.mt"]] <- PercentageFeatureSet(P14, pattern = "^MT-")
P15 <- CreateSeuratObject(counts=P15, project = 'P15', min.cells=3, min.features=200)
P15[["percent.mt"]] <- PercentageFeatureSet(P15, pattern = "^MT-")
P19 <- CreateSeuratObject(counts=P19, project = 'P19', min.cells=3, min.features=200)
P19[["percent.mt"]] <- PercentageFeatureSet(P19, pattern = "^MT-")
P20 <- CreateSeuratObject(counts=P20, project = 'P20', min.cells=3, min.features=200)
P20[["percent.mt"]] <- PercentageFeatureSet(P20, pattern = "^MT-")
P22 <- CreateSeuratObject(counts=P22, project = 'P22', min.cells=3, min.features=200)
P22[["percent.mt"]] <- PercentageFeatureSet(P22, pattern = "^MT-")
P23 <- CreateSeuratObject(counts=P23, project = 'P23', min.cells=3, min.features=200)
P23[["percent.mt"]] <- PercentageFeatureSet(P23, pattern = "^MT-")
P26 <- CreateSeuratObject(counts=P26, project = 'P26', min.cells=3, min.features=200)
P26[["percent.mt"]] <- PercentageFeatureSet(P26, pattern = "^MT-")
#save objects generated
save(P03, P04, P05, P06, P07, P08, P09, P10, P12, P13, P14, P15,
     P19, P20, P22, P23, P26, file="GSE205013.raw.seurat.RData")
print("Open and save raw count matrix:")
#Open and save raw count matrix:



exp.rawdataP03 <- as.matrix(P03@assays$RNA$counts)
exp.rawdataP04 <- as.matrix(P04@assays$RNA$counts)
exp.rawdataP05 <- as.matrix(P05@assays$RNA$counts)
exp.rawdataP06 <- as.matrix(P06@assays$RNA$counts)
exp.rawdataP07 <- as.matrix(P07@assays$RNA$counts)
exp.rawdataP08 <- as.matrix(P08@assays$RNA$counts)
exp.rawdataP09 <- as.matrix(P09@assays$RNA$counts)
exp.rawdataP10 <- as.matrix(P10@assays$RNA$counts)
exp.rawdataP12 <- as.matrix(P12@assays$RNA$counts)
exp.rawdataP13 <- as.matrix(P13@assays$RNA$counts)
exp.rawdataP14 <- as.matrix(P14@assays$RNA$counts)
exp.rawdataP15 <- as.matrix(P15@assays$RNA$counts)
exp.rawdataP19 <- as.matrix(P19@assays$RNA$counts)
exp.rawdataP20 <- as.matrix(P20@assays$RNA$counts)
exp.rawdataP22 <- as.matrix(P22@assays$RNA$counts)
exp.rawdataP23 <- as.matrix(P23@assays$RNA$counts)
exp.rawdataP26 <- as.matrix(P26@assays$RNA$counts)

saveRDS(exp.rawdataP03, file.path(wd,"exp.rawdataP03.rds"))
saveRDS(exp.rawdataP04, file.path(wd,"exp.rawdataP04.rds"))
saveRDS(exp.rawdataP05, file.path(wd,"exp.rawdataP05.rds"))
saveRDS(exp.rawdataP06, file.path(wd,"exp.rawdataP06.rds"))
saveRDS(exp.rawdataP07, file.path(wd,"exp.rawdataP07.rds"))
saveRDS(exp.rawdataP08, file.path(wd,"exp.rawdataP08.rds"))
saveRDS(exp.rawdataP09, file.path(wd,"exp.rawdataP09.rds"))
saveRDS(exp.rawdataP10, file.path(wd,"exp.rawdataP10.rds"))
saveRDS(exp.rawdataP12, file.path(wd,"exp.rawdataP12.rds"))
saveRDS(exp.rawdataP13, file.path(wd,"exp.rawdataP13.rds"))
saveRDS(exp.rawdataP14, file.path(wd,"exp.rawdataP14.rds"))
saveRDS(exp.rawdataP15, file.path(wd,"exp.rawdataP15.rds"))
saveRDS(exp.rawdataP19, file.path(wd,"exp.rawdataP19.rds"))
saveRDS(exp.rawdataP20, file.path(wd,"exp.rawdataP20.rds"))
saveRDS(exp.rawdataP22, file.path(wd,"exp.rawdataP22.rds"))
saveRDS(exp.rawdataP23, file.path(wd,"exp.rawdataP23.rds"))
saveRDS(exp.rawdataP26, file.path(wd,"exp.rawdataP26.rds"))

print("Generate the dataframe:")
#Generate the dataframe:
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

samples <- c(P03, P04, P05, P06, P07, P08, P09, P10, P12, P13, P14, P15, P19, P20, P22, P23, P26)

samples_chr <- c("P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P12", "P13", "P14", "P15", "P19", "P20", "P22", "P23", "P26")

prefilt_cells <- c()
i=1
for (samp in samples) {
  samp_chr <- samples_chr[i]
  # png(filename=file.path(resultsDir,paste0("QC/Violin.QC.",levels(samp$orig.ident),".png")),width=1000,height=800)
  # print(VlnPlot(samp, features = feats, pt.size = 0.05))
  # dev.off()
  prefilt_cells[[samp_chr]] <- ncol(samp)
  i=i+1
}

cells_df <- data.frame(prefilt_cells)
print('filtering')
#Filtering:
P03 <- subset(P03, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
# VlnPlot(P03, features = feats, pt.size = 0.1) 
P04 <- subset(P04, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P05 <- subset(P05, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P06 <- subset(P06, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P07 <- subset(P07, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P08 <- subset(P08, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P09 <- subset(P09, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P10 <- subset(P10, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P12 <- subset(P12, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P13 <- subset(P13, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P14 <- subset(P14, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P15 <- subset(P15, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P19 <- subset(P19, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P20 <- subset(P20, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P22 <- subset(P22, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P23 <- subset(P23, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)
P26 <- subset(P26, nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA>1500 & percent.mt < 15)


samples_post <- c(P03, P04, P05, P06, P07, P08, P09, P10, P12, P13, P14, P15, P19, P20, P22, P23, P26)
postfilt_cells <- c()
for (samp in samples) {
  #print(ncol(samp))
  postfilt_cells <- c(postfilt_cells, ncol(samp))
}

cells_df <- rbind(cells_df,postfilt_cells)


print('normal')
P03 <- NormalizeData(P03, normalization.method = "LogNormalize", scale.factor = 10000)
P04 <- NormalizeData(P04, normalization.method = "LogNormalize", scale.factor = 10000)
P05 <- NormalizeData(P05, normalization.method = "LogNormalize", scale.factor = 10000)
P06 <- NormalizeData(P06, normalization.method = "LogNormalize", scale.factor = 10000)
P07 <- NormalizeData(P07, normalization.method = "LogNormalize", scale.factor = 10000)
P08 <- NormalizeData(P08, normalization.method = "LogNormalize", scale.factor = 10000)
P09 <- NormalizeData(P09, normalization.method = "LogNormalize", scale.factor = 10000)
P10 <- NormalizeData(P10, normalization.method = "LogNormalize", scale.factor = 10000)
P12 <- NormalizeData(P12, normalization.method = "LogNormalize", scale.factor = 10000)
P13 <- NormalizeData(P13, normalization.method = "LogNormalize", scale.factor = 10000)
P14 <- NormalizeData(P14, normalization.method = "LogNormalize", scale.factor = 10000)
P15 <- NormalizeData(P15, normalization.method = "LogNormalize", scale.factor = 10000)
P19 <- NormalizeData(P19, normalization.method = "LogNormalize", scale.factor = 10000)
P20 <- NormalizeData(P20, normalization.method = "LogNormalize", scale.factor = 10000)
P22 <- NormalizeData(P22, normalization.method = "LogNormalize", scale.factor = 10000)
P23 <- NormalizeData(P23, normalization.method = "LogNormalize", scale.factor = 10000)
P26 <- NormalizeData(P26, normalization.method = "LogNormalize", scale.factor = 10000)

save(P03, P04, P05, P06, P07, P08, P09, P10, P12, P13, P14, P15,
     P19, P20, P22, P23, P26, file="GSE205013.norm.seurat.RData")














