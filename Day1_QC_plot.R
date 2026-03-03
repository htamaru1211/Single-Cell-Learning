library(Seurat)
data("pbmc_small")
VlnPlot(pbmc_small, features = c("nFeature_RNA", "nCount_RNA"))