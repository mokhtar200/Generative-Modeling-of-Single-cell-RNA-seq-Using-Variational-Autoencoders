# scripts/01_preprocessing.R
library(Seurat)
library(ggplot2)

# Load data (PBMC 3k from 10x)
data_dir <- "data/raw/"
pbmc <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = pbmc, project = "PBMC3K", min.cells = 3, min.features = 200)

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5)

# Normalization
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling and PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Save processed object
saveRDS(seurat_obj, file = "data/processed/seurat_pbmc.rds")
