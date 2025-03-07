setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output")
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)

scdata = ReadMtx(mtx = "~/BINF/yushi scrnaseq/hippocampus/sb output/LRP6_clean_RSEC_MolsPerCell_MEX/matrix.mtx.gz",
                   cells = "~/BINF/yushi scrnaseq/hippocampus/sb output/LRP6_clean_RSEC_MolsPerCell_MEX/barcodes.tsv.gz",
                   features = "~/BINF/yushi scrnaseq/hippocampus/sb output/LRP6_clean_RSEC_MolsPerCell_MEX/features.tsv.gz")

hippo = CreateSeuratObject(counts = scdata, project = "hippocampus")

# Calculate percentage (percent.mt) of mitochondrial RNA and store in metadata 
# [[ ]] operator adds a column to metadata
hippo[["percent.mt"]] <- PercentageFeatureSet(hippo, pattern = "^mt-")
hippo[["percent.ribo"]] <- PercentageFeatureSet(hippo, pattern = "^Rp[ls]")
head(hippo@meta.data)


VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

#find number of cells before filtering
nrow(hippo@meta.data)
#2167 cells

hippo = subset(hippo, subset = percent.mt < 5)
nrow(hippo@meta.data)
#167 cells
#low quality. 