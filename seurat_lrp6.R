#two standalone analyses: treated and WT
setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output/treated")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)
library(Seurat)

library(SingleR)
library(celldex)
#process treated

scdata =Seuratscdata = ReadMtx(mtx = "~/BINF/yushi scrnaseq/hippocampus/sb output/treated/treated_RSEC_MolsPerCell_MEX/matrix.mtx.gz",
                   cells = "~/BINF/yushi scrnaseq/hippocampus/sb output/treated/treated_RSEC_MolsPerCell_MEX/barcodes.tsv.gz",
                   features = "~/BINF/yushi scrnaseq/hippocampus/sb output/treated/treated_RSEC_MolsPerCell_MEX/features.tsv.gz")

hippo = CreateSeuratObject(counts = scdata, project = "hippocampus")

# Calculate percentage (percent.mt) of mitochondrial RNA and store in metadata 
# [[ ]] operator adds a column to metadata
hippo[["percent.mt"]] <- PercentageFeatureSet(hippo, pattern = "^mt-")
hippo[["percent.ribo"]] <- PercentageFeatureSet(hippo, pattern = "^Rp[ls]")
head(hippo@meta.data)


VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

#find number of cells before filtering
nrow(hippo@meta.data)
#989 cells

hippo = subset(hippo, subset = percent.mt < 25)
nrow(hippo@meta.data)
#932 cells if mt = 25%
#low quality. 




VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

hippo <- subset(hippo, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000)
nrow(hippo@meta.data)
#923

hippo <- NormalizeData(hippo)

hippo  <- FindVariableFeatures(hippo , selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(hippo), 10)
#sacling
all.genes = row.names(hippo)
hippo = ScaleData(hippo, features = all.genes)

#pca
hippo = RunPCA(hippo, features = VariableFeatures(object = hippo))
heat = DimHeatmap(hippo, dims = 1:20, cells = 500, balanced = TRUE)
elbow = ElbowPlot(hippo)
#clear elbow at 6

hippo = FindNeighbors(hippo, dims = 1:6)
resolution.range <- seq(from = 0, to = 1, by = 0.1)



# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  hippo<- FindClusters(hippo, resolution = res)
  
  # Find all markers for the clusters at this resolution
  hippo.markers <- FindAllMarkers(hippo, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(hippo.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  file_new = paste0("filt",file)
  write_xlsx(dff, file_new)
}

clustr = clustree(hippo)
#res = 0.6
hippo = RunUMAP(hippo, dims = 1:6)
DimPlot(hippo, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.5")

saveRDS(hippo,"hippo.rds")
#annotation using singleR

hippo = readRDS("hippo.rds")

ref = celldex::MouseRNAseqData()
View(as.data.frame(colData(ref)))

hippo.count = GetAssayData(hippo, layer = "counts")

pred = SingleR(test = hippo.count,
               ref = ref,
               labels = ref$label.main)

hippo$singler.labels = pred$labels[match(rownames(hippo@meta.data), rownames(pred))]
DimPlot(hippo, reduction = "umap", group.by = "singler.labels",label = TRUE, pt.size = 0.5,repel = TRUE)+ NoLegend()
saveRDS(hippo,"treated_anntd.rds")


#####################

#WT
setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output/WT")

scdata =Seuratscdata = ReadMtx(mtx = "~/BINF/yushi scrnaseq/hippocampus/sb output/WT/WT_clean_RSEC_MolsPerCell_MEX/matrix.mtx.gz",
                               cells = "~/BINF/yushi scrnaseq/hippocampus/sb output/WT/WT_clean_RSEC_MolsPerCell_MEX/barcodes.tsv.gz",
                               features = "~/BINF/yushi scrnaseq/hippocampus/sb output/WT/WT_clean_RSEC_MolsPerCell_MEX/features.tsv.gz")

hippo = CreateSeuratObject(counts = scdata, project = "hippocampus")

# Calculate percentage (percent.mt) of mitochondrial RNA and store in metadata 
# [[ ]] operator adds a column to metadata
hippo[["percent.mt"]] <- PercentageFeatureSet(hippo, pattern = "^mt-")
hippo[["percent.ribo"]] <- PercentageFeatureSet(hippo, pattern = "^Rp[ls]")
head(hippo@meta.data)


VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

#find number of cells before filtering
nrow(hippo@meta.data)
#1192 cells

hippo = subset(hippo, subset = percent.mt < 25)
nrow(hippo@meta.data)
#966 cells if mt = 25%
#low quality. 




VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

hippo <- subset(hippo, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000)
nrow(hippo@meta.data)
#938

hippo <- NormalizeData(hippo)

hippo  <- FindVariableFeatures(hippo , selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(hippo), 10)
#sacling
all.genes = row.names(hippo)
hippo = ScaleData(hippo, features = all.genes)

#pca
hippo = RunPCA(hippo, features = VariableFeatures(object = hippo))
heat = DimHeatmap(hippo, dims = 1:20, cells = 500, balanced = TRUE)
elbow = ElbowPlot(hippo)
#elbow at 11

hippo = FindNeighbors(hippo, dims = 1:11)
resolution.range <- seq(from = 0, to = 1, by = 0.1)



# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  hippo<- FindClusters(hippo, resolution = res)
  
  # Find all markers for the clusters at this resolution
  hippo.markers <- FindAllMarkers(hippo, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(hippo.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  file_new = paste0("filt",file)
  write_xlsx(dff, file_new)
}

clustr = clustree(hippo)
#res = 0.5
hippo = RunUMAP(hippo, dims = 1:11)
DimPlot(hippo, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.5")


ref = celldex::MouseRNAseqData()
View(as.data.frame(colData(ref)))

hippo.count = GetAssayData(hippo, layer = "counts")

pred = SingleR(test = hippo.count,
               ref = ref,
               labels = ref$label.main)

hippo$singler.labels = pred$labels[match(rownames(hippo@meta.data), rownames(pred))]
DimPlot(hippo, reduction = "umap", group.by = "singler.labels",label = TRUE, pt.size = 0.5,repel = TRUE)+ NoLegend()
saveRDS(hippo,"hippo_wt.rds")
