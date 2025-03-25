setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output/integration")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)
library(Seurat)
library(harmony)

library(SingleR)
library(celldex)

treated <- readRDS("~/BINF/yushi scrnaseq/hippocampus/seurat output/treated/treated_anntd.rds")
wt <- readRDS("~/BINF/yushi scrnaseq/hippocampus/seurat output/WT/hippo_wt.rds")
treated$condition = "treated"
wt$condition = "WT"

Combined = merge(treated, y = wt)
Combined = NormalizeData(Combined)
Combined = FindVariableFeatures(Combined)
Combined = ScaleData(Combined)
Combined = RunPCA(Combined)

Combined = RunHarmony(Combined, group.by.vars = "condition")
heat = DimHeatmap(Combined, dims = 1:20, cells = 500, balanced = TRUE)
elbow = ElbowPlot(Combined)
#elbow at 6

Combined = FindNeighbors(Combined, dims = 1:6)
resolution.range = seq(from=0, to = 1, by = 0.1)


# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  Combined<- FindClusters(Combined, resolution = res)
  
  # Find all markers for the clusters at this resolution
  Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(Combined.markers, file_name)
  
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

clustr = clustree(Combined)


Combined = RunUMAP(Combined, dims = 1:6)
DimPlot(Combined, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.4")

DimPlot(Combined, reduction = "umap", group.by = "condition")


ref = celldex::MouseRNAseqData()


Combined = JoinLayers(Combined)
Combined.count = LayerData(Combined, assay = "RNA", layer = "counts")


pred = SingleR(test = Combined.count,
               ref = ref,
               labels = ref$label.main)

Combined$singler.labels = pred$labels[match(rownames(Combined@meta.data), rownames(pred))]
DimPlot(Combined, reduction = "umap", group.by = "singler.labels",label = TRUE, pt.size = 0.5,repel = TRUE)+ NoLegend()
saveRDS(Combined,"integrated.rds")


