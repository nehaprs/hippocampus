setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output/integration/DGE")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)
library(harmony)

library(SingleR)
library(celldex)

#read Combined.rds

FeaturePlot(Combined, features = "Lrp6")

#Identify DGE MAST
Combined$celltype.condition = paste(Combined$singler.labels, Combined$condition, sep = "_")
table(Combined$celltype.condition)
Idents(Combined) = "celltype.condition"

#Loop for DGE

celltypes = unique(Combined$singler.labels)

#find cell types which have at least 2 treated cells and 2 WT cells
meta = Combined@meta.data
cell_counts = table(meta$singler.labels, meta$condition)

celltypes2 = rownames(cell_counts)[cell_counts[,"treated"] >= 2 & cell_counts[,"WT"] >= 2]


dge_results = list()


for(ct in celltypes){
  
  # Define identifiers for treated and WT
  ident_treated = paste(ct, "treated", sep = "_")
  ident_WT = paste(ct, "WT", sep = "_")
  
  
  # Run MAST
    dge = FindMarkers(object = Combined, ident.1 =ident_treated , ident.2 = ident_WT,
                    test.use = "MAST",min.cells.group = 0, logfc.threshold = 0.25, min.pct = 0.1)
    dge$GeneName = row.names(dge)
    dge_results[[ct]] = dge
    write_xlsx(dge, paste0("DGE_", ct,".xlsx"))
  
}


#DGE doe Lrp6


celltypes <- unique(Combined$singler.labels)
lrp6_results <- data.frame()

# Loop through each cell type
for (ct in celltypes) {
  
  # Subset cells of current cell type
  subset_cells <- subset(Combined, subset = singler.labels == ct)
  
  # Check counts per condition
  cond_counts <- table(subset_cells$condition)

  # Initialize placeholders
  res <- data.frame(celltype = ct, 
                    Lrp6_present = NA, 
                    p_val = NA, 
                    avg_log2FC = NA, 
                    pct_treated = NA, 
                    pct_WT = NA)

  # Check if Lrp6 is present in this subset
  if ("Lrp6" %in% rownames(subset_cells)) {
    res$Lrp6_present <- TRUE
    
    if (all(cond_counts[c("treated", "WT")] >= 2)) {
      Idents(subset_cells) <- subset_cells$condition
      
      # Perform differential expression for Lrp6
      de_res <- FindMarkers(subset_cells,
                            ident.1 = "treated",
                            ident.2 = "WT",
                            features = "Lrp6",
                            test.use = "MAST",
                            logfc.threshold = 0,
                            min.pct = 0)
      
      # Store relevant statistics
      res$p_val <- de_res$p_val
      res$avg_log2FC <- de_res$avg_log2FC
      res$pct_treated <- de_res$pct.1
      res$pct_WT <- de_res$pct.2
      
    } else {
      cat("Skipped DE for cell type:", ct, "- insufficient cells per condition.\n")
    }
    
  } else {
    res$Lrp6_present <- FALSE
    cat("Lrp6 not detected in cell type:", ct, "\n")
  }
  
  # Append results
  lrp6_results <- rbind(lrp6_results, res)
}

# Display final results clearly
print(lrp6_results)
