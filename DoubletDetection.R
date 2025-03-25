
#learned from https://biostatsquid.com/doubletfinder-tutorial/ 
#and https://github.com/chris-mcginnis-ucsf/DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(dplyr)
library(DoubletFinder)
library(Seurat)
library(writexl)
library(readxl)
setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output")

#run for WT
wt1 <- readRDS("~/BINF/yushi scrnaseq/hippocampus/seurat output/WT/hippo_wt.rds")

## pK Identification (no ground-truth
sweep.res.list_wt <- paramSweep(wt1, PCs = 1:11, sct = FALSE)
sweep.stats_wt <- summarizeSweep(sweep.res.list_wt, GT = FALSE)
bcmvn_wt<- find.pK(sweep.stats_wt)

## pK Identification (ground-truth)
sweep.res.list_wt <- paramSweep(wt1, PCs = 1:11, sct = FALSE)
sweep_stats_wt <- summarizeSweep(sweep.res.list_wt)
bcmvn <- find.pK(sweep_stats_wt) # computes a metric to find the optimal pK value 
# Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
optimal.pk <- bcmvn %>% 
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))

# Homotypic doublet proportion estimate
annotations <- wt1@meta.data$singler.labels # use the clusters as the user-defined cell types
homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets





'''
gt.calls <- wt1@meta.data[rownames(sweep.res.list_wt[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_wt <- summarizeSweep(sweep.res.list_wt, GT = TRUE, GT.calls = gt.calls)
bcmvn_wt <- find.pK(sweep.stats_wt)
'''