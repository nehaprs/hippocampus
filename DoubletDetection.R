
#learned from https://biostatsquid.com/doubletfinder-tutorial/ 
#and https://github.com/chris-mcginnis-ucsf/DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(dplyr)
library(tidyverse)
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

#assuming multiplet rate 2.5% based on estimation from https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/instruments/bd-rhapsody-single-cell-analysis-system/BD-Rhapsody-HT-System-brochure.pdf
nExp_poi <- round(0.025*nrow(wt1@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run DoubletFinder
wt1 <- doubletFinder(seu = wt1, 
                        PCs = 1:11, 
                        pK = optimal.pk, 
                        nExp = nExp_poi.adj)

# change name of metadata column with Singlet/Doublet information
colnames(wt1@meta.data)[grepl('DF.classifications.*', colnames(wt1@meta.data))] <- "doublet_finder"
double_finder_res <- wt1@meta.data['doublet_finder']
double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly




VlnPlot(wt1,  split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')


# Get doublets per sample
doublets_summary <- wt1@meta.data %>% 
  group_by( doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  #group_by(condition) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)

write_xlsx(doublets_summary,"doublets_summaryWT.xlsx")
saveRDS(wt1,"wt_doublets.rds")

############# treated

setwd("~/BINF/yushi scrnaseq/hippocampus/seurat output/treated")

#run rds of  treated. 
treated1 = treated

## pK Identification (no ground-truth
sweep.res.list_treated <- paramSweep(treated1, PCs = 1:6, sct = FALSE)
sweep.stats_treated <- summarizeSweep(sweep.res.list_treated, GT = FALSE)
bcmvn_treated<- find.pK(sweep.stats_treated)

## pK Identification (ground-truth)
sweep.res.list_treated <- paramSweep(treated1, PCs = 1:6, sct = FALSE)
sweep_stats_treated <- summarizeSweep(sweep.res.list_treated)
bcmvn <- find.pK(sweep_stats_treated) # computes a metric to find the optimal pK value 
# Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
optimal.pk <- bcmvn %>% 
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))

# Homotypic doublet proportion estimate
annotations <- treated1@meta.data$singler.labels # use the clusters as the user-defined cell types
homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets

#assuming multiplet rate 2.5% based on estimation from https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/instruments/bd-rhapsody-single-cell-analysis-system/BD-Rhapsody-HT-System-brochure.pdf
nExp_poi <- round(0.025*nrow(treated1@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run DoubletFinder
treated1 <- doubletFinder(seu = treated1, 
                     PCs = 1:6, 
                     pK = optimal.pk, 
                     nExp = nExp_poi.adj)

# change name of metadata column with Singlet/Doublet information
colnames(treated1@meta.data)[grepl('DF.classifications.*', colnames(treated1@meta.data))] <- "doublet_finder"
double_finder_res <- treated1@meta.data['doublet_finder']
double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly




VlnPlot(treated1,  split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')


# Get doublets per sample
doublets_summary <- treated1@meta.data %>% 
  group_by( doublet_finder) %>% 
  summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
  #group_by(condition) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(doublet_finder, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
  dplyr::select(-countT)

write_xlsx(doublets_summary,"doublets_summarytreated.xlsx")

