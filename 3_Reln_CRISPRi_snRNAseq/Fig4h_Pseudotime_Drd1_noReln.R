#Calculating pseudotime for Drd1/Oligodendrocyte subcluster
set.seed(1234)
#Load libraries
library(Seurat)
library(SeuratData)
library(monocle3)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)
library(tidyverse)

remotes::install_version("Seurat", version = "4.3.0.1")
#install packages
# devtools::install_github('satijalab/seurat-data')
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'))
# devtools::install_github('cole-trapnell-lab/monocle3')

SC0016 <- allRats
#Load Reln KD rn7 object
SC0016 <- readRDS("/Volumes/Day-Lab$/10X_SingleCell/2023_SC_0016/analysis/Seurat/R objects/allRats_souped_noDoub.rds")
SC0016$cellType <- Idents(SC0016)
DefaultAssay(SC0016) <- "RNA"
#Validate UMAP
DimPlot(object = SC0016,reduction = "umap")

#retrieve the counts from your Seurat object with:
counts <- as.data.frame((as.data.frame(GetAssayData(object = SC0016,
                                                    slot = "data",
                                                    assay = "RNA"))))

#Remove rows with Reln in them

counts <- counts[-which(rownames(counts) == "Reln"), ]

#And then subset your data based upon the genes (row names)

SC0016_sub <- subset(SC0016, features = rownames(counts))

#Subset based on gRNA
Idents(SC0016_sub)=SC0016_sub$target

lacZ <- subset(SC0016_sub, idents = "lacz")
Idents(lacZ)=lacZ$cellType

Reln <- subset(SC0016_sub, idents = "reln")
Idents(Reln)=Reln$cellType

Idents(SC0016_sub)=SC0016_sub$cellType

#Subcluster Drd1/oligos
#Subset for only Drd1/Oligos
D1 <- subset(SC0016_sub, subset = cellType == c("Drd1-MSN"))
#D1 <- D1_old

##D1
#Run standard dimensionality reduction and clustering workflow 
D1 <- ScaleData(D1,verbose = FALSE)
D1 <- FindVariableFeatures(D1)
D1 <- RunPCA(D1,npcs = 17,verbose = FALSE)
D1 <- RunUMAP(D1, reduction = "pca", dims = 1:17)
D1 <- FindNeighbors(D1, reduction = "pca", dims = 1:17)
D1 <- FindClusters(D1, resolution = 0.1)

#Visualize D1s on UMAP
DimPlot(D1,reduction = 'umap',label = TRUE) + NoLegend()

#Find genes of interest 
FeaturePlot(D1, features = c("Ebf1", "Htr4"))


#To run pseudotime, the Seurat4 object must be converted into a Monocle3 cds object
#To make cds object manually, we need to extract annotations, counts metadata, and barcodes
#SeuratWrappers contains documentation on a shorter way to do this, but it does not give as biologically believable results
#For example, when using SeuratWrappers, some Htr4 positive cells are assigned high pseudotimes whereas when creating the Monocle3 object manually, this does not occur

#Extract gene annotations
gene_annotation_D1 <- as.data.frame(rownames(D1@reductions[["pca"]]@feature.loadings),
                                         row.names = rownames(D1@reductions[["pca"]]@feature.loadings))
#Rename column
colnames(gene_annotation_D1) <- "gene_short_name"

#Extract barcodes
cell_metadata_D1 <- as.data.frame(D1@assays[["RNA"]]@counts@Dimnames[[2]],
                                       row.names = D1@assays[["RNA"]]@counts@Dimnames[[2]])
#Rename column
colnames(cell_metadata_D1) <- "barcode"

#Extract counts data
New_matrix_D1 <- D1@assays[["RNA"]]@counts
New_matrix_D1 <- New_matrix_D1[rownames(D1@reductions[["pca"]]@feature.loadings), ]
expression_matrix_D1 <- New_matrix_D1

#Create a cds object from Seurat using annotations, counts, and barcodes
cds_from_seurat_D1 <- new_cell_data_set(expression_matrix_D1,
                                             cell_metadata = cell_metadata_D1,
                                             gene_metadata = gene_annotation_D1)

#Monocle3 requires that partitions are calculated in conjunction with the clustering information
#However, we will not be using the partition data in this pseudotime analysis, so we are creating a filler partition element in the cds object to satisfy this requirement

#Create a dataframe that has the same number of rows as cds object - all values are "1" to stand in as partition data
recreate.partition_D1 <- c(rep(1, length(cds_from_seurat_D1@colData@rownames)))
#Add barcodes as names 
names(recreate.partition_D1) <- cds_from_seurat_D1@colData@rownames
recreate.partition_D1 <- as.factor(recreate.partition_D1)
#Add dummy dataframe to cds object 
cds_from_seurat_D1@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition_D1

#Add clustering data and UMAP dimensions from Seurat object to cds object so the UMAPs match
list_cluster_D1 <- D1@active.ident
names(list_cluster_D1) <- D1@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat_D1@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster_D1
cds_from_seurat_D1@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat_D1@int_colData@listData$reducedDims@listData[["UMAP"]] <-D1@reductions[["umap"]]@cell.embeddings
cds_from_seurat_D1@reduce_dim_aux$gene_loadings <- D1@reductions[["pca"]]@feature.loadings

#Perform pseudotime analysis
cds_from_seurat_D1 <- learn_graph(cds_from_seurat_D1, use_partition = F)

#Choose root cell/origin
cds_from_seurat_D1 <- order_cells(cds_from_seurat_D1, reduction_method = 'UMAP')

#Visualize on UMAP using Monocle3, should be the same as Seurat UMAP
plot_cells(cds_from_seurat_D1, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4)

#Visualize pseudotime on UMAP
plot_cells(cds_from_seurat_D1,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE,
)

#Write pseudotime data out as .csv
pseudo_D1 <- pseudotime(cds_from_seurat_D1)
write.csv(pseudo_D1, file = '/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/D1_OLD_rn7_pseudotime.csv', quote = FALSE)

#Save D1 object
saveRDS(D1, file = '/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/D1_noReln_rn7.rds')

D1_noReln <- readRDS('/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/D1_noReln_rn7.rds')
pseudo_csv <- read.csv('/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/D1_OLD_rn7_pseudotime.csv')
D1@meta.data$pseudo <- pseudo_csv$x
write.csv(D1@meta.data, file = '/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/D1_OLD_PseudoMeta.csv', quote = FALSE)


FeaturePlot(D1, "mCherry", order = T, split.by = "target", keep.scale = "all",
            min.cutoff = 0, col = c("lightgrey","#ed2324"))

FeaturePlot(D1, "pseudo", split.by = "target", order = T, keep.scale = "all",
            min.cutoff = 0, col = c("lightgrey","#ed2324"))

tmat <- t(expression_matrix_D1)
tmat <- as.data.frame(tmat)
D1_meta <- D1_noReln@meta.data

lacZ_meta <- subset(D1@meta.data, subset = target =="lacz")
Reln_meta <- subset(D1@meta.data, subset = target =="reln")
lacZ_tmat <- subset(tmat, subset = rownames(tmat) %in% rownames(lacZ_meta))
Reln_tmat <- subset(tmat, subset = rownames(tmat) %in% rownames(Reln_meta))

cor(tmat$mCherry, D1_meta$pseudo)
cor(lacZ_meta$pseudo, lacZ_tmat$mCherry)
cor(Reln_meta$pseudo, Reln_tmat$mCherry)

lacZ_all <- cbind(lacZ_meta, lacZ_tmat)
Reln_all <- cbind(Reln_meta, Reln_tmat)
write.csv(lacZ_all, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/lacZ_all.csv" )
write.csv(Reln_all, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 4 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/Reln_all.csv" )


VlnPlot(D1, "pseudo", split.by = "target")

lacZ_meta <- subset(D1@meta.data, subset = target =="lacz")
Reln_meta <- subset(D1@meta.data, subset = target =="reln")

write.csv(lacZ_meta, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 6 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/lacZ_meta_noReln.csv")
write.csv(Reln_meta, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 6 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/Reln_meta_noReln.csv")


#lacZ
lacZ_df <- as.data.frame(t(as.data.frame(GetAssayData(object = lacZ,
                                                    slot = "data",
                                                    assay = "RNA"))))
lacZ_df <- cbind(cell = rownames(lacZ_df), lacZ_df)
cell <- rownames(lacZ_meta)
lacZ_meta <- add_column(lacZ_meta, cell, .before = "orig.ident")
lacZ_df <- left_join(lacZ_df, lacZ_meta)
lacZ_df_sub <- data.frame(cell = lacZ_df$cell,
                GEM = lacZ_df$GEM,
                cluster = lacZ_df$seurat_clusters,
                Pseudo = lacZ_df$pseudo)
lacZ_df_sub <- na.omit(lacZ_df_sub)

Reln_df <- as.data.frame(t(as.data.frame(GetAssayData(object = Reln,
                                                      slot = "data",
                                                      assay = "RNA"))))
Reln_df <- cbind(cell = rownames(Reln_df), Reln_df)
cell <- rownames(Reln_meta)
Reln_meta <- add_column(Reln_meta, cell, .before = "orig.ident")
Reln_df <- left_join(Reln_df, Reln_meta)
Reln_df_sub <- data.frame(cell = Reln_df$cell,
                          GEM = Reln_df$GEM,
                          cluster = Reln_df$seurat_clusters,
                          Pseudo = Reln_df$pseudo)
Reln_df_sub <- na.omit(Reln_df_sub)

write.csv(lacZ_df_sub, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 6 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/lacZ_df_sub_noReln.csv")
write.csv(Reln_df_sub, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 6 - Reln KD snRNA-seq/Data&Code/Drd1s/Pseudotime/Reln_df_sub_noReln.csv")



# sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.3.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] shiny_1.8.0                 lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
# [5] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [9] tidyverse_2.0.0             dplyr_1.1.4                 magrittr_2.0.3              patchwork_1.2.0            
# [13] ggplot2_3.5.0               monocle3_1.3.5              SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [17] GenomicRanges_1.54.1        GenomeInfoDb_1.38.7         IRanges_2.36.0              S4Vectors_0.40.2           
# [21] MatrixGenerics_1.14.0       matrixStats_1.2.0           Biobase_2.62.0              BiocGenerics_0.48.1        
# [25] SeuratData_0.2.2.9001       Seurat_5.0.2                SeuratObject_5.0.1          sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.3.1           later_1.3.2             bitops_1.0-7            polyclip_1.10-6        
# [6] fastDummies_1.7.3       lifecycle_1.0.4         globals_0.16.3          lattice_0.22-5          MASS_7.3-60.0.1        
# [11] plotly_4.10.4           sass_0.4.8              jquerylib_0.1.4         httpuv_1.6.14           sctransform_0.4.1      
# [16] spam_2.10-0             spatstat.sparse_3.0-3   reticulate_1.35.0       cowplot_1.1.3           pbapply_1.7-2          
# [21] minqa_1.2.6             RColorBrewer_1.1-3      abind_1.4-5             zlibbioc_1.48.0         Rtsne_0.17             
# [26] RCurl_1.98-1.14         rappdirs_0.3.3          GenomeInfoDbData_1.2.11 ggrepel_0.9.5           irlba_2.3.5.1          
# [31] listenv_0.9.1           spatstat.utils_3.0-4    goftest_1.2-3           RSpectra_0.16-1         spatstat.random_3.2-3  
# [36] fitdistrplus_1.1-11     parallelly_1.37.1       leiden_0.4.3.1          codetools_0.2-19        DelayedArray_0.28.0    
# [41] tidyselect_1.2.0        farver_2.1.1            viridis_0.6.5           lme4_1.1-35.1           spatstat.explore_3.2-6 
# [46] jsonlite_1.8.8          ellipsis_0.3.2          progressr_0.14.0        ggridges_0.5.6          survival_3.5-8         
# [51] systemfonts_1.0.6       Matrix.utils_0.9.8      tools_4.3.1             ragg_1.2.7              ica_1.0-3              
# [56] Rcpp_1.0.12             glue_1.7.0              gridExtra_2.3           SparseArray_1.2.4       withr_3.0.0            
# [61] fastmap_1.1.1           boot_1.3-30             fansi_1.0.6             digest_0.6.34           timechange_0.3.0       
# [66] R6_2.5.1                mime_0.12               textshaping_0.3.7       colorspace_2.1-0        scattermore_1.2        
# [71] tensor_1.5              spatstat.data_3.0-4     utf8_1.2.4              generics_0.1.3          data.table_1.15.4      
# [76] httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.2.1          uwot_0.1.16             pkgconfig_2.0.3        
# [81] gtable_0.3.4            lmtest_0.9-40           XVector_0.42.0          htmltools_0.5.7         dotCall64_1.1-1        
# [86] scales_1.3.0            png_0.1-8               rstudioapi_0.15.0       tzdb_0.4.0              reshape2_1.4.4         
# [91] nlme_3.1-164            nloptr_2.0.3            proxy_0.4-27            zoo_1.8-12              cachem_1.0.8           
# [96] KernSmooth_2.23-22      parallel_4.3.1          miniUI_0.1.1.1          pillar_1.9.0            grid_4.3.1             
# [101] vctrs_0.6.5             RANN_2.6.1              promises_1.2.1          xtable_1.8-4            cluster_2.1.6          
# [106] cli_3.6.2               compiler_4.3.1          rlang_1.1.3             crayon_1.5.2            grr_0.9.5              
# [111] future.apply_1.11.1     labeling_0.4.3          plyr_1.8.9              stringi_1.8.3           viridisLite_0.4.2      
# [116] deldir_2.0-4            assertthat_0.2.1        munsell_0.5.0           lazyeval_0.2.2          spatstat.geom_3.2-9    
# [121] Matrix_1.6-5            RcppHNSW_0.6.0          hms_1.1.3               future_1.33.1           ROCR_1.0-11            
# [126] igraph_2.0.2            memoise_2.0.1           bslib_0.6.1 
