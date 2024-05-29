#Calculate gini coefficients for all clusters
library(reldist)
library(Seurat)
library(ggplot2)
library(scales)
library(ggrepel)
library(DescTools)
library(Libra)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexUpset)
library(Matrix.utils)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

### ------------- Gini Coefficient Analysis -------------------------------------------------

#Load the NAc_Combo object
Drd1 <- readRDS(file = "/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/D1_subset_of_MCNobj.rds")

DefaultAssay(Drd1) <- "RNA"
Reln <- FeaturePlot(Drd1, "Reln",  cols = c("#dbdbdb", "#ed1c24"))
Reln + coord_fixed(ratio = 1)  #will keep ratios consistent between UMAPs & featureplots

#Other FPs
Fos <- FeaturePlot(Drd1, "Fos", cols = c("#dbdbdb", "#ed1c24"))
Fos + coord_fixed(ratio = 1) 
Fosb <- FeaturePlot(Drd1, "Fosb", cols = c("#dbdbdb", "#ed1c24"))
Fosb + coord_fixed(ratio = 1)
Fosl2 <- FeaturePlot(Drd1, "Fosl2", cols = c("#dbdbdb", "#ed1c24"))
Fosl2 + coord_fixed(ratio = 1) 
Nr4a1 <- FeaturePlot(Drd1, "Nr4a1", cols = c("#dbdbdb", "#ed1c24"))
Nr4a1 + coord_fixed(ratio = 1) 
Nr4a3 <- FeaturePlot(Drd1, "Nr4a3", cols = c("#dbdbdb", "#ed1c24"))
Nr4a3 + coord_fixed(ratio = 1) 


#Calculate the cluster-specific average expression values for every gene
DefaultAssay(Drd1) <- "RNA"
CellAverages <- AverageExpression(Drd1)

#Pull the RNA values
RNA <- CellAverages$RNA

#Calculate the Gini coefficient, maximum/mean gene expression value, and ratio
Gini_calc <- apply(RNA, 1, Gini)
max_RNA <- apply(RNA, 1, max)
mean_RNA <- apply(RNA, 1, mean)
Ratio <- max_RNA/mean_RNA

#Create a dataframe consisting of the Gini coefficient, maximum/mean gene expression value, and ratio
Gini_dataframe <- data.frame(Gini_calc, max_RNA, mean_RNA, Ratio, RNA)
#Create a column of gene names
Gini_dataframe$Gene <- rownames(Gini_dataframe)


#Remove any NA Gini_calc values
Gini_dataframe_noNA <- Gini_dataframe[-which(is.na(Gini_dataframe$Gini_calc)),]
colnames(Gini_dataframe_noNA)[5] <- "0"
colnames(Gini_dataframe_noNA)[6] <- "1"
colnames(Gini_dataframe_noNA)[7] <- "2"
colnames(Gini_dataframe_noNA)[8] <- "3"
colnames(Gini_dataframe_noNA)[9] <- "4"


#Change all of the - to .
colnames(Gini_dataframe_noNA)[grep(pattern = "-",x = colnames(Gini_dataframe_noNA))] <- gsub(colnames(Gini_dataframe_noNA)[grep(pattern = "-",x = colnames(Gini_dataframe_noNA))],pattern = "-",replacement = ".")

#Write out the gini data.frame
write.table(x         = Gini_dataframe_noNA,
            file      = "/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Gini_Dataframe_Drd1.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep       = "\t",
            quote     = FALSE)
 

#Create a sample id column that is Dataset_Sex_Stim
Drd1$sample.id <- as.factor(paste(Drd1$Dataset,Drd1$Sex_Stim,sep = "_"))

#Get cell and sample metrics for aggregation
Drd1@meta.data$Combo_CellType <- as.character(Drd1@meta.data$seurat_clusters)
Drd1@meta.data$CellType <- as.factor(Drd1@meta.data$Combo_CellType)
cell_names <- purrr::set_names(levels(Drd1@meta.data$CellType))
cell_names
#"0" "1" "2" "3" "4"

#Number of clusters
cluster_total <- length(cell_names)
cluster_total
# [5]

# Named vector of sample names
Drd1$Sex_Stim <- factor(Drd1$Sex_Stim)
sample_ids <- purrr::set_names(levels(Drd1@meta.data$Sex_Stim))
sample_ids
# Fem_Coc    Fem_Sal   Male_Coc   Male_Sal 


# Total number of samples 
sample_total <- length(sample_ids)
sample_total
#4

#Figure out how many cells in each main group
table(Drd1@meta.data$Sex_Stim)
# Fem_Coc  Fem_Sal Male_Coc Male_Sal 
# 2259     1482     1887     1674 

##########Count aggregation to sample level
groups <- Drd1@meta.data[, c("CellType", "sample.id")]

# Aggregate across cluster-sample groups (raw counts, thus counts slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(Drd1@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

dim(count_aggr)
# [1]    40 30560

#Transpose count_aggr
count_aggr_t <- t(count_aggr)

#change the count_aggr_t
colnames(count_aggr_t) <- gsub(x = colnames(count_aggr_t),pattern = "-",replacement = ".")


# Create a data frame with the sample IDs, cluster IDs and condition
metadata <- data.frame(cluster_id = sub("_.*","", colnames(count_aggr_t)),
                       sample_id = colnames(count_aggr_t),
                       Sex.Stim = sub("^[^_]*_","", colnames(count_aggr_t)))
#Create extra columns in the metadata
metadata$dataset <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",1)))
metadata$Sex     <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",2)))
metadata$Stim    <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",3)))

counts_mat <- Drd1@assays$RNA@counts
Idents(Drd1) <- gsub(x = Idents(Drd1),pattern = "-",replacement = ".")
for(i in unique(metadata$cluster_id)){
  ############Make the DESEq object############
  metadata$CellType <- ifelse(metadata$cluster_id == i,
                              i,
                              "Other")
  dds <- DESeqDataSetFromMatrix(count_aggr_t, 
                                colData = metadata, 
                                design = ~ dataset + Stim + CellType)
  dds$CellType <- relevel(dds$CellType, ref = "Other")
  keep <- rowMeans(counts(dds)) > 5
  dds <- dds[keep,]
  print(paste("dds made for",i))
  ############Quality control############
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  #PCA for major celltype
  data <- plotPCA(vsd, intgroup = c("CellType"), returnData = TRUE)
  data$cluster_id <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",1))
  percentVar <- round(100 * attr(data, "percentVar"))
  pca.plot.CellType <- ggplot(data, aes(PC1, PC2, color = CellType)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  #PCA for cluster_id
  pca.plot.cluster_id <- ggplot(data, aes(PC1, PC2, color = cluster_id)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  ggsave(filename = paste0("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/",i,".pdf"),
         plot     = cowplot::plot_grid(plotlist = list(pca.plot.CellType,pca.plot.cluster_id),ncol = 2),
         width = 12,
         height = 7)
  print(paste("PCA plots made for",i))
  ############Get results############
  dds <- DESeq(dds, test="LRT", reduced = ~ dataset + Stim)
  res <- as.data.frame(results(dds))
  res$GeneName <- rownames(res)
  cells.1    <- WhichCells(Drd1,idents = i)
  cells.2    <- setdiff(x = row.names(Drd1@meta.data),y=cells.1)
  res$Pct_Expressing <- NA
  res$Pct_CellType_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100
  res$Pct_Other_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.2,drop = FALSE] > 0) / length(cells.2))*100
  write.table(file      = paste0("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/",i,".txt"),
              x         = res,
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE)
  print(paste("DEG lists written for",i))
  rm(dds,keep,vsd,data,percentVar,pca.plot.CellType,pca.plot.cluster_id,res)
}




colors <- hue_pal()(5)
names(colors) <- gsub(x = levels(Drd1$Combo_CellType),pattern = "-",replacement = ".")

Idents(Drd1) <- Drd1$Combo_CellType
#Make Gini plots
for(i in colnames(Gini_dataframe_noNA)[5:9]){
  print(i)
  #Pull the DEG list
  x <- read.table(file = paste0("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/",i,".txt"), sep = "\t",header = TRUE)
  #Calculate the percent of cells expressing 
  celltype_name <- gsub(pattern = "[.]",replacement = "-",x = i)
  cells.1    <- WhichCells(Drd1,idents = celltype_name)
  counts_mat <- Drd1@assays$RNA@counts
  x$Pct_Expressing <- NA
  x$Pct_Expressing <- (rowSums(x = counts_mat[x$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100
  #Merge the dataframes 
  y <- merge(x = x,
             y = Gini_dataframe_noNA,
             by.x = "GeneName",
             by.y = "Gene")
  y$logCluster <- log10(y[,i])
  y <- y[!is.infinite(y$logCluster),]
  

  #Make the plot
  Gini_plot <- ggplot(data = y,aes(x = logCluster,y = Gini_calc,size = Pct_Expressing)) +
    geom_point(color = "lightgrey",aes(size = Pct_Expressing)) +
    geom_point(data = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0)),
               color = colors[i]) +
    geom_text_repel(data = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))[order(subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))$Gini_calc,decreasing = TRUE)[1:25],],
                    label = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0 ))[order(subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))$Gini_calc,decreasing = TRUE)[1:25],"GeneName"],
                    show.legend = FALSE,
                    size = 6) +
    geom_density2d(color = "red",show.legend = FALSE) +
    ylim(c(0,1)) +
    theme_bw() + 
    NoGrid() +
    labs(x = paste0("log10(",i," Counts)"),
         y = "Gini coefficient",
         size = "% of Cells in Cluster\nExpressing Gene") +
    ggtitle(label = i) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot     = Gini_plot,
         filename = paste0("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Gini/",i,"_GiniPlot.pdf"),
         height   = 12,
         width    = 12)
  print(paste(i,"is complete"))
}

write.csv(y, "/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Gini/y.csv" )



#read in csvs & add cell type
D1_0 <- read.table("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/0.txt", header = TRUE)
D1_0$Cluster <- 0
D1_1 <- read.table("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/1.txt", header = TRUE)
D1_1$Cluster <- 1
D1_2 <- read.table("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/2.txt", header = TRUE)
D1_2$Cluster <- 2
D1_3 <- read.table("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/3.txt", header = TRUE)
D1_3$Cluster <- 3
D1_4 <- read.table("/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/Drd1_Subcluster_DESeq2/4.txt", header = TRUE)
D1_4$Cluster <- 4

D1_dds <- rbind(D1_0, D1_1, D1_2, D1_3, D1_4)
rm(D1_0, D1_1, D1_2, D1_3, D1_4)
Reln_only <- subset(D1_dds, subset = D1_dds$GeneName == "Reln")

#DotPlot of Reln Across Populations
ggplot(data = Reln_only,aes(x = GeneName, y = Cluster)) +
  geom_point(aes(color = Pct_CellType_Expressing, size = log2FoldChange)) +
  scale_color_gradient(low = "#999999", high = "dodger blue") +
  theme_bw() +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5))

#EES
EES_meta <- read.csv("/Volumes/Day-Lab$/Cathy/EES analysis of MCN2023 object/Drd1 all/metadata_with_EES_D1.csv")
Drd1@meta.data$EES <- EES_meta$EES
FeaturePlot(Drd1, "EES", order = T, min.cutoff = "q3")

EES_byCluster <- EES_meta %>% group_by(cluster) %>%
  mutate(EES_mean = mean(EES))


#Get reln 
EES_Reln <- as.data.frame(t(GetAssayData(object = Drd1,slot = "data",assay = "RNA")))
EES_Reln_df <- data.frame(X = rownames(EES_Reln),
                       Reln = EES_Reln$Reln)
EES_Reln_df <- left_join(EES_meta, EES_Reln_df)
EES_df <- data.frame(Cell = EES_Reln_df$X,
                          Cluster = EES_Reln_df$cluster,
                          Sample = paste0(EES_Reln_df$Sex_Stim, EES_Reln_df$Dataset),
                          EES = EES_Reln_df$EES,
                          Reln = EES_Reln_df$Reln)
write.csv(EES_df, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 1 - snRNA-seq/Data&Code/EES&Subcluster/EES_Reln_byCell.csv")
#####RepeatedMeasuresRegression#####
library(rmcorr)
#get sample size
n <- length(unique(EES_df$Sample))

EES_df0 <- subset(EES_df, subset = EES_df$Cluster == "Drd1-0")
EES_df1 <- subset(EES_df, subset = EES_df$Cluster == "Drd1-1")
EES_df2 <- subset(EES_df, subset = EES_df$Cluster == "Drd1-2")
EES_df3 <- subset(EES_df, subset = EES_df$Cluster == "Drd1-3")
EES_df4 <- subset(EES_df, subset = EES_df$Cluster == "Drd1-4")

## Calculate rmcorr using selected columns
Cluster0 <- rmcorr(participant = EES_df0$Sample,
                   measure1 = EES_df0$Reln,
                   measure2 = EES_df0$EES,
                   dataset = EES_df0,
                   CI.level = 0.95,
                   CIs = "bootstrap",
                   nreps = 100,
                   bstrap.out = FALSE)

Cluster1 <- rmcorr(participant = EES_df1$Sample,
                 measure1 = EES_df1$Reln,
                 measure2 = EES_df1$EES,
                 dataset = EES_df1,
                 CI.level = 0.95,
                 CIs = "bootstrap",
                 nreps = 100,
                 bstrap.out = FALSE)

Cluster2 <- rmcorr(participant = EES_df2$Sample,
                   measure1 = EES_df2$Reln,
                   measure2 = EES_df2$EES,
                   dataset = EES_df2,
                   CI.level = 0.95,
                   CIs = "bootstrap",
                   nreps = 100,
                   bstrap.out = FALSE)

Cluster3 <- rmcorr(participant = EES_df3$Sample,
                   measure1 = EES_df3$Reln,
                   measure2 = EES_df3$EES,
                   dataset = EES_df3,
                   CI.level = 0.95,
                   CIs = "bootstrap",
                   nreps = 100,
                   bstrap.out = FALSE)

Cluster4 <- rmcorr(participant = EES_df4$Sample,
                   measure1 = EES_df4$Reln,
                   measure2 = EES_df4$EES,
                   dataset = EES_df4,
                   CI.level = 0.95,
                   CIs = "bootstrap",
                   nreps = 100,
                   bstrap.out = FALSE)

All_D1 <- rmcorr(participant = EES_df$Sample,
                   measure1 = EES_df$Reln,
                   measure2 = EES_df$EES,
                   dataset = EES_df,
                   CI.level = 0.95,
                   CIs = "bootstrap",
                   nreps = 100,
                   bstrap.out = FALSE)
Cluster0
Cluster1
Cluster2
Cluster3
Cluster4
All_D1

#Is this the best way to adjust p-values??
p.vals <- list(Cluster0[["p"]], Cluster1[["p"]], Cluster2[["p"]], Cluster3[["p"]],
               Cluster4[["p"]], All_D1[["p"]])

p.vals.fdr <- p.adjust(p.vals, 
                       method = "fdr",
                       n = length(p.vals))

## And plot the data
rm_plot <- ggplot(data = EES_df, aes(x=Reln, y=EES)) +
  geom_point(aes(color = Cluster))+
  geom_smooth(span=0.25, se=FALSE, aes(color = Cluster))+ #LOWESS smoothing
  xlim(0,5.1)+
  ylim(-0.2,0.3)
rm_plot


rm_plot_all <- ggplot(data = EES_df, aes(x=Reln, y=EES)) +
  geom_point(aes(color = Cluster))+
  geom_smooth(span=0.25, se=FALSE)+
  xlim(0,5.1)+
  ylim(-0.2,0.3)

rm_plot_all
lines(lowess(EES_df$Reln, EES_df$EES))


rm_plot_all+rm_plot

#Trying to figure out how to do the Lowess plot in R for suppl
lowess(Reln ~ EES, data=EES_df, col = Cluster)
plotLowess(Reln ~ EES, data=EES_df, col = Cluster)

sessionInfo()

# > sessionInfo()
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
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-3          DESeq2_1.40.2               SummarizedExperiment_1.30.2
# [5] Biobase_2.60.0              MatrixGenerics_1.12.3       matrixStats_1.0.0           GenomicRanges_1.52.0       
# [9] GenomeInfoDb_1.36.3         IRanges_2.34.1              S4Vectors_0.38.2            BiocGenerics_0.46.0        
# [13] Matrix.utils_0.9.8          Matrix_1.6-1.1              ComplexUpset_1.3.3          Libra_1.7                  
# [17] nnls_1.5                    DescTools_0.99.52           scales_1.2.1                SeuratObject_4.1.4         
# [21] Seurat_4.4.0                lubridate_1.9.3             forcats_1.0.0               stringr_1.5.0              
# [25] purrr_1.0.2                 readr_2.1.4                 tibble_3.2.1                tidyverse_2.0.0            
# [29] class_7.3-22                MASS_7.3-60                 biomaRt_2.56.1              VennDetail_1.16.0          
# [33] BiocManager_1.30.22         ggrepel_0.9.3               ggplot2_3.4.3               tidyr_1.3.0                
# [37] dplyr_1.1.3                
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21        splines_4.3.1           later_1.3.1             bitops_1.0-7            filelock_1.0.2         
# [6] cellranger_1.1.0        polyclip_1.10-6         XML_3.99-0.14           lifecycle_1.0.3         globals_0.16.2         
# [11] lattice_0.21-9          magrittr_2.0.3          limma_3.56.2            plotly_4.10.2           rmarkdown_2.25         
# [16] yaml_2.3.7              httpuv_1.6.11           sctransform_0.4.0       sp_2.1-0                spatstat.sparse_3.0-2  
# [21] gld_2.6.6               reticulate_1.32.0       cowplot_1.1.1           pbapply_1.7-2           DBI_1.1.3              
# [26] abind_1.4-5             zlibbioc_1.46.0         expm_0.999-8            Rtsne_0.16              RCurl_1.98-1.12        
# [31] rappdirs_0.3.3          GenomeInfoDbData_1.2.10 irlba_2.3.5.1           listenv_0.9.0           spatstat.utils_3.0-3   
# [36] goftest_1.2-3           spatstat.random_3.1-6   fitdistrplus_1.1-11     parallelly_1.36.0       DelayedArray_0.26.7    
# [41] leiden_0.4.3            codetools_0.2-19        xml2_1.3.5              tidyselect_1.2.0        futile.logger_1.4.3    
# [46] farver_2.1.1            BiocFileCache_2.8.0     spatstat.explore_3.2-3  jsonlite_1.8.7          e1071_1.7-14           
# [51] ellipsis_0.3.2          progressr_0.14.0        ggridges_0.5.4          survival_3.5-7          systemfonts_1.0.5      
# [56] tools_4.3.1             progress_1.2.2          ragg_1.2.6              ica_1.0-3               Rcpp_1.0.11            
# [61] glue_1.6.2              gridExtra_2.3           xfun_0.40               withr_2.5.1             formatR_1.14           
# [66] fastmap_1.1.1           boot_1.3-28.1           fansi_1.0.5             digest_0.6.33           timechange_0.2.0       
# [71] R6_2.5.1                mime_0.12               textshaping_0.3.7       colorspace_2.1-0        scattermore_1.2        
# [76] tensor_1.5              spatstat.data_3.0-1     RSQLite_2.3.1           UpSetR_1.4.0            utf8_1.2.3             
# [81] generics_0.1.3          data.table_1.14.8       S4Arrays_1.0.6          prettyunits_1.2.0       httr_1.4.7             
# [86] htmlwidgets_1.6.2       uwot_0.1.16             pkgconfig_2.0.3         Exact_3.2               gtable_0.3.4           
# [91] blob_1.2.4              lmtest_0.9-40           XVector_0.40.0          htmltools_0.5.6.1       lmom_3.0               
# [96] png_0.1-8               knitr_1.44              lambda.r_1.2.4          rstudioapi_0.15.0       tzdb_0.4.0             
# [101] reshape2_1.4.4          nlme_3.1-163            curl_5.1.0              proxy_0.4-27            cachem_1.0.8           
# [106] zoo_1.8-12              rootSolve_1.8.2.4       KernSmooth_2.23-22      parallel_4.3.1          miniUI_0.1.1.1         
# [111] AnnotationDbi_1.62.2    pillar_1.9.0            grid_4.3.1              vctrs_0.6.3             RANN_2.6.1             
# [116] promises_1.2.1          dbplyr_2.3.4            xtable_1.8-4            cluster_2.1.4           evaluate_0.22          
# [121] isoband_0.2.7           VennDiagram_1.7.3       locfit_1.5-9.8          mvtnorm_1.2-3           cli_3.6.1              
# [126] compiler_4.3.1          futile.options_1.0.1    rlang_1.1.1             crayon_1.5.2            grr_0.9.5              
# [131] future.apply_1.11.0     labeling_0.4.3          plyr_1.8.9              stringi_1.7.12          BiocParallel_1.34.2    
# [136] viridisLite_0.4.2       deldir_1.0-9            munsell_0.5.0           Biostrings_2.68.1       lazyeval_0.2.2         
# [141] spatstat.geom_3.2-5     hms_1.1.3               patchwork_1.1.3         bit64_4.0.5             future_1.33.0          
# [146] KEGGREST_1.40.1         shiny_1.7.5             ROCR_1.0-11             igraph_1.5.1            memoise_2.0.1          
# [151] bit_4.0.5               readxl_1.4.3 
