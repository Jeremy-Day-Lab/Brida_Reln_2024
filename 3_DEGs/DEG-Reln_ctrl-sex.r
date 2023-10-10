####################################
##                                ##
## DIFFERENTIALLY EXPRESSED GENES ##
##      BETWEEN LACZ & RELN       ##
##        CONTROL FOR SEX         ##
##                                ##
####################################

# This workflow is based off of work from Lara Ianov

library(Seurat)
library(ggplot2)
library(Libra)
library(dplyr)
library(stringr)
library(ComplexUpset)
library(Matrix.utils)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

# Load combo object
allRats <- readRDS("allRats_souped_noDoub.rds")

# Make sure cellType is correct after integrating & DoubletFinder
# DESeq2 code below will error if cellType still includes NAs
allRats$cellType <- Idents(allRats)

# Group isn't relevant here but add it so that code downstream to create metadata works
allRats$group <- allRats$GEM
allRats$group[allRats$group == 2] <- 1
allRats$group[allRats$group == 3] <- 1
allRats$group[allRats$group == 4] <- 1
allRats$group[allRats$group == 5] <- 2
allRats$group[allRats$group == 6] <- 2
allRats$group[allRats$group == 7] <- 2
allRats$group[allRats$group == 8] <- 2
allRats$group <- as.factor(allRats$group)

# Create a sample ID column: dataset_sex_stim
allRats$sample.id <- as.factor(paste(allRats$GEM, allRats$sex_target, sep = "_"))
# Named vector of sample names
sample_ids <- purrr::set_names(levels(allRats$sample.id))

# Get cell and sample metrics for aggregation
allRats@meta.data$cellType <- as.factor(allRats@meta.data$cellType)
cluster_names <- purrr::set_names(levels(allRats@meta.data$cellType))

# Number of clusters
cluster_total <- length(cluster_names)
# Total number of samples
sample_total <- length(sample_ids)
# Number of cells in each group
table(allRats$sample.id)

###### Count aggregation to sample level

# Make list of all cells: vars = cellType, sample.id
groups <- allRats@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(allRats@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

# Original code below to create raw_counts_list fails with this data,
# because 2 clusters don't contain cells from all 8 datasets (Drd3-MSN missing 3 & 4; Pvalb-Int missing 5).
# This is a cumbersome workaround, but it works.

# Make lists for subsets:
keep_rows_astro <- c("Astrocyte_1_male_lacz", "Astrocyte_2_male_reln", "Astrocyte_3_female_lacz", "Astrocyte_4_female_reln", "Astrocyte_5_female_lacz", "Astrocyte_6_female_reln", "Astrocyte_7_male_reln", "Astrocyte_8_male_lacz")
keep_rows_Drd1 <- c("Drd1-MSN_1_male_lacz", "Drd1-MSN_2_male_reln", "Drd1-MSN_3_female_lacz", "Drd1-MSN_4_female_reln", "Drd1-MSN_5_female_lacz", "Drd1-MSN_6_female_reln", "Drd1-MSN_7_male_reln", "Drd1-MSN_8_male_lacz")
keep_rows_Drd2.1 <- c("Drd2-MSN-1_1_male_lacz", "Drd2-MSN-1_2_male_reln", "Drd2-MSN-1_3_female_lacz", "Drd2-MSN-1_4_female_reln", "Drd2-MSN-1_5_female_lacz", "Drd2-MSN-1_6_female_reln", "Drd2-MSN-1_7_male_reln", "Drd2-MSN-1_8_male_lacz")
keep_rows_Drd2.2 <- c("Drd2-MSN-2_1_male_lacz", "Drd2-MSN-2_2_male_reln", "Drd2-MSN-2_3_female_lacz", "Drd2-MSN-2_4_female_reln", "Drd2-MSN-2_5_female_lacz", "Drd2-MSN-2_6_female_reln", "Drd2-MSN-2_7_male_reln", "Drd2-MSN-2_8_male_lacz")
keep_rows_Drd3 <- c("Drd3-MSN_1_male_lacz", "Drd3-MSN_2_male_reln", "Drd3-MSN_3_female_lacz", "Drd3-MSN_4_female_reln", "Drd3-MSN_5_female_lacz", "Drd3-MSN_6_female_reln", "Drd3-MSN_7_male_reln", "Drd3-MSN_8_male_lacz")
keep_rows_grm8 <- c("Grm8-MSN_1_male_lacz", "Grm8-MSN_2_male_reln", "Grm8-MSN_3_female_lacz", "Grm8-MSN_4_female_reln", "Grm8-MSN_5_female_lacz", "Grm8-MSN_6_female_reln", "Grm8-MSN_7_male_reln", "Grm8-MSN_8_male_lacz")
keep_rows_gaba <- c("GABA-Undef_1_male_lacz", "GABA-Undef_2_male_reln", "GABA-Undef_3_female_lacz", "GABA-Undef_4_female_reln", "GABA-Undef_5_female_lacz", "GABA-Undef_6_female_reln", "GABA-Undef_7_male_reln", "GABA-Undef_8_male_lacz")
keep_rows_pvalb <- c("Pvalb-Int_1_male_lacz", "Pvalb-Int_2_male_reln", "Pvalb-Int_3_female_lacz", "Pvalb-Int_4_female_reln", "Pvalb-Int_5_female_lacz", "Pvalb-Int_6_female_reln", "Pvalb-Int_7_male_reln", "Pvalb-Int_8_male_lacz")
keep_rows_sst <- c("Sst-Int_1_male_lacz", "Sst-Int_2_male_reln", "Sst-Int_3_female_lacz", "Sst-Int_4_female_reln", "Sst-Int_5_female_lacz", "Sst-Int_6_female_reln", "Sst-Int_7_male_reln", "Sst-Int_8_male_lacz")
keep_rows_micro <- c("Microglia_1_male_lacz", "Microglia_2_male_reln", "Microglia_3_female_lacz", "Microglia_4_female_reln", "Microglia_5_female_lacz", "Microglia_6_female_reln", "Microglia_7_male_reln", "Microglia_8_male_lacz")
keep_rows_mural <- c("Mural_1_male_lacz", "Mural_2_male_reln", "Mural_3_female_lacz", "Mural_4_female_reln", "Mural_5_female_lacz", "Mural_6_female_reln", "Mural_7_male_reln", "Mural_8_male_lacz")
keep_rows_olig <- c("Olig_1_male_lacz", "Olig_2_male_reln", "Olig_3_female_lacz", "Olig_4_female_reln", "Olig_5_female_lacz", "Olig_6_female_reln", "Olig_7_male_reln", "Olig_8_male_lacz")
keep_rows_poly <- c("Polydend_1_male_lacz", "Polydend_2_male_reln", "Polydend_3_female_lacz", "Polydend_4_female_reln", "Polydend_5_female_lacz", "Polydend_6_female_reln", "Polydend_7_male_reln", "Polydend_8_male_lacz")

# Create subsets
subset_astro <- t(count_aggr[rownames(count_aggr) %in% keep_rows_astro, ])
subset_Drd1 <- t(count_aggr[rownames(count_aggr) %in% keep_rows_Drd1, ])
subset_Drd2.1 <- t(count_aggr[rownames(count_aggr) %in% keep_rows_Drd2.1, ])
subset_Drd2.2 <- t(count_aggr[rownames(count_aggr) %in% keep_rows_Drd2.2, ])
subset_Drd3 <- t(count_aggr[rownames(count_aggr) %in% keep_rows_Drd3, ])
subset_grm8 <- t(count_aggr[rownames(count_aggr) %in% keep_rows_grm8, ])
subset_gaba <- t(count_aggr[rownames(count_aggr) %in% keep_rows_gaba, ])
subset_pvalb <- t(count_aggr[rownames(count_aggr) %in% keep_rows_pvalb, ])
subset_sst <- t(count_aggr[rownames(count_aggr) %in% keep_rows_sst, ])
subset_micro <- t(count_aggr[rownames(count_aggr) %in% keep_rows_micro, ])
subset_mural <- t(count_aggr[rownames(count_aggr) %in% keep_rows_mural, ])
subset_olig <- t(count_aggr[rownames(count_aggr) %in% keep_rows_olig, ])
subset_poly <- t(count_aggr[rownames(count_aggr) %in% keep_rows_poly, ])

# Combine subsets into a list
raw_counts_list <- list(subset_astro, subset_Drd1, subset_Drd2.1, subset_Drd2.2, subset_Drd3, subset_grm8, subset_gaba, subset_pvalb, subset_sst, subset_micro, subset_mural, subset_olig, subset_poly)
# Add sublist names
names(raw_counts_list) <- c("Astrocyte", "Drd1.MSN", "Drd2.MSN.1", "Drd2.MSN.2", "Drd3.MSN", "Grm8.MSN", "GABA.Undef", "Pvalb.Int", "Sst.Int", "Microglia", "Mural", "Olig", "Polydend")

# Replace - with . in cluster names only
colnames(raw_counts_list[[2]]) <- c("Drd1.MSN_1_male_lacz","Drd1.MSN_2_male_reln","Drd1.MSN_3_female_lacz","Drd1.MSN_4_female_reln","Drd1.MSN_5_female_lacz","Drd1.MSN_6_female_reln","Drd1.MSN_7_male_reln","Drd1.MSN_8_male_lacz")
colnames(raw_counts_list[[3]]) <- c("Drd2.MSN.1_1_male_lacz","Drd2.MSN.1_2_male_reln","Drd2.MSN.1_3_female_lacz","Drd2.MSN.1_4_female_reln","Drd2.MSN.1_5_female_lacz","Drd2.MSN.1_6_female_reln","Drd2.MSN.1_7_male_reln","Drd2.MSN.1_8_male_lacz")
colnames(raw_counts_list[[4]]) <- c("Drd2.MSN.2_1_male_lacz","Drd2.MSN.2_2_male_reln","Drd2.MSN.2_3_female_lacz","Drd2.MSN.2_4_female_reln","Drd2.MSN.2_5_female_lacz","Drd2.MSN.2_6_female_reln","Drd2.MSN.2_7_male_reln","Drd2.MSN.2_8_male_lacz")
colnames(raw_counts_list[[5]]) <- c("Drd3.MSN_1_male_lacz","Drd3.MSN_2_male_reln","Drd3.MSN_5_female_lacz","Drd3.MSN_6_female_reln","Drd3.MSN_7_male_reln","Drd3.MSN_8_male_lacz")
colnames(raw_counts_list[[6]]) <- c("Grm8.MSN_1_male_lacz","Grm8.MSN_2_male_reln","Grm8.MSN_3_female_lacz","Grm8.MSN_4_female_reln","Grm8.MSN_5_female_lacz","Grm8.MSN_6_female_reln","Grm8.MSN_7_male_reln","Grm8.MSN_8_male_lacz")
colnames(raw_counts_list[[7]]) <- c("GABA.Undef_1_male_lacz","GABA.Undef_2_male_reln","GABA.Undef_3_female_lacz","GABA.Undef_4_female_reln","GABA.Undef_5_female_lacz","GABA.Undef_6_female_reln","GABA.Undef_7_male_reln","GABA.Undef_8_male_lacz")
colnames(raw_counts_list[[8]]) <- c("Pvalb.Int_1_male_lacz","Pvalb.Int_2_male_reln","Pvalb.Int_3_female_lacz","Pvalb.Int_4_female_reln","Pvalb.Int_6_female_reln","Pvalb.Int_7_male_reln","Pvalb.Int_8_male_lacz")
colnames(raw_counts_list[[9]]) <- c("Sst.Int_1_male_lacz","Sst.Int_2_male_reln","Sst.Int_3_female_lacz","Sst.Int_4_female_reln","Sst.Int_5_female_lacz","Sst.Int_6_female_reln","Sst.Int_7_male_reln","Sst.Int_8_male_lacz")

rm(list=ls(pattern="subset_"))
rm(list=ls(pattern="keep_rows_"))

#### This is the code that doesn't work with this data...
# Turn into a list, and split list into components for each cluster
# Transform so that in each cluster: rows = genes; cols = cluster_sampleID
# raw_counts_list <- split.data.frame(
#   count_aggr,
#   #  factor(cluster_names) # Changing list names from - to . is optional
#   factor(gsub("-", ".", cluster_names)) # Changing list names from - to . is optional
# ) %>%
#   lapply(function(x) {
#     magrittr::set_colnames(t(x), gsub("-", ".", rownames(x)))
#   })

######### START AGAIN HERE ###########

##### Sample-level metadata
get_sample_ids <- function(x){
  raw_counts_list[[x]] %>%
    colnames()
}

de_samples <- purrr::map(1:length(cluster_names), get_sample_ids) %>%
  unlist()

# Create metadata dataframe with sample IDs, cluster IDs, and condition
metadata <- data.frame(cluster_id = sub("_.*","", de_samples),
                       sample_id  = de_samples,
                       sex.target = sub("^[^_]*_","", de_samples))
#Create extra columns in the metadata
metadata$GEM     <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",1)))
metadata$sex       <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",2)))
metadata$target    <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",3)))

#################
# NEW FUNCTIONS #
#################

deter_direc <- function(x) {
  log_names <- c("avg_logFC", "avg_log2FC", "log2FC", "log2FoldChange")
  log_select <- x %>% select(any_of(log_names))
  ifelse(log_select <= 0, "down", ifelse(log_select >= 0, "up", "no_change"))
}

fromList <- function (input) {
  # Same as original UpSetR::fromList(), but modified as shown in https://github.com/hms-dbmi/UpSetR/issues/85
  # Thanks to @docmanny
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # Except now it conserves your original value names! (in this case gene names)
  row.names(data) <- elements
  return(data)
}

plot_output <- function(p, file_name, w_png=700, h_png=600, w_pdf=12, h_pdf=8, show_plot = TRUE){
  
  png(paste0(file_name,".png"), width = w_png, height = h_png)
  plot(eval(p))
  dev.off()
  
  pdf(paste0(file_name,".pdf"), width = w_pdf, height = h_pdf)
  plot(eval(p))
  dev.off()
  
  if (show_plot) {
    plot(eval(plot))
  }
}

# for pheatmaps as it uses grid system instead of ggplot
pheatmap_output <- function(x, file_name, w_png=900, h_png=700, w_pdf=12, h_pdf=8) {
  
  png(paste0(file_name,".png"), width = w_png, height = h_png)
  
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  
  pdf(paste0(file_name,".pdf"), width = w_pdf, height = h_pdf)
  
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

nested_lapply <- function(data, FUN) {
  lapply(data, function(sublist) { lapply(sublist, FUN) })
}

#####################################
# Create a DESEqDataSet from object #
#####################################

dds_all_cells <- mapply(FUN = function(x, y, z) {
  
  # Subset metadata by cell type being tested
  cluster_metadata <- subset(metadata, cluster_id == z)
  
  # Assign rownames of metadata to be sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  print(head(cluster_metadata))
  
  # Subset raw counts to cell type being tested
  counts <- raw_counts_list[[y]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  print(head(cluster_counts))
  
  # Check that all row names of metadata are the same and in the same order as cols of counts
  stopifnot(all(rownames(cluster_metadata) == colnames(cluster_counts)))
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sex + target)
  
  # Set reference to lacz
  dds$target <- relevel(dds$target, ref = "lacz")
  return(dds)
  
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = names(raw_counts_list))

####################################
# Remove genes with mean counts <5 #
####################################

dds_all_cells <- mapply(FUN = function(x, z) {
  
  #----- Counts Pre-filtering based on rowMeans -------
  # Here rowMeans is used, but can also use rowSums or min. per sample
  message(paste0("Number of genes before pre-filtering ", z, ": ",  nrow(counts(x))))
  keep <- rowMeans(counts(x)) > 5
  x <- x[keep,]
  message(paste0("Number of genes after pre-filtering ", z, ": ",  nrow(counts(x))))
  return(x)
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

# Check to see how many rows with counts > 5
mapply(FUN = function(x){
  message(paste0(nrow(counts(x))))
}, x = dds_all_cells)

###################
# QUALITY CONTROL #
###################

# Make output folder
dir.create("DESeq2_QC")

mapply(FUN = function(x, z) {
  vsd <- varianceStabilizingTransformation(x, blind = FALSE)
  vsd
  
  # PCA
  #-----------------Target-------------------
  data <- plotPCA(vsd, intgroup = c("target"), returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  pca.plot.target <- ggplot(data, aes(PC1, PC2, color = target)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  #-----------------Sex---------------------
  data <- plotPCA(vsd, intgroup = c("sex"), returnData = TRUE)
  
  pca.plot.sex <- ggplot(data, aes(PC1, PC2, color = sex)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  #-----------------Sex.Target------------------
  data <- plotPCA(vsd, intgroup = c("sex.target"), returnData = TRUE)
  
  pca.plot.sex.target <- ggplot(data, aes(PC1, PC2, color = sex.target)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  all_pca <- pca.plot.target + pca.plot.sex + pca.plot.sex.target
  file_name_pca <- paste0("DESeq2_QC/PCA_", z)
  plot_output(all_pca, file_name = file_name_pca, 
              w_png = 1800, w_pdf = 20, h_png = 400, h_pdf = 6, show_plot = FALSE)
  
  # Sample-to-samples distances heatmap:
  # Dist computes distance with Euclidean method
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  heatmap_dist <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors
  )
  
  file_name_heatmap <- paste0("DESeq2_QC/sample_to_sample_heatmap_", z)
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

####################
## CALCULATE DEGs ##
####################

# Make output folder
dir.create("DESeq2_DEGResults")

dds_all_cells <- mapply(FUN = function(x, z) {
  x <- DESeq(x, test="LRT", reduced = ~ sex)
  
  # Checking coeff by resultsNames(dds):
  resultsNames(x)
  
  # Dispersion plot
  file_name_dist <- paste0("DESeq2_DEGResults/DispPlot_", z, ".pdf")
  pdf(file_name_dist)
  plotDispEsts(x)
  dev.off()
  return(x)
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

# Calculate results
res_all_cells <- mapply(FUN = function(x) {
  res <- results(x) # ref level established earlier by relevel
  return(res)
}, x = dds_all_cells)

# Calculate percent cells expressing gene
Idents(allRats) <- gsub(Idents(allRats), pattern = "-", replacement = ".")
seuratCounts <- t(GetAssayData(object = allRats, slot = "data", assay = "RNA"))

for(i in names(res_all_cells)){
  
  print(paste("Beginning:", nrow(res_all_cells[[i]])))
  res_all_cells[[i]] <- as.data.frame(res_all_cells[[i]])
  
  # Create a gene name column
  res_all_cells[[i]]$gene <- row.names(res_all_cells[[i]])
  
  # Pull genes and calculate percent expressing
  genesPctExp <- as.data.frame(colMeans(seuratCounts[WhichCells(object = allRats, idents = i), res_all_cells[[i]]$gene]>0)*100)
  
  # Change column names
  genesPctExp$gene <- row.names(genesPctExp)
  colnames(genesPctExp)[1] <- "Pct_Expressing"
  
  # Add percent expressing to dataframe
  res_all_cells[[i]] <- merge(x = res_all_cells[[i]],
                              y = genesPctExp,
                              by = "gene")
  
  print(paste("End:", nrow(res_all_cells[[i]])))
}

# Make output folder
dir.create("DESeq2_Reln_Collapsed")

# Write out the files
mapply(FUN = function(x, z) {
  file_name <- paste0("DESeq2_Reln_Collapsed/DESeq2_res_", z, ".csv")
  write.csv(x, file=file_name, row.names = TRUE, quote = FALSE)
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

# Save counts
mapply(FUN = function(x, z) {
  file_name_norm <- paste0("DESeq2_Reln_Collapsed/normalized_counts_", z, ".csv")
  normalized_counts <- as.data.frame(counts(x, normalized = TRUE))
  write.csv(normalized_counts, file = file_name_norm, row.names = TRUE, quote = FALSE)
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

# Make table of number of DEGs
# NOTE: Change nrow to your correct number of names in dds_all_cells (i.e., # of clusters)
DEGs_DF <- as.data.frame(matrix(ncol = 2, nrow=13))
colnames(DEGs_DF) <- c("cellType", "number_of_DEGs")
rownames(DEGs_DF) <- names(dds_all_cells)
DEGs_DF$cellType  <- names(dds_all_cells)

for(i in names(dds_all_cells)){
  x <- as.data.frame(res_all_cells[[i]])
  DEGs_DF[i,"number_of_DEGs"] <- nrow(subset(x, subset=(padj <= 0.05)))
}

write.table(x         = DEGs_DF,
            file      = "DESeq2_Reln_Collapsed/number_of_DEGs_Per_cellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)


