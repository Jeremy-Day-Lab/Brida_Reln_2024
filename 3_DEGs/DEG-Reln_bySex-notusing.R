####################################
##                                ##
## DIFFERENTIALLY EXPRESSED GENES ##
##      BETWEEN LACZ & RELN       ##
##      BY SEX (LRT & WALD)       ##
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

allRats <- readRDS("allRats_souped_noDoub.rds")

# Make sure cellType is correct after integrating & DoubletFinder
# DESeq2 code below will error if cellType still includes NAs
allRats$cellType <- Idents(allRats)

# LRT requires a reducing factor, so control for FACS group here
# FACS group 1 = S1-4; FACS group 2 = S5-8

allRats$group <- allRats$GEM
allRats$group[allRats$group == 2] <- 1
allRats$group[allRats$group == 3] <- 1
allRats$group[allRats$group == 4] <- 1
allRats$group[allRats$group == 5] <- 2
allRats$group[allRats$group == 6] <- 2
allRats$group[allRats$group == 7] <- 2
allRats$group[allRats$group == 8] <- 2
allRats$group <- as.factor(allRats$group)


##################################################
## LIKELIHOOD RATIO TEST, CONTROLLING FOR GROUP ##
##################################################

###---------------- FEMALES ONLY ------------------------------------------------------------------

# Create subset of females only
# Have to remove Drd3 & Pvalb because clusters don't contain lacz & reln for both groups
females_only <- subset(allRats, subset = (sex == "female" & cellType != "Drd3-MSN" & cellType != "Pvalb-Int"))
females_only$cellType <- Idents(females_only)

# Create a sample ID column: group_sex_stim
females_only$sample.id <- as.factor(paste(females_only$group, females_only$sex_target, sep = "_"))
# Named vector of sample names
sample_ids <- purrr::set_names(levels(females_only$sample.id))

# Get cell and sample metrics for aggregation
females_only@meta.data$cellType <- as.factor(females_only@meta.data$cellType)
cluster_names <- purrr::set_names(levels(females_only@meta.data$cellType))

# Number of clusters
cluster_total <- length(cluster_names)
# Total number of samples
sample_total <- length(sample_ids)
# Number of cells in each group
table(females_only$sample.id)

###### Count aggregation to sample level

# Make list of all cells: vars = cellType, sample.id
groups <- females_only@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(females_only@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

# Turn into a list, and split list into components for each cluster
# Transform so that in each cluster: rows = genes; cols = cluster_sampleID
raw_counts_list <- split.data.frame(
  count_aggr,
  #  factor(cluster_names) # Changing list names from - to . is optional
  factor(gsub("-", ".", cluster_names)) # Changing list names from - to . is optional
) %>%
  lapply(function(x) {
    magrittr::set_colnames(t(x), gsub("-", ".", rownames(x)))
  })

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
metadata$group     <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",1)))
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
                                design = ~ group + target)
  
  # Set reference to lacz
  dds$target <- relevel(dds$target, ref = "lacz")
  return(dds)
  
# }, x = raw_counts_list, y = 1:length(raw_counts_list), z = names(raw_counts_list))
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = gsub("-", ".", names(raw_counts_list)))

################################
# Remove genes with < 5 counts #
################################

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
dir.create("DESeq2_QC_femLRT")

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
  file_name_pca <- paste0("DESeq2_QC_femLRT/PCA_", z)
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
  
  file_name_heatmap <- paste0("DESeq2_QC_femLRT/sample_to_sample_heatmap_", z)
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

####################
## CALCULATE DEGs ##
####################

# Make output folder
dir.create("DESeq2_DEGResults_femLRT")

dds_all_cells <- mapply(FUN = function(x, z) {
  x <- DESeq(x, test="LRT", reduced = ~ group)
  
  # Checking coeff by resultsNames(dds):
  resultsNames(x)
  
  # Dispersion plot
  file_name_dist <- paste0("DESeq2_DEGResults_femLRT/DispPlot_", z, ".pdf")
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
dir.create("DESeq2_Reln_Collapsed_femLRT")

# Write out the files
mapply(FUN = function(x, z) {
  file_name <- paste0("DESeq2_Reln_Collapsed_femLRT/DESeq2_res_", z, ".csv")
  write.csv(x, file=file_name, row.names = TRUE, quote = FALSE)
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

# Save counts
mapply(FUN = function(x, z) {
  file_name_norm <- paste0("DESeq2_Reln_Collapsed_femLRT/normalized_counts_", z, ".csv")
  normalized_counts <- as.data.frame(counts(x, normalized = TRUE))
  write.csv(normalized_counts, file = file_name_norm, row.names = TRUE, quote = FALSE)
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))


# Make table of number of DEGs
# NOTE: Change nrow to your correct number of names in dds_all_cells (i.e., # of clusters)
DEGs_DF <- as.data.frame(matrix(ncol = 2, nrow=11))
colnames(DEGs_DF) <- c("cellType", "number_of_DEGs")
rownames(DEGs_DF) <- names(dds_all_cells)
DEGs_DF$cellType  <- names(dds_all_cells)

for(i in names(dds_all_cells)){
  x <- as.data.frame(res_all_cells[[i]])
  DEGs_DF[i,"number_of_DEGs"] <- nrow(subset(x, subset=(padj <= 0.05)))
}

write.table(x         = DEGs_DF,
            file      = "DESeq2_Reln_Collapsed_femLRT/number_of_DEGs_Per_cellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

rm(sample_ids, cluster_names, cluster_total, sample_total, groups, count_aggr, raw_counts_list, 
   de_samples, metadata, dds_all_cells, DEGs_DF, genesPctExp, res_all_cells, seuratCounts)

### ---------------- MALES ONLY --------------------------------------------------------------------

# Create subset of males only
males_only <- subset(allRats, subset = (sex == "male"))
males_only$cellType <- Idents(males_only)

# Create a sample ID column: group_sex_stim
males_only$sample.id <- as.factor(paste(males_only$group, males_only$sex_target, sep = "_"))
# Named vector of sample names
sample_ids <- purrr::set_names(levels(males_only$sample.id))

# Get cell and sample metrics for aggregation
males_only@meta.data$cellType <- as.factor(males_only@meta.data$cellType)
cluster_names <- purrr::set_names(levels(males_only@meta.data$cellType))

# Number of clusters
cluster_total <- length(cluster_names)
# Total number of samples
sample_total <- length(sample_ids)
# Number of cells in each group
table(males_only$sample.id)

###### Count aggregation to sample level

# Make list of all cells: vars = cellType, sample.id
groups <- males_only@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(males_only@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 


# Turn into a list, and split list into components for each cluster
# Transform so that in each cluster: rows = genes; cols = cluster_sampleID
raw_counts_list <- split.data.frame(
  count_aggr,
  #  factor(cluster_names) # Changing list names from - to . is optional
  factor(gsub("-", ".", cluster_names)) # Changing list names from - to . is optional
) %>%
  lapply(function(x) {
    magrittr::set_colnames(t(x), gsub("-", ".", rownames(x)))
  })


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
metadata$group     <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",1)))
metadata$sex       <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",2)))
metadata$target    <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",3)))

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
                                design = ~ group + target)
  
  # Set reference to lacz
  dds$target <- relevel(dds$target, ref = "lacz")
  return(dds)
  
  # }, x = raw_counts_list, y = 1:length(raw_counts_list), z = names(raw_counts_list))
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = gsub("-", ".", names(raw_counts_list)))

################################
# Remove genes with < 5 counts #
################################

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
dir.create("DESeq2_QC_maleLRT")

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
  file_name_pca <- paste0("DESeq2_QC_maleLRT/PCA_", z)
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
  
  file_name_heatmap <- paste0("DESeq2_QC_maleLRT/sample_to_sample_heatmap_", z)
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

####################
## CALCULATE DEGs ##
####################

# Make output folder
dir.create("DESeq2_DEGResults_maleLRT")

dds_all_cells <- mapply(FUN = function(x, z) {
  x <- DESeq(x, test="LRT", reduced = ~ group)
  
  # Checking coeff by resultsNames(dds):
  resultsNames(x)
  
  # Dispersion plot
  file_name_dist <- paste0("DESeq2_DEGResults_maleLRT/DispPlot_", z, ".pdf")
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
dir.create("DESeq2_Reln_Collapsed_maleLRT")

# Write out the files
mapply(FUN = function(x, z) {
  file_name <- paste0("DESeq2_Reln_Collapsed_maleLRT/DESeq2_res_", z, ".csv")
  write.csv(x, file=file_name, row.names = TRUE, quote = FALSE)
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

# Save counts
mapply(FUN = function(x, z) {
  file_name_norm <- paste0("DESeq2_Reln_Collapsed_maleLRT/normalized_counts_", z, ".csv")
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
            file      = "DESeq2_Reln_Collapsed_maleLRT/number_of_DEGs_Per_cellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

rm(sample_ids, cluster_names, cluster_total, sample_total, groups, count_aggr, raw_counts_list, 
   de_samples, metadata, dds_all_cells, DEGs_DF, genesPctExp, res_all_cells, seuratCounts)


###################################
## WALD TEST, NO REDUCING FACTOR ##
###################################

###---------------- FEMALES ONLY ------------------------------------------------------------------

# Create subset of females only
# Have to remove Drd3 & Pvalb because clusters don't contain lacz & reln for both groups
females_only <- subset(allRats, subset = (sex == "female" & cellType != "Drd3-MSN" & cellType != "Pvalb-Int"))
females_only$cellType <- Idents(females_only)

# Create a sample ID column: group_sex_stim
females_only$sample.id <- as.factor(paste(females_only$group, females_only$sex_target, sep = "_"))
# Named vector of sample names
sample_ids <- purrr::set_names(levels(females_only$sample.id))

# Get cell and sample metrics for aggregation
females_only@meta.data$cellType <- as.factor(females_only@meta.data$cellType)
cluster_names <- purrr::set_names(levels(females_only@meta.data$cellType))

# Number of clusters
cluster_total <- length(cluster_names)
# Total number of samples
sample_total <- length(sample_ids)
# Number of cells in each group
table(females_only$sample.id)

###### Count aggregation to sample level

# Make list of all cells: vars = cellType, sample.id
groups <- females_only@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(females_only@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

# Turn into a list, and split list into components for each cluster
# Transform so that in each cluster: rows = genes; cols = cluster_sampleID
raw_counts_list <- split.data.frame(
  count_aggr,
  #  factor(cluster_names) # Changing list names from - to . is optional
  factor(gsub("-", ".", cluster_names)) # Changing list names from - to . is optional
) %>%
  lapply(function(x) {
    magrittr::set_colnames(t(x), gsub("-", ".", rownames(x)))
  })

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
metadata$group     <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",1)))
metadata$sex       <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",2)))
metadata$target    <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",3)))

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
                                design = ~ target)
  
  # Set reference to lacz
  dds$target <- relevel(dds$target, ref = "lacz")
  return(dds)
  
# }, x = raw_counts_list, y = 1:length(raw_counts_list), z = names(raw_counts_list))
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = gsub("-", ".", names(raw_counts_list)))

################################
# Remove genes with < 5 counts #
################################

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
dir.create("DESeq2_QC_femWald")

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
  file_name_pca <- paste0("DESeq2_QC_femWald/PCA_", z)
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
  
  file_name_heatmap <- paste0("DESeq2_QC_femWald/sample_to_sample_heatmap_", z)
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

####################
## CALCULATE DEGs ##
####################

# Make output folder
dir.create("DESeq2_DEGResults_femWald")

dds_all_cells <- mapply(FUN = function(x, z) {
  # LRT
#   x <- DESeq(x, test="LRT", reduced = ~ group)
  # Wald test, no reducing factor
  x <- DESeq(x)
  
  # Checking coeff by resultsNames(dds):
  resultsNames(x)
  
  # Dispersion plot
  file_name_dist <- paste0("DESeq2_DEGResults_femWald/DispPlot_", z, ".pdf")
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
dir.create("DESeq2_Reln_Collapsed_femWald")

# Write out the files
mapply(FUN = function(x, z) {
  file_name <- paste0("DESeq2_Reln_Collapsed_femWald/DESeq2_res_", z, ".csv")
  write.csv(x, file=file_name, row.names = TRUE, quote = FALSE)
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

# Save counts
mapply(FUN = function(x, z) {
  file_name_norm <- paste0("DESeq2_Reln_Collapsed_femWald/normalized_counts_", z, ".csv")
  normalized_counts <- as.data.frame(counts(x, normalized = TRUE))
  write.csv(normalized_counts, file = file_name_norm, row.names = TRUE, quote = FALSE)
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))


# Make table of number of DEGs
# NOTE: Change nrow to your correct number of names in dds_all_cells (i.e., # of clusters)
DEGs_DF <- as.data.frame(matrix(ncol = 2, nrow=11))
colnames(DEGs_DF) <- c("cellType", "number_of_DEGs")
rownames(DEGs_DF) <- names(dds_all_cells)
DEGs_DF$cellType  <- names(dds_all_cells)

for(i in names(dds_all_cells)){
  x <- as.data.frame(res_all_cells[[i]])
  DEGs_DF[i,"number_of_DEGs"] <- nrow(subset(x, subset=(padj <= 0.05)))
}

write.table(x         = DEGs_DF,
            file      = "DESeq2_Reln_Collapsed_femWald/number_of_DEGs_Per_cellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

rm(sample_ids, cluster_names, cluster_total, sample_total, groups, count_aggr, raw_counts_list, 
   de_samples, metadata, dds_all_cells, DEGs_DF, genesPctExp, res_all_cells, seuratCounts)


###---------------- MALES ONLY --------------------------------------------------------------------

# Create subset of males only
males_only <- subset(allRats, subset = (sex == "male"))
males_only$cellType <- Idents(males_only)

# Create a sample ID column: group_sex_stim
males_only$sample.id <- as.factor(paste(males_only$group, males_only$sex_target, sep = "_"))
# Named vector of sample names
sample_ids <- purrr::set_names(levels(males_only$sample.id))

# Get cell and sample metrics for aggregation
males_only@meta.data$cellType <- as.factor(males_only@meta.data$cellType)
cluster_names <- purrr::set_names(levels(males_only@meta.data$cellType))

# Number of clusters
cluster_total <- length(cluster_names)
# Total number of samples
sample_total <- length(sample_ids)
# Number of cells in each group
table(males_only$sample.id)

###### Count aggregation to sample level

# Make list of all cells: vars = cellType, sample.id
groups <- males_only@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(males_only@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

# Turn into a list, and split list into components for each cluster
# Transform so that in each cluster: rows = genes; cols = cluster_sampleID
raw_counts_list <- split.data.frame(
  count_aggr,
  #  factor(cluster_names) # Changing list names from - to . is optional
  factor(gsub("-", ".", cluster_names)) # Changing list names from - to . is optional
) %>%
  lapply(function(x) {
    magrittr::set_colnames(t(x), gsub("-", ".", rownames(x)))
  })

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
metadata$group     <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",1)))
metadata$sex       <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",2)))
metadata$target    <- as.factor(as.character(lapply(strsplit(metadata$sex.target,"_"),"[",3)))

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
                                design = ~ target)
  
  # Set reference to lacz
  dds$target <- relevel(dds$target, ref = "lacz")
  return(dds)
  
  # }, x = raw_counts_list, y = 1:length(raw_counts_list), z = names(raw_counts_list))
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = gsub("-", ".", names(raw_counts_list)))

################################
# Remove genes with < 5 counts #
################################

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
dir.create("DESeq2_QC_maleWald")

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
  file_name_pca <- paste0("DESeq2_QC_maleWald/PCA_", z)
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
  
  file_name_heatmap <- paste0("DESeq2_QC_maleWald/sample_to_sample_heatmap_", z)
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))

####################
## CALCULATE DEGs ##
####################

# Make output folder
dir.create("DESeq2_DEGResults_maleWald")

dds_all_cells <- mapply(FUN = function(x, z) {
  # LRT
#   x <- DESeq(x, test="LRT", reduced = ~ group)
  # Wald test, no reducing factor
  x <- DESeq(x)
  
  # Checking coeff by resultsNames(dds):
  resultsNames(x)
  
  # Dispersion plot
  file_name_dist <- paste0("DESeq2_DEGResults_maleWald/DispPlot_", z, ".pdf")
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
dir.create("DESeq2_Reln_Collapsed_maleWald")

# Write out the files
mapply(FUN = function(x, z) {
  file_name <- paste0("DESeq2_Reln_Collapsed_maleWald/DESeq2_res_", z, ".csv")
  write.csv(x, file=file_name, row.names = TRUE, quote = FALSE)
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

# Save counts
mapply(FUN = function(x, z) {
  file_name_norm <- paste0("DESeq2_Reln_Collapsed_maleWald/normalized_counts_", z, ".csv")
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
            file      = "DESeq2_Reln_Collapsed_maleWald/number_of_DEGs_Per_cellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

rm(sample_ids, cluster_names, cluster_total, sample_total, groups, count_aggr, raw_counts_list, 
   de_samples, metadata, dds_all_cells, DEGs_DF, genesPctExp, res_all_cells, seuratCounts)
