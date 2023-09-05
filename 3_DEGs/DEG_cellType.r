####################################
##                                ##
## DIFFERENTIALLY EXPRESSED GENES ##
##       BETWEEN CELL TYPES       ##
##                                ##
####################################

# Finds genes that are differentially expressed between cell types
# regardless of sample (GEM well), target (lacZ vs. Reln), or sex. (IS THIS CORRECT???)

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
allRats <- RenameIdents(allRats, "Mural?" = "Mural")


# Make sure cellType is correct after integrating & DoubletFinder
# DESeq2 code below will error if cellType still includes NAs
allRats$cellType <- Idents(allRats)

# Create a sample ID column: GEM_sex_target (as factor)
# In this experiment, each individual animal had own GEM well, so "sample" = animal
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
# Sum of counts for each sample within each cell type

# Make list of all cells: vars = cellType, sample.id
groups <- allRats@meta.data[, c("cellType", "sample.id")]

# Aggregate across cluster-sample groups using raw counts ("counts" slot)
# Matrix of raw counts: row names = cluster_sampleID; col names = genes
count_aggr <- Matrix.utils::aggregate.Matrix(t(allRats@assays$RNA@counts), 
                                             groupings = groups, fun = "sum")

# Transpose count_aggr
# row names = genes; col names = cluster_sampleID
count_aggr_t2 <- t(count_aggr2)

# Optional change, replacing - with .
# Colnames will now be this format: Drd1.MSN.1_3_fem_lacz
colnames(count_aggr_t2) <- gsub(x = colnames(count_aggr_t2), pattern = "-", replacement = ".")

# Create metadata dataframe
metadata <- data.frame(cluster_id     = sub("_.*","", colnames(count_aggr_t)), #takes first part of colname, before _ (e.g., Drd1.MSN.1)
                       cluster_sample = colnames(count_aggr_t),
                       sample_id      = sub("^[^_]*_","", colnames(count_aggr_t))) #takes second part of colname, after first _ (e.g., 3_fem_lacz)
metadata$GEM       <- as.factor(as.character(lapply(strsplit(metadata$sample_id,"_"),"[",1)))
metadata$sex       <- as.factor(as.character(lapply(strsplit(metadata$sample_id,"_"),"[",2)))
metadata$target    <- as.factor(as.character(lapply(strsplit(metadata$sample_id,"_"),"[",3)))
rownames(metadata) <- metadata$cluster_sample

####################
## CALCULATE DEGs ##
####################

counts_mat <- allRats@assays$RNA@counts
Idents(allRats) <- gsub(x = Idents(allRats), pattern = "-", replacement = ".")

# Make folders for output
dir.create("DESeq2_cellTypes")
dir.create("DESeq2_cellTypes/PCA")

for(i in unique(metadata$cluster_id)){

  ############ Make the DESEq object ############
  metadata$cellType <- ifelse(metadata$cluster_id == i,
                              i,
                              "Other")
  # Can't include both GEM & target in design because they are "identical"
  dds <- DESeqDataSetFromMatrix(count_aggr_t, 
                                colData = metadata, 
                                design = ~ GEM + cellType)
  # Set reference level = control group
  dds$cellType <- relevel(dds$cellType, ref = "Other")
  # Filter
  keep <- rowMeans(counts(dds)) > 5
  dds <- dds[keep,]
  print(paste("dds made for",i))

  ############ Quality control ############
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  # PCA for this cell type vs. all others
  data <- plotPCA(vsd, intgroup = c("cellType"), returnData = TRUE)
  data$cluster_id <- as.character(lapply(strsplit(x = data$name, split = "_"),"[",1))
  percentVar <- round(100 * attr(data, "percentVar"))
  pca.plot.cellType <- ggplot(data, aes(PC1, PC2, color = cellType)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)

  # PCA for cluster_id
  pca.plot.cluster_id <- ggplot(data, aes(PC1, PC2, color = cluster_id)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)

  ggsave(filename = paste0("DESeq2_cellTypes/PCA/",i,".pdf"),
    plot   = cowplot::plot_grid(plotlist = list(pca.plot.cellType, pca.plot.cluster_id), ncol = 2),
    width  = 12,
    height = 7)
  
  print(paste("PCA plots made for", i))

  ############ Get results ############
  dds <- DESeq(dds, test = "LRT", reduced = ~ GEM)
  res <- as.data.frame(results(dds))
  res$geneName <- rownames(res)
  # All cells of this type
  cells.1    <- WhichCells(allRats, idents = i)
  # All other cells
  cells.2    <- setdiff(x = row.names(allRats@meta.data), y = cells.1)
  # What is the point of this?
  res$Pct_Expressing <- NA
  # Percent of this cell type expressing the gene
  res$Pct_CellType_Expressing <- (rowSums(x = counts_mat[res$geneName, cells.1, drop = FALSE] > 0) / length(cells.1))*100
  # Percent of all other cell types expressing the gene
  res$Pct_Other_Expressing <- (rowSums(x = counts_mat[res$geneName, cells.2, drop = FALSE] > 0) / length(cells.2))*100
  write.table(file      = paste0("DESeq2_cellTypes/",i,".txt"),
              x         = res,
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE)
  print(paste("DEG lists written for", i))
  rm(dds, keep, vsd, data, percentVar, pca.plot.cellType, pca.plot.cluster_id, res)
}

