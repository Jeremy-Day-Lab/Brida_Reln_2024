###############################
##                           ##
## QC & OTHER VISUALIZATIONS ##
##                           ##
###############################

# Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Load Seurat object
allGroups <- readRDS("allGroups_souped_noDoub.rds")

# Run if necessary to make sure cellType matches renamed idents
allGroups$cellType <- Idents(allGroups)

## --------------------BARPLOTS FOR TREATMENT GROUP & SEX BY CELL TYPE----------------------------------------

# Treatment Group
T_count_table <- table(allGroups@meta.data$cellType, allGroupsRepeated@meta.data$target)
T_count_mtx <- as.data.frame.matrix(T_count_table)
T_count_mtx$cluster <- rownames(T_count_mtx)
T_melt_mtx <- melt(T_count_mtx)
T_melt_mtx$cluster <- factor(T_melt_mtx$cluster, levels = levels(allGroups@meta.data$cellType))

# Will be the same for all integration features; will only calculate once
cluster_size <- aggregate(value ~ cluster, data = T_melt_mtx, FUN = sum)
cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(allGroups@meta.data$cellType))

colnames(T_melt_mtx)[2] <- "Treatment"

# Sex
S_count_table <- table(allGroups@meta.data$cellType, allGroups@meta.data$sex)
S_count_mtx <- as.data.frame.matrix(S_count_table)
S_count_mtx$cluster <- rownames(S_count_mtx)
S_melt_mtx <- melt(S_count_mtx)
S_melt_mtx$cluster <- factor(S_melt_mtx$cluster)
S_melt_mtx$cluster <- factor(S_melt_mtx$cluster,levels = levels(allGroups@meta.data$cellType))
colnames(S_melt_mtx)[2] <- "Sex"

p1 <- ggplot(T_melt_mtx, aes(x = cluster, y = value, fill = Treatment)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("salmon", "bisque")) +
  ylab("") + xlab("Cluster") + theme(legend.position = "bottom") + ggtitle("Treatment Group")

p2 <- ggplot(S_melt_mtx, aes(x = cluster, y = value, fill = Sex)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("plum", "#1B9E77")) +
  ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle("Sex")

p1 + p2

## --------------------DOTPLOT: GENE EXPRESSION BY CELL TYPE--------------------------------------------------

# Set default assay to original data
DefaultAssay(allGroups) <- "RNA"

# Set genes to be plotted (these are the marker genes)
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", "Chat", "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", "Rgs5", "Mbp", "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", "Gad1", "Syt1")

# Create dotplot
DotPlot(allGroups, features = features, cols = "Spectral", dot.scale = 8) + RotatedAxis()

## --------------------HEATMAP: GENE EXPRESSION BY CELL TYPE--------------------------------------------------

# Set default assay
DefaultAssay(allGroupsRepeated) <- "integrated"

# Create heatmap
DoHeatmap(subset(allGroupsRepeated, downsample = 100), features = features, size = 3)






