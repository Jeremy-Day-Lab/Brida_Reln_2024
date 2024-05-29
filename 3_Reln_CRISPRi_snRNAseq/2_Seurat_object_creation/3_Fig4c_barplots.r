###############
## BAR PLOTS ##
###############

# Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Load Seurat object
allRats <- readRDS("allRats_souped_noDoub.rds")

# Run if necessary to make sure cellType matches renamed idents
allRats$cellType <- Idents(allRats)

## --------------------BARPLOTS FOR TREATMENT GROUP & SEX BY CELL TYPE----------------------------------------
# Modified from code written by Guy Twa

# Treatment Group by cluster
T_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$target)
T_count_mtx <- as.data.frame.matrix(T_count_table)
T_count_mtx$cluster <- rownames(T_count_mtx)
T_melt_mtx <- melt(T_count_mtx)
T_melt_mtx$cluster <- factor(T_melt_mtx$cluster, levels = levels(allRats@meta.data$cellType))

# Will be the same for all integration features; will only calculate once
cluster_size <- aggregate(value ~ cluster, data = T_melt_mtx, FUN = sum)
cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(allRats@meta.data$cellType))

colnames(T_melt_mtx)[2] <- "Treatment"

# Proportions of clusters in each GEM well (sample)
G_count_table2 <- table(allRats@meta.data$GEM, allRats@meta.data$cellType)
G_count_mtx2   <- as.data.frame.matrix(G_count_table2)
G_count_mtx2$GEM <- rownames(G_count_mtx2)
G_melt_mtx2    <- melt(G_count_mtx2)
colnames(G_melt_mtx2)[2] <- "Cluster"

# Plot clusters by GEM well (sample)
p1 <- ggplot(G_melt_mtx2, aes(x = GEM, y = value, fill = Cluster)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + theme_bw() + coord_flip() +
  ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                   legend.text = element_text(size = 16), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
  ggtitle("Cluster by GEM")


# ## Additional plots (not used)
# 
# # Sex by cluster
# S_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$sex)
# S_count_mtx <- as.data.frame.matrix(S_count_table)
# S_count_mtx$cluster <- rownames(S_count_mtx)
# S_melt_mtx <- melt(S_count_mtx)
# S_melt_mtx$cluster <- factor(S_melt_mtx$cluster)
# S_melt_mtx$cluster <- factor(S_melt_mtx$cluster,levels = levels(allRats@meta.data$cellType))
# colnames(S_melt_mtx)[2] <- "Sex"
# 
# # GEM well by cluster
# G_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$GEM)
# G_count_mtx   <- as.data.frame.matrix(G_count_table)
# G_count_mtx$cluster <- rownames(G_count_mtx)
# G_melt_mtx    <- melt(G_count_mtx)
# G_melt_mtx$cluster <- factor(G_melt_mtx$cluster)
# G_melt_mtx$cluster <- factor(G_melt_mtx$cluster,levels = levels(allRats@meta.data$cellType))
# colnames(G_melt_mtx)[2] <- "GEM"
# 
# # Plot treatment group
# # This is left-most plot, so includes Y-axis labels
# p1 <- ggplot(T_melt_mtx, aes(x = cluster, y = value, fill = Treatment)) + 
#   geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("#E6AB02", "#E7298A")) +
#   ylab("") + xlab("Cluster") + theme(legend.position = "bottom", axis.text.y = element_text(size = 16), legend.text = element_text(size = 16), 
#                                      plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
#   ggtitle("Treatment Group") + guides(fill = guide_legend(reverse = TRUE))
# 
# # Plot sex
# p2 <- ggplot(S_melt_mtx, aes(x = cluster, y = value, fill = Sex)) + 
#   geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("#7570B3", "#1B9E77")) +
#   ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
#                    legend.text = element_text(size = 16), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
#   ggtitle("Sex") + guides(fill = guide_legend(reverse = TRUE))
# 
# # Plot GEM well
# p3 <- ggplot(G_melt_mtx, aes(x = cluster, y = value, fill = GEM)) + 
#   geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + theme_bw() + coord_flip() + scale_fill_brewer(palette = "Dark2") +
#   ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
#                    legend.text = element_text(size = 16), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
#   ggtitle("GEM")
# 
# # Show 3 plots together
# p1 + p2 + p3

