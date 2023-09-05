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
allRats <- readRDS("allRats_souped_noDoub.rds")

# Run if necessary to make sure cellType matches renamed idents
allRats$cellType <- Idents(allRats)

## --------------------BARPLOTS FOR TREATMENT GROUP & SEX BY CELL TYPE----------------------------------------
# Modified from code written by Guy Twa!

# Treatment Group
T_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$target)
T_count_mtx <- as.data.frame.matrix(T_count_table)
T_count_mtx$cluster <- rownames(T_count_mtx)
T_melt_mtx <- melt(T_count_mtx)
T_melt_mtx$cluster <- factor(T_melt_mtx$cluster, levels = levels(allRats@meta.data$cellType))

# Will be the same for all integration features; will only calculate once
cluster_size <- aggregate(value ~ cluster, data = T_melt_mtx, FUN = sum)
cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(allRats@meta.data$cellType))

colnames(T_melt_mtx)[2] <- "Treatment"

# Sex
S_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$sex)
S_count_mtx <- as.data.frame.matrix(S_count_table)
S_count_mtx$cluster <- rownames(S_count_mtx)
S_melt_mtx <- melt(S_count_mtx)
S_melt_mtx$cluster <- factor(S_melt_mtx$cluster)
S_melt_mtx$cluster <- factor(S_melt_mtx$cluster,levels = levels(allRats@meta.data$cellType))
colnames(S_melt_mtx)[2] <- "Sex"

# GEM well
G_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$GEM)
G_count_mtx   <- as.data.frame.matrix(G_count_table)
G_count_mtx$cluster <- rownames(G_count_mtx)
G_melt_mtx    <- melt(G_count_mtx)
G_melt_mtx$cluster <- factor(G_melt_mtx$cluster)
G_melt_mtx$cluster <- factor(G_melt_mtx$cluster,levels = levels(allRats@meta.data$cellType))
colnames(G_melt_mtx)[2] <- "GEM"

# FACS group
F_count_table <- table(allRats@meta.data$cellType, allRats@meta.data$group)
F_count_mtx   <- as.data.frame.matrix(F_count_table)
F_count_mtx$cluster <- rownames(F_count_mtx)
F_melt_mtx    <- melt(F_count_mtx)
F_melt_mtx$cluster <- factor(F_melt_mtx$cluster)
F_melt_mtx$cluster <- factor(F_melt_mtx$cluster,levels = levels(allRats@meta.data$cellType))
colnames(F_melt_mtx)[2] <- "FACS_Group"

# Plot treatment group
# This is left-most plot, so includes Y-axis labels
p1 <- ggplot(T_melt_mtx, aes(x = cluster, y = value, fill = Treatment)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("#E6AB02", "#E7298A")) +
  ylab("") + xlab("Cluster") + theme(legend.position = "bottom", axis.text.y = element_text(size = 16), legend.text = element_text(size = 16), 
                                     plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
  ggtitle("Treatment Group") + guides(fill = guide_legend(reverse = TRUE))

# Plot sex
p2 <- ggplot(S_melt_mtx, aes(x = cluster, y = value, fill = Sex)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("#7570B3", "#1B9E77")) +
  ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                   legend.text = element_text(size = 16), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
  ggtitle("Sex") + guides(fill = guide_legend(reverse = TRUE))

# Plot GEM well
p3 <- ggplot(G_melt_mtx, aes(x = cluster, y = value, fill = GEM)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") + theme_bw() + coord_flip() + scale_fill_brewer(palette = "Dark2") +
  ylab("") + theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                   legend.text = element_text(size = 16), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
  ggtitle("GEM")

# Show 3 plots together
p1 + p2 + p3

# Plot FACS group (individually, so also has Y-axis labels)
plot_group <- ggplot(F_melt_mtx, aes(x = cluster, y = value, fill = FACS_Group)) + 
  geom_bar(position = "fill", stat = "identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c("#E6AB02", "#E7298A")) +
  ylab("") + xlab("Cluster") + theme(legend.position = "bottom", axis.text.y = element_text(size = 16), legend.text = element_text(size = 16), 
                                     plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_blank()) +
  ggtitle("FACS Group") + guides(fill = guide_legend(reverse = TRUE))


## --------------------DOTPLOT: GENE EXPRESSION BY CELL TYPE--------------------------------------------------

# Set default assay to original data
DefaultAssay(allRats) <- "RNA"

# Set genes to be plotted (these are the marker genes)
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", "Chat", "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", "Rgs5", "Mbp", "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", "Gad1", "Syt1")

# Create dotplot
DotPlot(allRats, features = features, cols = "Spectral", dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))

DotPlot(Drd1_lacz, features = "Reln", cols = "Spectral", dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))

## --------------------HEATMAP: GENE EXPRESSION BY CELL TYPE--------------------------------------------------

# Set default assay
DefaultAssay(allRats) <- "integrated"

# Create heatmap
DoHeatmap(subset(allRats, downsample = 100), features = features, size = 3)

## --------------------RELN EXPRESSION in DRD1 BY TREATMENT--------------------------------------------------

# Subset data
Drd1 <- subset(allRats, idents = c("Drd1-MSN"))
Drd1_Reln <- subset(Drd1, subset = (target == "reln"))
Drd1_lacz <- subset(Drd1, subset = (target == "lacz"))
Drd1_F <- subset(Drd1, subset = (sex == "female"))
Drd1_M <- subset(Drd1, subset = (sex == "male"))

# Correlation between Reln expression and mCherry expression
FeatureScatter(Drd1, "mCherry", "Reln") + NoLegend()
FeatureScatter(Drd1, "dCas9-KRAB-MeCP2", "Reln") + NoLegend()
FeatureScatter(Drd1_lacz, "mCherry", "Reln") + NoLegend()
FeatureScatter(Drd1_lacz, "dCas9-KRAB-MeCP2", "Reln") + NoLegend()
FeatureScatter(Drd1_Reln, "mCherry", "Reln") + NoLegend()
FeatureScatter(Drd1_Reln, "dCas9-KRAB-MeCP2", "Reln") + NoLegend()

# Violin plot of RELN, split by target (lacZ vs. Reln knockdown)
vln_drd1 <- VlnPlot(Drd1, features = "Reln", split.by = "target") + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Violin plot of Hs6st3 for each sex, split by target
vln_F <- VlnPlot(Drd1_F, features = "Hs6st3", split.by = "target") + ggtitle("Hs6st3 (Females)") + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
vln_M <- VlnPlot(Drd1_M, features = "Hs6st3", split.by = "target") + ggtitle("Hs6st3 (Males)") + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
vln_F + vln_M


####### UP-SET PLOT FOR DEGs

library(UpSetR)
library(ComplexUpset)

d1 <- read.csv("degs_drd1.csv", header = FALSE)
d2_1 <- read.csv("degs_drd2-1.csv", header = FALSE)
d2_2 <- read.csv("degs_drd2-2.csv", header = FALSE)
d3 <- read.csv("degs_drd3.csv", header = FALSE)
gaba <- read.csv("degs_gaba.csv", header = FALSE)
astro <- read.csv("degs_astro.csv", header = FALSE)
micro <- read.csv("degs_micro.csv", header = FALSE)
olig <- read.csv("degs_olig.csv", header = FALSE)
poly <- read.csv("degs_poly.csv", header = FALSE)

# Create list
degs_list <- list(Drd1 = d1$V1, Drd2_1 = d2_1$V1, Drd2_2 = d2_2$V1, Drd3 = d3$V1, GABA = gaba$V1, 
                  Astro = astro$V1, Micro = micro$V1, Olig = olig$V1, Poly = poly$V1)

# Create dataframe from list
# Code from UpSetR package
data <- fromList(degs_list)

# ComplexUpset code
# CAN sort intersection bars and sets
upset(data, colnames(data), intersections = list(
  c("Drd1"), c("Drd2_1"), c("Drd2_2"), c("Drd3"), c("GABA"), c("Drd1", "Drd2_1"), 
  c("Drd1", "Drd3"), c("Drd1", "Drd2_1", "Drd3"), c("Drd2_1", "Micro"), 
  c("Drd3", "Astro", "Olig"), c("Astro"), c("Micro"), c("Olig"), c("Poly")), 
  sort_intersections = FALSE, sort_sets = FALSE, 
  queries = list(
    upset_query(intersect = c("Drd1"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd2_1"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd2_2"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd3"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("GABA"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd1", "Drd2_1"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd1", "Drd3"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd1", "Drd2_1", "Drd3"), color = "red", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd2_1", "Micro"), color = "orange", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Drd3", "Astro", "Olig"), color = "orange", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Poly"), color = "blue", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Olig"), color = "blue", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Micro"), color = "blue", only_components = c("intersections_matrix")), 
    upset_query(intersect = c("Astro"), color = "blue", only_components = c("intersections_matrix"))))

