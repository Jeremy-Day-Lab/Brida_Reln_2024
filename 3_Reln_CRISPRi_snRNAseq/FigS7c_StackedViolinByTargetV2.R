#The point of this script is to make a stacked violin plot that compares the
#expression of marker genes between an experimental & control group

#Set working directory and seed
setwd("/Volumes/Day-Lab$/10X_SingleCell/2023_SC_0016/analysis/DataVis_Graphs/StackedViolin/")
set.seed(1234)

######Load libraries
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(cowplot)

#Read in object
SC0016 <- readRDS(file = "/Volumes/Day-Lab$/10X_SingleCell/2023_SC_0016/analysis/Seurat/R objects/allRats_souped_noDoub.rds")
SC0016$cellType <- Idents(SC0016)
DefaultAssay(SC0016) <- "RNA"

#Create metadata df to merge to add in gRNA info
Idents(SC0016) <- SC0016$target
SC0016_MetaTarget <- data.frame(SC0016@active.ident)
colnames(SC0016_MetaTarget)[1] = "Target"
Idents(SC0016) <- SC0016$cellType
SC0016_MetaType <- data.frame(SC0016@active.ident)
colnames(SC0016_MetaType)[1] = "CellType"
SC0016_MetaData <- cbind(cell = rownames(SC0016_MetaType), SC0016_MetaType, SC0016_MetaTarget)
lacZ_Meta <- subset(SC0016_MetaData, subset=SC0016_MetaData$Target == "lacz")
Reln_Meta <- subset(SC0016_MetaData, subset=SC0016_MetaData$Target == "reln")


# # Set genes to be plotted (these are the marker genes)
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", "Chat",
              "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", "Rgs5", "Mbp",
              "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", "Gad1", "Syt1")

#Create marker gene dataframe
Idents(SC0016) <- SC0016$target

#lacZ
SC0016_lacZ_df <- as.data.frame(t(as.data.frame(GetAssayData(object = SC0016,
                                                                     slot = "data",
                                                                     assay = "RNA")[features,
                                                                                    WhichCells(object = subset(x = subset(SC0016,idents = "lacz")))])))
SC0016_lacZ_df <- cbind(cell = rownames(SC0016_lacZ_df), SC0016_lacZ_df)
SC0016_lacZ_df <- left_join(SC0016_lacZ_df, lacZ_Meta)

#Reln
SC0016_Reln_df <- as.data.frame(t(as.data.frame(GetAssayData(object = SC0016,
                                                          slot = "data",
                                                          assay = "RNA")[features,
                                                                         WhichCells(object = subset(x = subset(SC0016,idents = "reln")))])))
SC0016_Reln_df <- cbind(cell = rownames(SC0016_Reln_df), SC0016_Reln_df)
SC0016_Reln_df <- left_join(SC0016_Reln_df, Reln_Meta)

SC0016_FullDF <- rbind(SC0016_lacZ_df, SC0016_Reln_df)

#Gotta Reshape this
SC0016_FullDF2 <- melt(SC0016_FullDF)

#####ViolinPlots#####
# ggplot(SC0016_FullDF, aes(CellType, Drd1, fill=Target)) +
#   geom_violin()

# Identity on x-axis
a <- ggplot(SC0016_FullDF2, aes(factor(CellType), value, fill = Target)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(variable), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Identity on x-axis") + xlab("Identity") + ylab("Expression Level")
a


b <- ggplot(SC0016_FullDF2, aes(value, factor(CellType), fill = Target)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(variable), scales = "free")  +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  ggtitle("Identity on y-axis") + xlab("Expression Level") + ylab("Identity")
b

color <- c("#999999", "#FF1549")

i <- ggplot(SC0016_FullDF2, aes(factor(variable), value, fill = Target)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(CellType), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Feature on x-axis") + xlab("Feature") + ylab("Expression Level")


i

