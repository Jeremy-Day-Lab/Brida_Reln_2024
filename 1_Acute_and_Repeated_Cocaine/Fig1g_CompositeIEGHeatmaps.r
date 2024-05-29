set.seed(1234)

######Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)


NAc_Combo <- readRDS(file = "/Volumes/Day-Lab$/KLB/Manuscripts/2023_Reelin/Acute&Repeated/NAc_Combo_Integrated.RDS")


# Set default assay
DefaultAssay(NAc_Combo) <- "RNA"

#Subset Drd1 MSNs
D1 <- subset(NAc_Combo, idents = c("Drd1-MSN-1", "Drd1-MSN-2"))

#Add metadata column to describe Reln status
poscells <- WhichCells(D1, expression = Reln > 0)
D1@meta.data$RelnStatus<- ifelse(colnames(D1) %in% poscells, "Pos", "Neg")

#Now make a column that describes treatment & reln status
D1$RelnTx<- ifelse(D1@meta.data$RelnStatus == "Neg" & D1@meta.data$Stim == "Saline", "SalNeg",
                   ifelse(D1@meta.data$RelnStatus == "Neg" & D1@meta.data$Stim == "Cocaine", "CokeNeg",
                          ifelse(D1@meta.data$RelnStatus == "Pos" & D1@meta.data$Stim == "Saline", "SalPos",
                                 "CokePos")))
                          


#Set this new metadata item as identity
Idents(D1) <- D1$RelnTx

#Make new Seurat object where all cells within an identity are averaged
D1_avg <- AverageExpression(D1,return.seurat = TRUE)


#Make IEG list
IEGs <- c("Arc", "Btg2", "Egr1", "Egr2", "Egr4", "Fos", "Fosb", "Fosl2", "Homer1",
          "Junb", "Nr4a1", "Nr4a3", "Scg2", "Tiparp", "Vgf")


#Make figure with heatmap of all key genes
DoHeatmap(object = D1_avg,features = IEGs,disp.min=0,disp.max=2,group.bar=TRUE,slot="data",label = TRUE,size = 3,draw.lines = FALSE, raster = FALSE) + 
  theme(axis.text.y = element_text(size = 7)) + 
  scale_fill_gradientn(colors = c("white", "#FF1549"))


#Normalize to row average
IEG_SalNeg <- as.data.frame(GetAssayData(object = D1_avg,slot = "data",assay = "RNA")[IEGs,WhichCells(object = subset(D1_avg,idents = "SalNeg"))])
IEG_SalPos <- as.data.frame(GetAssayData(object = D1_avg,slot = "data",assay = "RNA")[IEGs,WhichCells(object = subset(D1_avg,idents = "SalPos"))])
IEG_CocNeg <- as.data.frame(GetAssayData(object = D1_avg,slot = "data",assay = "RNA")[IEGs,WhichCells(object = subset(D1_avg,idents = "CokeNeg"))])
IEG_CocPos <- as.data.frame(GetAssayData(object = D1_avg,slot = "data",assay = "RNA")[IEGs,WhichCells(object = subset(D1_avg,idents = "CokePos"))])

IEG_df <- cbind(IEG_SalPos, IEG_SalNeg, IEG_CocNeg, IEG_CocPos)
IEG_df$Means <-apply(IEG_df,1,mean)
colnames(IEG_df) <- c("SalPos", "SalNeg", "CokeNeg", "CokePos", "Mean")
IEG_norm <- data.frame(SalPos = IEG_df$SalPos/IEG_df$Mean,
                          SalNeg = IEG_df$SalNeg/IEG_df$Mean,
                          CokeNeg = IEG_df$CokeNeg/IEG_df$Mean,
                          CokePos = IEG_df$CokePos/IEG_df$Mean)


#Merge with metadata
IEG_normT <- t(IEG_norm)
colnames(IEG_normT) <- c("Arc", "Btg2", "Egr1", "Egr2", "Egr4", "Fos", "Fosb", "Fosl2", "Homer1",
                         "Junb", "Nr4a1", "Nr4a3", "Scg2", "Tiparp", "Vgf")
D1_avg@meta.data <- cbind(D1_avg@meta.data, IEG_normT)

options(Seurat.object.assay.version = "v5")
norm <- CreateSeuratObject(counts = pbmc.counts)
class(pbmc[["RNA"]])

#Make figure with heatmap of all key genes

library(babelgene)
library(pheatmap)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
IEG_normTT <- t(IEG_normT)
write.csv(IEG_normTT, "/Users/kbrida/Library/CloudStorage/Box-Box/Brida2022_Reln/Figures/Figure 1 - snRNA-seq/Data&Code/EES&Subcluster/IEG_normTT.csv")

colrange <-  seq(0,3, by = 0.01)
colorpal <- c("#FFFFFF", "#FFFDFD",  "#fffafa", "#fff9f9", "#fff6f5", "#fff2f1", "#ffe8e7",
             "#ffe4e2", "#ffcac7", "#ffccc9", "#ffc8c6", "#ffbdba", "#ff9292", "#ff465a",
             "#FF1549")

mypalette <- colorRampPalette(colorpal)(300)


pheatmap(IEG_normTT,
         color= mypalette,
         cluster_cols=F, 
         cluster_rows=F,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=F, 
         border_color=NA, 
         fontsize_number=6.5,
         legend_breaks=c(seq(0,3, by = 0.5)))
