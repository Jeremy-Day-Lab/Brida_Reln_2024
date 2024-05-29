library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(SeuratObject)
library(babelgene)
library(pheatmap)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)

#load NAc_Combo cocaine dataset
NAc_Combo <- readRDS(file = "/Users/kbrida/Desktop/Manuscript/Data&Code/Acute&Repeated/NAc_Combo_Integrated.RDS")

#Load the Reln KD object
RelnKD <- readRDS(file = "/Volumes/Day-Lab$/10X_SingleCell/2023_SC_0016/analysis/Seurat/R objects/allRats_souped_noDoub.rds")

#Change the hyphen to a period
Idents(NAc_Combo) <- gsub(Idents(NAc_Combo),pattern = "-",replacement = ".")
Idents(RelnKD) <- gsub(Idents(RelnKD),pattern = "-",replacement = ".")

#Define marker genes for each cell type
#RelnKD
markers.RelnKD <- FindAllMarkers(RelnKD, only.pos = TRUE, min.pct = 0.25)
# Save significant markers for each cluster
top_KD <- subset(markers.RelnKD, subset = markers.RelnKD$p_val_adj<0.05)
#NAC_Combo
markers.NAc_Combo <- FindAllMarkers(NAc_Combo, only.pos = TRUE, min.pct = 0.25)
# Save significant markers for each cluster
top_combo <- subset(markers.NAc_Combo, subset = markers.NAc_Combo$p_val_adj<0.05)



######Reln_KD#####
Drd1_KD <- data_frame(gene = ifelse(top_KD$cluster == "Drd1.MSN", top_KD$gene, NA),
                     Drd1_KD = ifelse(top_KD$cluster == "Drd1.MSN", top_KD$avg_log2FC, NA))
Drd1_KD <- na.omit(Drd1_KD)


Drd2.1_KD <- data_frame(gene = ifelse(top_KD$cluster == "Drd2.MSN.1", top_KD$gene, NA),
                   Drd2.1_KD = ifelse(top_KD$cluster == "Drd2.MSN.1", top_KD$avg_log2FC, NA))
Drd2.1_KD <- na.omit(Drd2.1_KD)

Drd2.2_KD <- data_frame(gene = ifelse(top_KD$cluster == "Drd2.MSN.2", top_KD$gene, NA),
                     Drd2.2_KD = ifelse(top_KD$cluster == "Drd2.MSN.2", top_KD$avg_log2FC, NA))
Drd2.2_KD <- na.omit(Drd2.2_KD)


Drd3_KD <- data_frame(gene = ifelse(top_KD$cluster == "Drd3.MSN", top_KD$gene, NA),
                   Drd3_KD = ifelse(top_KD$cluster == "Drd3.MSN", top_KD$avg_log2FC, NA))
Drd3_KD <- na.omit(Drd3_KD)


GABA_KD <- data_frame(gene = ifelse(top_KD$cluster == "GABA.Undef", top_KD$gene, NA),
                      GABA_KD = ifelse(top_KD$cluster == "GABA.Undef", top_KD$avg_log2FC, NA))
GABA_KD <- na.omit(GABA_KD)


Grm8_KD <- data_frame(gene = ifelse(top_KD$cluster == "Grm8.MSN", top_KD$gene, NA),
                      Grm8_KD = ifelse(top_KD$cluster == "Grm8.MSN", top_KD$avg_log2FC, NA))
Grm8_KD <- na.omit(Grm8_KD)


Sst_KD <- data_frame(gene = ifelse(top_KD$cluster == "Sst.Int", top_KD$gene, NA),
                      Sst_KD = ifelse(top_KD$cluster == "Sst.Int", top_KD$avg_log2FC, NA))
Sst_KD <- na.omit(Sst_KD)


Pvalb_KD <- data_frame(gene = ifelse(top_KD$cluster == "Pvalb.Int", top_KD$gene, NA),
                     Pvalb_KD = ifelse(top_KD$cluster == "Pvalb.Int", top_KD$avg_log2FC, NA))
Pvalb_KD <- na.omit(Pvalb_KD)


Microglia_KD <- data_frame(gene = ifelse(top_KD$cluster == "Microglia", top_KD$gene, NA),
                       Microglia_KD = ifelse(top_KD$cluster == "Microglia", top_KD$avg_log2FC, NA))
Microglia_KD <- na.omit(Microglia_KD)


Astrocyte_KD <- data_frame(gene = ifelse(top_KD$cluster == "Astrocyte", top_KD$gene, NA),
                           Astrocyte_KD = ifelse(top_KD$cluster == "Astrocyte", top_KD$avg_log2FC, NA))
Astrocyte_KD <- na.omit(Astrocyte_KD)


Olig_KD <- data_frame(gene = ifelse(top_KD$cluster == "Olig", top_KD$gene, NA),
                           Olig_KD = ifelse(top_KD$cluster == "Olig", top_KD$avg_log2FC, NA))
Olig_KD <- na.omit(Olig_KD)


Polydend_KD <- data_frame(gene = ifelse(top_KD$cluster == "Polydend", top_KD$gene, NA),
                      Polydend_KD = ifelse(top_KD$cluster == "Polydend", top_KD$avg_log2FC, NA))
Polydend_KD <- na.omit(Polydend_KD)


Mural_KD <- data_frame(gene = ifelse(top_KD$cluster == "Mural", top_KD$gene, NA),
                          Mural_KD = ifelse(top_KD$cluster == "Mural", top_KD$avg_log2FC, NA))
Mural_KD <- na.omit(Mural_KD)



#Okay, need to generate a matrix that is the gene name x cell type with values being log2FC
#Then we can do the statistical testing

#"iterative" joining apologies to those who are offended by the clunkiness
test_join <- merge(Drd1_KD, Drd2.1_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Drd2.2_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Drd3_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Grm8_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, GABA_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Sst_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Pvalb_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Microglia_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Astrocyte_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Olig_KD, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Polydend_KD, all.x = TRUE, all.y = TRUE)
KD_DEGs <- merge(test_join, Mural_KD, all.x = TRUE, all.y = TRUE)
#replace NA with 0
KD_DEGs[is.na(KD_DEGs)] <- 0

#####NAc_Combo#####
Drd1.1_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Drd1.MSN.1", top_combo$gene, NA),
                      Drd1.1_Combo = ifelse(top_combo$cluster == "Drd1.MSN.1", top_combo$avg_log2FC, NA))
Drd1.1_Combo <- na.omit(Drd1.1_Combo)

Drd1.2_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Drd1.MSN.2", top_combo$gene, NA),
                          Drd1.2_Combo = ifelse(top_combo$cluster == "Drd1.MSN.2", top_combo$avg_log2FC, NA))
Drd1.2_Combo <- na.omit(Drd1.2_Combo)


Drd2.1_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Drd2.MSN.1", top_combo$gene, NA),
                        Drd2.1_Combo = ifelse(top_combo$cluster == "Drd2.MSN.1", top_combo$avg_log2FC, NA))
Drd2.1_Combo <- na.omit(Drd2.1_Combo)

Drd2.2_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Drd2.MSN.2", top_combo$gene, NA),
                        Drd2.2_Combo = ifelse(top_combo$cluster == "Drd2.MSN.2", top_combo$avg_log2FC, NA))
Drd2.2_Combo <- na.omit(Drd2.2_Combo)


Drd3_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Drd3.MSN", top_combo$gene, NA),
                      Drd3_Combo = ifelse(top_combo$cluster == "Drd3.MSN", top_combo$avg_log2FC, NA))
Drd3_Combo <- na.omit(Drd3_Combo)


GABA_Combo <- data_frame(gene = ifelse(top_combo$cluster == "GABAergic", top_combo$gene, NA),
                      GABA_Combo = ifelse(top_combo$cluster == "GABAergic", top_combo$avg_log2FC, NA))
GABA_Combo <- na.omit(GABA_Combo)


Grm8_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Grm8.MSN", top_combo$gene, NA),
                      Grm8_Combo = ifelse(top_combo$cluster == "Grm8.MSN", top_combo$avg_log2FC, NA))
Grm8_Combo <- na.omit(Grm8_Combo)


Sst_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Sst.Interneuron", top_combo$gene, NA),
                     Sst_Combo = ifelse(top_combo$cluster == "Sst.Interneuron", top_combo$avg_log2FC, NA))
Sst_Combo <- na.omit(Sst_Combo)


Pvalb_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Pvalb.Interneuron", top_combo$gene, NA),
                       Pvalb_Combo = ifelse(top_combo$cluster == "Pvalb.Interneuron", top_combo$avg_log2FC, NA))
Pvalb_Combo <- na.omit(Pvalb_Combo)


Microglia_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Microglia", top_combo$gene, NA),
                           Microglia_Combo = ifelse(top_combo$cluster == "Microglia", top_combo$avg_log2FC, NA))
Microglia_Combo <- na.omit(Microglia_Combo)


Astrocyte_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Astrocyte", top_combo$gene, NA),
                           Astrocyte_Combo = ifelse(top_combo$cluster == "Astrocyte", top_combo$avg_log2FC, NA))
Astrocyte_Combo <- na.omit(Astrocyte_Combo)
#Astrocyte_Combo <- as.list(Astrocyte_Combo)

Olig_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Olig.1", top_combo$gene, NA),
                      Olig_Combo = ifelse(top_combo$cluster == "Olig.1", top_combo$avg_log2FC, NA))
Olig_Combo <- na.omit(Olig_Combo)


Polydend_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Polydendrocyte", top_combo$gene, NA),
                          Polydend_Combo = ifelse(top_combo$cluster == "Polydendrocyte", top_combo$avg_log2FC, NA))
Polydend_Combo <- na.omit(Polydend_Combo)

Mural_Combo <- data_frame(gene = ifelse(top_combo$cluster == "Mural", top_combo$gene, NA),
                       Mural_Combo = ifelse(top_combo$cluster == "Mural", top_combo$avg_log2FC, NA))
Mural_Combo <- na.omit(Mural_Combo)


#Okay, need to generate a matrix that is the gene name x cell type with values being log2FC
#Then we can do the statistical testing

#"iterative" joining apologies to those who are offended by the clunkiness
test_join <- merge(Drd1.1_Combo, Drd1.2_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Drd2.1_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Drd2.2_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Drd3_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Grm8_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, GABA_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Sst_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Pvalb_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Microglia_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Astrocyte_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Olig_Combo, all.x = TRUE, all.y = TRUE)
test_join <- merge(test_join, Polydend_Combo, all.x = TRUE, all.y = TRUE)
Combo_DEGs <- merge(test_join, Mural_Combo, all.x = TRUE, all.y = TRUE)
#replace NA with 0
Combo_DEGs[is.na(Combo_DEGs)] <- 0

#Now merge the two dataframes
all_DEGs <- merge(x    = Combo_DEGs,
                        y    = KD_DEGs,
                        by.x = "gene",
                        by.y = "gene")


Cor_mat <- cor(all_DEGs[,grep(colnames(all_DEGs),pattern = "_Combo")],all_DEGs[,grep(colnames(all_DEGs),pattern = "_KD")])

colrange <-  seq(-.60,.60, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))

pheatmap(Cor_mat,
         color=colorpal,
         cluster_cols=F, 
         cluster_rows=F,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.6,.6, by = 0.3)))


sessionInfo()