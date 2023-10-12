######################################
## SEURAT - INITIAL QC & CLUSTERING ##
######################################

# Load libraries
library(Seurat)
library(SoupX)
library(patchwork)
library(dplyr)
library(ggplot2)

# Load data for SoupX
# Subfolders from CellRanger needed for SoupX: raw_feature_bc_matrix, filtered_feature_bc_matrix, analysis (clustering info)
male1.lacz.data <- load10X("../cellranger/completed/1_lacZ_S1_M1/outs")
male1.reln.data <- load10X("../cellranger/completed/2_Reln_S2_M1/outs")
fem1.lacz.data  <- load10X("../cellranger/completed/3_lacZ_S3_F1/outs")
fem1.reln.data  <- load10X("../cellranger/completed/4_Reln_S4_F1/outs")
fem2.lacz.data  <- load10X("../cellranger/completed/5_lacZ_S5_F2/outs")
fem2.reln.data  <- load10X("../cellranger/completed/6_Reln_S6_F2/outs")
male2.reln.data <- load10X("../cellranger/completed/7_Reln_S7_M2/outs")
male2.lacz.data <- load10X("../cellranger/completed/8_lacZ_S8_M2/outs")

## ---1. INITIAL QC------------------------------------------------------------------------------------------

# Run SoupX: estimate ambient mRNA contamination
male1.lacz.data <- autoEstCont(male1.lacz.data)
male1.reln.data <- autoEstCont(male1.reln.data)
fem1.lacz.data  <- autoEstCont(fem1.lacz.data)
fem1.reln.data  <- autoEstCont(fem1.reln.data) #failed
fem2.lacz.data  <- autoEstCont(fem2.lacz.data)
fem2.reln.data  <- autoEstCont(fem2.reln.data)
male2.reln.data <- autoEstCont(male2.reln.data)
male2.lacz.data <- autoEstCont(male2.lacz.data)

# Run SoupX: infer corrected expression matrix (clean the data of ambient mRNA)
# Rounding to integer is required for DESeq2 later (DEGs)
male1.lacz.clean <- adjustCounts(male1.lacz.data, roundToInt = TRUE)
male1.reln.clean <- adjustCounts(male1.reln.data, roundToInt = TRUE)
fem1.lacz.clean  <- adjustCounts(fem1.lacz.data, roundToInt = TRUE)
# fem1.reln.clean  <- adjustCounts(fem1.reln.data, roundToInt = TRUE)
fem2.lacz.clean  <- adjustCounts(fem2.lacz.data, roundToInt = TRUE)
fem2.reln.clean  <- adjustCounts(fem2.reln.data, roundToInt = TRUE)
male2.reln.clean <- adjustCounts(male2.reln.data, roundToInt = TRUE)
male2.lacz.clean <- adjustCounts(male2.lacz.data, roundToInt = TRUE)

# Create Seurat objects from SoupX output
# Arbitrary cutoffs so we can interrogate quality of every cell and can remove some during later QC
male1.lacz <- CreateSeuratObject(male1.lacz.clean, min.cells = 1, min.features = 1)
male1.reln <- CreateSeuratObject(male1.reln.clean, min.cells = 1, min.features = 1)
fem1.lacz  <- CreateSeuratObject(fem1.lacz.clean, min.cells = 1, min.features = 1)
# fem1.reln  <- CreateSeuratObject(fem1.reln.clean, min.cells = 1, min.features = 1)
fem2.lacz  <- CreateSeuratObject(fem2.lacz.clean, min.cells = 1, min.features = 1)
fem2.reln  <- CreateSeuratObject(fem2.reln.clean, min.cells = 1, min.features = 1)
male2.reln <- CreateSeuratObject(male2.reln.clean, min.cells = 1, min.features = 1)
male2.lacz <- CreateSeuratObject(male2.lacz.clean, min.cells = 1, min.features = 1)

# Import fem1.reln bypassing SoupX
fem1.reln.data <- Read10X("../cellranger/completed/4_Reln_S4_F1/outs/filtered_feature_bc_matrix")
fem1.reln <- CreateSeuratObject(fem1.reln.data, min.cells = 1, min.features = 1)

# Some mt genes in Rn7 reference don't have "Mt-" prefix (e.g., AY172581.x)
# Create vector directly from reference to contain all mt genes
Rn7_gtf <- as.data.frame(rtracklayer::import("../custom_reference/Rattus_norvegicus.mRatBN7.2.109_Chimera.gtf"))
MT_genes <- subset(Rn7_gtf,subset = (seqnames == "MT" & type == "gene"))$gene_name

# ID percentage of reads mapping to mt genes
male1.lacz <- PercentageFeatureSet(male1.lacz, 
                                   features = MT_genes[which(MT_genes %in% rownames(male1.lacz@assays$RNA@counts))], 
                                   col.name = "percent.mt")
male1.reln <- PercentageFeatureSet(male1.reln, 
                                   features = MT_genes[which(MT_genes %in% rownames(male1.reln@assays$RNA@counts))], 
                                   col.name = "percent.mt")
fem1.lacz  <- PercentageFeatureSet(fem1.lacz, 
                                  features = MT_genes[which(MT_genes %in% rownames(fem1.lacz@assays$RNA@counts))], 
                                  col.name = "percent.mt")
fem1.reln  <- PercentageFeatureSet(fem1.reln, 
                                   features = MT_genes[which(MT_genes %in% rownames(fem1.reln@assays$RNA@counts))], 
                                   col.name = "percent.mt")
fem2.lacz  <- PercentageFeatureSet(fem2.lacz, 
                                  features = MT_genes[which(MT_genes %in% rownames(fem2.lacz@assays$RNA@counts))], 
                                  col.name = "percent.mt")
fem2.reln  <- PercentageFeatureSet(fem2.reln, 
                                   features = MT_genes[which(MT_genes %in% rownames(fem2.reln@assays$RNA@counts))], 
                                   col.name = "percent.mt")
male2.reln <- PercentageFeatureSet(male2.reln, 
                                   features = MT_genes[which(MT_genes %in% rownames(male2.reln@assays$RNA@counts))], 
                                   col.name = "percent.mt")
male2.lacz <- PercentageFeatureSet(male2.lacz, 
                                  features = MT_genes[which(MT_genes %in% rownames(male2.lacz@assays$RNA@counts))], 
                                  col.name = "percent.mt")

#Subset data to have >200 features and <5% of reads mapping to mt genes
male1.lacz <- subset(male1.lacz, subset = nFeature_RNA > 200 & percent.mt < 5) #637 nuclei
male1.reln <- subset(male1.reln, subset = nFeature_RNA > 200 & percent.mt < 5) #640 nuclei
fem1.lacz  <- subset(fem1.lacz, subset = nFeature_RNA > 200 & percent.mt < 5) #568 nuclei
fem1.reln  <- subset(fem1.reln, subset = nFeature_RNA > 200 & percent.mt < 5) #513 nuclei
fem2.lacz  <- subset(fem2.lacz, subset = nFeature_RNA > 200 & percent.mt < 5) #922 nuclei
fem2.reln  <- subset(fem2.reln, subset = nFeature_RNA > 200 & percent.mt < 5) #1443 nuclei
male2.reln <- subset(male2.reln, subset = nFeature_RNA > 200 & percent.mt < 5) #598 nuclei
male2.lacz <- subset(male2.lacz, subset = nFeature_RNA > 200 & percent.mt < 5) #900 nuclei

# Add metadata
male1.lacz$target <- "lacz"
male1.reln$target <- "reln"
fem1.lacz$target  <- "lacz"
fem1.reln$target  <- "reln"
fem2.lacz$target  <- "lacz"
fem2.reln$target  <- "reln"
male2.reln$target <- "reln"
male2.lacz$target <- "lacz"

male1.lacz$sex <- "male"
male1.reln$sex <- "male"
fem1.lacz$sex  <- "female"
fem1.reln$sex  <- "female"
fem2.lacz$sex  <- "female"
fem2.reln$sex  <- "female"
male2.reln$sex <- "male"
male2.lacz$sex <- "male"

male1.lacz$sex_target <- "male_lacz"
male1.reln$sex_target <- "male_reln"
fem1.lacz$sex_target  <- "female_lacz"
fem1.reln$sex_target  <- "female_reln"
fem2.lacz$sex_target  <- "female_lacz"
fem2.reln$sex_target  <- "female_reln"
male2.reln$sex_target <- "male_reln"
male2.lacz$sex_target <- "male_lacz"

male1.lacz$GEM <- 1
male1.reln$GEM <- 2
fem1.lacz$GEM  <- 3
fem1.reln$GEM  <- 4
fem2.lacz$GEM  <- 5
fem2.reln$GEM  <- 6
male2.reln$GEM <- 7
male2.lacz$GEM <- 8

# Normalize data under default parameters
male1.lacz <- NormalizeData(male1.lacz)
male1.reln <- NormalizeData(male1.reln)
fem1.lacz  <- NormalizeData(fem1.lacz)
fem1.reln  <- NormalizeData(fem1.reln)
fem2.lacz  <- NormalizeData(fem2.lacz)
fem2.reln  <- NormalizeData(fem2.reln)
male2.reln <- NormalizeData(male2.reln)
male2.lacz <- NormalizeData(male2.lacz)

# ID features that show high cell-to-cell variation in expression
male1.lacz <- FindVariableFeatures(male1.lacz, selection.method = "vst", nfeatures = 3000)
male1.reln <- FindVariableFeatures(male1.reln, selection.method = "vst", nfeatures = 3000)
fem1.lacz  <- FindVariableFeatures(fem1.lacz, selection.method = "vst", nfeatures = 3000)
fem1.reln  <- FindVariableFeatures(fem1.reln, selection.method = "vst", nfeatures = 3000)
fem2.lacz  <- FindVariableFeatures(fem2.lacz, selection.method = "vst", nfeatures = 3000)
fem2.reln  <- FindVariableFeatures(fem2.reln, selection.method = "vst", nfeatures = 3000)
male2.reln <- FindVariableFeatures(male2.reln, selection.method = "vst", nfeatures = 3000)
male2.lacz <- FindVariableFeatures(male2.lacz, selection.method = "vst", nfeatures = 3000)

## ---2. INTEGRATE DATASETS----------------------------------------------------------------------------------

# Save individual dataset objects
saveRDS(male1.lacz, "male1.lacz.rds")
saveRDS(male1.reln, "male1.reln.rds")
saveRDS(fem1.lacz, "fem1.lacz.rds")
saveRDS(fem1.reln, "fem1.reln.rds")
saveRDS(fem2.lacz, "fem2.lacz.rds")
saveRDS(fem2.reln, "fem2.reln.rds")
saveRDS(male2.reln, "male2.reln.rds")
saveRDS(male2.lacz, "male2.lacz.rds")

# Integrate datasets
anchors.allRats <- FindIntegrationAnchors(object.list = list(male1.lacz, 
                                                             male1.reln, 
                                                             fem1.lacz, 
                                                             fem1.reln, 
                                                             fem2.lacz, 
                                                             fem2.reln, 
                                                             male2.reln, 
                                                             male2.lacz), 
                                          dims = 1:17)
allRats <- IntegrateData(anchorset = anchors.allRats, dims = 1:17)

## ---3. INITIAL DIMENSIONALITY REDUCTION & CLUSTERING-------------------------------------------------------

# Reset default assay to integrated
DefaultAssay(allRats) <- "integrated"

# Dimensionality reduction & clustering
allRats <- ScaleData(allRats)
allRats <- RunPCA(allRats, npcs = 17)
allRats <- RunUMAP(allRats, reduction = "pca", dims = 1:17)
allRats <- FindNeighbors(allRats, reduction = "pca", dims = 1:17)
allRats <- FindClusters(allRats, resolution = 0.3)

# Visualize UMAP clustering
DimPlot(object = allRats, reduction = "umap", label = TRUE) + NoLegend()

## ---4. IDENTIFY & LABEL CLUSTERS---------------------------------------------------------------------------

## Make stacked violin plot to confirm cluster identities by marker genes
# Vector of marker genes
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", 
              "Chat", "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", 
              "Rgs5", "Mbp", "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", 
              "Gad1", "Syt1")
# Reset default assay
DefaultAssay(allRats) <- "RNA"
# Stacked violin plot
stackedVln <- VlnPlot(allRats, features, stack = TRUE, sort = FALSE, flip = TRUE, fill.by = "ident") + 
                   theme(legend.position = "none")
stackedVln

# Rename clusters - change names below as needed, based on marker genes
allRats <- RenameIdents(object = allRats,
                          "0" = "Drd2-MSN-1",
                          "1" = "Olig",
                          "2" = "Astrocyte",
                          "3" = "Drd1-MSN-2",
                          "4" = "Drd1-MSN-1",
                          "5" = "Grm8-MSN",
                          "6" = "Microglia",
                          "7" = "GABA-Undef",
                          "8" = "Polydend",
                          "9" = "Unk-1",
                          "10" = "Drd2-MSN-2",
                          "11" = "Pvalb-Int",
                          "12" = "Mural",
                          "13" = "Sst-Int",
                          "14" = "Drd3-MSN",
                          "15" = "Unk-2")

# Re-order the levels
levels(allRats) <- c("Drd1-MSN-1", "Drd1-MSN-2", "Drd2-MSN-1", "Drd2-MSN-2", "Drd3-MSN", 
                     "Grm8-MSN", "GABA-Undef", "Pvalb-Int", "Sst-Int", "Astrocyte", 
                     "Microglia", "Mural", "Olig", "Polydend", "Unk-1", "Unk-2")

DefaultAssay(allRats) <- "integrated"
DimPlot(object = allRats, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

# Save labeled object + UMAP
saveRDS(allRats, "allRats_souped.rds")
