#################################################
##                                             ##
## DOUBLETFINDER - FINDING & REMOVING DOUBLETS ##
##                                             ##
#################################################

# Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(DoubletFinder)

# Load objects from initial QC
male1.lacz <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/male1.lacz.rds")
male1.reln <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/male1.reln.rds")
fem1.lacz  <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/fem1.lacz.rds")
fem1.reln  <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/fem1.reln.rds")
fem2.lacz  <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/fem2.lacz.rds")
fem2.reln  <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/fem2.reln.rds")
male2.reln <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/male2.reln.rds")
male2.lacz <- readRDS("0_expanded_custom_reference/individual_datasets_processed/3K_variable_features/male2.lacz.rds")
allRats <- readRDS("allRats_souped.rds")

## ---1. Set up combined obj for metadata retrieval----------------------------------------------------------

# Make column containing non-unique cell names.
# Seurat forces unique cell names by adding _ + integer. Splitting by _ gives original cell ID.
allRats$nonUniqueCellName <- as.character(lapply(strsplit(rownames(allRats@meta.data),split = "_"),"[",1))

# Make column that IDs forced unique cell names by writing out rownames
allRats$uniqueCellName <- row.names(allRats@meta.data)

# Add column with cell types
allRats$cellType <- Idents(allRats)

## ------------------------------------------------MALE 1 LACZ-------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
ML1.metadata <- subset(allRats, subset = (GEM == 1))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
male1.lacz$nonUniqueCellName <- rownames(male1.lacz@meta.data)
ML1 <- merge(x  = male1.lacz@meta.data, y = ML1.metadata, by = "nonUniqueCellName")
rownames(ML1) <- ML1$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(ML1 %>% distinct()) == nrow(ML1.metadata)
ML1 <- ML1[match(row.names(male1.lacz@meta.data), row.names(ML1)), ]
all(row.names(ML1) == row.names(male1.lacz@meta.data))
all(length(row.names(ML1)) == length(row.names(male1.lacz@meta.data)))

# Add cell type to metadata
male1.lacz <- AddMetaData(male1.lacz, metadata = ML1$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

male1.lacz <- ScaleData(male1.lacz)
male1.lacz <- RunPCA(male1.lacz, npcs = 17)
male1.lacz <- RunUMAP(male1.lacz, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
male1.lacz.sweep <- paramSweep_v3(male1.lacz, PCs = 1:17, sct = FALSE)

# Summarize sweep
male1.lacz.sweepStats <- summarizeSweep(male1.lacz.sweep, GT = FALSE)

# ID & plot pK
ML1_pk <- find.pK(male1.lacz.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
ML1_pk <- as.numeric(as.character(ML1_pk[which(ML1_pk$BCmetric == max(ML1_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
ML1.homotypic_props <- modelHomotypic(male1.lacz$cellType)
ML1.nExp_poi <- round(.04*nrow(male1.lacz@meta.data))
ML1.nExp_poi_adj <- round(ML1.nExp_poi*(1-ML1.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
male1.lacz <- doubletFinder_v3(male1.lacz, PCs = 1:17, pN = 0.25, pK = ML1_pk, nExp = ML1.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(male1.lacz@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "ML1"))
rm(male1.lacz.sweep, male1.lacz.sweepStats)

## ------------------------------------------------MALE 1 RELN-------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
MR1.metadata <- subset(allRats, subset = (GEM == 2))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
male1.reln$nonUniqueCellName <- rownames(male1.reln@meta.data)
MR1 <- merge(x  = male1.reln@meta.data, y  = MR1.metadata, by = "nonUniqueCellName")
rownames(MR1) <- MR1$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(MR1 %>% distinct()) == nrow(MR1.metadata)
MR1 <- MR1[match(row.names(male1.reln@meta.data), row.names(MR1)), ]
all(row.names(MR1) == row.names(male1.reln@meta.data)) 
all(length(row.names(MR1)) == length(row.names(male1.reln@meta.data)))

# Add cell type to metadata
male1.reln <- AddMetaData(male1.reln, metadata = MR1$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

male1.reln <- ScaleData(male1.reln)
male1.reln <- RunPCA(male1.reln, npcs = 17)
male1.reln <- RunUMAP(male1.reln, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
male1.reln.sweep <- paramSweep_v3(male1.reln, PCs = 1:17, sct = FALSE)

# Summarize sweep
male1.reln.sweepStats <- summarizeSweep(male1.reln.sweep, GT = FALSE)

# ID & plot pK
MR1_pk <- find.pK(male1.reln.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
MR1_pk <- as.numeric(as.character(MR1_pk[which(MR1_pk$BCmetric == max(MR1_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
MR1.homotypic_props <- modelHomotypic(male1.reln$cellType)
MR1.nExp_poi <- round(.04*nrow(male1.reln@meta.data))
MR1.nExp_poi_adj <- round(MR1.nExp_poi*(1-MR1.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
male1.reln <- doubletFinder_v3(male1.reln, PCs = 1:17, pN = 0.25, pK = MR1_pk, nExp = MR1.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(male1.reln@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "MR1"))
rm(male1.reln.sweep, male1.reln.sweepStats)

## ---------------------------------------------FEMALE 1 LACZ------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull subset metadata from combo obj: non-unique cell name, cell type, unique cell name forced by Seurat
FL1.metadata <- subset(allRats, subset = (GEM == 3))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]

# Subset obj should have rownames that correspond to non-unique cell names
# Make that a column; will be used to merge in the next command
fem1.lacz$nonUniqueCellName <- rownames(fem1.lacz@meta.data)

# Merge the two dataframes by non-unique cell name
FL1 <- merge(x  = fem1.lacz@meta.data, y  = FL1.metadata, by = "nonUniqueCellName")

# Make rownames the non-unique cell name
rownames(FL1) <- FL1$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

# Make sure nothing is duplicated when merging.
# If nothing is duplicated, # of unique rows should = # of rows when pulling subset cells from allRats metadata
nrow(FL1 %>% distinct()) == nrow(FL1.metadata)

# Now make sure FL is in the same order as original subset & same length
# We want FL to be in same order as subset, so put subset first & FL second.
FL1 <- FL1[match(row.names(fem1.lacz@meta.data), row.names(FL1)), ]
all(row.names(FL1) == row.names(fem1.lacz@meta.data))
all(length(row.names(FL1)) == length(row.names(fem1.lacz@meta.data)))

# Add cell type to metadata
fem1.lacz <- AddMetaData(fem1.lacz, metadata = FL1$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

# DoubletFinder requires each obj to be scaled + dimensionality reduction
# Use same # of PCs as in processing of integrated obj
fem1.lacz <- ScaleData(fem1.lacz)
fem1.lacz <- RunPCA(fem1.lacz, npcs = 17)
fem1.lacz <- RunUMAP(fem1.lacz, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# paramSweep_v3 tests multiple pN & pK parameters to ID optimal pK value. pK is main parameter that needs to be tuned
fem1.lacz.sweep <- paramSweep_v3(fem1.lacz, PCs = 1:17, sct = FALSE)

# Summarize sweep. GT = ground truth, which we don't have here
fem1.lacz.sweepStats <- summarizeSweep(fem1.lacz.sweep, GT = FALSE)

# ID & plot pK
FL1_pk <- find.pK(fem1.lacz.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
FL1_pk <- as.numeric(as.character(FL1_pk[which(FL1_pk$BCmetric == max(FL1_pk$BCmetric)), "pK"]))

# Model proportion of homotypic doublets. This calculates sum of squares of proportions of annotations
FL1.homotypic_props <- modelHomotypic(fem1.lacz$cellType)

# Total number of doublets (upper bound): multiply # of cells by estimated multiplet rate
# Estimated multiplet rate (0.04 for 8K cells loaded) on p. 18 of 10X Chromium Next GEM Single Cell 3' Reagent Kits v3.1 Dual Index manual
FL1.nExp_poi <- round(.04*nrow(fem1.lacz@meta.data))
# Number of heterotypic doublets only (lower bound): total doublets * (1 - proportion homotypic)
FL1.nExp_poi_adj <- round(FL1.nExp_poi*(1-FL1.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj because DoubletFinder only IDs homotypic doublets, so using nExp_poi may remove some singlets.
fem1.lacz <- doubletFinder_v3(fem1.lacz, PCs = 1:17, pN = 0.25, pK = FL1_pk, nExp = FL1.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(fem1.lacz@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove all objs except original subsets
# (SAVE SWEEP OUTPUTS TO FILE FIRST??)
rm(list = ls(pattern = "FL1"))
rm(fem1.lacz.sweep, fem1.lacz.sweepStats)

## -----------------------------------------------FEMALE 1 RELN------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
FR1.metadata <- subset(allRats, subset = (GEM == 4))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
fem1.reln$nonUniqueCellName <- rownames(fem1.reln@meta.data)
FR1 <- merge(x  = fem1.reln@meta.data, y  = FR1.metadata, by = "nonUniqueCellName")
rownames(FR1) <- FR1$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(FR1 %>% distinct()) == nrow(FR1.metadata)
FR1 <- FR1[match(row.names(fem1.reln@meta.data), row.names(FR1)), ]
all(row.names(FR1) == row.names(fem1.reln@meta.data)) 
all(length(row.names(FR1)) == length(row.names(fem1.reln@meta.data)))

# Add cell type to metadata
fem1.reln <- AddMetaData(fem1.reln, metadata = FR1$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

fem1.reln <- ScaleData(fem1.reln)
fem1.reln <- RunPCA(fem1.reln, npcs = 17)
fem1.reln <- RunUMAP(fem1.reln, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
fem1.reln.sweep <- paramSweep_v3(fem1.reln, PCs = 1:17, sct = FALSE)

# Summarize sweep
fem1.reln.sweepStats <- summarizeSweep(fem1.reln.sweep, GT = FALSE)

# ID & plot pK
FR1_pk <- find.pK(fem1.reln.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
FR1_pk <- as.numeric(as.character(FR1_pk[which(FR1_pk$BCmetric == max(FR1_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
FR1.homotypic_props <- modelHomotypic(fem1.reln$cellType)
FR1.nExp_poi <- round(.04*nrow(fem1.reln@meta.data))
FR1.nExp_poi_adj <- round(FR1.nExp_poi*(1-FR1.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
fem1.reln <- doubletFinder_v3(fem1.reln, PCs = 1:17, pN = 0.25, pK = FR1_pk, nExp = FR1.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(fem1.reln@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "FR1"))
rm(fem1.reln.sweep, fem1.reln.sweepStats)

## ---------------------------------------------FEMALE 2 LACZ------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
FL2.metadata <- subset(allRats, subset = (GEM == 5))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
fem2.lacz$nonUniqueCellName <- rownames(fem2.lacz@meta.data)
FL2 <- merge(x  = fem2.lacz@meta.data, y = FL2.metadata, by = "nonUniqueCellName")
rownames(FL2) <- FL2$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(FL2 %>% distinct()) == nrow(FL2.metadata)
FL2 <- FL2[match(row.names(fem2.lacz@meta.data), row.names(FL2)), ]
all(row.names(FL2) == row.names(fem2.lacz@meta.data))
all(length(row.names(FL2)) == length(row.names(fem2.lacz@meta.data)))

# Add cell type to metadata
fem2.lacz <- AddMetaData(fem2.lacz, metadata = FL2$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

fem2.lacz <- ScaleData(fem2.lacz)
fem2.lacz <- RunPCA(fem2.lacz, npcs = 17)
fem2.lacz <- RunUMAP(fem2.lacz, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
fem2.lacz.sweep <- paramSweep_v3(fem2.lacz, PCs = 1:17, sct = FALSE)

# Summarize sweep
fem2.lacz.sweepStats <- summarizeSweep(fem2.lacz.sweep, GT = FALSE)

# ID & plot pK
FL2_pk <- find.pK(fem2.lacz.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
FL2_pk <- as.numeric(as.character(FL2_pk[which(FL2_pk$BCmetric == max(FL2_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
FL2.homotypic_props <- modelHomotypic(fem2.lacz$cellType)
FL2.nExp_poi <- round(.04*nrow(fem2.lacz@meta.data))
FL2.nExp_poi_adj <- round(FL2.nExp_poi*(1-FL2.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
fem2.lacz <- doubletFinder_v3(fem2.lacz, PCs = 1:17, pN = 0.25, pK = FL2_pk, nExp = FL2.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(fem2.lacz@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "FL2"))
rm(fem2.lacz.sweep, fem2.lacz.sweepStats)

## -----------------------------------------------FEMALE 2 RELN------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
FR2.metadata <- subset(allRats, subset = (GEM == 6))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
fem2.reln$nonUniqueCellName <- rownames(fem2.reln@meta.data)
FR2 <- merge(x  = fem2.reln@meta.data, y  = FR2.metadata, by = "nonUniqueCellName")
rownames(FR2) <- FR2$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(FR2 %>% distinct()) == nrow(FR2.metadata)
FR2 <- FR2[match(row.names(fem2.reln@meta.data), row.names(FR2)), ]
all(row.names(FR2) == row.names(fem2.reln@meta.data)) 
all(length(row.names(FR2)) == length(row.names(fem2.reln@meta.data)))

# Add cell type to metadata
fem2.reln <- AddMetaData(fem2.reln, metadata = FR2$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

fem2.reln <- ScaleData(fem2.reln)
fem2.reln <- RunPCA(fem2.reln, npcs = 17)
fem2.reln <- RunUMAP(fem2.reln, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
fem2.reln.sweep <- paramSweep_v3(fem2.reln, PCs = 1:17, sct = FALSE)

# Summarize sweep
fem2.reln.sweepStats <- summarizeSweep(fem2.reln.sweep, GT = FALSE)

# ID & plot pK
FR2_pk <- find.pK(fem2.reln.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
FR2_pk <- as.numeric(as.character(FR2_pk[which(FR2_pk$BCmetric == max(FR2_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
FR2.homotypic_props <- modelHomotypic(fem2.reln$cellType)
FR2.nExp_poi <- round(.04*nrow(fem2.reln@meta.data))
FR2.nExp_poi_adj <- round(FR2.nExp_poi*(1-FR2.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
fem2.reln <- doubletFinder_v3(fem2.reln, PCs = 1:17, pN = 0.25, pK = FR2_pk, nExp = FR2.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(fem2.reln@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "FR2"))
rm(fem2.reln.sweep, fem2.reln.sweepStats)

## ------------------------------------------------MALE 2 RELN-------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
MR2.metadata <- subset(allRats, subset = (GEM == 7))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
male2.reln$nonUniqueCellName <- rownames(male2.reln@meta.data)
MR2 <- merge(x  = male2.reln@meta.data, y  = MR2.metadata, by = "nonUniqueCellName")
rownames(MR2) <- MR2$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(MR2 %>% distinct()) == nrow(MR2.metadata)
MR2 <- MR2[match(row.names(male2.reln@meta.data), row.names(MR2)), ]
all(row.names(MR2) == row.names(male2.reln@meta.data)) 
all(length(row.names(MR2)) == length(row.names(male2.reln@meta.data)))

# Add cell type to metadata
male2.reln <- AddMetaData(male2.reln, metadata = MR2$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

male2.reln <- ScaleData(male2.reln)
male2.reln <- RunPCA(male2.reln, npcs = 17)
male2.reln <- RunUMAP(male2.reln, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
male2.reln.sweep <- paramSweep_v3(male2.reln, PCs = 1:17, sct = FALSE)

# Summarize sweep
male2.reln.sweepStats <- summarizeSweep(male2.reln.sweep, GT = FALSE)

# ID & plot pK
MR2_pk <- find.pK(male2.reln.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
MR2_pk <- as.numeric(as.character(MR2_pk[which(MR2_pk$BCmetric == max(MR2_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
MR2.homotypic_props <- modelHomotypic(male2.reln$cellType)
MR2.nExp_poi <- round(.04*nrow(male2.reln@meta.data))
MR2.nExp_poi_adj <- round(MR2.nExp_poi*(1-MR2.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
male2.reln <- doubletFinder_v3(male2.reln, PCs = 1:17, pN = 0.25, pK = MR2_pk, nExp = MR2.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(male2.reln@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(list = ls(pattern = "MR2"))
rm(male2.reln.sweep, male2.reln.sweepStats)

## ------------------------------------------------MALE 2 LACZ-------------------------------------------------
## ---2. Pull metadata & merge-------------------------------------------------------------------------------

# Pull & merge metadata
ML2.metadata <- subset(allRats, subset = (GEM == 8))@meta.data[,c("nonUniqueCellName", "cellType", "uniqueCellName")]
male2.lacz$nonUniqueCellName <- rownames(male2.lacz@meta.data)
ML2 <- merge(x  = male2.lacz@meta.data, y = ML2.metadata, by = "nonUniqueCellName")
rownames(ML2) <- ML2$nonUniqueCellName

## ---3. Verify correct merging------------------------------------------------------------------------------

nrow(ML2 %>% distinct()) == nrow(ML2.metadata)
ML2 <- ML2[match(row.names(male2.lacz@meta.data), row.names(ML2)), ]
all(row.names(ML2) == row.names(male2.lacz@meta.data))
all(length(row.names(ML2)) == length(row.names(male2.lacz@meta.data)))

# Add cell type to metadata
male2.lacz <- AddMetaData(male2.lacz, metadata = ML2$cellType, col.name = "cellType")

## ---4. Scaling & dimensionality reduction------------------------------------------------------------------

male2.lacz <- ScaleData(male2.lacz)
male2.lacz <- RunPCA(male2.lacz, npcs = 17)
male2.lacz <- RunUMAP(male2.lacz, reduction = "pca", dims = 1:17)

## ---5. DoubletFinder: parameter selection------------------------------------------------------------------

# Optimize pN & pK parameters
male2.lacz.sweep <- paramSweep_v3(male2.lacz, PCs = 1:17, sct = FALSE)

# Summarize sweep
male2.lacz.sweepStats <- summarizeSweep(male2.lacz.sweep, GT = FALSE)

# ID & plot pK
ML2_pk <- find.pK(male2.lacz.sweepStats)

# ID pK value corresponding to max mean-variance-normalized-bimodality coefficient
ML2_pk <- as.numeric(as.character(ML2_pk[which(ML2_pk$BCmetric == max(ML2_pk$BCmetric)), "pK"]))

# Calculate expected number of doublets
ML2.homotypic_props <- modelHomotypic(male2.lacz$cellType)
ML2.nExp_poi <- round(.04*nrow(male2.lacz@meta.data))
ML2.nExp_poi_adj <- round(ML2.nExp_poi*(1-ML2.homotypic_props))

## ---6. DoubletFinder: find doublets------------------------------------------------------------------------

# Find doublets using less-strict nExp_poi_adj
male2.lacz <- doubletFinder_v3(male2.lacz, PCs = 1:17, pN = 0.25, pK = ML2_pk, nExp = ML2.nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

# Change names of dataframe so everything can be merged later
# MAKE SURE TO RENAME THE CORRECT COLUMNS - CHECK THE NUMBERS (11,12)
names(male2.lacz@meta.data)[c(11,12)] <- c("pANN","DF.Classification")

# Remove objs
rm(ML2, ML2.metadata, male2.lacz.sweep, male2.lacz.sweepStats)

####################################################
## CALCULATE PROPORTION OF DOUBLETS PER CELL TYPE ##
####################################################

# Merge all the individual object metadata
combo.meta <- rbind(male1.lacz@meta.data, 
                    male1.reln@meta.data, 
                    fem1.lacz@meta.data, 
                    fem1.reln@meta.data, 
                    fem2.lacz@meta.data, 
                    fem2.reln@meta.data, 
                    male2.reln@meta.data, 
                    male2.lacz@meta.data)

# Check that everything matches
all(combo.meta$nonUniqueCellName == allRats$nonUniqueCellName)
all(nrow(combo.meta) == nrow(allRats@meta.data))

# Add doublet classification to metadata
allRats <- AddMetaData(allRats, metadata = combo.meta$DF.Classification, col.name = "doubletClassification")

# Plot doublets vs. singlets in UMAP space
Idents(allRats) <- allRats$doubletClassification
DimPlot(object = allRats, reduction = "umap")

# Reset Idents to cell type
Idents(allRats) <- allRats$cellType

#Create new dataframe
proportion <- data.frame(cellType = NA, Num_Doublets = NA, Num_Singlets = NA)

for(i in levels(Idents(allRats))){
  counts <- data.frame(cellType = i,
                       Num_Doublets = nrow(subset(allRats@meta.data, subset = (cellType == i & doubletClassification == 'Doublet'))),
                       Num_Singlets = nrow(subset(allRats@meta.data, subset = (cellType == i & doubletClassification == 'Singlet'))))
  
  proportion <- rbind(proportion, counts)
}

# Add proportion values
proportion <- transform(proportion, propDoublets = Num_Doublets/(Num_Doublets + Num_Singlets))
proportion <- na.omit(proportion)

# Create new dataframe with only doublet proportions
doubletProportions <- data.frame(cellType = proportion$cellType, propDoublets = proportion$propDoublets)

# Create bar graph to visualize
bargraph <- ggplot(doubletProportions, aes(x = cellType, y = propDoublets)) + 
                        geom_bar(stat="identity", color = 'blue', fill = 'blue') + 
                        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
                              plot.title = element_text(hjust = 0.5)) + 
                        labs(title="Proportion of doublets in each cell cluster")

#############################################################
## REMOVING DOUBLETS & RECLUSTERING & FINDING MARKER GENES ##
#############################################################

# Subset to remove doublets
allRats.noDoublets <- subset(allRats, subset = (doubletClassification == 'Singlet'))

# Reset default assay
DefaultAssay(allRats.noDoublets) <- "integrated"

# Data scaling, dimensionality reduction, clustering
allRats.noDoublets <- ScaleData(allRats.noDoublets)
allRats.noDoublets <- RunPCA(allRats.noDoublets, npcs = 17)
allRats.noDoublets <- RunUMAP(allRats.noDoublets, reduction = "pca", dims = 1:17)
allRats.noDoublets <- FindNeighbors(allRats.noDoublets, reduction = "pca", dims = 1:17)
allRats.noDoublets <- FindClusters(allRats.noDoublets, resolution = 0.15)

# Plot
DimPlot(allRats.noDoublets, reduction = "umap", label = TRUE) + NoLegend()

# Make stacked violin plot to confirm cluster identities by marker genes
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", "Chat", "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", "Rgs5", "Mbp", "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", "Gad1", "Syt1")
DefaultAssay(allRats.noDoublets) <- "RNA"
stackedVln <- VlnPlot(allRats.noDoublets, features, stack = TRUE, sort = FALSE, flip = TRUE, fill.by = "ident") + 
                theme(legend.position = "none")

##################################
# HORIZONTAL STACKED VIOLIN PLOT #
##################################

library(cowplot)
library(sctree) #For creating dataframe from Seurat obj

DefaultAssay(allRats.noDoublets) <- "RNA"
# Vector of genes
features <- c("Drd1", "Ebf1", "Drd2", "Penk", "Drd3", "Grm8", "Elavl2", "Chat", 
              "Kit", "Sst", "Slc17a7", "Aqp4", "Gja1", "Arhgap15", "Rgs5", "Mbp", 
              "Opalin", "Pdgfra", "Ppp1r1b", "Foxp2", "Bcl11b", "Gad1", "Syt1")
# Create dataframe from Seurat object
allRats.noDoublets.DF <- as.data.frame(allRats.noDoublets, genes = features, fix_names = FALSE)

# Add cell ID and identity classes
allRats.noDoublets.DF$Cell <- rownames(allRats.noDoublets.DF)
allRats.noDoublets.DF$Idents <- Idents(allRats.noDoublets)

# Use melt to change data.frame format
allRats.noDoublets.DF <- reshape2::melt(allRats.noDoublets.DF, id.vars = c("Cell","Idents"), measure.vars = features, 
                                        variable.name = "Feat", value.name = "Expr")
head(allRats.noDoublets.DF, 10)

horizVln <- ggplot(allRats.noDoublets.DF, aes(factor(Feat), Expr, fill = Idents)) +
               geom_violin(scale = "width", adjust = 1, trim = TRUE) +
               scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                                  c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
               facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
               theme_cowplot(font_size = 14) +
               theme(legend.position = "none", panel.spacing = unit(0, "lines"),
                     plot.title = element_text(hjust = 0.5),
                     panel.background = element_rect(fill = NA, color = "black"),
                     strip.background = element_blank(),
                     strip.text = element_text(face = "bold"),
                     strip.text.y.left = element_text(size = 14, angle = 0),
                     axis.text.x = element_text(size = 14, angle = 60, hjust = 1)) +
               xlab("Gene") + ylab("Expression Level")

#######################
# ANNOTATE CELL TYPES #
#######################

# Rename clusters - change names below as needed, based on marker genes
allRats.noDoublets <- RenameIdents(object = allRats.noDoublets,
                          "0" = "Drd1-MSN",
                          "1" = "Drd2-MSN-1",
                          "2" = "Olig",
                          "3" = "Astrocyte",
                          "4" = "Grm8-MSN",
                          "5" = "Microglia",
                          "6" = "GABA-Undef",
                          "7" = "Polydend",
                          "8" = "Drd2-MSN-2",
                          "9" = "Pvalb-Int",
                          "10" = "Mural",
                          "11" = "Drd3-MSN",
                          "12" = "Sst-Int")

levels(allRats.noDoublets) <- c("Drd1-MSN", "Drd2-MSN-1", "Drd2-MSN-2", "Drd3-MSN", "Grm8-MSN", "GABA-Undef", "Pvalb-Int", "Sst-Int", "Astrocyte", "Microglia", "Mural", "Olig", "Polydend")

DimPlot(object = allRats.noDoublets, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = allRats.noDoublets, reduction = "umap") + NoLegend()

# Save labeled object + UMAP
saveRDS(allRats.noDoublets, "allRats_souped_noDoub.rds")
