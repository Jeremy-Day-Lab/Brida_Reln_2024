library(Seurat)
library(ggplot2)

### ---------- Create Drd1-MSN subset object ------------------------------------------------

# Load MCN integrated object
allRats <- readRds("NAc_Combo_Integrated.RDS")

# Create Drd1 subset
Drd1 <- subset(allRats, idents = c("Drd1-MSN-1", "Drd1-MSN-2"))

# Rerun standard processing for dimensionality reduction & clustering
DefaultAssay(Drd1) <- "integrated"
Drd1 <- ScaleData(Drd1)
Drd1 <- RunPCA(Drd1, npcs = 17)
Drd1 <- RunUMAP(Drd1, reduction = "pca", dims = 1:17)
Drd1 <- FindNeighbors(Drd1, reduction = "pca", dims = 1:17)
Drd1 <- FindClusters(Drd1, resolution = 0.1)

# Plot on UMAP
DimPlot(Drd1, reduction = "umap", label = TRUE) + NoLegend()

# Save subset object
saveRDS(Drd1, "D1_subset_of_MCNobj.rds")


### ---------- Reln and Composite IEG Expression -------------------------------------------

# Load MCN D1 subset object
Drd1 <- readRDS("D1_subset_of_MCNobj.rds")

# Plot on UMAP to visualize
DimPlot(Drd1, reduction = "umap", label = T) + coord_fixed() + NoLegend()

# Set IEG genes
gene.set <- c("Arc", "Fos", "Fosb", "Fosl2", "Egr1", "Egr2", "Egr4", "Nr4a1",
              "Nr4a3", "Homer1", "Tiparp", "Junb", "Scg2", "Btg2", "Vgf")
# Set Reln
reln <- "Reln"

# Get mean IEG expression
mean.exp <- colMeans(Drd1@assays$RNA@data[gene.set, ], na.rm = T)
# Get Reln expression
reln.exp <- Drd1@assays$RNA@data[reln, ]

# Add mean expression values to metadata in 'object@meta.data$gene.set.score'
if (all(names(mean.exp) == rownames(Drd1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'Drd1@meta.data':\n", 
      "adding gene set mean expression values in 'Drd1@meta.data$gene.set.score'")
  Drd1@meta.data$gene.set.score <- mean.exp
}

# Add Reln expression values to metadata in 'reln.exp'
if (all(names(reln.exp) == rownames(Drd1@meta.data))) {
  cat("Cell names order match in 'reln.exp' and 'Drd1@meta.data':\n", 
      "adding reln expression values in 'Drd1@meta.data$reln.exp'")
  Drd1@meta.data$Reln.exp <- reln.exp
}

# Plot mean (composite) expression
DefaultAssay(Drd1) <- "RNA"
FeaturePlot(Drd1, "gene.set.score", order = T, cols = c("lightgrey", "#FF1549")) + coord_fixed() + NoAxes()

# Write metadata to file
Drd1.metadata <- Drd1@meta.data
write.csv(Drd1.metadata, "Drd1_MCN_metadata.csv", quote = F)
