library(Seurat)

allRats <- readRDS("allRats_souped_noDoub.rds")
allRats$cellType <- Idents(allRats)

# Generate cellType subsets
astro <- subset(allRats, idents = "Astrocyte")
d1 <- subset(allRats, idents = "Drd1-MSN")
d2_1 <- subset(allRats, idents = "Drd2-MSN-1")
d2_2 <- subset(allRats, idents = "Drd2-MSN-2")
d3 <- subset(allRats, idents = "Drd3-MSN")
grm8 <- subset(allRats, idents = "Grm8-MSN")
gaba <- subset(allRats, idents = "GABA-Undef")
pvalb <- subset(allRats, idents = "Pvalb-Int")
sst <- subset(allRats, idents = "Sst-Int")
micro <- subset(allRats, idents = "Microglia")
mural <- subset(allRats, idents = "Mural")
olig <- subset(allRats, idents = "Olig")
poly <- subset(allRats, idents = "Polydend")

# Syntax for total read count of all genes in subset
sum(Matrix::colSums(D1, slot = "counts"))

# Syntax for read count for a given gene (e.g., Virus)
sum(FetchData(D1, vars = "Virus", slot = "counts"))


############################################################
# Export all barcodes (nonUniqueCellName) for each cluster #
############################################################

# Using subsets created above...
# Create "dictionary" of cluster name and subset object
subset_list <- list(c("astro", astro), c("d1", d1), c("d2_1", d2_1), c("d2_2", d2_2), 
                    c("d3", d3), c("grm8", grm8), c("gaba", gaba), c("pvalb", pvalb), 
                    c("sst", sst), c("micro", micro), c("mural", mural), c("olig", olig), 
                    c("poly", poly))

# Iterate over clusters
for (cluster in subset_list) {
  # Iterate over GEM wells
  for (x in 1:8) {
    # Removes Seurat error & loop exit upon hitting subset attempt with no cells
    if (sum(cluster[[2]]@meta.data[["GEM"]] == x) == 0) {
      # Skip to next GEM iteration
      next
    } else {
      # Get barcodes (nonUniqueCellName) for all cells in this GEM well
      data <- subset(cluster[[2]], subset = (GEM == x))@meta.data[,c("nonUniqueCellName")]
      # Write list of barcodes to CSV file
      outfilename <- paste(cluster[[1]], x, "barcodes.csv", sep = "_")
      write.table(data, file = outfilename, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}

