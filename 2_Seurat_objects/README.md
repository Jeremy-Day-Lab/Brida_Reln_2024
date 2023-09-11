# QC and Seurat Object Creation

The workflow of quality control (QC), Seurat object creation, and dataset integration in the two files here is summarized below:

### Seurat1_SoupX_integrate.r

1. Import data through **SoupX** and run SoupX to remove ambient RNA

2. Create initial **Seurat object** for each dataset

3. Filter out cells containing <200 genes or >5% mitochondrial genes

4. Normalize counts and identify cell-to-cell variable features in each dataset

5. Integrate the 8 datasets into a single integrated Seurat object

6. Dimensionality reduction and clustering

7. Assign cell types to clusters

### Seurat2_DoubletFinder.r

1. Run **DoubletFinder** on each of the 8 datasets

2. Combine metadata (includes DoubletFinder results) from all 8 datasets

3. Remove doublets (i.e., retain only singlets)

4. Re-run dimensionality reduction, clustering, and assigning cell types
