# Code for Brida et al. 2024
This repository contains data analysis code for "Reelin marks cocaine-activated striatal ensembles, promotes neuronal excitability, and regulates cocaine reward"

# Internal Project IDs

* Genomics for 2023_SC0016: 2023-JD-0062

# Contents

## 1. Acute and Repeated Cocaine
Analysis of integrated R object from [Phillips et al. 2023 (Mol. Cell. Neurosci.)](https://github.com/Jeremy-Day-Lab/Phillips_2023_NAc)

1. [Creating Drd1 subset object; Reln & composite IEG expression](1_Acute_and_Repeated_Cocaine/Fig1c_Coc_D1subset_IEG_Reln_expression.R) (Fig. 1c)
2. [Gini coefficient analysis](1_Acute_and_Repeated_Cocaine/Fig1d_GiniCalcs_Drd1SubclustersV2.R) (Fig. 1d)
3. [Composite IEG heatmaps](1_Acute_and_Repeated_Cocaine/Fig1g_CompositeIEGHeatmaps.r) (Fig. 1g)

## 2. Reln CRISPRi Bulk RNA-sequencing

1. [nf-core pipeline files](2_Reln_CRISPRi_bulk-RNAseq/1_nf-core_files)
2. [DEG testing](2_Reln_CRISPRi_bulk-RNAseq/2_DESeq2_CRISPRi_bulk.R)
3. [Manhattan plot](2_Reln_CRISPRi_bulk-RNAseq/Fig3g_gRNA_OffTarget_Manhattan.R) (Fig. 3g)

## 3. Reln CRISPRi single-nucleus RNA-sequencing

1. [CellRanger files](3_Reln_CRISPRi_snRNAseq/1_CellRanger)
2. Creating Seurat object
 1. [SoupX and assay integration](3_Reln_CRISPRi_snRNAseq/2_Seurat_object_creation/1_SoupX_and_integration.r)
 2. [DoubletFinder and cell-type annotation](3_Reln_CRISPRi_snRNAseq/2_Seurat_object_creation/DoubletFinder_and_obj_annotation.r)
 3. [Bar plots](3_Reln_CRISPRi_snRNAseq/2_Seurat_object_creation/3_Fig4c_barplots.r) (Fig. 4c)
3. [DEG testing](3_Reln_CRISPRi_snRNAseq/3_DEGs/README.md)
4. [Pseudotime](3_Reln_CRISPRi_snRNAseq/Fig4h_Pseudotime_Drd1_noReln.R) (Fig. 4h)
5. [Correlation matrix](3_Reln_CRISPRi_snRNAseq/FigS7a_CorrelationMatrix_Combined_MarkerGenesV2.R) (Fig. S7a)
6. [Stacked violin plot](3_Reln_CRISPRi_snRNAseq/FigS7c_StackedViolinByTargetV2.R) (Fig. S7c)
