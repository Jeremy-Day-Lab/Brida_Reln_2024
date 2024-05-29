#---DayLab DESeq2 template script - Jeremy Day ---
# Last update: 9.20.2018

#Based on: NBI DESeq2 template script - Lara Ianov Ph.D.

# Statement: The DESeq2 section of the current script is to serve as a guide to your own analysis, in this case, doing a simple
# comparison between two groups. More complicated designs will need further modifications. Please consult me and the DESeq2 - vignette("DESeq2")

# Package requirements (install these first if you have not already):
#   -DESeq2   (Bioconductor)
#   -pheatmap
#   -ggplot
#   -ggrepel
#   -RColorBrewer
#   -BiomaRt     (Bioconductor)
#   -calibrate
#   -reshape2
#   -tximport  (Bioconductor)

# Usage: To use this script you will need to generate one file: 
#             1) colData.csv
#                 -Provide a "colData" file specifying the design of your study. Current script assumes the grouping information is named "Condition";
#                 change it if you choose another name. The colData file should have two columns. The first column contains the the file names laid out
#                 exactly as they appear in the header of the fixed_raw_counts.txt file. The second column has a header (titled "Condition"), and
#                 contains the group assignments for each sample (e.g., "Vehicle", "Control", "SKF", etc).


############## Data import and DESeq Data Matrix construction ################

# Set working directory to file containing fixed_raw_counts.txt output from the NBI standard pipeline (e.g.: path to Raw_Counts/fixed_raw_counts.txt)
setwd("/Users/JDay/Desktop/2023-JD-0063")

# Make a new directory for this analysis in the working directory - change name as desired
dir.create("./DESeq2")

# Engage DESeq2 library
library(DESeq2)
library(tximport)


### Data prep to go from quant.sf file to DESeq2 object

tx2gene <- read.table("salmon_tx2gene.tsv", sep = '\t',header = FALSE)

head(tx2gene)

# Importing quant.sf file from secondary outputs within data:
# in this case all quant.sf files are in a single directory (and they were rename to contain sample names)

myFiles <- dir("./CRISPRi/", ".sf$", full.names = TRUE)
myFiles

# Adding names for columns:

myFiles_names <- c()

for (i in myFiles) {
  result <- gsub("_quant.sf","",i)
  result <- basename(result)
  myFiles_names[result] <- i
}

all(file.exists(myFiles_names))

# Making a log of the col names from full names:

dir.create("./results/", recursive = TRUE)

Log <- as.matrix(myFiles_names)

write.table(Log, file = "./results/Sample_names_tximport.txt", quote = FALSE, col.names = FALSE, sep = "\t")

### tximport ###

# With Ensembl genome set ignoreTxVersion=TRUE as the GTF does not have versions

txi <- tximport(myFiles_names, type = "salmon", tx2gene = tx2gene, txOut = FALSE, ignoreTxVersion=TRUE)
names(txi)

# NOTE txi$counts is the original with countsFromAbundance="no", DO NOT pass that into DESeq2 directly,
# In short use DESeqDataSetFromTximport in DESeq2 as shown below.
head(txi$counts)
head(txi$abundance)
#head(txi$length)

### Create DESeq2 object ###

colData_CRISPRi <- read.csv("colData_CRISPRi.csv",header=TRUE, row.names=1)


# tximport dds generation

dds <- DESeqDataSetFromTximport(txi,
                                colData = colData_CRISPRi,
                                design = ~ Treatment)


# Switch to new subdirectory
setwd("./DESeq2")

############## Exploring the data by visualization - GOOD QA-QC ################

# PCA plot
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("Treatment")) #to plot a simple PCA


# Heatmap and hierarchical clustering with Euclidean distances
sampleDists <- dist(t(assay(rld))) 
library(RColorBrewer)
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


############## Differential expression ################
# For differential expression with outlier detection as default, 
# where second argument in contrast is the 'control', as this is the denominator.
# Thus all results of the fold change direction will be relative to numerator (treatment).

# DESeq function
dds <- DESeq(dds)

# Dispersion plot
plotDispEsts(dds)


# Results function. Prints summary of differentially expressed genes in Group_1 vs Group_2 comparison.
res <- results(dds, contrast=c("Treatment", "Reln", "lacZ"), alpha = 0.05, pAdjustMethod = "BH")
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

# Order results by padj
resOrdered <- res[order(res$padj),]
write.table(resOrdered, file="DESeq2_CRISPRi.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = NA)



# Save sessionInfo
sink("sessionInfo.txt")
sessionInfo()
sink()
