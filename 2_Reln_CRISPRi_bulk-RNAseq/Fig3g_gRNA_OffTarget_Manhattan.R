#This script identifies overlap between gRNA off targets identified using the 
#online tool CasOFFinder and assess if there is overlap between off targets and
#DEGs identified in DESeq2. Any identified overlaps can be used in building a 
#Manhattan plot 

library(dplyr)
library(tidyr)
library(IRanges)
library(GenomicAlignments)

CRISPRi <- read.csv("/Volumes/Day-Lab$/KLB/Experiments/KLB059_CRISPR_BulkSequencing/csv/DESeq2_CRISPRi.csv")
CRISPRi$padj <- as.numeric(CRISPRi$padj)

#Get chromosome Info for off target manhattan plot
RN7chr <- read.csv("/Volumes/Day-Lab$/KLB/Experiments/KLB059_CRISPR_BulkSequencing/csv/mart_export_RN7_chr_lite.csv")
#rename column for ease of left join
colnames(RN7chr)[1] = "X"
#for whatever reason, the BioMart gene list I downloaded has duplicates, gotta get rid of em
RN7chr <- distinct(RN7chr, X, .keep_all = TRUE)
#we must convert chromosomes to numbers so their identities are not lost when 
#we convert from character to numeric class
RN7chr$Chromosome.scaffold.name[RN7chr$Chromosome.scaffold.name=="MT"]<-"0"
RN7chr$Chromosome.scaffold.name[RN7chr$Chromosome.scaffold.name=="X"]<-"21"
RN7chr$Chromosome.scaffold.name[RN7chr$Chromosome.scaffold.name=="Y"]<-"21"
#now convert to numeric class
#RN7chr$Chromosome.scaffold.name <- as.numeric(RN7chr$Chromosome.scaffold.name)
RN7chr$Chromosome.scaffold.name <- as.integer(RN7chr$Chromosome.scaffold.name)
#some chromosomes are poorly annotated, and do not have a number we will remove them
RN7chr <- na.omit(RN7chr)

#Join CRISPRi gene list with chromosome list
CRISPRi_chr <- left_join(CRISPRi, RN7chr)
CRISPRi_chr <- na.omit(CRISPRi_chr)

#subset by base mean >50
CRISPRi_chr <- subset(CRISPRi_chr, subset = CRISPRi_chr$baseMean > 50)

#read in csv generated from casOFFinder
RelnCasOFF <- read.csv("/Volumes/Day-Lab$/KLB/Experiments/KLB059_CRISPR_BulkSequencing/csv/Reln_CasOFF_Rn7.csv")
#read in refseq accession #
RefSeqChrom <- read.csv("/Volumes/Day-Lab$/KLB/Experiments/KLB059_CRISPR_BulkSequencing/csv/RefSeq_Chromosome.csv")
#Merge to get chromsome number instead of stupid refseq number
RelnCasOFF <- merge(RelnCasOFF, RefSeqChrom, by.x = "Chromosome", by.y = "RefSeq_Accession", all.x = TRUE, all.y = TRUE)
RelnCasOFF <- na.omit(RelnCasOFF)


CRISPRi_chr_sig <- subset(CRISPRi_chr, subset = CRISPRi_chr$padj < 0.05, basemean)

DEG_df <- data.frame(ENSORG = CRISPRi_chr_sig$X,
                     Gene = CRISPRi_chr_sig$name,
                     Chromosome = CRISPRi_chr_sig$Chromosome.scaffold.name,
                     GeneStart = CRISPRi_chr_sig$Gene.start..bp.,
                     GeneEnd = CRISPRi_chr_sig$Gene.end..bp.)

OffTarget_df <- data.frame(SeqString = RelnCasOFF$DNA,
                           Chromosome = RelnCasOFF$Chromosome.y,
                           SeqStart = RelnCasOFF$Position-2000, #range +/- 2kbp
                           SeqEnd = RelnCasOFF$Position+2000,
                           Direction = RelnCasOFF$Direction,
                           BulgeType = RelnCasOFF$X.Bulge.type,
                           BulgeSize = RelnCasOFF$Bulge.Size,
                           Mismatches = RelnCasOFF$Mismatches)

#OT_df_noBulge <- subset(OffTarget_df, subset = OffTarget_df$BulgeSize == 0)

DEG_ranges <- makeGRangesFromDataFrame(DEG_df,keep.extra.columns = TRUE)
OffTarget_ranges <- makeGRangesFromDataFrame(OffTarget_df,keep.extra.columns = TRUE)


table(!is.na(findOverlaps(DEG_ranges, OffTarget_ranges, select="arbitrary")))
overlap_count <- countOverlaps(DEG_ranges, OffTarget_ranges)
findOverlaps(DEG_ranges, OffTarget_ranges)

results_by_OT <- as.data.frame(subsetByOverlaps(OffTarget_ranges, DEG_ranges))
results_by_DEG <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges))

merged_results <- merge(results_by_DEG, CRISPRi_chr, by.x = "Gene", by.y = "name", all.x= TRUE, all.y = FALSE)

write.csv(results_by_DEG, "/Users/kbrida/Desktop/Manuscript/Data&Code/ManhattantPlot/all_results_by_DEG.csv")
write.csv(results_by_OT, "/Users/kbrida/Desktop/Manuscript/Data&Code/ManhattantPlot/all_results_by_OT.csv")

write.csv(merged_results, "/Users/kbrida/Desktop/Manuscript/Data&Code/ManhattantPlot/all_merged_results.csv")

#Want to identify how many mismatches there are for each gene & I am an R infant 
#so apologies super clunky way to ID how many mismatches for off target
OT_df_4mm <- subset(OffTarget_df, subset = OffTarget_df$Mismatches == 4)
OT_df_3mm <- subset(OffTarget_df, subset = OffTarget_df$Mismatches == 3)
OT_df_2mm <- subset(OffTarget_df, subset = OffTarget_df$Mismatches == 2)
OT_df_1mm <- subset(OffTarget_df, subset = OffTarget_df$Mismatches == 1)
OT_df_0mm <- subset(OffTarget_df, subset = OffTarget_df$Mismatches == 0)

OffTarget_ranges4 <- makeGRangesFromDataFrame(OT_df_4mm,keep.extra.columns = TRUE)
OffTarget_ranges3 <- makeGRangesFromDataFrame(OT_df_3mm,keep.extra.columns = TRUE)
OffTarget_ranges2 <- makeGRangesFromDataFrame(OT_df_2mm,keep.extra.columns = TRUE)
OffTarget_ranges1 <- makeGRangesFromDataFrame(OT_df_1mm,keep.extra.columns = TRUE)
OffTarget_ranges0 <- makeGRangesFromDataFrame(OT_df_0mm,keep.extra.columns = TRUE)

results_by_DEG4mm <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges4))
results_by_DEG3mm <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges3))
results_by_DEG2mm <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges2))
results_by_DEG1mm <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges1))
results_by_DEG0mm <- as.data.frame(subsetByOverlaps(DEG_ranges, OffTarget_ranges0))

results_by_DEG1mm$mismatches1 <- 1
results_by_DEG2mm$mismatches2 <- 2
results_by_DEG3mm$mismatches3 <- 3
results_by_DEG4mm$mismatches4 <- 4

mismatch_df <- merge(results_by_DEG1mm, results_by_DEG2mm, by.x ="Gene", by.y = "Gene", all.x = TRUE, all.y = TRUE)
mismatch_df <- merge(mismatch_df, results_by_DEG3mm, by.x ="Gene", by.y = "Gene", all.x = TRUE, all.y = TRUE)
mismatch_df <- merge(mismatch_df, results_by_DEG4mm, by.x ="Gene", by.y = "Gene", all.x = TRUE, all.y = TRUE)

mismatch_df2 <- data.frame(Gene = mismatch_df$Gene,
                           Mismtaches1 = mismatch_df$mismatches1,
                           Mismtaches2 = mismatch_df$mismatches2,
                           Mismtaches3 = mismatch_df$mismatches3,
                           Mismtaches4 = mismatch_df$mismatches4)

mismatch_df2 <- mismatch_df2 %>% replace(is.na(.), 5)
mismatch_df3 <- mismatch_df2 %>% rowwise() %>% mutate(min = min(c_across(!Gene)))
mismatch_tib <- as_tibble(mismatch_df3)
mismatches <- mismatch_df3[-c(2:5)]
colnames(mismatches)[2] <- "Mismatches"

#Now let's make a manhattan plot
#devtools::install_github("boxiangliu/manhattan")
library(manhattan)
#from vingette, requires column "chrom" for chromosome
#also requires column "pos" for gene start & column "y" for the y axis
#we will just make a simplified dataframe
#make chromosome column character class for paste function
CRISPRi_chr$chrom <- as.character(CRISPRi_chr$Chromosome.scaffold.name)
CRISPRi_man <- data.frame(chrom = paste("chr", CRISPRi_chr$chrom,sep=""),
                          pos = CRISPRi_chr$Gene.start..bp.,
                          pval = CRISPRi_chr$padj,
                          gene = CRISPRi_chr$name,
                          y = CRISPRi_chr$log2FoldChange)
#CRISPRi_man$label <- ifelse(CRISPRi_man$pval > 0.05, NA, CRISPRi_man$gene)
CRISPRi_man$Mismatches <- ifelse(CRISPRi_man$gene %in% mismatches$Gene, mismatches$Mismatches, NA)

CRISPRi_man$OffTarget <- ifelse(CRISPRi_man$pval < 0.05 & CRISPRi_man$gene %in% results_by_DEG$Gene == TRUE, "Significant",
                                ifelse(CRISPRi_man$pval > 0.05 & CRISPRi_man$gene %in% results_by_DEG$Gen == TRUE, "NS",
                                        NA))


#did not end up using this 
#man_highlight <- list(c("Dctn2", "Arel1", "Gabrg3", "Ppp4r3a", "Atp6v1e1", "Pkn3", "Prex1",
# "Gse1", "Ube2e3", "Bicd1", "Tle4", "Shroom3", "Nck2", "Reln"))
#man_highlight <- as.vector(man_highlight)


#highlight off-targets and gRNA target
CRISPRi_man$color <- ifelse(CRISPRi_man$Mismatches == 4, "#F5BABD",
                            ifelse(CRISPRi_man$Mismatches == 3, "#FC7382",
                                   ifelse(CRISPRi_man$Mismatches == 2, "#FF1549",
                                          ifelse(CRISPRi_man$Mismatches == 1, "#EA0652",
                                                 NA))))
#generate manhattan
manhattan(CRISPRi_man)
#this plot will likely appear as if nothing has worked because the non-highlighted genes
#far out-number the higlhighted ones. I am an R infant and so I export this plot and open 
#in illustrator. First release all clipping masks. 
#From there select one gray dot and then Select>Same>Appearance
#this will select all gray dots, send them to the back and repeat with black
#your highlighted genes will now be front and center and you can further edit colors
