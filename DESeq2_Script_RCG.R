# Load packages
library("clusterProfiler")
library("dplyr")
library("tidyr")
library("GOplot")
library("enrichplot")
library("fgsea")
library("data.table")
library("GO.db")
library("AnnotationForge")
library("AnnotationHub")
library("writexl")
library("cowplot")
library("ggrepel")
library("EnhancedVolcano")
library("amap")
library("DESeq2") # Differential expression analysis
library("tidyverse") # Misc. data manipulation and plotting
library("here") # Managing file paths
library("pheatmap") # Heatmap plot
library("apeglm") # For LFC shrinkage
library("knitr") # For table printing in Markdown file
#remotes::install_github("YuLab-SMU/clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages("IRanges")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IRanges")

sessionInfo()

theme_set(theme_bw())  # Set ggplot theme


#Set your working directory
setwd("~/R_dir/Pennycress RNASeq/Multi-comparisons")

#Load and format data ----

#input files
count_table_file <- "counts.txt"
metadata_file <- "metadata_file.csv"

#Load input data
raw_counts <- read.table(count_table_file,
                         sep = "\t", header = TRUE, skip = 1)
#Load metadata
metadata <- read_csv(metadata_file, col_types = cols(.default = "c"))

kable(head(metadata))

#Check column names
colnames(raw_counts)[7:8]
colnames(raw_counts)

#remove columns with metadata for each gene
raw_counts[1:5, 1:8]
counts <- raw_counts[, 7:ncol(raw_counts)]
rownames(counts) <- raw_counts$Geneid

#Sort both data frames alphabetically:
metadata <- metadata[order(metadata$SampleID), ]
counts <- counts[, order(colnames(counts))]

#check if the names are the same.. DESeq2 assumes that sample IDs in both tables match and are provided in the same order
metadata$SampleID
colnames(counts)

matching_names <- identical(metadata$SampleID, colnames(counts))
matching_names
if(matching_names == FALSE) stop("Sample ID in metadata and count matrix do not match!")


# Create the DESeq2 object ----
dds_raw <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ 1)

#How many genes have non-zero counts? 
dim(counts[rowSums(counts) > 0, ])

#How many genes have total counts of at least 10?
dim(counts[rowSums(counts) >= 10, ])

# Histogram of gene counts ----
theme_set(theme_bw())

summed_gene_counts <- data.frame(gene_count = rowSums(counts)) %>%
  rownames_to_column("gene_id")

ggplot(data = summed_gene_counts) +
  geom_histogram(aes(x = gene_count), binwidth = 10000) +
  scale_y_log10(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0))
#Zoom in...
ggplot(data = summed_gene_counts) +
  geom_histogram(aes(x = gene_count), binwidth = 1000) +
  scale_y_log10(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 200000), expand = c(0,0)) +
  theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
#How are the counts distributed across samples?
apply(X = counts, MARGIN = 2, FUN = sum)


#Run the PCA and prepare for plotting... 
#First, we normalize the count data to have even sampling across samples (with respect to library size) and approximately even variance
#VST produces transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors
#the point is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low
#for visualization or clustering – it is useful to work with transformed versions of the count data
#When blind equals TRUE (the default), the functions will re-estimate the dispersions using only an intercept. 
#This setting should be used in order to compare samples in a manner wholly unbiased by the information about experimental groups, 
#for example to perform sample QA (quality assurance)
vsd <- varianceStabilizingTransformation(dds_raw, blind = TRUE)


#PCA Plots ----
#Next, we run the PCA and retrieve the data to plot with ggplot2
# * All Conditions/Tissues/Genotypes ----
pcaData <- plotPCA(vsd,
                   ntop = 500,
                   intgroup = c("tissue", "genotype"),
                   returnData = TRUE)
pcaData <- plotPCA(vsd,
                   ntop = 500,
                   intgroup = c("condition", "genotype"),
                   returnData = TRUE)

#We extract the percentage of variance explained by different principal components, so we can later add this information to the plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar
#We create a plot title with the species name in italic using the somewhat bizarre expression() function
plot_title <- expression("PCA of " * italic(Thlaspi ~ arvense) * " transcript count data")

ggplot(pcaData,
       aes(x = PC1, y = PC2, color = tissue, shape = genotype)) +
  geom_point(size = 5) +
  #geom_label_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title) +
  theme(text = element_text(size = 20))




##But we also want to look at PCA of just root tissue and each condition
# * Roots ----
metadata_root <- metadata %>%
  filter(tissue == "roots")
head(metadata_root)

counts_root <- counts %>%
  dplyr::select(MCR1,  MCR2,  MCR3,  MCR4,  MCR5,
                MWR1,  MWR2,  MWR3,  MWR4, MWR5,
                MAR1, MAR2, MAR3, MAR4, MAR5,
                SCR1,   SCR2,   SCR3,   SCR4,   SCR5,
                SWR1,   SWR2,   SWR3,   SWR4,  SWR5,
                SAR1,  SAR2,  SAR3,  SAR4,  SAR5,
                STR1, STR2, STR3, STR4, STR5,
                SCR7, SCR8, SCR9, SCR10, SCR11)
head(counts_root)

metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]

#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ 1)
vsd.root <- varianceStabilizingTransformation(dds_raw_root, blind = TRUE)

pcaData_root <- plotPCA(vsd.root,
                        ntop = 500,
                        intgroup = c("condition", "genotype"),
                        returnData = TRUE)


#We extract the percentage of variance explained by different principal components, so we can later add this information to the plot
percentVar <- round(100 * attr(pcaData_root, "percentVar"))
percentVar
#We create a plot title with the species name in italic using the somewhat bizarre expression() function
plot_title <- expression("PCA of " * italic(Thlaspi ~ arvense) * " transcript count data in root tissue")

ggplot(pcaData_root,
       aes(x = PC1, y = PC2, color = condition, shape = genotype)) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title) +
  theme(text = element_text(size = 20))

#Also going to do it by genotype
#MN106
metadata_root <- metadata %>%
  filter(tissue == "roots", genotype == "MN106")
head(metadata_root)

counts_root <- counts %>%
  dplyr::select(MCR1,  MCR2,  MCR3,  MCR4,  MCR5,
                MWR1,  MWR2,  MWR3,  MWR4, MWR5,
                MAR1, MAR2, MAR3, MAR4, MAR5)

metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]

#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ 1)
vsd.root <- varianceStabilizingTransformation(dds_raw_root, blind = TRUE)

pcaData_root <- plotPCA(vsd.root,
                        ntop = 500,
                        intgroup = c("condition", "genotype"),
                        returnData = TRUE)


#We extract the percentage of variance explained by different principal components, so we can later add this information to the plot
percentVar <- round(100 * attr(pcaData_root, "percentVar"))
percentVar
#We create a plot title with the species name in italic using the somewhat bizarre expression() function
plot_title <- expression("PCA of " * italic(Thlaspi ~ arvense) * " transcript count data in root tissue")

ggplot(pcaData_root,
       aes(x = PC1, y = PC2, color = condition, shape = genotype)) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title) +
  theme(text = element_text(size = 20))

#SP32-10
metadata_root <- metadata %>%
  filter(tissue == "roots", genotype == "SP3210", condition == "control_early")
head(metadata_root)

counts_root <- counts %>%
  dplyr::select(#SCR1,   SCR2,   SCR3,   SCR4,   SCR5,
    #SWR1,   SWR2,   SWR3,   SWR4,  SWR5)
    #SAR1,  SAR2,  SAR3,  SAR4,  SAR5,
    #STR1, STR2, STR3, STR4, STR5)
    SCR7, SCR8, SCR9, SCR10, SCR11)
head(counts_root)

metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]

#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ 1)
vsd.root <- varianceStabilizingTransformation(dds_raw_root, blind = TRUE)

pcaData_root <- plotPCA(vsd.root,
                        ntop = 500,
                        intgroup = c("condition", "genotype"),
                        returnData = TRUE)


#We extract the percentage of variance explained by different principal components, so we can later add this information to the plot
percentVar <- round(100 * attr(pcaData_root, "percentVar"))
percentVar
#We create a plot title with the species name in italic using the somewhat bizarre expression() function
plot_title <- expression("PCA of " * italic(Thlaspi ~ arvense) * " transcript count data in root tissue")

ggplot(pcaData_root,
       aes(x = PC1, y = PC2, color = condition, shape = genotype)) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title) +
  theme(text = element_text(size = 20))



# A correlation heatmap of all samples ----
metadata_root <- metadata %>%
  filter(tissue == "roots" & condition == "control" |
           tissue == "roots" & condition == "waterlogged" |
           tissue == "roots" & condition == "recovery")
#remove sample MCR1 because it failed QC
metadata_root <- metadata_root[-6,]

counts_root <- counts %>%
  dplyr::select(MCR2,  MCR3,  MCR4,  MCR5,
                MWR1,  MWR2,  MWR3,  MWR4, MWR5,
                MAR1, MAR2, MAR3, MAR4, MAR5,
                SCR1,   SCR2,   SCR3,   SCR4,   SCR5,
                SWR1,   SWR2,   SWR3,   SWR4,  SWR5,
                SAR1,  SAR2,  SAR3,  SAR4,  SAR5)
head(counts_root)

metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]

#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ 1)

#Make a factor for condition (control and waterlogged)
dds_raw_root$condition <- factor(paste(dds_raw_root$condition))
#Make a factor for genotype (MN106 and SP32-10)
dds_raw_root$genotype <- factor(paste(dds_raw_root$genotype))

#check reference level - whichever one appears first is set as the reference
dds_raw_root$condition
dds_raw_root$genotype

#Set the design and run DESeq
design(dds_raw_root) <- ~ genotype + genotype:condition
dds <- DESeq(dds_raw_root)

normalized_counts <- counts(dds, normalized=TRUE)
#head(normalized_counts)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

write.table(normalized_counts, file="allsample.gene.counts.matrix.filter.DESeq2.normalized.root.xls",
            quote=F, sep="\t", row.names=T, col.names=T)

## log transfomation
rld <- rlog(dds, blind=FALSE)
rld <- vst(dds, blind=FALSE,fitType = "mean")
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(100)

# pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
write.table(pearson_cor, file="allsample.gene.counts.matrix.filter.DESeq2.normalized.pearson.root.txt",
            quote=F, sep="\t", row.names=T, col.names=T)
head(pearson_cor)

# hierachical clustering
hc <- hcluster(t(rlogMat), method="pearson")
library(gplots)
pdf("PearsonCorrelationRoots.pdf", pointsize=10)
jpeg(filename= "PearsonCorrelationRoots2.jpeg", width = 1000, height = 1000, pointsize=5, res=500, quality=95)
png(filename= "PearsonCorrelationRoots2.png", width = 1500, height = 1500, pointsize=5, res=500)
par(mar = c(0, 10, 0, 0), xpd = TRUE)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="Pearson correlation of root samples",
          ColSideColors = c(    # grouping row-variables into different
            rep("#009E73", 5),   # MN106 Recovery
            rep("#FFCC99", 4),    # MN106 Controls
            rep("#009E73", 5),     # MN106 Waterlogged
            rep("#56B4E9", 5),       # SP32-10 Recovery
            rep("#CC79A7", 5),         # SP32-10 Controls
            rep("#56B4E9", 5)),         #SP32-10 Waterlogged
          
          RowSideColors = c(    # grouping row-variables into different
            rep("#009E73", 5),   # MN106 Recovery
            rep("#FFCC99", 4),    # MN106 Controls
            rep("#009E73", 5),     # MN106 Waterlogged
            rep("#56B4E9", 5),       # SP32-10 Recovery
            rep("#CC79A7", 5),         # SP32-10 Controls
            rep("#56B4E9", 5)          #SP32-10 Waterlogged
          )         
          
)
legend("bottomright", inset = c(0, 0), cex=0.8, title = "Condition/Genotype",legend=c("MN106 Waterlogged/Recovery", "SP32-10 Waterlogged/Recovery", "MN106 Controls", "SP32-10 Controls"), 
       fill=c("#009E73","#56B4E9", "#FFCC99", "#CC79A7"))
dev.off()




#Now we can do the analysis - we'll do each tissue separately to keep things simple
# DESeq2 Analysis ----
#Chunk of code from above to filter out root samples - but we need to keep out the extra SP32 time points for now
#So we filter the dataframe by leaving out the early waterlogging time points
metadata_root <- metadata %>%
  filter(tissue == "roots" & condition == "control" |
           tissue == "roots" & condition == "waterlogged" |
           tissue == "roots" & condition == "recovery")

#remove sample MCR1 because it failed QC
metadata_root <- metadata_root[-6,]

counts_root <- counts %>%
  dplyr::select(MCR2,  MCR3,  MCR4,  MCR5,
                MWR1,  MWR2,  MWR3,  MWR4, MWR5,
                MAR1, MAR2, MAR3, MAR4, MAR5,
                SCR1,   SCR2,   SCR3,   SCR4,   SCR5,
                SWR1,   SWR2,   SWR3,   SWR4,  SWR5,
                SAR1,  SAR2,  SAR3,  SAR4,  SAR5)
head(counts_root)


metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]


#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ 1)

#Make a factor for condition (control and waterlogged)
dds_raw_root$condition <- factor(paste(dds_raw_root$condition))
#Make a factor for genotype (MN106 and SP32-10)
dds_raw_root$genotype <- factor(paste(dds_raw_root$genotype))


#check reference level - whichever one appears first is set as the reference
dds_raw_root$condition
dds_raw_root$genotype

#change the reference level if needed
#dds_raw_root$genotype = relevel(dds_raw_root$genotype, "SP3210")
#dds_raw_root$condition = relevel(dds_raw_root$condition, "waterlogged")

#Set the design and run DESeq
design(dds_raw_root) <- ~ genotype + genotype:condition
dds <- DESeq(dds_raw_root)


#Check the names of the comparisons/interactions
resultsNames(dds)


#if calling the results of dds, the default is the last level of the factor with the first level
#so to avoid this,
#here are all of our comparisons assigned to new, separate variables with FDR < 0.05 (FDR = alpha) and LFC=1
#We'll assign ALL of the DEGs to a separate variable - useful for downstream analyses like GSEA
# DEG Results ----

#MN106 waterlogged vs control
results.WvC.MN106 <- results(dds, name="genotypeMN106.conditionwaterlogged", alpha=0.05, lfcThreshold=1)
summary(results.WvC.MN106)
#SP32-10 waterlogged vs control
results.WvC.SP32 <- results(dds, name="genotypeSP3210.conditionwaterlogged", alpha=0.05, lfcThreshold=1)
summary(results.WvC.SP32)
#MN106 vs SP32-10 Waterlogged (so genes are up or down regulated in MN106 compared to SP32-10)
results.WvW <- results(dds, contrast=list("genotypeMN106.conditionwaterlogged","genotypeSP3210.conditionwaterlogged"), alpha=0.05, lfcThreshold=1)
head(results.WvW, 10)
summary(results.WvW)
results <- results(dds, alpha=0.05, lfcThreshold=1)
head(results, 10)
summary(results)
#Or we can reverse it and do SP32-10 vs MN106 (genes are up or down regulated in SP32-10 compared to MN106)
results.WvW2 <- results(dds, contrast=list("genotypeSP3210.conditionwaterlogged","genotypeMN106.conditionwaterlogged"), alpha=0.05, lfcThreshold=1)
summary(results.WvW2)
#MN106 Recovery vs Control
results.RvC.MN106 <- results(dds, name="genotypeMN106.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results.RvC.MN106)
#SP32-10 Recovery vs Control
results.RvC.SP32 <- results(dds, name="genotypeSP3210.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results.RvC.SP32)
#MN106 vs SP32-10 Recovery (with controls as reference in each)
results.RvR <- results(dds, contrast=list("genotypeMN106.conditionrecovery","genotypeSP3210.conditionrecovery"), alpha=0.05, lfcThreshold=1)
summary(results.RvR)
#MN106 Recovery vs Waterlogged
results.RvW.MN106 <- results(dds, contrast=list("genotypeMN106.conditionrecovery", "genotypeMN106.conditionwaterlogged"), alpha=0.05, lfcThreshold=1)
summary(results.RvW.MN106)
#SP32-10 Recovery vs Waterlogged
results.RvW.SP3210 <- results(dds, contrast=list("genotypeSP3210.conditionrecovery", "genotypeSP3210.conditionwaterlogged"), alpha=0.05, lfcThreshold=1)
summary(results.RvW.SP3210)
#SP3210 vs MN106 Controls
#If we want to compare the controls we have to relevel or redesign
metadata_root <- metadata %>%
  filter(tissue == "roots" & condition == "control" )
#remove sample MCR1 because it failed QC
metadata_root <- metadata_root[-6,]
head(metadata_root)
counts_root <- counts %>%
  dplyr::select(MCR2,  MCR3,  MCR4,  MCR5,
                SCR1,   SCR2,   SCR3,   SCR4,   SCR5)
head(counts_root)
metadata_root <- metadata_root[order(metadata_root$SampleID), ]
counts_root <- counts_root[, order(colnames(counts_root))]
#create the DESeq2 object
dds_raw_root <- DESeqDataSetFromMatrix(countData = counts_root,
                                       colData = metadata_root,
                                       design = ~ genotype)
design(dds_raw_root) <- ~ genotype
dds_2 <- DESeq(dds_raw_root) 
resultsNames(dds_2)
results.CvC <- results(dds_2, name="genotype_SP3210_vs_MN106", alpha=0.05, lfcThreshold=1)
summary(results.CvC)


#MN106 vs SP32-10 Recovery (with waterlogging as reference in each)
#check reference level - whichever one appears first is set as the reference
dds_raw_root$condition
dds_raw_root$genotype
#change the reference level
dds_raw_root$condition = relevel(dds_raw_root$condition, "waterlogged")
#Set the design and run DESeq
design(dds_raw_root) <- ~ genotype + genotype:condition
dds <- DESeq(dds_raw_root)
resultsNames(dds)
results.RvR.2 <- results(dds, contrast=list("genotypeMN106.conditionrecovery","genotypeSP3210.conditionrecovery"), alpha=0.05, lfcThreshold=1)
summary(results.RvR.2)



plotDispEsts(dds)

# LFC estimates ----
#We’ll also create an object with adjusted (shrunken) LFC estimates, which will be useful for visualization and ranking of genes by LFC
resLFC.WvC.MN106 <- lfcShrink(dds, coef="genotypeMN106.conditionwaterlogged", type="apeglm", lfcThreshold = 1)
resLFC.WvC.SP32 <- lfcShrink(dds, coef="genotypeSP3210.conditionwaterlogged", type="apeglm", lfcThreshold = 1)
resLFC.RvC.MN106 <- lfcShrink(dds, coef="genotypeMN106.conditionrecovery", type="apeglm", lfcThreshold = 1)
resLFC.RvC.SP32 <- lfcShrink(dds, coef="genotypeSP3210.conditionrecovery", type="apeglm", lfcThreshold = 1)
resLFC.CvC <- lfcShrink(dds_2, coef="genotype_SP3210_vs_MN106", type="apeglm", lfcThreshold = 1)

#We need a new design for comparing accessions with controls as reference so that there is only one coefficent
design(dds_raw_root) <- ~ genotype + condition + genotype:condition
dds_raw_root$condition = relevel(dds_raw_root$condition, "control")
dds_raw_root$genotype = relevel(dds_raw_root$genotype, "SP3210")
dds_LFC <- DESeq(dds_raw_root)
resultsNames(dds_LFC)
#Check it first
results_test <- results(dds_LFC, name="genotypeMN106.conditionwaterlogged", alpha=0.05, lfcThreshold=1)
head(results_test, 10)
summary(results_test)
results_test2 <- results(dds_LFC, name="genotypeMN106.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results_test2)
results.RvW.accessions <- results(dds_LFC, contrast=list("genotypeMN106.conditionrecovery", "genotypeMN106.conditionwaterlogged"), alpha=0.05, lfcThreshold=1)
summary(results.RvW.accessions)
#Looks correct, now let's do the lfcshrink
resLFC.WvW <- lfcShrink(dds_LFC, coef="genotypeMN106.conditionwaterlogged", type="apeglm", lfcThreshold = 1)
resLFC.RvR <- lfcShrink(dds_LFC, coef="genotypeMN106.conditionrecovery", type="apeglm", lfcThreshold = 1)

#For these two comparisons we can just relevel so that the reference is waterlogged
design(dds_raw_root) <- ~ genotype + genotype:condition
#dds_raw_root$genotype = relevel(dds_raw_root$genotype, "MN106")
dds_raw_root$condition = relevel(dds_raw_root$condition, "waterlogged")
dds_LFC2 <- DESeq(dds_raw_root)
resultsNames(dds_LFC2)
results_test3 <- results(dds_LFC2, name="genotypeSP3210.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results_test3)
results_test3df <- as.data.frame(results_test3)
resLFC.RvW.SP3210 <- lfcShrink(dds_LFC2, coef="genotypeSP3210.conditionrecovery", type="apeglm", lfcThreshold = 1)
results_test4 <- results(dds_LFC2, name="genotypeMN106.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results_test4)
results_test4df <- as.data.frame(results_test4)
resLFC.RvW.MN106 <- lfcShrink(dds_LFC2, coef="genotypeMN106.conditionrecovery", type="apeglm", lfcThreshold = 1)
#old code... But now we need to relevel the genotype too for MN106...
dds_raw_root$genotype = relevel(dds_raw_root$genotype, "SP3210")
dds_LFC3 <- DESeq(dds_raw_root)
resultsNames(dds_LFC3)
results_test4 <- results(dds_LFC3, name="genotypeMN106.conditionrecovery", alpha=0.05, lfcThreshold=1)
summary(results_test4)
resLFC.RvW.MN106 <- lfcShrink(dds_LFC3, coef="genotypeMN106.conditionrecovery", type="apeglm", lfcThreshold = 1)



# DEG Results - filtered ----

#DESeq2 also automatically filters genes with a low mean count in the sense that it does not include them in the multiple testing correction. 
#Therefore, in such cases, the p-value will not be NA, but the adjusted p-value will be.
#Because we have very low power to detect differential expression for such low-count genes, it is beneficial to remove them prior to the multiple testing correction: 
#that way, the correction becomes less severe for the remaining genes.

#change p-value as NA to 1
results.WvC.MN106$padj[is.na(results.WvC.MN106$padj)] <- 1
results.WvC.SP32$padj[is.na(results.WvC.SP32$padj)] <- 1
results.WvW$padj[is.na(results.WvW$padj)] <- 1
results.WvW2$padj[is.na(results.WvW2$padj)] <- 1
results.RvC.MN106$padj[is.na(results.RvC.MN106$padj)] <- 1
results.RvC.SP32$padj[is.na(results.RvC.SP32$padj)] <- 1
results.RvR$padj[is.na(results.RvR$padj)] <- 1
results.RvW.MN106$padj[is.na(results.RvW.MN106$padj)] <- 1
results.RvW.SP3210$padj[is.na(results.RvW.SP3210$padj)] <- 1
results.CvC$padj[is.na(results.CvC$padj)] <- 1
results.RvW.accessions$padj[is.na(results.RvW.accessions$padj)] <- 1

#order the results by adjusted pvalue and export the sorted DEG results (with NA padj values changed to 1) to a csv file
#Note- specifying alpha = 0.05 above does not exclude that data from the data.frame, so we need to subset here
#Will include a file for all ~26,000 genes and a file for the filtered DEGs based on FDR = .05 and LFC=1
results.WvC.MN106 <- results.WvC.MN106[order(results.WvC.MN106$padj),]
head(results.WvC.MN106)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Waterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.WvC.MN106), file=file, row.names=T)
results.WvC.MN106_padj <- subset(results.WvC.MN106, results.WvC.MN106$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Waterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.WvC.MN106_padj), file=file, row.names=T)
results.WvC.MN106_top20 <- results.WvC.MN106[1:20,]

results.WvC.SP32 <- results.WvC.SP32[order(results.WvC.SP32$padj),]
head(results.WvC.SP32)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.WvC.SP32), file=file,  row.names=T)
results.WvC.SP32_padj <- subset(results.WvC.SP32, results.WvC.SP32$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.WvC.SP32_padj), file=file,  row.names=T)

results.WvW <- results.WvW[order(results.WvW$padj),]
head(results.WvW)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Waterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.WvW), file=file,  row.names=T)
results.WvW_padj <- subset(results.WvW, results.WvW$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Waterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.WvW_padj ), file=file,  row.names=T)

results.RvC.MN106 <- results.RvC.MN106[order(results.RvC.MN106$padj),]
head(results.RvC.MN106)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Recovery_roots_ALL.csv"
write.csv(as.data.frame(results.RvC.MN106), file=file,  row.names=T)
results.RvC.MN106_padj <- subset(results.RvC.MN106, results.RvC.MN106$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Recovery_roots_filtered.csv"
write.csv(as.data.frame(results.RvC.MN106_padj), file=file,  row.names=T)

results.RvC.SP32 <- results.RvC.SP32[order(results.RvC.SP32$padj),]
head(results.RvC.SP32)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Recovery_roots_ALL.csv"
write.csv(as.data.frame(results.RvC.SP32), file=file,  row.names=T)
results.RvC.SP32_padj <- subset(results.RvC.SP32, results.RvC.SP32$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Recovery_roots_filtered.csv"
write.csv(as.data.frame(results.RvC.SP32_padj), file=file,  row.names=T)

results.RvR <- results.RvR[order(results.RvR$padj),]
head(results.RvR)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Recovery_roots_ALL.csv"
write.csv(as.data.frame(results.RvR), file=file,  row.names=T)
results.RvR_padj <- subset(results.RvR, results.RvR$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Recovery_roots_filtered.csv"
write.csv(as.data.frame(results.RvR_padj), file=file,  row.names=T)

results.RvW.MN106 <- results.RvW.MN106[order(results.RvW.MN106$padj),]
head(results.RvW.MN106)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_RecoveryvsWaterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.RvW.MN106), file=file,  row.names=T)
results.RvW.MN106_padj <- subset(results.RvW.MN106, results.RvW.MN106$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_RecoveryvsWaterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.RvW.MN106_padj), file=file,  row.names=T)

results.RvW.SP3210 <- results.RvW.SP3210[order(results.RvW.SP3210$padj),]
head(results.RvW.SP3210)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_RecoveryvsWaterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.RvW.SP3210), file=file,  row.names=T)
results.RvW.SP3210_padj <- subset(results.RvW.SP3210, results.RvW.SP3210$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_RecoveryvsWaterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.RvW.SP3210_padj), file=file,  row.names=T)

results.CvC <- results.CvC[order(results.CvC$padj),]
head(results.CvC)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/CvC_roots_ALL.csv"
write.csv(as.data.frame(results.CvC), file=file,  row.names=T)
results.CvC_padj <- subset(results.CvC, results.CvC$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/CvC_roots_filtered.csv"
write.csv(as.data.frame(results.CvC_padj), file=file,  row.names=T)

results.RvW.accessions <- results.RvW.accessions[order(results.RvW.accessions$padj),]
head(results.RvW.accessions)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/accessions_RecoveryvsWaterlogged_roots_ALL.csv"
write.csv(as.data.frame(results.RvW.accessions), file=file,  row.names=T)
results.RvW.accessions_padj <- subset(results.RvW.accessions, results.RvW.accessions$padj<=0.05)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/accessions_RecoveryvsWaterlogged_roots_filtered.csv"
write.csv(as.data.frame(results.RvW.accessions_padj), file=file,  row.names=T)

#UP-regulated ----

#Now we want to export separate files for upregulated genes and downregulated genes - this would be useful for ORA
# foldchange > 1 for upregulated genes
# log2Foldchange >=1 is the same as foldchange >=2
res_de_up_WvC.MN106_padj <- subset(results.WvC.MN106_padj, results.WvC.MN106_padj$log2FoldChange>=1)
summary(res_de_up_WvC.MN106_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Waterlogged_roots_UP.csv"
write.csv(as.data.frame(res_de_up_WvC.MN106_padj), file=file,  row.names=T)
res_de_up_WvC.MN106_top25 <- res_de_up_WvC.MN106_padj[1:25,]

res_de_up_WvC.SP32_padj <- subset(results.WvC.SP32_padj, results.WvC.SP32_padj$log2FoldChange>=1)
summary(res_de_up_WvC.SP32_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_UP.csv"
write.csv(as.data.frame(res_de_up_WvC.SP32_padj), file=file,  row.names=T)
res_de_up_WvC.SP32_top25 <- res_de_up_WvC.SP32_padj[1:25,]

res_de_up_WvW_padj <- subset(results.WvW_padj, results.WvW_padj$log2FoldChange>=1)
summary(res_de_up_WvW_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Waterlogged_roots_UP.csv"
write.csv(as.data.frame(res_de_up_WvW_padj ), file=file,  row.names=T)

res_de_up_RvC.MN106_padj <- subset(results.RvC.MN106_padj , results.RvC.MN106_padj$log2FoldChange>=1)
summary(res_de_up_RvC.MN106_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Recovery_roots_UP.csv"
write.csv(as.data.frame(res_de_up_RvC.MN106_padj ), file=file,  row.names=T)

res_de_up_RvC.SP32_padj <- subset(results.RvC.SP32_padj , results.RvC.SP32_padj$log2FoldChange>=1)
summary(res_de_up_RvC.SP32_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Recovery_roots_UP.csv"
write.csv(as.data.frame(res_de_up_RvC.SP32_padj ), file=file,  row.names=T)

res_de_up_RvR_padj <- subset(results.RvR_padj , results.RvR_padj$log2FoldChange>=1)
summary(res_de_up_RvR_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Recovery_roots_UP.csv"
write.csv(as.data.frame(res_de_up_RvR_padj ), file=file,  row.names=T)

res_de_up_RvW.MN106_padj <- subset(results.RvW.MN106_padj , results.RvW.MN106_padj$log2FoldChange>=1)
summary(res_de_up_RvW.MN106_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_RecoveryvsWaterlogged_roots_UP.csv"
write.csv(as.data.frame(res_de_up_RvW.MN106_padj), file=file,  row.names=T)

res_de_up_RvW.SP3210_padj <- subset(results.RvW.SP3210_padj, results.RvW.SP3210_padj$log2FoldChange>=1)
summary(res_de_up_RvW.SP3210_padj )
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_RecoveryvsWaterlogged_roots_UP.csv"
write.csv(as.data.frame(res_de_up_RvW.SP3210_padj), file=file,  row.names=T)


# DOWN-regulated ----
# foldchange > 1 for downregulated genes
res_de_dw_WvC.MN106_padj <- subset(results.WvC.MN106_padj , results.WvC.MN106_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_WvC.MN106_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Waterlogged_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_WvC.MN106_padj), file=file,  row.names=T)
res_de_dw_WvC.MN106_top25 <- res_de_dw_WvC.MN106_padj[1:25,]

res_de_dw_WvC.SP32_padj  <- subset(results.WvC.SP32_padj , results.WvC.SP32_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_WvC.SP32_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_WvC.SP32_padj), file=file,  row.names=T)
res_de_dw_WvC.SP32_top25 <- res_de_dw_WvC.SP32_padj[1:25,]

res_de_dw_WvW_padj <- subset(results.WvW_padj, results.WvW_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_WvW_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Waterlogged_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_WvW_padj), file=file,  row.names=T)

res_de_dw_RvC.MN106_padj<- subset(results.RvC.MN106_padj, results.RvC.MN106_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_RvC.MN106_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Recovery_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_RvC.MN106_padj), file=file,  row.names=T)

res_de_dw_RvC.SP32_padj<- subset(results.RvC.SP32_padj, results.RvC.SP32_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_RvC.SP32_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Recovery_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_RvC.SP32_padj), file=file,  row.names=T)

res_de_dw_RvR_padj<- subset(results.RvR_padj, results.RvR_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_RvR_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106vsSP32_Recovery_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_RvR_padj), file=file,  row.names=T)

res_de_dw_RvW.MN106_padj<- subset(results.RvW.MN106_padj, results.RvW.MN106_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_RvW.MN106_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_RecoveryvsWaterlogged_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_RvW.MN106_padj), file=file,  row.names=T)

res_de_dw_RvW.SP3210_padj<- subset(results.RvW.SP3210_padj, results.RvW.SP3210_padj$log2FoldChange<=(-1)*1)
summary(res_de_dw_RvW.SP3210_padj)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_RecoveryvsWaterlogged_roots_DOWN.csv"
write.csv(as.data.frame(res_de_dw_RvW.SP3210_padj), file=file,  row.names=T)




#Now we are going to visualize the results.
#We are first going to transform the data - rlog and vst and ntd
#We'll move foward with rlog
#Then we can plot the transformed data as a heatmap and to look at sample-to-sample distances
#This is probably better to run before down-stream analysis (along with the PCA plots)
#After this, we can start down-stream visualization of data
# rlog transform count data ----

rld <- rlog(dds, blind=FALSE)
head( assay(rld) )
#To show the effect of the transformation, we plot the first sample against the second, first simply using
#the log2 function (after adding 1, to avoid taking the log of zero), and then using the rlog-transformed values.
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
#we can see how genes with low counts seem to be excessively variable on the ordinary logarithmic scale, 
#while the rlog transform compresses differences for genes for which the data cannot provide good information anyway

#We can also compare to vsd transformation
vsd <- vst(dds, blind=FALSE)
plot( assay(vsd)[, 1:2], col="#00000020", pch=20, cex=0.3 )


sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition,
                                     rld$genotype, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(4, 10))

library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


#transform count data on the log2 scale normalized with respect to library size and other normalization factors
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(rld)
library("vsn")
meanSdPlot(assay(rld))

#To explore a count matrix, it is often instructive to look at it as a heatmap. 
#Below is how to produce a heatmap for various transformations of the data.
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","genotype")])
#ntd method
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#vsd method
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
#rld method
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[])


#Another use of the transformed data is sample clustering. 
#Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Now let's visualize our data
#MA-plot
plotMA(resLFC.WvC.MN106, ylim=c(-10,10))
#Make the plot interactive
idx <- identify(results.WvC.MN106$baseMean, results.WvC.MN106$log2FoldChange)
rownames(results.WvC.MN106)[idx]

# Dispersion plot ----
plotDispEsts(dds)

#histogram of p-values
hist(results.WvC.MN106$pvalue, breaks=20, col="grey" )



# Volcano plots ----
#we need to use the raw results without filtering
#we can use the max and min functions to figure out the xlim
#we need to manually enter the p-value!
#results.WvC.MN106 <- results(dds, name="genotypeMN106.conditionwaterlogged")
#max(results.WvC.MN106$log2FoldChange, na.rm=T)
#min(results.WvC.MN106$log2FoldChange, na.rm=T)
pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvC.MN106.pdf",wi = 8, he = 8)
print(EnhancedVolcano(results.WvC.MN106,
                      lab = rownames(results.WvC.MN106),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-12, 14),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = 'MN106 Waterlogged vs Control Roots',
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.008868,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()

EnhancedVolcano(results.WvC.MN106,
                lab = rownames(results.WvC.MN106),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-12, 14),
                #selectlab = names,
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = 'MN106 Waterlogged vs Control Roots',
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.008868,
                #pCutoffCol = 'padj',
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvC.MN106.png", device="png",
       width = 8,
       height = 7,)



pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvC.SP3210.pdf",wi = 8, he = 8)
print(EnhancedVolcano(results.WvC.SP32,
                      lab = rownames(results.WvC.SP32),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-9.5, 14),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = 'SP32-10 Waterlogged vs Control Roots',
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.004486,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.WvC.SP32,
                lab = rownames(results.WvC.SP32),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-9.5, 14),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = 'SP32-10 Waterlogged vs Control Roots',
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.004486,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvC.SP3210.png", device="png",
       width = 8,
       height = 7,)



#max(results.WvW$log2FoldChange, na.rm=T)
#min(results.WvW$log2FoldChange, na.rm=T)
pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvW.pdf",wi = 8, he = 8)
my_title <- 'Up and Down-regulated genes in MN106 vs SP32-10 Waterlogged Roots'
print(EnhancedVolcano(results.WvW,
                      lab = rownames(results.WvW),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-10.3, 10.62),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = str_wrap(my_title, 50),
                      #subtitle = 'Down-regulated genes are overexpressed in SP32-10 vs MN106 waterlogged Roots',
                      pCutoff = .000348,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()

EnhancedVolcano(results.WvW,
                lab = rownames(results.WvW),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-10.3, 10.62),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = str_wrap(my_title, 50),
                #subtitle = 'Down-regulated genes are overexpressed in SP32-10 vs MN106 waterlogged Roots',
                #pCutoffCol = 'padj',
                pCutoff = .000348,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.WvW.png", device="png",
       width = 8,
       height = 7,)



pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvC.MN106.pdf",wi = 8, he = 8)
print(EnhancedVolcano(results.RvC.MN106,
                      lab = rownames(results.RvC.MN106),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-8.33, 12.4),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = 'MN106 Recovery vs Control Roots',
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.006043,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.RvC.MN106,
                lab = rownames(results.RvC.MN106),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-8.33, 12.4),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = 'MN106 Recovery vs Control Roots',
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.006043,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvC.MN106.png", device="png",
       width = 8,
       height = 7,)


pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvC.SP32.pdf",wi = 8, he = 8)
print(EnhancedVolcano(results.RvC.SP32,
                      lab = rownames(results.RvC.SP32),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-11.9, 13),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = 'SP32-10 Recovery vs Control Roots',
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.00379,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.RvC.SP32,
                lab = rownames(results.RvC.SP32),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-11.9, 13),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = 'SP32-10 Recovery vs Control Roots',
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.00379,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvC.SP32.png", device="png",
       width = 8,
       height = 7,)



pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvR.pdf",wi = 8, he = 8)
my_title2 <- 'Up and Down-regulated genes in MN106 vs SP32-10 Recovery Roots'
print(EnhancedVolcano(results.RvR,
                      lab = rownames(results.RvR),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-9.4, 8.5),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = str_wrap(my_title2, 50),
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.00017,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.RvR,
                lab = rownames(results.RvR),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-9.4, 8.5),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = str_wrap(my_title2, 50),
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.00017,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvR.png", device="png",
       width = 8,
       height = 7,)      



pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvW.MN106.pdf",wi = 8, he = 8)
my_title3 <- 'Up and Down-regulated genes in MN106 Recovery vs Waterlogged Roots'
print(EnhancedVolcano(results.RvW.MN106,
                      lab = rownames(results.RvW.MN106),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-4.9, 7),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = str_wrap(my_title3, 50),
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 1.9e-6,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.RvW.MN106,
                lab = rownames(results.RvW.MN106),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-4.9, 7),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = str_wrap(my_title3, 50),
                #title = 'SP32-10 Roots - Recovery vs Control',
                #pCutoffCol = 'padj',
                pCutoff = 1.9e-6,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvW.MN106.png", device="png",
       width = 8,
       height = 7,)   


pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvW.SP3210.pdf",wi = 8, he = 8)
my_title4 <- 'Up and Down-regulated genes in SP32-10 Recovery vs Waterlogged Roots'
print(EnhancedVolcano(results.RvW.SP3210,
                      lab = rownames(results.RvW.SP3210),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-4.9, 5.9),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = str_wrap(my_title4, 50),
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 10e-5,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.RvW.SP3210,
                lab = rownames(results.RvW.SP3210),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-4.9, 5.9),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = str_wrap(my_title4, 50),
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 10e-5,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.RvW.SP3210.png", device="png",
       width = 8,
       height = 7,) 


pdf("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.CvC.pdf",wi = 8, he = 8)
my_title5 <- 'Up and Down-regulated genes in MN106 Control vs SP32-10 Control Roots'
print(EnhancedVolcano(results.CvC,
                      lab = rownames(results.CvC),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      #xlim = c(-4.9, 5.9),
                      #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                      title = str_wrap(my_title5, 50),
                      #title = 'SP32-10 Roots - Recovery vs Control',
                      pCutoff = 0.000298388,
                      FCcutoff = 2,
                      selectLab=""
))
dev.off()
EnhancedVolcano(results.CvC,
                lab = rownames(results.CvC),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-4.9, 5.9),
                #title = 'SP32-10 Roots - 3 Days Waterlogged vs Control',
                title = str_wrap(my_title5, 50),
                #title = 'SP32-10 Roots - Recovery vs Control',
                pCutoff = 0.000298388,
                FCcutoff = 2,
                selectLab="")
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/Volcano Plots/VAplot.CvC.png", device="png",
       width = 8,
       height = 7,)   





#Heatmaps ---- 
#need to use normalized values (rlog)

#rld is in large DESEQTransform format so we need to get into a more simplified format
rld_assay <- assay(rld)
head(rld_assay)

#Select top 20 genes and create as list
orthologs2 <- read.csv("Orthologs2.csv", header= TRUE)
head(orthologs2)
orthologs <- read.csv("Orthologs_FINAL_no-duplicates.csv", header=TRUE)


results.WvC.MN106_top20 <- results.WvC.MN106[1:20,]
results.WvC.MN106_top20 <- setNames(cbind(rownames(results.WvC.MN106_top20), results.WvC.MN106_top20, row.names = NULL),
                                    c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.WvC.MN106_top20)
class(results.WvC.MN106_top20)
results.WvC.MN106_top20 <- as.data.frame(results.WvC.MN106_top20)
Top20_AT <- inner_join(results.WvC.MN106_top20, orthologs2, by = "PennycressGeneID")
head(Top20_AT)
Top20_AT <- Top20_AT %>%
  dplyr::select(PennycressGeneID, ArabidopsisGeneID, Entrez_GeneID, ProteinFunction)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/WvC.MN106_top20.csv"
write.table(Top20_AT, file=file, row.names=T)

write.csv(Top20_AT, file = file, row.names = TRUE)
library(gridExtra)
grid.table(Top20_AT)

results.WvC.MN106_top20_list <- results.WvC.MN106_top20 %>%
  dplyr::select(PennycressGeneID)
results.WvC.MN106_top20_list <- list(results.WvC.MN106_top20_list$PennycressGeneID)
results.WvC.MN106_top20_list <- unlist(results.WvC.MN106_top20_list)

#other heatmap options besides top 20
results.WvC.MN106_padj.01 <- subset(results.WvC.MN106, results.WvC.MN106$padj<=0.000001)
results.WvC.MN106_padj.01 <- setNames(cbind(rownames(results.WvC.MN106_padj.01), results.WvC.MN106_padj.01, row.names = NULL),
                                      c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.WvC.MN106_padj.01 <- as.data.frame(results.WvC.MN106_padj.01)
padj.01 <- inner_join(results.WvC.MN106_padj.01, orthologs2, by = "PennycressGeneID")
padj.01 <- padj.01 %>%
  dplyr::select(PennycressGeneID)
padj.01_list <- list(padj.01$PennycressGeneID)
padj.01_list <- unlist(padj.01_list)

results.WvC.MN106_top50 <- results.WvC.MN106[1:50,]
results.WvC.MN106_top50 <- setNames(cbind(rownames(results.WvC.MN106_top50), results.WvC.MN106_top50, row.names = NULL),
                                    c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.WvC.MN106_top50)
results.WvC.MN106_top50 <- as.data.frame(results.WvC.MN106_top50)
results.WvC.MN106_top50_list <- results.WvC.MN106_top50 %>%
  dplyr::select(PennycressGeneID)
results.WvC.MN106_top50_list <- list(results.WvC.MN106_top50_list$PennycressGeneID)
results.WvC.MN106_top50_list <- unlist(results.WvC.MN106_top50_list)

#Top 25 up and down in MN106 wl vs controls
res_de_dw_WvC.MN106_top25 <- as.data.frame(res_de_dw_WvC.MN106_top25)
res_de_dw_WvC.MN106_top25 <- setNames(cbind(rownames(res_de_dw_WvC.MN106_top25), res_de_dw_WvC.MN106_top25, row.names = NULL),
                                      c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_up_WvC.MN106_top25 <- as.data.frame(res_de_up_WvC.MN106_top25)
res_de_up_WvC.MN106_top25 <- setNames(cbind(rownames(res_de_up_WvC.MN106_top25), res_de_up_WvC.MN106_top25, row.names = NULL),
                                      c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_dw_WvC.MN106_top25_list <- res_de_dw_WvC.MN106_top25 %>%
  dplyr::select(PennycressGeneID)
res_de_dw_WvC.MN106_top25_list <- list(res_de_dw_WvC.MN106_top25_list$PennycressGeneID)
res_de_dw_WvC.MN106_top25_list <- unlist(res_de_dw_WvC.MN106_top25_list)
res_de_up_WvC.MN106_top25_list <- res_de_up_WvC.MN106_top25 %>%
  dplyr::select(PennycressGeneID)
res_de_up_WvC.MN106_top25_list <- list(res_de_up_WvC.MN106_top25_list$PennycressGeneID)
res_de_up_WvC.MN106_top25_list <- unlist(res_de_up_WvC.MN106_top25_list)
WvC.MN106_top50 <- append(res_de_up_WvC.MN106_top25_list, res_de_dw_WvC.MN106_top25_list)

#Top 25 up and down in sp32-10 wl vs controls
res_de_dw_WvC.SP32_top25 <- as.data.frame(res_de_dw_WvC.SP32_top25)
res_de_dw_WvC.SP32_top25 <- setNames(cbind(rownames(res_de_dw_WvC.SP32_top25), res_de_dw_WvC.SP32_top25, row.names = NULL),
                                     c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_up_WvC.SP32_top25 <- as.data.frame(res_de_up_WvC.SP32_top25)
res_de_up_WvC.SP32_top25 <- setNames(cbind(rownames(res_de_up_WvC.SP32_top25), res_de_up_WvC.SP32_top25, row.names = NULL),
                                     c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_dw_WvC.SP32_top25_list <- res_de_dw_WvC.SP32_top25 %>%
  dplyr::select(PennycressGeneID)
res_de_dw_WvC.SP32_top25_list <- list(res_de_dw_WvC.SP32_top25_list$PennycressGeneID)
res_de_dw_WvC.SP32_top25_list <- unlist(res_de_dw_WvC.SP32_top25_list)
res_de_up_WvC.SP32_top25_list <- res_de_up_WvC.SP32_top25 %>%
  dplyr::select(PennycressGeneID)
res_de_up_WvC.SP32_top25_list <- list(res_de_up_WvC.SP32_top25_list$PennycressGeneID)
res_de_up_WvC.SP32_top25_list <- unlist(res_de_up_WvC.SP32_top25_list)
WvC.SP32_top50 <- append(res_de_up_WvC.SP32_top25_list, res_de_dw_WvC.SP32_top25_list)

#Change geneID to protein function
protein <- assay(rld)[WvC.MN106_top50, 6:14]
protein_1 <- setNames(cbind(rownames(protein), protein, row.names = NULL))
colnames(protein_1)[1] ="PennycressGeneID"
protein_1 <- as.data.frame(protein_1)
protein_2 <- inner_join(protein_1, orthologs, by = "PennycressGeneID")
protein_3 <- protein_2 %>%
  unite("PandG", c('Arabidopsis_Protein_Name', 'PennycressGeneID'), sep= "-", 
        remove = FALSE)
protein_3 <- protein_3 %>% remove_rownames %>% column_to_rownames(var="PandG")
protein_4 <- protein_3 %>%
  dplyr::select("MCR2", "MCR3", "MCR4", "MCR5",
                "MWR1", "MWR2", "MWR3", "MWR4", "MWR5")
protein_5 <- as.matrix(sapply(protein_4, as.numeric), rownames.force = TRUE)
protein_5 <- as.matrix(sapply(protein_4, rownames.force=TRUE, as.numeric))
protein_5 <- as.matrix(protein_4, rownames.force=TRUE)
class(protein_5)
print(sapply(protein_7, class))
transform(protein_5, MCR2 = as.numeric(MCR2))
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/MN106_Waterlogged_roots_proteins_top50.csv"
write.csv(protein_5, file=file, row.names=T)
protein_6 <- read.csv("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/old/MN106_Waterlogged_roots_proteins_top50.csv", row.names = 1)
protein_7 <- as.matrix(protein_6, rownames.force=TRUE)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/heatmaps/WvC_MN106_Roots_heatmap_top50.png",
    width=900, height=900)
pheatmap(protein_7, cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)
dev.off()


#Change geneID to protein function
protein <- assay(rld)[WvC.SP32_top50, 6:14]
protein_1 <- setNames(cbind(rownames(protein), protein, row.names = NULL))
colnames(protein_1)[1] ="PennycressGeneID"
protein_1 <- as.data.frame(protein_1)
protein_2 <- inner_join(protein_1, orthologs, by = "PennycressGeneID")
protein_3 <- protein_2 %>%
  unite("PandG", c('Arabidopsis_Protein_Name', 'PennycressGeneID'), sep= "-", 
        remove = FALSE)
protein_3 <- protein_3 %>% remove_rownames %>% column_to_rownames(var="PandG")
protein_4 <- protein_3 %>%
  dplyr::select("MCR2", "MCR3", "MCR4", "MCR5",
                "MWR1", "MWR2", "MWR3", "MWR4", "MWR5")
protein_5 <- as.matrix(sapply(protein_4, as.numeric), rownames.force = TRUE)
protein_5 <- as.matrix(sapply(protein_4, rownames.force=TRUE, as.numeric))
protein_5 <- as.matrix(protein_4, rownames.force=TRUE)
class(protein_5)
print(sapply(protein_7, class))
transform(protein_5, MCR2 = as.numeric(MCR2))
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_proteins_top50.csv"
write.csv(protein_5, file=file, row.names=T)
protein_6 <- read.csv("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/DEG results/SP32_Waterlogged_roots_top50.csv", row.names = 1)
protein_7 <- as.matrix(protein_6, rownames.force=TRUE)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/heatmaps/WvC_SP32_Roots_heatmap_top50.png",
    width=900, height=900)
pheatmap(protein_7, cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)
dev.off()


results.WvC.SP32_top20 <- results.WvC.SP32[1:20,]
results.WvC.SP32_top20 <- setNames(cbind(rownames(results.WvC.SP32_top20), results.WvC.SP32_top20, row.names = NULL),
                                   c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.WvC.SP32_top20)
results.WvC.SP32_top20 <- as.data.frame(results.WvC.SP32_top20)
results.WvC.SP32_top20_list <- results.WvC.SP32_top20 %>%
  dplyr::select(geneID)
results.WvC.SP32_top20_list <- list(results.WvC.SP32_top20_list$geneID)
results.WvC.SP32_top20_list <- unlist(results.WvC.SP32_top20_list)

results.WvW_top20 <- results.WvW[1:20,]
results.WvW_top20 <- setNames(cbind(rownames(results.WvW_top20), results.WvW_top20, row.names = NULL),
                              c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.WvW_top20)
results.WvW_top20 <- as.data.frame(results.WvW_top20)
results.WvW_top20_list <- results.WvW_top20 %>%
  dplyr::select(geneID)
results.WvW_top20_list <- list(results.WvW_top20_list$geneID)
results.WvW_top20_list <- unlist(results.WvW_top20_list)

results.WvW2_top20 <- results.WvW2[1:20,]
results.WvW2_top20 <- setNames(cbind(rownames(results.WvW2_top20), results.WvW2_top20, row.names = NULL),
                               c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.WvW2_top20)
results.WvW2_top20 <- as.data.frame(results.WvW2_top20)
results.WvW2_top20_list <- results.WvW2_top20 %>%
  dplyr::select(geneID)
results.WvW2_top20_list <- list(results.WvW2_top20_list$geneID)
results.WvW2_top20_list <- unlist(results.WvW2_top20_list)

results.RvC.MN106_top20 <- results.RvC.MN106[1:20,]
results.RvC.MN106_top20 <- setNames(cbind(rownames(results.RvC.MN106_top20), results.RvC.MN106_top20, row.names = NULL),
                                    c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.RvC.MN106_top20)
class(results.RvC.MN106_top20)
results.RvC.MN106_top20 <- as.data.frame(results.RvC.MN106_top20)
results.RvC.MN106_top20_list <- results.RvC.MN106_top20 %>%
  dplyr::select(geneID)
results.RvC.MN106_top20_list <- list(results.RvC.MN106_top20_list$geneID)
results.RvC.MN106_top20_list <- unlist(results.RvC.MN106_top20_list)

results.RvC.SP32_top20 <- results.RvC.SP32[1:20,]
results.RvC.SP32_top20 <- setNames(cbind(rownames(results.RvC.SP32_top20), results.RvC.SP32_top20, row.names = NULL),
                                   c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.RvC.SP32_top20)
class(results.RvC.SP32_top20)
results.RvC.SP32_top20 <- as.data.frame(results.RvC.SP32_top20)
results.RvC.SP32_top20_list <- results.RvC.SP32_top20 %>%
  dplyr::select(geneID)
results.RvC.SP32_top20_list <- list(results.RvC.SP32_top20_list$geneID)
results.RvC.SP32_top20_list <- unlist(results.RvC.SP32_top20_list)

results.RvR_top20 <- results.RvR[1:20,]
results.RvR_top20 <- setNames(cbind(rownames(results.RvR_top20), results.RvR_top20, row.names = NULL),
                              c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.RvR_top20)
class(results.RvR_top20)
results.RvR_top20 <- as.data.frame(results.RvR_top20)
results.RvR_top20_list <- results.RvR_top20 %>%
  dplyr::select(geneID)
results.RvR_top20_list <- list(results.RvR_top20_list$geneID)
results.RvR_top20_list <- unlist(results.RvR_top20_list)

results.RvW.MN106_top20 <- results.RvW.MN106[1:4,]
results.RvW.MN106_top20 <- setNames(cbind(rownames(results.RvW.MN106_top20), results.RvW.MN106_top20, row.names = NULL),
                                    c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
head(results.RvW.MN106_top20)
class(results.RvW.MN106_top20)
results.RvW.MN106_top20 <- as.data.frame(results.RvW.MN106_top20)
results.RvW.MN106_top20_list <- results.RvW.MN106_top20 %>%
  dplyr::select(geneID)
results.RvW.MN106_top20_list <- list(results.RvW.MN106_top20_list$geneID)
results.RvW.MN106_top20_list <- unlist(results.RvW.MN106_top20_list)


# Data frame with column annotations.
df <- as.data.frame(colData(dds)[,c("condition","genotype")])

#OPlot the heatmap - can specify list of top genes (row spot) and which samples (column spot) from rld assay
pheatmap(assay(rld)[results.WvC.MN106_top20_list, 6:14], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

pheatmap(assay(rld)[WvC.MN106_top50, 6:14], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)


pheatmap(assay(rld)[results.WvC.SP32_top20_list, 20:29], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

sample_list<- c("MCR2", "MCR3", "MCR4", "MCR5",
                "MWR1", "MWR2", "MWR3", "MWR4", "MWR5",
                "SCR1", "SCR2", "SCR3", "SCR4", "SCR5",
                "SWR1", "SWR2", "SWR3", "SWR4", "SWR5")
pheatmap(assay(rld)[results.WvW_top20_list, sample_list], cluster_row=T, cutree_cols = 4, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1.5, 0, 1.5), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

sample_list<- c("MWR1", "MWR2", "MWR3", "MWR4", "MWR5",
                "SWR1", "SWR2", "SWR3", "SWR4", "SWR5")
pheatmap(assay(rld)[results.WvW_top20_list, sample_list], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1.5, 0, 1.5), legend_labels = c("Low", "Medium", "High"), fontsize = 12)


pheatmap(assay(rld)[results.WvW2_top20_list, sample_list], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1.5, 0, 1.5), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

pheatmap(assay(rld)[results.RvC.MN106_top20_list, 1:9], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-2, 0, 2), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

pheatmap(assay(rld)[results.RvC.SP32_top20_list,15:24], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-2, 0, 2), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

sample_list<- c("MCR2", "MCR3", "MCR4", "MCR5",
                "MAR1", "MAR2", "MAR3", "MAR4", "MAR5",
                "SCR1", "SCR2", "SCR3", "SCR4", "SCR5",
                "SAR1", "SAR2", "SAR3", "SAR4", "SAR5")
pheatmap(assay(rld)[results.RvR_top20_list, sample_list], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-2, 0, 2), legend_labels = c("Low", "Medium", "High"), fontsize = 12)

sample_list<- c("MWR1", "MWR2", "MWR3", "MWR4", "MWR5",
                "MAR1", "MAR2", "MAR3", "MAR4", "MAR5")
pheatmap(assay(rld)[results.RvW.MN106_top20_list, sample_list], cluster_row=T, scale="row", annotation_col=df, legend = TRUE, display_numbers = FALSE, 
         legend_breaks= c(-1, 0, 1), legend_labels = c("Low", "Medium", "High"), fontsize = 12)



#GSEA ----
#Read in Arabidopsis ortholog file
orthologs <- read.csv("Orthologs_FINAL_no-duplicates.csv", header= TRUE)
head(orthologs)
#set Arabidopsis organism here
organism = "org.At.tair.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#GSEA with Arabidopsis Orthologs - MN106 Waterlogged vs Control Roots
results.WvC.MN106_ALL <- setNames(cbind(rownames(results.WvC.MN106), results.WvC.MN106, row.names = NULL),
                                  c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.WvC.MN106_ALL <- as.data.frame(results.WvC.MN106_ALL)
head(results.WvC.MN106_ALL)
results.WvC.MN106_ALL_GSEA <- results.WvC.MN106_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.WvC.MN106_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.WvC.MN106_ALL_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.WvC.MN106_df <- as.data.frame(resLFC.WvC.MN106)
resLFC.WvC.MN106_df  <- setNames(cbind(rownames(resLFC.WvC.MN106_df ), resLFC.WvC.MN106_df, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.WvC.MN106_df  <- as.data.frame(resLFC.WvC.MN106_df )
head(resLFC.WvC.MN106_df )
resLFC.WvC.MN106_df_GSEA <- resLFC.WvC.MN106_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.WvC.MN106_df_GSEA)
geneList <- inner_join(resLFC.WvC.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_WvC_MN106_Roots <- gseGO(geneList=geneList, 
                             ont ="ALL", 
                             keyType = "TAIR", 
                             #nPerm = 10000, 
                             minGSSize = 10, 
                             maxGSSize = 500, 
                             pvalueCutoff = 0.05,
                             verbose = TRUE,
                             OrgDb = organism, 
                             pAdjustMethod = "BH")


simplify <- clusterProfiler::simplify(gse_WvC_MN106_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_WvC_MN106_Roots_df <- as.data.frame(gse_WvC_MN106_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_roots_GSEA.csv"
write.csv(as.data.frame(gse_WvC_MN106_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC__MN106_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_WvC_MN106_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#tree plot
gse_WvC_MN106_Roots2 <- pairwise_termsim(gse_WvC_MN106_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_gsea_tree.png",
    width=900, height=600)
treeplot(
  gse_WvC_MN106_Roots2,
  showCategory = 30,
  nWords = 4,
  ncluster = 5,
  color = "p.adjust",
  xlim = c(0, 40),
  label_format = 18,
  fontsize = 6,
  split = NULL,
  hclust_method = "average",
)
dev.off()



#KEGG preparation
geneList <- inner_join(resLFC.WvC.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_WvC_MN106_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                    keyType = "ncbi-geneid",
                                    organism = "ath",
                                    exponent = 1,
                                    #nPerm = 10000, 
                                    minGSSize = 10, 
                                    maxGSSize = 500,
                                    eps = 1e-10,
                                    pvalueCutoff = .05,
                                    verbose = TRUE)

gse_WvC_MN106_Roots_KEGG_df <- as.data.frame(gse_WvC_MN106_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_WvC_MN106_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_WvC_MN106_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in MN106 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

browseKEGG(gse_WvC_MN106_Roots_KEGG, 'ath00010')
library("pathview")
ath00195 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00195",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00195", kegg.native = T)

ath00010 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00010",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00010", kegg.native = T)
ath00051 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00051",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00051", kegg.native = T)
ath00073 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00073",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00073", kegg.native = T)
ath00564 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00564",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00564", kegg.native = T)


#ORA for upregulated
res_de_up_WvC.MN106.GSEA <-  setNames(cbind(rownames(res_de_up_WvC.MN106_padj), res_de_up_WvC.MN106_padj, row.names = NULL),
                                      c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_up_WvC.MN106.GSEA <- as.data.frame(res_de_up_WvC.MN106.GSEA)
head(res_de_up_WvC.MN106.GSEA)
res_de_up_WvC.MN106.GSEA <- res_de_up_WvC.MN106.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_up_WvC.MN106.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=10)

#ORA for downregulated
res_de_dw_WvC.MN106.GSEA <-  setNames(cbind(rownames(res_de_dw_WvC.MN106_padj), res_de_dw_WvC.MN106_padj, row.names = NULL),
                                      c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_dw_WvC.MN106.GSEA <- as.data.frame(res_de_dw_WvC.MN106.GSEA)
head(res_de_dw_WvC.MN106.GSEA)
res_de_dw_WvC.MN106.GSEA <- res_de_dw_WvC.MN106.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_dw_WvC.MN106.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=10)

#ORA KEGG
ORA_WvC_MN106_Roots_KEGG <- enrichKEGG(gene = geneList3,
                                       #keyType = "ncbi-geneid",
                                       organism = "ath",
                                       pvalueCutoff = 1)

ORA_WvC_MN106_Roots_KEGG <- as.data.frame(ORA_WvC_MN106_Roots_KEGG)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_ORA_kegg.png",
    width=900, height=900)
barplot(ORA_WvC_MN106_Roots_KEGG, 
        drop = TRUE)
dotplot(ORA_WvC_MN106_Roots_KEGG)
dotplot(gse_WvC_MN106_Roots_KEGG, showCategory=10, font.size=20, label_format=30)
ggtitle("GSEA of Waterlogged Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()



## convert gene ID to Symbol
ath_symb <- setReadable(gse_WvC_MN106_Roots, 'org.At.tair.db', 'TAIR')

p1 <- cnetplot(gse_WvC_MN106_Roots, foldChange=geneList, node_label="category", 
               cex_label_category = 1.2)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(gse_WvC_MN106_Roots, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(gse_WvC_MN106_Roots, foldChange=geneList, showCategory = 2, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p3)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

#heatmap of GO terms
p1 <- heatplot(gse_WvC_MN106_Roots, showCategory=5)
p2 <- heatplot(gse_WvC_MN106_Roots, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

#tree plot
gse_WvC_MN106_Roots2 <- pairwise_termsim(gse_WvC_MN106_Roots)
p1 <- treeplot(gse_WvC_MN106_Roots2)
p2 <- treeplot(gse_WvC_MN106_Roots2, hclust_method = "average")
p2
aplot::plot_list(p1, p2, tag_levels='A')


gse_WvC_MN106_Roots <- pairwise_termsim(gse_WvC_MN106_Roots)
p1 <- emapplot(gse_WvC_MN106_Roots)
p2 <- emapplot(gse_WvC_MN106_Roots, cex_category=1.5)
p3 <- emapplot(gse_WvC_MN106_Roots, layout="kk")
p4 <- emapplot(gse_WvC_MN106_Roots, cex_category=1.5,layout="kk")
p4
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
cowplot::plot_grid(p1)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(gse_WvC_MN106_Roots, label_format = 80, showCategory = 6) + ggtitle("GSEA MN106 Waterlogged Roots") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15))
dev.off()


gse_WvC_MN106_Roots_ridge <- gse_WvC_MN106_Roots %>% filter(Description == "generation of precursor metabolites and energy")
gse_WvC_MN106_Roots_ridge

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_MN106_Roots_gsea_ridge_energy.png",
    width=800, height=900)
ridgeplot(gse_WvC_MN106_Roots_ridge, label_format = 80, showCategory = 1) + ggtitle("GSEA MN106 Waterlogged Roots") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15))
dev.off()


#GSEA with Arabidopsis Orthologs - SP32-10 Waterlogged vs Control Roots
results.WvC.SP32_ALL <- setNames(cbind(rownames(results.WvC.SP32), results.WvC.SP32, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.WvC.SP32_ALL <- as.data.frame(results.WvC.SP32_ALL)
head(results.WvC.SP32_ALL)
results.WvC.SP32_ALL_GSEA <- results.WvC.SP32_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.WvC.SP32_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.WvC.SP32_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.WvC.SP32_df <- as.data.frame(resLFC.WvC.SP32)
resLFC.WvC.SP32_df  <- setNames(cbind(rownames(resLFC.WvC.SP32_df ), resLFC.WvC.SP32_df, row.names = NULL),
                                c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.WvC.SP32_df  <- as.data.frame(resLFC.WvC.SP32_df )
head(resLFC.WvC.SP32_df )
resLFC.WvC.SP32_df_GSEA <- resLFC.WvC.SP32_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.WvC.SP32_df_GSEA)
geneList <- inner_join(resLFC.WvC.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)


set.seed(1234)
gse_WvC_SP32_Roots <- gseGO(geneList=geneList, 
                            ont ="ALL", 
                            keyType = "TAIR", 
                            #nPerm = 10000, 
                            minGSSize = 10, 
                            maxGSSize = 500, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE,
                            OrgDb = organism, 
                            pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_WvC_SP32_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_WvC_SP32_Roots_df <- as.data.frame(gse_WvC_SP32_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_roots_GSEA.csv"
write.csv(as.data.frame(gse_WvC_SP32_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_WvC_SP32_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in SP32-10") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in SP32-10") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(simplify, label_format = 80, showCategory = 20) + ggtitle("GSEA - 7 Days Waterlogged Roots SP32-10") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15))
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_ridge_ALL.jpg", width = 35, height = 23, unit = "cm", dpi = 300)
dev.off()

#tree plot
gse_WvC_SP32_Roots2 <- pairwise_termsim(gse_WvC_SP32_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_tree.png",
    width=1100, height=550)
treeplot(
  gse_WvC_SP32_Roots2,
  showCategory = 30,
  nWords = 4,
  ncluster = 5,
  color = "p.adjust",
  xlim = c(80, 0),
  label_format = 18,
  fontsize = 6,
  split = NULL,
  hclust_method = "average",
)
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.WvC.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_WvC_SP32_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                   keyType = "ncbi-geneid",
                                   organism = "ath",
                                   exponent = 1,
                                   #nPerm = 10000, 
                                   minGSSize = 10, 
                                   maxGSSize = 500,
                                   eps = 1e-10,
                                   pvalueCutoff = .05,
                                   verbose = TRUE)

gse_WvC_SP32_Roots_KEGG_df <- as.data.frame(gse_WvC_SP32_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_WvC_SP32_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_KEGG.png",
    width=900, height=1100)
dotplot(gse_WvC_SP32_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in SP32 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()



#ORA for upregulated
res_de_up_WvC.SP32.GSEA <-  setNames(cbind(rownames(res_de_up_WvC.SP32_padj), res_de_up_WvC.SP32_padj, row.names = NULL),
                                     c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_up_WvC.SP32.GSEA <- as.data.frame(res_de_up_WvC.SP32.GSEA)
head(res_de_up_WvC.SP32.GSEA)
res_de_up_WvC.SP32.GSEA <- res_de_up_WvC.SP32.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_up_WvC.SP32.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=10)

#ORA for downregulated
res_de_dw_WvC.SP32.GSEA <-  setNames(cbind(rownames(res_de_dw_WvC.SP32_padj), res_de_dw_WvC.SP32_padj, row.names = NULL),
                                     c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_dw_WvC.SP32.GSEA <- as.data.frame(res_de_dw_WvC.SP32.GSEA)
head(res_de_dw_WvC.SP32.GSEA)
res_de_dw_WvC.SP32.GSEA <- res_de_dw_WvC.SP32.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_dw_WvC.SP32.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=10)


gse_WvC_SP32_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                   keyType = "ncbi-geneid",
                                   organism = "ath",
                                   #nPerm = 10000, 
                                   minGSSize = 120, 
                                   #maxGSSize = 500, 
                                   pvalueCutoff = 1,
                                   verbose = TRUE)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvC_SP32_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_WvC_SP32_Roots_KEGG, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Waterlogged Roots in SP32-10") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#GSEA with Arabidopsis Orthologs - MN106 vs SP32-10 Waterlogged Roots
results.WvW_ALL <- setNames(cbind(rownames(results.WvW), results.WvW, row.names = NULL),
                            c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.WvW_ALL <- as.data.frame(results.WvW_ALL)
head(results.WvW_ALL)
results.WvW_ALL_GSEA <- results.WvW_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.WvW_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.WvW_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.WvW_df <- as.data.frame(resLFC.WvW)
resLFC.WvW_df  <- setNames(cbind(rownames(resLFC.WvW_df ), resLFC.WvW_df, row.names = NULL),
                           c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.WvW_df  <- as.data.frame(resLFC.WvW_df )
head(resLFC.WvW_df )
resLFC.WvW_df_GSEA <- resLFC.WvW_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.WvW_df_GSEA)
geneList <- inner_join(resLFC.WvW_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_WvW_Roots <- gseGO(geneList=geneList, 
                       ont ="ALL", 
                       keyType = "TAIR", 
                       #nPerm = 10000, 
                       minGSSize = 10, 
                       maxGSSize = 500, 
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_WvW_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_WvW_Roots_df <- as.data.frame(gse_WvW_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_roots_GSEA.csv"
write.csv(as.data.frame(gse_WvW_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_WvW_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 vs SP32-10 Waterlogged Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 vs SP32-10 Waterlogged Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
exclude <- c('GO:0036294', 'GO:0071453')
gse_WvW_Roots_drop <- dropGO(gse_WvW_Roots, level = BP, term = exclude)
ridgeplot(gse_WvW_Roots, label_format = 80, showCategory = 30) + theme(axis.text = element_text(size = 60, face="bold"))
dev.off()

#tree plot
gse_WvW_Roots2 <- pairwise_termsim(gse_WvW_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_tree.png",
    width=3000, height=900)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_tree.png",
    width=900, height=600)
treeplot(
  gse_WvW_Roots2,
  showCategory = 30,
  label_format = 50,
  nWords = 4,
  fontsize = 9,
  xlim = c(0, 40),
  hclust_method = "average",
)
dev.off()

treeplot(
  gse_WvW_Roots2,
  #showCategory = 30,
  hilight.params = list(hilight = TRUE,
                        align = "left"), 
  offset.params = list(bar_tree = rel(1),
                       tiplab = rel(1.2), 
                       extend = 0.3, 
                       hexpand = 0.3 
  ), 
  cluster.params = list(hclust_method = "average", 
                        label_format = 50,
                        nWords = 4,
                        fontsize = 6,
                        xlim = c(0, 40)))


#KEGG preparation
geneList <- inner_join(resLFC.WvW_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_WvW_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                              keyType = "ncbi-geneid",
                              organism = "ath",
                              exponent = 1,
                              #nPerm = 10000, 
                              minGSSize = 10, 
                              maxGSSize = 500,
                              eps = 1e-10,
                              pvalueCutoff = .05,
                              verbose = TRUE)

gse_WvW_Roots_KEGG_df <- as.data.frame(gse_WvW_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_WvW_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_WvW_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 vs SP32-10 Waterlogged Roots - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

library("pathview")
ath00195 <- pathview(gene.data  = geneList3,
                     pathway.id = "ath00195",
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList3)), cpd=1), out.suffix = "ath00195.WvW", kegg.native = T)



#ORA for upregulated
res_de_up_WvW.GSEA <-  setNames(cbind(rownames(res_de_up_WvW_padj), res_de_up_WvW_padj, row.names = NULL),
                                c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_up_WvW.GSEA <- as.data.frame(res_de_up_WvW.GSEA)
head(res_de_up_WvW.GSEA)
res_de_up_WvW.GSEA <- res_de_up_WvW.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_up_WvW.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
go_enrichdf <- as.data.frame(go_enrich)
dotplot(go_enrich, showCategory=10)

#ORA for downregulated
res_de_dw_WvW.GSEA <-  setNames(cbind(rownames(res_de_dw_WvW_padj), res_de_dw_WvW_padj, row.names = NULL),
                                c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
res_de_dw_WvW.GSEA <- as.data.frame(res_de_dw_WvW.GSEA)
head(res_de_dw_WvW.GSEA)
res_de_dw_WvW.GSEA <- res_de_dw_WvW.GSEA %>%
  dplyr::select(PennycressGeneID)
geneList <- inner_join(res_de_dw_WvW.GSEA, orthologs, by = "PennycressGeneID")
geneList2 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList2 %>%
  dplyr::select(ArabidopsisGeneID)
geneList3 <- geneList3 %>%
  pull(ArabidopsisGeneID)
geneList3 <- na.omit(geneList3)
head(geneList3)

go_enrich <- enrichGO(gene = geneList3,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
dotplot(go_enrich, showCategory=10)


#GSEA with Arabidopsis Orthologs - MN106 Recovery vs Control Roots
results.RvC.MN106_ALL <- setNames(cbind(rownames(results.RvC.MN106), results.RvC.MN106, row.names = NULL),
                                  c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvC.MN106_ALL <- as.data.frame(results.RvC.MN106_ALL)
head(results.RvC.MN106_ALL)
results.RvC.MN106_ALL_GSEA <- results.RvC.MN106_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvC.MN106_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvC.MN106_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.RvC.MN106_df <- as.data.frame(resLFC.RvC.MN106)
resLFC.RvC.MN106_df  <- setNames(cbind(rownames(resLFC.RvC.MN106_df ), resLFC.RvC.MN106_df, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvC.MN106_df  <- as.data.frame(resLFC.RvC.MN106_df )
head(resLFC.RvC.MN106_df )
resLFC.RvC.MN106_df_GSEA <- resLFC.RvC.MN106_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvC.MN106_df_GSEA)
geneList <- inner_join(resLFC.RvC.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvC_MN106_Roots <- gseGO(geneList=geneList, 
                             ont ="ALL", 
                             keyType = "TAIR", 
                             #nPerm = 10000, 
                             minGSSize = 10, 
                             maxGSSize = 500, 
                             pvalueCutoff = 0.05,
                             verbose = TRUE,
                             OrgDb = organism, 
                             pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvC_MN106_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvC_MN106_Roots_df <- as.data.frame(gse_RvC_MN106_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_MN106_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvC_MN106_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_MN106_Roots_gsea.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Control Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.RvC.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvC_MN106_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                    keyType = "ncbi-geneid",
                                    organism = "ath",
                                    exponent = 1,
                                    #nPerm = 10000, 
                                    minGSSize = 10, 
                                    maxGSSize = 500,
                                    eps = 1e-10,
                                    pvalueCutoff = .05,
                                    verbose = TRUE)

gse_RvC_MN106_Roots_KEGG_df <- as.data.frame(gse_RvC_MN106_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_MN106_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvC_MN106_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_MN106_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvC_MN106_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Control Roots in MN106 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#GSEA with Arabidopsis Orthologs - SP32-10 Recovery vs Control Roots
results.RvC.SP32_ALL <- setNames(cbind(rownames(results.RvC.SP32), results.RvC.SP32, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvC.SP32_ALL <- as.data.frame(results.RvC.SP32_ALL)
head(results.RvC.SP32_ALL)
results.RvC.SP32_ALL_GSEA <- results.RvC.SP32_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvC.SP32_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvC.SP32_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.RvC.SP32_df <- as.data.frame(resLFC.RvC.SP32)
resLFC.RvC.SP32_df  <- setNames(cbind(rownames(resLFC.RvC.SP32_df ), resLFC.RvC.SP32_df, row.names = NULL),
                                c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvC.SP32_df  <- as.data.frame(resLFC.RvC.SP32_df )
head(resLFC.RvC.SP32_df )
resLFC.RvC.SP32_df_GSEA <- resLFC.RvC.SP32_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvC.SP32_df_GSEA)
geneList <- inner_join(resLFC.RvC.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvC_SP32_Roots <- gseGO(geneList=geneList, 
                            ont ="ALL", 
                            keyType = "TAIR", 
                            #nPerm = 10000, 
                            minGSSize = 10, 
                            maxGSSize = 500, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE,
                            OrgDb = organism, 
                            pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvC_SP32_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvC_SP32_Roots_df <- as.data.frame(gse_RvC_SP32_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_SP32_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvC_SP32_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_SP32_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_RvC_SP32_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Control Roots in SP32") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.RvC.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvC_SP32_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                   keyType = "ncbi-geneid",
                                   organism = "ath",
                                   exponent = 1,
                                   #nPerm = 10000, 
                                   minGSSize = 10, 
                                   maxGSSize = 500,
                                   eps = 1e-10,
                                   pvalueCutoff = .05,
                                   verbose = TRUE)

gse_RvC_SP32_Roots_KEGG_df <- as.data.frame(gse_RvC_SP32_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_SP32_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvC_SP32_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvC_SP32_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvC_SP32_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Control Roots in SP32 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#GSEA with Arabidopsis Orthologs - MN106 vs SP32-10 Recovery Roots
results.RvR_ALL <- setNames(cbind(rownames(results.RvR), results.RvR, row.names = NULL),
                            c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvR_ALL <- as.data.frame(results.RvR_ALL)
head(results.RvR_ALL)
results.RvR_ALL_GSEA <- results.RvR_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvR_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvR_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.RvR_df <- as.data.frame(resLFC.RvR)
resLFC.RvR_df  <- setNames(cbind(rownames(resLFC.RvR_df ), resLFC.RvR_df, row.names = NULL),
                           c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvR_df  <- as.data.frame(resLFC.RvR_df )
head(resLFC.RvR_df )
resLFC.RvR_df_GSEA <- resLFC.RvR_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvR_df_GSEA)
geneList <- inner_join(resLFC.RvR_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvR_Roots <- gseGO(geneList=geneList, 
                       ont ="ALL", 
                       keyType = "TAIR", 
                       #nPerm = 10000, 
                       minGSSize = 10, 
                       maxGSSize = 500, 
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvR_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvR_Roots_df <- as.data.frame(gse_RvR_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_roots_GSEA.csv"
write.csv(as.data.frame(gse_RvR_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_RvR_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 vs SP32-10 Recovery Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 vs SP32-10 Recovery Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(gse_RvR_Roots, label_format = 80, showCategory = 20) + ggtitle("GSEA of MN106 vs SP32-10 Recovery Roots") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=20))
dev.off()

#tree plot
gse_RvR_Roots2 <- pairwise_termsim(gse_RvR_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_tree.png",
    width=1200, height=900)
treeplot(
  gse_RvR_Roots2,
  showCategory = 30,
  label_format = 30,
  nWords = 4,
  fontsize = 6,
  xlim = c(0, 40),
  hclust_method = "average",
)
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.RvR_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvR_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                              keyType = "ncbi-geneid",
                              organism = "ath",
                              exponent = 1,
                              #nPerm = 10000, 
                              minGSSize = 10, 
                              maxGSSize = 500,
                              eps = 1e-10,
                              pvalueCutoff = .05,
                              verbose = TRUE)

gse_RvR_Roots_KEGG_df <- as.data.frame(gse_RvR_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvR_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvR_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 compared to SP32-10 Recovery Roots - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#GSEA with Arabidopsis Orthologs - MN106 Recovery vs Waterlogged Roots
results.RvW.MN106_ALL <- setNames(cbind(rownames(results.RvW.MN106), results.RvW.MN106, row.names = NULL),
                                  c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvW.MN106_ALL <- as.data.frame(results.RvW.MN106_ALL)
head(results.RvW.MN106_ALL)
results.RvW.MN106_ALL_GSEA <- results.RvW.MN106_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvW.MN106_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvW.MN106_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.RvW.MN106_df <- as.data.frame(resLFC.RvW.MN106)
resLFC.RvW.MN106_df  <- setNames(cbind(rownames(resLFC.RvW.MN106_df ), resLFC.RvW.MN106_df, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvW.MN106_df  <- as.data.frame(resLFC.RvW.MN106_df )
head(resLFC.RvW.MN106_df )
resLFC.RvW.MN106_df_GSEA <- resLFC.RvW.MN106_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvW.MN106_df_GSEA)
geneList <- inner_join(resLFC.RvW.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvW_MN106_Roots <- gseGO(geneList=geneList, 
                             ont ="ALL", 
                             keyType = "TAIR", 
                             #nPerm = 10000, 
                             minGSSize = 10, 
                             maxGSSize = 500, 
                             pvalueCutoff = 0.05,
                             verbose = TRUE,
                             OrgDb = organism, 
                             pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvW_MN106_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvW_MN106_Roots_df <- as.data.frame(gse_RvW_MN106_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvW_MN106_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_RvW_MN106_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in MN106") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(gse_RvW_MN106_Roots, label_format = 80, showCategory = 20) + ggtitle("GSEA - Recovery vs Waterlogged Roots MN106") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15, face="bold"))
dev.off()

gse_RvW_MN106_Roots2 <- pairwise_termsim(gse_RvW_MN106_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_gsea_tree.png",
    width=1100, height=650)
treeplot(
  gse_RvW_MN106_Roots2,
  showCategory = 30,
  nWords = 4,
  ncluster = 5,
  color = "p.adjust",
  xlim = c(80, 0),
  label_format = 18,
  fontsize = 6,
  split = NULL,
  hclust_method = "average",
)
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.RvW.MN106_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvW_MN106_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                    keyType = "ncbi-geneid",
                                    organism = "ath",
                                    exponent = 1,
                                    #nPerm = 10000, 
                                    minGSSize = 10, 
                                    maxGSSize = 500,
                                    eps = 1e-10,
                                    pvalueCutoff = .05,
                                    verbose = TRUE)

gse_RvW_MN106_Roots_KEGG_df <- as.data.frame(gse_RvW_MN106_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvW_MN106_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_MN106_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvW_MN106_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in MN106 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#GSEA with Arabidopsis Orthologs - SP32-10 Recovery vs Waterlogged Roots
results.RvW.SP32_ALL <- setNames(cbind(rownames(results.RvW.SP3210), results.RvW.SP3210, row.names = NULL),
                                 c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvW.SP32_ALL <- as.data.frame(results.RvW.SP32_ALL)
head(results.RvW.SP32_ALL)
results.RvW.SP32_ALL_GSEA <- results.RvW.SP32_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvW.SP32_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvW.SP32_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

#LFC shrinkage
resLFC.RvW.SP32_df <- as.data.frame(resLFC.RvW.SP3210)
resLFC.RvW.SP32_df  <- setNames(cbind(rownames(resLFC.RvW.SP32_df ), resLFC.RvW.SP32_df, row.names = NULL),
                                c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvW.SP32_df  <- as.data.frame(resLFC.RvW.SP32_df )
head(resLFC.RvW.SP32_df )
resLFC.RvW.SP32_df_GSEA <- resLFC.RvW.SP32_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvW.SP32_df_GSEA)
geneList <- inner_join(resLFC.RvW.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvW_SP32_Roots <- gseGO(geneList=geneList, 
                            ont ="ALL", 
                            keyType = "TAIR", 
                            #nPerm = 10000, 
                            minGSSize = 10, 
                            maxGSSize = 500, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE,
                            OrgDb = organism, 
                            pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvW_SP32_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvW_SP32_Roots_df <- as.data.frame(gse_RvW_SP32_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvW_SP32_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_RvW_SP32_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in SP32") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in SP32-10") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(gse_RvW_SP32_Roots, label_format = 80, showCategory = 20) + ggtitle("GSEA - Recovery vs Waterlogged Roots SP32-10") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15, face="bold"))
dev.off()

gse_RvW_SP32_Roots2 <- pairwise_termsim(gse_RvW_SP32_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea_tree.png",
    width=1600, height=650)
treeplot(
  gse_RvW_SP32_Roots2,
  showCategory = 30,
  nWords = 4,
  ncluster = 5,
  color = "p.adjust",
  xlim = c(0, 40),
  label_format = 30,
  fontsize = 6,
  split = NULL,
  hclust_method = "average",
)
dev.off()

#KEGG preparation
geneList <- inner_join(resLFC.RvW.SP32_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvW_SP32_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                   keyType = "ncbi-geneid",
                                   organism = "ath",
                                   exponent = 1,
                                   #nPerm = 10000, 
                                   minGSSize = 10, 
                                   maxGSSize = 500,
                                   eps = 1e-10,
                                   pvalueCutoff = .05,
                                   verbose = TRUE)

gse_RvW_SP32_Roots_KEGG_df <- as.data.frame(gse_RvW_SP32_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvW_SP32_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvW_SP32_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in SP32 - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#GSEA with Arabidopsis Orthologs - Recovery vs Waterlogged Roots between accessions
results.RvW.accessions_ALL <- setNames(cbind(rownames(results.RvW.accessions), results.RvW.accessions, row.names = NULL),
                                       c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvW.accessions_ALL <- as.data.frame(results.RvW.accessions_ALL)
head(results.RvW.accessions_ALL)
results.RvW.accessions_ALL_GSEA <- results.RvW.accessions_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvW.accessions_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvW.accessions_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvW_accessions_Roots <- gseGO(geneList=geneList, 
                                  ont ="ALL", 
                                  keyType = "TAIR", 
                                  #nPerm = 10000, 
                                  minGSSize = 10, 
                                  maxGSSize = 500, 
                                  pvalueCutoff = 0.05,
                                  verbose = TRUE,
                                  OrgDb = organism, 
                                  pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvW_accessions_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvW_accessions_Roots_df <- as.data.frame(gse_RvW_accessions_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_accessions_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvW_accessions_Roots_df), file=file,  row.names=T)


png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_accessions_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots between accessions") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#KEGG preparation
geneList <- inner_join(results.RvW.accessions_ALL_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_RvW_accessions_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                                         keyType = "ncbi-geneid",
                                         organism = "ath",
                                         exponent = 1,
                                         #nPerm = 10000, 
                                         minGSSize = 10, 
                                         maxGSSize = 500,
                                         eps = 1e-10,
                                         pvalueCutoff = .05,
                                         verbose = TRUE)

gse_RvW_accessions_Roots_KEGG_df <- as.data.frame(gse_RvW_accessions_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_accessions_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_RvW_accessions_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_accessions_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_RvW_accessions_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots between accessions - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

#GSEA with Arabidopsis Orthologs - SP32-10 vs MN106 Control Roots
results.CvC_ALL <- setNames(cbind(rownames(results.CvC), results.CvC, row.names = NULL),
                            c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.CvC_ALL <- as.data.frame(results.CvC_ALL)
head(results.CvC_ALL)
results.CvC_ALL_GSEA <- results.CvC_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.CvC_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.CvC_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_RvW_SP32_Roots <- gseGO(geneList=geneList, 
                            ont ="ALL", 
                            keyType = "TAIR", 
                            #nPerm = 10000, 
                            minGSSize = 10, 
                            maxGSSize = 500, 
                            pvalueCutoff = 0.05,
                            verbose = TRUE,
                            OrgDb = organism, 
                            pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_RvW_SP32_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_RvW_SP32_Roots_df <- as.data.frame(gse_RvW_SP32_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_GSEA.csv"
write.csv(as.data.frame(gse_RvW_SP32_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_RvW_SP32_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in SP32") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvW_SP32_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of Recovery vs Waterlogged Roots in SP32-10") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#LFC shrinkage
resLFC.CvC_df <- as.data.frame(resLFC.CvC)
resLFC.CvC_df  <- setNames(cbind(rownames(resLFC.CvC_df ), resLFC.CvC_df, row.names = NULL),
                           c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.CvC_df  <- as.data.frame(resLFC.CvC_df )
head(resLFC.CvC_df )
resLFC.CvC_df_GSEA <- resLFC.CvC_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.CvC_df_GSEA)
geneList <- inner_join(resLFC.CvC_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList <- geneList %>%
  pull(log2FoldChange, ArabidopsisGeneID)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- na.omit(geneList)
head(geneList)

set.seed(1234)
gse_CvC_Roots <- gseGO(geneList=geneList, 
                       ont ="ALL", 
                       keyType = "TAIR", 
                       #nPerm = 10000, 
                       minGSSize = 10, 
                       maxGSSize = 500, 
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       OrgDb = organism, 
                       pAdjustMethod = "BH")

simplify <- clusterProfiler::simplify(gse_CvC_Roots, cutoff=0.7, by="p.adjust", select_fun=min)

gse_CvC_Roots_df <- as.data.frame(gse_CvC_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/gse_CvC_Roots_GSEA.csv"
write.csv(as.data.frame(gse_CvC_Roots_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/gse_CvC_Roots_gsea.png",
    width=900, height=900)
dotplot(gse_CvC_Roots, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of SP32 vs MN106 Control Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/gse_CvC_Roots_gsea_ALL.png",
    width=900, height=900)
dotplot(simplify, showCategory=10, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of SP32 vs MN106 Control Roots") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/gse_CvC_Roots_gsea_ridge_ALL.png",
    width=800, height=900)
ridgeplot(gse_CvC_Roots, label_format = 80, showCategory = 20) + ggtitle("GSEA of SP32 vs MN106 Control Roots") + theme(axis.text = element_text(size = 60, face="bold"), plot.title = element_text(size=15, face="bold"))
dev.off()

gse_CvC_Roots2 <- pairwise_termsim(gse_CvC_Roots)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/gse_CvC_Roots_gsea_tree.png",
    width=1600, height=650)
treeplot(
  gse_CvC_Roots2,
  showCategory = 30,
  nWords = 4,
  ncluster = 5,
  color = "p.adjust",
  xlim = c(0, 40),
  label_format = 30,
  fontsize = 6,
  split = NULL,
  hclust_method = "average",
)
dev.off()


#KEGG preparation
geneList <- inner_join(resLFC.CvC_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- as.data.frame(geneList)
geneList <- geneList %>%
  dplyr::select(log2FoldChange, ArabidopsisGeneID)
geneList2 <- geneList %>%
  pull(ArabidopsisGeneID)
eg2np <- bitr(geneList2, fromType='TAIR', toType='ENTREZID', OrgDb = organism)
colnames(eg2np)[1] ="ArabidopsisGeneID"
geneList <- as.data.frame(geneList)
geneList3 <- inner_join(geneList, eg2np, by = "ArabidopsisGeneID")
geneList3 <- geneList3 %>%
  dplyr::select(log2FoldChange, ENTREZID)
geneList3 <- geneList3 %>%
  pull(log2FoldChange, ENTREZID)
geneList3 <- sort(geneList3, decreasing = TRUE)
geneList3 <- na.omit(geneList3)
head(geneList3)

set.seed(1234)
gse_CvC_Roots_KEGG <- gseKEGG(geneList=geneList3, 
                              keyType = "ncbi-geneid",
                              organism = "ath",
                              exponent = 1,
                              #nPerm = 10000, 
                              minGSSize = 10, 
                              maxGSSize = 500,
                              eps = 1e-10,
                              pvalueCutoff = .05,
                              verbose = TRUE)

gse_CvC_Roots_KEGG_df <- as.data.frame(gse_CvC_Roots_KEGG)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/CvC_roots_GSEA_KEGG.csv"
write.csv(as.data.frame(gse_CvC_Roots_KEGG_df), file=file,  row.names=T)

png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/CvC_Roots_gsea_KEGG.png",
    width=900, height=900)
dotplot(gse_CvC_Roots_KEGG, showCategory=20, split=".sign", font.size=20, label_format=30) + 
  facet_grid(.~.sign) + ggtitle("GSEA of MN106 compared to SP32-10 Control Roots - KEGG") +
  theme(plot.title=element_text(size=26,face="bold", hjust=1), legend.text=element_text(size=24), legend.title=element_text(size=24), strip.text = element_text(size = 24, margin = margin()))
dev.off()


#Goplots####
library(GOplot)
head(gse_WvC_MN106_Roots)
file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/MN106_Waterlogged_roots_GSEA-ALL.csv"
write.csv(as.data.frame(gse_WvC_MN106_Roots), file=file, row.names=T)
#edited above file in excel - remove dashes and replace with commas
GOplot_gsea <- read.csv("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/MN106_Waterlogged_roots_GSEA-ALL.csv")

results.WvC.MN106_ALL <- results(dds, name="genotypeMN106.conditionwaterlogged")
results.WvC.MN106_ALL <- setNames(cbind(rownames(results.WvC.MN106_ALL), results.WvC.MN106_ALL, row.names = NULL),
                                  c("PennycressGeneID", "AveExpr", "logFC", "lfcSE", "stat", "P.Value", "adj.P.Val"))
results.WvC.MN106_ALL <- as.data.frame(results.WvC.MN106_ALL)
head(results.WvC.MN106_ALL)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.WvC.MN106_ALL, orthologs, by = "PennycressGeneID")
geneList <- na.omit(geneList)
head(geneList)
geneList <- geneList %>%
  dplyr::select(ArabidopsisGeneID, logFC, AveExpr, stat, P.Value, adj.P.Val)
geneList <- na.omit(geneList)
colnames(geneList)[1] ="ID"
colnames(geneList)[2] ="logFC"
circ <- circle_dat(GOplot_gsea, geneList)

GOCircle(circ)
GOBubble(circ, labels = 1, ID=FALSE)
reduced_circ <- reduce_overlap(circ, overlap = 0.6)
# ...and plot it
GOBubble(reduced_circ, labels = 1, ID=FALSE)


GOBar(subset(circ, category == 'BP'))

chord <- chord_dat(circ, GOplot_gsea$Genes, GOplot_gsea$Term)
head(chord)

l1 <- subset(circ, term == 'response to hypoxia', c(genes,logFC))
l2 <- subset(circ, term == 'plant-type secondary cell wall biogenesis', c(genes,logFC))
l3 <- subset(circ, term == 'translation', c(genes,logFC))
GOVenn(l1, l3, label = c('response to hypoxia', 'translation'))


#WvW
#file <- "~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_roots_GSEA-ALL.csv"
#write.csv(as.data.frame(gse_WvW_Roots), file=file, row.names=T)
#edited above file in excel - remove dashes and replace with commas
GOplot_gsea <- read.csv("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_roots_GSEA_ALL2.csv")


results.WvW <- setNames(cbind(rownames(results.WvW), results.WvW, row.names = NULL),
                        c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.WvW_ALL <- as.data.frame(results.WvW)
head(results.WvW_ALL)
results.WvW_ALL_GSEA <- results.WvW_ALL %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.WvW_ALL_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.WvW_ALL_GSEA, orthologs, by = "PennycressGeneID")
geneList <- na.omit(geneList)
head(geneList)
geneList <- geneList %>%
  dplyr::select(ArabidopsisGeneID, log2FoldChange)
colnames(geneList)[1] ="ID"
colnames(geneList)[2] ="logFC"
head(geneList)

#LFC shrinkage
resLFC.WvW_df <- as.data.frame(resLFC.WvW)
resLFC.WvW_df  <- setNames(cbind(rownames(resLFC.WvW_df ), resLFC.WvW_df, row.names = NULL),
                           c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.WvW_df  <- as.data.frame(resLFC.WvW_df)
head(resLFC.WvW_df )
resLFC.WvW_df_GSEA <- resLFC.WvW_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.WvW_df_GSEA)

geneList <- inner_join(resLFC.WvW_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- na.omit(geneList)
head(geneList)
geneList <- geneList %>%
  dplyr::select(ArabidopsisGeneID, log2FoldChange)
colnames(geneList)[1] ="ID"
colnames(geneList)[2] ="logFC"
head(geneList)


circ <- circle_dat(GOplot_gsea, geneList)


reduced <- c("GO:0009535", "GO:0036294", "GO:0009409", "GO:0009642", "GO:0010200", "GO:0050832", "GO:0009753", "GO:0031347", "GO:0009408", "GO:0006091")
GOCircle(circ, nsub = reduced)
GOCircle(circ)
GOBubble(circ, labels = 1, ID=FALSE)
reduced_circ <- reduce_overlap(circ, overlap = 0.8)
GOBubble(reduced_circ)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_Bubble_colors.png",
    width=1200, height=700)
#reduced_circ$jit<-with(reduced_circ, ifelse(term == "RNA modification" | term == "inorganic ion homeostasis", 0.5, 0.5))
Bubble <- GOBubble(reduced_circ, title = 'MN106 vs SP3210 Waterlogged Roots Enriched GO Terms - GSEA', labels = 10, ID=FALSE, table.legend = FALSE, colour = c("#009E73", "#E69F00", "#56B4E9")) + 
  #geom_text(aes(label= term), size=5, position=position_jitter(width=reduced_circ$jit,height=reduced_circ$jit)) +
  geom_text_repel(aes(label = term), size = 5.7) +
  theme(axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size=20), legend.text = element_text(size=20), title = element_text(size=20))
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/WvW_Roots_gsea_Bubble_colors.jpg", width = 35, height = 23, unit = "cm", dpi = 300)
dev.off()
#we can run this to get a table in grey for now
GOCircle(reduced_circ, nsub = 20)

GOBubble(reduced_circ, labels = 6, ID=FALSE)


#RvR
#edited gsea file in excel - remove dashes and replace with commas
GOplot_gsea <- read.csv("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_roots_GSEA_bubble.csv")


results.RvR <- setNames(cbind(rownames(results.RvR), results.RvR, row.names = NULL),
                        c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
results.RvR <- as.data.frame(results.RvR)
head(results.RvR)
results.RvR_GSEA <- results.RvR %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(results.RvR_GSEA)

#Convert Pennycress Gene IDs into Arabidopsis and remove NAs
geneList <- inner_join(results.RvR_GSEA, orthologs, by = "PennycressGeneID")
geneList <- na.omit(geneList)
head(geneList)
geneList <- geneList %>%
  dplyr::select(ArabidopsisGeneID, log2FoldChange)
colnames(geneList)[1] ="ID"
colnames(geneList)[2] ="logFC"
head(geneList)

#LFC shrinkage
resLFC.RvR_df <- as.data.frame(resLFC.RvR)
resLFC.RvR_df  <- setNames(cbind(rownames(resLFC.RvR_df ), resLFC.RvR_df, row.names = NULL),
                           c("PennycressGeneID", "baseMean", "log2FoldChange", "lfcSE", "svalue"))
resLFC.RvR_df  <- as.data.frame(resLFC.RvR_df)
head(resLFC.RvR_df )
resLFC.RvR_df_GSEA <- resLFC.RvR_df %>%
  dplyr::select(PennycressGeneID, log2FoldChange)
head(resLFC.RvR_df_GSEA)

geneList <- inner_join(resLFC.RvR_df_GSEA, orthologs, by = "PennycressGeneID")
class(geneList)
geneList <- na.omit(geneList)
head(geneList)
geneList <- geneList %>%
  dplyr::select(ArabidopsisGeneID, log2FoldChange)
colnames(geneList)[1] ="ID"
colnames(geneList)[2] ="logFC"
head(geneList)


circ <- circle_dat(GOplot_gsea, geneList)

GOCircle(circ)
reduced <- c("GO:0009535", "GO:0036294", "GO:0009409", "GO:0009642", "GO:0010200", "GO:0050832", "GO:0009753", "GO:0031347", "GO:0009408", "GO:0006091")
GOCircle(circ, nsub = reduced)
GOBubble(circ, labels = 1, ID=FALSE)
reduced_circ <- reduce_overlap(circ, overlap = 0.8)
png(file="~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_Bubble.png",
    width=1200, height=700)
#reduced_circ$jit<-with(reduced_circ, ifelse(term == "RNA modification" | term == "inorganic ion homeostasis", 0.5, 0.5))
GOBubble(reduced_circ, title = 'MN106 vs SP3210 Recovery Roots Enriched GO Terms - GSEA', labels = 10, ID=FALSE, table.legend = FALSE, colour = c("#009E73", "#E69F00", "#56B4E9")) +
  #geom_text(aes(label= term), size=5, position=position_jitter(width=reduced_circ$jit,height=reduced_circ$jit)) +
  geom_text_repel(aes(label = term), size = 5.7) +
  theme(axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size=20), legend.text = element_text(size=20), title = element_text(size=20))
ggsave("~/R_dir/Pennycress RNASeq/Multi-comparisons/Roots - Waterlogging & Recovery/GSEA results/RvR_Roots_gsea_Bubble_colors.jpg", width = 35, height = 23, unit = "cm", dpi = 300)
dev.off()

#We can run this to get a table in grey color for now
GOCircle(reduced_circ, nsub=20)

GOBubble(reduced_circ, labels = 4, ID=FALSE)

#Over-Representation Analysis with an OrgDb (Arabidopsis) ####

#set Arabidopsis organism here
organism = "org.At.tair.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#Input DEG list with corresponding Arabidopsis ortholog
gene <- read.csv("upregulated-MN106-recovery-roots.csv", header = TRUE)
gene <- inner_join(gene, orthologs, by = "PennycressGeneID")
gene <- gene %>%
  dplyr::select(ArabidopsisGeneID)
gene<- na.omit(gene)
head(gene)
class(gene)
gene <- as.vector(gene)
gene <- as.character(gene)

#Use enrichGO this time instead of enricher function
go_enrich <- enrichGO(gene = gene,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#graph the results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
dotplot(go_enrich)

linkage <- pairwise_termsim(go_enrich)
emapplot(linkage)
goplot(go_enrich, showCategory = 10)
#categorySize can be either 'pvalue' or 'geneNum'
#cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)


gene <- read.csv("upregulated-MN106-waterlogged-roots.csv", header = TRUE)
gene <- inner_join(gene, orthologs, by = "PennycressGeneID")
gene <- gene %>%
  dplyr::select(ArabidopsisGeneID)
gene<- na.omit(gene)
gene <- pull(gene)
head(gene)
class(gene)


#Use enrichGO this time instead of enricher function
go_enrich <- enrichGO(gene = gene,
                      #universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'TAIR',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#graph the results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
dotplot(go_enrich)

linkage <- pairwise_termsim(go_enrich)
emapplot(linkage)
goplot(go_enrich, showCategory = 10)
#categorySize can be either 'pvalue' or 'geneNum'
#cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)


