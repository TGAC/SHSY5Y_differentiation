#AUTHOR: dwright / haerty group / earlham institute
#EMAIL: wilfried.haerty@earlham.ac.uk
#DATE: 2020 12 03
#TITLE: DTE_analysis
#DESCRIPTION: Differential expression of differentiated and undifferentiated sh-sy5y human neuroblastoma cell lines using TALON annotation. Salmon quantification used (with 100 bootstraps) by counting transcriptome-mapped reads. Makes use of publicly available R packages, see relevant user manuals for details.    
# LICENCE: MIT

#Clear out R before starting anything
rm(list=ls())
#Set the working directory to my file
setwd("~/datafiles/DE")
getwd()

#Import libraries 
library(tximport) #import Salmon quantification data in appropriate format for EdgeR
library(plyr) #rename columns easily
library(edgeR) #DGE package
library(dplyr)

### DIY functions 

# 'Not in' function (a la grep -v)
'%!in%' <- function(x,y)!('%in%'(x,y))

##### 1
##### NANOPORE DATA. ISOFORM LEVEL QUANTIFICATION FROM SALMON USING TXIMPORT

### CATCHSALMON VERSION: DTE using Salmon counts

## 1) IMPORT WITH CATCHSALMON: use catchSalmon to read in the isoform count sample data directly from Salmon ready for EdgeR
iso.txi <- catchSalmon(c('D1','D2','D3','D4','D5','U1','U2','U3','U4','U5'), verbose = TRUE)

## 2) COUNT FILTERING: 

#A) remove from catchSalmon object the scaffold-only novel isoforms
scaffisos <- read.csv('novel_isoforms_scaffolds_toremove_N78_20200827.txt', sep = '\t', header = FALSE)$V1
scaffisos <- as.character(scaffisos)

#remove from counts
nrow(iso.txi$counts)
iso.txi$counts <- iso.txi$counts[!rownames(iso.txi$counts) %in% scaffisos, ]
nrow(iso.txi$counts)
#remove from annotation
nrow(iso.txi$annotation)
iso.txi$annotation <- iso.txi$annotation[!rownames(iso.txi$annotation) %in% scaffisos, ]
nrow(iso.txi$annotation)

#B) remove from catchSalmon object the zero-row counts (to do this - pre-calculate the overdispersion corrected counts) 
#create count as imported into DGElist and get 'zeros' from this so that read-filtering happens BEFORE creating DGElist...
iso.txi$countcheck <- iso.txi$counts/iso.txi$annotation$Overdispersion
nrow(iso.txi$countcheck)

# get list of isoforms with zero count (after overdispersion correction)
zerocount_isos <- iso.txi$counts[rowSums(iso.txi$counts[,1:ncol(iso.txi$counts)]) == 0,] # changed for counts or countcheck (overdispersion)
zeros_to_remove <- cbind(rownames(zerocount_isos))
zeros_to_remove <- as.vector(zeros_to_remove)

# remove zero-count rows from counts
nrow(iso.txi$counts)
iso.txi$counts <- iso.txi$counts[!rownames(iso.txi$counts) %in% zeros_to_remove, ]
nrow(iso.txi$counts)
#remove zero-count rows from annotation
#remove from annotation
nrow(iso.txi$annotation)
iso.txi$annotation <- iso.txi$annotation[!rownames(iso.txi$annotation) %in% zeros_to_remove, ]
nrow(iso.txi$annotation)

## 3) MAKE DGELIST from the filtered catchSalmon input (which uses bootstraps for overdispersion), add grouping data 
dge_catchSalmonVersion <- DGEList(counts=iso.txi$counts/iso.txi$annotation$Overdispersion, genes=iso.txi$annotation, group=rep(1:2,each=5)) #this takes the catchSalmon counts, adjusted with bootstrap-calculated overdispersion and makes DGEList with sample details
summary(dge_catchSalmonVersion)

## 4) add normalisation factors to dge-catchSalmonVersion object once filtered. catchSalmon imports lib sizes etc itself
dge_catchSalmonVersion <- calcNormFactors(dge_catchSalmonVersion, method="TMM")
plotMDS(dge_catchSalmonVersion)

# pseudoNormCountCheck -- think about this? Is this meaningful? 
eff.lib.size <- dge_catchSalmonVersion$samples$lib.size*dge_catchSalmonVersion$samples$norm.factors
eff.lib.size
normCounts <- cpm(dge_catchSalmonVersion)
pseudoNormCounts <- log2(normCounts + 1)
boxplot(pseudoNormCounts, col="gray", las=3, ylab = 'Pseudonormalised Counts')

# model creation for DTE
sampleType <- rep("D", ncol(dge_catchSalmonVersion)) #D = differentiated; U = undifferentiated
sampleType[grep("U", colnames(dge_catchSalmonVersion))] <- "U"
sampleType
sampleReplicate <- paste("S", rep(1:5, 2), sep="")
sampleReplicate
designMat <- model.matrix(~sampleReplicate + sampleType)
designMat

#Calculate dispersion statistics 
dge_catchSalmonVersion <- estimateGLMCommonDisp(dge_catchSalmonVersion, design=designMat)
dge_catchSalmonVersion <- estimateGLMTrendedDisp(dge_catchSalmonVersion, design=designMat)
dge_catchSalmonVersion <- estimateGLMTagwiseDisp(dge_catchSalmonVersion, design=designMat)
dge_catchSalmonVersion

#plot biological coefficient of variation 
plotBCV(dge_catchSalmonVersion)

#8 RUN DIFFERENTIAL EXPRESSION MODEL 
fit <- glmFit(dge_catchSalmonVersion, designMat)
fit
# check the column names of the model to make sure you get right column for 'coefficient' - the thing you're testing!
colnames(fit)
lrt <- glmLRT(fit, coef='sampleTypeU')
lrt
# check top results with topTags
edgeR_result <- topTags(lrt)
edgeR_result

# Smear plot results
deIsos <- decideTestsDGE(lrt, adjust.method = 'fdr', p=0.05)
deIsos <- rownames(lrt)[as.logical(deIsos)]
plotSmear(lrt, de.tags=deIsos)
abline(h=c(-1.5, 1.5), col=2)

lrt$table$logFC # list of log Fold-change from analysis
lrt$table$PValue # list of Pvalues from analysis 
lrt$table$FDR <- p.adjust(lrt$table$PValue, method = 'fdr') # manually correct for multiple testing 

# summary up and down regulated
summary(decideTests(lrt))

# Volcano plot results 
volcanoData <- cbind(lrt$table$logFC, -log10(lrt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)

# make file for GO enrichment analysis
DEout <- cbind(lrt$genes, lrt$table$logFC, lrt$table$PValue, lrt$table$FDR)
head(DEout)

#############################
##rename & rearrange columns for clarity 
#names(DEout) <- c('IsoformID','logFC','Pval','FDR','GeneID','TranscriptName','GeneName')
DEout <- cbind(rownames(DEout), data.frame(DEout, row.names=NULL))
names(DEout) <- c('transcript','length','effective_length','overdispersion','logFC','Pval','FDR')
head(DEout)
nrow(DEout)
### DE Main Output
write.table(DEout, 'cDNA_TALON_validated_isoform_diffexp_DEoutput_filter0reads_20201203.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
nrow(DEout)

### DE UPREG Threshold +1.5 & Threshold >0
# +1.5
upreg <- subset(DEout, DEout$logFC >= 1.5 & DEout$FDR < 0.05) 
nrow(upreg)
write.table(upreg, 'cDNA_TALON_validated_isoform_diffexp_UPREG_threshold1.5_filter0reads_20201203.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
# >0
upreg0 <- subset(DEout, DEout$logFC > 0 & DEout$FDR < 0.05) 
nrow(upreg0)
write.table(upreg0, 'cDNA_TALON_validated_isoform_diffexp_UPREG_threshold0_filter0reads_20201203.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)

### DE DOWNREG Threshold -1.5 & Threshold <0
# -1.5
dwnreg <- subset(DEout, DEout$logFC <= -1.5 & DEout$FDR < 0.05)
nrow(dwnreg)
write.table(dwnreg, 'cDNA_TALON_validated_isoform_diffexp_DWNREG_threshold1.5_filter0reads_20201203.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
# <0
dwnreg0 <- subset(DEout, DEout$logFC < 0 & DEout$FDR < 0.05)
nrow(dwnreg0)
write.table(dwnreg0, 'cDNA_TALON_isoform_diffexp_DWNREG_threshold0_filter0reads_20201203.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)

#print counts of each category
print(paste0("upreg 0 = ",nrow(upreg0)))
print(paste0("upreg 1.5 = ",nrow(upreg)))
print(paste0("downreg 0 = ",nrow(dwnreg0)))
print(paste0("downreg 1.5 = ",nrow(dwnreg)))

# Export the normalised counts (cpm(count) using TMM normalisation factors, see script above)
IsoformNormCounts = cbind(dge_catchSalmonVersion$genes, dge_catchSalmonVersion$counts)
IsoformNormCounts <- cbind(rownames(IsoformNormCounts), data.frame(IsoformNormCounts, row.names=NULL))
names(IsoformNormCounts) <- c('transcript','length','effective_length','overdispersion','D1','D2','D3','D4','D5','U1','U2','U3','U4','U5')
head(IsoformNormCounts)
write.table(IsoformNormCounts, 'Isoform_overdispersion-adjusted_counts_20201203.txt', sep = ' ', dec = '.', row.names = TRUE, col.names = TRUE, quote = FALSE)
