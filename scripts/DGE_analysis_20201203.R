#AUTHOR: dwright / haerty group / earlham institute
#EMAIL: wilfried.haerty@earlham.ac.uk
#DATE: 2020 12 03
#TITLE: DGE_analysis
#DESCRIPTION: Differential expression of differentiated and undifferentiated sh-sy5y human neuroblastoma cell lines using TALON annotation. Salmon quantification used (with 100 bootstraps) by counting transcriptome-mapped reads. Summarise to gene-level DGE. Makes use of publicly available R packages, see relevant user manuals for details.    
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
library(data.table)

### DIY functions 

# 'Not in' function (a la grep -v)
'%!in%' <- function(x,y)!('%in%'(x,y))

##### 1
##### NANOPORE DATA. ISOFORM LEVEL QUANTIFICATION FROM SALMON USING TXIMPORT

# Make an isoform-->gene table for summarising counts at gene level from isoform level quantification 
myannot <- annot2ids("~/datafiles/ANNOTATION/Talon_minID90_N5K3_purged_filtered_gencode_renamed_validated_NoSequins.gtf")
myannot_df <- as.data.frame(myannot)
write.table(myannot_df, 'tx2gene_conversion_20201203.txt', sep = '\t', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)

### TXIMPORT APPROACH: Perform DGE using Salmon transcript counts merged to gene level with tximport functionality. As using Nanopore with custom filtering of isoforms with zero-row counts (see methods), it is necessary to proceed with isoform-level importation to calculate which isoforms need removing, filtering the corresponding quant.sf files manually and then re-importing the filtered quant.sf transcript count files at gene-level using tximport. 

## 1) IMPORT WITH TXIMPORT: use catchSalmon to read in the isoform count sample data directly from Salmon ready for EdgeR. Find list of isoforms to remove due to zero-row counts (see methods)

# Import Quant.sf files using catchSalmon and calculate list of isoforms that need removing before gene-level merging of counts.
catchsalmon.iso.txi <- catchSalmon(c('D1','D2','D3','D4','D5','U1','U2','U3','U4','U5'), verbose = TRUE)
##A) remove from catchSalmon object the scaffold-only novel isoforms
scaffisos <- read.csv('novel_isoforms_scaffolds_toremove_N78_20200827.txt', sep = '\t', header = FALSE)$V1
scaffisos <- as.character(scaffisos)
#remove from counts
nrow(catchsalmon.iso.txi$counts)
catchsalmon.iso.txi$counts <- catchsalmon.iso.txi$counts[!rownames(catchsalmon.iso.txi$counts) %in% scaffisos, ]
nrow(catchsalmon.iso.txi$counts)
#remove from annotation
catchsalmon.iso.txi$annotation <- catchsalmon.iso.txi$annotation[!rownames(catchsalmon.iso.txi$annotation) %in% scaffisos, ]
nrow(catchsalmon.iso.txi$annotation)
##B) remove from catchSalmon object the zero-row counts
##create count as imported into DGElist and get 'zeros' from this so that read-filtering happens before creating DGElist - more appropriate with Nanopore read quantification
catchsalmon.iso.txi$countcheck <- catchsalmon.iso.txi$counts/catchsalmon.iso.txi$annotation$Overdispersion
nrow(catchsalmon.iso.txi$countcheck)
## get list of isoforms with zero-row counts
zerocount_isos <- catchsalmon.iso.txi$counts[rowSums(catchsalmon.iso.txi$counts[,1:ncol(catchsalmon.iso.txi$counts)]) == 0,] # no reads in any sample i.e. zero-row count
zeros_to_remove <- cbind(rownames(zerocount_isos))
zeros_to_remove <- as.vector(zeros_to_remove)
# Create list of zero-count isoforms that need to be removed from the quant.sf files manually/externally 
write.table(zeros_to_remove, 'zerosum_isoforms_to_remove_N110116.txt', row.names = FALSE, quote = FALSE, col.names = FALSE) # quote = FALSE VERY IMPORTANT HERE OR grep -Ff -v won't work!!!! 

##2) MANUALLY/EXTERNALLY FILTER THE QUANT.SF COUNT FILES, REMOVING ANY SCAFFOLD-ONLY ISOFORMS AND ISOFORMS WITH ZERO-COUNT ROWS (I.E. NO READS IN ANY SAMPLES)
# Note the zerosum_isoforms_to_remove_N110116.txt list was merged externally with novel_isoforms_scaffolds_toremove_N78_20200827.txt list to form a single file of isoforms to remove from the quant.sf files during the filtering step for simplicity. This final file of isoforms to be removed from the quant.sf files is named "isoforms_to_remove_ScaffoldLoci_and_DTEzerosumCounts_N110194_20201209.txt". Filtering was done externally by grepping out isoforms in the list with 'grep -Ff -v'. The filtered counts were then imported below WITHOUT bootstraps. Copies of the filtered counts are provided in ~/datafiles/DE/salmon_quants_filtered_no_bootstraps/   

# IMPORT THE MANUALLY FILTERED QUANT.SF FILES USING TXIMPORT AND GO STRAIGHT TO GENE-MERGING.PROCEED WITH USUAL DGE PIPELINE, NOTE THAT COUNTS ALREADY FILTERED
tx2gene <- read.table('tx2gene_conversion_20201203.txt', sep = '\t', dec = '.')
gen.tx <- tximport(c('salmon_quants_filtered_no_bootstraps/D1_quant.sf','salmon_quants_filtered_no_bootstraps/D2_quant.sf','salmon_quants_filtered_no_bootstraps/D3_quant.sf','salmon_quants_filtered_no_bootstraps/D4_quant.sf','salmon_quants_filtered_no_bootstraps/D5_quant.sf','salmon_quants_filtered_no_bootstraps/U1_quant.sf','salmon_quants_filtered_no_bootstraps/U2_quant.sf','salmon_quants_filtered_no_bootstraps/U3_quant.sf','salmon_quants_filtered_no_bootstraps/U4_quant.sf','salmon_quants_filtered_no_bootstraps/U5_quant.sf'), type = "salmon", txOut = FALSE, tx2gene = tx2gene)

## 3) MAKE DGELIST from the filtered catchSalmon input (which uses bootstraps for overdispersion), add grouping data as before 
dgeList_genes <- DGEList(counts=gen.tx$counts, group=rep(1:2,each=5), samples = c('D1','D2','D3','D4','D5','U1','U2','U3','U4','U5')) #this takes the catchSalmon counts, adjusted with bootstrap-calculated overdispersion and makes DGEList with sample details.  
summary(dgeList_genes)

## 4) add normalisation factors to dge-catchSalmonVersion object once filtered. catchSalmon imports lib sizes etc itself!
dgeList_genes <- calcNormFactors(dgeList_genes, method="TMM")
plotMDS(dgeList_genes, labels = dgeList_genes$samples$samples)

# pseudoNormCountCheck -- think about this? Is this meaningful? 
eff.lib.size <- dgeList_genes$samples$lib.size*dgeList_genes$samples$norm.factors
eff.lib.size
normCounts <- cpm(dgeList_genes)
pseudoNormCounts <- log2(normCounts + 1)
boxplot(pseudoNormCounts, col="gray", las=3, ylab = 'Pseudonormalised Counts')

# model creation for DTE
sampleType <- rep("D", ncol(dgeList_genes)) #D = differentiated; U = undifferentiated
sampleType[grep("U", dgeList_genes$samples$samples)] <- "U"
sampleType
sampleReplicate <- paste("S", rep(1:5, 2), sep="")
sampleReplicate
designMat <- model.matrix(~sampleReplicate + sampleType)
designMat

#Calculate dispersion statistics 
dgeList_genes <- estimateGLMCommonDisp(dgeList_genes, design=designMat)
dgeList_genes <- estimateGLMTrendedDisp(dgeList_genes, design=designMat)
dgeList_genes <- estimateGLMTagwiseDisp(dgeList_genes, design=designMat)
dgeList_genes

#plot biological coefficient of variation 
plotBCV(dgeList_genes)

#8 RUN DIFFERENTIAL EXPRESSION MODEL 
fit <- glmFit(dgeList_genes, designMat)
fit
# check the column names of the model to make sure you get right column for 'coefficient' - the thing you're testing!
colnames(fit)
lrt <- glmLRT(fit, coef='sampleTypeU')
lrt
# check top results with topTags
edgeR_result <- topTags(lrt)
edgeR_result

# Smear plot results
deGenes <- decideTestsDGE(lrt, adjust.method = 'fdr', p=0.05)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
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
DEout <- cbind(row.names(lrt$table), lrt$table$logFC, lrt$table$PValue, lrt$table$FDR)
head(DEout)
#rename for clarity & convert the numeric columns to number for filtering!!!
DEout <- data.frame(DEout)
names(DEout) <- c('gene','logFC','Pval','FDR')
DEout$logFC <- as.numeric(DEout$logFC)
DEout$Pval <- as.numeric(DEout$Pval)
DEout$FDR <- as.numeric(DEout$FDR)
head(DEout)
nrow(DEout)

### DE Main Output
write.table(DEout, 'cDNA_TALON_validated_gene_diffexp_DEoutput_FilterZeroSumReads_20201208.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
nrow(DEout)

#### There seem to be a lot of genes - so maybe look at threshold - what about FDR < 0.001 instead of 0.05? Precedant?

### DE UPREG Threshold +1.5 & Threshold >0
# +1.5
upreg <- subset(DEout, DEout$logFC >= 1.5 & DEout$FDR < 0.05) 
nrow(upreg)
write.table(upreg, 'cDNA_TALON_validated_gene_diffexp_UPREG_threshold1.5_filter0reads_20201209.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
# >0
upreg0 <- subset(DEout, DEout$logFC > 0 & DEout$FDR < 0.05) 
nrow(upreg0)
write.table(upreg0, 'cDNA_TALON_validated_gene_diffexp_UPREG_threshold0_filter0reads_20201209.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)

### DE DOWNREG Threshold -1.5 & Threshold <0
# -1.5
dwnreg <- subset(DEout, as.numeric(DEout$logFC) <= -1.5 & as.numeric(DEout$FDR) < 0.05)
nrow(dwnreg)
write.table(dwnreg, 'cDNA_TALON_validated_gene_diffexp_DWNREG_threshold1.5_filter0reads_20201209.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)
# <0
dwnreg0 <- subset(DEout, DEout$logFC < 0 & DEout$FDR < 0.05)
nrow(dwnreg0)
write.table(dwnreg0, 'cDNA_TALON_validated_gene_diffexp_DWNREG_threshold0_filter0reads_20201209.txt', sep = ' ', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)

print(paste0("upreg 0 = ",nrow(upreg0)))
print(paste0("upreg 1.5 = ",nrow(upreg)))
print(paste0("downreg 0 = ",nrow(dwnreg0)))
print(paste0("downreg 1.5 = ",nrow(dwnreg)))

# Export the normalised counts (cpm(count) using TMM normalisation factors, see script above)
GeneNormCounts = cbind(dgeList_genes$genes, pseudoNormCounts)
GeneNormCounts <- cbind(rownames(GeneNormCounts), data.frame(GeneNormCounts, row.names=NULL))
names(GeneNormCounts) <- c('gene','D1','D2','D3','D4','D5','U1','U2','U3','U4','U5')
head(GeneNormCounts)
write.table(GeneNormCounts, 'Gene_pseudonormalised_counts_20201209.txt', sep = ' ', dec = '.', row.names = TRUE, col.names = TRUE, quote = FALSE)

#############################

# Export the Normalised Counts for each gene: 
GeneTMMCPMCounts = cbind(dgeList_genes$genes, normCounts)
GeneTMMCPMCounts <- cbind(rownames(GeneTMMCPMCounts), data.frame(GeneTMMCPMCounts, row.names=NULL))
names(GeneTMMCPMCounts) <- c('gene','D1','D2','D3','D4','D5','U1','U2','U3','U4','U5')
head(GeneTMMCPMCounts)
write.table(GeneTMMCPMCounts, 'Gene_level_normalised_counts_CPMfromTMM_20210121.txt', sep = '\t', dec = '.', row.names = FALSE, col.names = TRUE, quote = FALSE)