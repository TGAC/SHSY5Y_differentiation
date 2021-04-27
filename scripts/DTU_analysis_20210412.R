#AUTHOR: dwright / haerty group / earlham institute
#EMAIL: wilfried.haerty@earlham.ac.uk
#DATE: 2021 04 12
#TITLE: DTU_analysis
#DESCRIPTION: Differential Transcript Usage (DTU) of differentiated and undifferentiated sh-sy5y human neuroblastoma cell lines using TALON re-annotated gene models. This script uses the IsoformSwitchAnalyzeR workflow by Vitting-Seerup & Sandelin Bioinformatics 35 (21): 4469–71, using functions and code available from https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html. See also Love et al 2018 F1000 for good overview of different approaches to DTU analysis    
# LICENCE: MIT

#Clear out R before starting
rm(list=ls())
#Set the working directory
setwd("~/datafiles/DTU/")
getwd()

#Import libraries 
library(plyr) #rename columns easily
library(tximport) #import and calculate input data correctly
library(IsoformSwitchAnalyzeR) # complete workflow for a DTU and functional consequences analysis by Vitting-Seerup & Sandelin
library(stats) #hypergeometric test of enrichment

### 1 Import data using ISAR (ensure correct format in dir - each quant.sf in own sub-dir, see Vitting-Seerup )
salmonQuant <- importIsoformExpression(parentDir = "salmon_quant_files_ONT/")

# check importation
head(salmonQuant$abundance, 5)
head(salmonQuant$counts, 5)

### 2 Create matrix of experimental design (see methods: differentiated vs undifferentiated cells)
sampleType <- rep("Differentiated", ncol(salmonQuant$abundance[-1])) #D = differentiated; U = undifferentiated
sampleType[grep("U", colnames(salmonQuant$abundance[-1]))] <- "Undifferentiated"
ProjectDesign <- data.frame(sampleID = colnames(salmonQuant$abundance)[-1], condition = sampleType) 
ProjectDesign$condition <- as.factor(ProjectDesign$condition)
ProjectDesign

### 3 Create an isoformswitchanalyzeRlist object with the other inputs (reference gtf and transcriptome .fa used for salmon quantification)
ONTSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = ProjectDesign,
  isoformExonAnnoation ="~/datafiles/ANNOTATION/Talon_minID90_N5K3_purged_filtered_gencode_renamed_validated_NoSequins.gtf",
  isoformNtFasta       ="~/datafiles/ANNOTATION/Talon_minID90_N5K3_purged_filtered_gencode_renamed_validated_NoSequins_transcriptome.fa",
  showProgress = TRUE)
ONTSwitchList

# extract the isoform IDs (and corresponding gene IDs) of the subset that survived filtering based on read counts (i.e. were assessed for DTU because at least 1 had reads!)
outfile <- cbind(ONTSwitchList$isoformFeatures$isoform_id, ONTSwitchList$isoformFeatures$gene_id)
write.table(outfile, 'assessed_isoformIDs_N98449_with_geneIDs_N32325.txt', sep = '\t', quote = FALSE)

### Run part 1 of analysis. This filters non-expressed isoforms, identifies isoform switches, annotates ORFs swtiches, and outputs both nuceotide and peptide fastas for part2 (function consequences). 
ONTSwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist = ONTSwitchList,
  pathToOutput = "~/datafiles/DTU//DTU_ISAR_out",
  switchTestMethod = "DEXSeq",
  outputSequences = FALSE, # already run, just rerun to get SwitchList results, normally set to TRUE for first time to output sequences
  prepareForWebServers = FALSE) # webservers option throws error

extractSwitchSummary(ONTSwitchList)

### Run external assessments of results and import for stage 2 of functional consequence prediction:
ONTSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = ONTSwitchList, 
  dIFcutoff                 = 0.1,   # Cutoff for defining switch size 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
  pathToCPATresultFile      = "DTU_ISAR_out/DTU_Dec2020_CPAT.txt",codingCutoff = 0.725,
  pathToPFAMresultFile      = "DTU_ISAR_out/DTU_Dec2020_Pfam.txt",
  pathToIUPred2AresultFile  = "DTU_ISAR_out/DTU_Dec2020_IUPRED2A.txt",
  pathToSignalPresultFile   = "DTU_ISAR_out/DTU_Dec2020_SignalP.txt",
  outputPlots               = FALSE,
  quiet = FALSE)

### Now explore the results. Two main post-result analysis streams: A) individual-gene-centric switching analysis and B) global genome-wide patterns

## A) individual gene level. E.g. extract top/most different/consequential switches 

# Extract top isoform switches (there are in total 104 here, so extract all rather than top #) 
gene_details <- extractTopSwitches(ONTSwitchList, filterForConsequences = TRUE, n=104)
write.table(gene_details, 'DTU_FunctionallyConsequentialSwitches_N104_20201201.txt', quote = FALSE)

#RBM5 - splicing factor identified by intersecting ENCODE 356 RBPs against list of DTUs (see methods)
switchPlot(ONTSwitchList, gene='RBM5', condition2 = 'Undifferentiated', condition1 = 'Differentiated')

## B) Global genome-wide approach to ID themes and patterns on genome-scale. There are 4 main result overviews for this:  

#1 genome-wide global summary statistics (two functions consequence summary and splicing summary)
extractSwitchSummary(
  ONTSwitchList,
  filterForConsequences = TRUE
) 

#2 genome-wide consequence summary & enrichment analysis
extractConsequenceSummary(
  ONTSwitchList,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)

conenr <- extractConsequenceEnrichment(
  ONTSwitchList,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = FALSE,
  returnResult = TRUE # if TRUE returns a data.frame with the results
)
write.table(conenr, 'DTU_ConsequenceEnrichment_20201202.txt', sep = '\t', quote = FALSE)

#3 genome-wide splicing enrichment analysis
extractSplicingSummary(
  ONTSwitchList,
  splicingToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)

splenr <- extractSplicingEnrichment(
  ONTSwitchList,
  returnResult = TRUE # if TRUE returns a data.frame with the results
)
write.table(conenr, 'DTU_SplicingEnrichment_20201202.txt', sep = '\t', quote = FALSE)

#4 Genome-wide changes in isoform usage - useful if expected differences are large
extractConsequenceGenomeWide(ONTSwitchList) 
extractSplicingGenomeWide(ONTSwitchList)

#######
# Hypergeometric test of enrichment of RBM5-binding between DTU-identified genes vs the background list 
#######

# for hypergeometric test we need: 
# m - number of marked elements (total linked to RBM5)
# N - total number of genes assessed (from processed gtf annotation) 
# n - this is N - m (i.e. proportion of non-marked elements, genes without RBM5 link)
# k - Size of selection (i.e. the DTU set of 104)
# x - the number marked in the selection (linked to RBM5 in the 104 DTUs)

#### RUNNING HYPERGEOMETRIC TEST ON INTERSECTS OF TWO CELLS LINES FOR CONFIDENCE IN RESULT
m2 <- 9882 # total linked to RBM5 interaction using ENCODE chipseq bedtools intersect from intersection of two cell line runs (ENCFF9xx and ENCFF1xx)
N2 <- 32325 # total assessed in the DTU analysis (from DTU ONTSwitchlist)
n2 = N2 - m2
k2 <- 104 # total found to be DTU 
x2 <- 69 # total of DTU found to be RBM5 associated

# Calculate Hypergeometric test p-value
intersects.p.value <-  phyper(q=x2-1, m=m2, n=n2, k=k2, lower.tail=FALSE)
intersects.p.value

fold.enrichment2 <- (x2 / k2) / (m2 / N2)
fold.enrichment2