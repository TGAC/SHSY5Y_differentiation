#AUTHOR: dwright / haerty group / earlham institute
#EMAIL: wilfried.haerty@earlham.ac.uk
#DATE: 2021 04 12
#TITLE: CACNA2D2_transcript_plotting
#DESCRIPTION: plot of CACNA2D2 transcripts showing additional exon identified by TALON
#LICENCE: MIT

#Clear out R 
rm(list=ls())
#Set the working directory
setwd("~/datafiles/CACNA2D2")
getwd()

# import libraries
library(IsoformSwitchAnalyzeR)
library(ggplot2)

# import gtf containing all gene information for CACNA2D2
CACNA2D2 <- importGTF("CACNA2D2_ENSG00000007402_entry.gtf")
head(CACNA2D2)

# Plot details of the gene
switchPlotTranscript(CACNA2D2, gene = 'CACNA2D2', rescaleTranscripts = TRUE, rectHegith = 0.2)

# ISOFORM LEVEL: read in the CPM TMM normalised counts for each isoform and plot:
isocounts3 <- read.csv("CACNA2D2_isoform_CPMcounts_20210121.txt", sep = '\t', header = TRUE)
# add a log2 +1 for plotting
isocounts3$logCPM1 <- log(isocounts3$CPM_TMM+1)

# plot
isocounts3$cell_state = factor(isocounts3$cell_state, levels=c("Undiff","Diff"))
ggplot(isocounts3, aes(x=isoform_id, y=logCPM1, fill=cell_state)) + 
  geom_boxplot() + theme_bw() + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + geom_point(aes(group=cell_state), position=position_jitterdodge(jitter.width=0.25, jitter.height = 0.0)) +theme(axis.text.x = element_text(angle = -40, hjust = -.005)) + xlab("Isoform ID") + ylab("log2(CPM+1)") + labs(fill=("Cell State"))

# GENE-LEVEL: plot
genecounts <- read.csv("CACNA2D2_gene_CPMcounts_20210121.txt", sep = '\t', header = TRUE)
# add a log(CPM+1) to improve visualisation of comparison between the cell states
genecounts$logCPM1 <- log2(genecounts$CPM_TMM+1)

genecounts$cell_state = factor(genecounts$cell_state, levels=c("Undifferentiated","Differentiated"))
ggplot(genecounts, aes(x=cell_state, y=logCPM1, fill = cell_state)) + 
  geom_boxplot() + theme_bw() + scale_fill_manual(values = c("#E69F00", "#56B4E9")) + geom_point(aes(group=cell_state), position=position_jitterdodge(jitter.width=0.1, jitter.height = 0.0)) + xlab("CACNA2D2 (ENSG00000007402)") + ylab("log(CPM+1)") + labs(fill=("Cell State")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + theme(axis.title.x = element_text(vjust=-0.9)) + theme(legend.position = "none")
