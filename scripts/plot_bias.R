########### 
# plot normalised read coverage for 10 read catasets 
# 2021 08 10
###########

###Clear out R before starting anything
rm(list=ls())
#Set the working directory to my file
setwd("~/datafiles2/SEQ_BIAS/")
getwd()

### Plot packages
library(ggplot2) #This is the released version

### 2021 07 26 both seq approaches, sequin ref_flat run: 
dtSRS <- read.table("noramlised_coverage_SRS_sequin_refflat_picardtools_20210726.csv", sep = ",", header = TRUE) 
dtONT <- read.table("noramlised_coverage_ONT_sequin_refflat_picardtools_20210726.csv", sep = ",", header = TRUE) 

### Mean & SD
line_types = c("short read"="longdash", "ONT"=2)
plot1 <- ggplot(data=dtSRS, aes(x=normalized_position)) + xlab("Normalised position") + ylab(expression(Mean~normalised~coverage~({"\u00b1"}~SD))) + theme_bw(base_size = 14) +
  # short reads
  geom_line(aes(y=mean, linetype = "short read")) + geom_ribbon(aes(ymin=SDlower, ymax=SDupper), alpha=0.3) +
  # ONT reads
  geom_line(data = dtONT, aes(y=meanONT, linetype = "ONT")) + geom_ribbon(aes(ymin=dtONT$SDlowerONT, ymax=dtONT$SDupperONT), alpha=0.5) + labs(linetype = "Sequencing")
plot1
