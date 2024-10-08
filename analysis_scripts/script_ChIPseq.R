#load libraries
library(DiffBind)
library(tidyverse)
library(rtracklayer)
library(GreyListChIP)
library(profileplyr)

#######################################################################################################
################################################# DiffBind ############################################
#######################################################################################################

# read in sample sheet 
samples <- read.csv("sample_sheet_LM_.csv")
## load data 
experiment <- dba(sampleSheet=samples)

## plot initial clustering 
plot(experiment)
#save
tiff(paste0('imr_BE_chip_initialcorrelation','.tiff'),height = 12, width = 12, units='cm', compression = "lzw", res = 300)
print(plot(experiment))
dev.off()

## retrieve black and grey list (instead of just subtracting the input from the samples, we calculate greylist regions from the input)
black <- dba.blacklist(experiment)
experiment.greylist <- dba.blacklist(black, Retrieve=DBA_GREYLIST)
# apply the black and grey lists
exp_filtered <- dba.blacklist(experiment, blacklist=DBA_BLACKLIST_HG38,
                              greylist=experiment.greylist)


# do occupancy analysis to choose best consensus peak
#BE
tiff(paste0('BE_chip_venn_repoverlap','.tiff'),height = 15, width = 19, units='cm', compression = "lzw", res = 300)
print(dba.plotVenn(exp_filtered, exp_filtered$masks$BE))
dev.off()
#IMR
tiff(paste0('imr_chip_venn_repoverlap','.tiff'),height = 15, width = 19, units='cm', compression = "lzw", res = 300)
print(dba.plotVenn(exp_filtered, exp_filtered$masks$IMR))
dev.off()

# create consensus peak of peaks in 3 out of the 4 replicates for each line
consensus_rep <- dba.peakset(exp_filtered, 
                             consensus=c(DBA_TISSUE),
                             minOverlap=0.75)

#create a new DBA object based on the two consensus peaks generated
exp_filtered_consensus <- dba(consensus_rep,
                              mask=consensus_rep$masks$Consensus,
                              minOverlap=1)
exp_filtered_consensus

#calculate an overall consensus peak to be using for counting
consensus_peaks <- dba.peakset(exp_filtered_consensus, bRetrieve=TRUE)

#export consensus for each line
consensus_BE <- dba.peakset(exp_filtered_consensus, 
                            exp_filtered_consensus$masks$BE,bRetrieve = TRUE,
                            writeFile = "consensus_BE.bed")

consensus_IMR <- dba.peakset(exp_filtered_consensus, 
                             exp_filtered_consensus$masks$IMR,bRetrieve = TRUE,
                             writeFile = "consensus_IMR.bed")

#plot Venn of the binding sites in IMR, BE and both
tiff(paste0('IMR_BE_chip_venn_overlap','.tiff'),height = 15, width = 19, units='cm', compression = "lzw", res = 300)
print(dba.plotVenn(consensus_rep, consensus_rep$masks$Consensus))
dev.off()


#count
exp_counted <- dba.count(exp_filtered, peaks=consensus_peaks, summits=250)
#
exp_counted


####################################### using the counted data ####################################################
#normalisation based on library size
norm_lib <- dba.normalize(exp_counted, normalize=DBA_NORM_LIB)

#contrast the normalised data
norm_lib_contrast <- dba.contrast(norm_lib,contrast=c("Tissue","BE","IMR"),
                                  reorderMeta = list(Tissue="BE"))

# run the analysis using DESEQ2
exp_analysed<- dba.analyze(norm_lib_contrast,method=DBA_DESEQ2)

##########################################PLOTS####################################

###############retrieve differentially bound sites
analysed.DB <- dba.report(exp_analysed)
analysed.DB

sum(analysed.DB$Fold>0)
sum(analysed.DB$Fold<0)

#MA plot
tiff(paste0('BEvsIMR_chip_counted_DESEq2_MAplot_LIBSIZE_norm','.tiff'),height = 12, width = 12, units='cm', compression = "lzw", res = 300)
dba.plotMA(exp_analysed, contrast=list(BE=exp_analysed$masks$BE), bNormalized=FALSE)
dev.off()


#PCA plot
tiff(paste0('BEvsIMR_chip_analysed_PCAplot_LIBSIZE_norm','.tiff'),height = 12, width = 12, units='cm', compression = "lzw", res = 300)
dba.plotPCA(exp_analysed,DBA_TISSUE,label=DBA_REPLICATE, vColors = c("#1A2C79", "#E40E4F"))
dev.off()

#volcano plots
tiff(paste0('BEvsIMR_chip_analysed_volcanoplot_LIBSIZE_norm','.tiff'),height = 12, width = 12, units='cm', compression = "lzw", res = 300)
dba.plotVolcano(exp_analysed)
dev.off()
#heatmap
hmap <- colorRampPalette(c('royalblue3', "ivory", 'darkorange'))(n = 13)
tiff(paste0('BEvsIMR_chip_analysed_heatmap_LIBSIZE_norm','.tiff'),height = 12, width = 12, units='cm', compression = "lzw", res = 300)
readscores <- dba.plotHeatmap(exp_analysed, contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap,colSideCols=c("#1A2C79","#E40E4F"))
dev.off()

#boxplot
my_colors <- c("#1A2C79", "#E40E4F")  # Specify your desired colors

tiff(paste0('BEvsIMR_chip_analysed_boxplot_LIBSIZE_norm_1','.tiff'),height = 14, width = 14, units='cm', compression = "lzw", res = 300)
dba.plotBox(exp_analysed, vColors = c("#1A2C79", "#E40E4F"))
dev.off()



#######################################SAVE DATA####################################

normCounts <- dba.peakset(norm_lib, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
write.csv(normCounts, file="IMR_BE_chip_normalized_libsize_counts.csv")

my.DB = dba.report(exp_analysed, th=1)
my.bd <- as.data.frame(my.DB)

write.csv(my.bd, "DiffBind_BEvsIMR_DESeq2_all.csv", row.names = FALSE)
write.table(my.bd, "DiffBind_BEvsIMR_DESeq2_all.txt", row.names = FALSE)
my.bd_ <- my.bd %>%
  subset(FDR > 0.1)
write.csv(my.bd_, "DiffBind_BEvsIMR_FDRabove0.1.csv", row.names = FALSE)
write.table(my.bd_, "DiffBind_BEvsIMR_FDRabove0.1.txt", row.names = FALSE)

dba.save(exp_analysed, "IMR_BE_experiment.analysis")
dba.save(exp_counted, "IMR_BE_experiment.counted")

my = dba.report(exp_analysed)
my <- as.data.frame(my)

my_BE <- my %>%
  subset(FDR < 0.01 & Fold>1)
my_IMR <- my %>%
  subset(FDR < 0.01 & Fold<(-1))

write.csv(my_BE, "DiffBind_BEvsIMR_BEup_FDR0.01_l2fc_1.csv", row.names = FALSE)
write.csv(my_IMR, "DiffBind_BEvsIMR_IMRup_FDR0.01_l2fc_1.csv", row.names = FALSE)

write.table(my_BE, "DiffBind_BEvsIMR_BEup_FDR0.01_l2fc_1.txt", row.names = FALSE)
write.table(my_IMR, "DiffBind_BEvsIMR_IMRup_FDR0.01_l2fc_1.txt", row.names = FALSE)
########################################################################################
