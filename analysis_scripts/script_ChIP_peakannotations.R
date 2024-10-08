# Annotation of ASCL1 peaks to proximal genes

#Load libraries
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(dplyr)

# Function to plot distance from TSS and region (promoter vs gene body , etc) of ASCL1 peaks
plot_ChIPseeker_multiple<-function(txdbX,filenames,plotname){
  peakFiles<-lapply(filenames,readPeakFile)
  peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdbX,tssRegion=c(-3000, 3000), verbose=FALSE)
  tiff(paste0('peakAnno_',plotname,'.tiff'),height = 8, width = 17.5, units='cm', compression = "lzw", res = 300)
  print(plotAnnoBar(peakAnnoList))
  dev.off()
  
  tiff(paste0('distToTSS_',plotname,'.tiff'),height = 5, width = 17.5, units='cm', compression = "lzw", res = 300)
  print(plotDistToTSS(peakAnnoList))
  dev.off()
}

########################## consensus peaks for each cell line, individually #########################################
# Run
my_list <- list('consensus_IMR.bed', 'consensus_BE.bed')
names(my_list) <- c("IMR", "BE")
plotname <- "IMRASCL1 and BEASCL1 anno"
txdb1<- TxDb.Hsapiens.UCSC.hg38.knownGene

plot_ChIPseeker_multiple(txdb1,my_list,plotname)

# Annotate peaks to nearest gene
peakFiles<-lapply(my_list,readPeakFile)
peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdb1,tssRegion=c(-3000, 3000), verbose=FALSE)

# Retrieve data for each cell line
IMR <- as.data.frame(peakAnnoList[1])
BE <- as.data.frame(peakAnnoList[2])

# Create a column with absolute distance values (for each cell line)
IMR$IMR.absolute_distToTSS <- ifelse(IMR$IMR.distanceToTSS<0, (-IMR$IMR.distanceToTSS), IMR$IMR.distanceToTSS)
BE$BE.absolute_distToTSS <- ifelse(BE$BE.distanceToTSS<0, (-BE$BE.distanceToTSS), BE$BE.distanceToTSS)

#Translate the EntrezID to Official Gene Symbol
x <- AnnotationDbi::select(org.Hs.eg.db, keys=IMR$IMR.geneId, 
            columns="SYMBOL", keytype="ENTREZID")
IMR$gene <- x$SYMBOL
y <- AnnotationDbi::select(org.Hs.eg.db, keys=BE$BE.geneId, 
            columns="SYMBOL", keytype="ENTREZID")
BE$gene <- y$SYMBOL

# save - this is all annotated peaks for each cell line
write.csv(IMR, "IMR_peakAnno_all.csv", row.names = FALSE)
write.csv(BE, "BE_peakAnno_all.csv", row.names = FALSE)

#Create a dataframe with a column listing all the genes associated to ASCL1 peaks in each cell line
#without duplicates
IMR_unique<-as.data.frame(unique(IMR$gene))
BE_unique<-as.data.frame(unique(BE$gene))
#save
write.csv(IMR_unique, "IMR_peakAnno_uniquelist_all.csv", row.names = FALSE)
write.csv(BE_unique, "BE_peakAnno_uniquelist_all.csv", row.names = FALSE)

#Shortlist to genes withing 50kb of the nearest peak
IMR_proxgenes<-subset(IMR,IMR.absolute_distToTSS<50000)
BE_proxgenes<-subset(BE,BE.absolute_distToTSS<50000)
#save
write.csv(IMR_proxgenes, "IMR_peakAnno_lessthan50kb.csv", row.names = FALSE)
write.csv(BE_proxgenes, "BE_peakAnno_lessthan50kb.csv", row.names = FALSE)

# Remove duplicates again to obtain a df of all proximal genes, up to 50kb
IMR_prox_unique<-as.data.frame(unique(IMR_proxgenes$gene))
BE_prox_unique<-as.data.frame(unique(BE_proxgenes$gene))
#save
write.csv(IMR_prox_unique, "IMR_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
write.csv(BE_prox_unique, "BE_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
###

# Prepare IMR_proxgenes and BE_proxgenes file for bed format
IMR_proxgenes_bed <- IMR_proxgenes[,c(1:3)]
colnames(IMR_proxgenes_bed) <- c("chrom", "start", "end")
write.table(IMR_proxgenes_bed ,"IMR_50kb_proxgenes_anno.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

BE_proxgenes_bed <- BE_proxgenes[,c(1:3)]
colnames(BE_proxgenes_bed) <- c("chrom", "start", "end")
write.table(BE_proxgenes_bed ,"BE_50kb_proxgenes_anno.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

#remove the input list of files for the function
rm(my_list)
rm(plotname)

######################### Differentially bound regions by ASCL1 between cell lines #########################################

# Run anno on differentially bound by ASCL1 sites between cell lines
my_list <- list('DiffBind_BEvsIMR_IMRup_FDR0.01_l2fc_1.txt','DiffBind_BEvsIMR_FDRabove0.1.txt','DiffBind_BEvsIMR_BEup_FDR0.01_l2fc_1.txt')
names(my_list) <- c("IMR_up","Common", "BE_up")
plotname <- "Differentially bound or common peaks anno"
txdb1<- TxDb.Hsapiens.UCSC.hg38.knownGene

plot_ChIPseeker_multiple(txdb1,my_list,plotname)

# Annotate peaks to nearest gene
peakFiles<-lapply(my_list,readPeakFile)
peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdb1,tssRegion=c(-3000, 3000), verbose=FALSE)

# Retrieve data for each cell line
IMR_up <- as.data.frame(peakAnnoList[1])
Common <- as.data.frame(peakAnnoList[2])
BE_up <- as.data.frame(peakAnnoList[3])

# Create a column with absolute distance values (for each cell line)
IMR_up$IMR_up.absolute_distToTSS <- abs(IMR_up$IMR_up.distanceToTSS)
Common$Common.absolute_distToTSS <- abs(Common$Common.distanceToTSS)
BE_up$BE_up.absolute_distToTSS <- abs(BE_up$BE_up.distanceToTSS)

#Translate the EntrezID to Official Gene Symbol
a <- AnnotationDbi::select(org.Hs.eg.db, keys=IMR_up$IMR_up.geneId, 
                           columns="SYMBOL", keytype="ENTREZID")
IMR_up$gene <- a$SYMBOL

b <- AnnotationDbi::select(org.Hs.eg.db, keys=Common$Common.geneId, 
                           columns="SYMBOL", keytype="ENTREZID")
Common$gene <- b$SYMBOL

c <- AnnotationDbi::select(org.Hs.eg.db, keys=BE_up$BE_up.geneId, 
                           columns="SYMBOL", keytype="ENTREZID")
BE_up$gene <- c$SYMBOL

# save - this is all annotated peaks for each cell line
write.csv(IMR_up, "IMR_up_peakAnno_all.csv", row.names = FALSE)
write.csv(Common, "Common_FDR0.1_peakAnno_all.csv", row.names = FALSE)
write.csv(BE_up, "BE_up_peakAnno_all.csv", row.names = FALSE)

#Create a dataframe with a column listing all the genes associated to ASCL1 peaks in each cell line
#without duplicates
IMR_up_unique<-as.data.frame(unique(IMR_up$gene))
Common_unique<-as.data.frame(unique(Common$gene))
BE_up_unique<-as.data.frame(unique(BE_up$gene))
#save
write.csv(IMR_up_unique, "IMR_up_peakAnno_uniquelist_all.csv", row.names = FALSE)
write.csv(Common_unique, "Common_FDR0.1_peakAnno_uniquelist_all.csv", row.names = FALSE)
write.csv(BE_up_unique, "BE_up_peakAnno_uniquelist_all.csv", row.names = FALSE)


# Prepare IMR_up, Common and BE_up file for bed format
IMR_up_bed <- IMR_up[,c(1:3)]
colnames(IMR_up_bed) <- c("chrom", "start", "end")
IMR_up_bed <- IMR_up_bed[order(IMR_up_bed$chrom, IMR_up_bed$start), ]
write.table(IMR_up_bed ,"IMR_up.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

Common_bed <- Common[,c(1:3)]
colnames(Common_bed) <- c("chrom", "start", "end")
Common_bed <- Common_bed[order(Common_bed$chrom, Common_bed$start), ]
write.table(Common_bed ,"Common_FDR0.1.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

BE_up_bed <- BE_up[,c(1:3)]
colnames(BE_up_bed) <- c("chrom", "start", "end")
BE_up_bed <- BE_up_bed[order(BE_up_bed$chrom, BE_up_bed$start), ]
write.table(BE_up_bed ,"BE_up.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

#Shortlist to genes withing 50kb of the nearest peak
IMR_up_proxgenes<-subset(IMR_up,IMR_up.absolute_distToTSS<50000)
Common_proxgenes<-subset(Common,Common.absolute_distToTSS<50000)
BE_up_proxgenes<-subset(BE_up,BE_up.absolute_distToTSS<50000)
#save
write.csv(IMR_up_proxgenes, "IMR_up_peakAnno_lessthan50kb.csv", row.names = FALSE)
write.csv(Common_proxgenes, "Common_peakAnno_lessthan50kb.csv", row.names = FALSE)
write.csv(BE_up_proxgenes, "BE_up_peakAnno_lessthan50kb.csv", row.names = FALSE)

# Remove duplicates again to obtain a df of all proximal genes, up to 50kb
IMR_up_prox_unique<-as.data.frame(unique(IMR_up_proxgenes$gene))
Common_prox_unique<-as.data.frame(unique(Common_proxgenes$gene))
BE_up_prox_unique<-as.data.frame(unique(BE_up_proxgenes$gene))
#save
write.csv(IMR_up_prox_unique, "IMR_up_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
write.csv(Common_prox_unique, "Common_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
write.csv(BE_up_prox_unique, "BE_up_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
###
#remove the input list of files for the function
rm(my_list)
rm(plotname)


########################## joint consensus peaks (for both IMR and BE) #########################################
# Run
my_list <- list('consensus_IMRandBE.bed')
names(my_list) <- c("consensus")
plotname <- "Consensus_IMRandBE"
txdb1<- TxDb.Hsapiens.UCSC.hg38.knownGene

plot_ChIPseeker_multiple(txdb1,my_list,plotname)

# Annotate peaks to nearest gene
peakFiles<-lapply(my_list,readPeakFile)
peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdb1,tssRegion=c(-3000, 3000), verbose=FALSE)

# Retrieve data 
consensus <- as.data.frame(peakAnnoList[1])

# Create a column with absolute distance values (for each cell line)
consensus$consensus.absolute_distToTSS <- abs(consensus$consensus.distanceToTSS)

#Translate the EntrezID to Official Gene Symbol
i <- AnnotationDbi::select(org.Hs.eg.db, keys=consensus$consensus.geneId, 
                           columns="SYMBOL", keytype="ENTREZID")
consensus$gene <- i$SYMBOL


# save - this is all annotated peaks
write.csv(consensus, "Consensus_IMRandBE_peakAnno_all.csv", row.names = FALSE)

#Create a dataframe with a column listing all the genes associated to ASCL1 peaks in both cell lines
#without duplicates
consensus_unique<-as.data.frame(unique(consensus$gene))
#save
write.csv(consensus_unique, "Consensus_IMRandBE_peakAnno_uniquelist_all.csv", row.names = FALSE)

#Shortlist to genes withing 50kb of the nearest peak
consensus_proxgenes<-subset(consensus,consensus.absolute_distToTSS<50000)
#save
write.csv(consensus_proxgenes, "Consensus_IMRandBE_peakAnno_lessthan50kb.csv", row.names = FALSE)

# Remove duplicates again to obtain a df of all proximal genes, up to 50kb
consensus_prox_unique<-as.data.frame(unique(consensus_proxgenes$gene))
#save
write.csv(consensus_prox_unique, "Consensus_IMRandBE_peakAnno_lessthan50kb_unique.csv", row.names = FALSE)
###

# Prepare files for bed format
consensus_proxgenes_bed <- consensus_proxgenes[,c(1:3)]
colnames(consensus_proxgenes_bed) <- c("chrom", "start", "end")
write.table(consensus_proxgenes_bed ,"Consensus_IMRandBE_50kb_proxgenes.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)

#remove the input list of files for the function
rm(my_list)
rm(plotname)
