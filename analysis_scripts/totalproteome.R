##ANALYSE TOTAL PROTEOME DATA USING QPLEX-ANALYZER PACKAGE##
#load libraries
library(ggplot2)
library(grid)
library(dplyr)
library(tidyverse)
library(qPLEXanalyzer)
library(UniProt.ws)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)

#Load data
metadata <- read_csv("metadata.csv")
intensities <- read_tsv("Lumos_02082022_PR1569_EVP_FP__PeptideGroups.txt") %>%
  dplyr::select(matches('Sequence|Modifications|Protein|Abundance'))

#check for NA values in each column
na_count <- sapply(intensities, function(y) sum(length(which(is.na(y)))))
na_count

#We have NA values in the Protein Accessions, 
#so one protein doesn't have an accession number. 
#Need to remove this row.

# Remove row with no accession number
intensities <- intensities %>%  
  filter(!is.na(intensities$"Master Protein Accessions"))

#Replace the column name of the TMT values with the actual sample names

repNam <- pull(metadata, SampleName, TMTLabel)
head(repNam)

colnames(intensities)  <- colnames(intensities) %>%
  str_remove_all("Abundance: F1: |, .*") %>%
  str_replace_all(repNam)
colnames(intensities)

#For the next function to work no repeating column names should be present, 
#so here we check if we do have these and then remove them or any columns 
#we don't need to tidy the data

# check if have duplicate column names in data
#duplicated(rlang:::chr_unserialise_unicode(rlang:::names2(intensities)))
# only keep the columns we need
intensities <- intensities %>%
  dplyr::select(matches('Sequence|Modifications|Protein|BE|IMR'))
colnames(intensities)


#There will be some peptides that are identified but do not have TMT values, 
#either for all samples or some  
#These are the missing values and are usually below 5%. 
#We can clean these up by either adding a small value to the IgG samples 
#(these will then be removed as unspecific interactors if they are higher in the IgG), 
#or add a small value to all the values in the dataset (after creating the MSnset_data)

# #add something to the NA values so they have a low value in case it 
#is just one sample missing it

#replace NA with 0.1
intensities <- intensities %>%
  mutate(IMR1 = replace_na(IMR1, 0.1)) %>%
  mutate(IMR2 = replace_na(IMR2, 0.1)) %>%
  mutate(IMR3 = replace_na(IMR3, 0.1)) %>%
  mutate(BE1 = replace_na(BE1, 0.1)) %>%
  mutate(BE2 = replace_na(BE2, 0.1)) %>%
  mutate(BE3 = replace_na(BE3, 0.1))

# Create MSnSet Object
colnames(intensities)
MSnset_data <- convertToMSnset(intensities,
                               metadata=metadata,
                               indExpData=c(8:13),
                               Sequences=1,
                               Accessions=5,
                               rmMissing = TRUE)
#Only keep unique peptides (that identify only one protein group/ protein)
MSnset_data <- MSnset_data[which(fData(MSnset_data)[,3]==1 &
                                   fData(MSnset_data)[,4]==1),]

# asign colours
colours <- c("#1A2C79","#E40E4F")
names_forcol <- c("BE", "IMR")
col <- setNames(colours, names_forcol)

#### QC
intensityPlot(
  MSnset_data,
  sampleColours = col,
  title = "Peptide intensity distribution",
  colourBy = "SampleGroup",
  transform = TRUE,
  xlab = "log2(intensity)",
  trFunc = log2xplus1
)

intensityBoxplot(
  MSnset_data,
  title = "Peptide intensity distribution",
  sampleColours = col,
  colourBy = "SampleGroup"
)

pcaPlot(
  MSnset_data,
  omitIgG = FALSE,
  sampleColours = col,
  transform = TRUE,  
  colourBy = "SampleGroup",
  title = "PCA - ASCL1 qPLEX-RIME BE vs IMR",
  labelColumn = "SampleName",
  labelsize = 4,
  pointsize = 4,
  x.nudge = 4,
  x.PC = 1
)

coveragePlot(MSnset_data,
             ProteinID = "P50553", 
             ProteinName = "ASCL1",
             fastaFile = "protein.faa")


#NORMALISATION:
#Mean/median scaling normalizeScaling: In this normalization 
#method the central tendencies (mean or median) of the samples 
#are aligned. The central tendency for each sample is computed 
#and log transformed. A scaling factor is determined by subtracting 
#from each central tendency the mean of all the central tendencies. 
#The raw intensities are then divided by the scaling factor to get normalized ones.

#normalise
MSnset_norm <- groupScaling(MSnset_data, 
                            scalingFunction= stats::median)

#sum peptide intensities to protein
human_anno_lm <-read.delim("human_anno_lm.csv", header=TRUE, sep=",")
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, 
                                     summarizationFunction = sum,
                                     annotation = human_anno_lm)

#Perform BEvsIMR analysis.
#Before doing this, altered the code in computeDiffStats of qPLEX-Analyzer to

#fittedContrasts <- treat(contrastsfit, lfc = 1, trend = trend, robust = robust)
#return(diffstats <- list(MSnSetObj = MSnSetObj, 
#                         fittedLM = fit, 
#                         fittedContrasts = fittedContrasts))

#To calculate significance based on a log2 fold change threshold of 1.

contrast_BEvsIMR <- c(BE_vs_IMR = "BE - IMR")
diffstats_BEvsIMR <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrast_BEvsIMR)
diffexp <- getContrastResults(diffstats_BEvsIMR, 
                              contrast = contrast_BEvsIMR,
                              controlGroup = NULL,
                              writeFile= TRUE)
#
BEvsIMR <- read.delim("BE_vs_IMR.txt", header=TRUE)

BEvsIMR$minus_log10_adj.P.Val <- -log10(BEvsIMR$adj.P.Val)

BEvsIMR <- BEvsIMR %>%
  mutate(Type=case_when(
    adj.P.Val <0.05 & log2FC>=1  ~ "Enriched in BE",
    adj.P.Val <0.05 & log2FC <=(-1) ~ "Enriched in IMR",
    TRUE ~ "Non-enriched"
  ))

BEvsIMR <- BEvsIMR[,-2]
write.csv(BEvsIMR, "BEvsIMR_proteome.csv", row.names = FALSE)

#######################################################################################################
##################################### Combine with qPLEX-RIME ########################################
#######################################################################################################
#read
prot<- read.delim("BEvsIMR_proteome.csv", header=TRUE, sep=",")
rime <- read.delim("BEvsIMR_filtered.csv", header=TRUE, sep=",")

#check whether column names align between df and adjust
names(prot)
names(rime)
#before merging the datasets, make sure the log2FC and adj.p.val are labeled 
#as proteome, or qPLEX-RIME
prot <- prot %>%
  rename(
    log2FC_prot=log2FC,
    adj.P.Val_prot=adj.P.Val
  )
rime <- rime %>%
  rename(
    log2FC_rime=log2FC,
    adj.P.Val_rime=adj.P.Val
  )

prot_rime <- merge(rime[,c("GeneSymbol", "log2FC_rime", "adj.P.Val_rime")],prot[,c("GeneSymbol", "log2FC_prot", "adj.P.Val_prot")],  by="GeneSymbol")
prot_rime <- prot_rime %>%
  mutate(type_rime=case_when(
    log2FC_rime>=1 & adj.P.Val_rime<0.05 ~ "BE",
    log2FC_rime<=(-1) & adj.P.Val_rime<0.05 ~ "IMR",
    TRUE ~ "Common"
  ))

##set colours and sizes
sCol <- c("IMR" = "#E40E4F", 
          "BE" = "#1A2C79", 
          "Common" = "gray")
sSize <- c("IMR" = 2, 
           "BE" = 2, 
           "Common" = 2)
#plot
rimevsprot<- ggplot() + 
  geom_point(data=prot_rime,shape=20,stroke=1,
             aes(log2FC_rime, log2FC_prot, fill=type_rime, size=type_rime))+
  theme_bw()+
  
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +
  
  scale_fill_manual(values=sCol, name="") +
  scale_color_manual(values=sCol, name="") +
  scale_size_manual(values=sSize, name="") +
  
  geom_text_repel(data = prot_rime %>%
                    subset(prot_rime$GeneSymbol %in% 
                             c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
                  aes(x=log2FC_rime, y=log2FC_prot, label= GeneSymbol), 
                  hjust= -0.15, 
                  vjust = -0.15, colour = "purple", size = 4.5, angle = 0,
                  fontface = 2, family = "Karla-Bold", max.overlaps = Inf, 
                  min.segment.length = 0, box.padding = 0.3, force_pull   = 0) + 
  
  geom_vline(xintercept=c(-0.8,0.8), linetype = "dotted", 
             linewidth = 1) +
  geom_hline(yintercept=c(-1,1), linetype = "dotted", 
             linewidth = 1) +
  labs(x = "log2FC_qPLEX-RIME", y = "log2FC_proteome", title = "ASCL1 BEvsIMR qPLEX-RIME vs proteome")
# plot #2
rimevsprot<- ggplot() + 
  geom_point(data=prot_rime,shape=20,stroke=1, size=1,fill="gray", color="gray",
             aes(log2FC_rime, log2FC_prot, fill=type_rime, size=type_rime))+
  theme_bw()+
  
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +
  
  #Non-enriched CDKs or Cyclins detected in the RIME
  geom_point(data = prot_rime %>%
               filter(log2FC_rime<(1) & GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")|
                        log2FC_rime>(-1) & GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
              aes(log2FC_rime, log2FC_prot), 
             shape = 21, size = 1, stroke = 1, 
             colour = "gray", fill = "gray") +
  
  #enriched in IMR points that are CDKs or Cyclins detected in the RIME
  geom_point(data = prot_rime %>%
               filter(log2FC_rime<=(-1)),
             #& GeneSymbol %in% 
             #  c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")), 
             aes(log2FC_rime, log2FC_prot), 
             shape = 21, size = 2, stroke = 1, 
             colour = "#E40E4F", fill = "#E40E4F") +
  
  #enriched in BE points that are CDKs or Cyclins detected in the RIME
  geom_point(data = prot_rime %>%
               filter(log2FC_rime>=(1)),
                      #& GeneSymbol %in% 
      #        c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")), 
             aes(log2FC_rime, log2FC_prot), 
             shape = 21, size = 2, stroke = 1, 
             colour = '#1A2C79', fill = '#1A2C79') +
  
  geom_text_repel(data = prot_rime %>%
                    subset(prot_rime$GeneSymbol %in% 
                             c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
                  aes(x=log2FC_rime, y=log2FC_prot, label= GeneSymbol), 
                  hjust= -0.15, 
                  vjust = -0.15, colour = "purple", size = 5, angle = 0,
                  fontface = 2, family = "Karla-Bold", max.overlaps = Inf, 
                  min.segment.length = 0, box.padding = 0.3, force_pull   = 0) + 
  
  geom_vline(xintercept=c(-1,1), linetype = "dotted", 
             linewidth = 1) +
  geom_hline(yintercept=c(-1,1), linetype = "dotted", 
             linewidth = 1) +
  labs(x = "log2FC_qPLEX-RIME", y = "log2FC_proteome", title = "ASCL1 BEvsIMR qPLEX-RIME vs proteome")
#
ggsave(
  filename = "BEvsIMR_RIMEandPROT_3.png",  
  plot = rimevsprot,                                      
  width = 8,                                        
  height = 6,                                       
  dpi = 300                                         
)
