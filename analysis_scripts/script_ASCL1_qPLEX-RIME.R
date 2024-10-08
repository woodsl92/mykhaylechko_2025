# Install CRAN packages
install.packages(c("ggplot2", "grid", "dplyr", "tidyverse", "ggrepel"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("qPLEXanalyzer", "UniProt.ws", "clusterProfiler", "org.Hs.eg.db"))

#Load libraries
library(ggplot2)
library(grid)
library(dplyr)
library(tidyverse)
suppressWarnings(library(qPLEXanalyzer))
library(UniProt.ws)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

#A metadata file containing "SampleName", "SampleGroup" "BioRep", "TechRep", "TMTLabel" was created.
#The peptide intensities txt file was loaded.

metadata <- read_csv("metadata.csv")
intensities <- read_tsv("Lumos_28042022_PR1541_EVP_ACN__PeptideGroups.txt") %>%
  dplyr::select(matches('Sequence|Modifications|Protein|Abundance'))

#Check for NA values in each column 
na_count <- sapply(intensities, function(y) sum(length(which(is.na(y)))))
na_count
#will give you a list with the counts for each column.

#We have NA values in the Protein Accessions, so one protein doesn't have an accession number. Need to remove this row.
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
#so here we check if we do have these and then remove them or any columns we don't need to tidy the data
# check if have duplicate column names in data
duplicated(rlang:::chr_unserialise_unicode(rlang:::names2(intensities)))
# only keep the columns we need
intensities <- intensities %>%
  dplyr::select(matches('Sequence|Modifications|Protein|BE|IMR'))
colnames(intensities)

#There will be some peptides that are identified but do not have TMT values, 
#either for all samples or some , e.g. IgG. These are the missing values 
#and are usually below 5%. We can clean these up by either adding a small value 
#to the IgG samples (these will then be removed as unspecific interactors 
#if they are higher in the IgG), or add a small value to all the values in the dataset 
#(after creating the MSnset_data)

#add something to the NA values in IgG so they have a value (added the smallest value in the data = 0.1)
intensities <- intensities %>%
  mutate(IMR_IgG_1 = replace_na(IMR_IgG_1, 0.1)) %>%
  mutate(IMR_IgG_2 = replace_na(IMR_IgG_2, 0.1)) %>%
  mutate(BE_IgG_1 = replace_na(BE_IgG_1, 0.1)) %>%
  mutate(BE_IgG_2 = replace_na(BE_IgG_2, 0.1))

#To create the object we need to continue the analysis using the code below. 
#Need to specify the file to use, the metadata, and which columns represent 
#the experimental values. As well as which have the sequences and accessions. 
#Here, we can remove any missing values, if needed.

# Create MSnSet Object
MSnset_data <- convertToMSnset(intensities,
                               metadata=metadata,
                               indExpData=c(9:22),
                               Sequences=1,
                               Accessions=5,
                               rmMissing = TRUE)

#Only keep unique peptides (that identify only one protein group/ protein)
MSnset_data <- MSnset_data[which(fData(MSnset_data)[,3]==1 &
                                   fData(MSnset_data)[,4]==1),]

# asign colours
colours <- c("#1A2C79","#797CD0","#E80566", "#FF95C4")
names_forcol <- c("BE", "BE_IgG", "IMR", "IMR_IgG")
col <- setNames(colours, names_forcol)


#Quality control plots.
#The intensityPlot function generates a peptide intensity distribution plot 
#that helps in identifying samples with outlier distributions. 

intensity_plot <- intensityPlot(
  MSnset_data,
  sampleColours = col,
  title = "Peptide intensity distribution",
  colourBy = "SampleGroup",
  transform = TRUE,
  xlab = "log2(intensity)"
)

# Save the plot
ggsave(
  filename = "peptide_intensity_distribution.png",  # Filename to save the plot
  plot = intensity_plot,                                      # Plot object
  width = 6,                                        # Width of the plot in inches
  height = 4,                                       # Height of the plot in inches
  dpi = 300                                         # Resolution (dots per inch)
)

#The intensities can also be viewed in the form of boxplots by intensityPlot. 
#The plot shows the distribution of peptides intensities for each sample.

intensity_box <- intensityBoxplot(
  MSnset_data,
  title = "Peptide intensity distribution",
  sampleColours = col,
  colourBy = "SampleGroup"
)

# Save the plot
ggsave(
  filename = "peptide_intensity_distribution_boxplot.png",  
  plot = intensity_box,                                      
  width = 6,                                        
  height = 4,                                       
  dpi = 300                                         
)

#PCA plot, before normalisation

pca<- pcaPlot(
  MSnset_data,
  omitIgG = FALSE,
  sampleColours = col,
  transform = TRUE,
  colourBy = "SampleGroup",
  title = "PCA - ASCL1 qPLEX-RIME BE vs IMR - before NORM",
  labelColumn = "SampleName",
  labelsize = 4,
  pointsize = 4,
  x.nudge = 4,
  x.PC = 1
)

ggsave(
  filename = "PCA_BEvsIMR_qPLEX-RIME_NOnorm.png",  
  plot = pca,                                      
  width = 6,                                        
  height = 4,                                       
  dpi = 300                                         
)

#Look at ASCL1 coverage in this experiment. This means how much of the 
#ASCL1 sequence is covered by the peptides detected in the experiment.

# Coverage Plot
#fastaFile character: fasta file of protein sequence

coverage <- coveragePlot(MSnset_data,
             ProteinID = "P50553", 
             ProteinName = "ASCL1",
             fastaFile = "protein.faa")

ggsave(
  filename = "ASCL1_coverage_plot_BEvsIMR_qPLEX-RIME.png",  
  plot = coverage,                                      
  width = 6,                                        
  height = 2,                                       
  dpi = 300                                         
)

#correlation plot
corrPlot(MSnset_data)
hierarchicalPlot(MSnset_data)


#######################################################################################################
####################################### normalisation of data #########################################
#######################################################################################################

#Information from the qPLEX-Analyzer package:

#Mean/median scaling normalizeScaling: In this normalization method 
#the central tendencies (mean or median) of the samples are aligned. 
#The central tendency for each sample is computed and log transformed. 
#A scaling factor is determined by subtracting from each central tendency 
#the mean of all the central tendencies. The raw intensities are then divided 
#by the scaling factor to get normalized ones.

#In qPLEX-RIME data, the IgG (or control samples) should be normalized 
#separately from the bait protein pull-down samples. As IgG samples 
#represent the low background intensity, their intensity distribution profile 
#is different from bait pull-downs. Hence, normalizing the two together would 
#result in over-correction of the IgG intensity resulting in inaccurate computation 
#of differences among groups. To this end we provide groupScaling, the additional 
#parameter groupingColumn defines a category for grouping the samples, scaling is 
#then carried out within each group independently.

#Normalise all ASCL1 samples together and IgG separately using median scaling. 
#For some reason, "median" was clashing with some other package? 
#so had to run this to add "stats::median".

MSnset_norm <- groupScaling(MSnset_data, 
                            scalingFunction= stats::median, 
                            groupingColumn = "Antibody")

# Re-make intensity plots after normalisation

norm_intensity <- intensityPlot(
  MSnset_norm,
  sampleColours = col,
  title = "Peptide intensity distribution - Normalised Data",
  colourBy = "SampleGroup",
  transform = TRUE,
  xlab = "log2(intensity)"
)

#Save the plot
ggsave(
  filename = "peptide_intensity_distribution_normalised.png",  # Filename to save the plot
  plot = norm_intensity,                                      # Plot object
  width = 6,                                        # Width of the plot in inches
  height = 4,                                       # Height of the plot in inches
  dpi = 300                                         # Resolution (dots per inch)
)

#boxplot
norm_intensity_box  <-intensityBoxplot(
  MSnset_norm,
  title = "Peptide intensity distribution",
  sampleColours = col,
  colourBy = "SampleGroup"
)
#Save the plot
ggsave(
  filename = "peptide_intensity_distribution_boxplot_normalised.png",  # Filename to save the plot
  plot = norm_intensity_box,                                      # Plot object
  width = 6,                                        # Width of the plot in inches
  height = 4,                                       # Height of the plot in inches
  dpi = 300                                         # Resolution (dots per inch)
)

#PCA plot of normalised data
norm_pca<- pcaPlot(
  MSnset_norm,
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

ggsave(
  filename = "PCA_BEvsIMR_qPLEX-RIME_norm.png",  
  plot = norm_pca,                                      
  width = 6,                                        
  height = 4,                                       
  dpi = 300                                         
)

## SUM PEPTIDE INTENSITIES TO GET PROTEIN INTENSITIES
#Create a dataframe with annotations of all the proteins in this data, 
#from the accessions. add gene name, descriptions and genesymbol

# UniProt.ws::select(the tax id, your proteins of interest, 
#what columns we want to retrieve, 
#the type of key used for our proteins"UniProtKB")

proteins <- unique(fData(MSnset_norm)$Accessions)
columns <- c("id", "protein_name", "gene_primary")
hs <- UniProt.ws::UniProt.ws(taxId = 9606)

human_anno_lm <- UniProt.ws::select(hs, proteins, columns,"UniProtKB") %>%
  dplyr::select(Accessions = "Entry", 
                Gene = "Entry.Name",
                Description = "Protein.names", 
                GeneSymbol = "Gene.Names..primary.") %>% 
  arrange(Accessions)

head(human_anno_lm)
#write
write.csv(human_anno_lm, "human_anno_lm.csv")

#read
human_anno_lm<-read.delim("human_anno_lm.csv", header=TRUE, sep=",")

#Obtain protein intensities by the sum of peptides.

MSnset_Pnorm <- summarizeIntensities(MSnset_norm, 
                                     summarizationFunction = sum,
                                     annotation = human_anno_lm)

#Check peptide and protein intensity in the samples of ASCL1 (bait)

ascl1_levels<-peptideIntensityPlot(MSnset_data,
                     combinedIntensities = MSnset_Pnorm,
                     ProteinID = "P50553", 
                     ProteinName = "ASCL1")

ggsave(
  filename = "ASCL1_levels_norm.png",  
  plot = ascl1_levels,                                      
  width = 6,                                        
  height = 4,                                       
  dpi = 300                                         
)

#Now compare ASCL1 BE vs IgG BE and ASCL1 IMR vs IgG IMR
contrasts_IgGcom <- c(BE_vs_BE_IgG = "BE - BE_IgG",
                      IMR_vs_IMR_IgG = "IMR - IMR_IgG")

diffstats_IgGcom <- computeDiffStats(MSnset_Pnorm,
                                     contrasts=contrasts_IgGcom)
diffexp_IgGcom <- list()
for (i in 1:length(contrasts_IgGcom))
  diffexp_IgGcom[[i]] <- getContrastResults(diffstats=diffstats_IgGcom,
                                            contrast=contrasts_IgGcom[i], 
                                            writeFile= TRUE)
#Using the graphs in the package:
#MA
for (i in 1:length(contrasts_IgGcom))
{
  toplot <- diffexp_IgGcom[[i]]$Accessions[diffexp_IgGcom[[i]]$adj.P.Val < 0.05]
  if(length(toplot) > 10)
    toplot <- toplot[1:10]
  print(maVolPlot(diffstats_IgGcom, contrast = contrasts_IgGcom[i], 
                  plotType="MA", 
                  title= contrasts_IgGcom[i], 
                  selectedGenes = toplot, 
                  fdrCutOff=0.01))
}
#Volcano
for (i in 1:length(contrasts_IgGcom))
{
  toplot <- diffexp_IgGcom[[i]]$Accessions[diffexp_IgGcom[[i]]$adj.P.Val < 0.05]
  if(length(toplot) > 10)
    toplot <- toplot[1:10]
  print(maVolPlot(diffstats_IgGcom, contrast = contrasts_IgGcom[i], 
                  plotType="Volcano", 
                  title= contrasts_IgGcom[i], 
                  selectedGenes = toplot, 
                  fdrCutOff=0.05,
                  lfcCutOff = 0.5))
}

#### read saved files and include significance groups and minus log10 adj p value column
#BE
BEvsIgG <- read.delim("BE_vs_BE_IgG.txt", header=TRUE)

BEvsIgG$minus_log10_adj.P.Val <- -log10(BEvsIgG$adj.P.Val)

BEvsIgG <- BEvsIgG %>%
  mutate(Significance=case_when(
    adj.P.Val <0.05 & log2FC>0.5  ~ "Enriched in ASCL1 pull-down",
    adj.P.Val <0.05 & log2FC <(-0.5) ~ "Enriched in IgG",
    TRUE ~ "Non-significant"
  ))

write.csv(BEvsIgG, "BEvsIgG.csv")
#same for IMR
IMRvsIgG <- read.delim("IMR_vs_IMR_IgG.txt", header=TRUE)

IMRvsIgG$minus_log10_adj.P.Val <- -log10(IMRvsIgG$adj.P.Val)

IMRvsIgG <- IMRvsIgG %>%
  mutate(Significance=case_when(
    adj.P.Val <0.05 & log2FC>0.5  ~ "Enriched in ASCL1 pull-down",
    adj.P.Val <0.05 & log2FC <(-0.5) ~ "Enriched in IgG",
    TRUE ~ "Non-significant"
  ))

write.csv(IMRvsIgG, "IMRvsIgG.csv")

### shortlist known interactors from literature
lit <- read.delim("ASCL1_interactors_lit.csv", header=TRUE, sep=",")
IMRvsIgG_lit <- merge(IMRvsIgG, lit, by="GeneSymbol")
BEvsIgG_lit <- merge(BEvsIgG, lit, by="GeneSymbol")

#plot ggplot

# IgG vs IMR - labeled wtih some known proteins
imr_igg<- ggplot(IMRvsIgG, aes(log2FC, minus_log10_adj.P.Val)) + 
  geom_point(data = IMRvsIgG, aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 1, stroke = 1, colour = "gray", fill = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +  # Increase right margin
  geom_point(data = IMRvsIgG %>%
               filter(log2FC>0.5) %>%
               filter(adj.P.Val <0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#E40E4F") +
  geom_point(data = IMRvsIgG %>%
               filter(log2FC<(-0.5)) %>%
               filter(adj.P.Val<0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#FF95C4") +
  geom_text_repel(data = IMRvsIgG_lit, aes(label= GeneSymbol), 
                  #hjust= 3, 
                  #vjust = -0.15, 
                  colour = "black", 
                  size = 5,
                  #,angle = 0, 
                  xlim = c(4,10), ylim = c(1,12),
                 # max.overlaps = Inf, min.segment.length = 0, box.padding = 1, force_pull   = 0, 
                 nudge_x = 4,
                 ) + 
  expand_limits(x = c(-6, 10)) +
  geom_vline(xintercept=c(-0.5,0.5), linetype = "dotted", 
             linewidth = 1) +
  #xlim(-6,6) +
  labs(x = "log2FC", y = "-log10(Adj.P.Val)", title = "ASCL1 IMR vs IgG IMR qPLEX-RIME")


# Render the plot
print(imr_igg)

#save
ggsave(
  filename = "IMRvsIgG_volcanoplot_litlabeled.png",  
  plot = imr_igg,                                      
  width = 9,                                        
  height = 7,                                       
  dpi = 300                                         
)

# IgG vs BE - labeled wtih some known proteins
be_igg<- ggplot(BEvsIgG, aes(log2FC, minus_log10_adj.P.Val)) + 
  geom_point(data = BEvsIgG, aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 1, stroke = 1, colour = "gray", fill = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +  # Increase right margin
  geom_point(data = BEvsIgG %>%
               filter(log2FC>0.5) %>%
               filter(adj.P.Val <0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#1A2C79") +
  geom_point(data = BEvsIgG %>%
               filter(log2FC<(-0.5)) %>%
               filter(adj.P.Val<0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#8ea6db") +
  geom_text_repel(data = BEvsIgG_lit, aes(label= GeneSymbol), 
                  #hjust= 3, 
                  #vjust = -0.15, 
                  colour = "black", 
                  size = 5,
                  #,angle = 0, 
                  xlim = c(4,10), ylim = c(1,12),
                  # max.overlaps = Inf, min.segment.length = 0, box.padding = 1, force_pull   = 0, 
                  nudge_x = 4,
  ) + 
  expand_limits(x = c(-6, 10)) +
  geom_vline(xintercept=c(-0.5,0.5), linetype = "dotted", 
             linewidth = 1) +
  #xlim(-6,6) +
  labs(x = "log2FC", y = "-log10(Adj.P.Val)", title = "ASCL1 BE vs IgG BE qPLEX-RIME")


# Render the plot
print(be_igg)

#save
ggsave(
  filename = "BEvsIgG_volcanoplot_litlabeled.png",  
  plot = be_igg,                                      
  width = 9,                                        
  height = 7,                                       
  dpi = 300                                         
)


######################### RE-analyse after removing IgG #############################

#Remove IgG columns and any peptides that are not enriched in the BE and IMR ASCL1 pull-downs

# Select only proteins with log2FC >0.5 and adj.p.val <0.05 in each line
ASCL1enriched <- lapply(diffexp_IgGcom, function(x) x$Accessions[which(x$log2FC > 0.5 & x$adj.P.Val <0.05)])
#Combine both lists to only keep proteins enriched over IgG in both lines
ASCL1enriched_combined <- Reduce(intersect,ASCL1enriched[1:2])

# Remove IgG samples
IgG_columns <- grep("IgG",colnames(MSnset_data))
MSnset_data_noIgG_columns <- MSnset_data[,-IgG_columns]

# keep only those peptides found as significant enriched in samples compared to IgG
tokeep <- fData(MSnset_data_noIgG_columns)$Accessions %in% ASCL1enriched_combined

MSnset_noIgG_data <- MSnset_data_noIgG_columns[tokeep,]

###
#Normalise all ASCL1 samples together using median scaling. 
MSnset_noIgG_data_norm <- groupScaling(MSnset_noIgG_data, 
                                       scalingFunction = stats::median, 
                                       groupingColumn = "Antibody")
#check normalised intensity plot
#set colours
colours2 <- c("#2C4789","#E40E4F")
names_forcol2 <- c("BE", "IMR")
col2 <- setNames(colours2, names_forcol2)
#
norm_no_igg <- intensityBoxplot(
  MSnset_noIgG_data_norm,
  title = "Peptide intensity distribution - Normalised",
  sampleColours = col2,
  colourBy = "SampleGroup"
)
#Save the plot
ggsave(
  filename = "peptide_intensity_distribution_boxplot_normalised_noIgG.png",  # Filename to save the plot
  plot = norm_no_igg,                                      # Plot object
  width = 6,                                        # Width of the plot in inches
  height = 4,                                       # Height of the plot in inches
  dpi = 300                                         # Resolution (dots per inch)
)


#Summarise peptide intensities into proteins
MSnset_noIgG_data_norm_P <- summarizeIntensities(MSnset_noIgG_data_norm, sum, human_anno_lm)

intensityPlot(
  MSnset_noIgG_data_norm_P,
  sampleColours = col2,
  title = "Protein intensity distribution - Normalised",
  colourBy = "SampleGroup",
  transform = TRUE,
  xlab = "log2(intensity)"
)

#ASCL1 levels, without IgG samples
ascl1_levels_noigg <- peptideIntensityPlot(MSnset_noIgG_data_norm,
                     combinedIntensities = MSnset_noIgG_data_norm_P,
                     ProteinID = "P50553", 
                     ProteinName = "ASCL1")

ggsave(
  filename = "ASCL1_levels_norm_noIgG.png",  
  plot = ascl1_levels_noigg,                                      
  width = 6,                                        
  height = 4,                                       
  dpi = 300                                         
)

############### Compare ASCL1 interactors in BE vs IMR cell lines overexpressing ASCL1 #####################

contrast_BEvsIMR <- c(BE_vs_IMR = "BE - IMR")
diffstats <- computeDiffStats(MSnSetObj=MSnset_noIgG_data_norm_P, contrasts=contrast_BEvsIMR)
diffexp <- getContrastResults(diffstats, 
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

write.csv(BEvsIMR, "BEvsIMR.csv")
#
BEvsIMR <- read.delim("BEvsIMR.csv", header=TRUE, sep=",")
########################## outside of R, filter the BEvsIMR proteins further ###########################

# First download information about all proteins in the datafram from UniProt,
#specifically the cellular location

#Remove proteins that are not localised in the nucleus (unless their location is unknown - keep these)

#read in the shortlisted proteins with nuclear or unknown cellular locations, from UniProt
nuc <- read.delim("NuclearProteins_UniProt.csv", header=TRUE, sep=",")

#merge this with BEvsIMR
BEvsIMR_nuc_filtered <- merge(BEvsIMR, nuc, by="Accessions")

#Then, download data from CRAPome database 
#This database lists the proteins pulled-down in negative controls of affinity purification experiments

#Download only mass spec affinity purification experiments performed on Human samples

#Remove proteins from the BEvsIMR dataframe that are detected in 50% or over 
#CRAPome listed experiments 


#load CRAPome filtered df
noCRAPome <- read.delim("CRAPome_filtered.csv", header=TRUE, sep=",")

#merge with BEvsIMR_nuc_filtered
#This will now be the fully filtered BEvsIMR protein list
BEvsIMR_filtered <- merge(BEvsIMR_nuc_filtered, noCRAPome, by="Accessions")

#tidy df
BEvsIMR_filtered <- BEvsIMR_filtered[,c(-2,-3)]

write.csv(BEvsIMR_filtered, "BEvsIMR_filtered.csv", row.names = FALSE)

#check how many in each group
count<- BEvsIMR_filtered %>%
  count(BEvsIMR_filtered$Type)
################################# PLOT BEvsIMR ###################################################
#PLOT
bevsimr<- ggplot(BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val)) + 
  geom_point(data = BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 1.5, stroke = 1, colour = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC>=1) %>%
               filter(adj.P.Val <0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = '#1A2C79') +
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC<=(-1)) %>%
               filter(adj.P.Val<0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#E40E4F") +
  geom_text_repel(data = BEvsIMR_filtered %>%
                    filter(log2FC>1), aes(label= GeneSymbol), hjust= -0.15, 
                  vjust = -0.15, colour = "black", size = 4, angle = 0,
                  fontface = 2, max.overlaps = Inf, min.segment.length = 0, box.padding = 0, force_pull   = 0) + 
  geom_text_repel(data = BEvsIMR_filtered %>%
                    filter(log2FC<(-1)), aes(label= GeneSymbol), hjust= -0.15, 
                  vjust = -0.15, colour = "black", size = 4, angle = 0,
                  fontface = 2,max.overlaps = Inf, min.segment.length = 0, box.padding = 0, force_pull   = 0) +
  geom_vline(xintercept=c(-1,1), linetype = "dotted", 
             linewidth = 1) +
  #xlim(-6,6) +
  labs(x = "log2FC", y = "-log10(Adj.P.Val)", title = "BE vs IMR ASCL1 qPLEX-RIME")

# plot
print(bevsimr)

#save
ggsave(
  filename = "BEvsIMR_volcano.png",  
  plot = bevsimr,                                      
  width = 9,                                        
  height = 7,                                       
  dpi = 300                                         
)
###############################################################################
#TRY OTHER PLOT SIZES, ETC
bevsimr <- ggplot(BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val)) + 
  geom_point(data = BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 1.5, stroke = 1, colour = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC >= 1) %>%
               filter(adj.P.Val < 0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = '#1A2C79') +
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC <= (-1)) %>%
               filter(adj.P.Val < 0.05), aes(log2FC, minus_log10_adj.P.Val), shape = 21, size = 3, stroke = 1, colour = "black", fill = "#E40E4F") +
  
  # Labels on the right side (log2FC >=1)
  geom_text_repel(data = BEvsIMR_filtered %>%
                    filter(log2FC >= 1), aes(label = GeneSymbol),
                  size = 4, fontface = 2, colour = "black",
                  box.padding = 0.5, point.padding = 0.5, 
                  xlim = c(0.1,4.8), 
              #    ylim = c(0,0),
              #    direction = "y",  # Keep labels on their side
                  hjust = 0, 
              #    nudge_x = 0.5,  # Nudge labels right of 0
                  max.overlaps = 100, 
                  min.segment.length = 0) +  
  
  # Labels on the left side (log2FC <= (-1))
  geom_text_repel(data = BEvsIMR_filtered %>%
                    filter(log2FC <= (-1)), aes(label = GeneSymbol),
                  size = 4, fontface = 2, colour = "black",
                  box.padding = 0.5, point.padding = 0.5, 
                 # xlim = c(-4,-0.1), 
                #  direction = "y",  # Keep labels on their side
                  hjust = 1, 
               #   nudge_x = -0.5,  # Nudge labels left of 0
                  max.overlaps = 100, 
                  min.segment.length = 0) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) +
  xlim(-3, 3) +  # Adjust x limits to provide more space
  labs(x = "log2FC", y = "-log10(Adj.P.Val)", title = "BE vs IMR ASCL1 qPLEX-RIME")

# Plot
print(bevsimr)

# Save
ggsave(
  filename = "BEvsIMR_volcano_final.png",  
  plot = bevsimr,                                      
  width = 9,                                        
  height = 7,                                       
  dpi = 300                                         
)
#####################################################################

#######################################################################################################
####################################### Gene Ontology Analysis ########################################
#######################################################################################################
data <- BEvsIMR_filtered
go<-compareCluster(GeneSymbol~Type, data=data[,c('GeneSymbol','Type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP",pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)

#
#plot result
plot <- dotplot(go, font.size = 15) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange')) +
  theme(
    legend.text = element_text(size = 16),       # Increase the size of the legend text
    legend.title = element_text(size = 18),      # Increase the size of the legend title
    legend.key.height = unit(1.5, "cm"),         # Increase the height of the legend keys
    legend.key.width = unit(1, "cm")             # Adjust the width of the legend keys if needed
  )
#save
ggsave(filename = "compareCluster_qPLEX-RIME_BEvsIMR_BP.png", plot = plot, width = 10.2, height = 10, units = "in", dpi = 300)
ggsave(filename = "compareCluster_qPLEX-RIME_BEvsIMR_BP.pdf", plot = plot, width = 8, height = 10, units = "in", dpi = 300)
write.csv(go,("compareCluster_qPLEX-RIME_BEvsIMR_BP.csv"))

########################## Plot only CDK and Cyclin proteins ############################
#read data
BEvsIMR_filtered <- read.delim("BEvsIMR_filtered.csv", header=TRUE, sep=",")
#plot ggplot
bevsimr_cdk<- ggplot(BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val)) + 
  geom_point(data = BEvsIMR_filtered, aes(log2FC, minus_log10_adj.P.Val), 
             shape = 20, size = 1, stroke = 0, colour = "gray") +
  theme_bw() +
  
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(1, 1.5, 1, 1, "cm")) +
  #enriched in BE points
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC>=1) %>%
               filter(adj.P.Val <0.05), 
             aes(log2FC, minus_log10_adj.P.Val), 
             shape = 20, size = 1, stroke = 0, 
             colour = '#1A2C79', fill = '#1A2C79') +
  #enriched in IMR points
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC<=(-1)) %>%
               filter(adj.P.Val<0.05), 
             aes(log2FC, minus_log10_adj.P.Val), 
             shape = 20, size = 1, stroke = 0, 
             colour = "#E40E4F", fill = "#E40E4F") +
  
  #Non-enriched CDKs or Cyclins detected in the RIME
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC<(1) & GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")|
                        log2FC>(-1) & GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
             #   filter(GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
             aes(log2FC, minus_log10_adj.P.Val), 
             shape = 21, size = 2, stroke = 1, 
             colour = "gray", fill = "gray") +
  
  #enriched in IMR points that are CDKs or Cyclins detected in the RIME
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC<=(-1)) %>%
               filter(adj.P.Val<0.05,GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")), 
               aes(log2FC, minus_log10_adj.P.Val), 
             shape = 21, size = 2, stroke = 1, 
             colour = "#E40E4F", fill = "#E40E4F") +
  
  #enriched in BE points that are CDKs or Cyclins detected in the RIME
  geom_point(data = BEvsIMR_filtered %>%
               filter(log2FC>=(1)) %>%
               filter(adj.P.Val<0.05,GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")), 
             #  filter(),
             aes(log2FC, minus_log10_adj.P.Val), 
             shape = 21, size = 2, stroke = 1, 
             colour = '#1A2C79', fill = '#1A2C79') +

  #label CDKs and CYclins only
  geom_text_repel(data = BEvsIMR_filtered %>%
                  filter(GeneSymbol %in% c("CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9", "CCNA2", "CCNF", "CCNK")),
                  aes(label= GeneSymbol), hjust= -0.15, 
                  vjust = -0.15, colour = "purple", 
                  size = 4.5, angle = 0,
                  fontface = 2, max.overlaps = Inf, 
                  min.segment.length = 0.2, box.padding = 0, 
                  force_pull   = 0, nudge_x = 0.4) + 
    geom_vline(xintercept=c(-1,1), linetype = "dotted", 
             linewidth = 1) +
  geom_hline(yintercept=c(1.3), linetype = "dotted", 
             linewidth = 1) +
  #xlim(-6,6) +
  labs(x = "log2FC", y = "-log10(Adj.P.Val)", title = "BE vs IMR ASCL1 qPLEX-RIME")

# plot
print(bevsimr_cdk)

#save
ggsave(
  filename = "BEvsIMR_CDKandCyclin_volcano.png",  
  plot = bevsimr_cdk,                                      
  width = 7,                                        
  height = 5,                                       
  dpi = 300                                         
)
ggsave(
  filename = "BEvsIMR_CDKandCyclin_volcano_1.png",  
  plot = bevsimr_cdk,                                      
  width = 5,                                        
  height = 4,                                       
  dpi = 300                                         
)
ggsave(
  filename = "BEvsIMR_CDKandCyclin_volcano_2.png",  
  plot = bevsimr_cdk,                                      
  width = 5,                                        
  height = 4,                                       
  dpi = 300                                         
)

