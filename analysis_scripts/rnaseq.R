#load libraries
library('ggplot2')
library('DESeq2')
library('apeglm')
library('tidyverse')
library('ggrepel')
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(Hmisc)
#######################################################################################################
###################################### Organise raw counts table ######################################
#######################################################################################################

# load raw counts
RNAseq_data<-read.delim('raw_counts.tsv',header=TRUE,stringsAsFactors = FALSE)
#view data
dim(RNAseq_data)
View(RNAseq_data)
names(RNAseq_data)
# select columns of interest (gene, BE and IMR cell lines)
RNAseq_data<-RNAseq_data[,grepl('gene|IMR|BE',names(RNAseq_data))==TRUE]

#load sample information sheet
sample_info<-read.csv('sample_information.csv',header=TRUE,stringsAsFactors = FALSE)
#factorise and reorder if needed 
sample_info$Cell_Line<-factor(sample_info$Cell_Line,levels=c('IMR','BE'))
sample_info$Condition<-factor(sample_info$Condition,levels=c('uninduced','induced'))
sample_info$Replicate<-factor(sample_info$Replicate)
#Create Sample_name column combining the information from the other columns (Cell_Line, Condition, Replicate)
sample_info$sample_name<-paste(sample_info$Cell_Line, sample_info$Condition,sample_info$Replicate,sep='_')

#mreorder the columns in the counts table
colnames(RNAseq_data)
RNAseq_data <- RNAseq_data[,c(1,13,15,17,19,21,12, 14, 16, 18, 20,3,5,7,9,11,2,4,6,8,10)]
#now that the order matches the order we want: IMR_uninduced, IMR_induced, BE_uninduced, BE_induced, 
#in the sample_info as well, we change the names
#to those in sample_info (tidier names)
names(RNAseq_data) <- c("gene", sample_info$sample_name)

#save the updated rawcounts file
write.csv(RNAseq_data, "rawcounts_updated.csv", row.names = FALSE)

#######################################################################################################
###################################### pre-filter the data ############################################
#######################################################################################################

count_df<-RNAseq_data# 1 = gene name column, rest = data columns of interest
names(count_df)
dim(count_df)
rownames(count_df)<-RNAseq_data$gene

#remove any gene with total counts being less than the number of samples (20)
count_df_filtered<-count_df[rowSums(count_df[,2:ncol(count_df)])>nrow(sample_info),]
#check new dimensions of df
dim(count_df_filtered)

#remove the column with gene names as now we have those as row names
count_df_filtered<-count_df_filtered[,-1]
dim(count_df_filtered)
#save
write.csv(count_df_filtered, file='filtered_counts.csv',row.names=TRUE,quote=FALSE)

#######################################################################################################
######################### Separate the cell lines in preparation for analysis ##########################
#######################################################################################################

#focus on one cell line at a time, so first select the information needed from the counts and sample info
#IMR-32
#select only IMR columns
count_df_IMR<-count_df_filtered[,grepl('IMR',names(count_df_filtered))==TRUE]
dim(count_df_IMR)
#create matrix
countsMatrix_IMR<-as.matrix(count_df_IMR)
#select only sample_info relating to IMR
sample_info_IMR<-subset(sample_info,(sample_info$sample_name %in% colnames(countsMatrix_IMR))==TRUE)
sample_info_IMR 
#save
write.csv(sample_info_IMR, file='sample_info_IMR.csv',row.names=FALSE,quote=FALSE)
write.csv(count_df_IMR, file='count_filtered_df_IMR.csv',row.names=TRUE,quote=FALSE)

#SK-N-BE(2)C
#select only BE columns
count_df_BE<-count_df_filtered[,grepl('BE',names(count_df_filtered))==TRUE]
dim(count_df_BE)
#create matrix
countsMatrix_BE<-as.matrix(count_df_BE)
#select only sample_info relating to IMR
sample_info_BE<-subset(sample_info,(sample_info$sample_name %in% colnames(countsMatrix_BE))==TRUE)
sample_info_BE 
#save
write.csv(sample_info_BE, file='sample_info_BE.csv',row.names=FALSE,quote=FALSE)
write.csv(count_df_BE, file='count_filtered_df_BE.csv',row.names=TRUE,quote=FALSE)


#######################################################################################################
###################################### Perform DESeq2 analysis ############################################
#######################################################################################################
#IMR-32
#use a model design of ~Replicate+Condition model
DESeq2data_IMR<-DESeqDataSetFromMatrix(countData = countsMatrix_IMR,colData = sample_info_IMR,design= ~Replicate + Condition)
DESeq2output_IMR <- DESeq(DESeq2data_IMR)
resultsNames(DESeq2output_IMR) 

#lfcShrink is a way of applying an adjustment to the log2FC to account for low counts likely to have more extreme log2FCs
DESeq2results_IMR <- lfcShrink(DESeq2output_IMR, coef='Condition_induced_vs_uninduced', type="apeglm") 

DESeq2results_IMR_ordered <- DESeq2results_IMR[order(DESeq2results_IMR$padj),]

#quick results overview
summary(DESeq2results_IMR_ordered)
sum(DESeq2results_IMR_ordered$padj < 0.05, na.rm=TRUE) 

#plot the spread of the results
png('MAplot_IMR_induced_vs_uninduced.png',width=450,height=400)
DESeq2::plotMA(DESeq2results_IMR , ylim=c(-4,4),alpha=0.05) #alpha = adjusted p-value threshold for colouring
dev.off()

#sort the results
DESeq2results_IMR_df<-as.data.frame(DESeq2results_IMR_ordered)
nrow(DESeq2results_IMR_df)
names(DESeq2results_IMR_df)
DESeq2results_IMR_df$gene<-rownames(DESeq2results_IMR_df)
DESeq2results_IMR_df<-DESeq2results_IMR_df[,c(6,1:5)]
View(DESeq2results_IMR_df)
#add fold change column
DESeq2results_IMR_df$FoldChange <- 2^(DESeq2results_IMR_df$log2FoldChange)

#save results as output table
write.csv(DESeq2results_IMR_df, file='IMR_induced_vs_uninduced_DESeq2output.csv',row.names=FALSE,quote=FALSE)
#IMR <- read.delim("IMR_induced_vs_uninduced_DESeq2output.csv", header=TRUE, sep=",")

#######################################################################################
#SK-N-BE(2)C
#use a model design of ~Replicate+Condition model
DESeq2data_BE<-DESeqDataSetFromMatrix(countData = countsMatrix_BE,colData = sample_info_BE,design= ~Replicate + Condition)
DESeq2output_BE <- DESeq(DESeq2data_BE)
resultsNames(DESeq2output_BE) #this gives you a list of relevant coefficients 

#lfcShrink is a way of applying an adjustment to the log2FC to account for low counts likely to have more extreme log2FCs
DESeq2results_BE <- lfcShrink(DESeq2output_BE, coef='Condition_induced_vs_uninduced', type="apeglm") 

DESeq2results_BE_ordered <- DESeq2results_BE[order(DESeq2results_BE$padj),]

#quick results overview
summary(DESeq2results_BE_ordered)
sum(DESeq2results_BE_ordered$padj < 0.05, na.rm=TRUE) 
#plot the spread of the results
png('MAplot_BE_induced_vs_uninduced.png',width=450,height=400)
DESeq2::plotMA(DESeq2results_BE , ylim=c(-4,4),alpha=0.05) #alpha = adjusted p-value threshold for colouring
dev.off()
#sort the results
DESeq2results_BE_df<-as.data.frame(DESeq2results_BE_ordered)
nrow(DESeq2results_BE_df)
names(DESeq2results_BE_df)
DESeq2results_BE_df$gene<-rownames(DESeq2results_BE_df)
DESeq2results_BE_df<-DESeq2results_BE_df[,c(6,1:5)]
View(DESeq2results_BE_df)
#add fold change column
DESeq2results_BE_df$FoldChange <- 2^(DESeq2results_BE_df$log2FoldChange)

#save results as output table
write.csv(DESeq2results_BE_df, file='BE_induced_vs_uninduced_DESeq2output.csv',row.names=FALSE,quote=FALSE)

#######################################################################################################
###################################### Visualisation of data ######################################
#######################################################################################################
#IMR-32 cell line

# Vst normalisation
vsd <- vst(DESeq2output_IMR, blind=FALSE) #blind=FALSE to prevent dispersions getting set too high for highly variable conditions
pca_data<-plotPCA(vsd, intgroup=c("Condition", "Replicate"),returnData=TRUE) 
# extract the % variance explained by PCA 1 and PCA 2 from pca_data, turn it to percentage and round in up
attributes(pca_data)
percentVar <- round(100 * attr(pca_data,which="percentVar")) 
#BE colours
data_colours<-c("#FF95C4","#E40E4F") 
#plot
f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = Condition,col = Condition, shape = Replicate)) + geom_point(size =3) + 
  scale_fill_manual(values=data_colours) + scale_colour_manual(values=data_colours) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(plot=f1,file=paste0("IMR_PCA_vst_norm_modelRandC.png"),width=5,height=4)
#

### Prepare data for plotting ####

DESeq2results_IMR_df_ <- DESeq2results_IMR_df %>%
  drop_na() %>% #drop any row that has NA
  mutate(diff_type=case_when( #create the categorical column. FIrst pick a name,then use tidyverse case_when function
    padj<0.05 & log2FoldChange>0.5 ~ "up",
    padj<0.05 & log2FoldChange<(-0.5) ~ "down",
    TRUE ~ "ns"
  )
  )
#save
write.csv(DESeq2results_IMR_df_,"DESeq2results_IMR_with_difftype.csv",row.names=FALSE,quote=FALSE)

#count how many in each group
table_IMR <- DESeq2results_IMR_df_ %>%
  count(diff_type) #get a count of that data
# We expect a similar number between up and down as that is the assumption of this analysis
write.csv(table_IMR,"DESeq2results_IMR_df_count.csv")

## volcano plot

#label only ASCL1
f1 <- ggplot(data=DESeq2results_IMR_df,mapping=aes(x=log2FoldChange, y=-log(padj),col=diff_type))+
  geom_point()+theme_classic()+scale_colour_manual(values=c("#FF95C4",'grey68',"#E40E4F"), name=NULL)+
  #coord_trans(ylim=c(0,10000))+
  geom_text_repel(data=DESeq2results_IMR_df %>%
                    filter(gene=="ASCL1"), aes(label=gene),size=8, col="black")+
theme(
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14)   
)
#save plot
ggsave(filename = "IMR_induced_vs_uninduced_ggplot_paper.png", plot = f1, width = 6, height = 4, units = "in", dpi = 300)

##########################################################################
# SK-N-BE(2)C
# Vst normalisation
vsd <- vst(DESeq2output_BE, blind=FALSE) #blind=FALSE to prevent dispersions getting set too high for highly variable conditions
pca_data<-plotPCA(vsd, intgroup=c("Condition", "Replicate"),returnData=TRUE) 
# extract the % variance explained by PCA 1 and PCA 2 from pca_data, turn it to percentage and round in up
attributes(pca_data)
percentVar <- round(100 * attr(pca_data,which="percentVar")) 
#BE colours
data_colours<-c("#8ea6db","#1A2C79") 
#plot
f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = Condition,col = Condition, shape = Replicate)) + geom_point(size =3) + 
  scale_fill_manual(values=data_colours) + scale_colour_manual(values=data_colours) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(plot=f1,file=paste0("BE_PCA_vst_norm_modelRandC.png"),width=5,height=4)


### Prepare data for plotting ####

DESeq2results_BE_df_ <- DESeq2results_BE_df %>%
  drop_na() %>% #drop any row that has NA
  mutate(diff_type=case_when( #create the categorical column. First pick a name,then use tidyverse case_when function
    padj<0.05 & log2FoldChange>0.5 ~ "up",
    padj<0.05 & log2FoldChange<(-0.5) ~ "down",
    TRUE ~ "ns"
  )
  )
#save
write.csv(DESeq2results_BE_df_, "DESeq2results_BE_with_difftype.csv", quote=FALSE, row.names = FALSE)
#count how many in each group
table_BE <- DESeq2results_BE_df_ %>%
  count(diff_type) #get a count of that data
# We expect a similar number between up and down as that is the assumption of this analysis
write.csv(table_BE,"DESeq2results_BE_df_count.csv")

# volcano plot
#label only ASCL1
b1 <- ggplot(data=DESeq2results_BE_df,mapping=aes(x=log2FoldChange, y=-log(padj),col=diff_type))+
  geom_point()+theme_classic()+
  scale_colour_manual(values=c("#8ea6db",'grey68',"#1A2C79"), name=NULL)+
  geom_text_repel(data=DESeq2results_BE_df %>%
                    filter(gene=="ASCL1"), aes(label=gene),size=8, col="black")+
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)   
  )

#save
ggsave(filename = "BE_induced_vs_uninduced_ggplot_paper.png", plot = b1, width = 6, height = 4, units = "in", dpi = 300)


#######################################################################################################
###################################### Gene Ontology Analysis ######################################
#######################################################################################################
#IMR-32 cell line
IMR <- read.delim("DESeq2results_IMR_with_difftype.csv", header=TRUE, sep=",")
#filter to keep top 1500 up / top 1500 down genes only, based on log2FC
IMR_top <- rbind(
  IMR %>% filter(diff_type == "up") %>% slice_max(order_by = log2FoldChange,n=1500),
  IMR %>% filter(diff_type == "down") %>% slice_min(order_by = log2FoldChange,n=1500)
)
#run compareCluster
go<-compareCluster(gene~diff_type, data=IMR_top[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='BP',pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
#plot result
plot<- dotplot(go, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange'))
#save
ggsave(filename = "compareCluster_IMR_top1500_BP.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
write.csv(go,("compareCluster_IMR_top1500_BP.csv"))

##################################################################
#SK-N-BE(2)C cell line
BE <- read.delim("DESeq2results_BE_with_difftype.csv", header=TRUE, sep=",")
#filter to keep top 1500 up / top 1500 down genes only, based on log2FC
BE_top <- rbind(
  BE %>% filter(diff_type == "up") %>% slice_max(order_by = log2FoldChange,n=1500),
  BE %>% filter(diff_type == "down") %>% slice_min(order_by = log2FoldChange,n=1500)
)
#run compareCluster
go_be<-compareCluster(gene~diff_type, data=BE_top[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='BP',pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
#plot result
plot_be<- dotplot(go_be, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange'))
#save
ggsave(filename = "compareCluster_BE_top1500_BP.png", plot = plot_be, width = 8, height = 6, units = "in", dpi = 300)
write.csv(go_be,("compareCluster_BE_top1500_BP.csv"))

#Filter all significantly changing genes
BE_sig <- rbind(
  BE %>% filter(diff_type == "up"),
  BE %>% filter(diff_type == "down")
)
#run compareCluster
go_sig<-compareCluster(gene~diff_type, data=BE_sig[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='BP',pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
#plot
plot_be_sig<-dotplot(go_sig, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange'))
#save
ggsave(filename = "compareCluster_BE_allsig_BP.png", plot = plot_be_sig, width = 8, height = 6, units = "in", dpi = 300)
write.csv(go_sig,("compareCluster_BE_allsig_BP.csv"))

#######################################################################################################
############ Gene Ontology Analysis and Heatmap of significant genes from both cell lines ##################
#######################################################################################################

################################################# Prepare data ####################################

####################
##calculate counts per million (CPM) for each cell line

#read filtered counts file
RNAseq_data <- read.delim("filtered_counts.csv", header=TRUE, sep=",")
#make sure column title for gene names is "gene"
RNAseq_data <- RNAseq_data%>%
  rename(
    gene=X
  )
#calculate CPM
RNAseq_CPM <- RNAseq_data

RNAseq_CPM[,2:ncol(RNAseq_CPM)] <-apply(RNAseq_CPM[,2:ncol(RNAseq_CPM)],2,
                                        function(x) (x/sum(x,na.rm = TRUE))*1000000)
#check it worked
colSums(RNAseq_CPM[,2:ncol(RNAseq_CPM)])
colSums(RNAseq_data[,2:ncol(RNAseq_data)])
#also make CPM files for each cell line
RNAseq_CPM_IMR <- RNAseq_CPM[,grepl('gene|IMR',names(RNAseq_CPM))==TRUE]
RNAseq_CPM_BE <- RNAseq_CPM[,grepl('gene|BE',names(RNAseq_CPM))==TRUE]

#save
write.csv(RNAseq_CPM,"RNAseq_CPM.csv", row.names = FALSE)
write.csv(RNAseq_CPM_IMR,"RNAseq_CPM_IMR.csv", row.names = FALSE)
write.csv(RNAseq_CPM_BE,"RNAseq_CPM_BE.csv", row.names = FALSE)

#########################
# Shortlist to only the set of genes significantly changing
#IMR-32
#read DESeq2 results dataframe with annotated significantly changing genes
IMR<- read.delim("DESeq2results_IMR_with_difftype.csv", header=TRUE, sep=",")
#filter for significant only
IMR_sig <- rbind(
  IMR %>% filter(diff_type == "up"),
  IMR %>% filter(diff_type == "down")
)
##########################################################
#SK-N-BE(2)C
#read DESeq2 results dataframe with annotated significantly changing genes
BE<- read.delim("DESeq2results_BE_with_difftype.csv", header=TRUE, sep=",")
#filter for significant only
BE_sig <- rbind(
  BE %>% filter(diff_type == "up"),
  BE %>% filter(diff_type == "down")
)
# join these two gene lists fully together
colnames(IMR_sig)

sig_set_all <- full_join(IMR_sig[,c("gene", "log2FoldChange", "lfcSE", "padj","FoldChange", "diff_type")], BE_sig[,c("gene","log2FoldChange", "lfcSE", "padj","FoldChange", "diff_type")], by="gene")
sig_set_all <- sig_set_all %>%
  dplyr::rename(
    log2FoldChange_IMR=log2FoldChange.x,
    log2FoldChange_BE=log2FoldChange.y,
    lfcSE_IMR = lfcSE.x,
    lfcSE_BE = lfcSE.y,
    padj_IMR = padj.x,
    padj_BE = padj.y,
    FoldChange_IMR = FoldChange.x,
    FoldChange_BE = FoldChange.y,
    diff_type_IMR = diff_type.x,
    diff_type_BE = diff_type.y
  )
# Create sugroups
sig_set_all <- sig_set_all %>%
  mutate(type=case_when( #create the categorical column. FIrst pick a name,then use tidyverse case_when function
    diff_type_IMR=="up" & diff_type_BE=="up" ~ "Both up",
    diff_type_IMR=="down" & diff_type_BE=="down" ~ "Both down",
    diff_type_IMR=="up" & is.na(diff_type_BE) ~ "IMR up",
    is.na(diff_type_IMR) & diff_type_BE=="up" ~ "BE up",
    diff_type_IMR=="down" & is.na(diff_type_BE) ~ "IMR down",
    is.na(diff_type_IMR) & diff_type_BE=="down" ~ "BE down",
    diff_type_IMR=="up" & diff_type_BE=="down" ~ "IMR up/BE down",
    diff_type_IMR=="down" & diff_type_BE=="up" ~ "IMR down/BE up",
    TRUE ~ "ns"
  )
  )

#Shortlist the CPM file to significant genes too

RNAseq_CPM_sig<- merge(RNAseq_CPM,sig_set_all[,c("gene","type")], by="gene")
RNAseq_CPM_IMR_sig<- merge(RNAseq_CPM_IMR,IMR_sig[,c("gene","diff_type")], by="gene")
RNAseq_CPM_BE_sig<- merge(RNAseq_CPM_BE,BE_sig[,c("gene","diff_type")], by="gene")

#save all
write.csv(IMR_sig,"sig_set_IMR.csv", row.names = FALSE)
write.csv(BE_sig,"sig_set_BE.csv", row.names = FALSE)
write.csv(sig_set_all,"sig_set_IMRandBE.csv", row.names = FALSE)

write.csv(RNAseq_CPM_IMR_sig,"CPM_sig_set_IMR.csv", row.names = FALSE)
write.csv(RNAseq_CPM_BE_sig,"CPM_sig_set_BE.csv", row.names = FALSE)
write.csv(RNAseq_CPM_sig,"CPM_sig_set_IMRandBE.csv", row.names = FALSE)

#############
#read DESeq2 dataframe with significantly changing genes in either cell line
sig_set_all <- read.delim("sig_set_IMRandBE.csv", header=TRUE, sep=",")
#check number of genes in each group
table_sig_all <- sig_set_all %>%
  count(type) #get a count of that data
# save
write.csv(table_sig_all,"DESeq2_sig_all_count.csv")

################################################# GENE ONTOLOGY ANALYSIS ####################################
#gene ontology analysis

#run compareCluster
go_both<-compareCluster(gene~type, data=sig_set_all[,c('gene','type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='BP',pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
#plot result
plot_both <- dotplot(go_both, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange')) +
  theme(
    legend.text = element_text(size = 16),       # Increase the size of the legend text
    legend.title = element_text(size = 18),      # Increase the size of the legend title
    legend.key.height = unit(1.5, "cm"),         # Increase the height of the legend keys
    legend.key.width = unit(1, "cm")             # Adjust the width of the legend keys if needed
  )
#save
ggsave(filename = "compareCluster_groups_bothlines_sigset_BP.png", plot = plot_both, width = 10, height = 14, units = "in", dpi = 300)
ggsave(filename = "compareCluster_groups_bothlines_sigset_BP.pdf", plot = plot_both, width = 10, height = 14, units = "in", dpi = 300)
write.csv(go_both,("compareCluster_groups_bothlines_sigset_BP.csv"))

################################################# HEATMAP ####################################
#Plot heatmap for these groups

#read CPM df with significantly changing genes in either cell line
RNAseq_CPM_sig <- read.delim("CPM_sig_set_IMRandBE.csv", header=TRUE, sep=",")

#scaled across both cell lines
scaled <- RNAseq_CPM_sig
for(r in 1:nrow(scaled)){
  scaled[r,2:(ncol(scaled)-1)] <- scale(as.numeric(scaled[r,2:(ncol(scaled)-1)]))
}
# heatmap with the pre-made groups and rows NOT clustered within those groups
#set colours + styles
col_fun= colorRamp2(c(-2, 0, 2), c('royalblue3', "ivory", 'darkorange'))
col_fun(seq(-3, 3))
ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
#heatmap
h1.1 <- Heatmap(as.matrix(scaled[,2:(ncol(scaled)-1)]),name= "z-score",cluster_rows=FALSE, cluster_columns=FALSE, # everything but second to last column bc its diff_type
                use_raster=FALSE, 
                row_dend_reorder=TRUE,
                show_row_dend=FALSE,
                show_column_dend=FALSE,
                show_column_names = FALSE,
                show_row_names=FALSE,
                column_split=c(rep('BE_induced',5),rep('BE_uninduced',5),rep('IMR_induced',5),rep('IMR_uninduced',5)),
                column_title = c("", "", "", ""),
                column_title_gp = gpar(fill = c("#FF95C4", "#E40E4F",'#797CD0', '#1A2C79'), col = "white", border = "black"), 
                col=col_fun,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 18),  # Font size of z-score label
                                            legend_direction = "horizontal",  # Make the legend horizontal
                                            title_position = "topleft"))   # Increase z-score label font size
h2.1 <- Heatmap(as.matrix(scaled[,'type']),cluster_columns=FALSE,cluster_rows =FALSE)
hlist2 <- h2.1+h1.1

#save
tiff(filename = "heatmap_IMRandBE_sig_groups_scaledtogether.tiff", 
     width = 7.5, height = 14, units = "in", res = 300, compression = "lzw")
draw(hlist2, row_split = scaled$type)
dev.off()
############

############
#######################################################################################################
################################ Creating Neuronal/Cell Division Gene Signatures ######################################
#######################################################################################################

#read enriched ontology terms from the 8 groups of genes 
#(significantly changing genes in at least one of the cell lines upon ASCL1 O/E)
goterms <- read.delim("compareCluster_groups_bothlines_sigset_BP.csv", header=TRUE, sep=",")
colnames(goterms)

####filter these gene ontologies to those relating to neuronal function/differentiation
# Define the words you want to include and exclude
include_neuronal <- c("neur", "synap", "axon")
exclude <- "negative"

# Filter rows that contain any of the include neuronal terms in the Description column
goterms_neuronal <- goterms[grepl(paste(include_neuronal, collapse = "|"), goterms$Description, ignore.case = TRUE), ]

# Exclude rows that contain "negative" in the Description column
goterms_neuronal <- goterms_neuronal[!grepl(exclude, goterms_neuronal$Description, ignore.case = TRUE), ]

#CHECK THROUGH THE LIST TO SEE IF IT MAKES SENSE

#save
write.csv(goterms_neuronal, "NeuronalGeneSignature_goterms.csv", row.names = FALSE, quote = FALSE)

#now save the gene names within those neuronal go terms.
#these will make up our Neuronal Gene Signature
data <- goterms_neuronal
# Separate geneID column
data_split <- data %>%
  separate_rows(geneID, sep = "/")
names(data_split)
# Group by geneID and summarize GO terms
result <- data_split %>%
  group_by(geneID) %>%
  summarize(ID = toString(unique(unlist(ID))),
            Description = toString(unique(unlist(Description))),
            GeneRatio = toString(unique(unlist(GeneRatio))),
            BgRatio = toString(unique(unlist(BgRatio))),
            pvalue = toString(unique(unlist(pvalue))),
            p.adjust = toString(unique(unlist(p.adjust))),
            qvalue = toString(unique(unlist(qvalue))),
            Count = toString(unique(unlist(Count))),
            type = toString(unique(unlist(type)))
  )

result <- result %>%
  rename(
    gene = geneID
  )
#save
write.csv(result, "NeuronalGeneSignature_genelist.csv", row.names = FALSE)

###now for Division gene signature

#filter these gene ontologies to those relating to cell division/mitosis
# Define the words you want to include and exclude
include_div <- c("mito", "cycl")
exclude <- "negative"

# Filter rows that contain any of the include neuronal terms in the Description column
goterms_div <- goterms[grepl(paste(include_div, collapse = "|"), goterms$Description, ignore.case = TRUE), ]

# Exclude rows that contain "negative" in the Description column
goterms_div <- goterms_div[!grepl(exclude, goterms_div$Description, ignore.case = TRUE), ]

#CHECK THROUGH THE LIST TO SEE IF IT MAKES SENSE
#in this case some terms are not cell cycle related. remove those manually

# Reset row names to sequential numbers, as currently not
rownames(goterms_div) <- NULL
# Assign new sequential row names
goterms_div <- goterms_div[seq_len(nrow(goterms_div)), ]

# Define which rows to exclude, based on go-terms (some mitochondrial terms came up)
rows_to_exclude <- c(14,17,40,42,43,45,50,51,53,54,59,60)
#exclude those rows
goterms_div <- goterms_div[-rows_to_exclude, ]
#save
write.csv(goterms_div, "CellDivisionGeneSignature_goterms.csv", row.names = FALSE, quote = FALSE)

#now save the gene names within those cell division go terms.
#these will make up our Cell Division Gene Signature
data1 <- goterms_div
# Separate geneID column
data_split1 <- data1 %>%
  separate_rows(geneID, sep = "/")
names(data_split1)
# Group by geneID and summarize GO terms
result1 <- data_split1 %>%
  group_by(geneID) %>%
  summarize(ID = toString(unique(unlist(ID))),
            Description = toString(unique(unlist(Description))),
            GeneRatio = toString(unique(unlist(GeneRatio))),
            BgRatio = toString(unique(unlist(BgRatio))),
            pvalue = toString(unique(unlist(pvalue))),
            p.adjust = toString(unique(unlist(p.adjust))),
            qvalue = toString(unique(unlist(qvalue))),
            Count = toString(unique(unlist(Count))),
            type = toString(unique(unlist(type)))
  )

result1 <- result1 %>%
  rename(
    gene = geneID
  )
#save
write.csv(result1, "CellDivisionGeneSignature_genelist.csv", row.names = FALSE)


#######################################################################################################
################# Plotting Neuronal/Cell Division Gene Signatures Average Expressions ######################################
#######################################################################################################

#To calculate average expressions of the gene signatures, need to also adjust for
#gene length, and calculate something resembling Transcript per Million (TPM)

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene lengths
gene_lengths <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"),
  mart = ensembl
)

# Calculate gene lengths
gene_lengths$length <- gene_lengths$end_position - gene_lengths$start_position
#change name of the GeneSymbol column to match the name in counts table
gene_lengths<-gene_lengths%>%
  rename(
    gene=external_gene_name
  )
#remove duplicates
gene_lengths_unique <- gene_lengths[!duplicated(gene_lengths$gene), ]

#read filtered counts file
counts <- read.delim("filtered_counts.csv", header=TRUE, sep=",")
#rename X to gene
counts<- counts%>%
  rename(
    gene=X
  )

#merge gene length with counts df
counts<- merge(gene_lengths_unique[,c("gene","length")], counts, by="gene")


#############Calculate ~TPM

#Divide values by gene lengths in kba (calculate RPK)
data_rpk <- counts
#
data_rpk[,3:ncol(data_rpk)] <- sapply(data_rpk[,3:ncol(data_rpk)], function(x) x/(data_rpk$length/1000))

#Scale per million
tpm <- data_rpk

tpm[,3:ncol(tpm)] <-apply(tpm[,3:ncol(tpm)],2,
                                        function(x) (x/sum(x,na.rm = TRUE))*1000000)

#remove lengths column now
tpm<- tpm[,-2]
#save tpm
write.csv(tpm, "TPM_IMRandBE.csv", row.names = FALSE, quote=FALSE)

################### Calculate and plot average expression of neuro and div signatures ######################
#read calculated TPM
tpm <- read.delim("TPM_IMRandBE.csv", header=TRUE, sep=",")

#read the list of neuronal or cell cycle genes
neuro <- read.delim("NeuronalGeneSignature_genelist.csv", header=TRUE, sep=",")
cycle <- read.delim("CellDivisionGeneSignature_genelist.csv", header=TRUE, sep=",")

# merge individually neuro/cycle with tpm
tpm_neuro <- merge(tpm, neuro_, by="gene")
tpm_cycle <- merge(tpm, cycle_, by="gene")

#make row names the gene name then remove the gene name column
rownames(tpm_neuro)<- tpm_neuro$gene
rownames(tpm_cycle)<- tpm_cycle$gene

tpm_neuro <- tpm_neuro[,-c(1)]
tpm_cycle <- tpm_cycle[,-c(1)]

#calculate mean tpm for neuro and cycle and create dataframe
mean_neuro <- sapply(tpm_neuro, mean)
mean_neuro <- as.data.frame(mean_neuro)
mean_neuro$Condition <- row.names(mean_neuro)

mean_cycle <- sapply(tpm_cycle, mean)
mean_cycle <- as.data.frame(mean_cycle)
mean_cycle$Condition <- row.names(mean_cycle)

####merge mean dataframes
means <- merge(mean_neuro, mean_cycle, by="Condition")

#plot
myplot<- ggplot(means, aes(x = mean_neuro, y = mean_cycle, fill = Condition)) +
  geom_point(shape = 21, size = 6) +
  
  # Uncomment and adjust geom_text_repel for labeling specific genes if needed
  # geom_text_repel(data = syn_prot %>% 
  #                 filter(log2FoldChange > 0.5) %>%
  #                 slice_min(padj, n = 5), 
  #                 aes(label = gene), 
  #                 size = 4, max.overlaps = 30) +
  
  labs(x = "Average Expression (TPM)\nNeuronal GO Terms", 
       y = "Average Expression (TPM)\nCell Division GO Terms") +
  
  scale_fill_manual(values = c("BE_uninduced_R1" = '#797CD0',
                               "BE_uninduced_R2" = '#797CD0',
                               "BE_uninduced_R3" = '#797CD0',
                               "BE_uninduced_R4" = '#797CD0',
                               "BE_uninduced_R5" = '#797CD0',
                               "BE_induced_R1" = "#1A2C79",
                               "BE_induced_R2" = "#1A2C79",
                               "BE_induced_R3" = "#1A2C79",
                               "BE_induced_R4" = "#1A2C79",
                               "BE_induced_R5" = "#1A2C79",
                               "IMR_uninduced_R1" = "#FF95C4",
                               "IMR_uninduced_R2" = "#FF95C4",
                               "IMR_uninduced_R3" = "#FF95C4",
                               "IMR_uninduced_R4" = "#FF95C4",
                               "IMR_uninduced_R5" = "#FF95C4",
                               "IMR_induced_R1" = "#E40E4F",
                               "IMR_induced_R2" = "#E40E4F",
                               "IMR_induced_R3" = "#E40E4F",
                               "IMR_induced_R4" = "#E40E4F",
                               "IMR_induced_R5" = "#E40E4F"), 
                    name = NULL, 
                    guide = "none") +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 16),  # Increase axis numbers size
    axis.title = element_text( size = 16),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )
# Save
ggsave("IMRvsBE_neuro_div_signatures_dotplot.tiff",plot=myplot, width = 5, height = 5, dpi = 300, compression = "lzw")


