# Load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db) 

#######################
#load genes proximal to IMR enriched ASCL1 peaks (up to 50kb)
IMR_enriched_genes <- read.delim("IMR_up_peakAnno_lessthan50kb_unique.csv", header=TRUE, sep=",")
colnames(IMR_enriched_genes) <- "gene"
#load RNAseq CPM data
CPM <- read.delim("RNAseq_CPM.csv", header=TRUE, sep=",")
#load significantly changing genes
sig <- read.delim("sig_set_IMRandBE.csv", header=TRUE, sep=",")

##shortlist CPM file to only sig genes
CPM_sig <- merge(CPM, sig[,c("gene", "type")], by="gene")
# shortlist to UP genes only
CPM_sig_up <- CPM_sig[grepl("up", CPM_sig$type, ignore.case = TRUE), ]
#shortlist to genes proximal to IMR enriched ASCL1 peaks
CPM_sig_up_IMRenr <- merge(CPM_sig_up, IMR_enriched_genes, by="gene")

#count how many in each group
table <- CPM_sig_up_IMRenr %>%
  count(type) 
# combine IMR_up with IMR_up/BE_down, and BE_up with IMR_down/BE_up
CPM_sig_up_IMRenr_ <- CPM_sig_up_IMRenr %>%
  mutate(type=case_when(
    type=="IMR up" ~ "IMR up",
    type=="BE up" ~ "BE up",
    type=="Both up" ~ "Both up",
    type=="IMR up/BE down" ~ "IMR up",
    type=="IMR down/BE up" ~ "BE up",
    TRUE ~ "ns"
  ))
#count how many in each group
table_updated <- CPM_sig_up_IMRenr_ %>%
  count(type) 
###

############################### PLOT HEATMAP #################
#scale across both cell lines
scaled <- CPM_sig_up_IMRenr_
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
h2.1 <- Heatmap(as.matrix(scaled[,'type']),cluster_columns=FALSE,cluster_rows =FALSE, show_row_names = FALSE,show_column_names = FALSE, show_heatmap_legend = FALSE)

hlist2 <- h2.1+h1.1

#save
tiff(filename = "heatmap_IMRenrpeaks_RNAseqUP_.tiff", 
     width = 7.5, height = 14, units = "in", res = 300, compression = "lzw")
draw(hlist2, row_split = scaled$type)
dev.off()
#save smaller
tiff(filename = "heatmap_IMRenrpeaks_RNAseqUP_3.tiff", 
     width = 5, height = 5, units = "in", res = 300, compression = "lzw")
draw(hlist2, row_split = scaled$type)
dev.off()


####################### GENE ONTOLOGY ######################
#run compareCluster
go_groups<-compareCluster(gene~type, data=CPM_sig_up_IMRenr_[,c('gene','type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='BP',pAdjustMethod='fdr', pvalueCutoff=0.01,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
#plot result
plot_groups <- dotplot(go_groups, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange')) +
  theme(
    legend.text = element_text(size = 16),       # Increase the size of the legend text
    legend.title = element_text(size = 18),      # Increase the size of the legend title
    legend.key.height = unit(1.5, "cm"),         # Increase the height of the legend keys
    legend.key.width = unit(1, "cm")             # Adjust the width of the legend keys if needed
  )
#save
ggsave(filename = "compareCluster_IMRenrpeaks_RNAseqUP_BP.png", plot = plot_groups, width = 8, height = 10, units = "in", dpi = 300)
ggsave(filename = "compareCluster_IMRenrpeaks_RNAseqUP_BP.pdf", plot = plot_groups, width = 8, height = 10, units = "in", dpi = 300)
write.csv(go_groups,("compareCluster_IMRenrpeaks_RNAseqUP_BP.csv"))
#save smaller
ggsave(filename = "compareCluster_IMRenrpeaks_RNAseqUP_BP_1.png", plot = plot_groups, width = 6, height = 8, units = "in", dpi = 300)

# Create the plot with horizontal orientation
plot_groups_h <- dotplot(go_groups, font.size = 14) +
  scale_fill_gradientn(colors = c('royalblue3', "ivory", 'darkorange')) +
  coord_flip() +  # Flip coordinates to make the plot horizontal
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
    legend.text = element_text(size = 16),       # Increase the size of the legend text
    legend.title = element_text(size = 18),      # Increase the size of the legend title
    legend.key.height = unit(1.5, "cm"),         # Increase the height of the legend keys
    legend.key.width = unit(1, "cm")             # Adjust the width of the legend keys if needed
  )
plot_groups_h
# Save the plot
ggsave(filename = "compareCluster_IMRenrpeaks_RNAseqUP_BP_horizontal_1.png", 
       plot = plot_groups_h, 
       width = 8,  # Adjust width to accommodate horizontal layout
       height = 5,  # Adjust height to accommodate horizontal layout
       units = "in", 
       dpi = 300)

#########################################################################
