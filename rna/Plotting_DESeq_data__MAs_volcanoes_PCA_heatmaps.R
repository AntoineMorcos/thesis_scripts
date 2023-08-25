# Clear environment
rm(list=ls())

# Load libraries
library(car)
library(multcomp)
library(reshape2)
library(tidyverse)
library(stats) 
library(dplyr)
library(tidyr)

# Libraries for this script
library(BiocManager) #install.packages("BiocManager")
library(DESeq2) #BiocManager::install("DESeq2")
library(ggrepel)
library(EnhancedVolcano) #BiocManager::install('EnhancedVolcano') #Nice volcano plots
library(pheatmap)
library("RColorBrewer")
library(ggrepel)

## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()

dev.off()






###########################################################################################################################
##                                                                                                                       ##
##                                                 Importing data                                                        ##
##                                                                                                                       ##
###########################################################################################################################


## This .Rdata is from the script "Metadata_and_counts_dfs_for_neutrophil_samples.R"

# Loading .Rdata saved
lnames <- load(file="data/Counts_and_TPM.RData")
lnames

# Loading .Rdata saved
lnames <- load(file="data/metadata_neutrophils_with_VIRGO.RData")
lnames


# This RData is from "Downstream_RNA-seq_DESeq2.R"
# Loading .Rdata saved 
lnames = load(file="data/DESeq2_essentials.RData")
lnames





###########################################################################################################################
##                                                                                                                       ##
##                                                         MA plots                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?plotMA

# Plot MA plot Ethnicity
png("plots/MA_plot_Ethnicity_same_axis.png")
par(mfrow=c(1,2))
plotMA(res_Ethnicity_sort, 
       ylim=c(-4,4), 
       main=paste0('White vs Black'))
plotMA(res_sLFC_Ethnicity_sort, 
       ylim=c(-4,4), 
       main=paste0('White vs Black shrinked'))
dev.off()







###########################################################################################################################
##                                                                                                                       ##
##                                             Volcano plots                                                             ##
##                                                                                                                       ##
###########################################################################################################################

# Tutorial for making nice Volcano plots https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# Info on function
?EnhancedVolcano
citation("EnhancedVolcano")



######## Plot - Over-ride colouring scheme with custom key-value pairs

# Create colour scheme for points of the DEGs by fold-change & p-values
Colour_scheme_volcano <- ifelse(
  res_sLFC_Ethnicity_sort$log2FoldChange < -1.5 & res_sLFC_Ethnicity_sort$padj < 0.05, 'royalblue',
  ifelse(res_sLFC_Ethnicity_sort$log2FoldChange > 1.5 & res_sLFC_Ethnicity_sort$padj < 0.05, 'red',
         'black'))
Colour_scheme_volcano[is.na(Colour_scheme_volcano)] <- 'black'
names(Colour_scheme_volcano)[Colour_scheme_volcano == 'red'] <- 'Upregulated'
names(Colour_scheme_volcano)[Colour_scheme_volcano == 'black'] <- 'No significant difference'
names(Colour_scheme_volcano)[Colour_scheme_volcano == 'royalblue'] <- 'Downregulated'


# Making a vector of the top 10 upregulated DEGs & top 10 downregulated DEGs
top_10_up_and_down <- c(top_10_genes_sLFC_up$Gene, top_10_genes_sLFC_down$Gene)
top_10_up_and_down 


# Plotting the volcano... only top 10 upregulated & downregulated labeled
?EnhancedVolcano
pdf("plots/Volcano_plot_top_10_up_and_down_DEGs_labeled.pdf", width=11, height=10)
EnhancedVolcano(res_sLFC_Ethnicity_sort, # final res from DESeq
                lab = rownames(res_sLFC_Ethnicity_sort),
                x = 'log2FoldChange',  # data to plot within res on the x-axis
                y = 'padj', # data to plot within res on the y-axis
                selectLab = c(top_10_up_and_down), #vector of the top 10 upregulated DEGs & top 10 downregulated DEGs
                labSize = 5, # font size of labels
                xlab = bquote(~Log[2]~ 'fold change'), # label x axis
                ylab = bquote("-" ~Log[10]~ '(FDR adjusted p-value)'), # label y axis
                title = NULL, subtitle=NULL,
                pCutoff = 0.05, FCcutoff = 1.5, #lines
                pointSize = 2,
                colCustom = Colour_scheme_volcano, # the colour scheme we made above with our ifelse statements
                colAlpha = 0.6, # Alpha - transparency of points
                legendPosition = 'none', # no legend - makes it messy. Previously was "left"
                drawConnectors = TRUE, widthConnectors = 0.4, colConnectors = 'grey', arrowheads = F, # Connector line between gene text and point options
                #lengthConnectors = 5,
                gridlines.major = T, gridlines.minor = F,
                border = 'partial', borderWidth = 1, borderColour = 'black') # Border options 
# no warning! :D All text labels plotted  
dev.off()
# this plotting function only allows labelling a small amount of DEGs otherwise it doesn't plot & says error






###########################################################################################################################
##                                                                                                                       ##
##                                           Normalising gene expression                                                 ##
##                                                                                                                       ##
###########################################################################################################################

# In normalisation options  you can either do vst or rlog
# Variance-stabilizing-transformation (vst) 
# Regularized-log-transformation (rld)


# For a fully unsupervised transformation, set blind = TRUE (which is the default)

# Normalisation needed for PCA etc
rld <- rlog(dds1, blind = TRUE) #?rlog
vst <- vst(dds1, blind = TRUE)

# I went with the rlog normalisation as final choice, but I did previously run all plots on both rlog & vst



# Saving rld & vst
save(rld, vst, file="data/rld_and_vst.RData")

## Loading .Rdata from the script "Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.R"
# lnames <- load(file="data/rld_and_vst.RData")
# lnames






###########################################################################################################################
##                                                                                                                       ##
##                                     Heatmaps of the correlations between samples                                      ##
##                                                                                                                       ##
###########################################################################################################################


############################################################### heatmap prep
# Info
?pheatmap
citation("pheatmap")

# metadata that we want in the heatmap
metadata_for_heatmap <- metadata_CVF%>%
  dplyr::select(BMI_category, Ethnicity, Alicia_group, clusters_hc_k4) %>% 
  dplyr::rename(Cluster=clusters_hc_k4, 
                Microbiome=Alicia_group) 


# Changing all spaces & dashes in the df to underscores https://stackoverflow.com/questions/20760547/removing-whitespace-from-a-whole-data-frame-in-r
metadata_for_heatmap <- as.data.frame(apply(metadata_for_heatmap,2,function(x)gsub('\\s+', '_',x)))
metadata_for_heatmap <- as.data.frame(apply(metadata_for_heatmap,2,function(x)gsub('-', '',x)))
metadata_for_heatmap <- as.data.frame(apply(metadata_for_heatmap,2,function(x)gsub('/', '_and_',x)))
metadata_for_heatmap <- as.data.frame(apply(metadata_for_heatmap,2,function(x)gsub('\\._', '_',x))) #https://stackoverflow.com/questions/20309876/r-how-to-replace-in-a-string
metadata_for_heatmap 

# Choosing colours for heatmap https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
metadata_for_heatmap
my_colours <- list(
  Cluster = c("1"="#FF226F", "2"="#3D6A89", "3"="#BEABFF", "4"="#30B2FF"),
  Microbiome=c(L_crispatus="#6ED088", L_iners="#99CCFF", L_iners_and_L_jensenii="blue", L_crispatus_and_L_iners_and_G_vaginalis="black", G_vaginalis_and_L_gasseri="red", G_vaginalis_and_A_vaginae="#FCA64C"),
  Ethnicity = c(White = "#F68F44", Black = "purple"),
  BMI_category = c(Healthy_weight="#00CC66", Overweight="yellow", Obese="#FF9933", Morbidly_obese="red"))
my_colours
###############################################################




############################################################### rld heatmap
#extract the rld matrix from the object
rld_mat_dds1 <- assay(rld)
head(rld_mat_dds1)

# Changing from a matrix to a df so we can change the rownames
rld_mat_dds1 <- as.data.frame(rld_mat_dds1)
head(rld_mat_dds1)


#Compute pairwise correlation values
rld_corr_dds1 <- stats::cor(rld_mat_dds1) # method="pearson" is the default (I ran a test)
head(rld_corr_dds1)
citation("stats")


# checks before heatmap
dim(metadata_for_heatmap) 
dim(rld_corr_dds1)
rownames(metadata_for_heatmap) 
rownames(rld_corr_dds1)


# Plot the heatmap FOR THESIS
pdf("plots/rlog_correlation.pdf", width=10, height=6)
par(mfrow=c(1,1))
pheatmap(rld_corr_dds1, annotation_col=metadata_for_heatmap, border_color=NA, 
         show_rownames=T, show_colnames=T, #Show IDs 
         annotation_colors = my_colours,
         legend=T, legend_breaks=c(seq(0,1,0.05)), 
         display_numbers=T, fontsize_number=6 #Prints the distance numbers over the cells
)
dev.off()
###############################################################






###########################################################################################################################
##                                                                                                                       ##
##                                            PCA plots of rlog/rld                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Plot PCA plot for paper- Ethnicity & Cluster clusters_hc_k4
PCA_for_paper_data  <- DESeq2::plotPCA(rld, intgroup=c("Ethnicity", "clusters_hc_k4"), returnData = TRUE)
PCA_for_paper_percentVar <- round(100 * attr(PCA_for_paper_data, "percentVar"))
PCA_for_paper <- ggplot(PCA_for_paper_data, aes(x=PC1, y=PC2, color=clusters_hc_k4, shape=Ethnicity)) +
  geom_point(size=4, stroke=1.5) +
  scale_color_manual(values = c("#FFA1CD", "turquoise", "#BEABFF", "#30B2FF")) +
  scale_shape_manual(values=c(3, 16), na.translate=F) +
  xlab(paste0("PC1: ", PCA_for_paper_percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", PCA_for_paper_percentVar[2], "% variance")) +
  guides(color=guide_legend(order=1), shape = guide_legend(order=2)) +
  labs(color="Cluster")
PCA_for_paper
ggsave("plots/PCA/PCA_for_paper.png", PCA_for_paper, height=4, width=5, dpi=700)

# Plot PCA plot for paper- Ethnicity & Cluster clusters_hc_k4 - diff colours
PCA_for_thesis <- ggplot(PCA_for_paper_data, aes(x=PC1, y=PC2, color=clusters_hc_k4, shape=Ethnicity)) +
  geom_point(size=4, stroke=1.5) +
  scale_color_manual(values = c("#FF226F", "#3D6A89", "#BEABFF", "#30B2FF")) +
  scale_shape_manual(values=c(3, 16), na.translate=F) +
  xlab(paste0("PC1: ", PCA_for_paper_percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", PCA_for_paper_percentVar[2], "% variance")) +
  guides(color=guide_legend(order=1), shape = guide_legend(order=2)) +
  labs(color="Cluster")
ggsave("plots/PCA/PCA_for_thesis.pdf", PCA_for_thesis, height=4, width=5)








###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting all data                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data
save.image(file = "data/ALL_objects_in_Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.RData", compress=F)


# # Loading .Rdata saved
# lnames = load(file="data/ALL_objects_in_Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.RData")
# lnames


