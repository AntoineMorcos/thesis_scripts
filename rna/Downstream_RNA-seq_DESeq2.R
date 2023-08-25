# Clear environment
rm(list=ls())

# Load libraries
library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(ggplot2)
library(dplyr)
library(tidyr)

# Libraries for this script
library(BiocManager) #install.packages("BiocManager")
library(DESeq2) #BiocManager::install("DESeq2")
library(devtools)
library(ggbiplot) #install_github("vqv/ggbiplot")
library(ggrepel)
library(sva) #BiocManager::install("sva")



#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()



# This R script uses the data outputted from featureCounts






###########################################################################################################################
##                                                                                                                       ##
##                                          Load in counts and metadata                                                  ##
##                                                                                                                       ##
###########################################################################################################################

## This .Rdata is from the script "Metadata_and_counts_dfs_for_neutrophil_samples.R"

# Loading .Rdata saved
lnames <- load(file="data/Counts_and_TPM.RData")
lnames

# Loading .Rdata saved
lnames <- load(file="data/metadata_neutrophils_with_VIRGO.RData")
lnames






###########################################################################################################################
##                                                                                                                       ##
##                              Begin DESeq2 - make dds from count data and metadata                                     ##
##                                                                                                                       ##
###########################################################################################################################

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start

# DESeq of RNA-seq CVF samples (Excluding blood samples)

### Check that sample names match in both files
all(colnames(counts_CVF_adjusted) %in% rownames(metadata_CVF))
all(colnames(counts_CVF_adjusted) == rownames(metadata_CVF))
# they match

# Checks
head(counts_CVF_adjusted)
head(metadata_CVF)


# Making a DESeq2 dataset object
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=counts_CVF_adjusted, colData=metadata_CVF,
                              design= ~ Ethnicity)





###########################################################################################################################
##                                                                                                                       ##
##                                            Look at initial dds object                                                 ##
##                               & filtering out genes w/ 0 counts in all samples                                        ##
##                                                                                                                       ##
###########################################################################################################################


# Look
dds
# class: DESeqDataSet 
# dim: 38557 10 
# metadata(1): version
# assays(1): counts
# rownames(38557): DDX11L1 WASH7P ... LOC105377244 LOC105369226
# rowData names(0):
#   colnames(10): ID_4 ID_5 ... ID_3 ID_10
# colData names(64): IDs IDs_2 ... clusters_km_k5 clusters_km_k6

# metadata column information
colData(dds)

# head of counts
head(assay(dds))
head(counts_CVF_adjusted) #same as above

# summary of counts per sample
head(summary(assay(dds)))
summary(counts_CVF_adjusted) #same as above

# number of reads per ID
colSums(assay(dds))




# Removing genes where all samples have 0 counts
dds_no0s <- dds[rowSums(counts(dds)) > 0, ]



# Checking 
nrow(dds) #38557
nrow(dds_no0s) #15982

colSums(assay(dds)) 
# After ComBat:
#     ID_4     ID_5     ID_9     ID_6     ID_7     ID_1     ID_2     ID_8     ID_3    ID_10 
# 18310664 16773793 15838403 15573435 15401082 12550658 15675338 16808650 14153179 18950090

colSums(assay(dds_no0s)) #sums are the same so no rows with genes were dropped 
#     ID_4     ID_5     ID_9     ID_6     ID_7     ID_1     ID_2     ID_8     ID_3    ID_10 
# 18310664 16773793 15838403 15573435 15401082 12550658 15675338 16808650 14153179 18950090







###########################################################################################################################
##                                                                                                                       ##
##                                       DESeq() algorithm on dds_no0s to make dds1                                      ##
##                                                                                                                       ##
###########################################################################################################################

# info
?DESeq


# Perform the DESeq2 differential gene expression command
set.seed(1)
dds1 <- DESeq(dds_no0s)





###########################################################################################################################
##                                                                                                                       ##
##                                                   Getting top DEGs                                                    ##
##                                                                                                                       ##
###########################################################################################################################

resultsNames(dds1) #resultsNames returns the names of the estimated effects (coefficents) of the model
#[1] "Intercept"               "Ethnicity_Black_vs_White"



# Make results table - base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values
# Results comparing specific conditions - Ethnicity
?DESeq2::results
res_Ethnicity <- results(dds1, contrast=c("Ethnicity", "White", "Black"), alpha=0.05, lfcThreshold = 0.58) 
#contrast is used to do White vs Black, rather than Black vs White as trying to change the reference factor the way the website suggests does not work https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#note-on-factor-levels

# Interpreting lfcThreshold:
# log2(0.58),     2^0.58=~1.5,    so the genes have a fold change of >|1.5|  
# https://support.bioconductor.org/p/62927/       https://support.bioconductor.org/p/9140109/



mcols(res_Ethnicity, use.names=T)
# DataFrame with 6 rows and 2 columns
# type            description
# <character>            <character>
# baseMean            intermediate mean of normalized c..
# log2FoldChange      results log2 fold change (ML..
# lfcSE               results standard error: Ethn..
# stat                results Wald statistic: Ethn..
# pvalue              results Wald test p-value: E..
# padj                results   BH adjusted p-values

summary(res_Ethnicity)
# out of 15176 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0.58 (up)    : 127, 0.84%
# LFC < -0.58 (down) : 73, 0.48%
# outliers [1]       : 1156, 7.6%
# low counts [2]     : 5908, 39%
# (mean count < 4)



# Order results by p values
res_Ethnicity_sort <- res_Ethnicity[order(res_Ethnicity$pvalue),]
head(res_Ethnicity_sort)


#List of the top 20 differentially expressed genes by most sig p-value
res_Ethnicity_sort@rownames[1:20] 



# How many adjusted p-values were less than 0.05?
sum(res_Ethnicity_sort$padj < 0.05, na.rm=T) 
# Now 200 with only ~Ethnicity and a min lfcThreshold = 0.58, i.e. actual fold change of 1.5






###########################################################################################################################
##                                                                                                                       ##
##                                          LFC estimates - Ethnicity                                                    ##
##                                                                                                                       ##
###########################################################################################################################

# library(apeglm) #BiocManager::install("apeglm")
library(ashr) #BiocManager::install("ashr")
?lfcShrink
?ash #"mplements Empirical Bayes shrinkage and false discovery rate methods based on unimodal prior distributions."
citation("ashr")

resultsNames(dds1)

# Log fold change shrinkage for visualization and ranking 
?DESeq2::lfcShrink
res_sLFC_Ethnicity <- lfcShrink(dds1, res=res_Ethnicity, contrast=c("Ethnicity", "White", "Black"), type="ashr")
# using 'ashr' for LFC shrinkage. If used in published research, please cite: Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. https://doi.org/10.1093/biostatistics/kxw041
# info on lfc shrink https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink


## When I tried using apeglm rather than ashr shrinkage
# Error in lfcShrink(dds1, res = res_Ethnicity, contrast = c("Ethnicity",  : type='apeglm' shrinkage only for use with 'coef'

# We can now play with the results table and for example order it by the smallest p value:
res_sLFC_Ethnicity_sort <- res_sLFC_Ethnicity[order(res_sLFC_Ethnicity$pvalue),]

# Looking
head(res_sLFC_Ethnicity_sort) #LFC
#log2 fold change (MMSE): Ethnicity White vs Black 
# Wald test p-value: Ethnicity White vs Black 
# DataFrame with 6 rows and 5 columns
#               baseMean log2FoldChange     lfcSE      pvalue        padj
#              <numeric>      <numeric> <numeric>   <numeric>   <numeric>
# SPRR2E       6713.3328       -15.9505   1.36948 8.01802e-31 7.15047e-27
# UCKL1-AS1     647.5065       -12.1325   1.40288 1.84526e-17 7.06467e-14
# CD300LD       680.7159       -12.9662   1.51744 2.37654e-17 7.06467e-14
# LOC105372803  326.7821       -11.4491   1.53128 1.70841e-13 3.80891e-10
# KLK10         204.3586        -9.5961   1.29628 9.38485e-13 1.67388e-09
# ELOVL7         80.7799        21.5165   3.19902 1.35293e-12 2.01091e-09


head(res_Ethnicity_sort) # not sLFC
# log2 fold change (MLE): Ethnicity White vs Black 
# Wald test p-value: Ethnicity White vs Black 
# DataFrame with 6 rows and 6 columns
#            baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#           <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# SPRR2E       6713.3328      -16.29889   1.36178 -11.54290 8.01802e-31 7.15047e-27
# UCKL1-AS1     647.5065      -12.52724   1.40504  -8.50314 1.84526e-17 7.06467e-14
# CD300LD       680.7159      -13.42511   1.51587  -8.47373 2.37654e-17 7.06467e-14
# LOC105372803  326.7821      -11.91635   1.53821  -7.36983 1.70841e-13 3.80891e-10
# KLK10         204.3586       -9.91343   1.30734  -7.13924 9.38485e-13 1.67388e-09
# ELOVL7         80.7799       23.31417   3.20706   7.08879 1.35293e-12 2.01091e-09




###########################################################################################################################
##                                                                                                                       ##
##                                         Export data for IPA - Ethnicity                                               ##
##                                                                                                                       ##
###########################################################################################################################


### Getting file ready for export
# Converting to a df
head(res_sLFC_Ethnicity_sort)
res_sLFC_Ethnicity_sort_df <- as.data.frame(res_sLFC_Ethnicity_sort)
head(res_sLFC_Ethnicity_sort_df)

# Adding a column of the Genes
res_sLFC_Ethnicity_sort_df$Gene <- rownames(res_sLFC_Ethnicity_sort_df)
str(res_sLFC_Ethnicity_sort_df)
head(res_sLFC_Ethnicity_sort_df)

# Rearranging column order
res_sLFC_Ethnicity_sort_df <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, everything())

# Rearranging column order, dropping unnecessary columns and keeping only sig adjusted p-values
dim(res_sLFC_Ethnicity_sort_df) #[1] 15982     6
res_sLFC_Ethnicity_sort_df_filtered <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, log2FoldChange, padj) %>% 
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>0.58)
head(res_sLFC_Ethnicity_sort_df_filtered)
dim(res_sLFC_Ethnicity_sort_df_filtered) #[1] 200   3 (was 404 3 before removing Outcome as a covariate & fixing the mislabelled ID)
summary(res_sLFC_Ethnicity_sort_df_filtered) #checks  on filtering


# Writing out the data as a txt without quotes
write.table(res_sLFC_Ethnicity_sort_df_filtered, "data/for_IPA/White_vs_Black_sLFC_ashr_1.5FC_0.05padj.txt", 
            sep="\t", row.names=F, col.names=T, quote=F)




###########################################################################################################################
##                                                                                                                       ##
##                 Checks - See if sLFC has made any sLFC<|0.58| , or if they are all happily >|0.58|                    ##
##                                                                                                                       ##
###########################################################################################################################

# Compare w/ summary 
summary(res_Ethnicity_sort)
summary(res_sLFC_Ethnicity_sort)
127+73 #200 DEGs

# This df in the previous section was filtered by filter(padj<0.05 & abs(log2FoldChange)>0.58)
dim(res_sLFC_Ethnicity_sort_df_filtered) #[1] 200   3 
# So still 200 DEGS with sLFC>|0.58| 

### These checking_up & checking_down are a variant of res_sLFC_Ethnicity_sort_df_filtered, where we look at up & down separately, not absolute value of 0.58
# Keeping only sig adjusted p-values & sLFC>0.58
dim(res_sLFC_Ethnicity_sort_df) #[1] 15982     6
checking_up <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, log2FoldChange, padj) %>% 
  dplyr::filter(padj<0.05 & log2FoldChange>0.58)

# Keeping only sig adjusted p-values & sLFC<-0.58
checking_down <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, log2FoldChange, padj) %>% 
  dplyr::filter(padj<0.05 & log2FoldChange<(-0.58))

#Checks
dim(checking_up) #[1] 127   3
summary(checking_up) # Min.   : 4.094      #so it's fine
dim(checking_down) #[1] 73  3
summary(checking_down) # Max.   :-15.951      #so it's fine

#200 DEGs all with sLFCs happily >|0.58| 




###########################################################################################################################
##                                                                                                                       ##
##                                                Top genes by sLFC                                                      ##
##                                                                                                                       ##
###########################################################################################################################
#?abs
#?order

# Sorted by absolute LFC
sort_sLFC_abs <- res_sLFC_Ethnicity_sort_df_filtered[order(abs(res_sLFC_Ethnicity_sort_df_filtered$log2FoldChange), decreasing=T),]
top_10_genes_sLFC_abs <- head(sort_sLFC_abs, 10)
top_10_genes_sLFC_abs

# Sorted by LFC upregulated 
sort_sLFC_up <- res_sLFC_Ethnicity_sort_df_filtered[order(res_sLFC_Ethnicity_sort_df_filtered$log2FoldChange, decreasing=T),]
top_10_genes_sLFC_up <- head(sort_sLFC_up, 10)
top_10_genes_sLFC_up

# Sorted by LFC downregulated 
sort_sLFC_down <- res_sLFC_Ethnicity_sort_df_filtered[order(res_sLFC_Ethnicity_sort_df_filtered$log2FoldChange, decreasing=F),]
top_10_genes_sLFC_down <- head(sort_sLFC_down, 10)
top_10_genes_sLFC_down


#Export to CSV
write_csv(top_10_genes_sLFC_up, "data/top_10_genes_sLFC_up.csv")
write_csv(top_10_genes_sLFC_down, "data/top_10_genes_sLFC_down.csv")




###########################################################################################################################
##                                                                                                                       ##
##       Export RData for comparing how many genes are in common when including or excluding the PTB ID                  ##
##                                                                                                                       ##
###########################################################################################################################

# Upregulated res df
dim(res_sLFC_Ethnicity_sort_df) #[1] 15982     6
up_ComBat <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, log2FoldChange, padj) %>% 
  dplyr::filter(padj<0.05 & log2FoldChange>0.58) %>% 
  dplyr::arrange(desc(log2FoldChange))

# Look
head(up_ComBat)
#                      Gene log2FoldChange         padj
# ELOVL7             ELOVL7       21.51647 2.010911e-09
# PEX19               PEX19       20.88268 4.323452e-08
# ZNF630             ZNF630       20.19653 1.174863e-07
# LOC101929356 LOC101929356       20.12108 1.236047e-07
# LINC02470       LINC02470       19.87564 1.815100e-07
# IL5RA               IL5RA       19.72674 2.223415e-07
dim(up_ComBat) # 127   3
summary(up_ComBat) #checks  on filtering

# Downregulated res df
dim(res_sLFC_Ethnicity_sort_df) #[1] 15982     6
down_ComBat <- res_sLFC_Ethnicity_sort_df %>% 
  dplyr::select(Gene, log2FoldChange, padj) %>% 
  dplyr::filter(padj<0.05 & log2FoldChange<0.58) %>% 
  dplyr::arrange(log2FoldChange)

# Look
head(down_ComBat)
#              Gene log2FoldChange         padj
# SPRR2E     SPRR2E      -15.95053 7.150468e-27
# SERPINB4 SERPINB4      -14.22344 4.945725e-07
# PGM3         PGM3      -14.16064 9.417740e-06
# ATL1         ATL1      -13.75252 4.008996e-05
# BCO2         BCO2      -13.60854 4.452679e-05
# SPATS2L   SPATS2L      -13.60097 5.719570e-05
dim(down_ComBat) # 73  3
summary(down_ComBat) #checks  on filtering 

#Renaming these for easier comprehension in the following script PTB__comparing_include_or_exclude_PTB.R
top_up_ComBat <- top_10_genes_sLFC_up
top_down_ComBat <- top_10_genes_sLFC_down

###Saving the data from this script
save(up_ComBat, down_ComBat, top_up_ComBat, top_down_ComBat, file="data/dds_set_up_comparison/ComBat.RData")


# Loading .Rdata saved 
# lnames = load(file="data/dds_set_up_comparison/ComBat.RData")
# lnames




###########################################################################################################################
##                                                                                                                       ##
##                                         Export DESeq2 data as .RData for other scripts                                ##
##                                                                                                                       ##
###########################################################################################################################

###Saving the data from this script
save(
     #dds
     dds1, 
     
     #differential expression white vs black results
     res_Ethnicity_sort, 
     res_sLFC_Ethnicity_sort, 
     res_sLFC_Ethnicity_sort_df_filtered,
     
     #top differentially expressed genes - for plotting single gene counts
     top_10_genes_sLFC_up, 
     top_10_genes_sLFC_down,
     
     file="data/DESeq2_essentials.RData")


## Loading .Rdata saved 
# lnames = load(file="data/DESeq2_essentials.RData")
# lnames

