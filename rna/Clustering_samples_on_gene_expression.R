# Clear environment
rm(list=ls())

# Load libraries
library(reshape2) # heatmaps
library(BiocManager) #install.packages("BiocManager")
library(tidyverse) # ggplot2, dplyr, tidyr etc
library(ggrepel) # labels on ggplots
library(stats) # for Wilcoxon or other stats tests
library(cluster) # has daisy function & K-means functions
library(DESeq2) #for assay()
library(pheatmap) #heatmaps
library(vegan) #Bray-Curtis
library(dendextend) #Dendograms coloured
library(purrr) # exploring k values in k-means clustering
 

## Set wd
# setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()





###########################################################################################################################
##                                                                                                                       ##
##                                                 Importing data                                                        ##
##                                                                                                                       ##
###########################################################################################################################


## Loading .Rdata from the script "Metadata_and_counts_dfs_for_neutrophil_samples.R"
lnames <- load(file="data/metadata_neutrophils_with_VIRGO.RData")
lnames


# Loading RData from "Downstream_RNA-seq_DESeq2.R"
lnames = load(file="data/DESeq2_essentials.RData")
lnames


# Loading RData from "Metagenome_VIRGO_basics.R"
lnames = load(file="data/metagenome/metagenome_abundances_VIRGO.RData")
lnames

# Loading .Rdata from the script "Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.R"
lnames <- load(file="data/rld_and_vst.RData")
lnames








###########################################################################################################################
##                                                                                                                       ##
##                                        K-means clustering on gene expression                                          ##
##                                                                                                                       ##
###########################################################################################################################


################################################################################## looking at df 
# Gene expression (rld)
rld_t[,1:5]
str(rld_t) 
# 'data.frame':	10 obs. of  15983 variables
# all num
##################################################################################




################################################################################## Elbow plot 
# Run models with varying value of k (centers)
?purrr::map_dbl
set.seed(1)
tot_withinss <- map_dbl(1:9,  function(k){
  model <- kmeans(x = rld_t, centers = k)
  model$tot.withinss #total within-cluster sum of square 
})

# Generate df containing both k and tot_withinss 
elbow_df <- data.frame(
  k = 1:9,
  tot_withinss = tot_withinss)

# Look
elbow_df

# Plot the elbow plot 
elbow_kmeans_plot <- ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() +
  scale_x_continuous(breaks = 1:9, minor_breaks = 1:9) +
  scale_y_continuous(labels = scales::comma, limits=c(0,350000)) + #https://stackoverflow.com/questions/14563989/force-r-to-stop-plotting-abbreviated-axis-labels-e-g-1e00-in-ggplot2
  labs(title="Elbow plot of k-means clustering of rlog gene expression", y="Total within-cluster sum of squares")
ggsave("plots/Clustering_samples/k_means/elbow_kmeans_plot.png", elbow_kmeans_plot, dpi=600, width=6, height=4)
# Maybe k=4? As the curve flattens a bit there... but there's no obvious elbow  
##################################################################################




################################################################################## Silhouette analysis 
# Info on k-means function
?pam

# Use map_dbl (in purrr)  to run many models with varying value of k
set.seed(1)
sil_width <- map_dbl(2:9,  function(k){
  model <- pam(x = rld_t, k = k)
  model$silinfo$avg.width
})

# Generate a data frame containing both k and sil_width
set.seed(1)
sil_df <- data.frame(
  k = 2:9,
  sil_width = sil_width
)
sil_df
#   k  sil_width
# 1 2 0.13199932
# 2 3 0.13551568
# 3 4 0.14021148
# 4 5 0.07287019
# 5 6 0.05954216
# 6 7 0.04294660
# 7 8 0.01304444
# 8 9 0.02176970


# Plot the relationship between k and sil_width to find the k with the highest average silhouette width
silhouette_averages_kmeans_plot <- ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() +
  scale_x_continuous(breaks = 2:9, minor_breaks = 2:9) + 
  scale_y_continuous(breaks=seq(-1, 1, 0.04), limits=c(0,0.16)) +
  labs(y="Average silhouette width")
ggsave("plots/Clustering_samples/k_means/silhouette_kmeans.pdf", silhouette_averages_kmeans_plot, width=3.5, height=2.5)
# k=4 has the highest average silhouette width
##################################################################################




################################################################################## k=2-9
# Get cluster assignments using k-means clustering for k=2 to k=9
set.seed(100) #(I originally did not set the seed but 100 comes close to the original run)
clusters_kmeans <- map_dfr(2:9,  function(k){
  model_kmeans <- kmeans(rld_t, centers = k)
  model_kmeans$cluster
})

# Generate a data frame containing both k and the cluster assignments
cluster_kmeans_df <- data.frame(k = 2:9, clusters_kmeans)
cluster_kmeans_df
#   k ID_4 ID_5 ID_9 ID_6 ID_7 ID_1 ID_2 ID_8 ID_3 ID_10
# 1 2    1    1    1    1    1    1    2    1    1     1
# 2 3    3    3    3    3    3    1    2    3    3     3
# 3 4    1    1    1    1    1    4    2    1    3     1
# 4 5    4    4    4    2    4    1    3    4    5     4
# 5 6    5    5    5    2    5    4    1    6    3     5
# 6 7    5    5    1    4    5    2    6    3    7     5
# 7 8    8    8    6    1    5    3    2    7    4     8
# 8 9    9    2    7    1    4    3    5    6    8     9
##################################################################################








###########################################################################################################################
##                                                                                                                       ##
##                                            Exporting cluster assignments                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Look
clusters_hc_cl 
cluster_kmeans_df
str(clusters_hc_cl)
str(cluster_kmeans_df)
rownames(clusters_hc_cl)
rownames(cluster_kmeans_df)


################################### Getting dfs in same format
# Hierachial clustering transpose
clusters_hc_cl_t <- as.data.frame(t(clusters_hc_cl))
clusters_hc_cl_t
clusters_hc_cl_t$IDs_2 <- rownames(clusters_hc_cl_t)

# K-means manipulation & transpose 
rownames(cluster_kmeans_df) <- paste("clusters_km_k", cluster_kmeans_df$k, sep="")
cluster_kmeans_df
cluster_kmeans_df$k <- NULL
cluster_kmeans_df_t <- as.data.frame(t(cluster_kmeans_df[2:5,])) # keep only k=3-6 in transposed df
cluster_kmeans_df_t
cluster_kmeans_df_t$IDs_2 <- rownames(cluster_kmeans_df_t)
###################################


################################### merge
# Look
clusters_hc_cl_t
cluster_kmeans_df_t

# Merge
clustering_IDs <- merge(clusters_hc_cl_t, cluster_kmeans_df_t, by="IDs_2")
clustering_IDs
#    IDs_2 clusters_hc_k4 clusters_hc_k5 clusters_km_k3 clusters_km_k4 clusters_km_k5 clusters_km_k6
# 1   ID_1              2              3              1              4              1              4
# 2  ID_10              1              1              3              1              4              5
# 3   ID_2              3              4              2              2              3              1
# 4   ID_3              4              5              3              3              5              3
# 5   ID_4              1              1              3              1              4              5
# 6   ID_5              1              1              3              1              4              5
# 7   ID_6              1              1              3              1              2              2
# 8   ID_7              1              2              3              1              4              5
# 9   ID_8              1              2              3              1              4              6
# 10  ID_9              1              1              3              1              4              5
###################################

################################### Change data type
str(clustering_IDs)
clustering_IDs[2:7] <- lapply(clustering_IDs[2:7], factor)
str(clustering_IDs)
clustering_IDs
#    IDs_2 clusters_hc_k4 clusters_hc_k5 clusters_km_k3 clusters_km_k4 clusters_km_k5 clusters_km_k6
# 1   ID_1              2              3              1              4              1              4
# 2  ID_10              1              1              3              1              4              5
# 3   ID_2              3              4              2              2              3              1
# 4   ID_3              4              5              3              3              5              3
# 5   ID_4              1              1              3              1              4              5
# 6   ID_5              1              1              3              1              4              5
# 7   ID_6              1              1              3              1              2              2
# 8   ID_7              1              2              3              1              4              5
# 9   ID_8              1              2              3              1              4              6
# 10  ID_9              1              1              3              1              4              5
###################################


################################### export to add to metadata_CVF
# Export 
save(clustering_IDs, file="data/clustering_IDs_by_rlog_gene_expression.RData")

# ## Loading .Rdata is from the script "Clustering_samples_on_gene_expression_and_other_data_layers.R"
# lnames = load(file="data/clustering_IDs_by_rlog_gene_expression.RData")
# lnames
################################### 


