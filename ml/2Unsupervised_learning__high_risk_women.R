# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(caret) #ML
library(randomForest)
library(rpart.plot) #Decision tree plots
library(reshape2) # heatmaps
library(BiocManager) 
library(ggrepel) # labels on ggplots
library(stats) # for Wilcoxon or other stats tests
library(cluster) # has daisy function & K-means functions
library(pheatmap) #heatmaps
library(dendextend) #Dendograms coloured
library(purrr) # exploring k values in k-means clustering
library(factoextra) #dendros - but didn't use this method in end
library(NbClust) #optimum number of clusters
library(gtsummary)   #nice tables
library(ConsensusClusterPlus) #BiocManager::install("ConsensusClusterPlus")
library(clValid) #Dunn index
library(fpc) #Cluster validation/tightness metrics

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/ML_clinical") #mac
setwd("C:/Users/alici/OneDrive - King's College London/PhD/Projects/ML_clinical") #pc
getwd()
#list.files()








###########################################################################################################################
##                                                                                                                       ##
##                                                   Load data                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Loading Rdata from the script "1Data_cleaning.R"
lnames <- load(file="data/1Data_cleaning_output.RData")
lnames




###########################################################################################################################
##                                                                                                                       ##
##                                                   df prep                                                             ##
##                                                                                                                       ##
###########################################################################################################################

############################################  Make dfs - NoOutcomes (data to cluster on) & HRClinical (all data for metadata etc)
# Looking at str
str(HR_clin)
str(outcomes_clin)
summary(outcomes_clin)
names(HR_clin)
names(outcomes_clin)

# Make NoOutcomes df - these will be the variables to cluster on
NoOutcomes <- HR_clin
NoOutcomes$PTB37 <- NULL

# Delete data in outcomes_clin that's in NoOutcomes (for a tidy merge)
outcomes_clin$Maternal.Infection <- NULL
outcomes_clin$Preeclampsia <- NULL

#Merge prep
dim(NoOutcomes) #1169   43
dim(outcomes_clin) #2059   20

# Merged data for metadata_for_heatmap
HRClinical <- merge(outcomes_clin, NoOutcomes, by="Participant.ID", all.x=F, all.y=T)
dim(HRClinical) #1169   
names(HRClinical)
row.names(HRClinical) <- paste0("ID_", HRClinical$Participant.ID, sep="")

# Delete dfs to avoid confusion
outcomes_clin <- NULL
LR_clin <- NULL

# Delete Participant.ID for clustering
NoOutcomes$Participant.ID <- NULL
str(NoOutcomes)
############################################ 


############################################ Delete duplicate data in NoOutcomes (df for clustering)
# In general I have deleted the simplified version of the duplicate options
NoOutcomes$BMI_category <- NULL
NoOutcomes$IMD_decile <- NULL

# CL
tapply(NoOutcomes$CL_minimum, NoOutcomes$Short_cervix, summary)
NoOutcomes$Short_cervix <- NULL

# Late miscarriages
tapply(NoOutcomes$Late_miscarriages_number, NoOutcomes$Previous.late.miscarriage, summary)
NoOutcomes$Previous.late.miscarriage <- NULL

# Uterine.abnormality
summary(NoOutcomes)[,12:16]
NoOutcomes$Bicornuate.uterus <- NULL
NoOutcomes$Double.cervix <- NULL
NoOutcomes$Intra.uterine.Septum <- NULL
NoOutcomes$Submucosal.fibroids <- NULL
table(NoOutcomes$Uterine.abnormality, useNA="always")

# Ethnicity
table(NoOutcomes$Ethnicity1, NoOutcomes$Ethnicity2, useNA="always")
NoOutcomes$Ethnicity2 <- NULL

############################################


############################################ Variables we'll cluster on
# Delete Centre in clustering (but want it in the demographic tables)
NoOutcomes$Centre <- NULL

# Look
names(NoOutcomes)
# [1] "Age"                          "BMI"                          "pH"                           "CL_minimum"                   "Early_miscarriages_number"   
# [6] "Late_miscarriages_number"     "IMD_rank"                     "Ethnicity1"                   "Previous.sPTB37"              "Previous.PPROM"              
# [11] "Previous.cervical.surgery"    "Uterine.abnormality"          "Cerclage"                     "Progesterone"                 "Preeclampsia"                
# [16] "Gestational.Diabetes"         "Primigravida"                 "Smoking"                      "Diagnosed.with.BV"            "History.of.BV"               
# [21] "History.of.UTIs.in.pregnancy" "History.of.GBS"               "Gender"                       "Domestic_violence"            "Recreational_drugs"          
# [26] "Amniocentesis"                "Pre.existing.hypertension"    "Asthma"                       "Type.1.diabetes"              "Type.2.diabetes"             
# [31] "Autoimmune.disease"           "Chronic.renal.disease"        "Chronic.viral.infection"  
head(NoOutcomes)
dim(NoOutcomes) #1169   33
str(NoOutcomes)

# Export for easy c&p
write.csv(t(names(NoOutcomes)), file="data/variables_for_hc_high_risk_women.csv", quote=F, row.names=F)
############################################





###########################################################################################################################
##                                                                                                                       ##
##                                                  Dissimilarity                                                        ##
##                                                                                                                       ##
###########################################################################################################################

############################################### Dissimilarity Matrix Calculation
# Info
?daisy

# Run daisy
set.seed(100)
Gower_NO_FullInfo <- daisy(NoOutcomes, metric = "gower") #Daisy/Gower scales the numeric data before calculating the dissimilarity 

# Make Gower_NO into a matrix
Gower_NO <- as.matrix(Gower_NO_FullInfo)
###############################################

############################################### sanity check for Dissimilarity Matrix Calculation
# Output most similar pair
NoOutcomes[
  which(Gower_NO == min(Gower_NO[Gower_NO != min(Gower_NO)]),
        arr.ind = TRUE)[1, ], ]

# Output most dissimilar pair
NoOutcomes[
  which(Gower_NO == max(Gower_NO[Gower_NO != max(Gower_NO)]),
        arr.ind = TRUE)[1, ], ]
#makes sense 
############################################### 


# Density-Based Spatial Clustering of Applications with Noise (DBSCAN) & k-means clustering only works with numeric data, so will stick to hc 






###########################################################################################################################
##                                                                                                                       ##
##                                       Hierarchical clustering - linkage                                               ##
##                                                                                                                       ##
###########################################################################################################################



############################################### Dendrograms - linkage methods
# Generate hclust for complete, single & average linkage methods
?hclust #Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
citation("stats")
hc_complete_Gower_NO <- hclust(as.dist(Gower_NO), method="complete")
hc_single_Gower_NO <- hclust(as.dist(Gower_NO), method="single")
hc_average_Gower_NO <- hclust(as.dist(Gower_NO), method="average")
# hc_ward_Gower_NO <- hclust(as.dist(Gower_NO), method="ward.D2") #not appropriate here as Ward does calculations in Euclidean space 

# Plot & label the 3 dendrograms looking at linkage 
pdf("plots/high_risk_women/Unsupervised_clustering/HR_Linkage_Gower.pdf", width=8, height=10)
par(mfrow = c(3,1), mar = c(5, 3, 5, 3)) #bottom, left, top, right)
plot(hc_complete_Gower_NO, main = 'Complete-linkage', labels=F)
plot(hc_single_Gower_NO, main = 'Single-linkage', labels=F)
plot(hc_average_Gower_NO, main = 'Average-linkage', labels=F) #labels=F to make IDs not printed, as they are squished at the bottom otherwise 
dev.off()

# Plot & label the 3 dendrograms looking at linkage - IDs labelled
pdf("plots/high_risk_women/Unsupervised_clustering/HR_Linkage_Gower_labels.pdf", width=10, height=12)
par(mfrow = c(3,1), mar = c(5, 3, 5, 3)) #bottom, left, top, right)
plot(hc_complete_Gower_NO, main = 'Complete-linkage')
plot(hc_single_Gower_NO, main = 'Single-linkage')
plot(hc_average_Gower_NO, main = 'Average-linkage')
dev.off()
# Complete looks best

# Individual dendro - complete linkage
pdf("plots/high_risk_women/Unsupervised_clustering/HR_Linkage_complete.pdf", width=8, height=3)
par(mfrow = c(1,1), mar = c(0, 5, 1, 0), ylab="Height") ##bottom, left, top, right)
plot(hc_complete_Gower_NO, labels=F, main=NULL)
dev.off()

# Individual dendro - average linkage
pdf("plots/high_risk_women/Unsupervised_clustering/HR_Linkage_average.pdf", width=8, height=3)
par(mfrow = c(1,1), mar = c(0, 5, 1, 0), ylab="Height") ##bottom, left, top, right)
plot(hc_average_Gower_NO, labels=F, main=NULL)
dev.off()

# Individual dendro - single linkage
pdf("plots/high_risk_women/Unsupervised_clustering/HR_Linkage_single.pdf", width=8, height=3)
par(mfrow = c(1,1), mar = c(0, 5, 1, 0), ylab="Height") ##bottom, left, top, right)
plot(hc_single_Gower_NO, labels=F, main=NULL)
dev.off()
############################################### 









###########################################################################################################################
##                                                                                                                       ##
##                                Silhouette analysis for optimal k for hc with different linkages                       ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?silhouette


############################################################### Silhouette widths for hc with "complete" method
# Set up an empty df
sil_df_HCcomplete <- data.frame(k_value=numeric(100), sil_width=numeric(100)) #whether you put 100 rows or any other number doesn't seem to impact df later as it just adds on extra rows

# Find average silhouette widths with varying value of k 
for (k in 2:100) {
  cluster_assignment <- cutree(hc_complete_Gower_NO, k=k, order_clusters_as_data = FALSE) #Assign clusters to IDs
  silhouette_values <- silhouette(cluster_assignment, dmatrix=Gower_NO) #Calculate silhouette widths
  sil_df_HCcomplete[k,1] <- k #fill column "k_value" in df
  sil_df_HCcomplete[k,2] <- summary(silhouette_values)$avg.width #fill column "sil_width" with the average silhouette widths in df
}
head(sil_df_HCcomplete)

# Delete k=0 for plotting
sil_df_HCcomplete <- sil_df_HCcomplete[-1,] #it does this because we say k=2-100, so it keeps the first row empty as there's no k=1, so the 1st row isn't filled

# Find best k
sil_df_HCcomplete[which.max(sil_df_HCcomplete$sil_width),] # 2 0.006050084

# Plot the relationship between k and sil_width
sil_width_HCcomplete <- ggplot(sil_df_HCcomplete, aes(x=k_value, y=sil_width)) +
  geom_line() +
  geom_label_repel(aes(label=ifelse(k_value==2, as.character(round(sil_width,3)), "")),
                   point.padding=0.3, box.padding = 0.5) +
  scale_x_continuous(minor_breaks = seq(0,100,10), breaks=seq(0,100,20),limits=c(1,100)) +
  scale_y_continuous(minor_breaks = seq(-0.8,0.1,0.1), breaks=seq(-0.8,0.1,0.2),limits=c(-0.6,0.1)) +
  labs(y="Mean silhouette width", x="k")
ggsave("plots/high_risk_women/Unsupervised_clustering/silhouette_hc/sil_width_HCcomplete.pdf", sil_width_HCcomplete, height=2.5, width=3.5)

# Checking - Silhouette for k=2
k2_complete <- cutree(hc_complete_Gower_NO, k=2, order_clusters_as_data = FALSE)
silhouette_k2_complete <- silhouette(k2_complete, dmatrix=Gower_NO)
plot(silhouette_k2_complete) #0.01
dev.off() 
###############################################################



############################################################### Silhouette widths for hc with "average" method
# Set up an empty df
sil_df_HCaverage <- data.frame(k_value=numeric(100), sil_width=numeric(100)) #whether you put 100 rows or any other number doesn't seem to impact df later as it just adds on extra rows

# Find average silhouette widths with varying value of k 
for (k in 2:100) {
  cluster_assignment <- cutree(hc_average_Gower_NO, k=k, order_clusters_as_data = FALSE) #Assign clusters to IDs
  silhouette_values <- silhouette(cluster_assignment, dmatrix=Gower_NO) #Calculate silhouette widths
  sil_df_HCaverage[k,1] <- k #fill column "k_value" in df
  sil_df_HCaverage[k,2] <- summary(silhouette_values)$avg.width #fill column "sil_width" with the average silhouette widths in df
}
head(sil_df_HCaverage)

# Delete k=0 for plotting
sil_df_HCaverage <- sil_df_HCaverage[-1,] #it does this because we say k=2-100, so it keeps the first row empty as there's no k=1, so the 1st row isn't filled

# Find best k
sil_df_HCaverage[which.max(sil_df_HCaverage$sil_width),] # 2 -0.1305285

# Plot the relationship between k and sil_width
sil_width_HCaverage <- ggplot(sil_df_HCaverage, aes(x=k_value, y=sil_width)) +
  geom_line() +
  geom_label_repel(aes(label=ifelse(k_value==2, as.character(round(sil_width,3)), "")),
                   point.padding=0.3, box.padding = 1) +
  scale_x_continuous(minor_breaks = seq(0,100,10), breaks=seq(0,100,20),limits=c(1,100)) +
  scale_y_continuous(minor_breaks = seq(-0.8,0.1,0.1), breaks=seq(-0.8,0.1,0.2),limits=c(-0.6,0.1)) +
  labs(y="Mean silhouette width", x="k")
ggsave("plots/high_risk_women/Unsupervised_clustering/silhouette_hc/sil_width_HCaverage.pdf", sil_width_HCaverage, height=2.5, width=3.5)

# Checking - Silhouette for k=2 
k2_average <- cutree(hc_average_Gower_NO, k=2, order_clusters_as_data = FALSE)
silhouette_k2_average <- silhouette(k2_average, dmatrix=Gower_NO)
plot(silhouette_k2_average) #-0.13
dev.off() 
###############################################################


############################################################### Selecting best k & linkage
# Note: Looks like the best option is complete linkage with k=2 

# Make vectors of cluster group (ordered by how they appear in dendrogram, rather than ID order)
final_cluster_assignment_HR_vector <- cutree(hc_complete_Gower_NO, k=2, order_clusters_as_data = FALSE) #pick which linkage method to use here
head(final_cluster_assignment_HR_vector)

# Make into a df for merging
final_cluster_assignment_HR <- as.data.frame(t(as.data.frame(t(final_cluster_assignment_HR_vector))))
final_cluster_assignment_HR <- final_cluster_assignment_HR %>% dplyr::rename("Cluster"="V1")
head(final_cluster_assignment_HR)
final_cluster_assignment_HR$Cluster <- paste0("Cluster_", final_cluster_assignment_HR$Cluster, sep="")
final_cluster_assignment_HR$Cluster <- factor(final_cluster_assignment_HR$Cluster)
final_cluster_assignment_HR$Participant.ID <- row.names(final_cluster_assignment_HR)
final_cluster_assignment_HR$Participant.ID <- gsub("ID_", "", final_cluster_assignment_HR$Participant.ID)
str(final_cluster_assignment_HR)
head(final_cluster_assignment_HR)

# Merge cluster assignment with HRClinical
HRClinical <- merge(HRClinical, final_cluster_assignment_HR, by="Participant.ID")
row.names(HRClinical) <- paste0("ID_", HRClinical$Participant.ID, sep="")
###############################################################







###########################################################################################################################
##                                                                                                                       ##
##                                        Optimal k for hc with complete linkage                                         ##
##                                                                                                                       ##
###########################################################################################################################


# Function to see optimum number of clusters - fails
?NbClust
NbClust(data = NoOutcomes, #matrix or dataset
        diss = as.dist(Gower_NO), #dissimilarity matrix to be used
        distance=NULL, #the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
        min.nc = 2, max.nc = 20, #min & max number of clusters
        method = "complete", #the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
        index = "all"#the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
        ) 
# Error in t(jeu) %*% jeu : 
#   requires numeric/complex matrix/vector arguments
# Seems like it only works when the "data" argument is all numeric


# Elbow method
?fviz_nbclust
fviz_nbclust(NoOutcomes, diss = as.dist(Gower_NO), FUNcluster=hcut, method = "wss", print.summary=T) +
  #geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method") #maybe k=3-6?
# Not sure what linkage this function is using for hclust......

# Silhouette method
fviz_nbclust(NoOutcomes, diss = as.dist(Gower_NO), FUNcluster=hcut, method = "silhouette", print.summary=T)+
  labs(subtitle = "Silhouette method") #k=2

# Gap statistic
set.seed(100)
fviz_nbclust(NoOutcomes, diss = as.dist(Gower_NO), FUNcluster=hcut, nstart = 25,  method = "gap_stat", 
             nboot = 50) + 
  labs(subtitle = "Gap statistic method") 
# Error in colMeans(x, na.rm = TRUE) : 'x' must be numeric
#So it doesn't work for mixed data types.
dev.off()













###########################################################################################################################
##                                                                                                                       ##
##                                         Dendrograms with colored bars                                                 ##
##                                                                                                                       ##
###########################################################################################################################

############################################### Make dendrogram object
# Create a dendrogram object from the hclust variable
dend_Gower_NO <- as.dendrogram(hc_complete_Gower_NO) #pick which linkage method to use here
#
# # Making dendrograms for each linkage method
# dend_Gower_NO_complete <- as.dendrogram(hc_complete_Gower_NO)
# dend_Gower_NO_average <- as.dendrogram(hc_average_Gower_NO)
# dend_Gower_NO_ward <- as.dendrogram(hc_ward_Gower_NO)
###############################################




############################################### Making colored bars for plot (barsForDendro)

# This is based on my previous script "Dendrograms with clinical bars - with early OTU matrix and late OTU matrix.R"
?colored_bars

# Selecting cols we want as bars
bars_clin <- HRClinical %>% 
  dplyr::select(Smoking, BMI_category, Ethnicity1, PTB37, Cluster)  %>% #This puts PTB37 in the first row in the heatmap
  dplyr::rename("Ethnicity"="Ethnicity1")
  
# Look
str(bars_clin)
head(bars_clin)
summary(bars_clin) #factors

# Changing "".s to "_"s
names(bars_clin) <- gsub(names(bars_clin), pattern = ".", replacement = "_", fixed=T)  
names(bars_clin)

# Duplicate df to keep one where the df is filled with actual data, not just color codes
barsForDendro <- bars_clin

# Changing all spaces & dashes in the df to underscores 
barsForDendro <- as.data.frame(apply(barsForDendro,2,function(x)gsub(' - ', '_',x)))
barsForDendro <- as.data.frame(apply(barsForDendro,2,function(x)gsub('\\s+', '_',x)))
head(barsForDendro)

# Look
str(barsForDendro)
head(barsForDendro)
summary(barsForDendro) #characters


# I will make the colours match the ones in my_colours when I used pheatmaps
my_colours

# PTB37
barsForDendro$PTB37[barsForDendro$PTB37=="Preterm"] <- "tomato"
barsForDendro$PTB37[barsForDendro$PTB37=="Term"] <- "#AAEAFF"

# Ethnicity
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="European"] <- "#FFE68A"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="African"] <- "#9976C9"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="African_Caribbean"] <- "#BE11F1"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="South_Asian"] <- "#00E67D"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="East_Asian"] <- "#1EB391"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="South_East_Asian"] <- "#097341"
barsForDendro$Ethnicity[barsForDendro$Ethnicity=="Other"] <- "#2001AE"

# BMI_category
barsForDendro$BMI_category[barsForDendro$BMI_category=="Underweight"] <- "#41D4FC"
barsForDendro$BMI_category[barsForDendro$BMI_category=="Healthy_weight"] <- "#00CC66"
barsForDendro$BMI_category[barsForDendro$BMI_category=="Overweight"] <- "yellow"
barsForDendro$BMI_category[barsForDendro$BMI_category=="Obese"] <- "#FF9933"
barsForDendro$BMI_category[barsForDendro$BMI_category=="Morbidly_obese"] <- "red"

# Smoking
barsForDendro$Smoking[barsForDendro$Smoking=="Never"] <- "#B7FF91"
barsForDendro$Smoking[barsForDendro$Smoking=="Ex_gave_up_before_pregnancy"] <- "#FDF9AE"
barsForDendro$Smoking[barsForDendro$Smoking=="Ex_gave_up_in_pregnancy"] <- "#FFCFAA"
barsForDendro$Smoking[barsForDendro$Smoking=="Current"] <- "#FE807E"

# Cluster
barsForDendro$Cluster[barsForDendro$Cluster=="Cluster_1"] <- "#C2A8FD" 
barsForDendro$Cluster[barsForDendro$Cluster=="Cluster_2"] <- "#FF519B" 


# Look
str(barsForDendro)
head(barsForDendro)

# Labels for bars for dendro
labels_barsForDendro <- gsub(names(barsForDendro), pattern = "_", replacement = " ", fixed=T)
labels_barsForDendro <- gsub(labels_barsForDendro, pattern = "PTB37", replacement = "PTB", fixed=T)
labels_barsForDendro
###############################################





############################################### Dendrograms with colored bars

# Plot dendrograms with colored bars
pdf("plots/high_risk_women/Unsupervised_clustering/Dendro_Gower_with_bars.pdf", width=15, height=8)
par(mar = c(12, 4, 3, 1), #bottom, left, top, and right
    cex=1, xpd = NA) # xpd allows content to go into outer margin
dend_Gower_NO %>% 
  #color_branches(k=4) %>%  #can use k= or h=
  set("labels_col", "white") %>% # Make IDs white font so they don't look messy
  set("labels_cex", 0.1) %>% 
  plot(horiz = F) #plot without flipping the x&y axes
title("Hierarchial clustering (complete linkage) of high risk women") 
colored_bars(colors = barsForDendro, #Colours I chose
             dend=dend_Gower_NO, 
             rowLabels = labels_barsForDendro, # Nice labels vector
             cex.rowLabels=0.7,
             horiz = F) #plot without flipping the x&y axes
dev.off()
############################################### 


############################################### Legend of what the colours mean in the dendrograms
# Can do later if needed. Because a legend was automatically made for Daisy_clustering_on_clinical_data__Gower.pdf & I could just copy this 
############################################### 




############################################### Checking dendros aligned correctly to clinical data
# Checking first ID on dendro
head(final_cluster_assignment_HR_vector, 1)
HRClinical %>% 
  filter(Participant.ID=="XXXX") %>% 
  dplyr::select(Participant.ID, Cluster, PTB37, Ethnicity1, BMI_category, Smoking)  



# Checking last ID on dendro
tail(final_cluster_assignment_HR_vector, 1)
HRClinical %>% 
  filter(Participant.ID=="XXXX") %>% 
  dplyr::select(Participant.ID, Cluster, PTB37, Ethnicity1, BMI_category, Smoking)  
############################################### 








###########################################################################################################################
##                                                                                                                       ##
##                                                   hc clusters & metadata                                              ##
##                                                                                                                       ##
###########################################################################################################################


#Look
names(HRClinical)
str(HRClinical)
head(HRClinical)
summary(HRClinical)

# Info
?gtsummary::tbl_summary # I had dependency issues with rlang, vctrs, cli, magrittr, xfun. Had to delete & reinstall most of these for gtsummary::tbl_summary() to work.
?gtsummary::add_p() 
?gt::gtsave




# Table of outcomes (factors which were NOT clustered on) 
Table_of_Outcomes_by_Cluster <- HRClinical %>% 
  dplyr::select(Gestation.at.delivery.wks.dec, PTB37, PTB_type, sPTB37, sPTB34,   #PTB/TB
                Chorioamnionitis, Apgar.Score.1.min, Apgar.Score.5.min, Major.congenital.abnormality, Neonatal.death, #baby health
                Customised.birthweight.centiles, SGA, Suspected.fetal.growth.restriction, #baby size
                Premature.Prelabour.Rupture.of.Membranes, Mode.of.delivery, Pregnancy.Outcome, Pregnancy.Outcome.Status, #less important ones
                Cluster) %>% 
  gtsummary::tbl_summary(     
    by = Cluster,                                               # stratify entire table by Cluster
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = all_categorical() ~ "categorical",                 # force all categorical levels to display
    missing_text = "Missing",                                   # how missing values should display
    label  = list(Gestation.at.delivery.wks.dec ~ "Gestation at delivery (weeks)", # display labels for column names
                  PTB37 ~ "PTB <37 weeks",
                  sPTB37 ~ "sPTB <37 weeks",
                  sPTB34 ~ "sPTB <34 weeks",
                  Customised.birthweight.centiles ~ "Customised birthweight centiles",
                  Apgar.Score.1.min ~ "Apgar score (1 min)",
                  Apgar.Score.5.min ~ "Apgar score (5 min)",
                  Major.congenital.abnormality ~ "Major congenital abnormality",
                  Neonatal.death ~ "Neonatal death",
                  Suspected.fetal.growth.restriction ~ "Suspected fetal growth restriction",
                  SGA ~ "Small for gestational age")
    ) %>%                              
  add_p() %>%
  # Export 
  as_gt() %>%
  gt::gtsave(filename = "data/tables_exported/Clusters_HR_PregOutcomes.tex")  # use extensions .html .tex .ltx .rtf
  






# Table of input variables (factors which were clustered on directly or indirectly) 
Table_of_hcFactors_by_Cluster <- HRClinical %>% 
  dplyr::select(Previous.sPTB37, Previous.PPROM, Previous.late.miscarriage, Previous.cervical.surgery, Uterine.abnormality, #make up Risk
                Early_miscarriages_number, Late_miscarriages_number, #number of miscarriages
                CL_minimum, 
                Cerclage, Progesterone, #intervention
                Preeclampsia, Pre.existing.hypertension, Gestational.Diabetes, Type.1.diabetes, Type.2.diabetes, Asthma, Autoimmune.disease, Chronic.renal.disease, Chronic.viral.infection, #medical conditions 
                Domestic_violence,
                Amniocentesis, #medical procedure
                pH, Diagnosed.with.BV, History.of.BV, History.of.UTIs.in.pregnancy, History.of.GBS, #bacterial
                Ethnicity1, Ethnicity2, Age, IMD_rank, IMD_decile, BMI, BMI_category, Primigravida, Gender, Smoking, Recreational_drugs, #demographics
                Cluster) %>% 
  tbl_summary(     
    by = Cluster,                                               # stratify entire table by Cluster
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(all_categorical() ~ "categorical",            # force all categorical levels to display
                  Late_miscarriages_number ~ "continuous"),     # force the autoconvert of INT to FACTOR for this variable https://www.pipinghotdata.com/posts/2021-07-14-polished-summary-tables-in-r-with-gtsummary/           
    label  = list(                                              # display labels for column names
      Previous.sPTB37 ~ "Previous sPTB <37 weeks",
      Previous.PPROM ~ "Previous PPROM",
      Previous.late.miscarriage ~ "Previous late miscarriage",
      Previous.cervical.surgery ~ "Previous cervical surgery",
      Uterine.abnormality ~ "Uterine abnormality",
      Early_miscarriages_number ~ "Number of early miscarriages",
      Late_miscarriages_number ~ "Number of late miscarriages",
      CL_minimum ~ "Minimum cervical length (mm)",
      Pre.existing.hypertension ~ "Pre-existing hypertension",
      Gestational.Diabetes ~ "Gestational diabetes",
      Type.1.diabetes ~ "Type 1 diabetes",
      Type.2.diabetes ~ "Type 2 diabetes",
      Autoimmune.disease ~ "Autoimmune disease",
      Domestic_violence ~ "Domestic violence",
      pH ~ "High vaginal pH",
      Diagnosed.with.BV ~ "Diagnosed with BV",
      History.of.BV ~ "History of BV",
      History.of.UTIs.in.pregnancy ~ "History of UTIs in pregnancy",
      History.of.GBS ~ "History of group B streptococcus (GBS) infection",
      Ethnicity1 ~ "Ethnicity",
      Gender ~ "Gender of baby",
      BMI_category ~ "Body mass index classification",
      IMD_decile ~ "Index of Multiple Deprivation (IMD) decile",
      Recreational_drugs ~ "Recreational drugs"
      ),
    missing_text = "Missing") %>%                               # how missing values should display
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% #fixing Fisher test error for Ethnicity1 https://github.com/ddsjoberg/gtsummary/issues/952
  # Export
  as_gt() %>%
  gt::gtsave(filename = "data/tables_exported/Clusters_HR_hcFactors.tex")  # use extensions .html .tex .ltx .rtf



# pH - Looking & Checking mean & Wilcoxon p work
tapply(HRClinical$pH, HRClinical$Cluster, summary)
ggpubr::compare_means(pH ~ Cluster, data=HRClinical, method = "wilcox.test", p.adjust.method = "fdr") # matches the table  
ggplot(data=HRClinical, aes(x=Cluster, y=pH)) + geom_boxplot()
dev.off()
Cluster_by_pH <- ggplot(data=HRClinical, aes(x=Cluster, y=pH)) + geom_violin()
ggsave("plots/high_risk_women/Unsupervised_clustering/Cluster_by_features/Cluster_by_pH.pdf", Cluster_by_pH, height=2.5, width=3)


# # Looking at sPTB & cervical surgery in Cluster_2
# table(HRClinical$Previous.cervical.surgery[HRClinical$Cluster=="Cluster_2"], 
#       HRClinical$PTB37[HRClinical$Cluster=="Cluster_2"], useNA="always") 
# # In women in Cluster 2: 56 had sPTB without cervical surgery, 19 had sPTB w/ surgery, 272 had surgery but not sPTB

# When is the test Chi-squared or Fisher:
# Default is "`chisq.test`" when expected cell counts >=5 and "`fisher.test`" when expected cell counts <5.
# https://github.com/ddsjoberg/gtsummary/blob/main/R/add_p.R








###########################################################################################################################
##                                                                                                                       ##
##                                                     Exporting                                                         ##
##                                                                                                                       ##
###########################################################################################################################

# Save main dfs
save(final_cluster_assignment_HR, #df of final Clusters for IDs (Gower, hc, complete, k=2)
     HRClinical, #df of clinical input, clinical outcomes & Cluster group
     file="data/2Unsupervised_learning_output_high_risk_women.RData")


# # Loading Rdata from the script "2Unsupervised_learning__high_risk_women.R"
# lnames <- load(file="data/2Unsupervised_learning_output_high_risk_women.RData")
# lnames




