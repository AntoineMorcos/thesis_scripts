# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats)
library(DESeq2) #for assay()
library(vegan) #for Hellinger
library(WGCNA)
library(GO.db)
library(org.Hs.eg.db) # Human genes

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()
#list.files()






###########################################################################################################################
##                                                                                                                       ##
##                                                  Load data                                                            ##
##                                                                                                                       ##
###########################################################################################################################

# Loading metagenome .Rdata created in the script "Metagenome_exploration.R"
lnames = load(file="data/metagenome_counts_and_RA.RData")
lnames

# Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
lnames = load(file="data/immune_CytokinesComplementsIgM_2023.RData")
lnames

## Loading metabolite .Rdata created in the script "Metabolites_exploration.R"
lnames = load(file="data/metabolites_df.RData")
lnames





###########################################################################################################################
##                                                                                                                       ##
##                                  Prep needed before making datExpr & datTraits                                        ##
##                                                                                                                       ##
###########################################################################################################################

########################################### datExpr prep
# Looking
dim(MG_Hel_t) #87 230
MG_Hel_t[1:5,1:5]
#summary(MG_Hel)
##########################################



########################################## datTraits prep

# Making ID1 from the Sample col
immu_MD$ID1 <- gsub("^","ID_", immu_MD$Participant.ID) # put ID_ at the beginning 
immu_MD$ID1
immu_MD %>% dplyr::select(ID1,Sample,Participant.ID)


# Setting rownames 
rownames(immu_MD) <- immu_MD$ID1

# Reorder df by ID - although it looks like it already is
immu_MD <- immu_MD %>% 
  dplyr::arrange(ID1)
immu_MD[1:5,1:5]
immu_MD$ID1

# Making a gestation in delivery weeks in decimals col & the same with visit
immu_MD$Gestation.at.delivery.wks.dec
immu_MD <- immu_MD %>%
  dplyr::rename(ga_birth_w_dec=Gestation.at.delivery.wks.dec) %>% 
  mutate(ga_visit_w_dec = (ga_visit_w + (ga_visit_d/7))) 
  

# Checks
plot(immu_MD$ga_birth_w_dec, immu_MD$ga_w)
plot(immu_MD$ga_visit_w_dec, immu_MD$ga_visit_w)
dev.off()
immu_MD %>% dplyr::select(ga_birth_w_dec, ga_w, ga_d, ga_visit_w_dec, ga_visit_w, ga_visit_d)


# Clinical quantitative data in immu_MD (not including metabolites & cytokines)
names(immu_MD)
clinical_quant <- immu_MD %>%
  dplyr::select(ID1, id,# for merges
                Age, BMI, IMD_rank,
                Short_cervix_value_this_visit, CL_absolute_small_value_at_any_visit, 
                ga_birth_w_dec, ga_visit_w_dec,
                medscinet_ph) 
head(clinical_quant)
str(clinical_quant)

# Change str
clinical_quant$id <- as.character(clinical_quant$id)
clinical_quant$Age <- as.numeric(clinical_quant$Age)
clinical_quant$CL_absolute_small_value_at_any_visit <- as.numeric(clinical_quant$CL_absolute_small_value_at_any_visit)
clinical_quant$Short_cervix_value_this_visit <- as.numeric(clinical_quant$Short_cervix_value_this_visit)
str(clinical_quant)


# Metabolite data - look
names(metabolites_df) 
str(metabolites_df) #all num

# Adding ID1 col for merge & move to beginning of df
metabolites_df2 <- metabolites_df
metabolites_df2$ID1 <- paste("ID_", rownames(metabolites_df2), sep="")
metabolites_df2 <- metabolites_df2 %>%
  dplyr::select(ID1, everything())
metabolites_df2[1:5,1:5]
metabolites_df2$Participant.ID <- NULL
rownames(metabolites_df2) <- metabolites_df2$ID1 
metabolites_df2[1:5,1:5]

# Merge prep checks
head(clinical_quant)
head(metabolites_df2)
dim(clinical_quant) #87  9
dim(metabolites_df2) #87 30 (was 44 mets when unknown metabolites & Propelene.glycol were included)
table(rownames(clinical_quant) == rownames(metabolites_df2), useNA="always") # All true

# Merge clinical_quant with metabolites_df2
clinical_and_metabolites <- merge(clinical_quant, metabolites_df2, by="ID1", all.x=T)
dim(clinical_and_metabolites) #87 39
head(clinical_and_metabolites)
#View(clinical_and_metabolites)



######## immu
# Looking at ID col
immu_small$Participant.ID
clinical_and_metabolites$id

# Looking at numbers in common and distinct
dim(immu_small) #81 women with cytokine data, but 87 women with metagenome 
table(immu_small$Participant.ID %in% clinical_and_metabolites$id) #TRUE: 81
table(clinical_and_metabolites$id %in% immu_small$Participant.ID) #FALSE:  6,  TRUE: 81

# Remove ID not in clinical_and_metabolites$id
immu_metagIDs <- as.data.frame(data.table::setDT(immu_small)[Participant.ID %chin% as.character(clinical_and_metabolites$id)]) 
dim(immu_metagIDs) # 81 women x 14 cytokines/complements/IgM
head(immu_metagIDs)
#######

# Merge prep checks
head(clinical_and_metabolites)
head(immu_metagIDs)
dim(clinical_and_metabolites) #87  39
dim(immu_metagIDs) #81 14
table(rownames(clinical_and_metabolites) == rownames(immu_metagIDs), useNA="always") # 6 FALSE as we don't have cytokine data on 6 women

# Make col to merge on have the same name
clinical_and_metabolites <- clinical_and_metabolites %>% dplyr::rename("Participant.ID" ="id")
head(clinical_and_metabolites)

# Merge clinical_and_metabolites with cytokine data
allTraits <- merge(clinical_and_metabolites, immu_metagIDs, by="Participant.ID", all.x=T)
dim(allTraits) #87 52
#View(allTraits) 



# Change Rownames
rownames(allTraits) <- allTraits$ID1
rownames(allTraits)

# Checks
head(allTraits)
allTraits$IMD_rank <- as.numeric(allTraits$IMD_rank)
str(allTraits) #All are num, as I've already changed them

# Remove extra ID column
allTraits$Participant.ID <- NULL

# Look
head(allTraits)

# Rename some cols
allTraits <- allTraits %>% 
  dplyr::rename(pH=medscinet_ph,
                Gestation_at_birth=ga_birth_w_dec,
                Gestation_at_visit=ga_visit_w_dec,
                CL_at_visit=Short_cervix_value_this_visit,
                CL_minimum=CL_absolute_small_value_at_any_visit)
head(allTraits)
names(allTraits)

##########################################





###########################################################################################################################
##                                                                                                                       ##
##                                WGCNA step 1 - Data input, cleaning and pre-processing                                 ##
##                                                                                                                       ##
###########################################################################################################################

# The WGCNA PDF tutorials, which were used to create this script, are here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html 


## 1.Data input, cleaning and pre-processing

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# In datExpr the rows should correspond to samples and the columns correspond to genes
datExpr <- MG_Hel_t

# Checks
datExpr[1:5, 1:5]
rownames(datExpr) #IDs
dim(datExpr) #[1]  87 230

# 1.b Checking data for excessive missing values and identification of outlier microarray samples
gsg = goodSamplesGenes(datExpr, verbose = 3) 
# Flagging genes and samples with too many missing values...
# ..step 1
# ..Excluding 3 genes from the calculation due to too many missing samples or zero variance.
# ..step 2
gsg$allOK #  FALSE


# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
#Removing genes: Eubacterium_saburreum, Granulicatella_elegans, Lactobacillus_ruminis

dim(datExpr) # 87 227 #3 species gone

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr), method = "average")

# Plot
pdf("plots/WGCNA/Hel/mg_sampleClustering.pdf", width=8, height=4)
par(cex = 0.6, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()



# 1.c Loading clinical trait data

# We now read in the trait data and match the samples for which they were measured to the expression samples.
# remove columns that hold information we do not need / choose columns we want #we did this earlier in the prep section
# 
# Looking
head(allTraits)
str(allTraits)
dim(allTraits) #87 51
names(allTraits)

 
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples_1 = rownames(datExpr)
traitRows = match(Samples_1, allTraits$ID1)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits.
# Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
?plotDendroAndColors
pdf("plots/WGCNA/Hel/mg_sampleClustering_datTraits.pdf", width=8, height=8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = NULL, #"Sample dendrogram and trait heatmap",
                    cex.colorLabels= 0.6, cex.dendroLabels=0.4, 
                    dendroLabels=F, #hide sample IDs from thesis
                    marAll=c(0,4,1,0)) #these margins work #bottom, left, top, right
dev.off()
# In the above plot white means a low value, red a high value, and grey a missing entry.





###########################################################################################################################
##                                                                                                                       ##
##                     WGCNA step 2b - Step-by-step network construction and module detection                            ##
##                                                                                                                       ##
###########################################################################################################################

# 2.b Step-by-step network construction and module detection
# 2.b.1 Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
powers

# Call the network topology analysis function
?pickSoftThreshold
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 2, networkType="signed") # looks like 12 is best sft
# pickSoftThreshold: will use block size 227.
# pickSoftThreshold: calculating connectivity for given powers...
# ..working on genes 1 through 227 of 227
#    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1  0.23300  4.360         0.7630 123.000   123.000 142.00
# 2      2  0.00423  0.211         0.8020  69.300    69.400  94.10
# 3      3  0.01390 -0.267         0.6660  41.000    40.400  65.60
# 4      4  0.37100 -0.977         0.8130  25.600    24.800  47.70
# 5      5  0.61700 -1.100         0.9360  16.800    16.100  36.40
# 6      6  0.64800 -1.120         0.8860  11.700    11.000  28.80
# 7      7  0.57600 -1.120         0.7510   8.550     7.980  23.30
# 8      8  0.51900 -1.110         0.6520   6.520     5.730  19.20
# 9      9  0.59700 -0.988         0.8130   5.160     4.570  16.10
# 10    10  0.68400 -0.945         0.8860   4.210     3.720  13.60
# 11    12  0.83000 -0.743         0.9000   3.020     2.570  10.00 #max here 
# 12    14  0.78300 -0.873         0.7900   2.320     1.690   8.72
# 13    16  0.82100 -0.911         0.8110   1.880     1.140   8.02
# 14    18  0.77600 -1.050         0.7240   1.580     0.867   7.49
# 15    20  0.78100 -1.090         0.7240   1.370     0.651   7.11
# 16    22  0.81000 -1.080         0.7560   1.210     0.484   6.89
# 17    24  0.76200 -1.160         0.6940   1.090     0.360   6.76
# 18    26  0.77200 -1.130         0.7070   0.998     0.282   6.64
# 19    28  0.17600 -2.090        -0.0375   0.922     0.218   6.53
# 20    30  0.18000 -2.110        -0.0321   0.860     0.169   6.43



# Plot the soft thresholding results - signed:
pdf("plots/WGCNA/Hel/mg_soft_thresholding_power_signed.pdf", width=8, height=4)
par(mfrow = c(1,2), #1x2 panels
    cex = 0.9) #scaling text size
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red") #cex=cex1 is making it error, so have removed it from this line...
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red") #cex=cex1 is making it error, so have removed it from this line...
dev.off()


# Plot the soft thresholding results - PART A only!!
pdf("plots/WGCNA/Hel/mg_soft_thresholding_power_signed_Rsquared.pdf", width=3.5, height=3.5)
par(mfrow = c(1,1), #1x1 panel
    cex = 0.8, mar = c(5,5,1,1)) #bottom, left, top, right 
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft-thresholding power",ylab=expression(paste("Scale free topology model fit (R"^"2", ")")),type="n",
     main = NULL)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red") #cex=cex1 is making it error, so have removed it from this line...
dev.off()
# superscript in axes labels https://stackoverflow.com/questions/10628547/use-superscripts-in-r-axis-labels


# Plot the soft thresholding results - PART B only!!
pdf("plots/WGCNA/Hel/mg_soft_thresholding_power_signed_mean_connectivity.pdf", width=3.5, height=3.5)
par(mfrow = c(1,1), #1x1 panel
    cex = 0.8, mar = c(5,5,1,1)) #bottom, left, top, right 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft-thresholding power",ylab="Mean connectivity", type="n",
     main = NULL)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red") #cex=cex1 is making it error, so have removed it from this line...
dev.off()






# 2.b.2 Co-expression similarity and adjacency
# We now calculate the adjacencies, using the selected soft thresholding power:
softPower = 12 # sft 
adjacency = adjacency(datExpr, power = softPower, type = "signed")



# 2.b.3 Topological Overlap Matrix (TOM)
# To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
# Turn adjacency into topological overlap
Sys.time()
?TOMsimilarity
TOM = TOMsimilarity(adjacency) 
dissTOM = 1-TOM
Sys.time() # <1min for metagenome




# 2.b.4 Clustering using TOM
# We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
# Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf("plots/WGCNA/Hel/mg_Gene_clustering.pdf", width=8, height=5)
par(mfrow = c(1,1), #1x1 panels
    cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(geneTree, xlab="", 
     sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



# Set the minimum module size
minModuleSize = 5 

# Module identification using dynamic tree cut:
?cutreeDynamic
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# ..cutHeight not given, setting it to 0.992  ===>  99% of the (truncated) height range in dendro.
# ..done.

# Look at number of modules and their membership sizes
table(dynamicMods) 
# The function returned modules from largest membership size to smallest. 
# Label 0 is reserved for unassigned genes. 


#We now plot the module assignment under the gene dendrogram:

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) #These are the unmerged modules

# Plot the dendrogram and colors underneath
pdf("plots/WGCNA/Hel/mg_Gene_clustering__unmerged_modules.pdf", width=12, height=6)
par(cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and first module colours")
dev.off()




# 2.b.5 Merging of modules whose expression profiles are very similar
# The Dynamic Tree Cut may identify modules whose expression profiles are very similar.
# It may be prudent to merge such modules since their genes are highly co-expressed.
# To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-WGCNA::cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf("plots/WGCNA/Hel/mg_ME_clustering.pdf", width=7, height=5)
par(cex = 0.9, #scaling text size
    mar = c(0.5,5,1,0)) #bottom, left, top, right
plot(METree, main = NULL,
     xlab = "", sub = "")
dev.off()


#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25 

# # Plot the result
# pdf("plots/WGCNA/Hel/mg_ME_clustering_with_diss_threshold.pdf", width=8, height=6)
# par(cex = 0.9, #scaling text size
#     mar = c(1,5,4,0)) #bottom, left, top, right
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# # Plot the cut line into the dendrogram
# abline(h=MEDissThres, col = "red")
# dev.off()
# no merging as all correlations are <0.75



# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors
table(mergedColors) # which colours
dim(table(mergedColors)) # number of colours 14

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
pdf("plots/WGCNA/Hel/mg_Gene_clustering__unmerged_and_merged_modules.pdf", width=12, height=9)
par(cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="Gene dendrogram comparing merged and unmerged modules with cut height of 0.25")
dev.off()


# Plot the dendrogram and colors underneath
?plotDendroAndColors
pdf("plots/WGCNA/Hel/mg_Gene_clustering__final_modules.pdf", width=8, height=4.5)
par(cex = 0.9)
plotDendroAndColors(geneTree, mergedColors, "Final module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = NULL,#"Gene dendrogram and final module colours"
                    marAll =c(0.5,4.5,0.5,0)) #bottom, left, top, right
dev.off()



#In the subsequent analysis, we will use the merged module colors in mergedColors. We save the relevant variables for use in subsequent parts of the tutorial:
## Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
table(moduleColors)
# moduleColors
# black        blue       brown       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#    15          36          19          17           6          13          10          15           8          17           6           6          42          17
dim(table(moduleColors)) # 14 before & after


# ################################################## What taxa to what module - can comment out so don't accidentally overwrite csv if we change parameters 
# # Looking at genes in 1 module
# names(datExpr)[moduleColors=="grey"]
# 
# # Make the df
# taxa_module_assignment <- as.data.frame(cbind(names(datExpr), moduleColors)) %>%
#   dplyr::rename(Taxa=V1,
#                 Module = moduleColors)
# 
# # Check
# taxa_module_assignment %>% filter(Module=="grey") #matching 
# 
# # Export to csv
# write_csv(taxa_module_assignment, "data/WGCNA/taxa_module_assignment_Hel_sft10_mm5_2023.csv")
# 
# # I checked they match the cytoscape data and they do 
# #diff taxa_module_assignment_Hel_sft10_mm5_2023.csv taxa_module_assignment_Hel_sft10_mm5_June2022.csv
# # 1c1
# # < Taxa,Module
# # ---
# #   > Taxa,Module_when_5minModuleSize
# 
# ##################################################






###########################################################################################################################
##                                                                                                                       ##
##                          WGCNA step 3 - Relating modules to external clinical trait                                   ##
##                                                                                                                       ##
###########################################################################################################################


#3 Relating modules to external clinical traits 

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
?moduleEigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Calculating ME correlations and p-values with datTraits
?WGCNA::cor
moduleTraitCor = WGCNA::cor(MEs, datTraits, use = "p") #the correlation values used in the heatmap
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) # the p-values used in the heatmap

#testing what the default cor is 
test_pear = WGCNA::cor(MEs, datTraits, use = "p", method = "pearson")
test_kend = WGCNA::cor(MEs, datTraits, use = "p", method = "kendall")
test_spear = WGCNA::cor(MEs, datTraits, use = "p", method = "spearman")
moduleTraitCor #pearson

#Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# # Display the correlation values within a heatmap plot
# pdf("plots/WGCNA/Hel/mg_WGCNA_minmod3_heatmap_cor.pdf", width=12, height=5)
# par(cex = 0.6, #scaling text size
#     mar = c(9, 10, 3, 2)) #bottom, left, top, right
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.3,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships (sft=12, mm=3)"))
# dev.off()



############################################################### Altered version of heatmap - where only sig values are shown in text
##This is based on Flavia's code

# Run the plot but maintaining text only for pvalues <0.05
moduleTraitPvalue_0.05 <- moduleTraitPvalue
moduleTraitPvalue_0.05[moduleTraitPvalue_0.05 >=0.05 ] <- NA
moduleTraitPvalue_0.05
moduleTraitPvalue_0.05 <- formatC(moduleTraitPvalue_0.05, format = "e", digits = 2)
moduleTraitPvalue_0.05[moduleTraitPvalue_0.05 == " NA" ] <- " "
moduleTraitPvalue_0.05[moduleTraitPvalue_0.05 == "NaN" ] <- " "
moduleTraitPvalue_0.05[1:5,1:5] # Checks

# making the p-values in brackets on a new line
p4heatmap <- moduleTraitPvalue_0.05
p4heatmap <- base::ifelse(p4heatmap==" ",p4heatmap,paste("\n(", p4heatmap, ")", sep = ""))
p4heatmap[1:5,1:5] # 
moduleTraitPvalue_0.05[1:5,1:5] # matching

## Getting correlations on the plot too, not just p-values
moduleTraitCor_0.05 <- moduleTraitCor # Duplicating
moduleTraitCor_0.05[moduleTraitPvalue_0.05 ==" " ] <- NA # All blanks in the pvalue matrix turn to NA in the correlation matrix
moduleTraitCor_0.05[1:5,1:5] # Checks
# moduleTraitPvalue_0.05[1:5,1:5] # Checks # Blank in same places 
moduleTraitCor_0.05 <- formatC(moduleTraitCor_0.05, format = "f", digits = 2)
moduleTraitCor_0.05[moduleTraitCor_0.05 == " NA" ] <- " "
moduleTraitCor_0.05[moduleTraitCor_0.05 == "NaN" ] <- " "
moduleTraitCor_0.05[1:5,1:5] # Checks

#Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix_0.05 = paste(moduleTraitCor_0.05, p4heatmap, sep = "")
dim(textMatrix_0.05) = dim(moduleTraitCor_0.05)

# Display the correlation values within a heatmap plot - with p-values
?labeledHeatmap
pdf("plots/WGCNA/Hel/mg_WGCNA_minmod5_heatmap_cor_p0.05.pdf", width=12, height=4.5)
par(cex = 0.6, #scaling text size
    mar = c(9, 10, 3, 2)) #bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_0.05,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
###############################################################
 




############################################################### Altered version of heatmap - where only sig values are shown in text & only some metadata is plotted

################## 1. Making a reduced version of datTraits for clearer heatmap
str(datTraits)
names(datTraits)
datTraits_small <- datTraits %>% dplyr::select(Age, BMI, CL_at_visit, CL_minimum, Gestation_at_birth, IMD_rank, pH, 
                                               Acetate, Alanine, Asparagine, Aspartate, Cadaverine, Formate, Glucose, Glutamate, Glutamine, Isoleucine, Lactate, Leucine, Lysine, 
                                               Maltose, Methionine, Phenylalanine, Putrescine, Succinate, Taurine, Threonine, Tryptophan, Tyramine, Tyrosine, Uracil, Valine, 
                                               C5, GCSF, IgM, IL10, IL13, IL17A, IL22, IL5)
# Columns not in heatmap 
setdiff(names(datTraits), names(datTraits_small)) # returns values that are only in the first vector
# Gestation_at_visit, Betaine, Choline, Ethanol, Pyruvate, C5a, IFNg, IL6, IL8, MCP1. 
##################

################## 2. Getting rid of grey module in heatmap
str(MEs)
names(MEs)
MEs_noGrey <- MEs
MEs_noGrey$MEgrey <- NULL
##################

##################  3. Calculating correlation & p-values
# Calculating ME correlations and p-values with datTraits_small
moduleTraitCor_smallTraits = WGCNA::cor(MEs_noGrey, datTraits_small, use = "p") #the correlation values used in the heatmap
moduleTraitPvalue_smallTraits = corPvalueStudent(moduleTraitCor_smallTraits, nSamples) # the p-values used in the heatmap

# Matrix of correlations and their p-values
textMatrix_smallTraits = paste(signif(moduleTraitCor_smallTraits, 2), "\n(",
                               signif(moduleTraitPvalue_smallTraits, 1), ")", sep = "")
dim(textMatrix_smallTraits) = dim(moduleTraitCor_smallTraits)
################## 


################## 4. Working with text overlay for plot - p-values
# Only keep test when pvalues <0.05
moduleTraitPvalue_smallTraits_0.05 <- moduleTraitPvalue_smallTraits
moduleTraitPvalue_smallTraits_0.05[moduleTraitPvalue_smallTraits_0.05 >=0.05 ] <- NA
moduleTraitPvalue_smallTraits_0.05
moduleTraitPvalue_smallTraits_0.05 <- formatC(moduleTraitPvalue_smallTraits_0.05, format = "e", digits = 2)
moduleTraitPvalue_smallTraits_0.05[moduleTraitPvalue_smallTraits_0.05 == " NA" ] <- " "
moduleTraitPvalue_smallTraits_0.05[moduleTraitPvalue_smallTraits_0.05 == "NaN" ] <- " "
moduleTraitPvalue_smallTraits_0.05[1:5,1:5] # Checks

# making the p-values in brackets on a new line
p4heatmap_smallTraits <- moduleTraitPvalue_smallTraits_0.05
p4heatmap_smallTraits <- base::ifelse(p4heatmap_smallTraits==" ",p4heatmap_smallTraits,paste("\n(", p4heatmap_smallTraits, ")", sep = ""))
p4heatmap_smallTraits[1:5,1:5] # 
moduleTraitPvalue_smallTraits_0.05[1:5,1:5] # matching
##################

################## 5. Working with text overlay for plot - correlations
## Getting correlations on the plot too, not just p-values
moduleTraitCor_smallTraits_0.05 <- moduleTraitCor_smallTraits # Duplicating
moduleTraitCor_smallTraits_0.05[moduleTraitPvalue_smallTraits_0.05 ==" " ] <- NA # All blanks in the pvalue matrix turn to NA in the correlation matrix
moduleTraitCor_smallTraits_0.05[1:5,1:5] # Checks
# moduleTraitPvalue_smallTraits_0.05[1:5,1:5] # Checks # Blank in same places 
moduleTraitCor_smallTraits_0.05 <- formatC(moduleTraitCor_smallTraits_0.05, format = "f", digits = 2)
moduleTraitCor_smallTraits_0.05[moduleTraitCor_smallTraits_0.05 == " NA" ] <- " "
moduleTraitCor_smallTraits_0.05[moduleTraitCor_smallTraits_0.05 == "NaN" ] <- " "
moduleTraitCor_smallTraits_0.05[1:5,1:5] # Checks

# Will display correlations and their p-values
textMatrix_smallTraits_0.05 = paste(moduleTraitCor_smallTraits_0.05, p4heatmap_smallTraits, sep = "")
dim(textMatrix_smallTraits_0.05) = dim(moduleTraitCor_smallTraits_0.05)
################## 


################## 6. Making prettier x-axis labels
# Making nice axis labels
x_axis_labels_smallTraits <- names(as.data.frame(moduleTraitCor_smallTraits_0.05))
x_axis_labels_smallTraits
x_axis_labels_smallTraits <- gsub("_"," ", x_axis_labels_smallTraits)
x_axis_labels_smallTraits <- gsub("IL","IL-", x_axis_labels_smallTraits)
x_axis_labels_smallTraits <- gsub("GCSF","G-CSF", x_axis_labels_smallTraits)
# x_axis_labels_smallTraits <- gsub("IFNg","IFN-γ", x_axis_labels_smallTraits) # not in heatmap anyway...
x_axis_labels_smallTraits <- gsub("MCP1","MCP-1", x_axis_labels_smallTraits)
x_axis_labels_smallTraits <- gsub("CL ","Cervical length ", x_axis_labels_smallTraits)
x_axis_labels_smallTraits <- gsub("Cervical length minimum","Minimum cervical length ", x_axis_labels_smallTraits)
#x_axis_labels_smallTraits <- gsub("."," ", x_axis_labels_smallTraits, fixed = TRUE) #without fixed argument it replaces everything lol
x_axis_labels_smallTraits
################## 

################## 7. Plotting
# Display the correlation values within a heatmap plot - with p-values
pdf("plots/WGCNA/Hel/mg_WGCNA_minmod5_heatmap_cor_p0.05_OnlySigTraitsLeft.pdf", width=10, height=4)
par(cex = 0.6, #scaling text size
    mar = c(7.5, 7.5, 0.5, 0.5)) #bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor_smallTraits,
               xLabels = x_axis_labels_smallTraits,
               yLabels = names(MEs_noGrey),
               ySymbols = names(MEs_noGrey),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_smallTraits_0.05,
               setStdMargins = FALSE,
               cex.text = 0.35,
               zlim = c(-1,1))
dev.off()
##################
###############################################################






###########################################################################################################################
##                                                                                                                       ##
##                               Exporting a gene network to external visualization software                             ##
##                                                                                                                       ##
###########################################################################################################################

## Exporting to Cytoscape

##################################################################### All modules (excluding grey)
# Select modules (WE WANT ALL apart from grey)
all_mods <- names(table(moduleColors))
all_mods
all_mods <- all_mods[!all_mods %in% c("grey")]
all_mods

# Select module genes
genes = names(datExpr)
Modules_in_all_mods = is.finite(match(moduleColors, all_mods))
Genes_in_all_mods = genes[Modules_in_all_mods]

# Looking
str(all_mods) #chr [1:13]
table(Modules_in_all_mods)
# Modules_in_all_mods
# FALSE  TRUE 
#    13   214 
length(Genes_in_all_mods) #214

# Select the corresponding Topological Overlap
modTOM_in_all_mods = TOM[Modules_in_all_mods, Modules_in_all_mods]
dimnames(modTOM_in_all_mods) = list(Genes_in_all_mods, Genes_in_all_mods)

# Export the network into edge and node list files Cytoscape can read
?exportNetworkToCytoscape
cyto_all_mods = exportNetworkToCytoscape(modTOM_in_all_mods,
                                         edgeFile = "Cytoscape/CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_2023.txt",
                                         nodeFile = "Cytoscape/CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_2023.txt",
                                         weighted = TRUE,
                                         threshold = 0.05, #adjacency threshold for including edges in the output - may be worth playing with this value
                                         nodeNames = Genes_in_all_mods,
                                         nodeAttr = moduleColors[Modules_in_all_mods])

# Looking
head(cyto_all_mods$nodeData)
head(cyto_all_mods$edgeData)
dim(cyto_all_mods$nodeData) #nodes
dim(cyto_all_mods$edgeData) #edges
# when 0.05 threshold: 201 nodes, 2251 edges
#####################################################################








###########################################################################################################################
##                                                                                                                       ##
##                                                     Save data                                                         ##
##                                                                                                                       ##
###########################################################################################################################

# Save as RData
save(
  # The 2 matrices for WGCNA analysis
  datExpr, datTraits,

  # Parameters
  sft, minModuleSize, MEDissThres,

  # Final module colours (& when using table() we can see how many genes per module)
  moduleColors,

  # Modoule eigengenes
  MEs,
  
  # Cytoscape dfs
  cyto_all_mods,

  file="data/WGCNA/WGCNA_on_metagenome_essentials__Hel_sft12_mm5.RData")


# ## Loading WGCNA Rdata created in the script "WGCNA_on_metagenome.R"
# # Loading .Rdata saved
# lnames = load(file="data/WGCNA/WGCNA_on_metagenome_essentials__Hel_sft12_mm5.RData")
# lnames



# Saving all data
save.image(file = "data/WGCNA/WGCNA_ALL_FINAL_DATA.RData", compress=F)

# ## Loading WGCNA Rdata created in the script "WGCNA_on_metagenome.R"
# # Loading .Rdata saved
# lnames = load(file="data/WGCNA/WGCNA_ALL_FINAL_DATA.RData")
# lnames


