# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(DESeq2) #for assay()
library(WGCNA) #BiocManager::install("WGCNA", force=T)
library(GO.db) #BiocManager::install("GO.db") 
library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db", force=T) #Human genes

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()


#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html Download the pdf tutorials from here





###########################################################################################################################
##                                                                                                                       ##
##                                            1. Load data & df prep                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# This .Rdata is from the script "Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.R"
# Loading .Rdata saved
lnames <- load(file="data/rld_and_vst.RData")
lnames

## This .Rdata is from the script "Metadata_and_counts_dfs_for_neutrophil_samples.R"
# Loading .Rdata saved
lnames <- load(file="data/metadata_neutrophils_with_VIRGO.RData")
lnames

# This RData is from "Metagenome_VIRGO_basics.R"
# Loading metagenome 
lnames = load(file="data/metagenome/metagenome_abundances_VIRGO.RData")
lnames



########################################### datExpr prep
# Gene expression 
rld_df <- as.data.frame(assay(rld))

# Look
rld_df[1:5,1:5]
head(rld_df)
dim(rld_df) # Genes where all IDs have 0 expression were already filtered out

# Re-naming the columns so the IDs match the metagenome etc 
names(rld_df) <- substring(names(rld_df), 1, 7)
names(rld_df) 
##########################################



########################################## datTraits prep
# Metadata of only quantitative data in metadata_CVF
names(metadata_CVF)
metadata_CVF_quant <- metadata_CVF %>% 
  dplyr::select(IDs_2, # for merge
                Gestation.at.delivery.wks.dec, Age, BMI) #quantitative data 
head(metadata_CVF_quant)


################# Selecting only interesting metagenome as otherwise the heatmap "WGCNA_correlation_heatmap_p0.05.pdf" is too large to be useful
# df of taxa present in 1 sample >=1% RA
VIRGO_t_1RA <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 1) 

# Look
VIRGO_t_1RA
names(VIRGO_t_1RA)
# [1] "Atopobium_vaginae"       "Bifidobacterium_breve"   "Gardnerella_vaginalis"   "Lactobacillus_crispatus" "Lactobacillus_gasseri"  
# [6] "Lactobacillus_iners"     "Lactobacillus_jensenii"  "Lactobacillus_vaginalis" "Prevotella_bivia"

# Select specific taxa we want to keep
selected_taxa <- VIRGO_t %>% 
  dplyr::select(Lactobacillus_acidophilus)  
selected_taxa 

# Bind these together
VIRGO_t_chosen <- cbind(VIRGO_t_1RA, selected_taxa)
VIRGO_t_chosen

# Putting df in alphabetical order https://stackoverflow.com/questions/29873293/dply-order-columns-alphabetically-in-r
VIRGO_t_chosen <- VIRGO_t_chosen %>% 
  dplyr::select(order(colnames(.)))
VIRGO_t_chosen
#################




# Merge prep checks
head(metadata_CVF_quant)
head(VIRGO_t_chosen)
VIRGO_t_chosen$IDs_2 <- rownames(VIRGO_t_chosen)
VIRGO_t_chosen$IDs_2 # missing the ID with no metagenome data
metadata_CVF_quant$IDs_2


# Merge metadata_CVF_quant with VIRGO_t_chosen 
allTraits <- merge(metadata_CVF_quant, VIRGO_t_chosen, by="IDs_2", all.x=T)
#View(allTraits)  success

# Fix Rownames
rownames(allTraits) <- allTraits$IDs_2

# Checks
head(allTraits)
str(allTraits)# All are num apart from IDs_2
##########################################





###########################################################################################################################
##                                                                                                                       ##
##                                WGCNA step 1 - Data input, cleaning and pre-processing                                 ##
##                                                                                                                       ##
###########################################################################################################################

## 1.Data input, cleaning and pre-processing

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# In datExpr the rows should correspond to samples and the columns correspond to genes
datExpr <- as.data.frame(t(rld_df))

# Checks
head(rld_df)
datExpr[1:5, 1:5]
rownames(datExpr)
dim(datExpr) #[1]  10  15982

# 1.b Checking data for excessive missing values and identification of outlier microarray samples
gsg = goodSamplesGenes(datExpr, verbose = 3)
# Flagging genes and samples with too many missing values...
# ..step 1
gsg$allOK #  TRUE 


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

dim(datExpr) #Still  10 15982

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr), method = "average")

# Plot
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_sampleClustering.pdf",  width=6, height=4)
par(cex = 0.6, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# 1.c Loading clinical trait data

# We now read in the trait data and match the samples for which they were measured to the expression samples.
# remove columns that hold information we do not need / choose columns we want #we did this earlier in the prep section

# Looking
head(allTraits)
str(allTraits)
dim(allTraits) 
names(allTraits)

# #Changing Col names
# names(allTraits)
# names(allTraits)[names(allTraits) == "ad_t1.ID"] <- "ID" 
# names(allTraits)


# Form a data frame analogous to expression data that will hold the clinical traits.
Samples_1 = rownames(datExpr)
traitRows = match(Samples_1, allTraits$IDs_2)
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
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_sampleClustering_datTraits.pdf",  width=8, height=6)
par(cex = 0.6)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = NULL,
                    marAll=c(1,10,1,2)) #these margins work #bottom, left, top, right  
dev.off()
# In the above plot white means a low value, red a high value, and grey a missing entry. 





###########################################################################################################################
##                                                                                                                       ##
##                     WGCNA step 2b - Step-by-step network construction and module detection                            ##
##                                                                                                                       ##
###########################################################################################################################

# 2.b Step-by-step network construction and module detection
# 2.b.1 Choosing the soft-thresholding power: analysis of network topology
# Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency [1]. 
# [1] = Bin Zhang, Steve Horvath, 2005, "A general framework for weighted gene co-expression network analysis"
# The authors of [1] have proposed to choose the soft thresholding power based on the criterion of approximate scale-free topology. 
# We refer the reader to that work for more details here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and aids the user in choosing a proper soft-thresholding power. 
# The user chooses a set of candidate powers (the function provides suitable default values), and the function returns a set of network indices that should be inspected

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
powers

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 2, networkType="signed")


# Plot the soft thresholding results:
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_soft_thresholding_power_signed.pdf",  width=9, height=5)
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
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_soft_thresholding_power_signed_Rsquared.pdf", width=3.5, height=3.5)
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
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_soft_thresholding_power_signed_mean_connectivity.pdf", width=3.5, height=3.5)
par(mfrow = c(1,1), #1x1 panel
    cex = 0.8, mar = c(5,5,1,1)) #bottom, left, top, right 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft-thresholding power",ylab="Mean connectivity", type="n",
     main = NULL)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red") #cex=cex1 is making it error, so have removed it from this line...
dev.off()



# In the tutorial they choose the power 6, which was the lowest power for which the scale-free topology fit index reaches 0.90.
# "If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers..." https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html#:~:text=If%20the,topology%20approximation
# Ours does not reach 0.8....
# "...0.8 for reasonable powers (less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) and the mean connectivity remains relatively high (in the hundreds or above), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest."

#Checking my data meets the conditions mentioned above
sft$fitIndices[c(1,2,5)] # Meets conditions mentioned above. Therefore... "data exhibit a strong driver that makes a subset of the samples globally different from the rest"

# "
# If the lack of scale-free topology fit turns out to be caused by an interesting biological variable that one does not want to remove (i.e., adjust the data for), the appropriate soft-thresholding power can be chosen based on the number of samples as in the table below. This table has been updated in December 2017 to make the resulting networks conservative.
# Number of samples             Unsigned and signed hybrid networks             Signed networks
#      Less than 20	                                              9	                         18
# "


# 2.b.2 Co-expression similarity and adjacency
# We now calculate the adjacencies, using the selected soft thresholding power:
softPower = 18 # based on WGCNA Q&A info C&P above
adjacency = adjacency(datExpr, power = softPower, type = "signed")

# # Save adjacency to use on Rosalind
# save(adjacency,
#      file="data/WGCNA/WGCNA_adjacency_with_VIRGO_in_datTraits.RData")

# # Perform a checksum to make sure I the adjacency in this script & the .RData match
# library(matter) #BiocManager::install("matter")
# checksum(matter(adjacency), algo = c("md5"))
# #"bec7799bce03eaf6492eb86e1689f888" # with set.seed(100) 18/06/21


# 2.b.3 Topological Overlap Matrix (TOM)
# To minimize effects of noise and spurious associat, ions, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
# Turn adjacency into topological overlap
Sys.time()
TOM = TOMsimilarity(adjacency) 
dissTOM = 1-TOM
Sys.time() # it takes about 25min to run



# # Loading .Rdata saved 
# lnames = load(file="data/WGCNA/WGCNA_dissTOM.RData") 
# lnames


# 2.b.4 Clustering using TOM
# We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. 
# Note that we use the function hclust that provides a much  faster hierarchical clustering routine than the standard hclust function.

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_Gene_clustering.pdf",  width=8, height=5)
par(mfrow = c(1,1), #1x1 panels
    cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



# Set the minimum module size:
minModuleSize = 30 

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# ..cutHeight not given, setting it to 0.987  ===>  99% of the (truncated) height range in dendro.
# ..done.

table(dynamicMods)
#The function returned 28 modules labeled 1–28 largest to smallest. Label 0 is reserved for unassigned genes. 


#We now plot the module assignment under the gene dendrogram:

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) #These are the unmerged modules
dim(table(dynamicColors)) #55 modules originally

# Plot the dendrogram and colors underneath
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_Gene_clustering__unmerged_modules.pdf",  width=8, height=5)
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
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_ME_clustering.pdf",  width=8, height=6)
par(cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25

# Plot the result
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_ME_clustering_with_diss_threshold.pdf",  width=8, height=6)
par(cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()



# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors
table(mergedColors) # which colours
dim(table(mergedColors)) # number of colours - 45

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath 
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_Gene_clustering__unmerged_and_merged_modules.pdf",  width=12, height=9)
par(cex = 0.9, #scaling text size
    mar = c(1,5,4,0)) #bottom, left, top, right
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="Gene dendrogram comparing merged and unmerged modules with cut height of 0.25")
dev.off()


# Plot the dendrogram and colors underneath
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_Gene_clustering__final_modules.pdf",  width=8, height=4.5)
par(cex = 0.9) 
plotDendroAndColors(geneTree, mergedColors, "Final module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = NULL,
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
# bisque4           black            blue           brown          brown4            cyan       darkgreen        darkgrey     darkmagenta 
#      69             662             870             954             113             418             232             283             303 
# darkolivegreen      darkorange     darkorange2         darkred   darkturquoise     floralwhite           green     greenyellow            grey 
#            134             192             182             242             273              80             723             484               6 
# grey60           ivory       lightcyan      lightcyan1 lightsteelblue1     lightyellow         magenta   mediumpurple3    midnightblue 
#    318            1034             321             272              93             247             582              97             323 
# navajowhite2          orange  palevioletred3            pink           plum2          purple       royalblue          salmon         sienna3 
#           35             200              36             626              55             552             246             448             132 
# skyblue        skyblue3       steelblue             tan        thistle1       turquoise          violet           white     yellowgreen 
#     183             109             170             482              44            1972             152             923             110
dim(table(moduleColors)) # 45



# Calculate dissimilarity of module eigengenes
MEDiss_final = 1-WGCNA::cor(MEs)

# Cluster module eigengenes
METree_final = hclust(as.dist(MEDiss_final), method = "average")

# Plot the result
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_ME_clustering_finalMEs.pdf",  width=8, height=5)
par(cex = 0.9,
    mar = c(0.5,5,1,0)) #bottom, left, top, right
plot(METree_final, main = NULL,
     xlab = "", sub = "")
dev.off()






###########################################################################################################################
##                                                                                                                       ##
##                          WGCNA step 3 - Relating modules to external clinical trait                                   ##
##                                                                                                                       ##
###########################################################################################################################


#3 Relating modules to external clinical traits 3.a Quantifying module–trait associations

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Excluding grey from heatmap
MEs_noGrey <- MEs
MEs_noGrey$MEgrey <- NULL

# Calculating ME correlations and p-values with datTraits
moduleTraitCor = WGCNA::cor(MEs_noGrey, datTraits, use = "p") #the correlation values used in the heatmap
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) # the p-values used in the heatmap

#Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot - no grey
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_WGCNA_heatmap.pdf",  width=6, height=9)
par(cex = 0.6, #scaling text size
    mar = c(9, 10, 3, 2)) #bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs_noGrey),
               ySymbols = names(MEs_noGrey),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



############################################################### Altered version of heatmap - where only sig values are shown

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


# Making nice axis labels
x_axis_labels1 <- names(as.data.frame(moduleTraitCor_0.05))
x_axis_labels1 <- gsub("_"," ", x_axis_labels1)
x_axis_labels1 <- gsub("."," ", x_axis_labels1, fixed = TRUE) #without fixed argument it replaces everything lol
x_axis_labels1 <- gsub(" wks dec","", x_axis_labels1, fixed = TRUE)
x_axis_labels1 


# Display the correlation values within a heatmap plot - with p-values
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_WGCNA_heatmap_p0.05.pdf",  width=8, height=9)
par(cex = 0.6, #scaling text size
    mar = c(8, 9, 1, 1)) #bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = x_axis_labels1,
               yLabels = names(MEs_noGrey),
               ySymbols = names(MEs_noGrey),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_0.05, 
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = NULL)
dev.off()
############################################################### 


############################################################### Altered version of heatmap - where only MEs we're interested in are shown
################# 1. working with the correlation matrix input for colours etc
# Looking at matrix used to plot heatmap above
moduleTraitCor
str(moduleTraitCor)
is.matrix(moduleTraitCor) # matrix type

# Make matrix into df
moduleTraitCor_df <- as.data.frame(moduleTraitCor)
str(moduleTraitCor_df)

# Make ME colours into a row
moduleTraitCor_df$MEs <- rownames(moduleTraitCor_df)
moduleTraitCor_df$MEs

# df of only MEs we are interested in https://www.datanovia.com/en/lessons/subset-data-frame-rows-in-r/
moduleTraitCor_selectedMEs_df <- moduleTraitCor_df %>% 
  filter(MEs %in% c("MEsienna3",       "MElightcyan",    "MEviolet",     "MEmagenta",      "MEdarkgreen", 
                    "MEplum2",        "MEdarkred",    "MEdarkmagenta",  "MEwhite", # MEmidnightblue removed because it is no longer sig with any metadata (as we have removed immaturity score)
                    "MEdarkorange",    "MEgreenyellow",  "MEsteelblue",  "MEgreen",        "MEskyblue3", 
                    "MEivory",         "MEbrown4",       "MEskyblue",    "MEturquoise",    "MEdarkturquoise")) # 20 colours

# Checks
dim(moduleTraitCor_selectedMEs_df) #19 14
moduleTraitCor_selectedMEs_df # matches

# Get rid of ME coloumn
moduleTraitCor_selectedMEs_df$MEs <- NULL 

# Look
names(moduleTraitCor_selectedMEs_df)
rownames(moduleTraitCor_selectedMEs_df)

# Make it into a matrix
moduleTraitCor_selectedMEs <- as.matrix(moduleTraitCor_selectedMEs_df)
str(moduleTraitCor_selectedMEs)
#################


################# 2. Working with text overlay for plot - correlations
# Now looking at the text we want to print on interesting modules
moduleTraitCor_0.05[,1:5] # we made this df for the heatmap above
str(moduleTraitCor_0.05)

# Make matrix into df
moduleTraitCor_0.05_df <- base::as.data.frame(moduleTraitCor_0.05)
str(moduleTraitCor_0.05_df)

# Make ME colours into a row
moduleTraitCor_0.05_df$MEs <- rownames(moduleTraitCor_0.05_df)
moduleTraitCor_0.05_df$MEs

# df of only MEs we are interested in https://www.datanovia.com/en/lessons/subset-data-frame-rows-in-r/
moduleTraitCor_0.05_selectedMEs_df <- moduleTraitCor_0.05_df %>% 
  filter(MEs %in% c("MEsienna3",       "MElightcyan",    "MEviolet",     "MEmagenta",      "MEdarkgreen", 
                    "MEplum2",        "MEdarkred",    "MEdarkmagenta",  "MEwhite", # MEmidnightblue removed because it is no longer sig with any metadata (as we have removed immaturity score)
                    "MEdarkorange",    "MEgreenyellow",  "MEsteelblue",  "MEgreen",        "MEskyblue3", 
                    "MEivory",         "MEbrown4",       "MEskyblue",    "MEturquoise",    "MEdarkturquoise")) # 20 colours

# Checks
dim(moduleTraitCor_0.05_selectedMEs_df) #19 14
moduleTraitCor_0.05_selectedMEs_df # Looks good 

# Get rid of ME coloumn
moduleTraitCor_0.05_selectedMEs_df$MEs <- NULL 

# Look
names(moduleTraitCor_0.05_selectedMEs_df)
rownames(moduleTraitCor_0.05_selectedMEs_df)

# Make it into a matrix
moduleTraitCor_0.05_selectedMEs <- as.matrix(moduleTraitCor_0.05_selectedMEs_df)
str(moduleTraitCor_0.05_selectedMEs)

# Look
moduleTraitCor_0.05_selectedMEs 
#################


################# 3. Working with text overlay for plot - p-values
# Now looking at the text we want to print on interesting modules
moduleTraitPvalue_0.05[,1:5] # we made this df for the heatmap above
str(moduleTraitPvalue_0.05)

# Make matrix into df
moduleTraitPvalue_0.05_df <- base::as.data.frame(moduleTraitPvalue_0.05)
str(moduleTraitPvalue_0.05_df)

# Make ME colours into a row
moduleTraitPvalue_0.05_df$MEs <- rownames(moduleTraitPvalue_0.05_df)
moduleTraitPvalue_0.05_df$MEs

# df of only MEs we are interested in 
moduleTraitPvalue_0.05_selectedMEs_df <- moduleTraitPvalue_0.05_df %>% 
  filter(MEs %in% c("MEsienna3",       "MElightcyan",    "MEviolet",     "MEmagenta",      "MEdarkgreen", 
                    "MEplum2",        "MEdarkred",    "MEdarkmagenta",  "MEwhite", # MEmidnightblue removed because it is no longer sig with any metadata (as we have removed immaturity score)
                    "MEdarkorange",    "MEgreenyellow",  "MEsteelblue",  "MEgreen",        "MEskyblue3", 
                    "MEivory",         "MEbrown4",       "MEskyblue",    "MEturquoise",    "MEdarkturquoise")) # 20 colours

# Checks
dim(moduleTraitPvalue_0.05_selectedMEs_df) #19 14
moduleTraitPvalue_0.05_selectedMEs_df # Looks good 

# Get rid of ME coloumn
moduleTraitPvalue_0.05_selectedMEs_df$MEs <- NULL 

# Look
names(moduleTraitPvalue_0.05_selectedMEs_df)
rownames(moduleTraitPvalue_0.05_selectedMEs_df)

# Make it into a matrix
moduleTraitPvalue_0.05_selectedMEs <- as.matrix(moduleTraitPvalue_0.05_selectedMEs_df)
str(moduleTraitPvalue_0.05_selectedMEs)

# Look
moduleTraitPvalue_0.05_selectedMEs 
################# 


################# 4. Merging cor & p for final text overlay 
# making the p-values in brackets on a new line
p4heatmap_selectedMEs <- moduleTraitPvalue_0.05_selectedMEs
p4heatmap_selectedMEs <- base::ifelse(p4heatmap_selectedMEs==" ",p4heatmap_selectedMEs,paste("\n(", p4heatmap_selectedMEs, ")", sep = ""))
p4heatmap_selectedMEs[1:5,] # 
moduleTraitPvalue_0.05_selectedMEs[1:5,] # matching

#Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix_0.05_selectedMEs = paste(moduleTraitCor_0.05_selectedMEs, p4heatmap_selectedMEs, sep = "")
textMatrix_0.05_selectedMEs
dim(textMatrix_0.05_selectedMEs) = dim(moduleTraitCor_0.05_selectedMEs)
#################

################# 5. Making nice axis labels
# Save as vector
x_axis_labels <- names(moduleTraitCor_selectedMEs_df)
x_axis_labels

# Making x axis labels without the underscores and dots
x_axis_labels <- gsub("_"," ", x_axis_labels)
x_axis_labels <- gsub("."," ", x_axis_labels, fixed = TRUE) #without fixed argument it replaces everything lol
x_axis_labels <- gsub(" wks dec","", x_axis_labels, fixed = TRUE)
x_axis_labels 
#################

################# 6. Plotting final heatmap
# Display the correlation values within a heatmap plot of only modules we're interested in - with p-values
# pdf:
pdf("plots/WGCNA/VIRGO_in_datTraits/rna_WGCNA_heatmap_p0.05_selectedMEs.pdf",  width=8, height=8)
par(cex = 0.8, #scaling text size
    mar = c(8, 9, 1, 1)) #bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor_selectedMEs,
               xLabels = x_axis_labels,
               yLabels = rownames(moduleTraitCor_selectedMEs_df),
               ySymbols = rownames(moduleTraitCor_selectedMEs_df),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_0.05_selectedMEs, 
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = NULL)
dev.off()
#################
############################################################### 





# 3.d Summary output of network analysis results
# We have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure. 
# We now merge this statistical information with gene annotation and write out a file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS Excel or Open Office Calc. 
# Our expression data are only annotated by probe ID [symbols in my case] names: the command
names(datExpr)
#will return all probe IDs included in the analysis. Similarly,
names(datExpr)[moduleColors=="lightcyan"]
#will return probe IDs belonging to the lightcyan module. 





# We now create a data frame holding the following information for all probes: gene symbol, Locus Link ID (Entrez code), module color, gene significance for Succinate, and module membership and p-values in all modules. 
# The modules will be ordered by their significance for Lactobacillus_crispatus, with the most significant ones to the left.

# Create the starting data frame
Geneinfo0 = data.frame(Genes = names(datExpr),
                       moduleColor = moduleColors #,
                       #geneTraitSignificance, # !need to run an above scatter with correct trait [Lactobacillus_crispatus for now] to look at so this geneTraitSignificance is correct
                       #GSPvalue
)
Geneinfo0[1:10, ]



# # Order modules by their significance for Lactobacillus_crispatus
# modOrder = order(-abs(WGCNA::cor(MEs, Lactobacillus_crispatus, use = "p")))
# 
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(Geneinfo0)
#   Geneinfo0 = data.frame(Geneinfo0, geneModuleMembership[, modOrder[mod]],
#                         MMPvalue[, modOrder[mod]])
#   names(Geneinfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }

############ Commented out reordering as it ruined my previous GO results because it made the allLLIDs (derived from Geneinfo) no longer be in the same order as datExpr
# # Order the genes in the Geneinfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(Geneinfo0$moduleColor, -abs(Geneinfo0$GS.Lactobacillus_crispatus))
Geneinfo = Geneinfo0#[geneOrder, ]

# Add rownames
rownames(Geneinfo) <- Geneinfo$Genes

# Adding the entrez ID 
Geneinfo$entrez <-  mapIds(org.Hs.eg.db,
                           keys=row.names(Geneinfo), 
                           column="ENTREZID", #the new column of the gene type ID we want
                           keytype="SYMBOL", #the ID in res already
                           multiVals="first")

# Checking
head(Geneinfo$entrez, 15)
# [1] "1"      "144568" "8086"   "79719"  "22848"  "28971"  "14"     "15"     "25980"  "16"     "132949" "60496"  "10157"  "26574"  "9625" 

# Reorder df 
Geneinfo <- Geneinfo %>% 
  dplyr::select(Genes, entrez, everything())

# Checking if the above code just reorders the columns
head(Geneinfo$entrez, 15)
# [1] "1"      "144568" "8086"   "79719"  "22848"  "28971"  "14"     "15"     "25980"  "16"     "132949" "60496"  "10157"  "26574"  "9625" 


# This data frame can be written into a text-format spreadsheet
write_csv(Geneinfo, "data/WGCNA/Geneinfo_Lactobacillus_crispatus_VIRGO.csv")









###########################################################################################################################
##                                                                                                                       ##
##     WGCNA step 4 - Interfacing network analysis with other data such as functional annotation and gene ontology       ##
##                                                                                                                       ##
###########################################################################################################################


# 4 Interfacing network analysis with other data such as functional annotation and gene ontology
# Our previous analysis has identified several modules that are highly associated with the trait of interest. 
# To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, whether they are significantly enriched in certain functional categories etc.



# 4.a Output gene lists for use with online software and services
# One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology and functional enrichment analysis suites such as David or AmiGO. 
# For example, we write out the LocusLinkID (entrez) codes for a specific module into a file:

## This bit is a substitute for reading in the gene annotation, as we already have gene labels (they had probe IDs)

# Checks
datExpr[1:10,1:10]
Geneinfo[1:10,]
dim(datExpr) # 9 15982       
dim(Geneinfo) #15982    3   # Cool - haven't dropped any genes

# Checking genes are in the correct order
datExpr[1:5,1:5]
Geneinfo[1:5,]
rownames(Geneinfo) == names(datExpr)
table(rownames(Geneinfo)== names(datExpr), useNA="always") # all matching


#  Locus Link IDs (i.e. entrez)
allLLIDs <- Geneinfo$entrez
length(allLLIDs) #15982
table(is.na(Geneinfo$entrez)) #295 genes without an entrez ID
head(allLLIDs, 15)
# [1] "1"      "144568" "8086"   "79719"  "22848"  "28971"  "14"     "15"     "25980"  "16"     "132949" "60496"  "10157"  "26574"  "9625"

# As background in the enrichment analysis, we will use all Genes in the analysis.
write.table(as.data.frame(allLLIDs), file = "data/WGCNA/gene_lists/entrez_IDs_ALL.txt",
            row.names = FALSE, col.names = FALSE)

# Checking
dim(as.data.frame(allLLIDs)) #[1] 15982     1



# 4.b Enrichment analysis directly within R
# The WGCNA package now contains a function to perform GO enrichment analysis using a simple, single step. 
# To run the function, Biconductor packages GO.db, AnnotationDBI, and the appropriate organism-specific annotation package(s) need to be installed before running this code. 
# The organism-specific packages have names of the form org.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc. 
# Please visit the Bioconductor main page at http://www.bioconductor.org to download and install the required packages. 
# (I have this installation step at beginning of the script)

# Calling the GO enrichment analysis function GOenrichmentAnalysis is very simple. 
# The function takes a vector of module labels, and the Entrez (a.k.a. Locus Link) codes for the genes whose labels are given.
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10) 
#NOTE: GOenrichmentAnalysis is deprecated. Please use function enrichmentAnalysis from R package anRichment, available from https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/
rownames(Geneinfo)
rownames(moduleColors)


# The function runs for awhile and returns a long list, the most interesting component of which is
tab = GOenr$bestPTerms[[4]]$enrichment # This is BP, MF & CC

# This is an enrichment table containing the 10 best terms for each module present in moduleColors. Names of the columns within the table can be accessed by
names(tab)

# We refer the reader to the help page of the function within R (available using ?GOenrichmentAnalysis at the R prompt) for details of what each column means. 
# Because the term definitions can be quite long, the table is a bit difficult to display on the screen,so we save as a csv:
write.table(tab, file = "data/WGCNA/GO/GOenrichmentAnalysis_top10_GO_in_all_modules_VIRGO.csv", sep = ",", quote = TRUE, row.names = FALSE)


##### Sig GO
summary(tab)
# all GO are sig using their nominal p-values
# no GO are sig using their Bonferoni p-values
##### 




# On the other hand, to quickly take a look at the results, one can also abridge the table a bit and display it directly on screen:
screenTab = tab[, c(1, 2, 5, 6, 7, 11, 12, 13)]
names(screenTab) # all looking fine

#  Round the numeric columns to 2 decimal places:
screenTab[, c(3, 4)] = signif(apply(screenTab[, c(3, 4)], 2, as.numeric), 2)

#  Truncate the the term name to at most 60 characters
screenTab[, 8] = substring(screenTab[, 8], 1, 60)

#  Shorten the column names:
colnames(screenTab) = c("module", "size", "p_val", "Bonf", "nInTerm", "GO_ID", "ont", "term_name")
rownames(screenTab) = NULL

# Finally, display the enrichment table:
screenTab



############################################### GO Subsets 
##### modules correlated with L. crispatus
screenTab %>% 
  filter(module=="lightcyan") 

##### modules correlated with Lactobacillus_iners
screenTab %>% 
  filter(module=="darkturquoise")
screenTab %>% 
  filter(module=="magenta")
screenTab %>% 
  filter(module=="violet")

##### modules correlated with anaerobic microbiota etc
screenTab %>% 
  filter(module=="white") 
screenTab %>% 
  filter(module=="darkmagenta")






###############################################  Better way to get GO
# The help page says: "NOTE: GOenrichmentAnalysis is deprecated" so trying their suggested replacement

# Installing
# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R")
# installAnRichment()
library(anRichment)

# The reference collection has to be created before calling enrichmentAnalysis. 
# Several reference collections are available internally. 
# Here we use the GO collection, accessible by calling the function buildGOcollection.
GOcollection <-  buildGOcollection(organism = "human")

# The next call evaluates the enrichment of the gene modules in the collection of GO terms.
?anRichmentMethods::enrichmentAnalysis
GOenrichment <- enrichmentAnalysis(
  classLabels = moduleColors, identifiers = allLLIDs,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 0.05,
  getFDR=T,
  thresholdType = "nominal", # "Bonferroni", "FDR", "nominal"
  getOverlapEntrez = F,
  getOverlapSymbols = T,
  ignoreLabels = "grey")
collectGarbage()

# The returned object is a list with several components:
names(GOenrichment)

# The enrichment results are summarized in the component enrichmentTable:
names(GOenrichment$enrichmentTable)

#The last column contains concatenated lists of overlap genes and would be diffcult to display 
# hence we shorten the last column for display purposes and display the first few rows:
table.display = GOenrichment$enrichmentTable
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 20,
                                            split = "|")
dim(table.display) # 7863   18 
summary(table.display)
head(table.display)
names(table.display)

# To check the GO in table.display is accurate let's make sure that APC2|ARFGEF2|KIF2C are in the bisque4 module
bisque4_module_genes <- names(datExpr)[moduleColors=="bisque4"] 
grep("APC2", bisque4_module_genes) #27
bisque4_module_genes[27] #"APC2"
grep("KIF2C", bisque4_module_genes) #1
bisque4_module_genes[1] #"KIF2C" # great 

# Export the full enrichment table 
write.csv(GOenrichment$enrichmentTable, file = "data/WGCNA/GO/enrichmentAnalysis_VIRGO.csv",
          row.names = FALSE)

# I look at the modules which have bonferroni sig p-values in the script WGCNA_plotting_selected_GO_VIRGO.R
############################################### 







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
  
  # Df with gene symbols, entrez IDs, L. crispatus & module membership info
  Geneinfo,
  
  # GO - 1st method
  GOenr, tab, screenTab, #tab & screenTab are reduced versions of GOenr  
  
  # GO - 2nd method
  GOenrichment, table.display, #table.display is a reduced version of GOenrichment 
  
  file="data/WGCNA/WGCNA_essentials_rld_with_VIRGO_in_datTraits.RData")

# ## This .Rdata is from the script "WGCNA_pipeline__neutrophils_VIRGO.R"
# # Loading .Rdata saved
# lnames = load(file="data/WGCNA/WGCNA_essentials_rld_with_VIRGO_in_datTraits.RData")
# lnames

# Saving all data
save.image(file = "data/WGCNA/ALL_WGCNA_objects.RData", compress=T)

