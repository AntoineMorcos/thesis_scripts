# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats)
library(DESeq2) #for assay()
library(vegan) #for Hellinger
library(corrplot)

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()
#list.files()






###########################################################################################################################
##                                                                                                                       ##
##                                             Load data & df prep                                                       ##
##                                                                                                                       ##
###########################################################################################################################

## These csvs were made by Flavia (I downloaded copies on 3rd Dec 21 and did not change the file names)
# Load metagenome counts data
MG_RA_t <- read.table("data/from_Flavia/Virgo.Percentage_NoBurlholderia.txt", sep="\t", header=T) #relative abundance here

# Load metadata
metadata_df <- read.csv("data/from_Flavia/Metadata_Version_December2021_TOPSPECIES.csv")


# Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
lnames = load(file="data/immune_CytokinesComplementsIgM_2023.RData")
lnames

## Loading metabolite .Rdata created in the script "Metabolites_exploration.R"
# Loading .Rdata saved
lnames = load(file="data/metabolites_df.RData")
lnames



###################################################### MG_RA_t prep
# Looking
dim(MG_RA_t) #89 231
MG_RA_t[1:5,1:5]
#summary(MG_RA_t)

# Change colname
MG_RA_t <- MG_RA_t %>% 
  dplyr::rename(Sample=PC)
MG_RA_t[1:5,1:5]

# Making ID1 from the Sample col
MG_RA_t$ID1 <- MG_RA_t$Sample
MG_RA_t$ID1 <- substring(MG_RA_t$ID1, 3) #keep characters after the 3rd character
MG_RA_t$ID1 <- gsub("_","", MG_RA_t$ID1) # getting rid of "_" in the AC women
MG_RA_t$ID1 <- gsub("^","ID_", MG_RA_t$ID1) # put ID_ at the beginning so it doesn't go mad when it is a col name
MG_RA_t$ID1
MG_RA_t %>% dplyr::select(ID1,Sample)
MG_RA_t$Sample <- NULL

# Reorder df by ID
MG_RA_t <- MG_RA_t %>% 
  dplyr::arrange(ID1)
MG_RA_t[1:5,1:5]
MG_RA_t$ID1

# Put ID1 col as the first col
MG_RA_t <- MG_RA_t %>% 
  dplyr::select(ID1, everything())
MG_RA_t[1:5,1:5]

# Remove IDs that should be excluded
dim(MG_RA_t) #89 231
MG_RA_t <- MG_RA_t %>% 
  dplyr::filter(!ID1 %in% c(IDs_to_exclude))
dim(MG_RA_t) #87 231, so 2 IDs dropped, which makes sense based on the metadata FF sent me

# Set rownames for transpose & delete ID1 col
rownames(MG_RA_t) <- MG_RA_t$ID1
MG_RA_t[1:5,1:5]
MG_RA_t$ID1 <- NULL
MG_RA_t[1:5,1:5]

# Transpose
MG_RA <- as.data.frame(t(MG_RA_t))
MG_RA[1:5,1:5]


# Looking
colSums(MG_RA_t)
rowSums(MG_RA_t) # all just under 100% 
######################################################




###################################################### MG_Hel_t prep
# Perform the Hellinger transformation
?decostand
MG_Hel_t <- decostand(MG_RA_t, method="hellinger")
MG_Hel_t[1:5,1:5]

# Transpose
MG_Hel <- as.data.frame(t(MG_Hel_t))
MG_Hel[1:5,1:5]

# Looking
colSums(MG_Hel_t)
rowSums(MG_Hel_t) 
######################################################




###########################################################################################################################
##                                                                                                                       ##
##                                               Exporting dfs                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Save as RData
save(
  # RA - where 100 is 100% (i.e. 1 in the matrix represents 1% not 100%)
  MG_RA_t, MG_RA,
  
  # Hellinger transformed data (counts -> RA -> Hellinger)
  MG_Hel_t, MG_Hel,
  
  file="data/metagenome_counts_and_RA.RData")


# # Loading metagenome Rdata created in the script "Metagenome_exploration.R"
# lnames = load(file="data/metagenome_counts_and_RA.RData")
# lnames





###########################################################################################################################
##                                                                                                                       ##
##                                           Looking at RA in general                                                    ##
##                                                                                                                       ##
###########################################################################################################################

# Looking at taxa of interest
summary(MG_RA_t$Lactobacillus_crispatus)
summary(MG_RA_t$Lactobacillus_iners)
summary(MG_RA_t$Lactobacillus_acidophilus)
summary(MG_RA_t$Gardnerella_vaginalis)







###########################################################################################################################
##                                                                                                                       ##
##                                     Looking at other specific taxa of interest                                        ##
##                                                                                                                       ##
###########################################################################################################################

# E. coli or Finegoldia_magna
MG_RA_t %>% dplyr::filter(Escherichia_coli>0 | Finegoldia_magna>0) %>% dplyr::select(Escherichia_coli, Finegoldia_magna)


# The species which were excluded in WGCNA (WGCNA_on_metagenome_in_Black_women_Hel.R)
# #Removing genes: Eubacterium_saburreum, Granulicatella_elegans, Lactobacillus_ruminis
MG_RA_t %>% dplyr::filter(Eubacterium_saburreum>0 | Granulicatella_elegans>0 | Lactobacillus_ruminis>0) %>% 
  dplyr::select(Eubacterium_saburreum, Granulicatella_elegans, Lactobacillus_ruminis) # no one has any >0
summary(MG_RA_t$Eubacterium_saburreum)
summary(MG_RA_t$Granulicatella_elegans)
summary(MG_RA_t$Lactobacillus_ruminis) # All IDs have 0 RA


# Looking at STIs 14/04/2023
MG_RA_t %>% filter(Chlamydia_trachomatis > 0) %>% select(Chlamydia_trachomatis) #1 ID 
MG_RA_t %>% filter(Neisseria_gonorrhoeae > 0) %>% select(Neisseria_gonorrhoeae) #1 ID 
MG_RA_t %>% filter(Mycoplasma_genitalium > 0) %>% select(Mycoplasma_genitalium) #1 ID








###########################################################################################################################
##                                                                                                                       ##
##                                           Correlation heatmap on RA                                                   ##
##                                                                                                                       ##
###########################################################################################################################

################################ INFO:
# Corrplot help from: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

# Looking at dim
dim(MG_RA_t) #87 230 # FIRST COLUMN IS ID HERE


# I'm thresholding the taxa to make the heatmap more readable. Also the corrplot won't generate if I try to include all taxa (just gives an error).
################################


################################ Looking at number of taxa at different minimum RA thresholds
# df with threshold of 1%
MG_RA_t_1 <- MG_RA_t %>% 
  select_if(~max(., na.rm = TRUE) >= 1) 
dim(MG_RA_t_1) #87 38 
names(MG_RA_t_1)

# df with threshold of 5%
MG_RA_t_5 <- MG_RA_t %>% 
  select_if(~max(., na.rm = TRUE) >= 5)
dim(MG_RA_t_5) #87 21 
names(MG_RA_t_5)
################################


################################ Taxa we want to keep regardless of whether it meets the threshold
# Select specific taxa we want to keep
selected_taxa <- MG_RA_t %>% 
  dplyr::select(Lactobacillus_acidophilus)  
selected_taxa 
################################ 


################################ Getting info on corrplot
# Functions we'll use
?corrplot
??corr.test

# Referencing
citation('corrplot') 
################################

dev.off()


################################################################### Heatmap with only certain taxa (taxa>=5% plus Lactobacillus_acidophilus)

######################## df prep
# Bind taxa with >=5% RA with selected taxa (i.e. Lactobacillus_acidophilus)
MG_RA_t_chosen5 <- merge(MG_RA_t_5, selected_taxa, by="row.names")
head(MG_RA_t_chosen5)
rownames(MG_RA_t_chosen5) <- MG_RA_t_chosen5$Row.names
head(MG_RA_t_chosen5)

# Putting df in alphabetical order 
MG_RA_t_chosen5 <- MG_RA_t_chosen5 %>% 
  dplyr::select(order(colnames(.))) 
MG_RA_t_chosen5
names(MG_RA_t_chosen5)

# Removing Sample/ID column
MG_RA_t_chosen5$Row.names <- NULL

# This will make a much more manageable heatmap now 
dim(MG_RA_t) #87 230 
dim(MG_RA_t_chosen5) #87 22
########################




########################  corrplot on selected taxa
# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_5 <- psych::corr.test(MG_RA_t_chosen5, adjust = "fdr")
corr_5 <- ct_5$r
p.mat_5 <- ct_5$p

# Finding n for captions
ct_5$n #87 all

####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_5 <- p.mat_5

# Make nominal p-values (lower.tri) into NAs
p_adj_5[lower.tri(p_adj_5)] <- NA # don't need this step but I include to verify
p_adj_5 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_5[lower.tri(p_adj_5)] = t(p_adj_5)[lower.tri(p_adj_5)]
p_adj_5 

# Checking symmetry 
isSymmetric(p_adj_5) #TRUE
#######

# Plot corr_5 with hclust
png("plots/metagenome_only/heatmap_metagenome_selectedTaxa5perc_allPs.png", res=600, units="in", width=7, height=7)
corrplot(corr_5, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



## Plot corr_5 with hclust (leave blank on non-significant coefficient & add significant correlation coefficients in text)
pdf("plots/metagenome_only/heatmap_metagenome_selectedTaxa5perc_FDRsig.pdf", width=7, height=7)
corrplot(corr_5, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_5, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.8, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.5, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
######################## 

################################################################### 





################################################################### Heatmap with only certain taxa (taxa>=1% plus Lactobacillus_acidophilus)

# Bind taxa with >=1% RA with selected taxa (i.e. Lactobacillus_acidophilus)
MG_RA_t_chosen1 <- merge(MG_RA_t_1, selected_taxa, by="row.names")
head(MG_RA_t_chosen1)
rownames(MG_RA_t_chosen1) <- MG_RA_t_chosen1$Row.names
head(MG_RA_t_chosen1)

# Putting df in alphabetical order 
MG_RA_t_chosen1 <- MG_RA_t_chosen1 %>% 
  dplyr::select(order(colnames(.))) 
MG_RA_t_chosen1
names(MG_RA_t_chosen1)

# Removing Sample/ID column
MG_RA_t_chosen1$Row.names <- NULL

# This will make a much more manageable heatmap now 
dim(MG_RA_t) #87 230 
dim(MG_RA_t_chosen1) #87 39 
########################




########################  corrplot on selected taxa
# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_1 <- psych::corr.test(MG_RA_t_chosen1, adjust = "fdr")
corr_1 <- ct_1$r
p.mat_1 <- ct_1$p

####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_1 <- p.mat_1

# Make nominal p-values (lower.tri) into NAs
p_adj_1[lower.tri(p_adj_1)] <- NA # don't need this step but I include to verify
p_adj_1 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_1[lower.tri(p_adj_1)] = t(p_adj_1)[lower.tri(p_adj_1)]
p_adj_1 

# Checking symmetry 
isSymmetric(p_adj_1) #TRUE
#######

# Plot corr_1 with hclust
png("plots/metagenome_only/heatmap_metagenome_selected_taxa_1perc_corrplot_pretty.png", res=600, units="in", width=12, height=12)
corrplot(corr_1, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 0.7, tl.col = 'black',
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



## Plot corr_1 with hclust (leave blank on non-significant coefficient & add significant correlation coefficients in text)
png("plots/metagenome_only/heatmap_metagenome_selected_taxa_1perc_corrplot_sig.png", res=600, units="in", width=12, height=12)
corrplot(corr_1, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_1, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.7, tl.col = 'black',
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.4, number.digits = 2) 
dev.off()
######################## 


###################################################################







###########################################################################################################################
##                                                                                                                       ##
##                                                  Immune vs taxa                                                       ##
##                                                                                                                       ##
###########################################################################################################################

################################################################### Prep & All together heatmap (same x&y axes)
# Looking
MG_RA_t_chosen5[1:10,1:7]
immu_small[1:10,]

# Merge prep - Making Participant.ID column
MG_pre_merge_immu <- MG_RA_t_chosen5
MG_pre_merge_immu$Participant.ID <- substring(rownames(MG_RA_t_chosen5), 4) #keep characters after the 4th character
MG_pre_merge_immu$Participant.ID
MG_pre_merge_immu[1:5, 21:23]
dim(immu_small) #81 14
dim(MG_pre_merge_immu) #87 23

# Merge 
?base::merge
immu_some_metag <- base::merge(immu_small, MG_pre_merge_immu, by="Participant.ID", 
                                all.x=F, all.y=T) # only keep immune women where we have metag (x), & keep all women with metag (y) 

# Tidy & check
dim(immu_some_metag) #87 36
row.names(immu_some_metag) <- immu_some_metag$Participant.ID
immu_some_metag$Participant.ID <- NULL

head(immu_some_metag)
str(immu_some_metag)


# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_immu_metag5 <- psych::corr.test(immu_some_metag, adjust = "fdr")
corr_immu_metag5 <- ct_immu_metag5$r
p.mat_immu_metag5 <- ct_immu_metag5$p #matrix isn't symmetrical which creates issues in the rectangle asymmetrical plots
# Probability values (Entries above the diagonal are adjusted for multiple tests.) 


####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_immu_metag5 <- p.mat_immu_metag5

# Make nominal p-values (lower.tri) into NAs
p_adj_immu_metag5[lower.tri(p_adj_immu_metag5)] <- NA # don't need this step but I include to verify
p_adj_immu_metag5 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_immu_metag5[lower.tri(p_adj_immu_metag5)] = t(p_adj_immu_metag5)[lower.tri(p_adj_immu_metag5)]
p_adj_immu_metag5 

# Checking symmetry 
isSymmetric(p_adj_immu_metag5) #TRUE
#######


# Plot corr_immu with hclust - square, no p-values
dev.off()
png("plots/metagenome_and_immune/heatmap_immune_metag_allPs.png", res=600, units="in", width=10, height=10)
corrplot(corr_immu_metag5, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



# Plot corr_immu_metag5 with hclust - leave blank on non-significant coefficient & add text of significant correlation coefficients
png("plots/metagenome_and_immune/heatmap_immune_metag_FDR_sig.png", res=600, units="in", width=10, height=10)
corrplot(corr_immu_metag5, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_immu_metag5 , sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.4, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
##################################################################


################################################################## Different x & y axes
# Setting up vectors
names_immu <- names(immu_small)[-1]
names_immu
names_MG_chosen <- names(MG_RA_t_chosen5)
names_MG_chosen

# I want cytokines on y & metag on x

# Filtering corr_immu_metag5 
str(corr_immu_metag5)
corr_immu_metag5_asym <- as.matrix(as.data.frame(corr_immu_metag5)[names_MG_chosen, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as intput for corrplot
corr_immu_metag5_asym 

# Filtering p_adj_immu_metag5
p_adj_immu_metag5_asym <- as.matrix(as.data.frame(p_adj_immu_metag5)[names_MG_chosen, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as intput for corrplot
p_adj_immu_metag5_asym

# Finding n for captions
n_immu_metag5_asym <- as.matrix(as.data.frame(ct_immu_metag5$n)[names_MG_chosen, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as input for corrplot
n_immu_metag5_asym
min(n_immu_metag5_asym)
max(n_immu_metag5_asym)
# 79-80


# Plot asymmetrical matrix immune vs metag 
?corrplot
pdf("plots/metagenome_and_immune/asym_heatmap_immune_metag_FDR_sig.pdf", width=6, height=7)
corrplot(corr_immu_metag5_asym, #order = "hclust" - cannot do hclust on rectangles as it cuts some variables to make it into a square matrix
         p.mat = p_adj_immu_metag5_asym, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.8, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.6, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
##################################################################






###########################################################################################################################
##                                                                                                                       ##
##                                               Metabolites vs taxa                                                     ##
##                                                                                                                       ##
###########################################################################################################################

################################################################### Prep & All together heatmap (same x&y axes)
# Looking
MG_RA_t_chosen5[1:10,1:7]
metabolites_df[1:10,]

# Merge prep - Making Participant.ID column
MG_pre_merge_metabolites <- MG_RA_t_chosen5
MG_pre_merge_metabolites$Participant.ID <- substring(rownames(MG_RA_t_chosen5), 4) #keep characters after the 4th character
MG_pre_merge_metabolites$Participant.ID
MG_pre_merge_metabolites[1:5, 21:23]
dim(metabolites_df) #87 30 - but only 84 women have metabolite data
dim(MG_pre_merge_metabolites) #87 23

# Merge 
?base::merge
metabolites_some_metag <- base::merge(metabolites_df, MG_pre_merge_metabolites, by="Participant.ID", 
                                      all.x=T, all.y=T) # keep all women

# Tidy & check
dim(metabolites_some_metag) #87 52
row.names(metabolites_some_metag) <- metabolites_some_metag$Participant.ID
metabolites_some_metag$Participant.ID <- NULL

head(metabolites_some_metag)
str(metabolites_some_metag)


# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_metabolites_metag5 <- psych::corr.test(metabolites_some_metag, adjust = "fdr")
corr_metabolites_metag5 <- ct_metabolites_metag5$r
p.mat_metabolites_metag5 <- ct_metabolites_metag5$p


####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_metabolites_metag5 <- p.mat_metabolites_metag5

# Make nominal p-values (lower.tri) into NAs
p_adj_metabolites_metag5[lower.tri(p_adj_metabolites_metag5)] <- NA # don't need this step but I include to verify
p_adj_metabolites_metag5 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_metabolites_metag5[lower.tri(p_adj_metabolites_metag5)] = t(p_adj_metabolites_metag5)[lower.tri(p_adj_metabolites_metag5)]
p_adj_metabolites_metag5 

# Checking symmetry 
isSymmetric(p_adj_metabolites_metag5) #TRUE
#######

# Plot corr_metabolites with hclust - square, no p-values
dev.off()
png("plots/metagenome_and_metabolites/heatmap_metabolites_metag_allPs.png", res=600, units="in", width=12, height=12)
corrplot(corr_metabolites_metag5, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



# Plot corr_metabolites_metag5 with hclust - leave blank on non-significant coefficient & add text of significant correlation coefficients
png("plots/metagenome_and_metabolites/heatmap_metabolites_metag_FDR_sig.png", res=600, units="in", width=12, height=12)
corrplot(corr_metabolites_metag5, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_metabolites_metag5, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.4, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
################################################################### 


################################################################### Different x & y axes
# Setting up vectors
names_metabolites <- names(metabolites_df)[-1]
names_metabolites



# Filtering corr_metabolites_metag5 
str(corr_metabolites_metag5)
corr_metabolites_metag5_asym <- as.matrix(as.data.frame(corr_metabolites_metag5)[names_metabolites, names_MG_chosen]) # needs to be a df to be filtered, but needs to be a matrix as intput for corrplot
corr_metabolites_metag5_asym 

# Filtering p_adj_metabolites_metag5
p_adj_metabolites_metag5_asym <- as.matrix(as.data.frame(p_adj_metabolites_metag5)[names_metabolites, names_MG_chosen]) # needs to be a df to be filtered, but needs to be a matrix as intput for corrplot
p_adj_metabolites_metag5_asym

# Finding n for captions
n_metabolites_metag5_asym <- as.matrix(as.data.frame(ct_metabolites_metag5$n)[names_metabolites, names_MG_chosen]) # needs to be a df to be filtered, but needs to be a matrix as input for corrplot
n_metabolites_metag5_asym
min(n_metabolites_metag5_asym)
max(n_metabolites_metag5_asym)
# 84

# Plot asymmetrical matrix metabolites vs metag 
?corrplot
pdf("plots/metagenome_and_metabolites/asym_heatmap_metabolites_metag_FDR_sig.pdf", width=7, height=9)
corrplot(corr_metabolites_metag5_asym, #order = "hclust" - cannot do hclust on rectangles as it cuts some variables to make it into a square matrix
         p.mat = p_adj_metabolites_metag5_asym, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.8, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.5, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
################################################################### 




