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



# Run this to get rid of any existing plots until it is null
dev.off()






###########################################################################################################################
##                                                                                                                       ##
##                                             Load in data                                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
lnames = load(file="data/immune_CytokinesComplementsIgM_2023.RData")
lnames


# Loading metagenome Rdata created in the script "Metagenome_exploration.R"
lnames = load(file="data/metagenome_counts_and_RA.RData")
lnames






###########################################################################################################################
##                                                                                                                       ##
##                                   df prep - Making cytokine w metabolites df                                          ##
##                                                                                                                       ##
###########################################################################################################################

################################### metabolites prep
# Look
names(immu_MD)
head(immu_MD)

# Metabolites only
metabolites_df <- immu_MD %>% dplyr::select(Participant.ID, Acetate:Tryptophan)

# Check
names(metabolites_df) # in new metadata the unknown metabolites have been dropped
str(metabolites_df) #all num  (apart from Participant.ID ofc)
head(metabolites_df)

# Putting df in alphabetical order 
metabolites_df <- metabolites_df %>%
  dplyr::select(order(colnames(.))) %>%
  dplyr::select(Participant.ID, everything())
head(metabolites_df)

####################### Dealing with "..._t"
# [17/08 15:41] Hadingham, Alicia - Hi James, what does it mean if a metabolite has a "_t" at the end of its name? I was wondering what the difference between Betaine and Betaine_t was
# [17/08 16:06] Mason, James - Tentative
# [17/08 20:32] Flaviani, Flavia - yeah I knew anbotu the t . Vicky mentioned to keep just the one without
# [18/08 10:30] Flaviani, Flavia - I have kept cadaverine as suggested by Vicky

#So will delete Betaine_t 
metabolites_df$Betaine_t <- NULL

# Rename Cadaverine_t
metabolites_df <- metabolites_df %>%
  dplyr::rename(Cadaverine=Cadaverine_t)

# Check
names(metabolites_df) # in new metadata the unknown metabolites have been dropped
head(metabolites_df)
dim(metabolites_df) #87 30
#######################

# Save as RData
save(metabolites_df, file="data/metabolites_df.RData")
# ## Loading metabolite .Rdata created in the script "Metabolites_exploration_in_Black_women.R"
# lnames = load(file="data/metabolites_df.RData")
# lnames
###################################


################################### immu prep
# Look
names(immu_small)
head(immu_small)

###################################



################################### merge
# Look
dim(metabolites_df) # 87  30
dim(immu_small) # 81 14

# Merge
metab_immu <- base::merge(immu_small, metabolites_df, by="Participant.ID",
                           all.x=F, all.y=T)# only keep cytokine women where we have metab (x), & keep all women with metab (y) 
# Removing the immu IDs not in metag, to keep it Black women only.

# Look
dim(metab_immu) #87 43
head(metab_immu)

# Heatmap prep - setting rownames by Participant.ID & getting rid of Participant.ID col 
rownames(metab_immu) <- metab_immu$Participant.ID
head(metab_immu)
metab_immu$Participant.ID <- NULL
head(metab_immu)
###################################







###########################################################################################################################
##                                                                                                                       ##
##                                           Heatmap of metabolites only                                                 ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?corrplot
?psych::corr.test
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_metabolites <- psych::corr.test(metabolites_df[,-1], adjust = "fdr") #[,-1] to get rid of Participant.ID 
corr_metabolites <- ct_metabolites$r
p.mat_metabolites <- ct_metabolites$p

# Finding n for captions
ct_metabolites$n 
min(ct_metabolites$n)
max(ct_metabolites$n)
#84 always

####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_metabolites <- p.mat_metabolites

# Make nominal p-values (lower.tri) into NAs
p_adj_metabolites[lower.tri(p_adj_metabolites)] <- NA # don't need this step but I include to verify
p_adj_metabolites 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_metabolites[lower.tri(p_adj_metabolites)] = t(p_adj_metabolites)[lower.tri(p_adj_metabolites)]
p_adj_metabolites 

# Checking symmetry 
isSymmetric(p_adj_metabolites) #TRUE
#######


# Plot corr_metabolites with hclust - square, no p-values
png("plots/metabolites_only/heatmap_metabolites_allPs.png", res=800, units="in", width=8, height=8)
corrplot(corr_metabolites, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



# Plot corr_metabolites with hclust - leave blank on non-significant coefficient & add text of significant correlation coefficients
pdf("plots/metabolites_only/heatmap_metabolites_FDR_sig.pdf", width=6, height=6)
corrplot(corr_metabolites, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_metabolites, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.7, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.3, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()








###########################################################################################################################
##                                                                                                                       ##
##                                          Heatmap of metabolites & immune                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?corrplot
?psych::corr.test
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html


################################################################### Prep & All together heatmap (same x&y axes)
# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_metab_immu <- psych::corr.test(metab_immu, adjust = "fdr") 
corr_metab_immu <- ct_metab_immu$r
p.mat_metab_immu <- ct_metab_immu$p



####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_metab_immu <- p.mat_metab_immu

# Make nominal p-values (lower.tri) into NAs
p_adj_metab_immu[lower.tri(p_adj_metab_immu)] <- NA # don't need this step but I include to verify
p_adj_metab_immu 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_metab_immu[lower.tri(p_adj_metab_immu)] = t(p_adj_metab_immu)[lower.tri(p_adj_metab_immu)]
p_adj_metab_immu 

# Checking symmetry 
isSymmetric(p_adj_metab_immu) #TRUE
#######

# Plot corr_metab_immu with hclust - square, no p-values
png("plots/metabolites_and_immune/heatmap_metabolites_and_immune_allPs.png", res=800, units="in", width=10, height=10)
corrplot(corr_metab_immu, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



# Plot corr_metab_immu with hclust - leave blank on non-significant coefficient & add text of significant correlation coefficients
png("plots/metabolites_and_immune/heatmap_metabolites_and_immune_FDR_sig.png", res=800, units="in", width=10, height=10)
corrplot(corr_metab_immu, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_metab_immu , sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.4, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
################################################################### 



################################################################### Different x & y axes
# Setting up vectors
names_immu <- names(immu_small)[-1]
names_immu
names_metab <- names(metabolites_df)[-1]
names_metab

# I want immune on y & metab on x

# Filtering corr_metab_immu 
str(corr_metab_immu)
corr_metab_immu_asym <- as.matrix(as.data.frame(corr_metab_immu)[names_metab, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as input for corrplot
corr_metab_immu_asym 

# Filtering p_adj_metab_immu
p_adj_metab_immu_asym <- as.matrix(as.data.frame(p_adj_metab_immu)[names_metab, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as input for corrplot
p_adj_metab_immu_asym

# Finding n for captions
n_metab_immu_asym <- as.matrix(as.data.frame(ct_metab_immu$n)[names_metab, names_immu]) # needs to be a df to be filtered, but needs to be a matrix as input for corrplot
n_metab_immu_asym
min(n_metab_immu_asym)
max(n_metab_immu_asym)
# 77-78

# Plot asymmetrical matrix immune vs metag 
?corrplot
pdf("plots/metabolites_and_immune/asym_heatmap_immune_metab_FDR_sig.pdf", width=4.5, height=7)
corrplot(corr_metab_immu_asym, #order = "hclust" - cannot do hclust on rectangles as it cuts some variables to make it into a square matrix
         p.mat = p_adj_metab_immu_asym, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.9, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.5, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()
################################################################### 


