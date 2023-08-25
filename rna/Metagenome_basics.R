# Clear environment
rm(list=ls())

# Load libraries
library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(dplyr)
library(tidyr)

library(microbiome) #BiocManager::install("microbiome")
library(knitr)
# library(phyloseq)
library(reshape2) #heatmaps
library(vegan) #diversity calculations
library(ggpmisc) #plotting with text labels
library(ggrepel) #plotting with text labels
library(ggpubr) #stats tests
library(corrplot)


#Versions of packages used
#sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()





###########################################################################################################################
##                                                                                                                       ##
##                                                 Importing data                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Loading .Rdata saved
lnames <- load(file="data/metadata_neutrophils_with_VIRGO.RData")
lnames

# The relative abundance files in %s as they add up to 100
# I'm calling it VIRGO_t just to be consistent with how mg_t and mg were
VIRGO_t <- read.table("data/metagenome/Flavia_4_Virgo_alicia/summary.Percentage.txt", header=T)
# I then did data cleaning/wrangling after importing this

# Importing counts
VIRGO_c_t <- read.table("data/metagenome/Flavia_4_Virgo_alicia/summary.Count.txt", header=T)
# I then did data cleaning/wrangling after importing this





###########################################################################################################################
##                                                                                                                       ##
##                                           df prep for stacked bars plots                                              ##
##                                                                                                                       ##
###########################################################################################################################

############ Deleting taxa columns where there is not at least 1 person with an abundance of >10% for that genus & add an others column

# df with threshold of 10%
VIRGO_t_filtered10 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 10) %>%     ##### At 10%
  mutate(Other = 100-(rowSums(.[,], na.rm=TRUE))) ###### Making an Others column for remaining taxa
dim(VIRGO_t_filtered10) #9 7 - same dim as metaphlan 

# df with threshold of 5%
VIRGO_t_filtered5 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 5) %>%     ##### At 5%
  mutate(Other = 100-(rowSums(.[,], na.rm=TRUE))) ###### Making an Others column for remaining taxa 
dim(VIRGO_t_filtered5) # 9 7 so it's 1 < metaphlan

# df with threshold of 1%
VIRGO_t_filtered1 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 1) %>%     ##### At 1%
  mutate(Other = 100-(rowSums(.[,], na.rm=TRUE))) ###### Making an Others column for remaining taxa
dim(VIRGO_t_filtered1) #9 10 - so it's 2 < metaphlan




            
# I will go with 1% for the stacked bars...

# The taxa at 1%:
names(VIRGO_t_filtered1)
# [1] "Atopobium_vaginae"       "Bifidobacterium_breve"   "Gardnerella_vaginalis"   "Lactobacillus_crispatus" "Lactobacillus_gasseri"   "Lactobacillus_iners"    
# [7] "Lactobacillus_jensenii"  "Lactobacillus_vaginalis" "Prevotella_bivia"        "Other"



# Fixing the Other column
str(VIRGO_t_filtered1)
VIRGO_t_filtered1$Other <- as.numeric(VIRGO_t_filtered1$Other)
str(VIRGO_t_filtered1)

# Adding the IDs as a proper column
VIRGO_t_filtered1$IDs_2 <- rownames(VIRGO_t_filtered1)
VIRGO_t_filtered1$IDs_2





#Only choosing interesting thing from the metadata
metadata_CVF_basic <- metadata_CVF %>% 
  dplyr::select(IDs, Participant.ID, IDs_1, IDs_2, Lane, Gestation.at.visit.wks.dec,
                Ethnicity, Age_category, BMI_category, Smoking, Intervention,
                Gestation.at.delivery.wks.dec, Outcome, 
                Did.she.receive.antibiotics., Antibiotic_indication, Maternal.Infection, Immaturity_score, 
                Ravel, PCoA_Erica, Alicia_group) # I made these 3 variables later on in this exact script lol
dim(metadata_CVF_basic)

#Data type - fixing / making better for plots
metadata_CVF_basic$IDs_2
metadata_CVF_basic$Participant.ID <- factor(metadata_CVF_basic$Participant.ID)
#metadata_CVF_basic$Intervention <- as.character(metadata_CVF_basic$Intervention)

# Merging the 2 dfs
VIRGO_f1 <- merge(VIRGO_t_filtered1, metadata_CVF_basic,  by="IDs_2")
dim(VIRGO_f1) #9 30
VIRGO_f1$IDs_2




#gathering Taxa into 1 column
VIRGO_f1_gather <- VIRGO_f1 %>%
  gather(key=Taxa, value=Abundance, 
         -IDs, -Participant.ID, -IDs_1, -IDs_2, -Lane, -Gestation.at.visit.wks.dec,
         -Ethnicity, -Age_category, -BMI_category, -Smoking, -Intervention,
         -Gestation.at.delivery.wks.dec, -Outcome, 
         -Did.she.receive.antibiotics., -Antibiotic_indication, -Maternal.Infection,
         -Immaturity_score, 
         -Ravel, -PCoA_Erica, -Alicia_group)

dim(VIRGO_f1_gather) #90 22
str(VIRGO_f1_gather)
names(VIRGO_f1_gather)
table(VIRGO_f1_gather$Taxa, useNA="always")

# Getting rid of underscores for prettier plotting
str(VIRGO_f1_gather$Taxa)
VIRGO_f1_gather$Taxa <- gsub("_"," ", VIRGO_f1_gather$Taxa)

#putting Other category to the end of the factors - https://forcats.tidyverse.org/reference/fct_relevel.html
VIRGO_f1_gather$Taxa <- as.factor(VIRGO_f1_gather$Taxa)
VIRGO_f1_gather$Taxa <- fct_relevel(VIRGO_f1_gather$Taxa, "Other", after = Inf)









###########################################################################################################################
##                                                                                                                       ##
##                        Plotting STACKED BARS of abundance at species level faceted by demographics                    ##
##                                                                                                                       ##
###########################################################################################################################

# Colours in metaphlan mg plots:
# scale_fill_manual(values = c("#FF0000", "yellow", "#00CC99", #Ato, Dia, Gard
#                              "#B2F5EB", "#02FEF4", "#99CCFF", #crisp, gas, iners
#                              "#0099FF", "#0003FF", "#F73A95", #jen, para, phage
#                              "#BF64FD", "#9A7BF9",  "black")) #Prev b, Prev t, Other

table(VIRGO_f1_gather$Taxa)



# plot faceted by Cluster group, which is also Ethnicity 
Cluster.labs <- c("Cluster 1", "Clusters 2-4")
names(Cluster.labs) <- c("White", "Black")
bar_Cluster <- ggplot(VIRGO_f1_gather, aes(fill=Taxa, y=Abundance, x=IDs_1)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~Ethnicity, scales = "free", space = "free",
             labeller=labeller(Ethnicity=Cluster.labs)) +
  labs(x="Participant", y="Relative abundance (%)") + 
  scale_fill_manual(values = c("#FF0000", "yellow", "#00CC99", #Ato, Dia, Gard
                               "#B2F5EB", "#02FEF4", "#99CCFF", #crisp, gas, iners
                               "#0099FF", "#000EA1", "#BF64FD", "black")) #jen, vag Prev b,  Other
ggsave("plots/metagenome/VIRGO/abundance_species/bar_Cluster_group.pdf", bar_Cluster, height=3.5, width=8)



# plot faceted by Ravel FOR THESIS
CST.labs <- c("CST I", "CST III", "CST IV", "Unclassified")
names(CST.labs) <- c("CST I", "CST III", "CST IV-B", "Unclassified")
bar_Ravel <- ggplot(VIRGO_f1_gather, aes(fill=Taxa, y=Abundance, x=IDs_1)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~Ravel, scales = "free", space = "free",
             labeller=labeller(Ravel=CST.labs)) +
  labs(x="Participant", y="Relative abundance (%)") +
  scale_fill_manual(values = c("#FF0000", "yellow", "#00CC99", #Ato, Dia, Gard
                               "#B2F5EB", "#02FEF4", "#99CCFF", #crisp, gas, iners
                               "#0099FF", "#000EA1", "#BF64FD", "black")) #jen, vag Prev b,  Other
ggsave("plots/metagenome/VIRGO/abundance_species/bar_Ravel.pdf", bar_Ravel, height=3.5, width=8)

# Working out taxa number in Unclassified CST person
VIRGO_f1_gather %>% 
  dplyr::filter(IDs_1=="3") %>% dplyr::select(IDs_2, Taxa, Abundance)






###########################################################################################################################
##                                                                                                                       ##
##                                        Looking at Lactobacillus_acidophilus                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Preview
VIRGO_t$Lactobacillus_acidophilus

# df of just acidophilus
acidophilus <- VIRGO_t %>% 
  dplyr::select(Lactobacillus_acidophilus)

# seeing who has > 0 acidophilus
acidophilus %>% 
  dplyr::filter(Lactobacillus_acidophilus>0)
#      Lactobacillus_acidophilus
# ID_4                0.09691910
# ID_2                0.00328423
# ID_3                0.00127836

# Add IDs
acidophilus$IDs_2 <- rownames(acidophilus)

# Merging
acidophilus_gg <- merge(acidophilus, metadata_CVF, by="IDs_2")

# Change str
acidophilus_gg$Participant.ID <- factor(acidophilus_gg$Participant.ID)

# Lactobacillus acidophilus for thesis 
library(ggtext) 
L_acidophilus_for_thesis <- ggplot(acidophilus_gg, aes(x=IDs_1, y=Lactobacillus_acidophilus))  +
  geom_point() +
  geom_text(label=round(acidophilus_gg$Lactobacillus_acidophilus ,3), vjust=-1) +
  labs(y="*L. acidophilus* relative abundance (%)", x="Participant")  +
  scale_color_manual(values = c("#F68F44", "purple")) +
  scale_y_continuous(breaks=seq(0, 0.16, 0.04), minor_breaks=seq(0, 0.16, 0.02), limits=c(0,0.12)) +
  theme(axis.title.y = element_markdown(size=9),
        ) #https://wilkelab.org/ggtext/articles/introduction.html
ggsave("plots/metagenome/VIRGO/abundance_species/L_acidophilus_for_thesis.pdf", L_acidophilus_for_thesis, height=2.5, width=8)








  
###########################################################################################################################
##                                                                                                                       ##
##                                     Correlation heatmap of taxa                                                       ##
##                                                                                                                       ##
###########################################################################################################################

# When I used VIRGO_t it is so big that it is not useful. So I'm going to threshold the taxa
dim(VIRGO_t) #9 135


# df with threshold of 0.5%
VIRGO_t_0.5 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 0.5) 
dim(VIRGO_t_0.5) #9 13 
names(VIRGO_t_0.5)

# df with threshold of 0.25%
VIRGO_t_0.25 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 0.25)
dim(VIRGO_t_0.25) #9 15 
names(VIRGO_t_0.25) # includes BVAB1

# df with threshold of 0.1%
VIRGO_t_0.1 <- VIRGO_t %>% 
  select_if(~max(., na.rm = TRUE) >= 0.1) 
dim(VIRGO_t_0.1) #9 27 
names(VIRGO_t_0.1) # Lactobacillus_acidophilus has a max of 0.09691910 so it is still not coming up!



##### Let's do 0.25% plus Lactobacillus_acidophilus
# Select specific taxa we want to keep
selected_taxa <- VIRGO_t %>% 
  dplyr::select(Lactobacillus_acidophilus)  
selected_taxa 

# Bind these together
VIRGO_t_chosen <- cbind(VIRGO_t_0.25, selected_taxa)
VIRGO_t_chosen

# Putting df in alphabetical order 
VIRGO_t_chosen <- VIRGO_t_chosen %>% 
  dplyr::select(order(colnames(.)))
VIRGO_t_chosen
names(VIRGO_t_chosen)


# This will make a much more manageable heatmap now 
dim(VIRGO_t) #9 135
dim(VIRGO_t_chosen) #9 16






################################### Heatmap of my data
#install.packages('corrplot')

# Info
?corrplot
?psych::corr.test
citation('corrplot')



# This little chunk is a c&p from above in case I delete the previous ggcorrplot section
# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_VIRGO <- psych::corr.test(VIRGO_t_chosen, adjust = "fdr") 
corr_VIRGO <- ct_VIRGO$r
p.mat_VIRGO <- ct_VIRGO$p


####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_VIRGO <- p.mat_VIRGO

# Make nominal p-values (lower.tri) into NAs
p_adj_VIRGO[lower.tri(p_adj_VIRGO)] <- NA # don't need this step but I include to verify
p_adj_VIRGO 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_VIRGO[lower.tri(p_adj_VIRGO)] = t(p_adj_VIRGO)[lower.tri(p_adj_VIRGO)]
p_adj_VIRGO 

# Checking symmetry 
isSymmetric(p_adj_VIRGO) #TRUE
#######


# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
## leave blank on non-significant coefficient
## add significant correlation coefficients
pdf("plots/metagenome/VIRGO/corrplot_taxa_sig.pdf", width=6, height=6)
corrplot(corr_VIRGO, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_VIRGO, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 0.7, tl.col = 'black',
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.5, number.digits = 2) 
dev.off()










###########################################################################################################################
##                                                                                                                       ##
##                                  Alpha (within sample) diversity on VIRGO counts                                      ##
##                                                                                                                       ##
###########################################################################################################################

# At first vegan was not working, but now that I've reinstalled it, it works fine
library(vegan) #install.packages("vegan")
# This is vegan 2.5-7
?vegan::diversity


############################################  Calculating alpha diversity....
# Shannon - alpha diversity
Shannon <- vegan::diversity(VIRGO_c_t, index="shannon")
Shannon
# Shannon_RA <- vegan::diversity(VIRGO_t, index="shannon")
# Shannon_RA

# Simpson - alpha diversity
Simpson <- vegan::diversity(VIRGO_c_t, index="simpson")
Simpson

# Inverse_Simpson - alpha diversity
Inverse_Simpson <- vegan::diversity(VIRGO_c_t, index="invsimpson")
Inverse_Simpson
############################################  




############################################ Export alpha data

# Combine alpha diversity scores
alpha <- as.data.frame(cbind(Shannon, Simpson, Inverse_Simpson))
alpha
str(alpha)

# Add IDs as a column
alpha$IDs_2 <- rownames(alpha)
alpha

### Saving alpha to add to metadata_CVF!
save(alpha, file="data/metagenome/alpha_diversity_df_VIRGO.RData")

## Loading .Rdata saved 
# lnames = load(file="data/metagenome/alpha_diversity_df_VIRGO.RData")
# lnames

############################################




############################################  Plotting prep

# Add in ethnicity data
for_alpha <- metadata_CVF %>% 
  dplyr::select(IDs_1, IDs_2, Ethnicity) 
for_alpha 

# Merge ethnicity with alpha data
alpha <- merge(alpha, for_alpha, by="IDs_2")
alpha

# Gather for plotting
alpha_gathered <- alpha %>% 
  gather(key=diversity_index, value=diversity_score, -IDs_1, -IDs_2, -Ethnicity)
head(alpha_gathered)

# Sorting out diversity_index
str(alpha_gathered)
alpha_gathered$diversity_index_copy <- alpha_gathered$diversity_index #making a copy of the column
alpha_gathered$diversity_index[alpha_gathered$diversity_index == "Inverse_Simpson"] <- "Inverse Simpson" #replacing
alpha_gathered$diversity_index <- factor(alpha_gathered$diversity_index, levels=c("Shannon", "Simpson", "Inverse Simpson")) # Ordering factor
table(alpha_gathered$diversity_index, alpha_gathered$diversity_index_copy, useNA="always")
alpha_gathered$diversity_index_copy <- NULL

############################################  




############################################  Plotting

# Boxplot of all alpha diversity FOR THESIS
labels_alpha <- c("Shannon's index", "Simpson's index", "Inverse Simpson's index")
names(labels_alpha) <- c("Shannon", "Simpson", "Inverse Simpson") #alpha_gathered$diversity_index
alpha_Cluster_group <- ggplot(alpha_gathered, aes(x=Ethnicity, y=diversity_score))  + 
  geom_boxplot() +
  #ylim(0,12.5) +
  geom_jitter(aes(colour=Ethnicity), alpha=0.6, width=0.3, height=0) + #hollow circles and small jitter width. Don't change height when jittering!
  scale_color_manual(values = c("#FF226F", "blue")) +
  facet_wrap(~diversity_index, strip.position="top", nrow=1,
             labeller=labeller(diversity_index=labels_alpha), scales="free") +
  labs(y="Alpha diversity", x = "Transcriptomics cluster group") +
  theme(legend.position="none") +
  ggpubr::stat_compare_means(comparisons = list( c("White", "Black")),
                             method = "wilcox.test", hide.ns = T, paired=F, label="p.format", vjust=-0.7) +
  scale_x_discrete(labels=c("White" = "Cluster 1", "Black" = "Clusters 2-4")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) # Add 10% spaces between the p-value labels and the plot border
ggsave("plots/metagenome/VIRGO/diversity/alpha_Cluster_group.pdf",alpha_Cluster_group, height=4, width=8)


############################################ 







###########################################################################################################################
##                                                                                                                       ##
##                                        PCoA on the metagenome (Bray-Curtis)                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Beta diversity = between sample diversity

# Calculate Bray-Curtis distance metric on Hellinger VIRGO data
?vegan::vegdist
Bray_Curtis_detailed <- vegdist(VIRGO_Hel, index = "bray", diag=T, upper=T)

# Making it into a matrix
str(Bray_Curtis_detailed)
Bray_matrix <- as.matrix(Bray_Curtis_detailed)
Bray_matrix

# PCoA is not included in vegan, so will use the ape package instead
library(ape)
?ape::pcoa
PCoA <- pcoa(Bray_matrix)


# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCoA <- pcoa(Bray_matrix, correction = "cailliez")

# Plot your results
biplot.pcoa(PCoA)



# Extract the plot scores from first two PCoA axes (if you need them):
PCoAaxes <- as.data.frame(PCoA$vectors[,c(1,2)])
PCoAaxes

# Add IDs to PCoAaxes for merging
PCoAaxes$IDs_2 <- base::rownames(PCoAaxes)

# Rename columns
PCoAaxes <- PCoAaxes %>% 
  dplyr::rename(PCoA1 = Axis.1,
                PCoA2 = Axis.2)
PCoAaxes

###Saving the PCoA axes to add to metadata_CVF!
save(PCoAaxes, file="data/metagenome/PCoAaxes_based_on_Bray-Curtis_VIRGO.RData")

## Loading .Rdata saved 
# lnames = load(file="data/metagenome/PCoAaxes_based_on_Bray-Curtis_VIRGO.RData")
# lnames

# Check that PCoA from this script matches the PCoA axes in the metadata_CVF df (i.e. to make sure we aren't using metaphlan or an older non-Hellinger version of PCoA)
metadata_CVF %>% dplyr::select(IDs_2, PCoA1, PCoA2) #matching


# Plot PCoA plot for paper by Alicia_group & Ethnicity
PCoA_for_paper <- ggplot(metadata_CVF, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Alicia_group, shape=Ethnicity), alpha=1, size=4, stroke=1.5, na.rm = T) +
  scale_color_manual(values = c("#6ED088", "#99CCFF", "blue", "black", "red", "#FCA64C"), na.translate=F) +
  scale_shape_manual(values=c(3,16), na.translate=F) +
  labs(subtitle="PCoA plot of the Bray-Curtis dissimilarity of the metagenome",
       colour="Dominant microbiota", x="PCoA axis 1", y="PCoA axis 2")
PCoA_for_paper
ggsave("plots/metagenome/VIRGO/diversity/PCoA_for_paper_VIRGO.png", PCoA_for_paper, height=4, width=6.25, dpi=700)

# Plot PCoA plot FOR THESIS
PCoA_for_thesis <- ggplot(metadata_CVF, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Alicia_group, shape=Ethnicity), alpha=1, size=3, stroke=1.5, na.rm = T) +
  scale_color_manual(values = c("#6ED088", "#99CCFF", "blue", "black", "red", "#FCA64C"), na.translate=F) +
  scale_shape_manual(values=c(3,16), na.translate=F) +
  labs(colour="Dominant microbiota", x="PCoA axis 1", y="PCoA axis 2") +
  guides(colour = guide_legend(order = 1)) 
ggsave("plots/metagenome/VIRGO/diversity/PCoA_for_thesis.pdf", PCoA_for_thesis, height=4, width=6.25)
#####################################################







###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting all data                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data
save.image(file = "data/ALL_objects_in_Metagenome_VIRGO_basics.RData", compress=F)


# # Loading .Rdata saved
# lnames = load(file="data/ALL_objects_in_Metagenome_VIRGO_basics.RData")
# lnames





