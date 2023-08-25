# Clear environment
rm(list=ls())

# Load libraries
library(reshape2) # heatmaps
library(BiocManager) #install.packages("BiocManager")
library(tidyverse) # ggplot2, dplyr, tidyr etc
library(ggrepel)
library(stats) # for Wilcoxon or other stats tests
library(ggpubr) # for stats test in the plots
library(lubridate) # for working with dates


## Set wd
# setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files("data/larger_cohort")






###########################################################################################################################
##                                                                                                                       ##
##                                                 Importing data                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Loading df
immune_complete_df <- read.csv("data/larger_cohort/Staining_sort_masterfile__merge_120821_monocytes_edited.csv", header=T)
# Note 04-02-22: Now running the analysis with monocytes split into 3 categories
# Did some data cleaning & added new cols to this df, which I've removed from here







###########################################################################################################################
##                                                                                                                       ##
##                                                Plotting set up                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Creating facet labels - https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
Cell_type_labs <- c("Neutrophils", "NK cells", "B cells", "T cells", "Patrolling monocytes", "Intermediate monocytes", "Classical monocytes")
names(Cell_type_labs) <- c("Neutrophils", "NKcells", "Bcells", "Tcells", "Patrolling_monocytes", "Intermediate_monocytes", "Classical_monocytes")


# Colour scheme of cell types
fills_of_neutrophil <- scale_fill_manual(values=c("tomato", "#35C250", "#9C73E7", "#FA70B4", "#3434C4", "#5E82FF", "#8BB6F4"))
colours_of_neutrophil <- scale_color_manual(values=c("tomato", "#35C250", "#9C73E7", "#FA70B4", "#3434C4", "#5E82FF", "#8BB6F4")) 






###########################################################################################################################
##                                                                                                                       ##
##                                        Making dfs of only 1 sample/ID                                                 ##
##                                                                                                                       ##
###########################################################################################################################

# IDs with multiple samples taken
table(as.character(multi_samples$Participant.ID))# 3 IDs

# Making immune with only 1 samples/ID - and selecting the earliest sample in IDs with multiple samples
immune_uniq <- immune %>% 
  dplyr::arrange(Participant.ID, Gest.at.cyto.wks.dec) %>% 
  dplyr::distinct(Participant.ID, .keep_all = T) #.keep_all = T keeps all the other columns

# Checks
dim(immune) #51 12
dim(immune_uniq) #46 12
table(immune$Participant.ID)
table(immune_uniq$Participant.ID) # all 1s 
# View(immune_uniq)


# Gather df - https://datacarpentry.org/R-ecology-lesson/03-dplyr.html#Reshaping_with_gather_and_spread
immune_gather_uniq <- immune_uniq %>%
  gather(key = "Cell_type", value = "Percentage", 
         -Participant.ID, -Ethnicity, -Outcome, -Gest.at.cyto.wks.dec, -INSIGHT_cat)

# Changing str
str(immune_gather_uniq)
table(immune_gather_uniq$Cell_type)
immune_gather_uniq$Cell_type <- factor(immune_gather_uniq$Cell_type, levels=c("Neutrophils", "NKcells", "Bcells", "Tcells", "Patrolling_monocytes", "Intermediate_monocytes", "Classical_monocytes"))
table(immune_gather_uniq$Cell_type)
str(immune_gather_uniq)






###########################################################################################################################
##                                                                                                                       ##
##                                  Plotting mean % of cell types across all samples                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Get summary df of info on cell type across all samples
?summarise
summary_immune_AllTogether <- immune_gather %>% 
  dplyr::group_by(Cell_type) %>%
  dplyr::summarise(n = n(),
                   Mean_percentage = mean(Percentage))
summary_immune_AllTogether 

summary(immune$Neutrophils) # matches summary_immune_AllTogether 
summary(immune$NKcells) # matches summary_immune_AllTogether 

# Add nominal col for bar plotting functionality
summary_immune_AllTogether$nominal <- ""

# Prettied Cell_type_pretty column
summary_immune_AllTogether$Cell_type
summary_immune_AllTogether$Cell_type_pretty <- c("Neutrophils", "NK cells", "B cells", "T cells", "Patrolling monocytes", "Intermediate monocytes", "Classical monocytes") # c&p from Cell_type_labs
summary_immune_AllTogether$Cell_type_pretty <- factor(summary_immune_AllTogether$Cell_type_pretty, levels=c("Neutrophils", "NK cells", "B cells", "T cells", "Patrolling monocytes", "Intermediate monocytes", "Classical monocytes"))
summary_immune_AllTogether
# # A tibble: 7 x 5
#   Cell_type                  n Mean_percentage nominal Cell_type_pretty      
#   <fct>                  <int>           <dbl> <chr>   <fct>                 
# 1 Neutrophils               51         80.0    ""      Neutrophils           
# 2 NKcells                   51          4.72   ""      NK cells              
# 3 Bcells                    51          1.44   ""      B cells               
# 4 Tcells                    51          0.992  ""      T cells               
# 5 Patrolling_monocytes      51          0.0588 ""      Patrolling monocytes  
# 6 Intermediate_monocytes    51          0.0165 ""      Intermediate monocytes
# 7 Classical_monocytes       51          0.0996 ""      Classical monocytes  
0.0588 + 0.0165 + 0.0996 #0.1749


# Stacked bar looking at cell proportions & sampling time categories FOR THESIS
Cell_types_AllTogether_bar <- ggplot(summary_immune_AllTogether, aes(x=nominal, fill=Cell_type_pretty, y=Mean_percentage)) + 
  geom_bar(stat="identity", position="stack") + 
  labs(x=NULL, y="Mean percentage of CD45+ cells (%)", fill="CD45+ cells") + #subtitle="Mean percentage of CD45+ cells across all samples", fill="CD45+ cells"
  fills_of_neutrophil +
  ylim(c(0,100)) +
  geom_text(aes(x=nominal, label=ifelse(Mean_percentage>4, paste0(round(Mean_percentage, 0), "%"),"")), #means label if it is above x% otherwise do not label
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  theme(axis.ticks = element_blank())
ggsave("plots/cell_proportions_in_larger_cohort/Cell_types_mean_percentages.pdf", Cell_types_AllTogether_bar, height=3, width=3.5)






###########################################################################################################################
##                                                                                                                       ##
##                                       Plotting cell types by sampling gestation                                       ##
##                                                                                                                       ##
###########################################################################################################################

# Look
str(immune_gather)
names(immune_gather)
head(immune_gather)
tapply(immune$Neutrophils, immune$INSIGHT_cat, mean)
tapply(immune$NKcells, immune$INSIGHT_cat, mean)

# Get summary df of info on cell type by INSIGHT_cat
?summarise
summary_immune_visit <- immune_gather %>% 
  dplyr::group_by(INSIGHT_cat, Cell_type) %>%
  dplyr::summarise(n = n(),
                   Mean_percentage = mean(Percentage))
summary_immune_visit # matches my tapply above 



# Scatter of cell proportions by gestation at visit time  - with LM for THESIS
Labs_LM_gestation <- c("Neutrophils *", "NK cells", "B cells", "T cells *", "Patrolling monocytes", "Intermediate monocytes", "Classical monocytes")
names(Labs_LM_gestation) <- c("Neutrophils", "NKcells", "Bcells", "Tcells", "Patrolling_monocytes", "Intermediate_monocytes", "Classical_monocytes")
Cell_types_time_scatter_lm <- ggplot(immune_gather, aes(x=Gest.at.cyto.wks.dec, y=Percentage)) + 
  geom_point(aes(color=Cell_type))  + 
  facet_wrap(~Cell_type, strip.position="top", nrow=1, labeller=labeller(Cell_type=Labs_LM_gestation), scales="free") + #, scales="free"
  labs(x="Sampling time (weeks' gestation)", y="Percentage of CD45+ cells (%)", color="CD45+ cells") +
  colours_of_neutrophil +
  theme(legend.position = "none", strip.text.x = element_text(size = 8)) +
  geom_smooth(method='lm') # for some reason I cannot get it to be black instead of blue using aes(color=c(replicate(7, "black")))
ggsave("plots/cell_proportions_in_larger_cohort/Cell_types_time_lm.pdf", Cell_types_time_scatter_lm, height=3, width=11)


################################################# Linear regression - all non-sig
?lm # first variable is the thing we're interested in, second is the group/variable that we're wondering if it impacts the first variable
names(immune)

# Neutrophils
lm_Neutrophils <- lm(immune$Neutrophils ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Neutrophils) #0.0303 sig

# NKcells
lm_NKcells <- lm(immune$NKcells ~ immune$Gest.at.cyto.wks.dec)
summary(lm_NKcells) #0.206

# Bcells
lm_Bcells <- lm(immune$Bcells ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Bcells) #0.241

# Tcells
lm_Tcells <- lm(immune$Tcells ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Tcells) #0.0165 sig

# Patrolling_monocytes
lm_Patrolling_monocytes <- lm(immune$Patrolling_monocytes ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Patrolling_monocytes) #0.883

# Intermediate_monocytes
lm_Intermediate_monocytes <- lm(immune$Intermediate_monocytes ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Intermediate_monocytes) #0.195

# Classical_monocytes
lm_Classical_monocytes <- lm(immune$Classical_monocytes ~ immune$Gest.at.cyto.wks.dec)
summary(lm_Classical_monocytes) #0.172


#################################################







###########################################################################################################################
##                                                                                                                       ##
##                                  Plotting cell types by Outcome & Ethnicity                                           ##
##                                                                                                                       ##
###########################################################################################################################

########################################################################## Ethnicity
# Boxplot of cell type percentages by Ethnicity FOR THESIS
set.seed(100)
Cell_types_Ethnicity_boxplot <- ggplot(subset(immune_gather_uniq, !is.na(Ethnicity)), aes(x=Ethnicity, y=Percentage))  + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(aes(colour=Ethnicity), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Cell_type, strip.position="top", nrow=1, labeller=labeller(Cell_type=Cell_type_labs), scales="free") +
  labs(y="Percentage of CD45+ cells (%)") +
  scale_color_manual(values=c("#F68F44", "purple")) +
  theme(legend.position="none", strip.text.x = element_text(size = 8)) #+
  # ggpubr::stat_compare_means(comparisons = list( c("White", "Black")),
  #                            method = "wilcox.test", hide.ns = T, paired=F, label="p.format", vjust=1.5) 
ggsave("plots/cell_proportions_in_larger_cohort/Cell_types_Ethnicity.pdf", Cell_types_Ethnicity_boxplot, height=4.5, width=11)


# Test to check that p-values on the plot is accurate
ggpubr::compare_means(Neutrophils ~ Ethnicity, data = immune_uniq, method = "wilcox.test") 

# Looking
tapply(immune_uniq$Neutrophil, immune_uniq$Ethnicity, summary)
##########################################################################



########################################################################## Outcome
# Boxplot of cell type percentages by Outcome FOR THESIS
set.seed(100)
Cell_types_Outcome_boxplot <- ggplot(immune_gather_uniq, aes(x=Outcome, y=Percentage))  + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(aes(colour=Outcome), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Cell_type, strip.position="top", nrow=1, labeller=labeller(Cell_type=Cell_type_labs), scales="free") +
  labs(y="Percentage of CD45+ cells (%)") +
  scale_color_manual(values = c("#00BFC4", "tomato", "grey50")) +
  theme(legend.position="none", strip.text.x = element_text(size = 8)) #+
  # ggpubr::stat_compare_means(comparisons = list( c("Term", "Preterm")),
  #                            method = "wilcox.test", hide.ns = T, paired=F, label="p.format", vjust=1.5) 
ggsave("plots/cell_proportions_in_larger_cohort/Cell_types_Outcome.pdf", Cell_types_Outcome_boxplot, height=4, width=11)



# Test to check that p-values on the plot is accurate
ggpubr::compare_means(Neutrophils ~ Outcome, data = immune_uniq, method = "wilcox.test") 

# Looking
tapply(immune_uniq$Neutrophils, immune_uniq$Outcome, summary)
tapply(immune_uniq$Bcells, immune_uniq$Outcome, summary)
##########################################################################









###########################################################################################################################
##                                                                                                                       ##
##                                                        n numbers                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Finding n numbers in the plots

########################################## 1 sample / ID
# Finding n numbers in the plots
table(immune_uniq$Ethnicity, useNA="always")
# White Black  <NA> 
# 27    15     4 

# Looking at the 4 missing ethnicity 
table(immune_uniq_all_demo$Ethnicity, immune_uniq_all_demo$Ethnicity_detailed, useNA="always") # the 4 are all Unclassified


table(immune_uniq$Outcome, useNA="always")
# Now:
# Term Preterm    <NA> 
#   40       6       0 
##########################################








###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting all data                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data
save.image(file = "data/ALL_objects_in_Immune_cell_types_in_larger_cohort.RData", compress=F)


# # Loading .Rdata from Immune_cell_types_in_larger_cohort
# lnames = load(file="data/ALL_objects_in_Immune_cell_types_in_larger_cohort.RData")
# lnames
