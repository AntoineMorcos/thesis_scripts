# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(WGCNA)
library(GO.db)
library(org.Hs.eg.db) # Human genes
library(emdbook) #lseq function

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()







###########################################################################################################################
##                                                                                                                       ##
##                                               Load data & df prep                                                     ##
##                                                                                                                       ##
###########################################################################################################################


# Loading .Rdata saved
lnames = load(file="data/WGCNA/WGCNA_essentials_rld_with_VIRGO_in_datTraits.RData")
lnames


# # GO - 1st method
# GOenr, tab, screenTab, #tab & screenTab are reduced versions of GOenr  
# 
# # GO - 2nd method
# GOenrichment, table.display, #table.display is a reduced version of GOenrichment




##### Method 1 - using the deprecated GOenrichmentAnalysis function from the WGCNA package
screenTab %>% 
  filter(module=="magenta") 


##### Method 2 - using the enrichmentAnalysis function from the anRichment package (WGCNA's recommended replacement [?enrichmentAnalysis])
head(table.display %>% 
       dplyr::select(class, dataSetID, dataSetName, pValue) %>% 
       dplyr::filter(class=="magenta"), 15)






###################################### New df of all GO

#### Making df for plotting based on data from method 2
str(table.display)

# Duplicating df to modify it
GO_All <- table.display

# GO_type - BP, MF, CC
GO_All$GO_type <- substring(GO_All$inGroups, 7, 8) #keep characters between 4 and 7 inclusive
table(GO_All$GO_type, GO_All$inGroups, useNA="always")
GO_All$inGroups <- NULL


# Trying to work out what classSize means
table(moduleColors)
table(GO_All$classSize, GO_All$class)

# Characters to factors
GO_All$class <- factor(GO_All$class)
GO_All$GO_type <- factor(GO_All$GO_type, levels=c("BP", "MF", "CC"))

######################################



###################################### New df of GO of Bonferroni sig results

# Look
names(GO_All)
dim(GO_All) #24224    18

# Filter df for magenta only (& arrange by pValue - even though I think it was already sorted like this anyway)
GO_bsig <- GO_All %>% 
  dplyr::filter(Bonferroni<=0.05) %>% 
  dplyr::arrange(class, pValue)

# Looking
dim(GO_bsig) #415  18
head(GO_bsig)
table(GO_bsig$class)

# Altering levels of GO_bsig$class to see module w/ sig GO more clearly 
GO_bsig$class <- as.character(GO_bsig$class)
GO_bsig$class <- factor(GO_bsig$class)
table(GO_bsig$class)
# black     darkgrey       grey60        ivory      magenta midnightblue         pink       purple 
# 103            6           27            3           78          100           89            9 

## Looking to see which modules have interesting things on the VIRGO heatmap (the heatmap with not many VIRGO taxa though)
# black          none
# darkgrey       none (not even in the metaphlan all taxa heatmap)
# grey60         none (not even in the metaphlan all taxa heatmap)
# ivory          positively w/ Lactobacillus_gasseri & Lactobacillus_vaginalis
# magenta        negatively w/ Gardnerella_vaginalis & positively with Lactobacillus_iners
# midnightblue   negatively w/ Immaturity score, but removed this score later     
# pink           none (not even in the metaphlan all taxa heatmap)
# purple         none (not even in the metaphlan all taxa heatmap)

######################################






###########################################################################################################################
##                                                                                                                       ##
##                                                   GO plot prep                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Reverse scales function https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# The x-axis scales I want in all the below GO plots so that they are comparable
x_axis_for_GO_plots <- scale_x_continuous(trans=reverselog_trans(10), minor_breaks=lseq(1e-01,1e-18,18), limits=c(1e-01,1e-14))
x_axis_for_GO_plots_e10 <- scale_x_continuous(trans=reverselog_trans(10), minor_breaks=lseq(1e-01,1e-10,10), limits=c(1e-01,1e-10))
# summary(GO_ivory_all$Bonferroni)
# summary(GO_magenta_all$Bonferroni)

# Faceting and theme I want in all the below GO plots 
faceting_for_GO_plots <- facet_grid(rows=vars(GO_type), space = "free_y", scales = "free") 
theme_for_GO_plots <- theme(legend.position="none",
        plot.title.position = "plot")





###########################################################################################################################
##                                                                                                                       ##
##                                             Exploring the magenta GO                                                  ##
##                                                                                                                       ##
###########################################################################################################################


###################################### df of GO of magenta module including only Bonferroni sig GO  
# Filter df for magenta only (& arrange by pValue - even though I think it was already sorted like this anyway)
GO_magenta_all <- GO_bsig %>% 
  dplyr::filter(class=="magenta") %>% 
  dplyr::arrange(pValue) %>% 
  dplyr::filter(dataSetName!="cellular_component" & dataSetName!="molecular_function" & dataSetName!="biological_process")

# Looking
dim(GO_magenta_all) # 75 18
head(GO_magenta_all)
names(GO_magenta_all)
summary(GO_magenta_all$Bonferroni)

# Sig GO terms
GO_magenta_all$dataSetName
######################################

###################################### magenta- biological process only
# BP only
GO_magenta_BP <- GO_magenta_all %>% filter(GO_type=="BP") #biological process
######################################


###################################### Looking
# Top 20 GO for Tan module 
GO_magenta_top <- GO_magenta_all[1:20,] # Top 20 GO terms

# Looking
dim(GO_magenta_top)
summary(GO_magenta_top$Bonferroni)



# All sig GO terms
GO_magenta_all %>% filter(GO_type=="BP") %>% dplyr::select(dataSetName)
GO_magenta_all %>% filter(GO_type=="MF") %>% dplyr::select(dataSetName)
GO_magenta_all %>% filter(GO_type=="CC") %>% dplyr::select(dataSetName)
# Potentially interesting ones:
# BP:
# neutrophil mediated immunity
# neutrophil activation
# neutrophil degranulation
# neutrophil activation involved in immune response
# granulocyte activation
# response to stimulus
# myeloid leukocyte mediated immunity
# leukocyte degranulation
# myeloid leukocyte activation
# myeloid cell activation involved in immune response
# leukocyte activation involved in immune response
# response to chemical
# cell activation involved in immune response
# regulated exocytosis
# cellular response to chemical stimulus
# cellular response to stimulus
# exocytosis
# cell activation
# leukocyte activation
# vesicle-mediated transport
# immune effector process
# leukocyte mediated immunity
# secretion
# secretion by cell
# response to hypoxia
# immune system process
# export from cell
# response to decreased oxygen levels
# response to oxygen levels
# MF:
# *nothing interesting*
# CC:
# tertiary granule
######################################



###################################### Plotting
# Plotting ALL sig p-values in magenta module
GO_magenta_all_plot <- ggplot(GO_magenta_all, aes(x=Bonferroni, y=reorder(dataSetName, -Bonferroni)))  + #order dataSetNames by Bonferroni p-values
  geom_point(colour="magenta") +
  labs(y=NULL, x="Bonferroni adjusted p-value", title="All statistically significant GO terms of the WGCNA magenta module") +
  faceting_for_GO_plots + 
  theme_for_GO_plots +
  x_axis_for_GO_plots 
ggsave("plots/WGCNA/VIRGO_in_datTraits/WGCNA_GO/GO_magenta_all_plot.pdf", GO_magenta_all_plot, height=15, width=6)


# Plotting ALL sig p-values for BP ONLY in magenta module - FOR THESIS
GO_magenta_BP$BP_full_name <- "Biological process GO"
GO_magenta_all_BP_plot <- ggplot(GO_magenta_BP, aes(x=Bonferroni, y=reorder(dataSetName, -Bonferroni)))  + #order dataSetNames by Bonferroni p-values
  geom_point(colour="magenta") +
  labs(y=NULL, x="Bonferroni adjusted p-value") +
  facet_grid(rows=vars(BP_full_name), space = "free_y", scales = "free")  + 
  theme_for_GO_plots +
  x_axis_for_GO_plots_e10 
ggsave("plots/WGCNA/VIRGO_in_datTraits/WGCNA_GO/GO_magenta_all_BP.pdf", GO_magenta_all_BP_plot, height=6, width=6)
######################################







###########################################################################################################################
##                                                                                                                       ##
##                                             Exploring the ivory GO                                                    ##
##                                                                                                                       ##
###########################################################################################################################


###################################### df of GO of ivory module including only Bonferroni sig GO  
# Filter df for ivory only (& arrange by pValue - even though I think it was already sorted like this anyway)
GO_ivory_all <- GO_bsig %>% 
  dplyr::filter(class=="ivory") %>% 
  dplyr::arrange(pValue) %>% 
  dplyr::filter(dataSetName!="cellular_component" & dataSetName!="molecular_function" & dataSetName!="biological_process")

# Looking
dim(GO_ivory_all) # 3 18
head(GO_ivory_all)
names(GO_ivory_all)
summary(GO_ivory_all$Bonferroni)

# Sig GO terms
GO_ivory_all$dataSetName #"type I interferon signaling pathway"    "cellular response to type I interferon" "response to type I interferon" 
######################################



# There's only 3 sig GO terms



###################################### Plotting
# Plotting ALL sig p-values in ivory module - FOR THESIS
GO_ivory_all$BP_full_name <- "GO"
GO_ivory_all_plot <- ggplot(GO_ivory_all, aes(x=Bonferroni, y=reorder(dataSetName, -Bonferroni)))  + #order dataSetNames by Bonferroni p-values
  geom_point(colour="black", fill="ivory", shape=21) + #ivory points are too difficult to see
  labs(y=NULL, x="Bonferroni adjusted p-value")  +
  facet_grid(rows=vars(BP_full_name), space = "free_y", scales = "free") + 
  theme(legend.position="none",
        plot.title.position = "plot") + #  axis.text.y = element_text(size = 10)
  x_axis_for_GO_plots_e10 
ggsave("plots/WGCNA/VIRGO_in_datTraits/WGCNA_GO/GO_ivory_all.pdf", GO_ivory_all_plot, height=1.2, width=5) #pdf
#######################################


    




###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting all data                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data
save.image(file = "data/ALL_objects_in_WGCNA_plotting_selected_GO_VIRGO.RData", compress=F)


# # Loading .Rdata saved
# lnames = load(file="data/ALL_objects_in_WGCNA_plotting_selected_GO_VIRGO.RData")
# lnames



