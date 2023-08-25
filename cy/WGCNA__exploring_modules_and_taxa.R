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


# This script looks at how the taxa are divided into modules based on the min module size. 
# This data comes from the end of step 2 in WGCNA.



###########################################################################################################################
##                                                                                                                       ##
##                                               Load data & df prep                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Load csvs - Hel matrix with all IDs
mod_assign <- read.csv("data/WGCNA/taxa_module_assignment_Hel_sft10_mm5_July2022.csv")


# Rename columns
names(mod_assign)
mod_assign <- mod_assign %>% 
  dplyr::rename(Module_assigned=Module_when_5minModuleSize)

 



###########################################################################################################################
##                                                                                                                       ##
##                                         Looking at taxa module assignment                                             ##
##                                                                                                                       ##
###########################################################################################################################

######################################################## Looking at most interesting taxa w/ our current bio knowledge
### Lactobacillus_crispatus
mod_assign %>% filter(Taxa=="Lactobacillus_crispatus") # greenyellow
mod_assign %>% filter(Module_assigned=="greenyellow")

### Lactobacillus_iners 
mod_assign %>% filter(Taxa=="Lactobacillus_iners") #salmon
mod_assign %>% filter(Module_assigned=="salmon") 

### Lactobacillus_acidophilus 
mod_assign %>% filter(Taxa=="Lactobacillus_acidophilus")
mod_assign %>% filter(Module_assigned=="green") 

### Gardnerella_vaginalis
mod_assign %>% filter(Taxa=="Gardnerella_vaginalis")
mod_assign %>% filter(Module_assigned=="red") 

# Note: Previous analysis on raw & Hellinger Transformed raw counts showed that the min mod size doesn't seem to impact if there are high correlations b/ MEs & metadata

########################################################




######################################################## Looking at grey module in each mm condition
# Grey (the unassigned taxa module) taxa
mod_assign %>% filter(Module_assigned=="grey")

# Lazier c&p
grey <- mod_assign %>% filter(Module_assigned=="grey")
vec_grey <- gsub("_"," ", grey$Taxa, fixed = TRUE)
vec_grey
#Chryseobacterium gleum, Clostridium botulinum, Corynebacterium urealyticum, Granulicatella adiacens, Lactobacillus amylolyticus, Lactobacillus amylovorus, Lactobacillus coleohominis, Lactobacillus kefiranofaciens, Mycobacterium parascrofulaceum, Pasteurella multocida, Sphingobacterium spiritivorum, Staphylococcus hominis, Staphylococcus simulans 

########################################################








###########################################################################################################################
##                                                                                                                       ##
##                                                Looking at module size                                                 ##
##                                                                                                                       ##
###########################################################################################################################


## Loading WGCNA Rdata created in the script "WGCNA_on_metagenome.R"
# Loading .Rdata saved
lnames = load(file="data/WGCNA/WGCNA_ALL_FINAL_DATA.RData")

# Number of taxa per module
table(moduleColors)
# black        blue       brown       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#    15          36          19          17           6          13          10          15           8          17           6           6          42          17

# Module size stats
module_size <- as.data.frame(table(moduleColors)) %>% dplyr::rename(Size=Freq, Module=moduleColors)
module_size
summary(module_size)






###########################################################################################################################
##                                                                                                                       ##
##                                 Getting stats on module size (for thesis results)                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Look
head(module_size)
dim(module_size)

# Exclude grey
module_size_exclGrey <- module_size %>% 
  dplyr::filter(Module != "grey")
module_size_exclGrey

# Getting stats of valid modules
summary(module_size_exclGrey)







###########################################################################################################################
##                                                                                                                       ##
##                                  Table of modules & what taxa are in them                                             ##
##                                                                                                                       ##
###########################################################################################################################

# Look
table(moduleColors)
dim(table(moduleColors)) #14 #this & above match WGCNA script 

# Make the df
taxa_in_modules <- as.data.frame(cbind(names(datExpr), moduleColors)) %>%
  dplyr::rename(Taxa=V1,
                Module = moduleColors) %>% 
  arrange(Taxa)
head(taxa_in_modules)


# Nicer df format 
taxa_in_modules_nice <- taxa_in_modules %>% group_by(Module) %>% 
  dplyr::mutate(Taxa = paste(Taxa, collapse=", ")) %>%
  dplyr::distinct(Module, Taxa) %>% 
  dplyr::select(Module, Taxa)%>% 
  arrange(Module)

# Checking
taxa_in_modules %>% filter(Module=="grey") 
taxa_in_modules_nice %>% filter(Module=="grey") #matching 


# Putting grey as the last row in the df
taxa_in_modules_nice 
taxa_in_modules_nice <- taxa_in_modules_nice[order(taxa_in_modules_nice$Module %in% "grey"), ]
taxa_in_modules_nice 

# Change _ to space
taxa_in_modules_nice$Taxa <- gsub("_"," ", taxa_in_modules_nice$Taxa, fixed=T)
taxa_in_modules_nice

# Export to csv
write_csv(taxa_in_modules_nice, "data/WGCNA/Taxa_module_assignment__thesis_format.csv")


# Convert to latex for c&p
# Export
?knitr::kable
taxa_in_modules_nice %>%
  knitr::kable(caption = "Taxa clustered into the wgcna modules",   
               format = "latex",
               #digits=3 , 
               booktabs = TRUE) 
# Note: I had to convert to longtable not table as it is v. long.


