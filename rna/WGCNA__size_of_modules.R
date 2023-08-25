# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(DESeq2) #for assay()
library(WGCNA)
library(GO.db)
library(org.Hs.eg.db) # Human genes

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()




###########################################################################################################################
##                                                                                                                       ##
##                                                     Load data                                                         ##
##                                                                                                                       ##
###########################################################################################################################

## This .Rdata is from the script "WGCNA_pipeline__neutrophils_VIRGO.R"
# Loading .Rdata saved
lnames = load(file="data/WGCNA/WGCNA_essentials_rld_with_VIRGO_in_datTraits.RData")
lnames





###########################################################################################################################
##                                                                                                                       ##
##                                         WGCNA module size for csv                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Look at module size
table(moduleColors)

# Save as df
module_size <- as.data.frame(table(moduleColors))
module_size

# Rename colnames
names(module_size)
names(module_size) <- c("module_colour", "module_size")
names(module_size)
head(module_size)

# Looking
summary(module_size$module_size)

# Exluding the grey module
modules <- module_size %>% 
  dplyr::filter(module_colour!="grey")

# Check
dim(module_size) #45
dim(modules) #44

# Looking
summary(modules$module_size)
boxplot(modules$module_size)
dev.off()






###########################################################################################################################
##                                                                                                                       ##
##                                 Getting stats on module size (for thesis results)                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Number of taxa per module
table(moduleColors)
# black        blue       brown       green greenyellow        grey     magenta        pink      purple         red      salmon         tan   turquoise      yellow 
#    15          36          19          17           6          13          10          15           8          17           6           6          42          17

# Module size stats
module_size2 <- as.data.frame(table(moduleColors)) %>% dplyr::rename(Size=Freq, Module=moduleColors)
module_size2
summary(module_size2)

# Look
head(module_size2)
dim(module_size2)

# Exclude grey
module_size2_exclGrey <- module_size2 %>% 
  dplyr::filter(Module != "grey")
module_size2_exclGrey

# Getting stats of valid modules
summary(module_size2_exclGrey)


