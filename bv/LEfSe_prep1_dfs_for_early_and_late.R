## Clear environment
rm(list=ls())

## Load libraries
library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(ggplot2)
library(dplyr)
library(tidyr)

## Set wd
setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
#setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
getwd()
#list.files()


## Loading data in: 
#lnames = load(file="data_BV/main_dfs.RData")
#lnames 
# This data only has women where we have 2 metabolite samples per person. So instead I will be going back to the df csv file itself, so not to reduce sample sizes unnecessarily. 


##Loading data in:
data1 <- read.csv("data_BV/Final_Db_04_07_2019_noIUD_Alicia_edit_CSV.csv", 
                     header=TRUE, na.strings=c("NA", " ", ""), fileEncoding = 'UTF-8-BOM')



##Checks
table(data1$participant_id)
dim(data1) #648 3176

#changing row names to Sample_IDs
rownames(data1) <- data1$Sample_ID
data1[1:5, 1:4]




####Making a new BV categories column, which is shorter to take up less room on the plots
table(data1$Bacterial.Vaginosis..categories) 
data1$BV_categories <- data1$Bacterial.Vaginosis..categories #making a copy of the df
data1$BV_categories <- as.character(data1$BV_categories)
data1$BV_categories[data1$BV_categories == "Bacterial vaginosis (3)"] <- "BV"
data1$BV_categories[data1$BV_categories == "Intermediate (2)"] <- "Intermediate"
data1$BV_categories[data1$BV_categories == "Normal (0-1)"] <- "Normal"
data1$BV_categories <- factor(data1$BV_categories, levels=c("Normal", "Intermediate", "BV"))
table(data1$BV_categories, data1$Bacterial.Vaginosis..categories)



###Adding new column for long ethnicity
table(data1$Ethnicity)
data1$Ethnicity_long <- data1$Ethnicity
data1$Ethnicity_long <- as.character(data1$Ethnicity_long)
data1$Ethnicity_long[data1$Ethnicity_long == "Black"] <- "Black ethnicity"
data1$Ethnicity_long[data1$Ethnicity_long == "White"] <- "White ethnicity"
data1$Ethnicity_long[data1$Ethnicity_long == "Other"] <- "Other ethnicity"
data1$Ethnicity_long <- factor(data1$Ethnicity_long, levels=c("White ethnicity", "Black ethnicity", "Other ethnicity"))
table(data1$Ethnicity, data1$Ethnicity_long)

#Shorterned ethnicity
data1$Ethnicity_short <- data1$Ethnicity
data1$Ethnicity_short <- factor(data1$Ethnicity_short, levels=c("White", "Black", "Other"))
table(data1$Ethnicity, data1$Ethnicity_short)
data1$Ethnicity <- NULL




###Subdivide data set by the 2 timepoints.
table(data1$weeks.at.visit.group)
data1[1:5, 1:5]

data_early <- data1 %>% 
  filter(weeks.at.visit.group == "10-15_weeks") 
dim(data_early)

data_late <- data1 %>% 
  filter(weeks.at.visit.group == "16-24_weeks")
dim(data_late)




#changing row names to Sample_IDs
rownames(data_early) <- data_early$Sample_ID
data_early[1:5, 1:4]
dim(data_early) #325 3180

#changing row names to Sample_IDs
rownames(data_late) <- data_late$Sample_ID
data_late[1:5, 1:4]
dim(data_late) #323 3180



###Saving the data from this script
save(data1, data_early, data_late, 
     file="data_BV/dfs_specifically_for_LEfSe.RData")


#Loading .Rdata saved 
#lnames = load(file="data_BV/dfs_specifically_for_LEfSe.RData")
#lnames

