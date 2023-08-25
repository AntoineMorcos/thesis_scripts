# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(reshape2) # heatmaps
library(BiocManager) 
library(ggrepel) # labels on ggplots
library(stats) # for Wilcoxon or other stats tests
library(gtsummary)   #nice tables

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

# Loading .Rdata from Immune_cell_types_in_larger_cohort
lnames = load(file="data/ALL_objects_in_Immune_cell_types_in_larger_cohort.RData")
lnames





###########################################################################################################################
##                                                                                                                       ##
##                                                     df prep                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(immune_uniq_all_demo)
# View(immune_uniq_all_demo)


# Filter data of only cols I want in demo table so I can check the str
for_demo <- immune_uniq_all_demo %>% 
  dplyr::mutate(Gestation.at.delivery.wks.dec = ((Gestationatdeliveryw*7) + Gestationatdeliveryd) / 7) %>%  
  #number of weeks x7 plus number of days, divided by 7 again to put it back into weeks 
  dplyr::select(Participant.ID, Outcome,
                Gestation.at.delivery.wks.dec,
                Lowriskatenrolment,
                Ageatregistration, BMI, BMI_category, Ethnicity_detailed, Smoking, 
                Didshereceiveprogesterone,
                Gest.at.cyto.wks.dec) %>% 
  dplyr::rename(Risk=Lowriskatenrolment,
                Progesterone=Didshereceiveprogesterone,
                Age=Ageatregistration)

# Look               
str(for_demo)


# Smoking
summary(for_demo$Smoking)
for_demo$Smoking <- factor(for_demo$Smoking, levels=c("Never", "Ex - gave up before pregnancy", "Ex - gave up in pregnancy", "Current"))
summary(for_demo$Smoking)


# Risk
table(for_demo$Risk, useNA="always")
for_demo$Risk[for_demo$Risk=="Yes"] <- "Low"
for_demo$Risk[for_demo$Risk=="No"] <- "High"
for_demo$Risk <- factor(for_demo$Risk, levels=c("Low", "High"))
table(for_demo$Risk, useNA="always")

# Ethnicity_detailed
table(for_demo$Ethnicity_detailed, useNA="always")
for_demo$Ethnicity_detailed[for_demo$Ethnicity_detailed=="Unclassified"] <- "Other"
for_demo$Ethnicity_detailed[for_demo$Ethnicity_detailed=="Afro-Caribbean"] <- "African Caribbean"
for_demo$Ethnicity_detailed <- factor(for_demo$Ethnicity_detailed, levels=c("European", "African", "African Caribbean", "Other")) 
table(for_demo$Ethnicity_detailed, useNA="always")


# Progesterone
table(for_demo$Progesterone, useNA="always")
for_demo$Progesterone[for_demo$Progesterone==""] <- NA
for_demo$Progesterone <- factor(for_demo$Progesterone) 
table(for_demo$Progesterone, useNA="always")


# Age - Integer to numeric
for_demo$Age 
for_demo$Age <- as.numeric(for_demo$Age)
for_demo$Age 

# Final check
str(for_demo)
summary(for_demo)
names(for_demo)

# Participant.ID
for_demo$Participant.ID <- NULL 








###########################################################################################################################
##                                                                                                                       ##
##                                                     Demo table                                                        ##
##                                                                                                                       ##
###########################################################################################################################


# Info
?gtsummary::tbl_summary # I had dependency issues with rlang, vctrs, cli, magrittr, xfun. Had to delete & reinstall most of these for gtsummary::tbl_summary() to work.
?gtsummary::add_p() 
?gt::gtsave


# Export demo table - columns by Outcome status
Table_by_Outcome <- for_demo %>% gtsummary::tbl_summary(     
  by = Outcome,                                               
  statistic = list(all_continuous() ~ "{mean} ({sd})",        
                   all_categorical() ~ "{n} / {N} ({p}%)"),   
  digits = all_continuous() ~ 1,                              
  type   = list(all_categorical() ~ "categorical"),    
  missing_text = "Missing",                             
  label  = list(
    BMI ~ "BMI (kg/m^2)",
    BMI_category ~ "BMI classification",
    Ethnicity_detailed ~ "Ethnicity",
    Gestation.at.delivery.wks.dec ~ "Gestation at delivery (weeks)",
    Gest.at.cyto.wks.dec ~ "Gestation at cytobrush sampling (weeks)")) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% 
  # Export 
  as_gt() %>%
  gt::gtsave(filename = "data/Demographics_n46_Outcome.tex")  # use extensions .html .tex .ltx .rtf


# What % were ptb
6/46 # 0.1304348


