#clear environment
rm(list=ls())

getwd()
#list.files()

library(multcomp)
library(reshape2)
library(tidyverse)
library(tidyverse)
library(reshape2) # heatmaps
library(BiocManager) 
library(ggrepel) # labels on ggplots
library(stats) # for Wilcoxon or other stats tests
library(gtsummary)   #nice tables

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
getwd()
#list.files()





###########################################################################################################################
##                                                                                                                       ##
##                                                     Load data                                                         ##
##                                                                                                                       ##
###########################################################################################################################

##Loading .Rdata saved 
lnames = load(file="data_BV/main_dfs.RData") 
lnames






###########################################################################################################################
##                                                                                                                       ##
##                                                     df prep                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(data1)[1:210]
#View(data1)
table(data_early$participant_id==data_late$participant_id)
table(data_early$Term_outcome)
table(data_early$Ethnicity_long, data_early$Ethnicity_short)
table(data_early$Reported.history.of.domestic.violence)


# Looking at gestational age at delivery, checking they match
test <- data_early %>% 
  dplyr::mutate(Gestation.at.delivery.wks.dec = ((GA.del.weeks.1*7) + GA.del.days) / 7,
                Gestation.at.delivery.wks.dec2 = Gestation.at.delivery..days/7 ) %>% 
  dplyr::select(participant_id,
                Gestation.at.delivery.wks.dec, Gestation.at.delivery.wks.dec2,
                sPTB_37w, sPTB_34w)
table(test$Gestation.at.delivery.wks.dec==test$Gestation.at.delivery.wks.dec)



# Filter data of only cols I want in demo table so I can check the str
for_demo <- data_early %>% 
  dplyr::mutate(Gestation.at.delivery.wks.dec = ((GA.del.weeks.1*7) + GA.del.days) / 7) %>%  
  #number of weeks x7 plus number of days, divided by 7 again to put it back into weeks 
  dplyr::select(participant_id, sPTB_37w, 
                Gestation.at.delivery.wks.dec, sPTB_34w,
                risk.cat,
                Previous.PTB, Previous.PPROM, Previous.late.miscarriage, Previous.cervical.surgery, UterineAbnormality, 
                cl.min,
                Pre.existing.hypertension, 
                Type.1.diabetes, Type.2.diabetes,
                Asthma, Autoimmune.disease,
                Reported.history.of.domestic.violence,
                History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy,
                Past.or.present.history.of.GBS,
                Ethnicity_long,
                Age, BMI, BMI_category,  
                Primigravida,
                Smoking, Past.or.present.history.of.recreational.drug.use) %>% 
  dplyr::rename(Participant.ID=participant_id,
                Risk=risk.cat)

# Look               
str(for_demo)
head(for_demo)


# Smoking
summary(for_demo$Smoking)
for_demo$Smoking <- factor(for_demo$Smoking, levels=c("never", "ex - gave up before pregnancy", "ex - gave up in pregnancy", "current"))
summary(for_demo$Smoking)

# Ethnicity_long
table(for_demo$Ethnicity_long, useNA="always")

# Primigravida
table(for_demo$Primigravida, useNA="always")
for_demo$Primigravida[for_demo$Primigravida=="Yes"] <- "Primigravida"
for_demo$Primigravida[for_demo$Primigravida=="No"] <- "Multigravida"
table(for_demo$Primigravida, useNA="always")

# Reported.history.of.domestic.violence
table(for_demo$Reported.history.of.domestic.violence, useNA="always")
for_demo$Reported.history.of.domestic.violence[for_demo$Reported.history.of.domestic.violence=="0"] <- "No"
for_demo$Reported.history.of.domestic.violence[for_demo$Reported.history.of.domestic.violence=="1"] <- "Yes"
table(for_demo$Reported.history.of.domestic.violence, useNA="always")

# sPTB_37w
table(for_demo$sPTB_37w, useNA="always")
for_demo$sPTB_37w <- as.character(for_demo$sPTB_37w)
for_demo$sPTB_37w[for_demo$sPTB_37w=="Yes"] <- "Preterm"
for_demo$sPTB_37w[for_demo$sPTB_37w=="No"] <- "Term"
for_demo$sPTB_37w <- factor(for_demo$sPTB_37w, levels=c("Term", "Preterm"))
table(for_demo$sPTB_37w, useNA="always")



# Integer to numeric
str(for_demo)
for_demo$Age 
for_demo$Age <- as.numeric(for_demo$Age)
for_demo$Age 
for_demo$cl.min 
for_demo$cl.min <- as.numeric(for_demo$cl.min)
for_demo$cl.min 

# Changing multiple columns from character data to factor data
str(for_demo)
cha_names_outcomes <- names(for_demo %>% dplyr::select(Previous.PTB:UterineAbnormality, 
                                                       Pre.existing.hypertension:Past.or.present.history.of.GBS,
                                                       Primigravida,
                                                       Past.or.present.history.of.recreational.drug.use))
for_demo[cha_names_outcomes] <- lapply(for_demo[cha_names_outcomes], factor)
str(for_demo)

# Final check
str(for_demo)
summary(for_demo)
names(for_demo)

# Participant.ID
for_demo$Participant.ID <- NULL#as.character(for_demo$Participant.ID)








###########################################################################################################################
##                                                                                                                       ##
##                                                     Demo table                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?gtsummary::tbl_summary # I had dependency issues with rlang, vctrs, cli, magrittr, xfun. Had to delete & reinstall most of these for gtsummary::tbl_summary() to work.
?gtsummary::add_p() 
?gt::gtsave



set.seed(100)
# Export demo table - columns by Outcome status
Table_by_Outcome <- for_demo %>% gtsummary::tbl_summary(     
  by = sPTB_37w,                                                # stratify entire table by this variable
  statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                   all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
  digits = all_continuous() ~ 1,                              # rounding for continuous columns
  type   = list(all_categorical() ~ "categorical"),            # force all categorical levels to display
  missing_text = "Missing",                                   # how missing values should display
  label  = list(
    Risk ~ "PTB risk",
    Previous.PTB ~ "Previous PTB",
    Previous.PPROM ~ "Previous PPROM",
    Previous.late.miscarriage ~ "Previous late miscarriage",
    Previous.cervical.surgery ~ "Previous cervical surgery",
    UterineAbnormality ~ "Uterine abnormality",
    cl.min ~ "Minimum cervical length (mm)",
    Pre.existing.hypertension ~ "Pre-existing hypertension",
    Type.1.diabetes ~ "Type 1 diabetes",
    Type.2.diabetes ~ "Type 2 diabetes",
    Autoimmune.disease ~ "Autoimmune disease",
    Reported.history.of.domestic.violence ~ "Reported history of domestic violence",
    History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy ~ "History of UTIs in pregnancy",
    Past.or.present.history.of.GBS ~ "History of GBS infection",
    Ethnicity_long ~ "Ethnicity",
    BMI ~ "BMI (kg/m2)",
    BMI_category ~ "BMI classification",
    Past.or.present.history.of.recreational.drug.use ~ "Recreational drugs")) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% #fixing Fisher test error for cols with many categories https://github.com/ddsjoberg/gtsummary/issues/952
  # Export 
  as_gt() %>%
  gt::gtsave(filename = "data_BV/Demographics_n302_sPTB_37w.tex")  # use extensions .html .tex .ltx .rtf


# What % were ptb
46/302 # 0.1523179


