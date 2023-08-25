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
#setwd("~/OneDrive - King's College London/PhD/Projects") #mac
setwd("C:/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()
#list.files()







###########################################################################################################################
##                                                                                                                       ##
##                                             Load in data                                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
lnames = load(file="data/immune_CytokinesComplementsIgM_2023.RData")
lnames



 




###########################################################################################################################
##                                                                                                                       ##
##                                                     df prep                                                           ##
##                                                                                                                       ##
###########################################################################################################################

#Look
names(immu_MD)
# View(immu_MD)

# Checking when these samples were taken
test <- immu_MD %>% 
  dplyr::mutate(Gestation.at.visit.wks.dec = ((ga_visit_w*7) + ga_visit_d) / 7)  #number of weeks x7 plus number of days, divided by 7 again to put it back into weeks
summary(test$Gestation.at.visit.wks.dec)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   12.14   12.71   12.86   13.43   15.86 

# Filter data so I can check the str
data_for_table <- immu_MD %>% 
  dplyr::select(Gestation.at.delivery.wks.dec, sPTB34, 
                Risk, 
                Previous.spontaneous.PTB...37.weeks., Previous.PPROM...37.weeks., Previous.late.miscarriage..16.23.6.weeks., Previous.cervical.surgery..e.g..LLETZ..Cone., Uterine.abnormality, #make up Risk
                CL_absolute_small_value_at_any_visit, Short_cervix_value_this_visit,
                Did.she.receive.progesterone.,  #intervention
                Preeclampsia, Pre.existing.hypertension, Gestational.Diabetes, Type.2.diabetes, Asthma, #medical conditions
                medscinet_ph, Diagnosed.with.BV, Past.or.present.history.of.BV, History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy, Past.or.present.history.of.GBS, #bacterial
                Past.or.present.history.of.Domestic.Violence,
                Age, BMI, BMI_category, IMD_rank, IMD_decile,
                Ethnicity, Primigravida, Gender, Smoking, Past.or.present.history.of.recreational.drug.use, #demographics
                Apgar.Score.1.min, Apgar.Score.5.min,  
                Chorioamnionitis,
                sPTB37) 
# No women with: Autoimmune.disease, Type.1.diabetes, Suspected.fetal.growth.restriction, SGA
str(data_for_table)


# Making outcome data more clear for table
table(data_for_table$sPTB37)
data_for_table$sPTB37[data_for_table$sPTB37=="Yes"] <- "Preterm"
data_for_table$sPTB37[data_for_table$sPTB37=="No"] <- "Term"
table(data_for_table$sPTB37)
data_for_table$sPTB37 <- factor(data_for_table$sPTB37, levels=c("Term", "Preterm"))



#Unknowns to NAs - Gender
table(data_for_table$Gender, useNA = "always")
data_for_table$Gender[data_for_table$Gender=="Unknown"] <- NA
table(data_for_table$Gender, useNA = "always")

#Unknowns to NAs - Past.or.present.history.of.Domestic.Violence
table(data_for_table$Past.or.present.history.of.Domestic.Violence, useNA = "always")
data_for_table$Past.or.present.history.of.Domestic.Violence[data_for_table$Past.or.present.history.of.Domestic.Violence=="Unknown"] <- NA
table(data_for_table$Past.or.present.history.of.Domestic.Violence, useNA = "always")



# Character data to factor data
str(data_for_table)
cha_names_outcomes <- names(data_for_table %>% dplyr::select(sPTB34:Uterine.abnormality, 
                                                             Did.she.receive.progesterone.:Asthma,
                                                             Diagnosed.with.BV:Past.or.present.history.of.Domestic.Violence,
                                                             Ethnicity:Past.or.present.history.of.recreational.drug.use,
                                                             Chorioamnionitis))
data_for_table[cha_names_outcomes] <- lapply(data_for_table[cha_names_outcomes], factor)
str(data_for_table)
summary(data_for_table)

# Integer to numeric
# CL_absolute_small_value_at_any_visit
data_for_table$CL_absolute_small_value_at_any_visit 
data_for_table$CL_absolute_small_value_at_any_visit <- as.numeric(data_for_table$CL_absolute_small_value_at_any_visit)
data_for_table$CL_absolute_small_value_at_any_visit 

# Short_cervix_value_this_visit
data_for_table$Short_cervix_value_this_visit 
data_for_table$Short_cervix_value_this_visit <- as.numeric(data_for_table$Short_cervix_value_this_visit) #"Not_measured" automatically to NAs
data_for_table$Short_cervix_value_this_visit 

# Age
data_for_table$Age 
data_for_table$Age <- as.numeric(data_for_table$Age)
data_for_table$Age 

# Final check
str(data_for_table)
summary(data_for_table)
names(data_for_table)








###########################################################################################################################
##                                                                                                                       ##
##                                                     Demo table                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?gtsummary::tbl_summary # I had dependency issues with rlang, vctrs, cli, magrittr, xfun. Had to delete & reinstall most of these for gtsummary::tbl_summary() to work.
?gtsummary::add_p() 
?gt::gtsave


# Export demo table - columns by sPTB37 status
Table_by_sPTB37 <- data_for_table %>% gtsummary::tbl_summary(     
    by = sPTB37,                                                # stratify entire table by this variable
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(all_categorical() ~ "categorical",            # force all categorical levels to display
                  Apgar.Score.5.min ~ "continuous",
                  IMD_decile ~ "continuous",
                  medscinet_ph ~ "continuous"),     # force the autoconvert of INT to FACTOR for this variable https://www.pipinghotdata.com/posts/2021-07-14-polished-summary-tables-in-r-with-gtsummary/           
    missing_text = "Missing",                                   # how missing values should display
    label  = list(
      Previous.spontaneous.PTB...37.weeks. ~ "Previous sPTB <37 weeks",
      Previous.PPROM...37.weeks. ~ "Previous PPROM",
      Previous.late.miscarriage..16.23.6.weeks. ~ "Previous late miscarriage",
      Previous.cervical.surgery..e.g..LLETZ..Cone. ~ "Previous cervical surgery",
      Uterine.abnormality ~ "Uterine abnormality",
      Pre.existing.hypertension ~ "Pre-existing hypertension",
      Gestational.Diabetes ~ "Gestational diabetes",
      Type.2.diabetes ~ "Type 2 diabetes",
      Past.or.present.history.of.Domestic.Violence ~ "Domestic violence",
      medscinet_ph ~ "High vaginal pH",
      Diagnosed.with.BV ~ "Diagnosed with BV",
      Past.or.present.history.of.BV ~ "History of BV",
      History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy ~ "History of UTIs in pregnancy",
      Past.or.present.history.of.GBS ~ "History of GBS infection",
      Gender ~ "Gender of baby",
      BMI ~ "BMI (kg/m^2)",
      BMI_category ~ "BMI classification",
      IMD_decile ~ "IMD decile",
      IMD_rank ~ "IMD rank",
      Past.or.present.history.of.recreational.drug.use ~ "Recreational drugs",
      Gestation.at.delivery.wks.dec ~ "Gestation at delivery (weeks)",
      sPTB34 ~ "sPTB <34 weeks",
      Apgar.Score.1.min ~ "Apgar score (1 min)",
      Apgar.Score.5.min ~ "Apgar score (5 min)",
      Did.she.receive.progesterone. ~ "Progesterone",
      CL_absolute_small_value_at_any_visit ~ "Minimum cervical length across pregnancy",
      Short_cervix_value_this_visit ~ "Cervical length at visit")) %>%
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% #fixing Fisher test error for cols with many categories https://github.com/ddsjoberg/gtsummary/issues/952
  # Export 
  as_gt() %>%
  gt::gtsave(filename = "data/Demographics_n87_sPTB37.tex")  # use extensions .html .tex .ltx .rtf


# Checking
data_for_table %>% gtsummary::tbl_summary(     
  by = sPTB37,                                                # stratify entire table by this variable
  statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                   all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
  digits = all_continuous() ~ 1,                              # rounding for continuous columns
  type   = list(all_categorical() ~ "categorical",            # force all categorical levels to display
                Apgar.Score.5.min ~ "continuous",
                IMD_decile ~ "continuous",
                medscinet_ph ~ "continuous"),     # force the autoconvert of INT to FACTOR for this variable https://www.pipinghotdata.com/posts/2021-07-14-polished-summary-tables-in-r-with-gtsummary/           
  missing_text = "Missing") %>%
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9))
59+28


# What % were ptb
28/87
