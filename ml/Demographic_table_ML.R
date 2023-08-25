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
#setwd("~/OneDrive - King's College London/PhD/Projects/ML_clinical") #mac
setwd("C:/Users/alici/OneDrive - King's College London/PhD/Projects/ML_clinical") #pc
getwd()
#list.files()







###########################################################################################################################
##                                                                                                                       ##
##                                                   Load data                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Loading Rdata from the script "2Unsupervised_learning__high_risk_women.R"
lnames <- load(file="data/2Unsupervised_learning_output_high_risk_women.RData")
lnames






###########################################################################################################################
##                                                                                                                       ##
##                                                     Demo table                                                        ##
##                                                                                                                       ##
###########################################################################################################################

#Look
names(HRClinical)
str(HRClinical)
tapply(HRClinical$Gestation.at.delivery.wks.dec, HRClinical$PTB37, summary)

# Info
?gtsummary::tbl_summary # I had dependency issues with rlang, vctrs, cli, magrittr, xfun. Had to delete & reinstall most of these for gtsummary::tbl_summary() to work.
?gtsummary::add_p() 
?gt::gtsave



# Demo table - columns by PTB37 status
Table_by_PTB37 <- HRClinical %>% 
  dplyr::select(Gestation.at.delivery.wks.dec, #sPTB37, sPTB34,
                Previous.sPTB37, Previous.PPROM, Previous.late.miscarriage, Previous.cervical.surgery, Uterine.abnormality, #make up Risk
                Early_miscarriages_number, Late_miscarriages_number, #number of miscarriages
                CL_minimum, 
                Cerclage, Progesterone, #intervention
                Preeclampsia, Pre.existing.hypertension, Gestational.Diabetes, Type.1.diabetes, Type.2.diabetes, Autoimmune.disease, Asthma, #medical conditions 
                Domestic_violence,
                Amniocentesis, #medical procedure
                pH, Diagnosed.with.BV, History.of.BV, History.of.UTIs.in.pregnancy, History.of.GBS, #bacterial
                Ethnicity1, Age, Centre, IMD_rank, IMD_decile, BMI, BMI_category, Primigravida, Gender, Smoking, Recreational_drugs, #Ethnicity2, #demographics
                Apgar.Score.1.min, Apgar.Score.5.min,  #Customised.birthweight.centiles,
                Suspected.fetal.growth.restriction, SGA, Chorioamnionitis, Major.congenital.abnormality, Neonatal.death,  
                PTB37) %>% 
  gtsummary::tbl_summary(     
    by = PTB37,                                                # stratify entire table by this variable
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(all_categorical() ~ "categorical",            # force all categorical levels to display
                  Late_miscarriages_number ~ "continuous"),     # force the autoconvert of INT to FACTOR for this variable https://www.pipinghotdata.com/posts/2021-07-14-polished-summary-tables-in-r-with-gtsummary/           
    missing_text = "Missing",                                   # how missing values should display
    label  = list(
      Previous.sPTB37 ~ "Previous sPTB <37 weeks",
      Previous.PPROM ~ "Previous PPROM",
      Previous.late.miscarriage ~ "Previous late miscarriage",
      Previous.cervical.surgery ~ "Previous cervical surgery",
      Uterine.abnormality ~ "Uterine abnormality",
      Early_miscarriages_number ~ "Number of early miscarriages",
      Late_miscarriages_number ~ "Number of late miscarriages",
      CL_minimum ~ "Minimum cervical length (mm)",
      Pre.existing.hypertension ~ "Pre-existing hypertension",
      Gestational.Diabetes ~ "Gestational diabetes",
      Type.1.diabetes ~ "Type 1 diabetes",
      Type.2.diabetes ~ "Type 2 diabetes",
      Autoimmune.disease ~ "Autoimmune disease",
      Domestic_violence ~ "Domestic violence",
      pH ~ "High vaginal pH",
      Diagnosed.with.BV ~ "Diagnosed with BV",
      History.of.BV ~ "History of BV",
      History.of.UTIs.in.pregnancy ~ "History of UTIs in pregnancy",
      History.of.GBS ~ "History of group B streptococcus (GBS) infection",
      Ethnicity1 ~ "Ethnicity",
      # Ethnicity2 ~ "Ethnicity2",
      Gender ~ "Gender of baby",
      BMI ~ "Body mass index (kg/m^2)",
      BMI_category ~ "Body mass index classification",
      IMD_decile ~ "Index of Multiple Deprivation (IMD) decile",
      Recreational_drugs ~ "Recreational drugs",
      Gestation.at.delivery.wks.dec ~ "Gestation at delivery (weeks)",
      # sPTB37 ~ "sPTB <37 weeks",
      # sPTB34 ~ "sPTB <34 weeks",
      # Customised.birthweight.centiles ~ "Customised birthweight centiles",
      Apgar.Score.1.min ~ "Apgar score (1 min)",
      Apgar.Score.5.min ~ "Apgar score (5 min)",
      Major.congenital.abnormality ~ "Major congenital abnormality",
      Neonatal.death ~ "Neonatal death",
      Suspected.fetal.growth.restriction ~ "Suspected fetal growth restriction",
      SGA ~ "Small for gestational age")) %>%                                                         
  add_p(test.args = all_tests("fisher.test") ~ list(workspace=2e9)) %>% #fixing Fisher test error for Ethnicity1 https://github.com/ddsjoberg/gtsummary/issues/952
  # Export 
  as_gt() %>%
  gt::gtsave(filename = "data/tables_exported/HR_demog_PTB37.tex")  # use extensions .html .tex .ltx .rtf




