# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(caret) #ML
library(randomForest)
library(rpart.plot) #Decision tree plots
library(e1071)
library(pROC) #ROC curve
library(ROCR)
library(knitr) #latex export

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

# Loading Rdata from the script "3Data_preprocessing_for_supervised__high_risk_women.R"
lnames <- load(file="data/3Data_preprocessing_output_high_risk_women.RData")
lnames







###########################################################################################################################
##                                                                                                                       ##
##                                                       Info                                                            ##
##                                                                                                                       ##
###########################################################################################################################

# For RF we will be using the true data, rather than the normalised data, as RF does not need normalised data

# Looking at median, as this is what is used in 1 imputation method
summary(train_HR$CL_minimum) #median 33.7mm 






###########################################################################################################################
##                                                                                                                       ##
##                                  Different variants of the features/predictors                                        ##
##                                                                                                                       ##
###########################################################################################################################

###################################### Notes on Risk

# High-risk women (10–24 weeks of gestation) were defined by the following criteria: 
#one or more of prior sPTB or late miscarriage between 16–36+6 weeks" gestation,
#previous cervical surgery, 
#uterine anomaly or 
#incidental finding of a cervical length < 25 mm on transvaginal ultrasound scan.

# Previous.sPTB37, Previous.PPROM, Previous.late.miscarriage, Previous.cervical.surgery, Short_cervix, .................

# History > Risk Factors > "Uterine Abnormality and/or other relevant history:"	& the options were:
# None
# Bicornuate uterus
# Double cervix
# Intra-uterine Septum
# Submucosal fibroids
# Other
######################################


###################################### General
# Look
dim(train_HR) #820  45
names(train_HR)

# All predictors
predi_All <- names(train_HR)
predi_All
############
########################## 

###################################### p1 = identical features to hc features (variants of the features/predictors) 
## p1 parameters:
# Including BMI,                  Ethnicity1,     Uterine.abnormality                               CL_minimum      IMD_rank
# Excluding: BMI_category,        Ethnicity2,     [features that make up "Uterine.abnormality"]     Short_cervix    IMD_decile

# p1 vector - select variables we DON'T WANT
p1_predictors <- predi_All[!predi_All %in% c("PTB37", "BMI_category", "Ethnicity2", "IMD_decile", "Short_cervix", 
                                             "Bicornuate.uterus", "Double.cervix", "Intra.uterine.Septum", "Submucosal.fibroids",
                                             "Previous.late.miscarriage",
                                             "Centre",
                                             "Cluster")]

# Looking
p1_predictors
# [1] "Age"                          "BMI"                          "pH"                           "CL_minimum"                  
# [5] "Early_miscarriages_number"    "Late_miscarriages_number"     "IMD_rank"                     "Ethnicity1"                  
# [9] "Previous.sPTB37"              "Previous.PPROM"               "Previous.cervical.surgery"    "Uterine.abnormality"         
# [13] "Cerclage"                     "Progesterone"                 "Preeclampsia"                 "Gestational.Diabetes"        
# [17] "Primigravida"                 "Smoking"                      "Diagnosed.with.BV"            "History.of.BV"               
# [21] "History.of.UTIs.in.pregnancy" "History.of.GBS"               "Gender"                       "Domestic_violence"           
# [25] "Recreational_drugs"           "Amniocentesis"                "Pre.existing.hypertension"    "Asthma"                      
# [29] "Type.1.diabetes"              "Type.2.diabetes"              "Autoimmune.disease"           "Chronic.renal.disease"       
# [33] "Chronic.viral.infection"   
length(p1_predictors) #33

# Formatting for input to model
p1 <- paste(p1_predictors, collapse = " + " ) # adds + signs between exploratory variables
PTB37p1 <- as.formula ( paste ( "PTB37" , p1, sep = " ~ " ))
######################################



###################################### p2 - same as hc features but with the addition of Cluster from hc (variants of the features/predictors)
## p2 parameters:
# Including BMI,                  Ethnicity1,     Uterine.abnormality                               CL_minimum      IMD_rank
# Excluding: BMI_category,        Ethnicity2,     [features that make up "Uterine.abnormality"]     Short_cervix    IMD_decile

# p2 vector - select variables we DON'T WANT
p2_predictors <- predi_All[!predi_All %in% c("PTB37", "BMI_category", "Ethnicity2", "IMD_decile", "Short_cervix", 
                                  "Bicornuate.uterus", "Double.cervix", "Intra.uterine.Septum", "Submucosal.fibroids",
                                  "Previous.late.miscarriage",
                                  "Centre")]
 
# Looking
p2_predictors
# [1] "Age"                          "BMI"                          "pH"                           "CL_minimum"                  
# [5] "Early_miscarriages_number"    "Late_miscarriages_number"     "IMD_rank"                     "Ethnicity1"                  
# [9] "Previous.sPTB37"              "Previous.PPROM"               "Previous.cervical.surgery"    "Uterine.abnormality"         
# [13] "Cerclage"                     "Progesterone"                 "Preeclampsia"                 "Gestational.Diabetes"        
# [17] "Primigravida"                 "Smoking"                      "Diagnosed.with.BV"            "History.of.BV"               
# [21] "History.of.UTIs.in.pregnancy" "History.of.GBS"               "Gender"                       "Domestic_violence"           
# [25] "Recreational_drugs"           "Amniocentesis"                "Pre.existing.hypertension"    "Asthma"                      
# [29] "Type.1.diabetes"              "Type.2.diabetes"              "Autoimmune.disease"           "Chronic.renal.disease"       
# [33] "Chronic.viral.infection"      "Cluster" 
length(p2_predictors) #34

# Formatting for input to model
p2 <- paste(p2_predictors, collapse = " + " ) # adds + signs between exploratory variables
PTB37p2 <- as.formula ( paste ( "PTB37" , p2, sep = " ~ " ))
######################################


###################################### p3 - same as hc features but without CL_minimum (variants of the features/predictors)
## p3 parameters:
# Including BMI,                 Ethnicity1,     Uterine.abnormality                                 IMD_rank
# Excluding: BMI_category,       Ethnicity2,     [features that make up "Uterine.abnormality"]       IMD_decile

# p3 vector - select variables we DON'T WANT
p3_predictors <- predi_All[!predi_All %in% c("PTB37", "BMI_category", "Ethnicity2", "IMD_decile", "Short_cervix", 
                                             "Bicornuate.uterus", "Double.cervix", "Intra.uterine.Septum", "Submucosal.fibroids",
                                             "Previous.late.miscarriage",
                                             "Centre",
                                             "Cluster",
                                             "CL_minimum")]

# Looking
p3_predictors
#  [1] "Age"                          "BMI"                          "pH"                           "Early_miscarriages_number"   
# [5] "Late_miscarriages_number"     "IMD_rank"                     "Ethnicity1"                   "Previous.sPTB37"             
# [9] "Previous.PPROM"               "Previous.cervical.surgery"    "Uterine.abnormality"          "Cerclage"                    
# [13] "Progesterone"                 "Preeclampsia"                 "Gestational.Diabetes"         "Primigravida"                
# [17] "Smoking"                      "Diagnosed.with.BV"            "History.of.BV"                "History.of.UTIs.in.pregnancy"
# [21] "History.of.GBS"               "Gender"                       "Domestic_violence"            "Recreational_drugs"          
# [25] "Amniocentesis"                "Pre.existing.hypertension"    "Asthma"                       "Type.1.diabetes"             
# [29] "Type.2.diabetes"              "Autoimmune.disease"           "Chronic.renal.disease"        "Chronic.viral.infection"     
length(p3_predictors) #32

# Formatting for input to model
p3 <- paste(p3_predictors, collapse = " + " ) # adds + signs between exploratory variables
PTB37p3 <- as.formula ( paste ( "PTB37" , p3, sep = " ~ " ))
######################################


###################################### p4 - only features definitely available at first antenatal visit (variants of the features/predictors)
## p4 parameters:
# Including BMI,                 Ethnicity1,     Uterine.abnormality                                 IMD_rank
# Excluding: BMI_category,       Ethnicity2,     [features that make up "Uterine.abnormality"]       IMD_decile  ... et al. ()


# Predictor variant 4: all features in hc excluding cervical length, and data on cerclage, progesterone, 
#         gestational diabetes, preeclampsia, gender of fetus, amniocentesis, high vaginal pH, and incidental diagnosis of BV

# p4 vector - select variables we DON'T WANT
p4_predictors <- predi_All[!predi_All %in% c("PTB37", "BMI_category", "Ethnicity2", "IMD_decile", "Short_cervix", 
                                             "Bicornuate.uterus", "Double.cervix", "Intra.uterine.Septum", "Submucosal.fibroids",
                                             "Previous.late.miscarriage",
                                             "Centre",
                                             "Cluster",
                                             "CL_minimum",
                                             "Cerclage", "Progesterone", "Gestational.Diabetes", "Preeclampsia", "Gender",
                                             "Amniocentesis", "pH", "Diagnosed.with.BV")]

# Looking
p4_predictors
# [1] "Age"                          "BMI"                          "Early_miscarriages_number"    "Late_miscarriages_number"    
# [5] "IMD_rank"                     "Ethnicity1"                   "Previous.sPTB37"              "Previous.PPROM"              
# [9] "Previous.cervical.surgery"    "Uterine.abnormality"          "Primigravida"                 "Smoking"                     
# [13] "History.of.BV"                "History.of.UTIs.in.pregnancy" "History.of.GBS"               "Domestic_violence"           
# [17] "Recreational_drugs"           "Pre.existing.hypertension"    "Asthma"                       "Type.1.diabetes"             
# [21] "Type.2.diabetes"              "Autoimmune.disease"           "Chronic.renal.disease"        "Chronic.viral.infection"  
length(p4_predictors) #24

# Formatting for input to model
p4 <- paste(p4_predictors, collapse = " + " ) # adds + signs between exploratory variables
PTB37p4 <- as.formula ( paste ( "PTB37" , p4, sep = " ~ " ))
######################################







###########################################################################################################################
##                                                                                                                       ##
##                                                RF imputation                                                          ##
##                                                                                                                       ##
###########################################################################################################################

# RF imputation on missing data from https://github.com/StatQuest/random_forest_demo/blob/master/random_forest_demo.R "Breiman says 4 to 6 iterations is usually good enough." 

# Info
?rfImpute

######################################## p1
# Impute any missing values in the training set using proximities
set.seed(100)
train_HR_RFimputed_p1 <- rfImpute(PTB37p1, data = train_HR, iter=6, ntree=500) #6 iterations #chose best b/ 4-6

# Look
summary(train_HR_RFimputed_p1)
summary(train_HR)
########################################

######################################## p2
# Impute any missing values in the training set using proximities
set.seed(100)
train_HR_RFimputed_p2 <- rfImpute(PTB37p2, data = train_HR, iter=4, ntree=500) #4 iterations #chose best b/ 4-6

# Look
summary(train_HR_RFimputed_p2)
summary(train_HR)
########################################

######################################## p3
# Impute any missing values in the training set using proximities
set.seed(100)
train_HR_RFimputed_p3 <- rfImpute(PTB37p3, data = train_HR, iter=6, ntree=500) #6 iterations #chose best b/ 4-6

# Look
summary(train_HR_RFimputed_p3)
summary(train_HR)
########################################

######################################## p4
# Impute any missing values in the training set using proximities
set.seed(100)
train_HR_RFimputed_p4 <- rfImpute(PTB37p4, data = train_HR, iter=4, ntree=500) #6 iterations #chose best b/ 4-6

# Look
summary(train_HR_RFimputed_p4)
summary(train_HR)
########################################


# Export this df (as the imputation doesn't seem to always produce a consistent result after opening & closing R)
save(train_HR_RFimputed_p1, train_HR_RFimputed_p2, train_HR_RFimputed_p3, train_HR,
     file=paste0("data/training_data_with_RF_imputation_", as.character(Sys.Date()), ".RData", sep=""))

# Export this df (as the imputation doesn't seem to always produce a consistent result after opening & closing R)
save(train_HR_RFimputed_p4, 
     file=paste0("data/training_data_with_RF_imputation_p4_", as.character(Sys.Date()), ".RData", sep=""))


# # Loading imputed data from the script "2Unsupervised_learning__high_risk_women.R"
# lnames <- load(file="data/training_data_with_RF_imputation_2022-12-20.RData")
# lnames








###########################################################################################################################
##                                                                                                                       ##
##                                                 Other RF prep                                                         ##
##                                                                                                                       ##
###########################################################################################################################

######################################## Looking 
# Info
?randomForest

## Working out which axis is which on the confusion matrix
table(train_HR$PTB37)
# Term Preterm 
#  663     157

########################################



# It appears there is no need for cross-validation in Random Forests. Sources:
# https://stackoverflow.com/questions/19760169/how-to-perform-random-forest-cross-validation-in-r
# https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#ooberr



######################################## dfs needed for tuneRF()

# Prep for tune for RFp1imp (median/mode imputation)
train_HR_PTB37p1 <- train_HR %>% 
  dplyr::select(PTB37, all_of(p1_predictors))
names(train_HR_PTB37p1)

# Prep for tune for RFp1impRf (RF imputation)
train_HR_RFimputed_PTB37p1 <- train_HR_RFimputed_p1 %>% 
  dplyr::select(PTB37, all_of(p1_predictors))
names(train_HR_RFimputed_PTB37p1)

# Prep for tune for RFp2imp (median/mode imputation)
train_HR_PTB37p2 <- train_HR %>% 
  dplyr::select(PTB37, all_of(p2_predictors))
names(train_HR_PTB37p2)

# Prep for tune for RFp2impRf (RF imputation)
train_HR_RFimputed_PTB37p2 <- train_HR_RFimputed_p2 %>% 
  dplyr::select(PTB37, all_of(p2_predictors))
names(train_HR_RFimputed_PTB37p2)

# Prep for tune for RFp3imp (median/mode imputation)
train_HR_PTB37p3 <- train_HR %>% 
  dplyr::select(PTB37, all_of(p3_predictors))
names(train_HR_PTB37p3)

# Prep for tune for RFp3impRf (RF imputation)
train_HR_RFimputed_PTB37p3 <- train_HR_RFimputed_p3 %>% 
  dplyr::select(PTB37, all_of(p3_predictors))
names(train_HR_RFimputed_PTB37p3)

# Prep for tune for RFp4imp (median/mode imputation)
train_HR_PTB37p4 <- train_HR %>% 
  dplyr::select(PTB37, all_of(p4_predictors))
names(train_HR_PTB37p4)

# Prep for tune for RFp4impRf (RF imputation)
train_HR_RFimputed_PTB37p4 <- train_HR_RFimputed_p4 %>% 
  dplyr::select(PTB37, all_of(p4_predictors))
names(train_HR_RFimputed_PTB37p4)
########################################






###########################################################################################################################
##                                                                                                                       ##
##                                          Running Random Forest models                                                 ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?randomForest
?na.roughfix #Impute Missing Values by median/mode.
?rfImpute #https://stat.ethz.ch/pipermail/r-help/2003-August/037278.html
?tuneRF
?caret::train





######################################## RF when imputing NAs to median/mode - p1 (i.e. ignoring Cluster)
# Imputing NAs to median/mode - p1 
set.seed(100)
RFp1imp <- randomForest(PTB37p1, data=train_HR,  ntree=500, mtry=4, importance=T, proximity=T, na.action = na.roughfix) #imputing
RFp1imp
# OPTIMAL mtry = 4
# OOB estimate of  error rate: 19.63%
# Confusion matrix:
#         Term Preterm class.error
# Term     654       9  0.01357466
# Preterm  152       5  0.96815287
########################################



######################################## RF when imputing NAs using RF - p1 (i.e. ignoring Cluster)
# Imputing NAs using RF - p1 
set.seed(100)
RFp1impRf <- randomForest(PTB37p1, data=train_HR_RFimputed_p1, mtry=4, ntree=500, importance=T, proximity=T) 
RFp1impRf
# OPTIMAL mtry = 4
#OOB estimate of  error rate: 20%
# Confusion matrix:
#         Term Preterm class.error
# Term     650      13  0.01960784
# Preterm  151       6  0.96178344
########################################



######################################## RF when imputing NAs to median/mode - p2
# Imputing NAs to median/mode - p2
set.seed(100)
RFp2imp <- randomForest(PTB37p2, data=train_HR,  ntree=500, mtry=4, importance=T, proximity=T, na.action = na.roughfix) #imputing
RFp2imp
# OPTIMAL mtry = 4
# OOB estimate of  error rate: 19.88%
# Confusion matrix:
#          Term Preterm class.error
# Term     652      11  0.01659125
# Preterm  152       5  0.96815287
# ########################################



######################################## RF when imputing NAs using RF - p2
# Imputing NAs using RF - p2
set.seed(100)
RFp2impRf <- randomForest(PTB37p2, data=train_HR_RFimputed_p2, mtry=4,  ntree=500, importance=T, proximity=T) 
RFp2impRf
# OPTIMAL mtry = 4
#         OOB estimate of  error rate: 19.39%
# Confusion matrix:
#          Term Preterm class.error
# Term     654       9  0.01357466
# Preterm  150       7  0.95541401
#######################################



######################################## RF when imputing NAs to median/mode - p3 (i.e. ignoring CL)
# Imputing NAs to median/mode - p3
set.seed(100)
RFp3imp <- randomForest(PTB37p3, data=train_HR,  ntree=500, mtry=4, importance=T, proximity=T, na.action = na.roughfix) #imputing
RFp3imp
# OPTIMAL mtry = 4
# OOB estimate of  error rate: 19.39%
# Confusion matrix:
#   Term Preterm class.error
# Term     659       4 0.006033183
# Preterm  155       2 0.987261146
########################################



######################################## RF when imputing NAs using RF - p3 (i.e. ignoring CL)
# Imputing NAs using RF - p3 
set.seed(100)
RFp3impRf <- randomForest(PTB37p3, data=train_HR_RFimputed_p3, mtry=6, ntree=500, importance=T, proximity=T) 
RFp3impRf
# OPTIMAL mtry = 6
# OOB estimate of  error rate: 19.88%
# Confusion matrix:
#         Term Preterm class.error
# Term     653      10  0.01508296
# Preterm  153       4  0.97452229
########################################




######################################## RF when imputing NAs to median/mode - p4 (i.e. data at 1st visit)
# Imputing NAs to median/mode - p4
set.seed(100)
RFp4imp <- randomForest(PTB37p4, data=train_HR,  ntree=500, mtry=6, importance=T, proximity=T, na.action = na.roughfix) #imputing
RFp4imp
# OPTIMAL mtry = 6
# OOB estimate of  error rate: 19.63%
# Confusion matrix:
#   Term Preterm class.error
# Term     653      10  0.01508296
# Preterm  151       6  0.96178344
########################################



######################################## RF when imputing NAs using RF - p4 (i.e. data at 1st visit)
# Imputing NAs using RF - p4 
set.seed(100)
RFp4impRf <- randomForest(PTB37p4, data=train_HR_RFimputed_p4, mtry=6, ntree=500, importance=T, proximity=T) 
RFp4impRf
# OPTIMAL mtry = 6. 
#  OOB estimate of  error rate: 19.39%
# Confusion matrix:
#   Term Preterm class.error
# Term     654       9  0.01357466
# Preterm  150       7  0.95541401
########################################



######################################## Export RF models
# Export this df (as the imputation doesn't seem to always produce a consistent result after opening & closing R)
save(RFp1imp, RFp1impRf, RFp2imp, RFp2impRf, RFp3imp, RFp3impRf, 
     file=paste0("data/RF_models_", as.character(Sys.Date()), ".RData", sep=""))

save(RFp4imp, RFp4impRf, 
     file=paste0("data/RF_models_p4_", as.character(Sys.Date()), ".RData", sep=""))

# # Loading imputed data from the script "4Supervised_learning__high_risk_women.R"
# lnames <- load(file="data/RF_models_2022-12-20.RData")
# lnames
########################################







###########################################################################################################################
##                                                                                                                       ##
##                                                     ROC curves                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?roc
?auc
citation("pROC")



# RFp1imp - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp1imp.pdf", width=4, height=4)
ROC_RFp1imp_train <- roc(train_HR$PTB37, RFp1imp$votes[,2])
plot(ROC_RFp1imp_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp1imp_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp1imp_train)
dev.off()

# RFp1impRf - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp1impRf.pdf", width=4, height=4)
ROC_RFp1impRf_train <- roc(train_HR$PTB37, RFp1impRf$votes[,2])
plot(ROC_RFp1impRf_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp1impRf_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp1impRf_train)
dev.off()

# RFp2imp - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp2imp.pdf", width=4, height=4)
ROC_RFp2imp_train <- roc(train_HR$PTB37, RFp2imp$votes[,2])
plot(ROC_RFp2imp_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp2imp_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp2imp_train)
dev.off()

# RFp2impRf - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp2impRf.pdf", width=4, height=4)
ROC_RFp2impRf_train <- roc(train_HR$PTB37, RFp2impRf$votes[,2])
plot(ROC_RFp2impRf_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp2impRf_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp2impRf_train)
dev.off()

# RFp3imp - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp3imp.pdf", width=4, height=4)
ROC_RFp3imp_train <- roc(train_HR$PTB37, RFp3imp$votes[,2])
plot(ROC_RFp3imp_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp3imp_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp3imp_train)
dev.off()

# RFp3impRf - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp3impRf.pdf", width=4, height=4)
ROC_RFp3impRf_train <- roc(train_HR$PTB37, RFp3impRf$votes[,2])
plot(ROC_RFp3impRf_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp3impRf_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp3impRf_train)
dev.off()

# RFp4imp - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp4imp.pdf", width=4, height=4)
ROC_RFp4imp_train <- roc(train_HR$PTB37, RFp4imp$votes[,2])
plot(ROC_RFp4imp_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp4imp_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp4imp_train)
dev.off()

# RFp4impRf - ROC curve on training set
pdf("plots/high_risk_women/Supervised/ROC/ROC_train_RFp4impRf.pdf", width=4, height=4)
ROC_RFp4impRf_train <- roc(train_HR$PTB37, RFp4impRf$votes[,2])
plot(ROC_RFp4impRf_train)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_RFp4impRf_train),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_RFp4impRf_train)
dev.off()




###########################################################################################################################
##                                                                                                                       ##
##                                              Variable importance plots                                                ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?varImpPlot
?randomForest::importance 
#type	 = either 1 or 2, specifying the type of importance measure (1=mean decrease in accuracy, 2=mean decrease in node impurity).

######################################## Plots of what variables are most important in the RFs

# Variable Importance Plot - RFp1imp
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p1imp.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp1imp, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed with median/mode)")
varImp_p1imp <- randomForest::importance(RFp1imp, type=2, scale=F)
str(varImp_p1imp)
rownames(varImp_p1imp) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p1imp) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p1imp) <- gsub(".", " ", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p1imp) <- gsub("_", " ", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p1imp) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p1imp)
dotchart(tail(sort(varImp_p1imp[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp1impRf
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p1impRf.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp1impRf, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed by RF)")
varImp_p1impRf <- randomForest::importance(RFp1impRf, type=2, scale=F)
str(varImp_p1impRf)
rownames(varImp_p1impRf) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p1impRf), fixed=T)
rownames(varImp_p1impRf) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p1impRf), fixed=T)
rownames(varImp_p1impRf) <- gsub(".", " ", rownames(varImp_p1impRf), fixed=T)
rownames(varImp_p1impRf) <- gsub("_", " ", rownames(varImp_p1impRf), fixed=T)
rownames(varImp_p1impRf) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p1impRf), fixed=T)

rownames(varImp_p1impRf)
dotchart(tail(sort(varImp_p1impRf[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp2imp
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p2imp.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp2imp, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed with median/mode)")
varImp_p2imp <- randomForest::importance(RFp2imp, type=2, scale=F)
str(varImp_p2imp)
rownames(varImp_p2imp) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p2imp), fixed=T)
rownames(varImp_p2imp) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p2imp), fixed=T)
rownames(varImp_p2imp) <- gsub(".", " ", rownames(varImp_p2imp), fixed=T)
rownames(varImp_p2imp) <- gsub("_", " ", rownames(varImp_p2imp), fixed=T)
rownames(varImp_p1imp) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p1imp), fixed=T)
rownames(varImp_p2imp)
dotchart(tail(sort(varImp_p2imp[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp2impRf
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p2impRf.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp2impRf, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed by RF)")
varImp_p2impRf <- randomForest::importance(RFp2impRf, type=2, scale=F)
str(varImp_p2impRf)
rownames(varImp_p2impRf) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p2impRf), fixed=T)
rownames(varImp_p2impRf) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p2impRf), fixed=T)
rownames(varImp_p2impRf) <- gsub(".", " ", rownames(varImp_p2impRf), fixed=T)
rownames(varImp_p2impRf) <- gsub("_", " ", rownames(varImp_p2impRf), fixed=T)
rownames(varImp_p1impRf) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p1impRf), fixed=T)
rownames(varImp_p2impRf)
dotchart(tail(sort(varImp_p2impRf[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp3imp
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p3imp.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp3imp, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed with median/mode)")
varImp_p3imp <- randomForest::importance(RFp3imp, type=2, scale=F)
str(varImp_p3imp)
rownames(varImp_p3imp) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p3imp), fixed=T)
rownames(varImp_p3imp) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p3imp), fixed=T)
rownames(varImp_p3imp) <- gsub(".", " ", rownames(varImp_p3imp), fixed=T)
rownames(varImp_p3imp) <- gsub("_", " ", rownames(varImp_p3imp), fixed=T)
rownames(varImp_p3imp) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p3imp), fixed=T)
rownames(varImp_p3imp)
dotchart(tail(sort(varImp_p3imp[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp3impRf
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p3impRf.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp3impRf, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed by RF)")
varImp_p3impRf <- randomForest::importance(RFp3impRf, type=2, scale=F)
str(varImp_p3impRf)
rownames(varImp_p3impRf) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p3impRf), fixed=T)
rownames(varImp_p3impRf) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p3impRf), fixed=T)
rownames(varImp_p3impRf) <- gsub(".", " ", rownames(varImp_p3impRf), fixed=T)
rownames(varImp_p3impRf) <- gsub("_", " ", rownames(varImp_p3impRf), fixed=T)
rownames(varImp_p3impRf) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p3impRf), fixed=T)
rownames(varImp_p3impRf)
dotchart(tail(sort(varImp_p3impRf[,1]), n=10), xlim=c(0,35), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp4imp
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p4imp.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp4imp, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed with median/mode)")
varImp_p4imp <- randomForest::importance(RFp4imp, type=2, scale=F)
str(varImp_p4imp)
rownames(varImp_p4imp) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p4imp), fixed=T)
rownames(varImp_p4imp) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p4imp), fixed=T)
rownames(varImp_p4imp) <- gsub(".", " ", rownames(varImp_p4imp), fixed=T)
rownames(varImp_p4imp) <- gsub("_", " ", rownames(varImp_p4imp), fixed=T)
rownames(varImp_p4imp) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p4imp), fixed=T)
rownames(varImp_p4imp)
dotchart(tail(sort(varImp_p4imp[,1]), n=10), xlim=c(0,45), xlab="Mean decrease in gini")
dev.off()

# Variable Importance Plot - RFp4impRf
pdf("plots/high_risk_women/Supervised/variable_importance/varImp_p4impRf.pdf", width=4, height=4)
par(mar = c(4.1, 0.5, 0.5, 0.5)) #bottom, left, top, right)
# varImpPlot(RFp4impRf, sort = T, n.var= 10, main=NULL, type=2)
# title("Variable Importance\n(missing data imputed by RF)")
varImp_p4impRf <- randomForest::importance(RFp4impRf, type=2, scale=F)
str(varImp_p4impRf)
rownames(varImp_p4impRf) <- gsub("Ethnicity1", "Ethnicity", rownames(varImp_p4impRf), fixed=T)
rownames(varImp_p4impRf) <- gsub("CL_minimum", "Cervical length", rownames(varImp_p4impRf), fixed=T)
rownames(varImp_p4impRf) <- gsub(".", " ", rownames(varImp_p4impRf), fixed=T)
rownames(varImp_p4impRf) <- gsub("_", " ", rownames(varImp_p4impRf), fixed=T)
rownames(varImp_p4impRf) <- gsub("Previous sPTB37", "Previous sPTB", rownames(varImp_p4impRf), fixed=T)
rownames(varImp_p4impRf)
dotchart(tail(sort(varImp_p4impRf[,1]), n=10), xlim=c(0,45), xlab="Mean decrease in gini")
dev.off()

######################################## 


# ######################################## Variable Importance TABLE of RFp3imp
# # Variable Importance Table of all predictors - RFp3imp
# VarImp_RFp3imp_Accuracy <- data.frame(importance(RFp3imp, type=1)) #1 for model accuracy (MeanDecreaseAccuracy)
# VarImp_RFp3imp_Gini <- data.frame(importance(RFp3imp, type=2)) #2 for node impurity (MeanDecreaseGini)
# VarImp_RFp3imp <- merge(VarImp_RFp3imp_Gini, VarImp_RFp3imp_Accuracy, by="row.names", all=T) #all 32 predictors
# VarImp_RFp3imp <- VarImp_RFp3imp %>% dplyr::rename("Predictors"="Row.names") #rename col
# VarImp_RFp3imp <- VarImp_RFp3imp[order(VarImp_RFp3imp$MeanDecreaseGini,decreasing = T),] #reorder df
# VarImp_RFp3imp 
# ######################################## 












###########################################################################################################################
##                                                                                                                       ##
##                                           Optimise probability threshold                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?caret::confusionMatrix


######################################################################### RFp1imp
################################ Confusion matrices for RFp1imp
# Generate data frame of predictions
pred_RFp1imp <- data.frame(predict(RFp1imp, type="prob"),
                          actual=train_HR$PTB37,
                          thresh0.5=predict(RFp1imp))
head(pred_RFp1imp)

# Add prediction if we set probability threshold of 0.3 (instead of 0.5) for classifying a prediction as "Preterm"
pred_RFp1imp$thresh0.3 <- factor(ifelse(pred_RFp1imp$Preterm > 0.3, "Preterm", "Term"),
                                 levels=c("Term", "Preterm"))
pred_RFp1imp$thresh0.2 <- factor(ifelse(pred_RFp1imp$Preterm > 0.2, "Preterm", "Term"),
                                 levels=c("Term", "Preterm"))
pred_RFp1imp$thresh0.15 <- factor(ifelse(pred_RFp1imp$Preterm > 0.15, "Preterm", "Term"),
                                 levels=c("Term", "Preterm"))
pred_RFp1imp$thresh0.1 <- factor(ifelse(pred_RFp1imp$Preterm > 0.1, "Preterm", "Term"),
                                 levels=c("Term", "Preterm"))
head(pred_RFp1imp)
str(pred_RFp1imp)


# Look at confusion matrix probability threshold of 0.5 (standard threshold)
confusionMatrix(data=pred_RFp1imp$thresh0.5, reference=pred_RFp1imp$actual, positive="Preterm")

# Look at confusion matrix probability threshold of 0.3
confusionMatrix(data=pred_RFp1imp$thresh0.3, reference=pred_RFp1imp$actual, positive="Preterm")

# Look at confusion matrix probability threshold of 0.2
confusionMatrix(data=pred_RFp1imp$thresh0.2, reference=pred_RFp1imp$actual, positive="Preterm")
#            Reference
# Prediction Term Preterm
#    Term     445      62
#    Preterm  218      95
# ....

# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_RFp1imp$thresh0.15, reference=pred_RFp1imp$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     351      35
#    Preterm  312     122

# Look at confusion matrix probability threshold of 0.1
confusionMatrix(data=pred_RFp1imp$thresh0.1, reference=pred_RFp1imp$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     241      16
#    Preterm  422     141
# ....
################################

################################ Metrics for different probability threshold for RFp1imp
# Find best probability threshold RFp1imp
prob_RFp1imp <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_RFp1imp$Preterm >= Threshold, "Preterm", "Term"))

  conf <- confusionMatrix(data=predicted_response, reference=pred_RFp1imp$actual, positive="Preterm")

  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293

  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1imp) <- NULL #otherwise kable exports rownames
prob_RFp1imp # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1imp %>%
  knitr::kable(caption = "Metrics when optimising rf with P1 with median/mode imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 
# write_csv(prob_RFp1imp, "data/RF_metrics/HR_train_prob_RFp1imp.csv")

################################
#########################################################################





######################################################################### RFp1impRf
################################ Confusion matrices for RFp1impRf
# Generate data frame of predictions
pred_RFp1impRf <- data.frame(predict(RFp1impRf, type="prob"),
                           actual=train_HR$PTB37,
                           thresh0.5=predict(RFp1impRf))
head(pred_RFp1impRf)

# Add prediction if we set probability threshold of 0.3 (instead of 0.5) for classifying a prediction as "Preterm"
pred_RFp1impRf$thresh0.3 <- factor(ifelse(pred_RFp1impRf$Preterm > 0.3, "Preterm", "Term"),
                                   levels=c("Term", "Preterm"))
pred_RFp1impRf$thresh0.2 <- factor(ifelse(pred_RFp1impRf$Preterm > 0.2, "Preterm", "Term"),
                                   levels=c("Term", "Preterm"))
pred_RFp1impRf$thresh0.15 <- factor(ifelse(pred_RFp1impRf$Preterm > 0.15, "Preterm", "Term"),
                                   levels=c("Term", "Preterm"))
pred_RFp1impRf$thresh0.1 <- factor(ifelse(pred_RFp1impRf$Preterm > 0.1, "Preterm", "Term"),
                                   levels=c("Term", "Preterm"))
head(pred_RFp1impRf)
str(pred_RFp1impRf)


# Look at confusion matrix probability threshold of 0.5 (standard threshold)
confusionMatrix(data=pred_RFp1impRf$thresh0.5, reference=pred_RFp1impRf$actual, positive="Preterm")

# Look at confusion matrix probability threshold of 0.3
confusionMatrix(data=pred_RFp1impRf$thresh0.3, reference=pred_RFp1impRf$actual, positive="Preterm")

# Look at confusion matrix probability threshold of 0.2
confusionMatrix(data=pred_RFp1impRf$thresh0.2, reference=pred_RFp1impRf$actual, positive="Preterm")

# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_RFp1impRf$thresh0.15, reference=pred_RFp1impRf$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     359      38
#    Preterm  304     119

# Look at confusion matrix probability threshold of 0.1
confusionMatrix(data=pred_RFp1impRf$thresh0.1, reference=pred_RFp1impRf$actual, positive="Preterm")
################################

################################ Metrics for different probability threshold for RFp1impRf
# Find best probability threshold RFp1impRf
prob_RFp1impRf <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_RFp1impRf$Preterm >= Threshold, "Preterm", "Term"))
  
  conf <- confusionMatrix(data=predicted_response, reference=pred_RFp1impRf$actual, positive="Preterm")
  
  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293
  
  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1impRf) <- NULL #otherwise kable exports rownames
prob_RFp1impRf # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1impRf %>%
  knitr::kable(caption = "Metrics when optimising rf with P1 with RF imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 

################################
#########################################################################








###########################################################################################################################
##                                                                                                                       ##
##                                               Test set - RF imputation                                                ##
##                                                                                                                       ##
###########################################################################################################################

# Look
dim(test_HR) #568  47


######################################## Imputing test data with RF
# Info
?rfImpute

################# p1
# Impute any missing values in the testing set using proximities
set.seed(100)
test_HR_teRF_p1 <- rfImpute(PTB37p1, data = test_HR, iter=6, ntree=500) #6 iterations #chose best b/ 4-6

# Look
summary(test_HR_teRF_p1)
summary(test_HR)
#################

################# p2
# Impute any missing values in the testing set using proximities
set.seed(100)
test_HR_teRF_p2 <- rfImpute(PTB37p2, data = test_HR, iter=5, ntree=500) #5 iterations #chose best b/ 4-6

# Look
summary(test_HR_teRF_p2)
summary(test_HR)
#################

################# p3
# Impute any missing values in the testing set using proximities
set.seed(100)
test_HR_teRF_p3 <- rfImpute(PTB37p3, data = test_HR, iter=5, ntree=500) #5 iterations #chose best b/ 4-6

# Look
summary(test_HR_teRF_p3)
summary(test_HR)
#################

################# p4
# Impute any missing values in the testing set using proximities
set.seed(100)
test_HR_teRF_p4 <- rfImpute(PTB37p4, data = test_HR, iter=6, ntree=500) #6 iterations #chose best b/ 4-6

# Look
summary(test_HR_teRF_p4)
summary(test_HR)
#################
########################################

######################################## Imputing test data with median/mode
set.seed(100)
test_HR_teMM <- na.roughfix(test_HR)
summary(test_HR_teMM)
######################################## 






###########################################################################################################################
##                                                                                                                       ##
##                                               Test set - ROCs                                                         ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?roc

######################################## ROC curve for RFp1imp on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp1imp <- data.frame(predict(RFp1imp, newdata=test_HR_teMM, type="prob"),
                                actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp1imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp1imp.pdf", width=4, height=4)
ROC_testMM_RFp1imp <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp1imp$Preterm)
plot(ROC_testMM_RFp1imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp1imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp1imp)
dev.off()
########################################

######################################## ROC curve for RFp2imp on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp2imp <- data.frame(predict(RFp2imp, newdata=test_HR_teMM, type="prob"),
                                  actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp2imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp2imp.pdf", width=4, height=4)
ROC_testMM_RFp2imp <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp2imp$Preterm)
plot(ROC_testMM_RFp2imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp2imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp2imp)
dev.off()
########################################

######################################## ROC curve for RFp3imp on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp3imp <- data.frame(predict(RFp3imp, newdata=test_HR_teMM, type="prob"),
                                  actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp3imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp3imp.pdf", width=4, height=4)
ROC_testMM_RFp3imp <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp3imp$Preterm)
plot(ROC_testMM_RFp3imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp3imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp3imp)
dev.off()
########################################

######################################## ROC curve for RFp4imp on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp4imp <- data.frame(predict(RFp4imp, newdata=test_HR_teMM, type="prob"),
                                  actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp4imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp4imp.pdf", width=4, height=4)
ROC_testMM_RFp4imp <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp4imp$Preterm)
plot(ROC_testMM_RFp4imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp4imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp4imp)
dev.off()
########################################


######################################## ROC curve for RFp1imp on test_HR_teRF_p1 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp1imp <- data.frame(predict(RFp1imp, newdata=test_HR_teRF_p1, type="prob"),
                                  actual=test_HR_teRF_p1$PTB)

# Look
head(test_HR_teRF_p1 %>% dplyr::select(PTB37))
head(pred_testRF_RFp1imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp1imp.pdf", width=4, height=4)
ROC_testRF_RFp1imp <- roc(response=test_HR_teRF_p1$PTB, predictor=pred_testRF_RFp1imp$Preterm)
plot(ROC_testRF_RFp1imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp1imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp1imp)
dev.off()
########################################

######################################## ROC curve for RFp2imp on test_HR_teRF_p2 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp2imp <- data.frame(predict(RFp2imp, newdata=test_HR_teRF_p2, type="prob"),
                                  actual=test_HR_teRF_p2$PTB)

# Look
head(test_HR_teRF_p2 %>% dplyr::select(PTB37))
head(pred_testRF_RFp2imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp2imp.pdf", width=4, height=4)
ROC_testRF_RFp2imp <- roc(response=test_HR_teRF_p2$PTB, predictor=pred_testRF_RFp2imp$Preterm)
plot(ROC_testRF_RFp2imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp2imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp2imp)
dev.off()
########################################

######################################## ROC curve for RFp3imp on test_HR_teRF_p3 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp3imp <- data.frame(predict(RFp3imp, newdata=test_HR_teRF_p3, type="prob"),
                                  actual=test_HR_teRF_p3$PTB)

# Look
head(test_HR_teRF_p3 %>% dplyr::select(PTB37))
head(pred_testRF_RFp3imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp3imp.pdf", width=4, height=4)
ROC_testRF_RFp3imp <- roc(response=test_HR_teRF_p3$PTB, predictor=pred_testRF_RFp3imp$Preterm)
plot(ROC_testRF_RFp3imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp3imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp3imp)
dev.off()
########################################

######################################## ROC curve for RFp4imp on test_HR_teRF_p4 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp4imp <- data.frame(predict(RFp4imp, newdata=test_HR_teRF_p4, type="prob"),
                                  actual=test_HR_teRF_p4$PTB)

# Look
head(test_HR_teRF_p4 %>% dplyr::select(PTB37))
head(pred_testRF_RFp4imp)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp4imp.pdf", width=4, height=4)
ROC_testRF_RFp4imp <- roc(response=test_HR_teRF_p4$PTB, predictor=pred_testRF_RFp4imp$Preterm)
plot(ROC_testRF_RFp4imp)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp4imp),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp4imp)
dev.off()
########################################



######################################## ROC curve for RFp1impRf on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp1impRf <- data.frame(predict(RFp1impRf, newdata=test_HR_teMM, type="prob"),
                                    actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp1impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp1impRf.pdf", width=4, height=4)
ROC_testMM_RFp1impRf <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp1impRf$Preterm)
plot(ROC_testMM_RFp1impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp1impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp1impRf)
dev.off()
########################################

######################################## ROC curve for RFp2impRf on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp2impRf <- data.frame(predict(RFp2impRf, newdata=test_HR_teMM, type="prob"),
                                    actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp2impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp2impRf.pdf", width=4, height=4)
ROC_testMM_RFp2impRf <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp2impRf$Preterm)
plot(ROC_testMM_RFp2impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp2impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp2impRf)
dev.off()
########################################

######################################## ROC curve for RFp3impRf on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp3impRf <- data.frame(predict(RFp3impRf, newdata=test_HR_teMM, type="prob"),
                                    actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp3impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp3impRf.pdf", width=4, height=4)
ROC_testMM_RFp3impRf <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp3impRf$Preterm)
plot(ROC_testMM_RFp3impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp3impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp3impRf)
dev.off()
########################################

######################################## ROC curve for RFp4impRf on test_HR_teMM (median/mode imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testMM_RFp4impRf <- data.frame(predict(RFp4impRf, newdata=test_HR_teMM, type="prob"),
                                    actual=test_HR_teMM$PTB)

# Look
head(test_HR_teMM %>% dplyr::select(PTB37))
head(pred_testMM_RFp4impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testMM_RFp4impRf.pdf", width=4, height=4)
ROC_testMM_RFp4impRf <- roc(response=test_HR_teMM$PTB, predictor=pred_testMM_RFp4impRf$Preterm)
plot(ROC_testMM_RFp4impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testMM_RFp4impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testMM_RFp4impRf)
dev.off()
########################################




######################################## ROC curve for RFp1impRf on test_HR_teRF_p1 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp1impRf <- data.frame(predict(RFp1impRf, newdata=test_HR_teRF_p1, type="prob"),
                                    actual=test_HR_teRF_p1$PTB)

# Look
head(test_HR_teRF_p1 %>% dplyr::select(PTB37))
head(pred_testRF_RFp1impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp1impRf.pdf", width=4, height=4)
ROC_testRF_RFp1impRf <- roc(response=test_HR_teRF_p1$PTB, predictor=pred_testRF_RFp1impRf$Preterm)
plot(ROC_testRF_RFp1impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp1impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp1impRf)
dev.off()
########################################

######################################## ROC curve for RFp2impRf on test_HR_teRF_p2 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp2impRf <- data.frame(predict(RFp2impRf, newdata=test_HR_teRF_p2, type="prob"),
                                    actual=test_HR_teRF_p2$PTB)

# Look
head(test_HR_teRF_p2 %>% dplyr::select(PTB37))
head(pred_testRF_RFp2impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp2impRf.pdf", width=4, height=4)
ROC_testRF_RFp2impRf <- roc(response=test_HR_teRF_p2$PTB, predictor=pred_testRF_RFp2impRf$Preterm)
plot(ROC_testRF_RFp2impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp2impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp2impRf)
dev.off()
########################################

######################################## ROC curve for RFp3impRf on test_HR_teRF_p3 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp3impRf <- data.frame(predict(RFp3impRf, newdata=test_HR_teRF_p3, type="prob"),
                                    actual=test_HR_teRF_p3$PTB)

# Look
head(test_HR_teRF_p3 %>% dplyr::select(PTB37))
head(pred_testRF_RFp3impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp3impRf.pdf", width=4, height=4)
ROC_testRF_RFp3impRf <- roc(response=test_HR_teRF_p3$PTB, predictor=pred_testRF_RFp3impRf$Preterm)
plot(ROC_testRF_RFp3impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp3impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp3impRf)
dev.off()
########################################

######################################## ROC curve for RFp4impRf on test_HR_teRF_p4 (RF imputation on missing data in test set)
# Generate data frame for thresholding & ROCs
set.seed(100)
pred_testRF_RFp4impRf <- data.frame(predict(RFp4impRf, newdata=test_HR_teRF_p4, type="prob"),
                                    actual=test_HR_teRF_p4$PTB)

# Look
head(test_HR_teRF_p4 %>% dplyr::select(PTB37))
head(pred_testRF_RFp4impRf)

# Plot
pdf("plots/high_risk_women/Supervised/ROC/ROC_testRF_RFp4impRf.pdf", width=4, height=4)
ROC_testRF_RFp4impRf <- roc(response=test_HR_teRF_p4$PTB, predictor=pred_testRF_RFp4impRf$Preterm)
plot(ROC_testRF_RFp4impRf)
text(0.2, 0.5, paste0("AUC = ", round(auc(ROC_testRF_RFp4impRf),3), sep="")) #text x&y positioning is from 0 to 1
auc(ROC_testRF_RFp4impRf)
dev.off()
########################################







###########################################################################################################################
##                                                                                                                       ##
##                                        Test set - threshold optimisation                                              ##
##                                                                                                                       ##
###########################################################################################################################


######################################################################### RFp1imp on testMM
################################ Confusion matrices for RFp1imp on testMM
# df of using this RF model on the test data (with MM imputation)
head(pred_testMM_RFp1imp)

# Add prediction if we set probability threshold of 0.15 (instead of 0.5) for classifying a prediction as "sPTB"
pred_testMM_RFp1imp$thresh0.15 <- factor(ifelse(pred_testMM_RFp1imp$Preterm > 0.15, "Preterm", "Term"),
                                         levels=c("Term", "Preterm"))
head(pred_testMM_RFp1imp)
str(pred_testMM_RFp1imp)


# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_testMM_RFp1imp$thresh0.15, reference=pred_testMM_RFp1imp$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     145      12
#    Preterm  138      54

################################

################################ Metrics for different probability threshold for RFp1imp on testMM
# Find best probability threshold RFp1imp
prob_RFp1imp_testMM <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_testMM_RFp1imp$Preterm >= Threshold, "Preterm", "Term"))
  
  conf <- confusionMatrix(data=predicted_response, reference=pred_testMM_RFp1imp$actual, positive="Preterm")
  
  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293
  
  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1imp_testMM) <- NULL #otherwise kable exports rownames
prob_RFp1imp_testMM # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1imp_testMM %>%
  knitr::kable(caption = "Metrics when optimising the classification threshold in the P1 rf model with median/mode imputation, on the training set where missing data was imputed with median/mode imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 

################################
#########################################################################


######################################################################### RFp1imp on testRF
################################ Confusion matrices for RFp1imp on testRF
# df of using this RF model on the test data (with RF imputation)
head(pred_testRF_RFp1imp)

# Add prediction if we set probability threshold of 0.15 (instead of 0.5) for classifying a prediction as "sPTB"
pred_testRF_RFp1imp$thresh0.15 <- factor(ifelse(pred_testRF_RFp1imp$Preterm > 0.15, "Preterm", "Term"),
                                         levels=c("Term", "Preterm"))
head(pred_testRF_RFp1imp)
str(pred_testRF_RFp1imp)


# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_testRF_RFp1imp$thresh0.15, reference=pred_testRF_RFp1imp$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     137      12
#    Preterm  146      54

################################

################################ Metrics for different probability threshold for RFp1imp on testRF
# Find best probability threshold RFp1imp
prob_RFp1imp_testRF <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_testRF_RFp1imp$Preterm >= Threshold, "Preterm", "Term"))
  
  conf <- confusionMatrix(data=predicted_response, reference=pred_testRF_RFp1imp$actual, positive="Preterm")
  
  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293
  
  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1imp_testRF) <- NULL #otherwise kable exports rownames
prob_RFp1imp_testRF # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1imp_testRF %>%
  knitr::kable(caption = "Metrics when optimising the classification threshold in the P1 rf model with median/mode imputation, on the training set where missing data was imputed with rf imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 

################################
#########################################################################



######################################################################### RFp1impRf on testMM
################################ Confusion matrices for RFp1impRf on testMM
# df of using this RF model on the test data (with MM imputation)
head(pred_testMM_RFp1impRf)

# Add prediction if we set probability threshold of 0.15 (instead of 0.5) for classifying a prediction as "sPTB"
pred_testMM_RFp1impRf$thresh0.15 <- factor(ifelse(pred_testMM_RFp1impRf$Preterm > 0.15, "Preterm", "Term"),
                                           levels=c("Term", "Preterm"))
head(pred_testMM_RFp1impRf)
str(pred_testMM_RFp1impRf)


# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_testMM_RFp1impRf$thresh0.15, reference=pred_testMM_RFp1impRf$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     145      14
#    Preterm  138      52

################################

################################ Metrics for different probability threshold for RFp1impRf on testMM
# Find best probability threshold RFp1impRf
prob_RFp1impRf_testMM <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_testMM_RFp1impRf$Preterm >= Threshold, "Preterm", "Term"))
  
  conf <- confusionMatrix(data=predicted_response, reference=pred_testMM_RFp1impRf$actual, positive="Preterm")
  
  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293
  
  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1impRf_testMM) <- NULL #otherwise kable exports rownames
prob_RFp1impRf_testMM # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1impRf_testMM %>%
  knitr::kable(caption = "Metrics when optimising the classification threshold in the P1 rf model with rf imputation, on the training set where missing data was imputed with median/mode imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 

################################
#########################################################################



######################################################################### RFp1impRf on testRF
################################ Confusion matrices for RFp1impRf on testRF
# df of using this RF model on the test data (with RF imputation)
head(pred_testRF_RFp1impRf)

# Add prediction if we set probability threshold of 0.3 (instead of 0.5) for classifying a prediction as "sPTB"
pred_testRF_RFp1impRf$thresh0.15 <- factor(ifelse(pred_testRF_RFp1impRf$Preterm > 0.15, "Preterm", "Term"),
                                           levels=c("Term", "Preterm"))
head(pred_testRF_RFp1impRf)
str(pred_testRF_RFp1impRf)


# Look at confusion matrix probability threshold of 0.15
confusionMatrix(data=pred_testRF_RFp1impRf$thresh0.15, reference=pred_testRF_RFp1impRf$actual, positive="Preterm")
#             Reference
# Prediction Term Preterm
#    Term     142      14
#    Preterm  141      52

################################

################################ Metrics for different probability threshold for RFp1impRf on testRF
# Find best probability threshold RFp1impRf
prob_RFp1impRf_testRF <- map_dfr(seq(0.05,0.5,0.05),  function(Threshold){
  predicted_response <- as.factor(ifelse(pred_testRF_RFp1impRf$Preterm >= Threshold, "Preterm", "Term"))
  
  conf <- confusionMatrix(data=predicted_response, reference=pred_testRF_RFp1impRf$actual, positive="Preterm")
  
  Accuracy <- conf$overall[1]
  Sensitivity <- conf$byClass[1]
  Specificity <- conf$byClass[2]
  PPV <- conf$byClass[3]
  NPV <- conf$byClass[4]
  #G_mean <- sqrt(as.numeric(conf$overall[1])*as.numeric(conf$byClass[2])) #https://towardsdatascience.com/optimal-threshold-for-imbalanced-classification-5884e870c293
  
  data.frame(Threshold, Accuracy, Sensitivity, Specificity, PPV, NPV) #, G_mean
})
rownames(prob_RFp1impRf_testRF) <- NULL #otherwise kable exports rownames
prob_RFp1impRf_testRF # see what is best
# Threshold  Accuracy Sensitivity Specificity       PPV       NPV    

# Export
?knitr::kable
prob_RFp1impRf_testRF %>%
  knitr::kable(caption = "Metrics when optimising the classification threshold in the P1 rf model with rf imputation, on the training set where missing data was imputed with rf imputation",   
               format = "latex",
               digits=3 , booktabs = TRUE) 

################################
#########################################################################






