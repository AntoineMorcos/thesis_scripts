# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(caret)


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

# Loading Rdata from the script "1Data_cleaning.R"
lnames <- load(file="data/1Data_cleaning_output.RData")
lnames

# Loading Rdata from the script "2Unsupervised_learning.R"
lnames <- load(file="data/2Unsupervised_learning_output_high_risk_women.RData")
lnames






###########################################################################################################################
##                                                                                                                       ##
##                                          df prep before partitioning                                                  ##
##                                                                                                                       ##
###########################################################################################################################

# Look
dim(HR_clin) #1169   45
names(HR_clin)
dim(final_cluster_assignment_HR) #1169    2
head(final_cluster_assignment_HR)


# Pre split clinical data in high risk women
PS_HR <- base::merge(HR_clin, final_cluster_assignment_HR, by="Participant.ID")
dim(PS_HR) #1169   46

# Fix rownames
row.names(PS_HR) <- paste0("ID_", PS_HR$Participant.ID, sep="")
head(PS_HR)

# Delete ID column
PS_HR$Participant.ID <- NULL




###########################################################################################################################
##                                                                                                                       ##
##                                                Data Partitioning                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Ratios planned - 70 train :  30 test
?caret::createDataPartition

# Set seed to always give the same IDs in each set
set.seed(100)



###################################### Partition_TrainTest
# Partition data into train & test
set.seed(100)
Partition_TrainTest <- createDataPartition(y=PS_HR$PTB37, #the class label, caret ensures an even split of classes
                                           p=0.7, #proportion of samples assigned to train
                                           list=FALSE #logical - should the results be in a list (TRUE) or a matrix 
) 

# Look
head(Partition_TrainTest)

# Make training set
train_HR <- PS_HR[Partition_TrainTest,] #take the corresponding rows for training

# Make test set
test_HR  <- PS_HR[-Partition_TrainTest,] #take the corresponding rows for test, by removing training rows

# Looking
nrow(train_HR) #820
nrow(test_HR) #349
1169-820-349
######################################





###################################### IDs check for future
train_IDs <- data.frame(ID=row.names(train_HR), set="train")
test_IDs <- data.frame(ID=row.names(test_HR), set="test")

HR_test_train <- rbind(train_IDs, test_IDs)

write_csv(HR_test_train, file=paste0("data/partitioning_check/HR_train_and_test_IDs_", as.character(Sys.Date()), ".csv"))
# cd "/c/Users/alici/OneDrive - King's College London/PhD/Projects/ML_clinical/data/partitioning_check"
# md5sum *
# diff file1 file2
######################################








###########################################################################################################################
##                                                                                                                       ##
##                                                Data normalising                                                        ##
##                                                                                                                       ##
###########################################################################################################################


###################################### Near-Zero Variance 
# Info
?caret::nzv #nzv is the original version of the function.
?caret::nearZeroVar #newer version

# Which data is numeric?
str(train_HR)
train_HR_numONLY <- train_HR %>% dplyr::select(Age, BMI, pH, CL_minimum, Early_miscarriages_number, Late_miscarriages_number, IMD_rank, IMD_decile)
test_HR_numONLY <- test_HR %>% dplyr::select(Age, BMI, pH, CL_minimum, Early_miscarriages_number, Late_miscarriages_number, IMD_rank, IMD_decile)

# Cols with Near-Zero Variance
nearZeroVar(train_HR_numONLY, names=T) # no numeric data with Near-Zero Variance
nearZeroVar(test_HR_numONLY, names=T) # no numeric data with Near-Zero Variance
######################################


###################################### Correlation between features
# Info
?cor
?caret::findCorrelation
?psych::corr.test

# Calculate correlation matrix on the numeric predictors
calculateCor2 <- psych::corr.test(train_HR_numONLY, adjust = "fdr") 
calculateCor2 # IMD rank & decile are highly correlated of course but will just pick one for the ML. Otherwise nothing is highly correlated out of the numeric predictors so I don't need to remove any. 
# Call:psych::corr.test(x = train_HR_numONLY, adjust = "fdr")
# Correlation matrix 
#                             Age   BMI    pH CL_minimum Early_miscarriages_number Late_miscarriages_number IMD_rank IMD_decile
# Age                        1.00 -0.01 -0.14       0.07                      0.12                     0.10     0.09       0.09
# BMI                       -0.01  1.00  0.11      -0.06                      0.13                     0.33    -0.12      -0.12
# pH                        -0.14  0.11  1.00      -0.14                      0.01                     0.15    -0.01      -0.01
# CL_minimum                 0.07 -0.06 -0.14       1.00                     -0.04                    -0.27     0.05       0.04
# Early_miscarriages_number  0.12  0.13  0.01      -0.04                      1.00                     0.14    -0.04      -0.04
# Late_miscarriages_number   0.10  0.33  0.15      -0.27                      0.14                     1.00    -0.13      -0.12
# IMD_rank                   0.09 -0.12 -0.01       0.05                     -0.04                    -0.13     1.00       0.99
# IMD_decile                 0.09 -0.12 -0.01       0.04                     -0.04                    -0.12     0.99       1.00
#
# Sample Size 
#                           Age BMI  pH CL_minimum Early_miscarriages_number Late_miscarriages_number IMD_rank IMD_decile
# Age                       820 817 699        760                       819                      819      724        724
# BMI                       817 817 697        758                       816                      816      722        722
# pH                        699 697 699        671                       698                      698      611        611
# CL_minimum                760 758 671        760                       759                      759      670        670
# Early_miscarriages_number 819 816 698        759                       819                      819      723        723
# Late_miscarriages_number  819 816 698        759                       819                      819      723        723
# IMD_rank                  724 722 611        670                       723                      723      724        724
# IMD_decile                724 722 611        670                       723                      723      724        724
#
# Probability values (Entries above the diagonal are adjusted for multiple tests.) 
#                             Age  BMI   pH CL_minimum Early_miscarriages_number Late_miscarriages_number IMD_rank IMD_decile
# Age                       0.00 0.77 0.00       0.11                      0.00                     0.01     0.02       0.02
# BMI                       0.74 0.00 0.01       0.19                      0.00                     0.00     0.00       0.00
# pH                        0.00 0.00 0.00       0.00                      0.78                     0.00     0.77       0.77
# CL_minimum                0.07 0.13 0.00       0.00                      0.36                     0.00     0.28       0.32
# Early_miscarriages_number 0.00 0.00 0.78       0.31                      0.00                     0.00     0.32       0.36
# Late_miscarriages_number  0.00 0.00 0.00       0.00                      0.00                     0.00     0.00       0.00
# IMD_rank                  0.01 0.00 0.75       0.20                      0.25                     0.00     0.00       0.00
# IMD_decile                0.01 0.00 0.74       0.25                      0.31                     0.00     0.00       0.00
######################################


###################################### Skewness & Scaling
# Info
?preProcess

# Looking at skewness
histogram(train_HR$Age) #normal distribution
histogram(train_HR$BMI) #a bit right-skewed
histogram(train_HR$pH) #right-skewed
histogram(train_HR$CL_minimum) # a bit left-skewed
histogram(train_HR$Early_miscarriages_number) #right-skewed
histogram(train_HR$Late_miscarriages_number) #right-skewed
histogram(train_HR$IMD_rank) #a bit right-skewed?
histogram(train_HR$IMD_decile) #a bit right-skewed?
dev.off()

# Calculate preProcess
set.seed(100)
?preProcess
calculatePreProcess <- preProcess(train_HR, #a matrix or df. Non-numeric predictors are allowed but will be ignored.
                                  method = c("center", "scale", "BoxCox"), #perform preprocessing
                                  na.remove = T) #a logical; should missing values be removed from the calculations?
# I didn't include "corr","nzv" as we did these above & they weren't a problem here

# Look
calculatePreProcess
# Created from 439 samples and 45 variables
# 
# Pre-processing:
#   - Box-Cox transformation (5)
# - centered (8)
# - ignored (37)
# - scaled (8)
# 
# Lambda estimates for Box-Cox transformation:
#   1.3, -1, -2, 0.5, 0.4


# Apply preprocessing to train_HR
set.seed(100)
train_HR_norm <- predict(calculatePreProcess, train_HR) 
summary(train_HR_norm[,2:10])

# Apply preprocessing to test_HR
set.seed(100)
test_HR_norm <- predict(calculatePreProcess, test_HR) 
summary(test_HR_norm[,2:10])
######################################






###########################################################################################################################
##                                                                                                                       ##
##                                                 Exporting data                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Saving data
save(train_HR_norm, test_HR_norm, # partitioned, normalised data (scaled, centered, BoxCox)
     train_HR, test_HR, # partitioned, but non-normalised data
     file="data/3Data_preprocessing_output_high_risk_women.RData")



# # Loading Rdata from the script "3Data_preprocessing_for_supervised__high_risk_women.R"
# lnames <- load(file="data/3Data_preprocessing_output_high_risk_women.RData")
# lnames




