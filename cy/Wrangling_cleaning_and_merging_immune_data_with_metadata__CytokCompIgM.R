# Clear environment
rm(list=ls())

# Load libraries
library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(dplyr)
library(tidyr)
library(BiocManager) #install.packages("BiocManager")
library(devtools)
library(ggbiplot) #install_github("vqv/ggbiplot")
library(ggrepel)
library(sva) #BiocManager::install("sva")
library(data.table) #filtering dfs
library(corrplot)#install.packages('corrplot')
library(ggpubr) #adj ps on ggplot 
library(rstatix) #adj ps on ggplot

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()
#list.files()

#Versions of packages used
sessionInfo()






###########################################################################################################################
##                                                                                                                       ##
##                                             Load in data                                                              ##
##                                                                                                                       ##
###########################################################################################################################

# 30/06/2022 Someone in the lab incorrectly labelled/matched samples to the the barcodes/IDs, so will now redo with the real data 
cytok_raw <- read.csv("data/from_Amirah/cytokine_results_with_ID_all_plates.csv", na.strings = c("", "NaN", "NA"))
comp_raw <- read.csv("data/from_Amirah/complement_results_with_ID_all_plates.csv", na.strings = c("", "NaN", "NA"))
IgM_raw <- read.csv("data/from_Amirah/IgM_results_with_ID_all_plates.csv", na.strings = c("", "NaN", "NA"))

# Import Flavia's new metadata with metabolites, clinical data, etc.
metadataFlavia <- read.csv("data/from_Flavia/Metadata_final_May2022.csv") 

# Import Flavia's 2023 metadata with corrected pH and IMD scores
NEWmetadataFlavia <- read.csv("data/from_Flavia/metadata_complete_January2023_withPcoaaxis_imuSmall_AliciaClusters.csv") 

# Import metadata I downloaded from medscinet for IDs XXXX
metadataHUGE <- read.csv("data/medscinet_metadata_for_IDs_included_here.csv", sep="|", na.strings = c("", "NaN", "NA")) 






###########################################################################################################################
##                                                                                                                       ##
##                                               IDs to exclude                                                          ##
##                                                                                                                       ##
###########################################################################################################################

# Vector of IDs to exclude
IDs_to_exclude <- as.character(c(XXXX, XXXX, XXXX, XXXX)) #iPTB IDs removed (IDs are kept confidential here in this public repository)
IDs_to_exclude 







###########################################################################################################################
##                                                                                                                       ##
##                                           Data cleaning immune dfs                                                    ##
##                                                                                                                       ##
###########################################################################################################################

####################################### cytok
# Only keeping relevant columns
names(cytok_raw)
cytok_only <- cytok_raw[,c(1,29:38, 41:50)]
names(cytok_only)

# Getting rid of any controls or measurements which don't belong to an ID & Renaming ID colummn
table(is.na(cytok_only$insight_id))
cytok_only <- cytok_only %>% 
  dplyr::rename(Participant.ID=insight_id) %>% 
  dplyr::filter(Participant.ID!="NA")
table(is.na(cytok_only$Participant.ID)) #NAs gone 

# Get rid of rows/IDs where all measurements are NAs 
measurments_cytok <- names(cytok_only[-1])
dim(cytok_only) #111  21
cytok_only <- cytok_only %>%
  dplyr::filter(!if_all(c(measurments_cytok), is.na))
dim(cytok_only) #110  21 #


# Making <X equal to 0 as per recommendations
# Changing <Xs in the df to blanks 
cytok_only[75:84,1:5]
int1 <- as.data.frame(apply(cytok_only,2,function(x)gsub('<', '_', x, fixed=T))) #converts < to _ (we need this intermediate step otherwise \\<.* doesn't work)
int1[75:84,1:5]
cytok_clean <- as.data.frame(apply(int1,2,function(x)gsub('\\_.*', '0', x))) #deletes all the characters in the column names after the "_" character
cytok_clean[75:84,1:5]

# Changing NAs to 0s as per recomendations
cytok_clean[is.na(cytok_clean)] = "0" #converts NAs to 0s 
cytok_clean[75:84,1:5]

##### Only 1 of all IDs:
# Looking
table(table(cytok_clean$Participant.ID))
table(cytok_clean$Participant.ID) # 2 of XXXX
cytok_clean %>% dplyr::filter(Participant.ID=="XXXX")
table(duplicated(cytok_clean)) # yes 1 duplicate row
dim(cytok_clean) #110  21
# Delete duplicate rows
cytok_clean <- unique(cytok_clean)
dim(cytok_clean) #109  21
table(duplicated(cytok_clean)) # All F - i.e. 1 row / ID
##### 
#######################################


####################################### comp
# Only keeping relevant columns
names(comp_raw)
comp_only <- comp_raw %>% dplyr::select(insight_id, C5.1, C5a.1, C5.2, C5a.2)
names(comp_only)

# Getting rid of any controls or measurements which don't belong to an ID & Renaming ID colummn
table(is.na(comp_only$insight_id))
comp_only <- comp_only %>% 
  dplyr::rename(Participant.ID=insight_id) %>% 
  dplyr::filter(Participant.ID!="NA")
table(is.na(comp_only$Participant.ID)) #NAs gone 

# Get rid of rows/IDs where all measurements are NAs 
measurments_comp <- names(comp_only[-1])
dim(comp_only) # 111   5 
comp_only <- comp_only %>%
  dplyr::filter(!if_all(c(measurments_comp), is.na))
dim(comp_only) # 110   5


# Making <X equal to 0 as it's too small for the machine to detect it as per recommendations
# Changing <s in the df to blanks 
comp_only[75:84,1:5]
int1 <- as.data.frame(apply(comp_only,2,function(x)gsub('<', '_', x, fixed=T))) #converts < to _ (we need this intermediate step otherwise \\<.* doesn't work)
int1[75:84,1:5]
comp_clean <- as.data.frame(apply(int1,2,function(x)gsub('\\_.*', '0', x))) #deletes all the characters in the column names after the "_" character
comp_clean[75:84,1:5]


# Changing NAs to 0s as per recommendations
comp_clean[is.na(comp_clean)] = "0" #converts NAs to 0s 
comp_clean[75:84,1:5]

##### Only 1 of all IDs:
# Looking
table(table(comp_clean$Participant.ID))
table(comp_clean$Participant.ID) # 2 of XXXX - like cytok was too
comp_clean %>% dplyr::filter(Participant.ID=="XXXX")
table(duplicated(comp_clean)) # yes 1 duplicate row
dim(comp_clean) #110  5
# Delete duplicate rows
comp_clean <- unique(comp_clean)
dim(comp_clean) #109  5
table(duplicated(comp_clean)) # All F - i.e. 1 row / ID
#####
#######################################


####################################### IgM
# Only keeping relevant columns
names(IgM_raw)
IgM_only <- IgM_raw %>% dplyr::select(insight_id, IgM.1, IgM.2)
names(IgM_only)

# Getting rid of any controls or measurements which don't belong to an ID & Renaming ID colummn
table(is.na(IgM_only$insight_id))
IgM_only <- IgM_only %>% 
  dplyr::rename(Participant.ID=insight_id) %>% 
  dplyr::filter(Participant.ID!="NA")
table(is.na(IgM_only$Participant.ID)) #NAs gone 

# Get rid of rows/IDs where all measurements are NAs 
measurments_IgM <- names(IgM_only[-1])
dim(IgM_only) # 110   3
IgM_only <- IgM_only %>%
  dplyr::filter(Participant.ID!="XXXX" & Participant.ID!="XXXX") #!if_all(c(measurments_IgM), is.na) 
dim(IgM_only) #108   3

# There's no "<" symbols in IgM_only so I don't need to convert them like I did in cytok & comp
IgM_clean <- IgM_only

# Changing NAs to 0s (proper NAs have been removed from IgM_only df)
IgM_clean[is.na(IgM_clean)] = "0" #converts NAs to 0s 
IgM_clean[75:84,]

##### Only 1 of all IDs:
# Looking
table(table(IgM_clean$Participant.ID))
table(IgM_clean$Participant.ID) 
IgM_clean %>% dplyr::filter(Participant.ID=="XXXX")
# not a duplicate, but a different measurement
table(duplicated(IgM_clean)) # no duplicate rows, but 1 duplicate ID 
dim(IgM_clean) #108 3
# Delete duplicate rows
IgM_clean <- IgM_clean[!duplicated(IgM_clean$Participant.ID), ] # Apply duplicated
IgM_clean %>% dplyr::filter(Participant.ID=="XXXX") # Great - have kept the measurement I was advised to keep
dim(IgM_clean) #107   3
table(duplicated(IgM_clean)) # All F - i.e. 1 row / ID
#####
#######################################








###########################################################################################################################
##                                                                                                                       ##
##                           Merging immune data into 1 df & adding mean measurement cols                                ##
##                                                                                                                       ##
###########################################################################################################################

#################### Merging dfs to make immu_1s2s
# Look
dim(cytok_clean) #109  21
dim(comp_clean) #109   5
dim(IgM_clean) #107   3
str(cytok_clean)
str(comp_clean)
str(IgM_clean)
IgM_clean$Participant.ID <- as.character(IgM_clean$Participant.ID) #fix str for merge

# Merging dfs
int_cytok_comp <- base::merge(cytok_clean, comp_clean, by="Participant.ID", all=T)
dim(int_cytok_comp) #109  25
immu_1s2s <- base::merge(int_cytok_comp, IgM_clean, by="Participant.ID", all=T)
dim(immu_1s2s) #110  27 
#View(immu_1s2s) #XXXX has IgM only, no cytok/comp
####################



#################### Tidying
# Reorder columns alphabetically, but with Participant.ID at the front
immu_1s2s[1:10,1:10]
immu_1s2s <- immu_1s2s %>% 
  dplyr::select(order(colnames(.))) %>% 
  dplyr::select(Participant.ID, everything())
immu_1s2s[1:10,1:10]
####################




#################### Changing str to numeric 
# Looking
str(immu_1s2s) #all chr
table(is.na(immu_1s2s))

# Character to numeric - I had to remove non-sample measurements as they had a > symbol in them & I had to read NaNs as NAs when read.csv to make the as.numeric not coerce any values into NAs
immu_1s2s_num <- immu_1s2s %>% 
  mutate_at(c(2:length(immu_1s2s)), as.numeric)

# Checks
table(is.na(immu_1s2s_num))
str(immu_1s2s_num)
####################



# Calculating mean cytokines
names(immu_1s2s_num)
immu_big_ALLIDs <- immu_1s2s_num %>% 
  dplyr::mutate(C5= (C5.1 + C5.2)/2,
                C5a= (C5a.1 + C5a.2)/2,
                GCSF= (GCSF.1 + GCSF.2)/2,
                IFNg= (IFNg.1 + IFNg.2)/2,
                IgM= (IgM.1 + IgM.2)/2,
                IL10= (IL10.1 + IL10.2)/2,
                IL13= (IL13.1 + IL13.2)/2,
                IL17A= (IL17A.1 + IL17A.2)/2,
                IL22= (IL22.1 + IL22.2)/2,
                IL5= (IL5.1 + IL5.2)/2,
                IL6= (IL6.1 + IL6.2)/2,
                IL8= (IL8.1 + IL8.2)/2,
                MCP1= (MCP1.1 + MCP1.2)/2
  )

# Checking
immu_big_ALLIDs %>% 
  dplyr::select(Participant.ID, GCSF, GCSF.1, GCSF.2) %>% head()

# Reorder columns alphabetically, but with Participant.ID at the front
immu_big_ALLIDs[1:10,1:10]
immu_big_ALLIDs <- immu_big_ALLIDs %>% 
  dplyr::select(order(colnames(.))) %>% 
  dplyr::select(Participant.ID, everything())
immu_big_ALLIDs[1:10,1:10]

# Setting rownames
row.names(immu_big_ALLIDs) <- immu_big_ALLIDs$Participant.ID

# Keeping only columns of interest (means of cytokines, getting rid of the duplicate results)
names(immu_big_ALLIDs)
immu_small_ALLIDs <- as.data.frame(immu_big_ALLIDs[,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38)])
#immu_small_ALLIDs 
names(immu_big_ALLIDs)
names(immu_small_ALLIDs)
immu_small_ALLIDs[1:5,]









###########################################################################################################################
##                                                                                                                       ##
##                                             Add metadata to immu data                                                 ##
##                                                                                                                       ##
###########################################################################################################################

############################################ Looking
# Making a common ID column
#View(metadataFlavia)
str(immu_small_ALLIDs$Participant.ID)
metadataFlavia$Participant.ID <- as.character(metadataFlavia$id)
str(metadataFlavia$Participant.ID)


# Looking at the women in common b/ the dfs
dim(immu_small_ALLIDs) #110 women of cytokines data
dim(metadataFlavia) #87 women of metabolite and/or metagenome data
table(immu_small_ALLIDs$Participant.ID %in% metadataFlavia$Participant.ID) #FALSE: 29,  TRUE: 81
table(metadataFlavia$Participant.ID %in% immu_small_ALLIDs$Participant.ID) #FALSE:  6,  TRUE: 81
81+29+6 # 116 women in total, but we won't use all these women in final analysis as we later remove women without metag data
############################################



############################################ Filtering immu IDs from metadataHUGE the medscinet export
# IDs of all women in immu_small_ALLIDs
IDs_vector_immu <- immu_small_ALLIDs$Participant.ID


# Need to add metadata of women who are in cytokines but not metadataFlavia before merge to avoid NAs!
names(metadataHUGE) 


# Selecting only our immu IDs & renaming some columns
metadataHUGE_immuIDs <- metadataHUGE %>% 
  dplyr::filter(Participant.ID %in% c(IDs_vector_immu)) %>% 
  dplyr::rename("sPTB37" = "Spontaneous.onset.of.labour.resulting.in.delivery..37.40",
                "sPTB34" = "Spontaneous.onset.of.labour.resulting.in.delivery..34.40") %>%
  dplyr::select(Participant.ID, Ethnicity, sPTB34, sPTB37, Low.risk.at.enrolment, # Selecting all columns, but putting them in this order
                Premature.Prelabour.Rupture.of.Membranes, Onset.of.labour, Pregnancy.Outcome.Status, Preeclampsia, Smoking, 
                everything()) %>% 
  dplyr::distinct(Participant.ID, .keep_all = TRUE) #Before applying distinct there were multiple rows per ID for some IDs. But contents of this reduced col df is same for each duplicate row (same Date.of.discharge, bmi, etc.). So that's why I used distinct()

# Checks
#View(metadataHUGE_immuIDs)
dim(metadataHUGE_immuIDs) #110 119
names(metadataHUGE_immuIDs)
############################################



############################################ Merging both metadata dfs
# Looking
?data.table::rbindlist
?base::merge
dim(metadataHUGE_immuIDs)
dim(metadataFlavia)


# Merge 
?base::merge
metadata_FF_MSN <- base::merge(metadataFlavia, metadataHUGE_immuIDs, by="Participant.ID", 
                                all.x=T, all.y=F, # keep all Flavia metadata (x), & only keep medscinet data where we have Flavia/metag data (y)
                                suffixes = c(".Flavia", ".medscinet")) #Add suffixes to cols where the colnames are the same

# Look
dim(metadata_FF_MSN) #87 204
table(metadata_FF_MSN$eth, useNA="always") 
table(table(metadata_FF_MSN$Participant.ID)) # good - no duplicate IDs
############################################ 



############################################ Looking at sPTB conflicting data

# Flavia says to use her sPTB, as medscinet DB is sometimes wrong. Her values are all checked with Natasha.

# # Looking at sPTB outcomes
# table(metadata_FF_MSN$sPTB34.medscinet, metadata_FF_MSN$sPTB34.Flavia, useNA="always")
# table(metadata_FF_MSN$sPTB37.medscinet, metadata_FF_MSN$sPTB37.Flavia, useNA="always")

 
# Making sentence case to make them in the same format
metadata_FF_MSN$sPTB37.Flavia
metadata_FF_MSN$sPTB37.Flavia[metadata_FF_MSN$sPTB37.Flavia=="yes"] <- "Yes"
metadata_FF_MSN$sPTB37.Flavia[metadata_FF_MSN$sPTB37.Flavia=="no"] <- "No"
table(metadata_FF_MSN$sPTB37.Flavia, useNA="always")
metadata_FF_MSN$sPTB34.Flavia
metadata_FF_MSN$sPTB34.Flavia[metadata_FF_MSN$sPTB34.Flavia=="yes"] <- "Yes"
metadata_FF_MSN$sPTB34.Flavia[metadata_FF_MSN$sPTB34.Flavia=="no"] <- "No"
table(metadata_FF_MSN$sPTB34.Flavia, useNA="always")

# # Seeing which IDs differ
# metadata_FF_MSN$Test1 <- ifelse(metadata_FF_MSN$sPTB37.medscinet==metadata_FF_MSN$sPTB37.Flavia,"Match","DIFFERENT")
# metadata_FF_MSN %>% dplyr::filter(Test1!="Match")%>% dplyr::select(Participant.ID, Test1, sPTB37.medscinet, sPTB37.Flavia) 
# metadata_FF_MSN %>% dplyr::select(Participant.ID, Sample, Test1, sPTB37.medscinet, sPTB37.Flavia) 
# #metadata_FF_MSN$Test1 <- NULL
# 
#
# ############################################



############################################ Tidying
names(metadata_FF_MSN)

# Deleting columns
metadata_FF_MSN$sPTB34.medscinet <- NULL
metadata_FF_MSN$sPTB37.medscinet <- NULL
metadata_FF_MSN$Ethnicity <- NULL
metadata_FF_MSN$BMI <- NULL
metadata_FF_MSN$match_id <- NULL
metadata_FF_MSN$Age_group <- NULL
metadata_FF_MSN$BMI_group<- NULL
metadata_FF_MSN$Ravel_updated_csv <- NULL
metadata_FF_MSN$Preeclampsia.1 <- NULL
metadata_FF_MSN$Smoking <- NULL



# Gestation at delivery decimalised & renaming columns
metadata_FF_MSN <- metadata_FF_MSN %>% 
  dplyr::mutate(Gestation.at.delivery.wks.dec = ((ga_w*7) + ga_d) / 7) %>%  #number of weeks x7 plus number of days, divided by 7 again to put it back into weeks 
    dplyr::rename("sPTB34"="sPTB34.Flavia",
                  "sPTB37"="sPTB37.Flavia",
                  "sPTB34"="sPTB34.Flavia",
                  "BMI"="bmi",
                  "Age"="age",
                  "Ethnicity"="eth",
                  "Smoking"="smoking") 
    

#Making BMI_category column 
metadata_FF_MSN$BMI_category <- NA
metadata_FF_MSN$BMI_category[metadata_FF_MSN$BMI<18.5] <- "Underweight"
metadata_FF_MSN$BMI_category[metadata_FF_MSN$BMI>=18.5 & metadata_FF_MSN$BMI<25] <- "Healthy weight"
metadata_FF_MSN$BMI_category[metadata_FF_MSN$BMI>=25 & metadata_FF_MSN$BMI<30] <- "Overweight"
metadata_FF_MSN$BMI_category[metadata_FF_MSN$BMI>=30 & metadata_FF_MSN$BMI<40] <- "Obese"
metadata_FF_MSN$BMI_category[metadata_FF_MSN$BMI>=40 ] <- "Morbidly obese"
str(metadata_FF_MSN$BMI_category)
metadata_FF_MSN$BMI_category <- factor(metadata_FF_MSN$BMI_category, levels=c("Underweight", "Healthy weight", "Overweight", "Obese", "Morbidly obese")) 
table(metadata_FF_MSN$BMI_category, useNA="always")
tapply(metadata_FF_MSN$BMI, metadata_FF_MSN$BMI_category, summary)


#Making Age_category column 
summary(metadata_FF_MSN$Age)
metadata_FF_MSN$Age_category <- NA
metadata_FF_MSN$Age_category[metadata_FF_MSN$Age>=20 & metadata_FF_MSN$Age<25] <- "20-24"                          
metadata_FF_MSN$Age_category[metadata_FF_MSN$Age>=25 & metadata_FF_MSN$Age<30] <- "25-29"
metadata_FF_MSN$Age_category[metadata_FF_MSN$Age>=30 & metadata_FF_MSN$Age<35] <- "30-34"
metadata_FF_MSN$Age_category[metadata_FF_MSN$Age>=35 & metadata_FF_MSN$Age<40] <- "35-39"
metadata_FF_MSN$Age_category[metadata_FF_MSN$Age>=40 ] <- "40+"
metadata_FF_MSN$Age_category
metadata_FF_MSN$Age_category <- factor(metadata_FF_MSN$Age_category, levels=c("20-24", "25-29", "30-34", "35-39", "40+")) 
table(metadata_FF_MSN$Age_category, useNA="always")
tapply(metadata_FF_MSN$Age, metadata_FF_MSN$Age_category, summary)
############################################ 



############################################ Adding IMD index and corrected pH from 2023 data
# Look
names(NEWmetadataFlavia)

# Just columns we want to add
metadata2023 <- NEWmetadataFlavia %>% 
  dplyr::select(id, IMD_rank, IMD_decile, medscinet_ph) %>% 
  dplyr::rename(Participant.ID=id)
head(metadata2023)
dim(metadata2023)

# Remove pHs pre merge from 1 df
metadata_FF_MSN$medscinet_ph <- NULL
metadata_FF_MSN[1:5,1:5]


# All metadata
metadata_ALL <- base::merge(metadata2023, metadata_FF_MSN, by="Participant.ID", 
                            all.x=T, all.y=T, # keep all
                            suffixes = c(".2023", ".old")) #Add suffixes to cols where the colnames are the same

# Look
metadata_ALL[1:5,1:5]
dim(metadata_ALL) #[1]  87 199
############################################ 





############################################ Merging immu_small_ALLIDs with metadata 
# Merge the immu_small_ALLIDs data with metadata_ALL
?base::merge
immu_MD <- base::merge(immu_small_ALLIDs, metadata_ALL, by="Participant.ID", 
                       all.x=F, all.y=T) # only keep immune data where we have metadata data (x) & keep all metadata (y)
dim(immu_MD) #87 women x 212 cols of data
#View(immu_MD)

# Set rownames
rownames(immu_MD) <- immu_MD$Participant.ID
head(immu_MD)

# Order cols
names(immu_MD)
immu_MD <- immu_MD %>% 
  dplyr::select(Participant.ID, sPTB34, sPTB37, Gestation.at.delivery.wks.dec, 
                Ethnicity, Risk, BMI, BMI_category, Age, Age_category, Smoking,
                CL_absolute_small_value_at_any_visit, Short_cervix_value_this_visit,
                IMD_rank, IMD_decile, medscinet_ph,
                everything())
head(names(immu_MD), 20)
############################################


############################################ Checking for iPTBs
# Making sure excluded IDs aren't present
table(immu_MD$Participant.ID %in% IDs_to_exclude)

# Checking for iPTBs - none here
metadata_ALL %>% dplyr::filter(Gestation.at.delivery.wks.dec<37 & sPTB37=="No") %>% dplyr::select(Participant.ID, sPTB37, Gestation.at.delivery.wks.dec, Preeclampsia, Maternal.Infection, Onset.of.labour, Other.reason.details) #Now we have no <37w that aren't sPTB (i.e. all iPTB have been removed)
############################################







###########################################################################################################################
##                                                                                                                       ##
##                                        immu data of only women with metag                                             ##
##                                                                                                                       ##
###########################################################################################################################

# IDs of all women in immu_MD
IDs_vector_metagImmune <- immu_MD$Participant.ID

# Selecting only our immu IDs 
immu_big <- immu_big_ALLIDs %>% 
  dplyr::filter(Participant.ID %in% c(IDs_vector_metagImmune))
immu_small <- immu_small_ALLIDs %>% 
  dplyr::filter(Participant.ID %in% c(IDs_vector_metagImmune))

# Checks
dim(immu_big_ALLIDs) #110 women with immu
dim(immu_small_ALLIDs)
dim(immu_MD) #87 women with metag
dim(immu_big) # 81 women with immu that also have metag
dim(immu_small)







###########################################################################################################################
##                                                                                                                       ##
##                                                   Exporting                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Exporting cytokines/complements/IgM data as RData
save(immu_big,  #means of immune markers + their duplicate measurements + INSIGHT ID (81 women with immu that also have metag)
     immu_small, #means of immune markers + INSIGHT ID (81 women with immu that also have metag)
     immu_MD, #immune markers df with all metadata - not tidy (87 women with metag)
     immu_small_ALLIDs, immu_big_ALLIDs, #immune markers of the full 106 women (including White women, & may contain iPTBs)
     file="data/immune_CytokinesComplementsIgM_2023.RData")

# # Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
# lnames = load(file="data/immune_CytokinesComplementsIgM.RData")
# lnames








