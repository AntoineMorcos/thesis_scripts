#N of 312&317 in old df,  305& 313 in the previously used unmatched sample pair df. Now we have 302 women with 604 samples.

#clear environment
rm(list=ls())

#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc on windows
getwd()
#list.files()


library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(calibrate)

############################################################################################
##                                                                                        ##
##               Only keep women where we have an early & a late sample                   ##
##                                                                                        ##
############################################################################################

##Loading data in:
raw_data <- read.csv("data_BV/Final_Db_04_07_2019_noIUD_Alicia_edit_CSV.csv", 
                     header=TRUE, na.strings=c("NA", " ", ""), fileEncoding = 'UTF-8-BOM')


###Deleting people where we do not have their data at both sample timing groups 
###i.e. Erasing data where we only have 1 sample per person, this we will mean all remaining people will have 2 samples (10-15 & 16-24).

##Look at df - some IDs have 2 entries, others only 1
table(raw_data$participant_id) 
dim(raw_data)

##Create df of just rows where we have 2 entries (i.e. 2 rows for each participant_id). 
data1 <- raw_data[raw_data$participant_id %in% raw_data$participant_id[duplicated(raw_data$participant_id)],]

##Checks
table(data1$participant_id)
dim(data1)


############################################################################################
##                                                                                        ##
##                   Adding new columns & modifying existing ones                         ##
##                                                                                        ##
############################################################################################

#adding a gestation in delivery weeks in decimals col
data1 <- data1 %>% 
  mutate(Gestation.at.delivery.wks.dec = Gestation.at.delivery..days / 7) 
plot(data1$Gestation.at.delivery..days, data1$Gestation.at.delivery.wks.dec)


#changing row names to Sample_IDs
rownames(data1) <- data1$Sample_ID
data1[1:5, 1:4]

#Changing term classifications in column Term_outcome
table(data1$Term_outcome)
data1$Term_outcome <- as.character(data1$Term_outcome)
data1$Term_outcome[data1$Term_outcome == "PT =<33"] <- "PT <34"
data1$Term_outcome[data1$Term_outcome == "36 >= PT >=34"] <- "PT 34 to <37"
data1$Term_outcome[data1$Term_outcome == "T >= 37"] <- "T >=37"
data1$Term_outcome <- factor(data1$Term_outcome, levels = c("PT <34", "PT 34 to <37", "T >=37"))
table(data1$Term_outcome)

tapply(data1$Gestation.at.delivery.wks.dec, data1$Term_outcome, summary)

#Making new sPTB binary data
table(data1$Sptb.34w) 
data1$sPTB_34w<- data1$Sptb.34w
data1$sPTB_34w <- as.character(data1$sPTB_34w)
data1$sPTB_34w[data1$sPTB_34w == "No"] <- "No"
data1$sPTB_34w[data1$sPTB_34w == "No "] <- "No"
data1$sPTB_34w[data1$sPTB_34w == "Yes"] <- "Yes"
data1$sPTB_34w <- factor(data1$sPTB_34w, levels=c("No", "Yes"))
table(data1$sPTB_34w, data1$Sptb.34w) # all fine

table(data1$Sptb.37w) 
data1$sPTB_37w<- data1$Sptb.37w
data1$sPTB_37w <- as.character(data1$sPTB_37w)
data1$sPTB_37w[data1$sPTB_37w == "No"] <- "No"
data1$sPTB_37w[data1$sPTB_37w == "Yes"] <- "Yes"
data1$sPTB_37w <- factor(data1$sPTB_37w, levels=c("No", "Yes"))
table(data1$sPTB_37w, data1$Sptb.37w)



####Making a new BV categories column, which is shorter to take up less room on the plots
table(data1$Bacterial.Vaginosis..categories) 
data1$BV_categories <- data1$Bacterial.Vaginosis..categories #making a copy of the df
data1$BV_categories <- as.character(data1$BV_categories)
data1$BV_categories[data1$BV_categories == "Bacterial vaginosis (3)"] <- "BV"
data1$BV_categories[data1$BV_categories == "Intermediate (2)"] <- "Intermediate"
data1$BV_categories[data1$BV_categories == "Normal (0-1)"] <- "Normal"
data1$BV_categories <- factor(data1$BV_categories, levels=c("Normal", "Intermediate", "BV"))
table(data1$BV_categories, data1$Bacterial.Vaginosis..categories)


#making BV grading a factor
str(data1$BV.grading)
data1$BV.grading <-  as.factor(data1$BV.grading)


#making IDs characters
str(data1$Sample_ID)
str(data1$participant_id)
data1$participant_id <- as.character(data1$participant_id)
summary(data1$participant_id)



#making short cervix a factor
table(data1$short.cervix) 
data1$short.cervix_copy<- data1$short.cervix
data1$short.cervix <- as.character(data1$short.cervix)
data1$short.cervix[data1$short.cervix == "0"] <- "No"
data1$short.cervix[data1$short.cervix == "1"] <- "Yes"
data1$short.cervix <- factor(data1$short.cervix, levels=c("No", "Yes"))
table(data1$short.cervix, data1$short.cervix_copy)
data1$short.cervix_copy <- NULL

##Changing names of the ABCDEU to make it clear what they mean:
data1$PCoA_groups <- data1$status #I previously used $group.pcoa.Erica, but have been advised by Flavia to use $status instead (25/08/2020)
table(data1$PCoA_groups)
data1$PCoA_groups <- as.character(data1$PCoA_groups)
data1$PCoA_groups[data1$PCoA_groups == "A"] <- "A) L. crispatus"
data1$PCoA_groups[data1$PCoA_groups == "B"] <- "B) L. gasseri"
data1$PCoA_groups[data1$PCoA_groups == "C"] <- "C) L. iners"
data1$PCoA_groups[data1$PCoA_groups == "D"] <- "D) Mixed dysbiotic"
data1$PCoA_groups[data1$PCoA_groups == "E"] <- "E) L. crispatus/gasseri"
#data1$PCoA_groups[data1$PCoA_groups == "U"] <- "U) Unclear"
data1$PCoA_groups <- as.factor(data1$PCoA_groups)
table(data1$PCoA_groups, data1$status, useNA="always")



#Making a new column based on conditions of another column - BMI
data1$BMI_category <- "sgdfklgskldfj"
data1$BMI_category[data1$BMI<18.5] <- "Underweight"
data1$BMI_category[data1$BMI>=18.5 & data1$BMI<25] <- "Healthy weight"
data1$BMI_category[data1$BMI>=25 & data1$BMI<30] <- "Overweight"
data1$BMI_category[data1$BMI>=30 & data1$BMI<40] <- "Obese"
data1$BMI_category[data1$BMI>=40 ] <- "Morbidly obese"

str(data1$BMI_category)
data1$BMI_category <- factor(data1$BMI_category, levels=c("Underweight", "Healthy weight", "Overweight", "Obese", "Morbidly obese")) 

tapply(data1$BMI, data1$BMI_category, summary)
table(data1$BMI_category, data1$BMI_cat)



#Making a new column based on conditions of another column - mat age
summary(data1$Age)
plot(data1$Age)
data1$Age_category <- "dsjkfjaklfj"
data1$Age_category[data1$Age<20] <- "<20"
data1$Age_category[data1$Age>=20 & data1$Age<25] <- "20-24"
data1$Age_category[data1$Age>=25 & data1$Age<30] <- "25-29"
data1$Age_category[data1$Age>=30 & data1$Age<35] <- "30-34"
data1$Age_category[data1$Age>=35 & data1$Age<40] <- "35-39"
data1$Age_category[data1$Age>=40 & data1$Age<45] <- "40-44"
data1$Age_category[data1$Age>=45] <- "45+"

str(data1$Age_category)
data1$Age_category <- factor(data1$Age_category, levels=c("<20", "20-24", "25-29", "30-34", "35-39", "40-44", "45+"))

tapply(data1$Age, data1$Age_category, summary)



#reordering factor
data1$Smoking <- factor(data1$Smoking, levels=c("current", "ex - gave up in pregnancy","ex - gave up before pregnancy", "never"))


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




#reordering factor levels
table(data1$risk.cat)
data1$risk.cat <- factor(data1$risk.cat , levels=c("Low", "High"))
table(data1$risk.cat)



#Dummy column for when we want to use boxplots or otherwise need it
data1$study <- "INSIGHT"



# Making the OTU data numeric
dim(data1)
str(data1$Otu00001_Lactobacillus_crispatus)
data1[1:10, 3175:3176] #ends at 3175
data1[1:10, 235:237] #starts at 235
data1[235:3175] <- lapply(data1[235:3175], as.numeric)



##changing the dots to NAs in the SNP genotypes
##Check they're all characters
str(data1$rs17333103)
str(data1$rs1983649)
str(data1$rs2664581)
str(data1$rs6032040)
##Dots to NAs
data1 <- data1 %>% 
  mutate(rs17333103 = na_if(rs17333103, ".")) %>%
  mutate(rs1983649 = na_if(rs1983649, ".")) %>%
  mutate(rs2664581 = na_if(rs2664581, ".")) %>%
  mutate(rs6032040 = na_if(rs6032040, "."))

##Changing aesthetics of SNP rs17333103
data1$rs17333103_copy <- data1$rs17333103
table(data1$rs17333103)
data1$rs17333103 <- as.character(data1$rs17333103)
data1$rs17333103[data1$rs17333103 == "c:c"] <- "CC"
data1$rs17333103[data1$rs17333103 == "t:c"] <- "CT"
data1$rs17333103[data1$rs17333103 == "t:t"] <- "TT"
data1$rs17333103 <- factor(data1$rs17333103, levels=c("CC", "CT", "TT"))
table(data1$rs17333103, data1$rs17333103_copy, useNA="always")
data1$rs17333103_copy <- NULL

##Changing aesthetics of SNP rs1983649
data1$rs1983649_copy <- data1$rs1983649
table(data1$rs1983649)
data1$rs1983649 <- as.character(data1$rs1983649)
data1$rs1983649[data1$rs1983649 == "a:a"] <- "AA"
data1$rs1983649[data1$rs1983649 == "t:a"] <- "AT"
data1$rs1983649[data1$rs1983649 == "t:t"] <- "TT"
data1$rs1983649 <- factor(data1$rs1983649, levels=c("TT", "AT", "AA"))
table(data1$rs1983649, data1$rs1983649_copy, useNA="always")
data1$rs1983649_copy <- NULL

##Changing aesthetics of SNP rs2664581
data1$rs2664581_copy <- data1$rs2664581
table(data1$rs2664581)
data1$rs2664581 <- as.character(data1$rs2664581)
data1$rs2664581[data1$rs2664581 == "a:a"] <- "AA"
data1$rs2664581[data1$rs2664581 == "c:a"] <- "AC"
data1$rs2664581[data1$rs2664581 == "c:c"] <- "CC"
data1$rs2664581 <- factor(data1$rs2664581, levels=c("CC", "AC", "AA"))
table(data1$rs2664581, data1$rs2664581_copy, useNA="always")
data1$rs2664581_copy <- NULL

##Changing aesthetics of SNP rs6032040
data1$rs6032040_copy <- data1$rs6032040
table(data1$rs6032040)
data1$rs6032040 <- as.character(data1$rs6032040)
data1$rs6032040[data1$rs6032040 == "a:a"] <- "AA"
data1$rs6032040[data1$rs6032040 == "t:a"] <- "AT"
data1$rs6032040[data1$rs6032040 == "t:t"] <- "TT"
data1$rs6032040 <- factor(data1$rs6032040, levels=c("TT","AT", "AA"))
table(data1$rs6032040, data1$rs6032040_copy, useNA="always")
data1$rs6032040_copy <- NULL






############################################################################################
##                                                                                        ##
##                           Making data_early and data_late                              ##
##                                                                                        ##
############################################################################################

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

#changing row names to Sample_IDs
rownames(data_late) <- data_late$Sample_ID
data_late[1:5, 1:4]


##Checks that both samples/timepoints of the women have the same SNP data
library(arsenal) #https://cran.r-project.org/web/packages/arsenal/vignettes/comparedf.html
SNPs_early <- data_early %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040)
SNPs_late <- data_late %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040)
#Compare
summary(comparedf(SNPs_early, SNPs_late, by="participant_id")) # they all match 

##Checks that both samples/timepoints of the women have the same gestational age data
gestation_early <- data_early %>% 
  select(participant_id, Gestation.at.delivery.wks.dec)
gestation_late <- data_late %>% 
  select(participant_id, Gestation.at.delivery.wks.dec)
#Compare
summary(comparedf(gestation_early, gestation_late, by="participant_id")) # they all match 



############################################################################################
##                                                                                        ##
##                                     Saving Rdata                                       ##
##                                                                                        ##
############################################################################################

###Saving the data from this script
save(data1, data_early, data_late, 
     file="data_BV/main_dfs.RData")


#Loading .Rdata saved 
#lnames = load(file="data_BV/main_dfs.RData")
#lnames



