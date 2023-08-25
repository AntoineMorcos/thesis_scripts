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
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
getwd()
#list.files()


# This R script comes before the LEfSe shell script 


###########################################################################################################################
##                                                                                                                       ##
##                                                  Loading data                                                         ##
##                                                                                                                       ##
###########################################################################################################################

#Loading .Rdata saved 
lnames = load(file="data_BV/dfs_specifically_for_LEfSe.RData")
lnames

# This data came from the script "LEfSe_prep1_dfs_for_early_and_late.R"

###########################################################################################################################
##                                                                                                                       ##
##                                              Dropping samples                                                          ##
##                                                                                                                       ##
###########################################################################################################################

## Drop samples with missing OTU data
dim(data_late) #[1]  323 3180
data_late_filtered <- data_late %>% 
  drop_na("Otu00001_Lactobacillus_crispatus") 
dim(data_late_filtered) #[1]  313 3180


## Drop samples with missing BV status
dim(data_late_filtered)
data_late_filtered <- data_late_filtered %>% 
  drop_na("Bacterial.Vaginosis..categories") 
dim(data_late_filtered) #[1]  304 3180





###########################################################################################################################
##                                                                                                                       ##
##                                  Clinical data needed for LEfSe                                                       ##
##                                                                                                                       ##
###########################################################################################################################

## Checks & getting N numbers 
table(data_late_filtered$Ethnicity_short, data_late_filtered$BV_categories, useNA="always")
#      Normal Intermediate  BV <NA>
#White    165           27  19    0
#Black     40           16   8    0
#Other     22            5   2    0
#<NA>       0            0   0    0


table(data_late_filtered$BV_categories, useNA="always")
#      Normal Intermediate           BV         <NA>  
#         227           48           29            0 


## make df w/ no unnecessary columns
clinical_data_for_lefse <- data_late_filtered %>% 
  select(Sample_ID, BV_categories, Ethnicity_short)

## More checks
clinical_data_for_lefse[1:5,]
table(clinical_data_for_lefse$BV_categories)


###########################################################################################################################
##                                                                                                                       ##
##                                     OTUs into the format like hmp_aerobiosis_small                                    ##
##                                                                                                                       ##
###########################################################################################################################

OTU_list <- read.csv("data_BV/OTU_full_taxa_list.csv", header=TRUE, na.strings=c("NA", " ", ""), fileEncoding = 'UTF-8-BOM')

#Checks
str(OTU_list)

#Making a new column with taxa data added on to OTU names
OTU_list$OTU_with_taxa <- paste(OTU_list$OTU, OTU_list$Species, sep="_")

#Making spaces into underscores
OTU_list$OTU_with_taxa <- str_replace_all(OTU_list$OTU_with_taxa, c(" "="_"))
OTU_list$Species <- str_remove_all(OTU_list$Species, c(" "="_"))

#Checks
str(OTU_list)
names(OTU_list)
table(OTU_list$Kingdom) # all Bacteria

# Add an indication of what classification to make things easier on lefse plots
head(OTU_list$Phylum)
OTU_list$Phylum <- paste("p_", OTU_list$Phylum, sep="")
head(OTU_list$Phylum)
head(OTU_list$Class)
OTU_list$Class <- paste("c_", OTU_list$Class, sep="")
head(OTU_list$Class)
head(OTU_list$Order)
OTU_list$Order <- paste("o_", OTU_list$Order, sep="")
head(OTU_list$Order)
head(OTU_list$Family)
OTU_list$Family <- paste("f_", OTU_list$Family, sep="")
head(OTU_list$Family)
head(OTU_list$Genus)
OTU_list$Genus <- paste("g_", OTU_list$Genus, sep="")
head(OTU_list$Genus)
head(OTU_list$Species)
OTU_list$Species <- paste("s_", OTU_list$Species, sep="")
head(OTU_list$Species)
str(OTU_list)


## Making a column to be used in LEfSe - with species
OTU_list$taxonomy_long_list_for_lefse <- paste(OTU_list$Kingdom, OTU_list$Phylum, OTU_list$Class, OTU_list$Order, 
                                               OTU_list$Family, OTU_list$Genus, OTU_list$Species, OTU_list$OTU_with_taxa,
                                               sep=" | ")

## Checks
head(OTU_list$taxonomy_long_list_for_lefse)


###########################################################################################################################
##                                                                                                                       ##
##                               Relative abundance of OTUs                                                        ##
##                                                                                                                       ##
###########################################################################################################################

## Making OTUs into a relative abundance matrix https://rdrr.io/cran/funrar/man/make_relative.html#heading-3
library(funrar) #install.packages("funrar")

## Raw OTUs
#data_late_filtered[1:10, 3175:3176] #ends at 3175
#data_late_filtered[1:10, 235:237] #starts at 235
raw_OTUs_late <- data_late_filtered[,c(235:3175)]

# checks
dim(raw_OTUs_late)
names(raw_OTUs_late[1:6])
names(raw_OTUs_late[2936:2941])
raw_OTUs_late[1:5,1:5]

# changing the structure to a matrix
raw_OTUs_late_matrix <- as.matrix(raw_OTUs_late)
dim(raw_OTUs_late_matrix)
raw_OTUs_late_matrix[1:5,1:5]

# make relative abundance matrix command
AM_OTUs_late <- make_relative(raw_OTUs_late_matrix)
# preview abundance matrix
AM_OTUs_late[1:5,1:5]

## Making AM into a df
AM_OTUs_late_df <- as.data.frame(AM_OTUs_late)
AM_OTUs_late_df[1:5,1:5]






###########################################################################################################################
##                                                                                                                       ##
##                          Replacing OTU names with the long ones so LEfSe will work                                    ##
##                                                                                                                       ##
###########################################################################################################################

##Swapping old OTU labels to new long LEfSe appropriate ones
##Transposing the df
AM_OTUs_late_df_t_matrix <- t(AM_OTUs_late_df)
AM_OTUs_late_df_t_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
AM_OTUs_late_df_t <- as.data.frame(AM_OTUs_late_df_t_matrix)
str(AM_OTUs_late_df_t)
AM_OTUs_late_df_t[1:5,1:5]

## Adding column from OTU_list
AM_OTUs_late_df_t$taxonomy_long_list_for_lefse <- OTU_list$taxonomy_long_list_for_lefse

## Move taxonomy_long_list_for_lefse to the beginning of the df
AM_OTUs_late_df_t <- AM_OTUs_late_df_t %>% 
  select(taxonomy_long_list_for_lefse, everything())

## Checks
AM_OTUs_late_df_t[1:5,1:5]
str(AM_OTUs_late_df_t)



########################################################## Checking OTUs are identical
#write to csv to check they are the same (in excel)
write.csv(AM_OTUs_late_df_t, "data_BV/LEfSe/AM_OTUs_late_df_t.csv")
#looks fine


####Checking OTUs are identical (in R)......
## Making a row ID column for merging with clinical data
AM_OTUs_late_df_t$OTU_Flavia_label <- rownames(AM_OTUs_late_df_t)

## Move ID to the beginning of the df
AM_OTUs_late_df_t <- AM_OTUs_late_df_t %>% 
  select(OTU_Flavia_label, everything())

##Keeping only the Otu part of the name of the long LEfSe label
AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number <- gsub(".*(Otu)", "\\1", AM_OTUs_late_df_t$taxonomy_long_list_for_lefse)
head(AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number)

##Test if the OTU columns are identical 
identical(AM_OTUs_late_df_t[['OTU_Flavia_label']], 
          AM_OTUs_late_df_t[['taxonomy_long_list_for_lefse__OTU_number']]) #different...
AM_OTUs_late_df_t$OTU_Flavia_label[!(AM_OTUs_late_df_t$OTU_Flavia_label %in% AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number)]
AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number[!(AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number %in% AM_OTUs_late_df_t$OTU_Flavia_label)]
# looks like the only diff is in Flavia's label she uses .s not -s 

## Move ID to the beginning of the df
AM_OTUs_late_df_t <- AM_OTUs_late_df_t %>% 
  select(taxonomy_long_list_for_lefse, taxonomy_long_list_for_lefse__OTU_number, OTU_Flavia_label,
         everything())

#write to csv to view in excel
write_csv(AM_OTUs_late_df_t, "data_BV/LEfSe/AM_OTUs_late_df_t_many_OTU_cols.csv")
########################################################## 

# OTUs match so now we can delete non-LEfSe ones
AM_OTUs_late_df_t$taxonomy_long_list_for_lefse__OTU_number <- NULL
AM_OTUs_late_df_t$OTU_Flavia_label <- NULL
AM_OTUs_late_df_t$OTU_Flavia_label <- NULL

# Change rownames to lefse ones
rownames(AM_OTUs_late_df_t) <- AM_OTUs_late_df_t$taxonomy_long_list_for_lefse
#View(AM_OTUs_late_df_t)

#Now delete final OTU name ready for transposing
AM_OTUs_late_df_t$taxonomy_long_list_for_lefse <- NULL
#View(AM_OTUs_late_df_t)

##Transposing the df back again
AM_OTUs_late_df_final_matrix<- t(AM_OTUs_late_df_t)
AM_OTUs_late_df_final_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
AM_OTUs_late_df_final <- as.data.frame(AM_OTUs_late_df_final_matrix)
str(AM_OTUs_late_df_final)
AM_OTUs_late_df_final[1:5,1:5]



## Making a row ID column for merging with clinical data
AM_OTUs_late_df_final$Sample_ID <- rownames(AM_OTUs_late_df_final)

## Move ID to the beginning of the df
AM_OTUs_late_df_final <- AM_OTUs_late_df_final %>% 
  select(Sample_ID, everything())

## Checks
AM_OTUs_late_df_final[1:5,1:5]
#View(AM_OTUs_late_df_final)




###########################################################################################################################
##                                                                                                                       ##
##                                      Merging and getting file ready for LEfSe                                         ##
##                                                                                                                       ##
###########################################################################################################################

## Checks
clinical_data_for_lefse[1:5,]
AM_OTUs_late_df_final[1:5,1:5]

dim(clinical_data_for_lefse)
dim(AM_OTUs_late_df_final)

str(clinical_data_for_lefse)
str(AM_OTUs_late_df_final)

## late df merge - clinical BV and AM OTU data
##Merge
lefse_late_df_unordered_t <- merge(clinical_data_for_lefse, AM_OTUs_late_df_final, by="Sample_ID")

##Moving the order of columns t
lefse_late_df_t <- lefse_late_df_unordered_t %>% 
  select(BV_categories, Ethnicity_short, everything())

lefse_late_df_t[1:5,1:5]
#View(lefse_late_df_t)
############### I start from this point in the below section when I don't look at samples with an intermediate status


##Transposing the df
lefse_late_df_matrix <- t(lefse_late_df_t)
lefse_late_df_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df <- as.data.frame(lefse_late_df_matrix)
str(lefse_late_df)
lefse_late_df[1:5,1:5]
#View(lefse_late_df)

##Moving the data stored in the rownames into an actual column
lefse_late_df$row_labels <- rownames(lefse_late_df)

##Moving the column row_labels to the beginning 
lefse_late_df_final <- lefse_late_df %>% 
  select(row_labels, everything())

lefse_late_df_final[1:5,1:5]
#View(lefse_late_df_final)

# Looking for N numbers for lefse plots - late
table(lefse_late_df_t$BV_categories, useNA = "always")
# Normal Intermediate           BV         <NA> 
#   227           48           29            0 

table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
#               White Black Other <NA>
# Normal         165    40    22    0
# Intermediate    27    16     5    0
# BV              19     8     2    0
# <NA>             0     0     0    0






###########################################################################################################################
##                                                                                                                       ##
##                                 Working on BV vs Normal (Intermediates are gone)                                      ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering out all intermediates (for just BV vs Normal analysis)
lefse_late_df_t_BV_vs_normal <- lefse_late_df_t %>% 
  filter(BV_categories!="Intermediate")

# Checks:
#View(lefse_late_df_t_BV_vs_normal)
table(lefse_late_df_t_BV_vs_normal$BV_categories, useNA = "always")
dim(lefse_late_df_t_BV_vs_normal) # 256 2946

##############################################################


##Transposing the df
lefse_late_df_BV_vs_normal_matrix <- t(lefse_late_df_t_BV_vs_normal)
lefse_late_df_BV_vs_normal_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_BV_vs_normal <- as.data.frame(lefse_late_df_BV_vs_normal_matrix)
str(lefse_late_df_BV_vs_normal)
lefse_late_df_BV_vs_normal[1:5,1:5]
#View(lefse_late_df_BV_vs_normal)

##Moving the data stored in the rownames into an actual column
lefse_late_df_BV_vs_normal$row_labels <- rownames(lefse_late_df_BV_vs_normal)

##Moving the column row_labels to the beginning 
lefse_late_df_BV_vs_normal_final <- lefse_late_df_BV_vs_normal %>% 
  select(row_labels, everything())

lefse_late_df_BV_vs_normal_final[1:5,1:5]

#View(lefse_late_df_BV_vs_normal_final)

###df with only "BV_categories", Sample_ID  and OTU data
lefse_late_df_final_BV_vs_normal <- lefse_late_df_BV_vs_normal_final[-c(1,2,4), ] #delete row 1&2&4




###########################################################################################################################
##                                                                                                                       ##
##                                          Working on Int vs Normal (BV are gone)                                       ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering out all BV (for just intermediates vs normal analysis)
lefse_late_df_t_int_vs_normal <- lefse_late_df_t %>% 
  filter(BV_categories!="BV")

# Checks:
#View(lefse_late_df_t_int_vs_normal)
table(lefse_late_df_t_int_vs_normal$BV_categories, useNA = "always")
dim(lefse_late_df_t_int_vs_normal) # 275 2946

##############################################################


##Transposing the df
lefse_late_df_int_vs_normal_matrix <- t(lefse_late_df_t_int_vs_normal)
lefse_late_df_int_vs_normal_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_int_vs_normal <- as.data.frame(lefse_late_df_int_vs_normal_matrix)
str(lefse_late_df_int_vs_normal)
lefse_late_df_int_vs_normal[1:5,1:5]
#View(lefse_late_df_int_vs_normal)

##Moving the data stored in the rownames into an actual column
lefse_late_df_int_vs_normal$row_labels <- rownames(lefse_late_df_int_vs_normal)

##Moving the column row_labels to the beginning 
lefse_late_df_int_vs_normal_final <- lefse_late_df_int_vs_normal %>% 
  select(row_labels, everything())

lefse_late_df_int_vs_normal_final[1:5,1:5]
#View(lefse_late_df_int_vs_normal_final)

###df with only "BV_categories", Sample_ID  and OTU data
lefse_late_df_final_int_vs_normal <- lefse_late_df_int_vs_normal_final[-c(1,2,4), ] #delete row 1&2&4



###########################################################################################################################
##                                                                                                                       ##
##                                 Working on BV vs Intermediates (Normals are gone)                                      ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering out all intermediates (for just BV vs Normal analysis)
lefse_late_df_t_BV_vs_int <- lefse_late_df_t %>% 
  filter(BV_categories!="Normal")

# Checks:
#View(lefse_late_df_t_BV_vs_int)
table(lefse_late_df_t_BV_vs_int$BV_categories, useNA = "always")
dim(lefse_late_df_t_BV_vs_int) # 77 2946

##############################################################


##Transposing the df
lefse_late_df_BV_vs_int_matrix <- t(lefse_late_df_t_BV_vs_int)
lefse_late_df_BV_vs_int_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_BV_vs_int <- as.data.frame(lefse_late_df_BV_vs_int_matrix)
str(lefse_late_df_BV_vs_int)
lefse_late_df_BV_vs_int[1:5,1:5]
#View(lefse_late_df_BV_vs_int)

##Moving the data stored in the rownames into an actual column
lefse_late_df_BV_vs_int$row_labels <- rownames(lefse_late_df_BV_vs_int)


##Moving the column row_labels to the beginning 
lefse_late_df_BV_vs_int_final <- lefse_late_df_BV_vs_int %>% 
  select(row_labels, everything())

lefse_late_df_BV_vs_int_final[1:5,1:5]
#View(lefse_late_df_BV_vs_int_final)

###df with only "BV_categories", Sample_ID  and OTU data
lefse_late_df_final_BV_vs_int <- lefse_late_df_BV_vs_int_final[-c(1,2,4), ] #delete row 1&2&4





###########################################################################################################################
##                                                                                                                       ##
##                      White women - Working on BV vs Normal (Intermediates are gone)                                   ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering out all intermediates (for just BV vs Normal analysis)
lefse_late_df_t_White_BV_vs_normal <- lefse_late_df_t %>% 
  filter(BV_categories!="Intermediate") %>% 
  filter(Ethnicity_short=="White")

# Checks:
#View(lefse_late_df_t_White_BV_vs_normal)
table(lefse_late_df_t_White_BV_vs_normal$BV_categories, lefse_late_df_t_White_BV_vs_normal$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t_White_BV_vs_normal) # 184 2946

##############################################################


##Transposing the df
lefse_late_df_White_BV_vs_normal_matrix <- t(lefse_late_df_t_White_BV_vs_normal)
lefse_late_df_White_BV_vs_normal_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_White_BV_vs_normal <- as.data.frame(lefse_late_df_White_BV_vs_normal_matrix)
str(lefse_late_df_White_BV_vs_normal)
lefse_late_df_White_BV_vs_normal[1:5,1:5]
#View(lefse_late_df_White_BV_vs_normal)

##Moving the data stored in the rownames into an actual column
lefse_late_df_White_BV_vs_normal$row_labels <- rownames(lefse_late_df_White_BV_vs_normal)

##Moving the column row_labels to the beginning 
lefse_late_df_White_BV_vs_normal_final <- lefse_late_df_White_BV_vs_normal %>% 
  select(row_labels, everything())

lefse_late_df_White_BV_vs_normal_final[1:5,1:5]
#View(lefse_late_df_White_BV_vs_normal_final)

###df with only "BV_categories", Sample_ID  and OTU data
lefse_late_df_final_White_BV_vs_normal <- lefse_late_df_White_BV_vs_normal_final[-c(1,2,4), ] #delete row 1&2&4
#View(lefse_late_df_final_White_BV_vs_normal)








###########################################################################################################################
##                                                                                                                       ##
##                      Black women - Working on BV vs Normal (Intermediates are gone)                                   ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering out all intermediates (for just BV vs Normal analysis)
lefse_late_df_t_Black_BV_vs_normal <- lefse_late_df_t %>% 
  filter(BV_categories!="Intermediate") %>% 
  filter(Ethnicity_short=="Black")

# Checks:
#View(lefse_late_df_t_Black_BV_vs_normal)
table(lefse_late_df_t_Black_BV_vs_normal$BV_categories, lefse_late_df_t_Black_BV_vs_normal$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t_Black_BV_vs_normal) # 48 2946

##############################################################


##Transposing the df
lefse_late_df_Black_BV_vs_normal_matrix <- t(lefse_late_df_t_Black_BV_vs_normal)
lefse_late_df_Black_BV_vs_normal_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_Black_BV_vs_normal <- as.data.frame(lefse_late_df_Black_BV_vs_normal_matrix)
str(lefse_late_df_Black_BV_vs_normal)
lefse_late_df_Black_BV_vs_normal[1:5,1:5]
#View(lefse_late_df_Black_BV_vs_normal)

##Moving the data stored in the rownames into an actual column
lefse_late_df_Black_BV_vs_normal$row_labels <- rownames(lefse_late_df_Black_BV_vs_normal)

##Moving the column row_labels to the beginning 
lefse_late_df_Black_BV_vs_normal_final <- lefse_late_df_Black_BV_vs_normal %>% 
  select(row_labels, everything())

lefse_late_df_Black_BV_vs_normal_final[1:5,1:5]
#View(lefse_late_df_Black_BV_vs_normal_final)

###df with only "BV_categories", Sample_ID  and OTU data
lefse_late_df_final_Black_BV_vs_normal <- lefse_late_df_Black_BV_vs_normal_final[-c(1,2,4), ] #delete row 1&2&4
#View(lefse_late_df_final_Black_BV_vs_normal)






###########################################################################################################################
##                                                                                                                       ##
##                                      Working on ethnicity impact in women w/ BV                                         ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering 
lefse_late_df_t_ethnicity_in_BV <- lefse_late_df_t %>% 
  filter(BV_categories=="BV") %>% 
  filter(Ethnicity_short=="White" | Ethnicity_short=="Black")


# Checks:
#View(lefse_late_df_t_ethnicity_in_BV)
table(lefse_late_df_t_ethnicity_in_BV$BV_categories, lefse_late_df_t_ethnicity_in_BV$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t_ethnicity_in_BV) # 27 2946

##############################################################


##Transposing the df
lefse_late_df_ethnicity_in_BV_matrix <- t(lefse_late_df_t_ethnicity_in_BV)
lefse_late_df_ethnicity_in_BV_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_ethnicity_in_BV <- as.data.frame(lefse_late_df_ethnicity_in_BV_matrix)
str(lefse_late_df_ethnicity_in_BV)
lefse_late_df_ethnicity_in_BV[1:5,1:5]
#View(lefse_late_df_ethnicity_in_BV)

##Moving the data stored in the rownames into an actual column
lefse_late_df_ethnicity_in_BV$row_labels <- rownames(lefse_late_df_ethnicity_in_BV)

##Moving the column row_labels to the beginning 
lefse_late_df_ethnicity_in_BV_final <- lefse_late_df_ethnicity_in_BV %>% 
  select(row_labels, everything())

lefse_late_df_ethnicity_in_BV_final[1:5,1:5]
#View(lefse_late_df_ethnicity_in_BV_final)

###df with only "Ethnicity_short", Sample_ID  and OTU data
lefse_late_df_final_ethnicity_in_BV <- lefse_late_df_ethnicity_in_BV_final[-c(1,2,3), ] #delete row 1&2&3
#View(lefse_late_df_final_ethnicity_in_BV)





###########################################################################################################################
##                                                                                                                       ##
##                                      Working on ethnicity impact in women w/ Intermediate                             ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering 
lefse_late_df_t_ethnicity_in_Intermediate <- lefse_late_df_t %>% 
  filter(BV_categories=="Intermediate") %>% 
  filter(Ethnicity_short=="White" | Ethnicity_short=="Black")


# Checks:
#View(lefse_late_df_t_ethnicity_in_Intermediate)
table(lefse_late_df_t_ethnicity_in_Intermediate$BV_categories, lefse_late_df_t_ethnicity_in_Intermediate$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t_ethnicity_in_Intermediate) # 43 2946

##############################################################


##Transposing the df
lefse_late_df_ethnicity_in_Intermediate_matrix <- t(lefse_late_df_t_ethnicity_in_Intermediate)
lefse_late_df_ethnicity_in_Intermediate_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_ethnicity_in_Intermediate <- as.data.frame(lefse_late_df_ethnicity_in_Intermediate_matrix)
str(lefse_late_df_ethnicity_in_Intermediate)
lefse_late_df_ethnicity_in_Intermediate[1:5,1:5]
#View(lefse_late_df_ethnicity_in_Intermediate)

##Moving the data stored in the rownames into an actual column
lefse_late_df_ethnicity_in_Intermediate$row_labels <- rownames(lefse_late_df_ethnicity_in_Intermediate)

##Moving the column row_labels to the beginning 
lefse_late_df_ethnicity_in_Intermediate_final <- lefse_late_df_ethnicity_in_Intermediate %>% 
  select(row_labels, everything())

lefse_late_df_ethnicity_in_Intermediate_final[1:5,1:5]
#View(lefse_late_df_ethnicity_in_Intermediate_final)

###df with only "Ethnicity_short", Sample_ID  and OTU data
lefse_late_df_final_ethnicity_in_Intermediate <- lefse_late_df_ethnicity_in_Intermediate_final[-c(1,2,3), ] #delete row 1&2&3
#View(lefse_late_df_final_ethnicity_in_Intermediate)







###########################################################################################################################
##                                                                                                                       ##
##                                      Working on ethnicity impact in women w/ Normal                                   ##
##                                                                                                                       ##
###########################################################################################################################

#Starting with this df made in previous section .....
lefse_late_df_t[1:5,1:5]

############################################################## filtering
# Getting stats
#View(lefse_late_df_t)
table(lefse_late_df_t$BV_categories, lefse_late_df_t$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t) # 304 2946

####Filtering 
lefse_late_df_t_ethnicity_in_Normal <- lefse_late_df_t %>% 
  filter(BV_categories=="Normal") %>% 
  filter(Ethnicity_short=="White" | Ethnicity_short=="Black")


# Checks:
#View(lefse_late_df_t_ethnicity_in_Normal)
table(lefse_late_df_t_ethnicity_in_Normal$BV_categories, lefse_late_df_t_ethnicity_in_Normal$Ethnicity_short, useNA = "always")
dim(lefse_late_df_t_ethnicity_in_Normal) # 205 2946

##############################################################


##Transposing the df
lefse_late_df_ethnicity_in_Normal_matrix <- t(lefse_late_df_t_ethnicity_in_Normal)
lefse_late_df_ethnicity_in_Normal_matrix[1:5,1:5]

##Making it to a df as transposing made it a matrix somehow
lefse_late_df_ethnicity_in_Normal <- as.data.frame(lefse_late_df_ethnicity_in_Normal_matrix)
str(lefse_late_df_ethnicity_in_Normal)
lefse_late_df_ethnicity_in_Normal[1:5,1:5]
#View(lefse_late_df_ethnicity_in_Normal)

##Moving the data stored in the rownames into an actual column
lefse_late_df_ethnicity_in_Normal$row_labels <- rownames(lefse_late_df_ethnicity_in_Normal)

##Moving the column row_labels to the beginning 
lefse_late_df_ethnicity_in_Normal_final <- lefse_late_df_ethnicity_in_Normal %>% 
  select(row_labels, everything())

lefse_late_df_ethnicity_in_Normal_final[1:5,1:5]
#View(lefse_late_df_ethnicity_in_Normal_final)

###df with only "Ethnicity_short", Sample_ID  and OTU data
lefse_late_df_final_ethnicity_in_Normal <- lefse_late_df_ethnicity_in_Normal_final[-c(1,2,3), ] #delete row 1&2&3
#View(lefse_late_df_final_ethnicity_in_Normal)







###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting files                                                        ##
##                                                                                                                       ##
###########################################################################################################################


###### DFs where we are doing a pairwise comparison of the 3 BV categories (& therefore we are excluding 1 category per TXT)

##Writing out the BV_vs_normal df as a txt without quotes
write.table(lefse_late_df_final_BV_vs_normal,"data_BV/LEfSe/BV_status___late__BV_vs_normal.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)

##Writing out the int_vs_normal df as a txt without quotes
write.table(lefse_late_df_final_int_vs_normal,"data_BV/LEfSe/BV_status___late__int_vs_normal.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)

##Writing out the BV_vs_int df as a txt without quotes
write.table(lefse_late_df_final_BV_vs_int,"data_BV/LEfSe/BV_status___late__BV_vs_int.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)



###### BV vs normal, stratified by ethnicity

##Writing out the Black_BV_vs_normal df as a txt without quotes
write.table(lefse_late_df_final_Black_BV_vs_normal,"data_BV/LEfSe/BV_status___late__BV_vs_normal__in_Black_women.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)

##Writing out the White_BV_vs_normal df as a txt without quotes
write.table(lefse_late_df_final_White_BV_vs_normal,"data_BV/LEfSe/BV_status___late__BV_vs_normal__in_White_women.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)



###### Ethnicity impact, stratified by BV status

##Writing out the ethnicity_in_BV df as a txt without quotes
write.table(lefse_late_df_final_ethnicity_in_BV,"data_BV/LEfSe/Ethnicity__late__women_w_BV.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)

##Writing out the ethnicity_in_Intermediate df as a txt without quotes
write.table(lefse_late_df_final_ethnicity_in_Intermediate,"data_BV/LEfSe/Ethnicity__late__women_w_Intermediate.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)

##Writing out the ethnicity_in_Normal df as a txt without quotes
write.table(lefse_late_df_final_ethnicity_in_Normal,"data_BV/LEfSe/Ethnicity__late__women_w_Normal.txt", 
            sep="\t", row.names=F, col.names=F, quote=F)











