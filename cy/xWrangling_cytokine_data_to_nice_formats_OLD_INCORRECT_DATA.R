###########################################################################################################################
##                                                                                                                       ##
##                                               Load data & df prep                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Clear environment
rm(list=ls())

## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()


# Downloaded these csvs from DATA folder in Amirah's Luminex folder
list.files("data/from_Amirah/Cytokines_QuiteRaw/")

library(tidyverse)




######################################################### clean plate1
# plate1 - unclean file
?read.csv
plate1_unclean <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/Amirah Cytokine plate1 140422.csv", 
                           skip=364, header=T) #"DataType:	Result" is line 364
#View(plate1_unclean)
names(plate1_unclean)
plate1 <- plate1_unclean[1:96,1:13] #96 is the last line of the "DataType:	Result" section
plate1_unclean[96:97,] #verify last line
tail(plate1)
head(plate1)

# Getting rid of "." in colnames to match the format of plate2
names(plate1)
names(plate1) <- sub(".", "", names(plate1),  fixed = TRUE)
names(plate1)
#########################################################




######################################################### clean plate2
# plate2 - file has been cleaned-ish
plate2 <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/Amirah Cytokine Plate 2 200422.csv") #whole spreadsheet is only "DataType:	Result"
plate2 <- plate2 %>% 
  dplyr::rename(Sample=X)
tail(plate2)
head(plate2)
#########################################################




######################################################### clean plate3
# plate3 - unclean file
plate3_unclean <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/Amirah Cytokine Plate 3 210422.csv", 
                           skip=370, header=T) #DataType:	Result is line 370
#View(plate3_unclean)
plate3 <- plate3_unclean[1:96,] #96 is the last line of the "DataType:	Result" section
plate3_unclean[96:97,] #verify last line
tail(plate3)
head(plate3)


# Getting rid of "." in colnames to match the format of plate2
names(plate3)
names(plate3) <- sub(".", "", names(plate3),  fixed = TRUE)
names(plate3)
#########################################################




######################################################### Adding plate column
# Adding plate column
plate1$plate <- "plate1"
plate2$plate <- "plate2"
plate3$plate <- "plate3"
#########################################################






###########################################################################################################################
##                                                                                                                       ##
##                                        Combine plates 1-3 (in long format)                                            ##
##                                                                                                                       ##
###########################################################################################################################

# Seeing if colnames match
names(plate1)
names(plate2)
names(plate3)


# Options
?merge # temp1 <- merge(plate1, plate2, by="Sample", all=T) #tricky as duplicates so will try another method. Also Sample numbers are not the same IDs across plates.
?rbind


# Merging the dfs
merged_long <- data.table::rbindlist(list(plate1, plate2, plate3), #bind all plates under same columns
                                     use.names=T, #use.names	    TRUE binds by matching column name
                                     fill=T)      #fill	          TRUE fills missing columns with NAs

# From data table to df
head(merged_long)
merged_long_df <- as.data.frame(merged_long)
head(merged_long_df)


# Sort df
head(merged_long_df)
merged_long_sorted <- merged_long_df %>%
  dplyr::arrange(Sample, plate) %>% 
  dplyr::select(Sample, plate, everything())
head(merged_long_sorted)




# Checks merge worked correctly 
merged_long_sorted %>% dplyr::filter(Sample == "Unknown9") 
plate1 %>% dplyr::filter(Sample == "Unknown9") %>%  dplyr::select(Sample, plate, everything())
plate2 %>% dplyr::filter(Sample == "Unknown9") %>%  dplyr::select(Sample, plate, everything())
plate3 %>% dplyr::filter(Sample == "Unknown9") %>%  dplyr::select(Sample, plate, everything())


# Exporting
write_csv(merged_long_sorted, file = "data/cytokines_exports/cytokine_results_longformat_plates123.csv")







###########################################################################################################################
##                                                                                                                       ##
##               Reshape df (in wide format) for the measurements taken in duplicate, for each plate                     ##
##                                                                                                                       ##
###########################################################################################################################

######################################################### plate1
# Ordering dfs
head(plate1)
plate1 <- plate1 %>% 
  dplyr::select(Sample, plate, everything())
head(plate1)

# Adding kind of a dummy variable of measurement 1 & 2 (as each measurement was taken in duplicate)
plate1$Sample
table(plate1$Sample)
plate1$Measurement_number <- as.character(rep(c(1,2)))
plate1 %>% dplyr::select(Sample, Measurement_number)

# Making the format wide for this df str 
?stats::reshape
plate1_wide <- stats::reshape(plate1, idvar=c("Sample", "plate"), timevar = "Measurement_number", direction="wide")

# Importing metadata for plate1
plate1_IDs <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/plate1_IDs.csv", header=T) # this is the plate1 tab from "Luminex_cytokine_run.xlsx"
names(plate1_IDs)
plate1_IDs$Luminex.sample.number

# Making IDs for merge
plate1_wide$Luminex.sample.number <- gsub("Unknown", "", plate1_wide$Sample)

# Check pre-merge
str(plate1_wide$Luminex.sample.number)
str(plate1_IDs$Luminex.sample.number)
plate1_IDs$Luminex.sample.number <- as.character(plate1_IDs$Luminex.sample.number)
str(plate1_IDs$Luminex.sample.number)

# Merging plate1 data with its metadata
?merge 
plate1_ALLinfo <- merge(plate1_IDs, plate1_wide, by="Luminex.sample.number", all=T) 
names(plate1_ALLinfo)

# Ordering cols in dfs
head(plate1_ALLinfo)
plate1_ALLinfo_sorted <- plate1_ALLinfo %>% 
  dplyr::select(insight_id, Luminex.sample.number, Sample, plate, everything()) %>% 
  dplyr::arrange(insight_id)
head(plate1_ALLinfo_sorted)

# Checking data is consistent
plate1_ALLinfo_sorted %>% dplyr::filter(Sample == "Unknown9")
plate1 %>% dplyr::filter(Sample == "Unknown9")
plate1_IDs %>% dplyr::filter(Luminex.sample.number == "9")
# data matches 

# Exporting
write_csv(plate1_ALLinfo_sorted, file = "data/cytokines_exports/separate_plates/cytokine_results_wideformat_plate1.csv")
#########################################################




######################################################### plate2
# Ordering dfs
head(plate2)
plate2 <- plate2 %>% 
  dplyr::select(Sample, plate, everything())
head(plate2)

# Adding kind of a dummy variable of measurement 1 & 2 (as each measurement was taken in duplicate)
plate2$Sample
table(plate2$Sample)
plate2$Measurement_number <- as.character(rep(c(1,2)))
plate2 %>% dplyr::select(Sample, Measurement_number)

# Making the format wide for this df str 
?stats::reshape
plate2_wide <- stats::reshape(plate2, idvar=c("Sample", "plate"), timevar = "Measurement_number", direction="wide")

# Importing metadata for plate2
plate2_IDs <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/plate2_IDs.csv", header=T) # this is the plate2 tab from "Luminex_cytokine_run.xlsx"
names(plate2_IDs)
plate2_IDs$Luminex.sample.number

# Making IDs for merge
plate2_wide$Luminex.sample.number <- gsub("Unknown", "", plate2_wide$Sample)

# Check pre-merge
str(plate2_wide$Luminex.sample.number)
str(plate2_IDs$Luminex.sample.number)
plate2_IDs$Luminex.sample.number <- as.character(plate2_IDs$Luminex.sample.number)
str(plate2_IDs$Luminex.sample.number)

# Merging plate2 data with its metadata
?merge 
plate2_ALLinfo <- merge(plate2_IDs, plate2_wide, by="Luminex.sample.number", all=T) 
names(plate2_ALLinfo)

# Ordering cols in dfs
head(plate2_ALLinfo)
plate2_ALLinfo_sorted <- plate2_ALLinfo %>% 
  dplyr::select(insight_id, Luminex.sample.number, Sample, plate, everything()) %>% 
  dplyr::arrange(insight_id)
head(plate2_ALLinfo_sorted)

# Checking data is consistent
plate2_ALLinfo_sorted %>% dplyr::filter(Sample == "Unknown9")
plate2 %>% dplyr::filter(Sample == "Unknown9")
plate2_IDs %>% dplyr::filter(Luminex.sample.number == "9")
# data matches 

# Exporting
write_csv(plate2_ALLinfo_sorted, file = "data/cytokines_exports/separate_plates/cytokine_results_wideformat_plate2.csv")
#########################################################



######################################################### plate3
# Ordering dfs
head(plate3)
plate3 <- plate3 %>% 
  dplyr::select(Sample, plate, everything())
head(plate3)

# Adding kind of a dummy variable of measurement 1 & 2 (as each measurement was taken in duplicate)
plate3$Sample
table(plate3$Sample)
plate3$Measurement_number <- as.character(rep(c(1,2)))
plate3 %>% dplyr::select(Sample, Measurement_number)

# Making the format wide for this df str 
?stats::reshape
plate3_wide <- stats::reshape(plate3, idvar=c("Sample", "plate"), timevar = "Measurement_number", direction="wide")

# Importing metadata for plate3
plate3_IDs <- read.csv("data/from_Amirah/Cytokines_QuiteRaw/plate3_IDs.csv", header=T) # this is the plate3 tab from "Luminex_cytokine_run.xlsx"
names(plate3_IDs)
plate3_IDs$Luminex.sample.number

# Making IDs for merge
plate3_wide$Luminex.sample.number <- gsub("Unknown", "", plate3_wide$Sample)

# Check pre-merge
str(plate3_wide$Luminex.sample.number)
str(plate3_IDs$Luminex.sample.number)
plate3_IDs$Luminex.sample.number <- as.character(plate3_IDs$Luminex.sample.number)
str(plate3_IDs$Luminex.sample.number)

# Merging plate3 data with its metadata
?merge 
plate3_ALLinfo <- merge(plate3_IDs, plate3_wide, by="Luminex.sample.number", all=T) 
names(plate3_ALLinfo)

# Ordering cols in dfs
head(plate3_ALLinfo)
plate3_ALLinfo_sorted <- plate3_ALLinfo %>% 
  dplyr::select(insight_id, Luminex.sample.number, Sample, plate, everything()) %>% 
  dplyr::arrange(insight_id)
head(plate3_ALLinfo_sorted)

# Checking data is consistent
plate3_ALLinfo_sorted %>% dplyr::filter(Sample == "Unknown9")
plate3 %>% dplyr::filter(Sample == "Unknown9")
plate3_IDs %>% dplyr::filter(Luminex.sample.number == "9")
# data matches 

# Exporting
write_csv(plate3_ALLinfo_sorted, file = "data/cytokines_exports/separate_plates/cytokine_results_wideformat_plate3.csv")
#########################################################









###########################################################################################################################
##                                                                                                                       ##
##                                         Combine plates 1-3 (in wide format)                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Looking
dim(plate1_ALLinfo_sorted)
dim(plate2_ALLinfo_sorted) #this has less columns as "Amirah Cytokine Plate 2 200422.csv" had less cols than the others
dim(plate3_ALLinfo_sorted)
names(plate1_ALLinfo_sorted)
names(plate2_ALLinfo_sorted)
names(plate3_ALLinfo_sorted)


# Options
?rbind
?data.table::rbindlist

# Merging the dfs
merged_wide <- data.table::rbindlist(list(plate1_ALLinfo_sorted, plate2_ALLinfo_sorted, plate3_ALLinfo_sorted), #bind all plates under same columns
                                     use.names=T, #use.names	    TRUE binds by matching column name
                                     fill=T)      #fill	          TRUE fills missing columns with NAs

# Looking
dim(merged_wide)
names(merged_wide)

# From data table to df
head(merged_wide)
merged_wide_df <- as.data.frame(merged_wide)
head(merged_wide_df)
table(merged_wide_df$insight_id)

# Sort df
head(merged_wide_df)
merged_wide_sorted <- merged_wide_df %>%
  dplyr::arrange(insight_id) %>% 
  dplyr::select(insight_id, Sample, plate, everything())
head(merged_wide_sorted)
names(merged_wide_sorted)


# Checks merge worked correctly 
merged_wide_sorted %>% dplyr::filter(Sample == "Unknown9") 
merged_long_sorted %>% dplyr::filter(Sample == "Unknown9") 
# All matching 



# Exporting
write_csv(merged_wide_sorted, file = "data/cytokines_exports/cytokine_results_wideformat_plates123.csv")


