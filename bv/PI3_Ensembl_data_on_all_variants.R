#clear environment
rm(list=ls())

getwd()
#list.files()

library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(ggplot2)
library(dplyr)
library(tidyr)


##Loading .Rdata saved 
setwd("~/OneDrive - King's College London/PhD/Projects/SNPs") #mac
#setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/SNPs") #pc
lnames





############################################################################################
##                                                                                        ##
##                  Looking at SNPs info from Ensembl                                     ##
##                                                                                        ##
############################################################################################

#data from http://www.ensembl.org/Homo_sapiens/Gene/Variation_Gene/Table?db=core;g=ENSG00000124102;r=20:45174902-45176544;t=ENST00000243924

#loading in df
PI3_variants <- read.csv("data_SNPs/PI3_variant_table_exported_from_Ensembl.csv", header=T, na.strings=c("NA", " ", ""))

#Numbers of different types of variants
table(PI3_variants$Class) #481 SNPs

#df of just SNPs (no other variant types)
PI3_SNPs <- PI3_variants %>% 
  filter(Class=="SNP")

#df of SNPs of interest
PI3_variant_table_exported_from_Ensembl___just_4_SNPs_we_looked_at <- PI3_SNPs %>% 
  filter(Variant.ID=="rs17333103" | Variant.ID=="rs1983649" | Variant.ID=="rs6032040" | Variant.ID=="rs2664581" )
write_csv(PI3_variant_table_exported_from_Ensembl___just_4_SNPs_we_looked_at, "data_SNPs/PI3_variant_table_exported_from_Ensembl___just_4_SNPs_we_looked_at.csv")

#There's multiple rows for some SNPs 
#(e.g. here there's 3 rows for rs1983649)
PI3_SNPs %>% 
  filter(Variant.ID=="rs1983649")

# filter for 1 row per SNP
PI3_SNPs_unq <- PI3_SNPs %>%
  distinct(Variant.ID, .keep_all=T) 

#check for duplicates
PI3_SNPs_unq %>% 
  filter(Variant.ID=="rs1983649")

#how many SNPs are there in this gene?
dim(PI3_SNPs_unq) #435 SNPs










