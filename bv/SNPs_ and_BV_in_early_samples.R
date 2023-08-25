#clear environment
rm(list=ls())

getwd()
#list.files()

library(car)
library(multcomp)
library(reshape2)
library(tidyverse)
library(stats) 
library(tidyr)
library(pracma) 
library(scales)
library(MASS)
library(rcompanion)

##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #pc
lnames = load(file="main_dfs.RData") 
#setwd("~/OneDrive - King's College London/PhD/Projects/SNPs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/SNPs") #pc
lnames



############################################################################################
##                                                                                        ##
##                                      df prep                                           ##
##                                                                                        ##
############################################################################################

data_early[1:5, 1:5]
dim(data_early)
table(data_early$BV_categories, data_early$BV.only)
table(data_early$BV_categories, data_early$BV.Int)


##Making SNPs_early dfs
SNPs_early <- data_early %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040, 
         BV_categories)


#reordering factor
table(SNPs_early$BV_categories)
SNPs_early$BV_categories <- factor(SNPs_early$BV_categories, levels=c("Normal", "Intermediate", "BV"))


##get rid of samples where we have NAs for ALL SNPs 
dim(SNPs_early) #302   
SNPs_early_copy <- SNPs_early #make df copy
SNPs_early <- SNPs_early[!with(SNPs_early, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(SNPs_early) #275       


##get rid of samples where we have NAs for BV 
SNPs_early_copy <- SNPs_early #make df copy
SNPs_early <- SNPs_early[!with(SNPs_early, is.na(BV_categories) ),]
dim(SNPs_early) #260

head(SNPs_early)


#as.character for merge/gathering
SNPs_early$rs17333103 <- as.character(SNPs_early$rs17333103)
SNPs_early$rs1983649 <- as.character(SNPs_early$rs1983649)
SNPs_early$rs2664581 <- as.character(SNPs_early$rs2664581)
SNPs_early$rs6032040 <- as.character(SNPs_early$rs6032040)


#Change genotype formating in dataframe ready for summary df
SNPs_early_gathered <- SNPs_early %>% 
  gather(key=SNP, value=Genotype, 
         -participant_id, -BV_categories) #-BV.only, -BV.Int, -BV.grading

str(SNPs_early_gathered)



# Remove Genotype NAs in SNPs_early_gathered to avoid them getting some of the %/proportion out of the SNP genotypes
table(SNPs_early_gathered$Genotype, useNA="always") #8 NAs
SNPs_early_gathered <- SNPs_early_gathered[!is.na(SNPs_early_gathered$Genotype),]
table(SNPs_early_gathered$Genotype, useNA="always") #0 NAs

#check
str(SNPs_early_gathered)

#Changing chr to factors
SNPs_early_gathered$Genotype <- factor(SNPs_early_gathered$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA"))
SNPs_early_gathered$SNP <- factor(SNPs_early_gathered$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))


#Look
summary(SNPs_early_gathered)
str(SNPs_early_gathered)







############################################################################################
##                                                                                        ##
##                                Plotting SNPs vs BV                                     ##
##                                                                                        ##
############################################################################################

# Getting N numbers for plots
dim(SNPs_early) #260
summary(SNPs_early)# no NAs for BV
table(SNPs_early$rs17333103, useNA="always") #257
table(SNPs_early$rs1983649, useNA="always") #258
table(SNPs_early$rs6032040 , useNA="always") #260
table(SNPs_early$rs2664581, useNA="always") #257



##Get summary df of INSIGHT so we can add n & %s to plots
SNPs_early_BV_categories_summary <- group_by(SNPs_early_gathered, SNP, BV_categories, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(BV_categories, SNP) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
SNPs_early_BV_categories_summary


#plot faceted by SNPs - with n and % labels, don't make a column of NA in BV_categories for each SNP
SNPs_BV_categories_early_plot <- ggplot(SNPs_early_BV_categories_summary, aes(x=BV_categories, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~SNP, ncol = 4) +
  labs(x="BV status", y="Proportion") +
  theme(axis.text.x=element_text(size=5)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=BV_categories, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")), #means label if it is above 3% otherwise do not label!
            colour = "white",  position=position_stack(vjust=0.5), size=1.5)
ggsave("plots/BV_early/BV_SNPs_early.pdf", SNPs_BV_categories_early_plot, height=3, width=6.5)

table(data_early$BV_categories, data_early$rs17333103) #consistent







############################################################################################
##                                                                                        ##
##      Stats tests we can do to test for diff in BV proportions among genotypes          ##
##                                                                                        ##
############################################################################################

names(SNPs_early)
head(SNPs_early)

# 2023 note: Nothing is sig when using the appropriate test - BV early


############################################## rs17333103_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs17333103_BV_categories <- table(SNPs_early$rs17333103, SNPs_early$BV_categories)
forX2_rs17333103_BV_categories
#     Normal Intermediate  BV
#CC    143           26  21
#CT     46            9   5
#TT      3            4   0

#Chi-square test
X2results_rs17333103_BV_categories <- chisq.test(forX2_rs17333103_BV_categories) 
X2results_rs17333103_BV_categories #p-value = 0.03343

#Looking at expected counts 
round(X2results_rs17333103_BV_categories$expected,2) # Expected counts
2/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs17333103_BV_categories) # p-value = 0.09424

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs17333103_BV_categories,
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################







############################################## rs1983649_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs1983649_BV_categories <- table(SNPs_early$rs1983649, SNPs_early$BV_categories)
forX2_rs1983649_BV_categories

#Chi-square test
X2results_rs1983649_BV_categories <- chisq.test(forX2_rs1983649_BV_categories) 
X2results_rs1983649_BV_categories #p-value = 0.5506

#Looking at expected counts 
round(X2results_rs1983649_BV_categories$expected,2) # Expected counts
0/9 #therefore DOES meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(forX2_rs1983649_BV_categories) # p-value = 0.541

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs1983649_BV_categories,
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################




############################################## rs6032040_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs6032040_BV_categories <- table(SNPs_early$rs6032040, SNPs_early$BV_categories)
forX2_rs6032040_BV_categories

#Chi-square test
X2results_rs6032040_BV_categories <- chisq.test(forX2_rs6032040_BV_categories) 
X2results_rs6032040_BV_categories #p-value = 0.7568

#Looking at expected counts 
round(X2results_rs6032040_BV_categories$expected,2) # Expected counts
3/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs6032040_BV_categories) # p-value = 0.936

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs6032040_BV_categories,
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################



############################################## rs2664581_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs2664581_BV_categories <- table(SNPs_early$rs2664581, SNPs_early$BV_categories)
forX2_rs2664581_BV_categories

#Chi-square test
X2results_rs2664581_BV_categories <- chisq.test(forX2_rs2664581_BV_categories) 
X2results_rs2664581_BV_categories #p-value = 0.03046

#Looking at expected counts 
round(X2results_rs2664581_BV_categories$expected,2) # Expected counts
2/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs2664581_BV_categories) # p-value = 0.08093

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs2664581_BV_categories,
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################




