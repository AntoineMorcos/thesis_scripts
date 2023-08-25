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

data_late[1:5, 1:5]
dim(data_late)
table(data_late$BV_categories, data_late$BV.only)
table(data_late$BV_categories, data_late$BV.Int)


##Making SNPs_late dfs
SNPs_late <- data_late %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040, 
         BV_categories)


#reordering factor
table(SNPs_late$BV_categories)
SNPs_late$BV_categories <- factor(SNPs_late$BV_categories, levels=c("Normal", "Intermediate", "BV"))


##get rid of samples where we have NAs for ALL SNPs 
dim(SNPs_late)
SNPs_late_copy <- SNPs_late #make df copy
SNPs_late <- SNPs_late[!with(SNPs_late, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(SNPs_late) #275   9 #dropped 27 samples

##get rid of samples where we have NAs for BV 
SNPs_late_copy <- SNPs_late #make df copy
SNPs_late <- SNPs_late[!with(SNPs_late, is.na(BV_categories) ),]
dim(SNPs_late) #269

head(SNPs_late)


#as.character for merge/gathering 
SNPs_late$rs17333103 <- as.character(SNPs_late$rs17333103)
SNPs_late$rs1983649 <- as.character(SNPs_late$rs1983649)
SNPs_late$rs2664581 <- as.character(SNPs_late$rs2664581)
SNPs_late$rs6032040 <- as.character(SNPs_late$rs6032040)


#Change genotype formating in dataframe ready for summary df
SNPs_late_gathered <- SNPs_late %>% 
  gather(key=SNP, value=Genotype, 
         -participant_id, -BV_categories)

str(SNPs_late_gathered)



# Remove Genotype NAs in SNPs_late_gathered to avoid them getting some of the %/proportion out of the SNP genotypes
table(SNPs_late_gathered$Genotype, useNA="always") #8 NAs
SNPs_late_gathered <- SNPs_late_gathered[!is.na(SNPs_late_gathered$Genotype),]
table(SNPs_late_gathered$Genotype, useNA="always") #0 NAs



#Changing chr to factors
SNPs_late_gathered$Genotype <- factor(SNPs_late_gathered$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA"))
SNPs_late_gathered$SNP <- factor(SNPs_late_gathered$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))

#Look
summary(SNPs_late_gathered)
str(SNPs_late_gathered)





############################################################################################
##                                                                                        ##
##                                Plotting SNPs vs BV                                     ##
##                                                                                        ##
############################################################################################

# Getting N numbers for plots
dim(SNPs_late) #269
summary(SNPs_late)# no NAs for BV
table(SNPs_late$rs17333103, useNA="always") #266
table(SNPs_late$rs1983649, useNA="always") #267
table(SNPs_late$rs6032040 , useNA="always") #269
table(SNPs_late$rs2664581, useNA="always") #266



##Get summary df of INSIGHT so we can add n & %s to plots
SNPs_late_BV_categories_summary <- group_by(SNPs_late_gathered, SNP, BV_categories, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(BV_categories, SNP) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))

str(SNPs_late_BV_categories_summary)


#plot faceted by SNPs - with n and % labels, don't make a column of NA in BV_categories for each SNP
SNPs_BV_categories_late_plot <- ggplot(SNPs_late_BV_categories_summary, aes(x=BV_categories, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~SNP, ncol = 4) +
  labs(x="BV status", y="Proportion") +
  theme(axis.text.x=element_text(size=5)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=BV_categories, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")), #means label if it is above 3% otherwise do not label!
            colour = "white",  position=position_stack(vjust=0.5), size=1.5)
ggsave("plots/BV_late/BV_SNPs_late.pdf", SNPs_BV_categories_late_plot, height=3, width=6.5)

table(data_late$BV_categories, data_late$rs17333103) #consistent








############################################################################################
##                                                                                        ##
##      Stats tests we can do to test for diff in BV proportions among genotypes          ##
##                                                                                        ##
############################################################################################

# 2023 note: Nothing is sig when using the appropriate test - BV late


############################################## rs17333103_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs17333103_BV_categories <- table(SNPs_late$rs17333103, SNPs_late$BV_categories)
forX2_rs17333103_BV_categories

#Chi-square test
X2results_rs17333103_BV_categories <- chisq.test(forX2_rs17333103_BV_categories) 
X2results_rs17333103_BV_categories #p-value = 0.8506

#Looking at expected counts 
round(X2results_rs17333103_BV_categories$expected,2) # Expected counts
3/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs17333103_BV_categories) # p-value = 0.8854

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs17333103_BV_categories, 
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################

############################################## rs1983649_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs1983649_BV_categories <- table(SNPs_late$rs1983649, SNPs_late$BV_categories)
forX2_rs1983649_BV_categories

#Chi-square test
X2results_rs1983649_BV_categories <- chisq.test(forX2_rs1983649_BV_categories) 
X2results_rs1983649_BV_categories #p-value = 0.4832

#Looking at expected counts 
round(X2results_rs1983649_BV_categories$expected, 2) # Expected counts
1/9 #therefore DOES meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(forX2_rs1983649_BV_categories) # p-value = 0.4977

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs1983649_BV_categories, 
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################


############################################## rs6032040_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs6032040_BV_categories <- table(SNPs_late$rs6032040, SNPs_late$BV_categories)
forX2_rs6032040_BV_categories

#Chi-square test
X2results_rs6032040_BV_categories <- chisq.test(forX2_rs6032040_BV_categories) 
X2results_rs6032040_BV_categories #p-value = 0.8212

#Looking at expected counts 
round(X2results_rs6032040_BV_categories$expected,2) # Expected counts
3/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs6032040_BV_categories) # p-value = 0.9639

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs6032040_BV_categories, 
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################



############################################## rs2664581_BV_categories
#Making df for Chi-square test and Fisher's exact test
forX2_rs2664581_BV_categories <- table(SNPs_late$rs2664581, SNPs_late$BV_categories)
forX2_rs2664581_BV_categories

#Chi-square test
X2results_rs2664581_BV_categories <- chisq.test(forX2_rs2664581_BV_categories) 
X2results_rs2664581_BV_categories #p-value = 0.7091

#Looking at expected counts 
round(X2results_rs2664581_BV_categories$expected,2) # Expected counts
2/9 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(forX2_rs2664581_BV_categories) # p-value = 0.6749

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(forX2_rs2664581_BV_categories, 
                            fisher=T, chisq=T, gtest=F, 
                            method="fdr", digits=4)

##############################################



