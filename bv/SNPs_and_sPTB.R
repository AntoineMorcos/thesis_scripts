#clear environment
rm(list=ls())

library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(dplyr)
library(MASS)
library(tidyr)
library(beanplot)
library(pracma) 
library(rcompanion) #https://rcompanion.org/rcompanion/b_07.html


##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #pc
lnames = load(file="main_dfs.RData") 
#setwd("~/OneDrive - King's College London/PhD/Projects/SNPs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/SNPs") #pc
lnames




############################################################################################
##                                                                                        ##
##                             SNPs and sPTB df prep                                      ##
##                                                                                        ##
############################################################################################

data_early[1:5, 1:5]
dim(data_early)


##Making SNPs - here it is arbitrary whether we use data_early or data_late as we are not using data that changes by timepoint
SNPs <- data_early %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040,
         sPTB_34w, sPTB_37w, Ethnicity_long, Ethnicity_short)

#check
str(SNPs)


##get rid of samples where we have NAs for ALL SNPs 
dim(SNPs) #302
SNPs_copy <- SNPs #make df copy
SNPs <- SNPs[!with(SNPs, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(SNPs) #275


#Compare dfs to check NA drop worked well 
subset(SNPs_copy, !(participant_id %in% SNPs$participant_id)) 
subset(SNPs, !(participant_id %in% SNPs_copy$participant_id)) 

#delete df copy
SNPs_copy <- NULL




#as.character so merge works
SNPs$rs17333103 <- as.character(SNPs$rs17333103)
SNPs$rs1983649 <- as.character(SNPs$rs1983649)
SNPs$rs2664581 <- as.character(SNPs$rs2664581)
SNPs$rs6032040 <- as.character(SNPs$rs6032040)

# Gather df for plotting
?gather
?pivot_longer
SNPs2 <- pivot_longer(SNPs,
                      cols=rs17333103:rs6032040, 
                      names_to="SNP",
                      values_to="Genotype")
head(SNPs2)


# Remove Genotype NAs in SNPs2 to avoid them getting some of the %/proportion out of the SNP genotypes
table(SNPs2$Genotype, useNA="always") #8 NAs
SNPs2 <- SNPs2[!is.na(SNPs2$Genotype),]
table(SNPs2$Genotype, useNA="always") #0 NAs


#check
str(SNPs2)

#Changing chr to factors
SNPs2$Genotype <- factor(SNPs2$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA"))
SNPs2$SNP <- factor(SNPs2$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))
summary(SNPs2)





############################################################################################
##                                                                                        ##
##                            sPTB vs elafin SNP genotypes                                ##
##                                                                                        ##
############################################################################################

# Getting N numbers for sPTB plots
dim(SNPs) #275
summary(SNPs)# no NAs for sPTB
table(SNPs$rs17333103, useNA="always") #272
table(SNPs$rs1983649, useNA="always") #273
table(SNPs$rs6032040 , useNA="always") #275
table(SNPs$rs2664581, useNA="always") #272



################################################################### 34
#Get summary df of INSIGHT so we can add %s 
sPTB_34w_SNPs <- group_by(SNPs2, SNP, sPTB_34w, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(sPTB_34w, SNP) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
sPTB_34w_SNPs

#plot faceted by SNPs 34 - with n and % labels
sPTB_34w_SNPs_plot_labelled <- ggplot(sPTB_34w_SNPs, aes(x=sPTB_34w, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~SNP, ncol = 4) +
  labs(x="sPTB <34 weeks' gestation", y="Proportion") +
  theme(axis.text.x=element_text(size=7)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=sPTB_34w, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")), #means label if it is above 3% otherwise do not label!
            colour = "white",  position=position_stack(vjust=0.5), size=1.9)
ggsave("plots/sPTB_binary/sPTB_34w_SNPs.pdf", sPTB_34w_SNPs_plot_labelled, height=3, width=6)

############################################## rs17333103_34w
#Making df for Chi-square test
rs17333103_34w <- table(SNPs$rs17333103, SNPs$sPTB_34w)
rs17333103_34w

#Chi-square test
X2results_rs17333103_34w <- chisq.test(rs17333103_34w) 
X2results_rs17333103_34w #p-value = 0.8042

#Looking at expected counts
round(X2results_rs17333103_34w$expected,2) # Expected counts
2/6 #therefore does NOT meet criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs17333103_34w) # p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_34w
#Making df for Chi-square test
rs1983649_34w <- table(SNPs$rs1983649, SNPs$sPTB_34w)
rs1983649_34w

#Chi-square test
X2results_rs1983649_34w <- chisq.test(rs1983649_34w) 
X2results_rs1983649_34w #p-value = 0.4844

#Looking at expected counts
round(X2results_rs1983649_34w$expected,2) # Expected counts
1/6 #therefore [marginally] meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs1983649_34w) #p-value = 0.5118

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_34w
#Making df for Chi-square test
rs2664581_34w <- table(SNPs$rs2664581, SNPs$sPTB_34w)
rs2664581_34w

#Chi-square test
X2results_rs2664581_34w <- chisq.test(rs2664581_34w) 
X2results_rs2664581_34w #p-value = 0.6848

#Looking at expected counts
round(X2results_rs2664581_34w$expected,2) # Expected counts
2/6 #therefore meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs2664581_34w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs6032040_34w
#Making df for Chi-square test
rs6032040_34w <- table(SNPs$rs6032040, SNPs$sPTB_34w)
rs6032040_34w

#Chi-square test
X2results_rs6032040_34w <- chisq.test(rs6032040_34w) 
X2results_rs6032040_34w #p-value = 0.8846

#Looking at expected counts
round(X2results_rs6032040_34w$expected,2) # Expected counts
3/6 #therefore DOES NOT meet criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs6032040_34w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################
################################################################### 








################################################################### 37
#Get summary df of INSIGHT so we can add %s 
sPTB_37w_SNPs <- group_by(SNPs2, SNP, sPTB_37w, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(sPTB_37w, SNP) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
sPTB_37w_SNPs

#plot faceted by SNPs 37 - with n and % labels
sPTB_37w_SNPs_plot_labelled <- ggplot(sPTB_37w_SNPs, aes(x=sPTB_37w, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~SNP, ncol = 4) +
  labs(x="sPTB <37 weeks' gestation", y="Proportion") +
  theme(axis.text.x=element_text(size=7)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=sPTB_37w, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")),
            colour = "white",  position=position_stack(vjust=0.5), size=1.9)
ggsave("plots/sPTB_binary/sPTB_37w_SNPs.pdf", sPTB_37w_SNPs_plot_labelled, height=3, width=6)

############################################## rs17333103_37w
#Making df for Chi-square test
rs17333103_37w <- table(SNPs$rs17333103, SNPs$sPTB_37w)
rs17333103_37w

#Chi-square test
X2results_rs17333103_37w <- chisq.test(rs17333103_37w) 
X2results_rs17333103_37w #p-value = 0.432

#Looking at expected counts
round(X2results_rs17333103_37w$expected,2) # Expected counts
1/6 #therefore meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs17333103_37w) #p-value = 0.3394

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_37w
#Making df for Chi-square test
rs1983649_37w <- table(SNPs$rs1983649, SNPs$sPTB_37w)
rs1983649_37w

#Chi-square test
X2results_rs1983649_37w <- chisq.test(rs1983649_37w) 
X2results_rs1983649_37w #p-value = 0.1044

#Looking at expected counts
round(X2results_rs1983649_37w$expected,2) # Expected counts
0/6 #therefore meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs1983649_37w) #p-value = 0.2612

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_37w
#Making df for Chi-square test
rs2664581_37w <- table(SNPs$rs2664581, SNPs$sPTB_37w)
rs2664581_37w

#Chi-square test
X2results_rs2664581_37w <- chisq.test(rs2664581_37w) 
X2results_rs2664581_37w #p-value = 0.5041

#Looking at expected counts
round(X2results_rs2664581_37w$expected,2) # Expected counts
1/6 #therefore meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs2664581_37w) #p-value = 0.694

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs6032040_37w
#Making df for Chi-square test
rs6032040_37w <- table(SNPs$rs6032040, SNPs$sPTB_37w)
rs6032040_37w

#Chi-square test
X2results_rs6032040_37w <- chisq.test(rs6032040_37w) 
X2results_rs6032040_37w #p-value = 0.7202

#Looking at expected counts
round(X2results_rs6032040_37w$expected,2) # Expected counts
2/6 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs6032040_37w) #p-value = 0.8879

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################
###################################################################








############################################################################################
##                                                                                        ##
##             Genotypes by ethnicity and binary sPTB - WHITE BLACK ONLY                  ##
##                                                                                        ##
############################################################################################

# Getting N numbers for sPTB plots
dim(SNPs) #275
SNPs_whiteBlack <- SNPs %>% filter(Ethnicity_short!="Other")
table(SNPs_whiteBlack$Ethnicity_short)
dim(SNPs_whiteBlack) #247
summary(SNPs_whiteBlack)# no NAs for sPTB
table(SNPs_whiteBlack$rs17333103, useNA="always") #244
table(SNPs_whiteBlack$rs1983649, useNA="always") #245
table(SNPs_whiteBlack$rs6032040 , useNA="always") #247
table(SNPs_whiteBlack$rs2664581, useNA="always") #244



######### 34 #Get summary df of INSIGHT so we can add %s 
sPTB_34w_SNPs_Ethn <- group_by(SNPs2, SNP, sPTB_34w, Genotype, Ethnicity_long) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(sPTB_34w, SNP, Ethnicity_long) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
sPTB_34w_SNPs_Ethn


######### 37 #Get summary df of INSIGHT so we can add %s 
sPTB_37w_SNPs_Ethn <- group_by(SNPs2, SNP, sPTB_37w, Genotype, Ethnicity_long) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(sPTB_37w, SNP, Ethnicity_long) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
sPTB_37w_SNPs_Ethn


######### 34
sPTB_34w_SNPs_Ethn_white_black <- sPTB_34w_SNPs_Ethn %>% 
  filter(Ethnicity_long!="Other ethnicity")

#plot faceted by SNPs 34 and ethnicity - with n and % labels white_black
sPTB_34w_SNPs_Ethn_white_black_plot_labelled <- ggplot(sPTB_34w_SNPs_Ethn_white_black, aes(x=sPTB_34w, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid(Ethnicity_long ~ SNP) +
  labs(x="sPTB <34 weeks' gestation", y="Proportion") +
  theme(axis.text.x=element_text(size=7)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=sPTB_34w, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")),
            colour = "white",  position=position_stack(vjust=0.5), size=1.55)
ggsave("plots/sPTB_binary/sPTB_34w_SNPs_Ethn.pdf", sPTB_34w_SNPs_Ethn_white_black_plot_labelled, height=4, width=5)
#########

######### 37
sPTB_37w_SNPs_Ethn_white_black <- sPTB_37w_SNPs_Ethn %>% 
  filter(Ethnicity_long!="Other ethnicity")

#plot faceted by SNPs 37 and ethnicity - with n and % labels white_black
sPTB_37w_SNPs_Ethn_white_black_plot_labelled <- ggplot(sPTB_37w_SNPs_Ethn_white_black, aes(x=sPTB_37w, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid(Ethnicity_long ~ SNP) +
  labs(x="sPTB <37 weeks' gestation", y="Proportion") +
  theme(axis.text.x=element_text(size=7)) +
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=sPTB_37w, label=ifelse(percentage>3, paste0(percentage, "%, n=", n),"")), #means label if it is above 3% otherwise do not label
            colour = "white",  position=position_stack(vjust=0.5), size=1.55)
ggsave("plots/sPTB_binary/sPTB_37w_SNPs_Ethn.pdf", sPTB_37w_SNPs_Ethn_white_black_plot_labelled, height=4, width=6)
#########












############################################################################################
##                                                                                        ##
##                   Stats tests of SNPs and sPTB in White ethnicity only                 ##
##                                                                                        ##
############################################################################################

#all ethnicity df
table(SNPs$Ethnicity_long, useNA="always")

#only white
SNPs_white <- SNPs %>% 
  filter(Ethnicity_long=="White ethnicity")

#check
table(SNPs_white$Ethnicity_long, useNA="always")



############################################## rs17333103_34w
#Making df for Chi-square test
rs17333103_34w <- table(SNPs_white$rs17333103, SNPs_white$sPTB_34w)
rs17333103_34w #same n numbers as plot

#Chi-square test
X2results_rs17333103_34w <- chisq.test(rs17333103_34w) 
X2results_rs17333103_34w #p-value = 0.4017

#Looking at expected counts
round(X2results_rs17333103_34w$expected,2) # Expected counts
3/6 #therefore does  meet criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs17333103_34w) # p-value = 0.3131

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_34w
#Making df for Chi-square test
rs1983649_34w <- table(SNPs_white$rs1983649, SNPs_white$sPTB_34w)
rs1983649_34w #same n numbers as plot

#Chi-square test
X2results_rs1983649_34w <- chisq.test(rs1983649_34w) 
X2results_rs1983649_34w #p-value = 0.07541

#Looking at expected counts
round(X2results_rs1983649_34w$expected,2) # Expected counts
3/6 #therefore does NOT meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs1983649_34w) #p-value = 0.08529 

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
#  Comparison p.Fisher p.adj.Fisher p.Chisq p.adj.Chisq
#1    AA : AT  0.62700       0.6270  0.7515      0.7515
#2    AA : TT  0.06061       0.1818  0.1063      0.3189
#3    AT : TT  0.14360       0.2154  0.2871      0.4306
##############################################

############################################## rs6032040_34w
#Making df for Chi-square test
rs6032040_34w <- table(SNPs_white$rs6032040, SNPs_white$sPTB_34w)
rs6032040_34w #same n numbers as plot

#Chi-square test
X2results_rs6032040_34w <- chisq.test(rs6032040_34w) 
X2results_rs6032040_34w #p-value =  0.7961

#Looking at expected counts
round(X2results_rs6032040_34w$expected,2) # Expected counts
4/6 #therefore DOES NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs6032040_34w) #p-value = 0.7061

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_34w
#Making df for Chi-square test
rs2664581_34w <- table(SNPs_white$rs2664581, SNPs_white$sPTB_34w)
rs2664581_34w #same n numbers as plot

#Chi-square test
X2results_rs2664581_34w <- chisq.test(rs2664581_34w) 
X2results_rs2664581_34w #p-value = 0.5081

#Looking at expected counts
round(X2results_rs2664581_34w$expected,2) # Expected counts
3/6 #therefore DOES NOT meet criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs2664581_34w) #p-value = 0.4631

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################






############################################## rs17333103_37w
#Making df for Chi-square test
rs17333103_37w <- table(SNPs_white$rs17333103, SNPs_white$sPTB_37w)
rs17333103_37w #same n numbers as plot

#Chi-square test
X2results_rs17333103_37w <- chisq.test(rs17333103_37w) 
X2results_rs17333103_37w #p-value = 0.6958

#Looking at expected counts
round(X2results_rs17333103_37w$expected,2) # Expected counts
3/6 #therefore DOES NOT meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs17333103_37w) #p-value = 0.8685

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_37w
#Making df for Chi-square test
rs1983649_37w <- table(SNPs_white$rs1983649, SNPs_white$sPTB_37w)
rs1983649_37w #same n numbers as plot

#Chi-square test
X2results_rs1983649_37w <- chisq.test(rs1983649_37w) 
X2results_rs1983649_37w #p-value =  0.4354

#Looking at expected counts
round(X2results_rs1983649_37w$expected,2) # Expected counts
1/6 #therefore meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs1983649_37w) #p-value = 0.4622

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs6032040_37w
#Making df for Chi-square test
rs6032040_37w <- table(SNPs_white$rs6032040, SNPs_white$sPTB_37w)
rs6032040_37w #same n numbers as plot

#Chi-square test
X2results_rs6032040_37w <- chisq.test(rs6032040_37w) 
X2results_rs6032040_37w #p-value = 0.731

#Looking at expected counts
round(X2results_rs6032040_37w$expected,2) # Expected counts
2/6 #therefore DOES NOT meet criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(rs6032040_37w) #p-value = 0.7015

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_37w
#Making df for Chi-square test
rs2664581_37w <- table(SNPs_white$rs2664581, SNPs_white$sPTB_37w)
rs2664581_37w 

#Chi-square test
X2results_rs2664581_37w <- chisq.test(rs2664581_37w) 
X2results_rs2664581_37w #p-value = 0.7222

#Looking at expected counts
round(X2results_rs2664581_37w$expected,2) # Expected counts
2/6 #therefore DOES NOT meets criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(rs2664581_37w) #p-value = 0.8712

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################










############################################################################################
##                                                                                        ##
##                   Stats tests of SNPs and sPTB in Black ethnicity only                 ##
##                                                                                        ##
############################################################################################

#all ethnicity df
table(SNPs$Ethnicity_long, useNA="always")

#only Black ethnicity
SNPs_black <- SNPs %>% 
  filter(Ethnicity_long=="Black ethnicity")

#check
table(SNPs_black$Ethnicity_long, useNA="always")



############################################## rs17333103_34w
#Making df for Chi-square test
rs17333103_34w <- table(SNPs_black$rs17333103, SNPs_black$sPTB_34w)
rs17333103_34w #same n numbers as plot

#Fisher's exact test
fisher.test(rs17333103_34w) # p-value = 0.5939

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_34w
#Making df for Chi-square test
rs1983649_34w <- table(SNPs_black$rs1983649, SNPs_black$sPTB_34w)
rs1983649_34w #same n numbers as plot

#Fisher's exact test
fisher.test(rs1983649_34w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs6032040_34w
#Making df for Chi-square test
rs6032040_34w <- table(SNPs_black$rs6032040, SNPs_black$sPTB_34w)
rs6032040_34w #same n numbers as plot

#Fisher's exact test
fisher.test(rs6032040_34w) #p-value = 0.7005

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_34w
#Making df for Chi-square test
rs2664581_34w <- table(SNPs_black$rs2664581, SNPs_black$sPTB_34w)
rs2664581_34w #same n numbers as plot

#Fisher's exact test
fisher.test(rs2664581_34w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################






############################################## rs17333103_37w
#Making df for Chi-square test
rs17333103_37w <- table(SNPs_black$rs17333103, SNPs_black$sPTB_37w)
rs17333103_37w #same n numbers as plot

#Fisher's exact test
fisher.test(rs17333103_37w) #p-value = 0.7122

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs17333103_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs1983649_37w
#Making df for Chi-square test
rs1983649_37w <- table(SNPs_black$rs1983649, SNPs_black$sPTB_37w)
rs1983649_37w #same n numbers as plot

#Fisher's exact test
fisher.test(rs1983649_37w) #p-value = 0.893

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs1983649_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs6032040_37w
#Making df for Chi-square test
rs6032040_37w <- table(SNPs_black$rs6032040, SNPs_black$sPTB_37w)
rs6032040_37w #same n numbers as plot

#Fisher's exact test
fisher.test(rs6032040_37w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs6032040_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################

############################################## rs2664581_37w
#Making df for Chi-square test
rs2664581_37w <- table(SNPs_black$rs2664581, SNPs_black$sPTB_37w)
rs2664581_37w 

#Fisher's exact test
fisher.test(rs2664581_37w) #p-value = 1

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(rs2664581_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################






