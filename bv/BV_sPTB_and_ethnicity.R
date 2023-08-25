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
library(tidyr)
library(rstatix)
library(ggpubr)
library(MASS)
library(rcompanion) #https://rcompanion.org/rcompanion/b_07.html


#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc


##Loading .Rdata saved 
lnames = load(file="data_BV/main_dfs.RData") 


############################################################################################
##                                                                                        ##
##                                         Checks                                         ##
##                                                                                        ##
############################################################################################

#early
#n numbers
table(data_early$sPTB_37w, data_early$BV_categories, useNA="always")
table(data_early$sPTB_34w, data_early$BV_categories, useNA="always")
table(data_early$Term_outcome, data_early$BV_categories, useNA="always")

tapply((data1$Gestation.at.delivery..days/7), data1$Term_outcome, summary) #yes, categories are correctly labelled

#late
#n numbers
table(data_late$sPTB_37w, data_late$BV_categories, useNA="always")
table(data_late$sPTB_34w, data_late$BV_categories, useNA="always")
table(data_late$Term_outcome, data_late$BV_categories, useNA="always")


#by ethnicity...
#n numbers
table(data_early$Ethnicity_long, data_early$BV_categories, useNA="always")
table(data_late$Ethnicity_long, data_late$BV_categories, useNA="always")

#BV numbers
table(data_early$BV_categories, useNA="always")
table(data_late$BV_categories, useNA="always")
table(data1$BV_categories, data1$weeks.at.visit.group, useNA="always")







############################################################################################
##                                                                                        ##
##                                      Plotting prep                                     ##
##                                                                                        ##
############################################################################################

#change the order of the factor
data1$BV_categories <- factor(data1$BV_categories, levels=c("BV", "Intermediate", "Normal"))
table(data1$BV_categories, data1$Bacterial.Vaginosis..categories)


#df of only variables we want
data1_sml <- data1 %>% 
  dplyr::select(participant_id, weeks.at.visit.group, BV_categories, sPTB_34w, sPTB_37w, Ethnicity_long, Ethnicity_short, Term_outcome)

# Remove samples with missing BV status
dim(data1_sml) #604   8
data1_sml <- subset(data1_sml, !is.na(BV_categories))
dim(data1_sml) #578   8

#Change dataframe for plot
data1_sml_gath <- data1_sml %>% 
  gather(key=sPTB_cutoff, value=sPTB, 
         -participant_id, -weeks.at.visit.group, -BV_categories, -Ethnicity_long, -Ethnicity_short, -Term_outcome)

colnames(data1_sml_gath)

#Making characters to factors
str(data1_sml_gath)
data1_sml_gath$sPTB <- factor(data1_sml_gath$sPTB, levels=c("No", "Yes"))
data1_sml_gath$sPTB_cutoff <- factor(data1_sml_gath$sPTB_cutoff, levels=c("sPTB_34w", "sPTB_37w"))
str(data1_sml_gath)


#Creating facet labels 
weeks.at.visit.group.labs <- c("Sample from 10-15 weeks' gestation", "Sample from 16-24 weeks' gestation")
names(weeks.at.visit.group.labs) <- c("10-15_weeks", "16-24_weeks")
sPTB_cutoff.labs <- c("sPTB <34 weeks' gestation", "sPTB <37 weeks' gestation")
names(sPTB_cutoff.labs) <- c("sPTB_34w", "sPTB_37w")

#checks
table(data_late$sPTB_37w, data_late$Bacterial.Vaginosis..categories) 




############################################################################################
##                                                                                        ##
##                                BV & ethnicity  - stat tests                            ##
##                                                                                        ##
############################################################################################

############################################## early_BV_Ethnicity
# Remove women labeled as Other ethnicity, as not enough N 
table(data_early$Ethnicity_short, useNA="always")
data_early_whiteBlack <- data_early %>% dplyr::filter(Ethnicity_short=="White" | Ethnicity_short=="Black")
data_early_whiteBlack$Ethnicity_short <- droplevels(data_early_whiteBlack$Ethnicity_short) #remove "Other" level
table(data_early_whiteBlack$Ethnicity_short, useNA="always")
table(data_early_whiteBlack$BV_categories, useNA="always")
dim(data_early_whiteBlack) #272 
272-17 #255

#n numbers
table(data_early_whiteBlack$Ethnicity_short, data_early_whiteBlack$BV_categories)

#Making df for stat test
table_early_BV_Ethnicity <- table(data_early_whiteBlack$Ethnicity_short, data_early_whiteBlack$BV_categories)
table_early_BV_Ethnicity

#Chi-square test
?chisq.test
X2results_early_BV_Ethnicity <- chisq.test(table_early_BV_Ethnicity) 
X2results_early_BV_Ethnicity #p-value = 3.79e-06

#Looking at expected counts 
round(X2results_early_BV_Ethnicity$expected,2) # Expected counts
6/6 #it meets criteria for chi-square as all are >5 count

#Fisher's exact test
fisher.test(table_early_BV_Ethnicity) #p-value = 6.313e-06

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_Ethnicity,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)

##############################################


############################################## late_BV_Ethnicity
# Remove women labeled as Other ethnicity, as not enough N 
table(data_late$Ethnicity_short, useNA="always")
data_late_whiteBlack <- data_late %>% dplyr::filter(Ethnicity_short=="White" | Ethnicity_short=="Black")
data_late_whiteBlack$Ethnicity_short <- droplevels(data_late_whiteBlack$Ethnicity_short) #remove "Other" level
table(data_late_whiteBlack$Ethnicity_short, useNA="always")
table(data_late_whiteBlack$BV_categories, useNA="always")
dim(data_late_whiteBlack) #272 
272-7 #265

#n numbers
table(data_late_whiteBlack$Ethnicity_short, data_late_whiteBlack$BV_categories)

#Making df for stat test
table_late_BV_Ethnicity <- table(data_late_whiteBlack$Ethnicity_short, data_late_whiteBlack$BV_categories)
table_late_BV_Ethnicity

#Chi-square test
X2results_late_BV_Ethnicity <- chisq.test(table_late_BV_Ethnicity) 
X2results_late_BV_Ethnicity #p-value = 0.02403

#Looking at expected counts 
round(X2results_late_BV_Ethnicity$expected,2) # Expected counts
6/6 #therefore does meet criteria for chi-square as it's <20%

#Fisher's exact test
fisher.test(table_late_BV_Ethnicity) #p-value = 0.02446

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_Ethnicity,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)

##############################################





############################################################################################
##                                                                                        ##
##                                  BV & ethnicity  - Plots                               ##
##                                                                                        ##
############################################################################################

# Remove women labeled as Other ethnicity, as not enough N 
table(data1_sml$Ethnicity_short, useNA="always")
data1_sml_whiteBlack <- data1_sml %>% dplyr::filter(Ethnicity_short=="White" | Ethnicity_short=="Black")
table(data1_sml_whiteBlack$Ethnicity_short, useNA="always")

#Now summary df so we can add %s 
BV_ethnicity_summary <- dplyr::group_by(data1_sml_whiteBlack, weeks.at.visit.group,  Ethnicity_short, BV_categories) %>% 
  dplyr::summarise(n = length(participant_id)) %>% 
  ungroup %>% dplyr::group_by(weeks.at.visit.group, Ethnicity_short) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
str(BV_ethnicity_summary)
BV_ethnicity_summary 



#plot of BV in early and late by ethnicity
?annotate
BV_ethnicity_summary_facet <- ggplot(BV_ethnicity_summary, aes(x=Ethnicity_short, fill=BV_categories, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid( ~ weeks.at.visit.group, 
              labeller=labeller(weeks.at.visit.group=weeks.at.visit.group.labs)) + 
  labs(x="Ethnicity", y="Proportion") +
  guides(fill=guide_legend(title="BV status")) +
  geom_text(aes(x=Ethnicity_short, label=paste0(percentage, "%, n=", n)),
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  annotate(geom = 'text', label=("*"), 
           size=6, x = -Inf, y = Inf, hjust = -10.2, vjust = 1.2) + #positive v is down, negative h is down
  theme(strip.text.x = element_text(size = 7.5)) +
  scale_fill_manual(values = c("red", "#FF9933", "#00CC33"), na.value="grey50") #red, orange, green 
ggsave("plots/BV_sPTB_eth/BV_ethnicity.pdf", BV_ethnicity_summary_facet, height=4, width=6)

# Checking
table(data_early$Ethnicity_long, data_early$Bacterial.Vaginosis..categories)







############################################################################################
##                                                                                        ##
##                             Stat tests - BV & sPTB                                     ##
##                                                                                        ##
############################################################################################
# Using Fisher's exact tests as too many tables did not meet the criteria for chi-square tests

#n numbers
table(data_early$BV_categories, useNA="always")
dim(data_early)
302-18
table(data_late$BV_categories, useNA="always")
dim(data_late)
302-8

# Any missing sPTB34/37 outcomes? No
table(data_late$sPTB_34w, useNA="always")
table(data_late$sPTB_37w, useNA="always")




############################################## early_BV_categories_34w
#Making df for stat test
table_early_BV_categories_34w <- table(data_early$BV_categories, data_early$sPTB_34w)
table_early_BV_categories_34w

#Chi-square test
X2results_early_BV_categories_34w <- chisq.test(table_early_BV_categories_34w) 
X2results_early_BV_categories_34w #p-value = 0.3867

#Looking at expected counts 
round(X2results_early_BV_categories_34w$expected,2) # Expected counts
2/6 #therefore does NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(table_early_BV_categories_34w) #p-value = 0.4728

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## late_BV_categories_34w
#Making df for stat test
table_late_BV_categories_34w <- table(data_late$BV_categories, data_late$sPTB_34w)
table_late_BV_categories_34w

#Chi-square test
X2results_late_BV_categories_34w <- chisq.test(table_late_BV_categories_34w) 
X2results_late_BV_categories_34w #p-value = 0.1819

#Looking at expected counts 
round(X2results_late_BV_categories_34w$expected,2) # Expected counts
2/6 #therefore does NOT meets criteria for chi-square as it's >20%

#Fisher's exact test
fisher.test(table_late_BV_categories_34w) #p-value = 0.1333 

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## early_BV_categories_37w
#Making df for stat test
table_early_BV_categories_37w <- table(data_early$BV_categories, data_early$sPTB_37w)
table_early_BV_categories_37w

#Chi-square test
X2results_early_BV_categories_37w <- chisq.test(table_early_BV_categories_37w) 
X2results_early_BV_categories_37w #p-value = 0.6782

#Looking at expected counts 
round(X2results_early_BV_categories_37w$expected,2) # Expected counts
1/6 #therefore does meets criteria for chi-square as it's <20%, but still quite borderline

#Fisher's exact test
fisher.test(table_early_BV_categories_37w) #p-value = 0.608

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################



############################################## late_BV_categories_37w
#Making df for stat test
table_late_BV_categories_37w <- table(data_late$BV_categories, data_late$sPTB_37w)
table_late_BV_categories_37w

#Chi-square test
X2results_late_BV_categories_37w <- chisq.test(table_late_BV_categories_37w) 
X2results_late_BV_categories_37w #p-value = 0.3466

#Looking at expected counts 
round(X2results_late_BV_categories_37w$expected,2) # Expected counts
1/6 #therefore does meets criteria for chi-square as it's <20%, but still quite borderline

#Fisher's exact test
fisher.test(table_late_BV_categories_37w) #p-value = 0.3035

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################





############################################################################################
##                                                                                        ##
##                      BV & sPTB <34 or BV & sPTB <37  - Plots                           ##
##                                                                                        ##
############################################################################################

dim(data1_sml_gath)

#DF of only <34 data
data1_sml_gath_34 <- data1_sml_gath %>% 
  filter(sPTB_cutoff=="sPTB_34w")

#DF of only <37 data
data1_sml_gath_37 <- data1_sml_gath %>% 
  filter(sPTB_cutoff=="sPTB_37w")


################################ 34
#DF of only <34 data
head(data1_sml_gath_34)
dim(data1_sml_gath_34)


#Now summary df so we can add %s of only <34 data 
data1_summary_34 <- dplyr::group_by(data1_sml_gath_34, weeks.at.visit.group, BV_categories, sPTB) %>% 
  dplyr::summarise(n = length(participant_id)) %>% 
  ungroup %>% dplyr::group_by(weeks.at.visit.group, sPTB) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
data1_summary_34


#plot of BV & sPTB <34 
sPTB_34_BV <- ggplot(data1_summary_34, aes(x=sPTB, fill=BV_categories, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid( ~ weeks.at.visit.group, 
              labeller=labeller(weeks.at.visit.group=weeks.at.visit.group.labs)) + 
  labs(x="sPTB <34 weeks' gestation", y="Proportion") +
  guides(fill=guide_legend(title="BV status")) +
  geom_text(aes(x=sPTB, label=paste0(percentage, "%, n=", n)),
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  theme(strip.text.x = element_text(size = 7)) +
  scale_fill_manual(values = c("red", "#FF9933", "#00CC33"), na.value="grey50") #red, orange, green 
ggsave("plots/BV_sPTB_eth/BV_sPTB_34.pdf", sPTB_34_BV, height=3.5, width=5.5)


table(data_early$BV_categories, data_early$sPTB_34w)
table(data_late$BV_categories, data_late$sPTB_34w)
################################




################################ 37
#DF of only <37 data
head(data1_sml_gath_37)

#Now summary df so we can add %s of only <37 data 
data1_summary_37 <- dplyr::group_by(data1_sml_gath_37, weeks.at.visit.group, BV_categories, sPTB) %>% 
  dplyr::summarise(n = length(participant_id)) %>% 
  ungroup %>% dplyr::group_by(weeks.at.visit.group, sPTB) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))

#plot of BV & sPTB <37 
sPTB_37_BV <- ggplot(data1_summary_37, aes(x=sPTB, fill=BV_categories, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid( ~ weeks.at.visit.group, 
              labeller=labeller(weeks.at.visit.group=weeks.at.visit.group.labs)) + 
  labs(x="sPTB <37 weeks' gestation", y="Proportion") +
  guides(fill=guide_legend(title="BV status")) +
  geom_text(aes(x=sPTB, label=paste0(percentage, "%, n=", n)),
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  theme(strip.text.x = element_text(size = 7)) +
  scale_fill_manual(values = c("red", "#FF9933", "#00CC33"), na.value="grey50") #red, orange, green 
ggsave("plots/BV_sPTB_eth/BV_sPTB_37.pdf", sPTB_37_BV, height=3.5, width=5.5)


table(data_early$BV_categories, data_early$sPTB_37w)
table(data_late$BV_categories, data_late$sPTB_37w)
################################








############################################################################################
##                                                                                        ##
##                   Ethnicity_long & BV categories vs sPTB in  - Plots                   ##
##                                                                                        ##
############################################################################################


#Now summary df so we can add %s of only <34 data 
data1_summary_34_eth <- dplyr::group_by(data1_sml_gath_34, weeks.at.visit.group, BV_categories, Ethnicity_long, sPTB) %>% 
  dplyr::summarise(n = length(participant_id)) %>% 
  ungroup %>% dplyr::group_by(weeks.at.visit.group, Ethnicity_long, sPTB) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))

#Now summary df so we can add %s of only <37 data 
data1_summary_37_eth <- dplyr::group_by(data1_sml_gath_37, weeks.at.visit.group, BV_categories, Ethnicity_long, sPTB) %>% 
  dplyr::summarise(n = length(participant_id)) %>% 
  ungroup %>% dplyr::group_by(weeks.at.visit.group, Ethnicity_long, sPTB) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))


#checks
table(data_late$Ethnicity_long, data_late$Bacterial.Vaginosis..categories)




######### Now ONLY WHITE & BLACK, NO OTHER

###Now df with only white & black of only <34 data 
data1_summary_34_eth_white_black <- data1_summary_34_eth %>% 
  filter(Ethnicity_long=="White ethnicity" | Ethnicity_long=="Black ethnicity")
data1_summary_34_eth_white_black
str(data1_summary_34_eth_white_black)


#plot of BV in early and late with Ethnicity_long of only <34 data - with n and % labels
Ethnicity_long_BV_facet_cutoff_time_34w_white_black <- ggplot(data1_summary_34_eth_white_black, aes(x=sPTB, fill=BV_categories, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid(Ethnicity_long ~ weeks.at.visit.group, 
             labeller=labeller(weeks.at.visit.group=weeks.at.visit.group.labs)) + 
  labs(x="sPTB <34 weeks' gestation", y="Proportion") +
  guides(fill=guide_legend(title="BV status")) +
  geom_text(aes(x=sPTB, label=paste0(percentage, "%, n=", n)),
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  theme(strip.text.x = element_text(size = 7.5),
        plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("red", "#FF9933", "#00CC33"), na.value="grey50") #red, orange, green 
ggsave("plots/BV_sPTB_eth/BV_sPTB_34_ethnicity.pdf", Ethnicity_long_BV_facet_cutoff_time_34w_white_black, height=6, width=6)





###Now df with only white & black of only <37 data 
data1_summary_37_eth_white_black <- data1_summary_37_eth %>% 
  filter(Ethnicity_long=="White ethnicity" | Ethnicity_long=="Black ethnicity")


#plot of BV in early and late with Ethnicity_long of only <37 data - with n and % labels
Ethnicity_long_BV_facet_cutoff_time_37w_white_black <- ggplot(data1_summary_37_eth_white_black, aes(x=sPTB, fill=BV_categories, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid(Ethnicity_long ~ weeks.at.visit.group, 
             labeller=labeller(weeks.at.visit.group=weeks.at.visit.group.labs)) + 
  labs(x="sPTB <37 weeks' gestation", y="Proportion") +
  guides(fill=guide_legend(title="BV status")) +
  geom_text(aes(x=sPTB, label=paste0(percentage, "%, n=", n)),
            colour = "white",  position=position_stack(vjust=0.5), size=3) + 
  theme(strip.text.x = element_text(size = 7.5),
        plot.title = element_text(size=12)) +
  scale_fill_manual(values = c("red", "#FF9933", "#00CC33"), na.value="grey50") #red, orange, green 
ggsave("plots/BV_sPTB_eth/BV_sPTB_37_ethnicity.pdf", Ethnicity_long_BV_facet_cutoff_time_37w_white_black, height=6, width=6)
















############################################################################################
##                                                                                        ##
##               Stats test for BV & sPTB when stratified by ethnicity white early        ##
##                                                                                        ##
############################################################################################
# Using Fisher's exact tests, as too many tables did not meet the criteria for chi-square tests

#filter only white ethnicity
data_early_white <- data_early %>% 
  filter(Ethnicity_short=="White")

#checks
table(data_early_white$Ethnicity_short)
table(data_early$Ethnicity_short)



############################################## early_BV_categories_34w
#Making df for stat test
table_early_BV_categories_34w <- table(data_early_white$BV_categories, data_early_white$sPTB_34w)
table_early_BV_categories_34w #yes

#Fisher's exact test
fisher.test(table_early_BV_categories_34w) #p-value = 0.2841

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## early_BV_categories_37w
#Making df for stat test
table_early_BV_categories_37w <- table(data_early_white$BV_categories, data_early_white$sPTB_37w)
table_early_BV_categories_37w

#Fisher's exact test
fisher.test(table_early_BV_categories_37w) #p-value = 0.2022

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################






############################################################################################
##                                                                                        ##
##               Stats test for BV & sPTB when stratified by ethnicity white late         ##
##                                                                                        ##
############################################################################################
#Using Fisher's exact tests, as too many tables did not meet the criteria for chi-square tests

#filter only white ethnicity
data_late_white <- data_late %>% 
  filter(Ethnicity_short=="White")

#checks
table(data_late_white$Ethnicity_short)
table(data_late$Ethnicity_short)



############################################## late_BV_categories_34w
#Making df for stat test
table_late_BV_categories_34w <- table(data_late_white$BV_categories, data_late_white$sPTB_34w)
table_late_BV_categories_34w #yes

#Fisher's exact test
fisher.test(table_late_BV_categories_34w) #p-value = 0.07466

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## late_BV_categories_37w
#Making df for stat test
table_late_BV_categories_37w <- table(data_late_white$BV_categories, data_late_white$sPTB_37w)
table_late_BV_categories_37w

#Fisher's exact test
fisher.test(table_late_BV_categories_37w) #p-value = 0.2339

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################










############################################################################################
##                                                                                        ##
##               Stats test for BV & sPTB when stratified by ethnicity black early        ##
##                                                                                        ##
############################################################################################
#Using Fisher's exact tests, as too many tables did not meet the criteria for chi-square tests

#filter only black ethnicity
data_early_black <- data_early %>% 
  filter(Ethnicity_short=="Black")

#checks
table(data_early_black$Ethnicity_short)
table(data_early$Ethnicity_short)



############################################## early_BV_categories_34w
#Making df for stat test
table_early_BV_categories_34w <- table(data_early_black$BV_categories, data_early_black$sPTB_34w)
table_early_BV_categories_34w #yes

#Fisher's exact test
fisher.test(table_early_BV_categories_34w) #p-value = 0.3857

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## early_BV_categories_37w
#Making df for stat test
table_early_BV_categories_37w <- table(data_early_black$BV_categories, data_early_black$sPTB_37w)
table_early_BV_categories_37w

#Fisher's exact test
fisher.test(table_early_BV_categories_37w) #p-value = 0.5863

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_early_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################






############################################################################################
##                                                                                        ##
##               Stats test for BV & sPTB when stratified by ethnicity black late         ##
##                                                                                        ##
############################################################################################
#Using Fisher's exact tests, as too many tables did not meet the criteria for chi-square tests

#filter only black ethnicity
data_late_black <- data_late %>% 
  filter(Ethnicity_short=="Black")

#checks
table(data_late_black$Ethnicity_short)
table(data_late$Ethnicity_short)



############################################## late_BV_categories_34w
#Making df for stat test
table_late_BV_categories_34w <- table(data_late_black$BV_categories, data_late_black$sPTB_34w)
table_late_BV_categories_34w #yes

#Fisher's exact test
fisher.test(table_late_BV_categories_34w) #p-value = 0.3003

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_34w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################


############################################## late_BV_categories_37w
#Making df for stat test
table_late_BV_categories_37w <- table(data_late_black$BV_categories, data_late_black$sPTB_37w)
table_late_BV_categories_37w

#Fisher's exact test
fisher.test(table_late_BV_categories_37w) #p-value = 0.2717

#Fisher & X2 post-hoc pairwise testing
pairwiseNominalIndependence(table_late_BV_categories_37w,
                            fisher=T, chisq=T, gtest=F,
                            digits=4)
##############################################





