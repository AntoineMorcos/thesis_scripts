#clear environment
rm(list=ls())


library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(calibrate)
library(ggbiplot)

#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc






############################################################################################
##                                                                                        ##
##                          PCA of Mets in data_early                                     ##
##                                                                                        ##
############################################################################################

##Loading .Rdata saved 
lnames = load(file="data_BV/main_dfs.RData") 


## df of just metabolites (Excluding Propylene_Glycol AS IT'S AN INGREDIENT USED IN THE INTERNAL SCAN)
metabolites_early <- data_early[,c(206:208, 210:234)]
head(metabolites_early)
summary(metabolites_early) #no NAs
dim(metabolites_early) #302  28

#Performing PCA
?prcomp
met_early_PCA <- prcomp(metabolites_early, center=T, scale.=T)
summary(met_early_PCA)



############################################################################################
##                                                                                        ##
##                                PCA plots for thesis                                    ##
##                                                                                        ##
############################################################################################
?ggbiplot # doesn't seem like there's a way to not plot the NAs


#PCA plot w/ BV 
PCA_met_early_group_BV <- ggbiplot(met_early_PCA, ellipse=F, var.axes=F,  group=data_early$BV_categories, alpha=0.6) + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3,4)) +
  labs(color="BV status") +
  scale_color_manual(values = c("#00CC00","#FF9900","tomato"), na.value = "grey50")
ggsave("plots/PCA_met_early/_PCA_early_BV.pdf", PCA_met_early_group_BV , height=3.5, width=4.5)


#PCA plot w/ ethnicity
PCA_met_early_group_ethnicity <- ggbiplot(met_early_PCA, ellipse=F, var.axes=F,  group=data_early$Ethnicity_long, alpha=0.6) + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3,4)) +
  labs(color="Ethnicity") +
  scale_color_manual(values = c("purple", "#0099FF", "#FF3399"), na.value="grey50") 
ggsave("plots/PCA_met_early/_PCA_early_ethnicity.pdf", PCA_met_early_group_ethnicity, height=3.5, width=4.5)




