#clear environment
rm(list=ls())


library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(calibrate)
library(devtools)
library(ggbiplot) #install_github("vqv/ggbiplot")


#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc






############################################################################################
##                                                                                        ##
##                           PCA of Mets in data_late                                     ##
##                                                                                        ##
############################################################################################

##Loading .Rdata saved 
lnames = load(file="data_BV/main_dfs.RData") 


## df of just metabolites (Excluding Propylene_Glycol AS IT'S AN INGREDIENT USED IN THE INTERNAL SCAN)
metabolites_late <- data_late[,c(206:208, 210:234)]
head(metabolites_late)
summary(metabolites_late) #no NAs
dim(metabolites_late) #302  28

#Performing PCA
met_late_PCA <- prcomp(metabolites_late, center=T, scale.=T)
summary(met_late_PCA)




############################################################################################
##                                                                                        ##
##                                PCA plots for thesis                                    ##
##                                                                                        ##
############################################################################################
?ggbiplot # doesn't seem to be a way to not plot the NAs


#PCA plot w/ BV - best
PCA_met_late_group_BV <- ggbiplot(met_late_PCA, ellipse=F, var.axes=F,  group=data_late$BV_categories, alpha=0.6) + 
  scale_y_continuous(breaks=c(-6,-5,-4,-3,-2,-1,0,1,2,3,4)) +
  labs(color="BV status") +
  scale_color_manual(values = c("#00CC00","#FF9900","tomato"), na.value = "grey50")
ggsave("plots/PCA_met_late/_PCA_late_BV.pdf", PCA_met_late_group_BV , height=3.5, width=4.5)

#PCA plot w/ ethnicity - best
PCA_met_late_group_ethnicity <- ggbiplot(met_late_PCA, ellipse=F, var.axes=F,  group=data_late$Ethnicity_long, alpha=0.6) + 
  scale_y_continuous(breaks=c(-6,-5,-4,-3,-2,-1,0,1,2,3,4)) +
  labs(color="Ethnicity") +
  scale_color_manual(values = c("purple", "#0099FF", "#FF3399"), na.value="grey50") 
ggsave("plots/PCA_met_late/_PCA_late_ethnicity.pdf", PCA_met_late_group_ethnicity, height=3.5, width=4.5)


