#clear environment
rm(list=ls())


library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(funrar) #install.packages("funrar")


##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
lnames = load(file="data_BV/main_dfs.RData") 
lnames




############################################################################################
##                                                                                        ##
##                 Early samples (10-15 weeks) making abundance matrix                    ##
##                                                                                        ##
############################################################################################



########################################Making matrix of raw OTU numbers from data_early
#find limits of OTUs in main df
data_early[1:5, 3175:3176] #ends at 3175
data_early[1:5, 235:237] #starts at 235

##Making a OTU df
raw_OTUs_early <- data_early[,c(235:3175)]
raw_OTUs_early[1:5,1:5]

###Deleting rows where there is no OTU data 
dim(raw_OTUs_early) #302
raw_OTUs_early <- raw_OTUs_early %>% 
  drop_na("Otu00001_Lactobacillus_crispatus") 
dim(raw_OTUs_early) #N=283 for women with OTU data


#checks
dim(raw_OTUs_early)
names(raw_OTUs_early[1:6])
names(raw_OTUs_early[2936:2941])
raw_OTUs_early[1:5,1:5]

#checking data in numerical
str(raw_OTUs_early[1:6])

#making a matrix
dim(raw_OTUs_early)
raw_OTUs_early[1:5,1:4]
raw_OTUs_early_matrix <- as.matrix(raw_OTUs_early)
dim(raw_OTUs_early_matrix)
raw_OTUs_early_matrix[1:5,1:4]
########################################



##Looking at OTU read data - number of reads per sample across the raw_OTUs 
#Early read numbers
OTUs_row_sum_early <- rowSums(raw_OTUs_early[,])
summary(OTUs_row_sum_early)
boxplot(OTUs_row_sum_early)





######################################## Making relative abundance matrix
## Making a relative abundance matrix of the OTUs https://rdrr.io/cran/funrar/man/make_relative.html#heading-3
?funrar::make_relative

#preview raw matrix
raw_OTUs_early_matrix[1:10,1:4]

#make relative abundance matrix command
AM_OTUs_early <- make_relative(raw_OTUs_early_matrix)

#preview abundance matrix
AM_OTUs_early[1:10,1:4]
########################################


######################################## Saving as RData
###Saving the data from this script
save(AM_OTUs_early,
     file="data_BV/AM_OTUs_early.RData")
########################################




############################################################################################
##                                                                                        ##
##                    Early samples (10-15 weeks) prep for plots                          ##
##                                                                                        ##
############################################################################################

#save df under new name to not change original df
df_early <- as.data.frame(AM_OTUs_early)

##Get rid of the Otu00000_ part of the OTU columns to make aggregating easier
names(df_early) <- substring(names(df_early), 10) #keep characters from the 10th character onwards

##Getting rid of species part of the name
names(df_early) <- gsub("\\_.*","", names(df_early))

str(df_early)

##merging columns with the same names
df_early_merged <- t(rowsum(t(df_early), group = colnames(df_early), na.rm = T))

##Checks and conversion to a df
dim(df_early_merged)
str(df_early_merged)
df_early_merged <- as.data.frame(df_early_merged)

##Add sample ID column
df_early_merged$Sample_ID <- rownames(df_early_merged) 
df_early_merged[1:10, 73:75]

############ Deleting genus columns where there is not at least 1 person with an abundance of >10% for that genus & add an others column
dim(df_early_merged) #283  75

#df_early_merged_filtered <- df_early_merged %>% 
#  select_if(~max(., na.rm = TRUE) >= 0.05) %>%     #####At 5%
#  mutate(Other = 1-(rowSums(.[1:20], na.rm=TRUE))) #20 for number of rows (21) minus the sample row

df_early_merged_filtered <- df_early_merged %>% 
  select_if(~max(., na.rm = TRUE) >= 0.1) %>%     #####At 10%
  mutate(Other = 1-(rowSums(.[1:13], na.rm=TRUE))) #########################1:13 is all rows after the filtering from the previous line of 10% minus the sample ID row

dim(df_early_merged_filtered) #283  15
str(df_early_merged_filtered)

##Making the other col numeric
df_early_merged_filtered$Other <- as.numeric(df_early_merged_filtered$Other)
  

##Making df of just clinical cols of interest
clinical_early <- data_early %>% 
  select(Sample_ID, participant_id, weeks.at.visit.group, risk.cat, Ethnicity_long, Ethnicity_short, group.pcoa.Erica,
         Age, BMI, BMI_category, Smoking, BV_categories, sPTB_34w, sPTB_37w)


#Changing sPTB_34 so it is understandable when using facet_grid
table(clinical_early$sPTB_34w) 
clinical_early$sPTB_34w_copy <- clinical_early$sPTB_34w
clinical_early$sPTB_34w <- as.character(clinical_early$sPTB_34w)
clinical_early$sPTB_34w[clinical_early$sPTB_34w == "No"] <- "Not sPTB <34"
clinical_early$sPTB_34w[clinical_early$sPTB_34w == "Yes"] <- "sPTB <34"
clinical_early$sPTB_34w <- factor(clinical_early$sPTB_34w, levels=c("Not sPTB <34", "sPTB <34"))
table(clinical_early$sPTB_34w, clinical_early$sPTB_34w_copy)
clinical_early$sPTB_34w_copy <- NULL

#Changing sPTB_37w so it is understandable when using facet_grid
table(clinical_early$sPTB_37w) 
clinical_early$sPTB_37w_copy <- clinical_early$sPTB_37w
clinical_early$sPTB_37w <- as.character(clinical_early$sPTB_37w)
clinical_early$sPTB_37w[clinical_early$sPTB_37w == "No"] <- "Not sPTB <37"
clinical_early$sPTB_37w[clinical_early$sPTB_37w == "Yes"] <- "sPTB <37"
clinical_early$sPTB_37w <- factor(clinical_early$sPTB_37w, levels=c("Not sPTB <37", "sPTB <37"))
table(clinical_early$sPTB_37w, clinical_early$sPTB_37w_copy)
clinical_early$sPTB_37w_copy <- NULL

#Changing risk.cat so it is understandable when using facet_grid
table(clinical_early$risk.cat) 
clinical_early$risk.cat_copy <- clinical_early$risk.cat
clinical_early$risk.cat <- as.character(clinical_early$risk.cat)
clinical_early$risk.cat[clinical_early$risk.cat == "Low"] <- "Low risk"
clinical_early$risk.cat[clinical_early$risk.cat == "High"] <- "High risk"
clinical_early$risk.cat <- factor(clinical_early$risk.cat, levels=c("Low risk", "High risk"))
table(clinical_early$risk.cat, clinical_early$risk.cat_copy)
clinical_early$risk.cat_copy <- NULL


#merging dfs
early_spread <- merge(clinical_early, df_early_merged_filtered, by="Sample_ID")


#Making a new column to use for plotting
early_spread_sort <- early_spread %>% 
  mutate(Label_gg=str_c(early_spread$Ethnicity_long, ", ", early_spread$risk.cat, ", Age ", early_spread$Age, ", ID:", early_spread$Sample_ID)) %>% 
  arrange(Ethnicity_long, risk.cat, Age)

# Make new col for ordering the plots
early_spread_sort$Lacto_copy <- early_spread_sort$Lactobacillus

# Look
early_spread_sort 
str(early_spread_sort)


##looking at n numbers  for plot
dim(early_spread_sort) #N=283
table(early_spread_sort$Ethnicity_short, useNA="always")
table(early_spread_sort$BV_categories, useNA="always")
283-17
table(early_spread_sort$risk.cat, useNA="always")
table(early_spread_sort$group.pcoa.Erica, useNA="always")
table(early_spread_sort$sPTB_34w, useNA="always")
table(early_spread_sort$sPTB_37w, useNA="always")


#gathering genera into 1 column
?gather
?pivot_longer
early_gather <- pivot_longer(early_spread_sort, 
                             cols=Aerococcus:Other, 
                             names_to="Genus", values_to="Abundance")

dim(early_gather)
str(early_gather)

#putting Other category to the end of the factors 
early_gather$Genus<- fct_relevel(early_gather$Genus, "Other", after = Inf)








############################################################################################
##                                                                                        ##
##                     Plotting Early samples (10-15 weeks) FACETING                      ##
##                                                                                        ##
############################################################################################

# Ordering plots by Lactobacillus abundance
?reorder 

######################### BV
##plot faceted by BV - blank x axis
Genus_early_plot_facet_BV_blank <- ggplot(subset(early_gather, !is.na(BV_categories)), aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~BV_categories, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm")) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_BV.pdf", Genus_early_plot_facet_BV_blank, height=2.3, width=14)
######################### 


######################### Ethnicity_long
##plot faceted by ethnicity - blank x axis
Genus_early_plot_facet_ethnicity_blank <- ggplot(early_gather, aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~Ethnicity_long, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm")) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_ethnicity.pdf", Genus_early_plot_facet_ethnicity_blank, height=2.3, width=14)
######################### 



######################### Risk
##plot faceted by risk - blank x axis
Genus_early_plot_facet_risk_blank <- ggplot(early_gather, aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~risk.cat, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm")) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_risk.pdf", Genus_early_plot_facet_risk_blank, height=2.3, width=14)
######################### 



######################### sPTB_34w 
##plot faceted by sPTB_34w - blank x axis
Genus_early_plot_facet_sPTB_34w_blank <- ggplot(early_gather, aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~sPTB_34w, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm")) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_sPTB_34w.pdf", Genus_early_plot_facet_sPTB_34w_blank, height=2.3, width=14)
######################### 


######################### sPTB_37w 
##plot faceted by sPTB_37w - blank x axis
Genus_early_plot_facet_sPTB_37w_blank <- ggplot(early_gather, aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~sPTB_37w, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm")) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_sPTB_37w.pdf", Genus_early_plot_facet_sPTB_37w_blank, height=2.3, width=14)
######################### 







############################################################################################
##                                                                                        ##
##        Plotting early samples   FACET GRID        white black only                     ##
##                                                                                        ##
############################################################################################

# White black only
early_gather_white_black <- early_gather %>% 
  filter(Ethnicity_long != "Other ethnicity") 
early_gather_white_black

# N numbers
dim(early_spread_sort) #N=283
table(early_spread_sort$BV_categories, early_spread_sort$Ethnicity_short, useNA="always") #
table(early_spread_sort$BV_categories, useNA="always") 
145+22+13+27+20+12 #239
283-25-2-1-13-3 #239




######################### BV + Ethnicity
##plot faceted by BV - blank x axis
Genus_early_plot_facet_BV_eth_blank2 <- ggplot(subset(early_gather_white_black, !is.na(BV_categories)), aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~Ethnicity_long+BV_categories, scales = "free", nrow=2) +  #  with ,  but then it looks ugly
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values = c("#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#0066FF", "#6633FF",#Gard, Lacto, Mega, Mob, Prev
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/early/Genus_early_BV_eth2.pdf", Genus_early_plot_facet_BV_eth_blank2, height=5, width=10)
######################### 

