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
library(funrar) #install.packages("funrar")

##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
lnames = load(file="data_BV/main_dfs.RData") 
lnames




############################################################################################
##                                                                                        ##
##                 Late samples (16-24 weeks) making abundance matrix                    ##
##                                                                                        ##
############################################################################################



########################################Making matrix of raw OTU numbers from data_late
#find limits of OTUs in main df
data_late[1:5, 3175:3176] #ends at 3175
data_late[1:5, 235:237] #starts at 235

##Making a OTU df
raw_OTUs_late <- data_late[,c(235:3175)]
raw_OTUs_late[1:5,1:5]

###Deleting rows where there is no OTU data 
dim(raw_OTUs_late) #302
raw_OTUs_late <- raw_OTUs_late %>% 
  drop_na("Otu00001_Lactobacillus_crispatus") 
dim(raw_OTUs_late) #293


#checks
dim(raw_OTUs_late)
names(raw_OTUs_late[1:6])
names(raw_OTUs_late[2936:2941])
raw_OTUs_late[1:5,1:5]

#checking data in numerical
str(raw_OTUs_late[1:6])

#making a matrix
dim(raw_OTUs_late)
raw_OTUs_late[1:5,1:4]
raw_OTUs_late_matrix <- as.matrix(raw_OTUs_late)
dim(raw_OTUs_late_matrix)
raw_OTUs_late_matrix[1:5,1:4]
########################################



##Looking at OTU read data - number of reads per sample across the raw_OTUs 
#Late read numbers
OTUs_row_sum_late <- rowSums(raw_OTUs_late[,])
summary(OTUs_row_sum_late)
boxplot(OTUs_row_sum_late)



######################################## Making relative abundance matrix
## Making a relative abundance matrix of the OTUs https://rdrr.io/cran/funrar/man/make_relative.html#heading-3
?funrar::make_relative

#preview raw matrix
raw_OTUs_late_matrix[1:10,1:4]

#make relative abundance matrix command
AM_OTUs_late <- make_relative(raw_OTUs_late_matrix)

#preview abundance matrix
AM_OTUs_late[1:10,1:4]
########################################


######################################## Saving as RData
###Saving the data from this script
save(AM_OTUs_late,
     file="data_BV/AM_OTUs_late.RData")
########################################




############################################################################################
##                                                                                        ##
##                    Late samples (10-15 weeks) prep for plots                          ##
##                                                                                        ##
############################################################################################

#save df under new name to not change original df
df_late <- as.data.frame(AM_OTUs_late)

##Get rid of the Otu00000_ part of the OTU columns to make aggregating easier
names(df_late) <- substring(names(df_late), 10) #keep characters from the 10th character onwards

##Getting rid of species part of the name
names(df_late) <- gsub("\\_.*","", names(df_late))

str(df_late)

##merging columns with the same names
df_late_merged <- t(rowsum(t(df_late), group = colnames(df_late), na.rm = T))

##Checks and conversion to a df
dim(df_late_merged)
str(df_late_merged)
df_late_merged <- as.data.frame(df_late_merged)

##Add sample ID column
df_late_merged$Sample_ID <- rownames(df_late_merged) 
df_late_merged[1:10, 73:75]

############Deleting genus columns where there is not at least 1 person with an abundance of >10% for that genus & add an others column
dim(df_late_merged) #293  75

#df_late_merged_filtered <- df_late_merged %>% 
#  select_if(~max(., na.rm = TRUE) >= 0.05) %>%     #####At 5%
#  mutate(Other = 1-(rowSums(.[1:20], na.rm=TRUE))) #20 for number of rows (21) minus the sample row

df_late_merged_filtered <- df_late_merged %>% 
  select_if(~max(., na.rm = TRUE) >= 0.1) %>%     #####At 10%
  mutate(Other = 1-(rowSums(.[1:15], na.rm=TRUE)))  #########################1:15 is all rows after the filtering from the previous line of 10% minus the sample ID row

dim(df_late_merged_filtered) #293  17
str(df_late_merged_filtered)

##Making the other col numeric
df_late_merged_filtered$Other <- as.numeric(df_late_merged_filtered$Other)


##Making df of just clinical things of interest
clinical_late <- data_late %>% 
  select(Sample_ID, participant_id, weeks.at.visit.group, risk.cat, Ethnicity_long, Ethnicity_short, group.pcoa.Erica,
         Age, BMI, BMI_category, Smoking, BV_categories, sPTB_34w, sPTB_37w)


#Changing sPTB_34 so it is understandable when using facet_grid
table(clinical_late$sPTB_34w) 
clinical_late$sPTB_34w_copy <- clinical_late$sPTB_34w
clinical_late$sPTB_34w <- as.character(clinical_late$sPTB_34w)
clinical_late$sPTB_34w[clinical_late$sPTB_34w == "No"] <- "Not sPTB <34"
clinical_late$sPTB_34w[clinical_late$sPTB_34w == "Yes"] <- "sPTB <34"
clinical_late$sPTB_34w <- factor(clinical_late$sPTB_34w, levels=c("Not sPTB <34", "sPTB <34"))
table(clinical_late$sPTB_34w, clinical_late$sPTB_34w_copy)
clinical_late$sPTB_34w_copy <- NULL

#Changing sPTB_34 so it is understandable when using facet_grid
table(clinical_late$sPTB_37w) 
clinical_late$sPTB_37w_copy <- clinical_late$sPTB_37w
clinical_late$sPTB_37w <- as.character(clinical_late$sPTB_37w)
clinical_late$sPTB_37w[clinical_late$sPTB_37w == "No"] <- "Not sPTB <37"
clinical_late$sPTB_37w[clinical_late$sPTB_37w == "Yes"] <- "sPTB <37"
clinical_late$sPTB_37w <- factor(clinical_late$sPTB_37w, levels=c("Not sPTB <37", "sPTB <37"))
table(clinical_late$sPTB_37w, clinical_late$sPTB_37w_copy)
clinical_late$sPTB_37w_copy <- NULL

#Changing sPTB_34 so it is understandable when using facet_grid
table(clinical_late$risk.cat) 
clinical_late$risk.cat_copy <- clinical_late$risk.cat
clinical_late$risk.cat <- as.character(clinical_late$risk.cat)
clinical_late$risk.cat[clinical_late$risk.cat == "Low"] <- "Low risk"
clinical_late$risk.cat[clinical_late$risk.cat == "High"] <- "High risk"
clinical_late$risk.cat <- factor(clinical_late$risk.cat, levels=c("Low risk", "High risk"))
table(clinical_late$risk.cat, clinical_late$risk.cat_copy)
clinical_late$risk.cat_copy <- NULL


#merging dfs
late_spread <- merge(clinical_late, df_late_merged_filtered, by="Sample_ID")

##looking at n numbers
table(late_spread$Ethnicity_short, useNA="always")
table(late_spread$BV_categories, useNA="always")
table(late_spread$risk.cat, useNA="always")
table(late_spread$group.pcoa.Erica, useNA="always")
table(late_spread$sPTB_34w, useNA="always")
table(late_spread$sPTB_37w, useNA="always")

#Making a new column to use for plotting
late_spread_sort <- late_spread %>% 
  mutate(Label_gg=str_c(late_spread$Ethnicity_long, ", ", late_spread$risk.cat, ", Age ", late_spread$Age, ", ID:", late_spread$Sample_ID)) %>% 
  arrange(Ethnicity_long, risk.cat, Age)


# Make new col for ordering the plots
late_spread_sort$Lacto_copy <- late_spread_sort$Lactobacillus

# Look
late_spread_sort 
str(late_spread_sort)
names(late_spread_sort)

##looking at n numbers  for plot
dim(late_spread_sort) #N=293
table(late_spread_sort$Ethnicity_short, useNA="always")
table(late_spread_sort$BV_categories, useNA="always")
293-8
table(late_spread_sort$risk.cat, useNA="always")
table(late_spread_sort$group.pcoa.Erica, useNA="always")
table(late_spread_sort$sPTB_34w, useNA="always")
table(late_spread_sort$sPTB_37w, useNA="always")


#gathering genera into 1 column
?gather
?pivot_longer
late_gather <- pivot_longer(late_spread_sort, 
                            cols=Actinomyces:Other, 
                            names_to="Genus", values_to="Abundance")

dim(late_gather)
str(late_gather)


#putting Other category to the end of the factors 
late_gather$Genus<- fct_relevel(late_gather$Genus, "Other", after = Inf)




############################################################################################
##                                                                                        ##
##                     Plotting late samples (16-24 weeks) FACETING                      ##
##                                                                                        ##
############################################################################################


######################### Ethnicity_long
##plot faceted by ethnicity - blank x axis
Genus_late_ethnicity_blank <- ggplot(late_gather, aes(fill=Genus, y=Abundance, x=Label_gg)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~Ethnicity_long, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance", title="Genus composition in samples at 16-24 weeks' gestation, faceted by ethnicity") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6.5), 
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.25, "cm")) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_ethnicity.pdf", Genus_late_ethnicity_blank, height=2.3, width=14, dpi=600)
######################### 




######################### Risk
##plot faceted by risk - blank x axis
Genus_late_risk_blank <- ggplot(late_gather, aes(fill=Genus, y=Abundance, x=Label_gg)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~risk.cat, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance", title="Genus composition in samples at 16-24 weeks' gestation, faceted by risk category") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6.5), 
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.25, "cm")) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_risk.pdf", Genus_late_risk_blank, height=2.3, width=14, dpi=600)
######################### 




######################### BV
##plot faceted by BV - blank x axis
Genus_late_BV_blank <- ggplot(subset(late_gather, !is.na(BV_categories)), aes(fill=Genus, y=Abundance, x=Label_gg)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~BV_categories, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance", title="Genus composition in samples at 16-24 weeks' gestation, faceted by bacterial vaginosis status") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6.5), 
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.25, "cm")) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_BV.pdf", Genus_late_BV_blank, height=2.3, width=14, dpi=600)
######################### 




######################### sPTB
##plot faceted by sPTB_34w - blank x axis
Genus_late_sPTB_34w_blank <- ggplot(late_gather, aes(fill=Genus, y=Abundance, x=Label_gg)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~sPTB_34w, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance", title="Genus composition in samples at 16-24 weeks' gestation, faceted by sPTB <34 weeks' gestation outcome") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6.5), 
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.25, "cm")) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_sPTB_34w.pdf", Genus_late_sPTB_34w_blank, height=2.3, width=14, dpi=600)

##plot faceted by sPTB_37w - blank x axis
Genus_late_sPTB_37w_blank <- ggplot(late_gather, aes(fill=Genus, y=Abundance, x=Label_gg)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~sPTB_37w, scales = "free", space = "free") +
  labs(x="Samples IDs", y="Relative abundance", title="Genus composition in samples at 16-24 weeks' gestation, faceted by sPTB <37 weeks' gestation outcome") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6.5), 
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.25, "cm")) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_sPTB_37w.pdf", Genus_late_sPTB_37w_blank, height=2.3, width=14, dpi=600)
######################### 












############################################################################################
##                                                                                        ##
##          Plotting late samples (16-24 weeks) FACET GRID  white black only              ##
##                                                                                        ##
############################################################################################

# white black only
late_gather_white_black <- late_gather %>% 
  filter(Ethnicity_long != "Other ethnicity") 

# N numbers
dim(late_spread_sort) #N=293
table(late_spread_sort$BV_categories, late_spread_sort$Ethnicity_short, useNA="always") #
table(late_spread_sort$BV_categories, useNA="always") 
154+25+16+39+15+8 #257
293-22-4-2-1-7 #257


##plot faceted by BV - blank x axis
Genus_late_BV_eth_blank2 <- ggplot(subset(late_gather_white_black, !is.na(BV_categories)), aes(fill=Genus, y=Abundance, x=reorder(Label_gg, -Lacto_copy))) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~Ethnicity_long+BV_categories, scales = "free", nrow=2) +  #  with ,  but then it looks ugly
  labs(x="Samples IDs", y="Relative abundance") +
  guides(fill=guide_legend(title="Genera")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values = c("#663333", "#990000", "#FF0000", "#FF9933", "yellow", "#66FF66", #ACTINO, Aero, Ato, Bifo, Dia, Fuso
                               "#00CC99", "#99CCFF", "#0099FF", "#6633FF", "#9933FF", "#CC99FF", #Gard, Lacto, Mega,  Prev, #SCAR, SHUT #Mob="#0066FF", Mob not in the lates
                               "#9999FF","#FF99FF", "#FF3399", "black")) #Sne, Strep, Ten, Other
ggsave("plots/Genus_abundance/late/Genus_late_BV_eth2.pdf", Genus_late_BV_eth_blank2, height=5, width=10)
######################### 

