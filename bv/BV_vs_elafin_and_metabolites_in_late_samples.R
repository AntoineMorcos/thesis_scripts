#clear environment
rm(list=ls())

getwd()
#list.files()

library(car)
library(multcomp)
library(reshape2)
library(tidyverse)
library(stats) 
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)



#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc










############################################################################################
##                                                                                        ##
##                                     df prep                                            ##
##                                                                                        ##
############################################################################################

## Loading .Rdata saved 
lnames = load(file="data_BV/main_dfs.RData") 
lnames

##Changing BV factor order 
table(data_late$BV_categories)
data_late$BV_categories <- factor(data_late$BV_categories, levels=c("Normal", "Intermediate", "BV")) 
table(data_late$BV_categories)


## df of just metabolites (Excluding Propylene_Glycol AS IT'S AN INGREDIENT USED IN THE INTERNAL SCAN)
metabolites_late <- data_late %>% 
  dplyr::select(participant_id, BV_categories, BV.grading, sPTB_34w, sPTB_37w, Ethnicity_long, Ethnicity_short,
         Leucine,	Isoleucine,	Valine,	Lactate,	Alanine,	Acetate,	Methionine,	
         Glutamate,	Glutamine,	Succinate,	Aspartate,	Asparagine,	Lysine,	Betaine,	Choline,	
         Carnitine,	Threonine,	Formate,	Tyrosine,	Phenylalanine,	Tryptophan,	Serine,	Taurine,	
         Glucose,	Uracil,	Free_EDTA,	CaEDTA,	MgEDTA)

##look
head(metabolites_late)
str(metabolites_late)

###Deleting rows where there is no BV data 
dim(metabolites_late) #302
metabolites_late <- metabolites_late %>% 
  drop_na("BV_categories") 
dim(metabolites_late) #294



################################################## boxplot prep
#gathering
?gather
?pivot_longer
mets_gathered_late <- pivot_longer(metabolites_late, 
                             cols=Leucine:MgEDTA, 
                             names_to="Metabolites", values_to="Quantity")

#looking
head(mets_gathered_late)
str(mets_gathered_late)

#Changing from character to factor
mets_gathered_late$Metabolites <- as.factor(mets_gathered_late$Metabolites)
str(mets_gathered_late)
##################################################







############################################################################################
##                                                                                        ##
##                                   BV vs metabolites                                    ##
##                                                                                        ##
############################################################################################

# Looking
head(mets_gathered_late)
?rstatix::pairwise_wilcox_test
?rstatix::adjust_pvalue
?rstatix::add_significance
?rstatix::add_xy_position
?stat_pvalue_manual

# df of stats for plot - BV_categories fdr
stats_fdr_BV_metabolites_late <- mets_gathered_late %>%
  group_by(Metabolites) %>%
  pairwise_wilcox_test(Quantity ~ BV_categories, p.adjust.method = "fdr") %>%
  add_significance() %>% 
  add_xy_position(x = "BV_categories", scales="free", step.increase=0.3)

# Look
stats_fdr_BV_metabolites_late
table(stats_fdr_BV_metabolites_late$p.adj.signif) #an idea of how many are sig
str(stats_fdr_BV_metabolites_late)
summary(metabolites_late)
tapply(stats_fdr_BV_metabolites_late$p.adj, stats_fdr_BV_metabolites_late$p.adj.signif, summary)


# Checks on p-values 
head(metabolites_late)
ggpubr::compare_means(Acetate ~ BV_categories, data=metabolites_late, method = "wilcox.test", p.adjust.method = "fdr") # matches 


# Boxplot of Metabolites by BV_categories 
late_metabolites_BV_boxplots_thesis <- ggplot(subset(mets_gathered_late, !is.na(BV_categories)), aes(x=BV_categories, y=Quantity)) + 
  geom_boxplot(aes(color=BV_categories)) + 
  facet_wrap(~Metabolites, nrow=4, scales="free_y") +
  labs(x="Bacterial vaginosis status", y="Metabolite concentration") + 
  theme(legend.position = "none", 
        axis.text.x=element_text(size=7),
        strip.text.x = element_text(size = 11)) +
  scale_color_manual(values = c("#00CC00", "#FF9900", "tomato"), na.value = "grey50") +
  stat_pvalue_manual(stats_fdr_BV_metabolites_late, 
                     hide.ns = T, label = "{p.adj.signif}", #fdr
                     size=4) + # font size
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) # Add 10% spaces between the top of the y and the plot border (& 0% space to bottom of y)
ggsave("plots/BV_metabolites_elafin_late/late_metabolites_BV_boxplots_thesis.pdf", late_metabolites_BV_boxplots_thesis, height=7, width=13)

# Ns for plot
stats_fdr_BV_metabolites_late
222+45+27 # 294
dim(metabolites_late) #294

# Looking for write up - look at median
tapply(metabolites_late$Taurine, metabolites_late$BV_categories, summary)
tapply(metabolites_late$Serine, metabolites_late$BV_categories, summary)
tapply(metabolites_late$Succinate, metabolites_late$BV_categories, summary)





############################################################################################
##                                                                                        ##
##                                      BV vs elafin                                      ##
##                                                                                        ##
############################################################################################

#################################################################### All ethnicities together
#Look
table(data_late$BV_categories)
?compare_means

# Find p-value
compare_means(subst.elafin20.x ~ BV_categories, data = data_late, method = "wilcox.test") #ns

# Elafin vs BV for thesis
BV_elafin_late <- ggplot(subset(data_late, !is.na(BV_categories) & !is.na(subst.elafin20.x)), 
                         aes(x=BV_categories, y=subst.elafin20.x)) +
  geom_violin(aes(fill=BV_categories)) + #, outlier.alpha =0.4
  labs(x="BV status", y="Elafin concentration (pg/\u00B5L)") +
  scale_fill_manual(values = c("#00CC00", "#FF9900", "tomato"), na.value="grey50") +
  theme(legend.position = "none") 
ggsave("plots/BV_metabolites_elafin_late/late_BV_elafin.pdf", BV_elafin_late, height=3, width=3)

# Find ns
dim(data_late) #302 IDs
n_elafin_BV <- data_late %>% 
  dplyr::select(Sample_ID, subst.elafin20.x, BV_categories)
n_elafin_BV <- subset(n_elafin_BV, !is.na(BV_categories) & !is.na(subst.elafin20.x)) 
dim(n_elafin_BV) #281 for N
#################################################################### 




#################################################################### stratified by ethnicity
#Look
table(data_late$BV_categories)
str(data_late)

# Cols we want 
BV_elafin_late <- data_late %>% 
  select(participant_id, BV_categories, subst.elafin20.x, 
         Ethnicity_short, Ethnicity_long, sPTB_34w, sPTB_37w,
         BV_categories) 

#Look
table(BV_elafin_late$BV_categories, useNA="always")
str(BV_elafin_late)
summary(BV_elafin_late)
head(BV_elafin_late)
dim(BV_elafin_late) #302

#drop samples without elafin & BV
BV_elafin_late <- BV_elafin_late %>% 
  drop_na(subst.elafin20.x) %>% 
  drop_na(BV_categories)
dim(BV_elafin_late) #281


#drop Other eth
BV_elafin_late_WB <- BV_elafin_late %>% 
  dplyr::filter(Ethnicity_short != "Other") 
dim(BV_elafin_late_WB) #254
table(BV_elafin_late_WB$Ethnicity_short)

#N numbers in plot
dim(BV_elafin_late_WB) #254
table(BV_elafin_late_WB$Ethnicity_short, BV_elafin_late_WB$BV_categories, useNA="always")
table(BV_elafin_late_WB$Ethnicity_short)

# df of stats for plot 
stats_fdr_BV_elafin_late <- BV_elafin_late_WB %>%
  group_by(Ethnicity_long) %>%
  pairwise_wilcox_test(subst.elafin20.x ~ BV_categories, p.adjust.method = "fdr") %>%
  add_significance() %>% 
  add_xy_position(x = "BV_categories", scales="fixed", step.increase=0.3)
stats_fdr_BV_elafin_late 
table(stats_fdr_BV_elafin_late$p.adj.signif) 

# Elafin vs BV 
BV_elafin_late_eth <- ggplot(BV_elafin_late_WB, aes(x=BV_categories, y=subst.elafin20.x)) +
  geom_violin(aes(fill=BV_categories)) + #, outlier.alpha =0.4
  facet_wrap(~Ethnicity_long, nrow=1) + 
  labs(x="BV status", y="Elafin concentration (pg/\u00B5L)") +
  scale_fill_manual(values = c("#00CC00", "#FF9900", "tomato"), na.value="grey50") +
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 11)) +
  stat_pvalue_manual(stats_fdr_BV_elafin_late, 
                     hide.ns = T, label = "{p.adj.signif}", #fdr
                     size=4) + # font size
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) # Add 10% spaces between the top of the y and the plot border (& 0% space to bottom of y)
ggsave("plots/BV_metabolites_elafin_late/late_BV_elafin_eth.pdf", BV_elafin_late_eth, height=3, width=6)
####################################################################
