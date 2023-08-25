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


##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
lnames = load(file="data_BV/main_dfs.RData") 
#setwd("~/OneDrive - King's College London/PhD/Projects/SNPs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/SNPs") #pc
lnames





############################################################################################
##                                                                                        ##
##                                 INSIGHT SNPs df prep                                   ##
##                                                                                        ##
############################################################################################

# Looking
data_early[1:5, 1:5]
dim(data_early)


## Making SNPs - here it is arbitrary whether we use data_early or data_late as we are not using data that changes by timepoint
SNPs <- data_early %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040,
         Ethnicity_short)

# Check
str(SNPs)


# Get rid of samples where we have NAs for ALL SNPs 
dim(SNPs)
SNPs_copy <- SNPs #302
SNPs <- SNPs[!with(SNPs, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(SNPs) # 275

# Looking
summary(SNPs)
# participant_id     rs17333103 rs1983649  rs2664581  rs6032040 Ethnicity_short
# Length:275         CC  :202   TT  : 58   CC  :  7   TT:184    White:188      
# Class :character   CT  : 63   AT  :120   AC  : 68   AT: 87    Black: 59      
# Mode  :character   TT  :  7   AA  : 95   AA  :197   AA:  4    Other: 28      
#                    NA's:  3   NA's:  2   NA's:  3   


#Compare dfs to check NA drop worked well 
subset(SNPs_copy, !(participant_id %in% SNPs$participant_id)) 
subset(SNPs, !(participant_id %in% SNPs_copy$participant_id)) 

#delete df copy
SNPs_copy <- NULL




##data shaping prep
str(SNPs)
SNPs$rs17333103 <- as.character(SNPs$rs17333103)
SNPs$rs1983649 <- as.character(SNPs$rs1983649)
SNPs$rs2664581 <- as.character(SNPs$rs2664581)
SNPs$rs6032040 <- as.character(SNPs$rs6032040)

#Change dataframe for plot
SNPs2 <- SNPs %>% 
  gather(key=SNP, value=Genotype, 
         -participant_id, -Ethnicity_short)

str(SNPs2)

#Changing chr to factors
SNPs2$Genotype <- factor(SNPs2$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA"))
SNPs2$SNP <- factor(SNPs2$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))

#Adding col
SNPs2$Ensembl_or_INSIGHT_data <- "INSIGHT_data"


#Looking
head(SNPs2)


# Remove NAs in SNPs2 to avoid them getting some of the %/proportion out of the SNP genotypes
table(SNPs2$Genotype, useNA="always") #8 NAs
SNPs2 <- SNPs2[!is.na(SNPs2$Genotype),]
table(SNPs2$Genotype, useNA="always") #0 NAs


#Now get summary df of INSIGHT so we can add %s 
SNPs2_summary_basic <- group_by(SNPs2, SNP, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(SNP) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))

# Looking
SNPs2_summary_basic


#Now get summary df of INSIGHT so we can add %s , by ethnicity
SNPs2_summary <- group_by(SNPs2, Ethnicity_short, SNP, Genotype) %>% 
  summarize(n = length(participant_id)) %>% 
  ungroup %>% group_by(SNP, Ethnicity_short) %>% 
  mutate(Proportion = n / sum(n)) %>% 
  mutate(percentage=round(Proportion*100))
SNPs2_summary



#Adding column to basic df so we can have all and separated by ethnicity in one df
SNPs2_summary_basic$Ethnicity_short <- "All"

str(SNPs2_summary_basic)
str(SNPs2_summary)


#Making ethnicity to characters to stop warning when merging
SNPs2_summary$Ethnicity_short <- as.character(SNPs2_summary$Ethnicity_short)

SNPs_ethnicity <- rbind(SNPs2_summary, SNPs2_summary_basic)
str(SNPs_ethnicity)


#ordering factors
SNPs_ethnicity$Ethnicity_short <- factor(SNPs_ethnicity$Ethnicity_short, levels=c("All", "White", "Black", "Other"))
table(SNPs_ethnicity$Ethnicity_short)












############################################################################################
##                                                                                        ##
##                                  Ensembl data prep                                     ##
##                                                                                        ##
############################################################################################

#Now to excel to manually add in genotype theoretical frequencies from Ensembl....
theo <- read.csv("data_SNPs/Ensembl_elafin_SNP_proportions.csv", header=T, fileEncoding="UTF-8-BOM") #first column is called "ï..SNP" without the fileEncoding argument  

str(theo)

#ordering factors
theo$Genotype <- factor(theo$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA"))
theo$Ethnicity <- factor(theo$Ethnicity, levels=c("All", "European", "African"))

#making percentage column
theo <- theo %>% 
  mutate(percentage=round(Proportion*100))








############################################################################################
##                                                                                        ##
##                                 INSIGHT vs Ensembl                                     ##
##                                                                                        ##
############################################################################################

#Filtering out Other ethnicity
SNPs_ethnicity_all_white_black <- SNPs_ethnicity %>% 
  filter(Ethnicity_short!="Other")

#df for plot faceted by SNP and whether its INSIGHT or Ensembl data
#converting from tibble to proper df
SNPs_ethnicity_all_white_black_df <- as.data.frame(SNPs_ethnicity_all_white_black)

#Add study column
SNPs_ethnicity_all_white_black_df$study <- "INSIGHT data"
theo$study <- "Ensembl"

#Delete column we don't want in this plot
SNPs_ethnicity_all_white_black_df$n <- NULL

#renaming ethnicity column so col names match
names(SNPs_ethnicity_all_white_black_df)[names(SNPs_ethnicity_all_white_black_df) == "Ethnicity_short"] <- "Ethnicity"

#merge
Ensembl_INSIGHT_df <- rbind(SNPs_ethnicity_all_white_black_df, theo)
names(Ensembl_INSIGHT_df)


# Plot of proportions
Ensembl_INSIGHT_plot <- ggplot(Ensembl_INSIGHT_df, aes(x=Ethnicity, fill=Genotype, y=Proportion)) + 
  geom_bar(stat="identity") + 
  facet_grid(SNP ~ study, 
             scales = "free", space = "free") +
  labs(x="Ethnicity", y="Proportion") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=8),
        strip.text.x = element_text(size = 9)) + #facet labels
  guides(fill=guide_legend(title="Genotype")) +
  geom_text(aes(x=Ethnicity, label=ifelse(percentage>8, paste0(percentage, "%"),"")), #means label if it is above 8%, otherwise do not label
            colour = "white",  position=position_stack(vjust=0.5), size=3) 
ggsave("plots/genotype_proportions/SNPs_Ensembl_INSIGHT.pdf", Ensembl_INSIGHT_plot, height=6, width=6)






