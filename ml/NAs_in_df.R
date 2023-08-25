# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(ggrepel)
library(stats) 
library(caret)


#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/ML_clinical") #mac
setwd("C:/Users/alici/OneDrive - King's College London/PhD/Projects/ML_clinical") #pc
getwd()
#list.files()








###########################################################################################################################
##                                                                                                                       ##
##                                                   Load data                                                           ##
##                                                                                                                       ##
###########################################################################################################################

# Loading Rdata from the script "1Data_cleaning.R"
lnames <- load(file="data/1Data_cleaning_output.RData")
lnames

# Look
dim(clin_input) #2056   48

# Delete ID column
clin_input$Participant.ID <- NULL










###########################################################################################################################
##                                                                                                                       ##
##                                            High risk women only (HR_clin)                                             ##
##                                                                                                                       ##
###########################################################################################################################

# Look
dim(HR_clin) #1169   45
summary(HR_clin)

# Number of NAs_vector 
NAs_HR_vector <- colSums(is.na(HR_clin))
NAs_HR_vector 

# Make into df to get %
NAs_HR_df <- data.frame(Number_of_NAs=NAs_HR_vector)
NAs_HR_df

# Add columns
NAs_HR_df <- NAs_HR_df %>% 
  dplyr::mutate(Percentage_of_NAs=(Number_of_NAs/dim(HR_clin)[1])*100)
NAs_HR_df$Variable <- row.names(NAs_HR_df)
NAs_HR_df


# Remove rows to avoid confusion/duplication
NAs_HR_df <- NAs_HR_df %>% 
  dplyr::filter(!Variable %in% c("Ethnicity2", "BMI_category", "Maternal.Infection", "Short_cervix", "IMD_decile", #duplicates & not used variables
                                 "Participant.ID", "PTB37" )) #no NAs due to filtering
NAs_HR_df


# Change "."s to "_"s  & rename rows
NAs_HR_df$Variable <- gsub('.', ' ', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('_', ' ', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('Ethnicity1', 'Ethnicity', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('Previous sPTB37', 'Previous sPTB', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('Gestational Diabetes', 'Gestational diabetes', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('CL minimum', 'Cervical length', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable <- gsub('Pre existing hypertension', 'Pre-existing hypertension', NAs_HR_df$Variable, fixed = T)
NAs_HR_df$Variable
NAs_HR_df

# Str to factor
str(NAs_HR_df)
NAs_HR_df$Variable <- factor(NAs_HR_df$Variable, levels = NAs_HR_df$Variable)
str(NAs_HR_df)



# Plot percentage of data that is missing in lollipop plot
NAs_plot_HR <- ggplot(NAs_HR_df, aes(x=Percentage_of_NAs, y=reorder(Variable, -Percentage_of_NAs))) +
  geom_point(aes(colour=Percentage_of_NAs), size=3.5) +
  geom_segment(aes(x=0, xend=Percentage_of_NAs, yend=Variable, colour=Percentage_of_NAs), size=2) + 
  geom_text(size=2, color="black", #fontface = "bold",
            aes(label=paste0(round(Percentage_of_NAs, digits=0)#, "%", sep=""
            ))) +
  scale_colour_gradient2(low="#3DCF74", high="red", mid = "yellow", midpoint = 50) +
  xlim(0,30) +
  labs(y=NULL, x="Missing data (%)") +
  theme(legend.position = "none") #rect=element_blank()
  
ggsave("plots/basic/NAs_plot_HR.pdf", NAs_plot_HR, height=6.5, width=5)


# Checking
summary(HR_clin$CL_minimum) #86NAs
86/1169 #~7% #matching the plot 






