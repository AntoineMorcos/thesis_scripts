# Clear environment
rm(list=ls())

# Load libraries
library(tidyverse)
library(stats) 
library(lubridate)

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

# Load medscinet export
INSIGHT_et_al <- read.csv("data/medscinet_export_for_ML_2022-11-22.csv", na.strings = c("", "NA"), sep="|")






###########################################################################################################################
##                                                                                                                       ##
##                                          Keeping only INSIGHT women                                                   ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(INSIGHT_et_al)
dim(INSIGHT_et_al) #14555   125

# Only keep INSIGHT women
raw_INSIGHT <- INSIGHT_et_al %>% 
  dplyr::filter(Enrolled.into.INSIGHT == "Yes" | Enrolled.into.INSIGHT.pH.study =="Yes" | Enrolled.into.INSIGHT.COVID.19.study == "Yes")

# Look at overlap enrollment in INSIGHT & its subdivisions
dim(raw_INSIGHT) #6853  125
table(raw_INSIGHT$Enrolled.into.INSIGHT, raw_INSIGHT$Enrolled.into.INSIGHT.pH.study, useNA="always")
table(raw_INSIGHT$Enrolled.into.INSIGHT, raw_INSIGHT$Enrolled.into.INSIGHT.COVID.19.study, useNA="always")

# Look at Enrolled.into.INSIGHT vs INSIGHT 
table(raw_INSIGHT$Enrolled.into.INSIGHT, raw_INSIGHT$INSIGHT, useNA="always")





###########################################################################################################################
##                                                                                                                       ##
##                                              Looking at raw                                                           ##
##                                                                                                                       ##
###########################################################################################################################


########################################################### Overall df 
# How many unique IDs?
n_distinct(raw_INSIGHT$Participant.ID) #2084 IDs


# Looking
dim(raw_INSIGHT) #6853  125
names(raw_INSIGHT)

###########################################################


########################################################### Looking at data columns 
# Looking
table(raw_INSIGHT$Current.medications..at.time.of.enrolment....Other, useNA="always")
#table(raw_INSIGHT$Current.medications..at.time.of.enrolment....Other.details, useNA="always") #fills up console
#table(raw_INSIGHT$Other.reason.details, useNA="always") #onset of labour I think 
table(raw_INSIGHT$Maternal.pyrexia, useNA="always")
table(raw_INSIGHT$Maternal.Death, useNA="always")
table(raw_INSIGHT$Neonatal.death, raw_INSIGHT$Maternal.Death, useNA="always")
table(raw_INSIGHT$Labour.augmented, useNA="always")
table(raw_INSIGHT$Maternal.Infection, useNA="always") # this is just as a reason for if labour was induced not infection in general
summary(raw_INSIGHT$WCC.highest.value.recorded, useNA="always")
table(raw_INSIGHT$Ethnicity, useNA="always") #Renamed as Previous.PCN.ethnicity in this export
table(raw_INSIGHT$Amniocentesis, useNA="always")
table(raw_INSIGHT$X10Q.code.reading, useNA="always") #lots of random letters not just numbers... may leave out ffn for now
table(raw_INSIGHT$Centre, useNA="always") 
table(raw_INSIGHT$Autoimmune.disease, useNA="always") 

# pH
summary(raw_INSIGHT$Low.in.the.vagina.pH, useNA="always")
summary(raw_INSIGHT$High.in.the.vagina.pH, useNA="always")

# CL
summary(raw_INSIGHT$Transvaginal.cervical.length.mean) #2165 NAs
table(raw_INSIGHT$Cervical.length..25mm.this.pregnancy, useNA="always") #0NAs

# Look at possible duplicate cols 
table(raw_INSIGHT$Pregnancy.Outcome, raw_INSIGHT$Pregnancy.Outcome.Status, useNA="always") #different info in cols 

# Miscarriage numbers
summary(raw_INSIGHT$No..of.pregnancies.ending.13..0.weeks.or.earlier)
summary(raw_INSIGHT$No..of.pregnancies.ending.14..0..23..6.weeks)

# Look at stitch/cerclage info
table(raw_INSIGHT$Cerclage...Type, useNA="always")
table(raw_INSIGHT$Cerclage...Indication, useNA="always")
table(raw_INSIGHT$Cerclage...GA.inserted..d., useNA="always")
table(raw_INSIGHT$Cerclage...GA.inserted..w., useNA="always")
table(raw_INSIGHT$Vaginal.suture.in.situ, useNA="always")
table(raw_INSIGHT$Abdominal.cerclage.in.situ, useNA="always")
# Intervention info from Rachel: "Vaginal stitch in Situ means they have a stitch in place at that visit.
# Abdominal stitch is a high cervical stitch but usually done outside of pregnancy in women that have had previous losses"

# Ethnicity info - this has changed since the last export
dim(raw_INSIGHT)
eth_test <- raw_INSIGHT %>% dplyr::select(Previous.PCN.ethnicity, Ethnicity, Country.of.Birth..if.Other.or.Mixed.Multiple.Ethnicity.please.also.describe.here., 
                                          Partner.s.Ethnicity, Other.or.Mixed.Multiple.Ethnicity.details..partner.)
summary(eth_test) #all NAs for partner, only "Previous.PCN.ethnicity" contains data
eth_test <- NULL
###########################################################

########################################################### Make Date (of visit) to correct data type
# Date of visit
table(raw_INSIGHT$Date.of.visit, useNA="always")
str(raw_INSIGHT$Date.of.visit)
raw_INSIGHT$Date <- substring(raw_INSIGHT$Date.of.visit, 1,10) #keep characters between 1-10 inclusive
raw_INSIGHT$Date <- as.Date(raw_INSIGHT$Date, format="%d/%m/%Y")
raw_INSIGHT %>% dplyr::select(Date.of.visit, Date)
str(raw_INSIGHT$Date)
summary(raw_INSIGHT$Date) 
###########################################################


################################ Preeclampsia & Autoimmune.disease
# Raw data
table(raw_INSIGHT$Preeclampsia, useNA="always") #After Nov medscinet export I only kept the correct preeclampsia column

# Sent Delphine ML_clinical\data\data_to_send_people\Preeclampsia_discrepancies.csv & Autoimmune_disease_discrepancies.csv
# Notes from meeting 27/10/22 with Delphine:
# • Some women have existing hypertension, not preeclampsia, which is why they were induced
# • ID X	has preeclampsia
# • Known.autoimmune.disease not found in individual women's medscinet profiles. So will ignore this col
# • Preeclampsia.1 is more reliable. As "Preeclampsia" in "Delivery Details" was detailing if preeclampsia impacted delivery. Whereas "Preeclampsia" in "End Report" [Preeclampsia.1] is just whether the woman had preeclampsia at some stage in pregnancy

# Checking that the women were corrected as per Delphine's spreadsheet at C:\Users\alici\OneDrive - King's College London\PhD\Projects\ML_clinical\data\data_from_others\Delphine_Preeclampsia_discrepancies.xlsx
pd <- read.csv("data/data_to_send_people/Preeclampsia_discrepancies.csv")
str(pd)
pd$Participant.ID <- as.character(pd$Participant.ID)
unique(raw_INSIGHT %>% dplyr::select(Participant.ID, Preeclampsia) %>% dplyr::filter(Participant.ID %in% pd$Participant.ID))
#these all match the spreadsheet Delphine sent me thankfully, so I know my latest medscinet export has the corrected data too


################################





###########################################################################################################################
##                                                                                                                       ##
##                                         Adding columns for models                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Add Gestation.at.delivery.wks.dec col
medscinet1 <- raw_INSIGHT %>% 
  dplyr::mutate(Gestation.at.delivery.wks.dec = ((Gestation.at.delivery..w.*7) + Gestation.at.delivery..d.) / 7) #number of weeks x7 plus number of days, divided by 7 again to put it back into weeks 
  

# Making BMI_category column 
# BMI category definitions as per https://digital.nhs.uk/data-and-information/publications/statistical/statistics-on-obesity-physical-activity-and-diet/statistics-on-obesity-physical-activity-and-diet-england-2019/part-3-adult-obesity & https://www.nhs.uk/common-health-questions/lifestyle/what-is-the-body-mass-index-bmi/
medscinet1$BMI_category <- NA
medscinet1$BMI_category[medscinet1$BMI<18.5] <- "Underweight"
medscinet1$BMI_category[medscinet1$BMI>=18.5 & medscinet1$BMI<25] <- "Healthy weight"
medscinet1$BMI_category[medscinet1$BMI>=25 & medscinet1$BMI<30] <- "Overweight"
medscinet1$BMI_category[medscinet1$BMI>=30 & medscinet1$BMI<40] <- "Obese"
medscinet1$BMI_category[medscinet1$BMI>=40 ] <- "Morbidly obese"
str(medscinet1$BMI_category)
medscinet1$BMI_category <- factor(medscinet1$BMI_category, levels=c("Underweight", "Healthy weight", "Overweight", "Obese", "Morbidly obese"), ordered=T) 
table(medscinet1$BMI_category, useNA="always")
tapply(medscinet1$BMI, medscinet1$BMI_category, summary)


################################ Simplify ethnicity columns
table(medscinet1$Previous.PCN.ethnicity, useNA="always")
medscinet1$Ethnicity1 <- medscinet1$Previous.PCN.ethnicity
medscinet1$Ethnicity2 <- medscinet1$Previous.PCN.ethnicity

# Ethnicity1
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="0 - Unclassified (other)"] <- "Other"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="1 - European"] <- "European"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="2 - Indian"] <- "South Asian"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="3 - Pakistani"] <- "South Asian"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="4 - Bangladeshi"] <- "South Asian"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="5 - AfroCaribbean"] <- "African Caribbean"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="6 - African"] <- "African"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="7 - Middle Eastern"] <- "Other"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="8 - Far East Asian"] <- "East Asian"
medscinet1$Ethnicity1[medscinet1$Ethnicity1=="9 - South East Asian"] <- "South East Asian"
medscinet1$Ethnicity1 <- factor(medscinet1$Ethnicity1, levels=c("European", "African", "African Caribbean", "South Asian", "East Asian", "South East Asian", "Other")) 
table(medscinet1$Previous.PCN.ethnicity, medscinet1$Ethnicity1, useNA="always")
table(medscinet1$Ethnicity1, useNA="always")
#Note: When I had low & high risk women together I had Middle Eastern as its own category, but now that I'm doing the analysis stratified by risk, I have moved Middle Eastern into "Other" as the numbers are too few
#Before re-exporting data in Nov:
# > table(LR_clin$PTB37, LR_clin$Ethnicity1)
#         European African African Caribbean South Asian East Asian South East Asian Middle Eastern Other
# Term         529     100                28          36         20               20              5    56
# Preterm       25      16                 5           1          0                4              0     3
#
#> table(HR_clin$PTB37, HR_clin$Ethnicity1) 
# 
#         European African African Caribbean South Asian East Asian South East Asian Middle Eastern Other
# Term         643     130                60          20         18                7              4    32
# Preterm      125      46                16          13          2                5              1    11


# Ethnicity2
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="0 - Unclassified (other)"] <- "Other"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="1 - European"] <- "White"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="2 - Indian"] <- "Asian"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="3 - Pakistani"] <- "Asian"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="4 - Bangladeshi"] <- "Asian"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="5 - AfroCaribbean"] <- "Black"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="6 - African"] <- "Black"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="7 - Middle Eastern"] <- "Other"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="8 - Far East Asian"] <- "Asian"
medscinet1$Ethnicity2[medscinet1$Ethnicity2=="9 - South East Asian"] <- "Asian"
table(medscinet1$Previous.PCN.ethnicity, medscinet1$Ethnicity2, useNA="always")
table(medscinet1$Ethnicity2, useNA="always")
medscinet1$Ethnicity2 <- factor(medscinet1$Ethnicity2, levels=c("White", "Black", "Asian", "Other")) 
################################


################################ Making Centre col better 
table(medscinet1$Centre, useNA="always")
medscinet1$Centre[medscinet1$Centre=="STH"] <- "St Thomas' Hospital"
medscinet1$Centre[medscinet1$Centre=="PMH"] <- "Poole - St Mary's Maternity Hospital"
medscinet1$Centre[medscinet1$Centre=="MAN"] <- "Manchester - Saint Mary's Hospital"
medscinet1$Centre[medscinet1$Centre=="WM"] <- "West Middlesex University Hospital"
medscinet1$Centre <- factor(medscinet1$Centre, levels=c("St Thomas' Hospital", "Poole - St Mary's Maternity Hospital", 
                                                        "Manchester - Saint Mary's Hospital", "West Middlesex University Hospital"))
table(medscinet1$Centre, useNA="always")
################################



###########################################################################################################################
##                                                                                                                       ##
##                           Selecting feature/predictor candidates & labels/outcomes for models                         ##
##                                                                                                                       ##
###########################################################################################################################

# Looking
names(medscinet1)[1:30]
names(medscinet1)[31:60]
names(medscinet1)[61:90]
names(medscinet1)[91:130]

# Selecting only certain columns - we can filter more out later
clin_unclean <- medscinet1 %>% 
  dplyr::select(
    
    # ID
    Participant.ID, 
    
    # PTB Outcome
    Gestation.at.delivery.wks.dec,
    Spontaneous.onset.of.labour.resulting.in.delivery..37.40,
    Spontaneous.onset.of.labour.resulting.in.delivery..34.40, # early entries are NAs for <34
    Premature.Prelabour.Rupture.of.Membranes,
    
    # Pregnancy outcome
    Mode.of.delivery,
    Pregnancy.Outcome,
    Pregnancy.Outcome.Status,
    Onset.of.labour, 
    
    # Pregnancy complications - other (helpful for removing iPTBs later?)
    Preeclampsia,
    Gestational.Diabetes,
    
    # Other outcomes - Relevant, but removed as they are outcomes we are not looking at
    Chorioamnionitis,
    #Pre.labour.ruptured.membranes, #we're not looking at term PROM only PPROM (& I've included Premature.Prelabour.Rupture.of.Membranes above)
    Apgar.Score.1.min,
    Apgar.Score.5.min,
    Major.congenital.abnormality,
    Neonatal.death,
    Maternal.Death,
    Customised.birthweight.centiles,
    Suspected.fetal.growth.restriction,
    SGA,
    
    # Maternal demographics
    Age.at.registration, 
    BMI, BMI_category, # will need to pick 1 of 2 for model not both
    Low.risk.at.enrolment,
    Ethnicity1, Ethnicity2, # will need to pick 1 of 2 for model not both
    Smoking,
    Lower.super.output.area, #for IMD later
    Centre, #hospital
    Primigravida,
    Gender,
    
    
    # Previous history - High/Low risk definitions
    Previous.spontaneous.PTB...37.weeks., 
    Previous.PPROM...37.weeks.,
    Previous.late.miscarriage..16.23.6.weeks.,
    Previous.cervical.surgery..e.g..LLETZ..Cone.,
    Cervical.length..25mm.this.pregnancy,
    Transvaginal.cervical.length.mean,
    
    # Risk factors & possible risk factors
    Uterine.abnormality,
    History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy,
    Past.or.present.history.of.GBS,
    Past.or.present.history.of.BV,
    Past.or.present.history.of.Domestic.Violence,
    Past.or.present.history.of.recreational.drug.use,
    Bicornuate.uterus,
    Double.cervix,
    Intra.uterine.Septum,
    Submucosal.fibroids,
    Amniocentesis,
    
    # Number of previous miscarriages https://www.tommys.org/baby-loss-support/miscarriage-information-and-support/causes-miscarriage https://www.tommys.org/baby-loss-support/miscarriage-information-and-support/types-of-miscarriage
    No..of.pregnancies.ending.13..0.weeks.or.earlier, # Early miscarriages 
    No..of.pregnancies.ending.14..0..23..6.weeks, # Late miscarriages
    
    # Interventions
    Did.she.receive.progesterone.,
    Cerclage...Type,
    Vaginal.suture.in.situ,
    Abdominal.cerclage.in.situ,
    
    # Health status 
    Pre.existing.hypertension,
    Asthma,
    Type.1.diabetes,
    Type.2.diabetes,
    Autoimmune.disease,
    Chronic.renal.disease,
    Chronic.viral.infection,
    
    # # Medication - removed because if these are flagged, it's probably the health problem, not the medicine that's involved with PTB
    # Antihypertensives,
    # Steroids,
    # Immunosuppressive.agents,
    # Antibiotics,

    # pH
    Low.in.the.vagina.pH,
    High.in.the.vagina.pH,
    Date, #to remove pHs that were measured incorrectly
    
    # Reproductive infections
    Diagnosed.with.BV, #Diagnosed.with.BV...Was.she.treated, Antibiotic, Antibiotic...Indication
    
    # Quantitative measurements in labour
    CRP.highest.value.recorded,
    WCC.highest.value.recorded
    )


# Looking 
dim(medscinet1) #6853  130
dim(clin_unclean) #6853   67
names(clin_unclean)


# Renaming some cols
clin_unclean <- clin_unclean %>% 
    dplyr::rename("Age"="Age.at.registration",
                  "sPTB37"="Spontaneous.onset.of.labour.resulting.in.delivery..37.40",
                  "sPTB34"="Spontaneous.onset.of.labour.resulting.in.delivery..34.40",
                  "Progesterone"="Did.she.receive.progesterone.",
                  "CRP"="CRP.highest.value.recorded",
                  "WCC"="WCC.highest.value.recorded",
                  "Domestic_violence"="Past.or.present.history.of.Domestic.Violence",
                  "Recreational_drugs"="Past.or.present.history.of.recreational.drug.use",
                  "Previous.sPTB37"="Previous.spontaneous.PTB...37.weeks.",
                  "Previous.PPROM"="Previous.PPROM...37.weeks.",
                  "Previous.late.miscarriage"="Previous.late.miscarriage..16.23.6.weeks.",
                  "Previous.cervical.surgery"="Previous.cervical.surgery..e.g..LLETZ..Cone.",
                  "Short_cervix_at_enrolment"="Cervical.length..25mm.this.pregnancy",
                  "Cervical_length"="Transvaginal.cervical.length.mean",
                  "History.of.UTIs.in.pregnancy"="History.of.2.or.more..proven..recurrent.UTIs.in.pregnancy",
                  "History.of.GBS"="Past.or.present.history.of.GBS",
                  "History.of.BV"="Past.or.present.history.of.BV",
                  "Early_miscarriages_number"="No..of.pregnancies.ending.13..0.weeks.or.earlier", 
                  "Late_miscarriages_number"="No..of.pregnancies.ending.14..0..23..6.weeks",
                  "LSOA"="Lower.super.output.area"
                  )

# Making IDs from int to characters
clin_unclean$Participant.ID <- as.character(clin_unclean$Participant.ID)









###########################################################################################################################
##                                                                                                                       ##
##                                         Choosing pH & removing faulty data                                            ##
##                                                                                                                       ##
###########################################################################################################################

######################### Only keep valid CRP/WCC
#Look
summary(clin_unclean$CRP) #some 999 values
summary(clin_unclean$WCC) #some 999 values

# Get rid of non-valid CRP/WCC
clin_unclean$CRP[clin_unclean$CRP>998] <- NA
clin_unclean$WCC[clin_unclean$WCC>998] <- NA

# Look
summary(clin_unclean$CRP) 
summary(clin_unclean$WCC)
#########################

######################### Only keep valid pH scores
# Look
summary(clin_unclean$High.in.the.vagina.pH) # the pH we used previously in the 16S paper
summary(clin_unclean$Low.in.the.vagina.pH)

# Renaming pH & deleting others
clin_unclean$pH <- clin_unclean$High.in.the.vagina.pH
clin_unclean$High.in.the.vagina.pH <- NULL
clin_unclean$Low.in.the.vagina.pH <- NULL
#########################






###########################################################################################################################
##                                                                                                                       ##
##                                        Excluding IDs & dealing with weird data                                        ##
##                                                                                                                       ##
###########################################################################################################################

################################################## Missing delivery timing
# Looking at Gestation.at.delivery.wks.dec
summary(clin_unclean$Gestation.at.delivery.wks.dec) #41 NAs

# Getting rid of IDs where we don't know when labour was in gestation
dim(clin_unclean) #6853   66
clin_unclean <- clin_unclean %>% 
  dplyr::filter(Gestation.at.delivery.wks.dec != "NA")
dim(clin_unclean) #6812   66
##################################################

################################################## Excluding early miscarriages
# Looking at Gestation.at.delivery.wks.dec
summary(clin_unclean$Gestation.at.delivery.wks.dec)

# Getting rid of IDs which were early miscarriages (definition used here is births <14 weeks)
dim(clin_unclean) #6812   66
clin_unclean <- clin_unclean %>% 
  dplyr::filter(Gestation.at.delivery.wks.dec >= 14)
dim(clin_unclean) #6809   65
summary(clin_unclean$Gestation.at.delivery.wks.dec)
##################################################



# ################################################## sPTB<34, but apparently not sPTB<37 - now fixed on medscinet
# # Looking at outcome data
# table(clin_unclean$sPTB37, useNA="always") #0 NAs - due to medscinet filtering
# table(clin_unclean$sPTB34, useNA="always") #0 NAs
# table(clin_unclean$sPTB37, clin_unclean$sPTB34, useNA="always") #It was 30 entries who had sPTB<34, but not sPTB<37..... But now I've fixed these w/ Delphine
# table(clin_unclean$sPTB37, clin_unclean$Premature.Prelabour.Rupture.of.Membranes, useNA="always")
# 
# 
# # Look at 7 entries who had sPTB<34, but not sPTB<37
# names(clin_unclean)
# weird_sPTBs <- clin_unclean %>% 
#   filter(sPTB37=="No" & sPTB34=="Yes") %>% 
#   dplyr::select(Participant.ID, Gestation.at.delivery.wks.dec, sPTB37, sPTB34, Pregnancy.Outcome, Pregnancy.Outcome.Status)
# weird_sPTBs #all <34
# weird_sPTBs_IDs <- unique(weird_sPTBs$Participant.ID)
# weird_sPTBs_IDs 
# # All Gestation.at.delivery.wks.dec<34
# # Delphine now has manually edit the DB to make them sPTB<37 too.
# 
# ######### Edit weird_sPTBs data in clin_unclean
# # Edit the weird_sPTBs_IDs
# clin_unclean$sPTB37old <- clin_unclean$sPTB37
# clin_unclean$sPTB37
# clin_unclean$sPTB37[which(clin_unclean$sPTB34 == "Yes")] <- "Yes"
# 
# # Checks
# table(clin_unclean$sPTB37, clin_unclean$sPTB37old, useNA="always") #worked
# clin_unclean %>% 
#   dplyr::filter(Participant.ID %in% weird_sPTBs_IDs) %>% 
#   dplyr::select(Participant.ID, Gestation.at.delivery.wks.dec, sPTB37, sPTB37old, sPTB34)
# 
# # Delete checking column
# clin_unclean$sPTB37old <- NULL
# ######### 
# ##################################################


# ################################################## Preeclampsia_delivery & Preeclampsia.1
# # Preeclampsia is in "End Report" & "Delivery Details"
# table(clin_unclean$Preeclampsia_delivery, clin_unclean$Preeclampsia.1, useNA="always")
# 
# # Preeclampsia discrepancy df 1
# weird_Preeclampsia1 <- clin_unclean %>% 
#   filter(Preeclampsia_delivery=="No" & Preeclampsia.1=="Yes") %>% 
#   dplyr::select(Participant.ID, Preeclampsia_delivery, Preeclampsia.1)
# weird_Preeclampsia1
# weird_Preeclampsia1_uniq <- unique(weird_Preeclampsia1) #$Participant.ID
# weird_Preeclampsia1_uniq 
# 
# # Preeclampsia discrepancy df 2
# weird_Preeclampsia2 <- clin_unclean %>% 
#   filter(Preeclampsia_delivery=="Yes" & Preeclampsia.1=="No") %>% 
#   dplyr::select(Participant.ID, Preeclampsia_delivery, Preeclampsia.1)
# weird_Preeclampsia2
# weird_Preeclampsia2_uniq <- unique(weird_Preeclampsia2) #$Participant.ID
# weird_Preeclampsia2_uniq 
# 
# # Combine
# weird_Preeclampsia_ALL <- rbind(weird_Preeclampsia1_uniq, weird_Preeclampsia2_uniq)
# 
# # Write out for midwives
# write.csv(weird_Preeclampsia_ALL, "data/data_to_send_people/Preeclampsia_discrepancies.csv", row.names = F)
# 
# # Delete duplicate data after chat with Delphine
# clin_unclean$Preeclampsia_delivery <- NULL # & kept Preeclampsia
# 
# clin_unclean <- clin_unclean %>% 
#   dplyr::rename("Preeclampsia"="Preeclampsia.1")
# ##################################################


# ################################################## Autoimmune.disease & Known.autoimmune.disease
# # Look
# table(clin_unclean$Autoimmune.disease, clin_unclean$Known.autoimmune.disease, useNA="always")
# 
# # Autoimmune.disease discrepancy df 1
# weird_Autoimmune.disease1 <- clin_unclean %>% 
#   filter(Autoimmune.disease=="No" & Known.autoimmune.disease=="Yes") %>% 
#   dplyr::select(Participant.ID, Autoimmune.disease, Known.autoimmune.disease)
# weird_Autoimmune.disease1
# weird_Autoimmune.disease1_uniq <- unique(weird_Autoimmune.disease1) #$Participant.ID
# weird_Autoimmune.disease1_uniq 
# 
# # Autoimmune.disease discrepancy df 2
# weird_Autoimmune.disease2 <- clin_unclean %>% 
#   filter(Autoimmune.disease=="Yes" & Known.autoimmune.disease=="No") %>% 
#   dplyr::select(Participant.ID, Autoimmune.disease, Known.autoimmune.disease)
# weird_Autoimmune.disease2
# weird_Autoimmune.disease2_uniq <- unique(weird_Autoimmune.disease2) #$Participant.ID
# weird_Autoimmune.disease2_uniq 
# 
# # Combine
# weird_Autoimmune.disease_ALL <- rbind(weird_Autoimmune.disease1_uniq, weird_Autoimmune.disease2_uniq)
# 
# # Write out for midwives
# write.csv(weird_Autoimmune.disease_ALL, "data/data_to_send_people/Autoimmune_disease_discrepancies.csv", row.names = F)
# 
# # Delete duplicate data after chat with Delphine
# clin_unclean$Known.autoimmune.disease <- NULL # & kept Autoimmune.disease
# ##################################################


################################################## Looking at rows with missing dates
# Looking at rows with missing dates
missing_dates <- clin_unclean[is.na(clin_unclean$Date),]
# View(missing_dates)
missing_dates$pH #all NAs
missing_dates$Cervical_length #all NAs
# So missing date data doesn't impact CL & pH calculations
##################################################







###########################################################################################################################
##                                                                                                                       ##
##                                       Deleting pHs that were taken incorrectly                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Rachel "pH values for samples taken before April 7th 2014 are invalid"

# So pHs taken before 2014-04-07 need to be removed

# Look
dim(clin_unclean) #6809   66
summary(clin_unclean$pH) #2581 NAs

# # Test to Remove pH measurements taken before 2014-04-07
# test <- clin_unclean 
# test$validpH <- "s"
# test$validpH[test$Date<"2014-04-07"] <- "no"
# test$validpH[test$Date>="2014-04-07"] <- "Y"
# table(test$validpH)
# tapply(test$pH, test$validpH, summary)
# test$pH[test$Date<"2014-04-07"] <- NA
# tapply(test$pH, test$validpH, summary)
# summary(test$pH) #2915NAs
# summary(test$Date)
# tapply(test$Date, test$validpH, summary) #correct 

# Remove pH measurements taken before 2014-04-07
summary(clin_unclean$pH) #2581 NAs
clin_unclean$pH[clin_unclean$Date<"2014-04-07"] <- NA
summary(clin_unclean$pH) #2912 NAs
summary(clin_unclean$Date)

# Removing dates from df now as we only want 1 row per ID soon
clin_unclean$Date <- NULL







###########################################################################################################################
##                                                                                                                       ##
##                                       Finding why multiple rows per ID                                                ##
##                                                                                                                       ##
###########################################################################################################################

# Looking
dim(clin_unclean) #6809   64
str(clin_unclean$Participant.ID)

# How many unique IDs?
n_distinct(clin_unclean$Participant.ID) #2056 (Was 1989 IDs 2 months ago)

########################################### IDs that occur >=2 times
# Looking at extra IDs - this doesn't include all instances of the duplicate IDs though, just the extra instances
?duplicated
extra_IDs <- clin_unclean$Participant.ID[duplicated(clin_unclean$Participant.ID)]
head(extra_IDs)
length(extra_IDs) 

# Create df where we have >=2 rows for each Participant.ID 
dup_IDs_df <- clin_unclean[clin_unclean$Participant.ID %in% clin_unclean$Participant.ID[duplicated(clin_unclean$Participant.ID)],]

# Looking
dim(dup_IDs_df) #6408   64
table(dup_IDs_df$Participant.ID)
#View(dup_IDs_df)
# "Cerclage...Type" can differ in 1 ID e.g. ID X. Other IDs are identical across all the columns I kept
# Others will be Cervical_length & pH 
###########################################

########################################### Finding which data differs between IDs occuring >=2 times (It is only the col "Cerclage...Type")
# Remove rows that are complete duplicates of other rows 
dup_IDs_diff_data <- distinct(dup_IDs_df)
dim(dup_IDs_diff_data) #5404   64
#View(dup_IDs_diff_data)

# Create df where we have >=2 rows for each Participant.ID 
dup_IDs_diff_data <- dup_IDs_diff_data[dup_IDs_diff_data$Participant.ID %in% dup_IDs_diff_data$Participant.ID[duplicated(dup_IDs_diff_data$Participant.ID)],]
dim(dup_IDs_diff_data) #5103   64

################### Let's see if it is just "Cerclage...Type" which differs between the same IDs
# Make df without "Cerclage...Type"
dup_IDs_diff_data_noCerInf <- dup_IDs_diff_data
dup_IDs_diff_data_noCerInf$Cerclage...Type <- NULL

# # Remove rows that are complete duplicates of other rows 
# dim(dup_IDs_diff_data_noCerInf) 
# dup_IDs_diff_data_noCerInf <- distinct(dup_IDs_diff_data_noCerInf)
# dim(dup_IDs_diff_data_noCerInf) 
# table(dup_IDs_diff_data_noCerInf$Participant.ID) # old note: 1 of each ID.. Therefore only different column is "Cerclage...Type"!! 
# #07.09.22 This was true before I added in pH & Cervical_length which are often taken at every visit

# Looking more at this column in whole df
table(clin_unclean$Cerclage...Type, useNA="always")
###################
########################################### 


############################ Working on Intervention
# Looking 
names(clin_unclean)
table(clin_unclean$Cerclage...Type, clin_unclean$Progesterone, useNA="always")
table(clin_unclean$Cerclage...Type, clin_unclean$Vaginal.suture.in.situ, useNA="always")
table(clin_unclean$Cerclage...Type, clin_unclean$Abdominal.cerclage.in.situ, useNA="always")      
str(clin_unclean$Cerclage...Type)
table(clin_unclean$Cerclage...Type, useNA="always") 

# Add a generic Cerclage column
clin_unclean$Cerclage_type <- clin_unclean$Cerclage...Type #making a copy of the column
clin_unclean$Cerclage_type[clin_unclean$Cerclage_type == "Transabdominal cerclage"] <- "Cerclage"
clin_unclean$Cerclage_type[clin_unclean$Cerclage_type == "Vaginal cerclage - high (Shirodkar)"] <- "Cerclage"
clin_unclean$Cerclage_type[clin_unclean$Cerclage_type == "Vaginal cerclage - low (McDonald's)"] <- "Cerclage"
clin_unclean$Cerclage_type[clin_unclean$Cerclage_type == "Vaginal cerclage - unknown"] <- "Cerclage"
clin_unclean$Cerclage_type[clin_unclean$Cerclage_type == "Arabin pessary"] <- "Arabin_pessary" 
table(clin_unclean$Cerclage...Type, clin_unclean$Cerclage_type, useNA="always")
clin_unclean$Cerclage...Type <- NULL
############################

#Also, pH & Cervical_length are often taken at multiple visits







###########################################################################################################################
##                                                                                                                       ##
##                                                 1 row per ID                                                          ##
##                                                                                                                       ##
###########################################################################################################################

# Remove rows that are complete duplicates of other rows 
dim(clin_unclean) # 6809   64
clin_unclean <- distinct(clin_unclean)
dim(clin_unclean) # 5771   64

# How many unique IDs?
n_distinct(clin_unclean$Participant.ID) #2056 IDs



############################################ Making binary Cerclage column
# 30.09.2022 Now I want a binary cerclage col & a binary progesterone col, instead of generic intervention
# Looking
table(clin_unclean$Cerclage_type, clin_unclean$Vaginal.suture.in.situ, useNA="always")
table(clin_unclean$Cerclage_type, clin_unclean$Abdominal.cerclage.in.situ, useNA="always")
table(clin_unclean$Abdominal.cerclage.in.situ, clin_unclean$Vaginal.suture.in.situ, useNA="always")

# Make new Cerclage based on multiple columns
?case_when
clin_unclean2 <- clin_unclean %>%
  mutate(Cerclage = case_when(Vaginal.suture.in.situ=="Yes" | Abdominal.cerclage.in.situ=="Yes"  ~ "Yes",
                              Cerclage_type == "Cerclage" | Cerclage_type == "Arabin_pessary"    ~ "Yes",
                              is.na(Vaginal.suture.in.situ) & is.na(Abdominal.cerclage.in.situ) ~ NA_character_,
                              TRUE ~ "No"))
table(clin_unclean2$Cerclage, useNA="always")


# Looking
table(clin_unclean2$Cerclage, clin_unclean2$Vaginal.suture.in.situ, useNA="always") # none missing (top right box)
table(clin_unclean2$Cerclage, clin_unclean2$Abdominal.cerclage.in.situ, useNA="always") # none missing (top right box)
table(clin_unclean2$Cerclage, clin_unclean2$Cerclage_type, useNA="always") # none missing (top right box)
############################################ 


############################################ Eliminating IDs duplicates by intervention
# Delete cols we won't use - this will mean distinct() will now be able to do 1ID/row, apart from pH & Cervical_length
clin_unclean2$Cerclage_type <- NULL
clin_unclean2$Vaginal.suture.in.situ <- NULL
clin_unclean2$Abdominal.cerclage.in.situ <- NULL
names(clin_unclean2)

# Remove rows that are complete duplicates of other rows 
dim(clin_unclean2) # 5771   62
clin_unclean2 <- distinct(clin_unclean2)
dim(clin_unclean2) # 5726   62
n_distinct(clin_unclean2$Participant.ID) #2056 distinct IDs
############################################ 


############################################ From CL & pH at each visit -> 1 entry/ID
# Look
str(clin_unclean2$Cervical_length)
str(clin_unclean2$pH)
names(clin_unclean2)

# Info
?mutate
?transmute

# Do mean pH and min CL by ID
clin_unclean2 <- clin_unclean2 %>% 
  group_by(Participant.ID) %>% #we want to perform the calulations by ID
  mutate(pH_mean=mean(pH, na.rm = TRUE), #this is all high in the vagina pH
         CL_minimum=min(Cervical_length, na.rm = TRUE)) #makes errors but this is just because there is missing data

# Checking
checking_CL_pH_maths <- clin_unclean2 %>% dplyr::select(Participant.ID, pH_mean, pH, CL_minimum, Cervical_length)
#View(checking_CL_pH_maths) #it worked 
summary(clin_unclean2$pH_mean) #now max is only 6.1
summary(clin_unclean2$CL_minimum)

# Change to df
clin_unclean2<- as.data.frame(clin_unclean2)

# Chang Inf and NaNs to NAs
clin_unclean2$CL_minimum[clin_unclean2$CL_minimum==Inf] <- NA
clin_unclean2$pH_mean[clin_unclean2$pH_mean=="NaN"] <- NA

# Erase extra columns
clin_unclean2$pH <- NULL
clin_unclean2$Cervical_length <- NULL

# Remove rows that are complete duplicates of other rows 
dim(clin_unclean2) # 5726   62
clin_unclean2 <- distinct(clin_unclean2)
dim(clin_unclean2) # 2132   62
n_distinct(clin_unclean2$Participant.ID) #2056 distinct IDs
############################################


############################################ Finding the problem IDs
# Create df where we have >=2 rows for each Participant.ID 
IDs_tricky <- clin_unclean2[clin_unclean2$Participant.ID %in% clin_unclean2$Participant.ID[duplicated(clin_unclean2$Participant.ID)],]
# IDs_tricky %>% dplyr::select(Participant.ID, Cerclage) #Cerclage is the problem still. 
# IDs_tricky %>% dplyr::select(Participant.ID, Cerclage) %>% arrange(desc(Cerclage)) #ordered by "Yes", "No", NA
dim(IDs_tricky) # 151  62
IDs_tricky_vector <- unique(IDs_tricky$Participant.ID)

# Arrange df so Cerclage is in the order "Yes", "No", NA
clin_unclean3 <- clin_unclean2 %>% 
  arrange(desc(Cerclage)) %>% 
  dplyr::select(Participant.ID, Cerclage, everything())
clin_unclean3 <- as.data.frame(clin_unclean3)
dim(clin_unclean3) #2056   62

# 1 row per ID for df (Take first instance of the ID from the top-bottom of df) 
clin_unique <- clin_unclean3[!duplicated(clin_unclean3$Participant.ID), ]    # Apply duplicated
dim(clin_unique) #2056   62
n_distinct(clin_unclean2$Participant.ID) #2059 distinct IDs 

# Arrange by ID to check
clin_unique <- clin_unique %>% arrange(Participant.ID)
clin_unique %>% dplyr::select(Participant.ID, Cerclage) %>% filter(Participant.ID %in% IDs_tricky_vector) #matches the above, yes overiding no & no overiding NA
############################################ 



############################################ Make short cervix col
# # Sanity check for CL & short cervix 
# boxplot(CL_minimum ~ Short_cervix_at_enrolment, data=clin_unclean2)
# dev.off()
# tapply(clin_unclean2$CL_minimum, clin_unclean2$Short_cervix_at_enrolment, summary) # trend in the right direction, but still not okay for many IDs
# # $No
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# #     0.0    28.7    33.7    32.3    38.0    62.7     898 
# # 
# # $Yes
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# #    0.00   21.00   23.70   23.99   27.00   39.30       1 # Find out what exactly binary Short_cervix_at_enrolment is based on. Delphine says it's a short cervix at recruitment which is why the sanity check failed

# Make column for short cervix during pregnancy
clin_unique$Short_cervix <- NA
clin_unique$Short_cervix[clin_unique$CL_minimum<25] <- "Yes"
clin_unique$Short_cervix[clin_unique$CL_minimum>=25] <- "No"

# Check
head(as.character(clin_unique$CL_minimum), 20)
head(clin_unique$Short_cervix, 20)
tapply(clin_unique$CL_minimum, clin_unique$Short_cervix, summary)
table(clin_unique$Short_cervix, useNA="always") #935 NAs
############################################ 






###########################################################################################################################
##                                                                                                                       ##
##                              Labelling PTBs - sPTB & iPTBs - but not guaranteed accuracy                              ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(clin_unique)
tapply(clin_unique$Gestation.at.delivery.wks.dec, clin_unique$sPTB37, summary)
dim(clin_unique) # 2056   63

# ######################### Looking at women labelled as sPTB37 but born >=37 weeks - now no one as Debbie fixed it
# # Identify fake data (IDs classified as sPTB but actually born >=37 weeks)
# fake_sPTB37_df <- clin_unique %>% 
#   dplyr::filter(sPTB37=="Yes" & Gestation.at.delivery.wks.dec>=37) %>% 
#   dplyr::select(Participant.ID, Gestation.at.delivery.wks.dec, sPTB37)
# 
# # Looking
# fake_sPTB37_df
# tapply(fake_sPTB37_df$Gestation.at.delivery.wks.dec, fake_sPTB37_df$sPTB37, summary)
# fake_sPTB37_IDs <- fake_sPTB37_df$Participant.ID
# fake_sPTB37_IDs 
# #on the DB they are all from PMH centre (Poole maternity hospital)
# 
# # Write out for midwives
# write.csv(fake_sPTB37_df, "data/data_to_send_people/labelled_as_sPTB37s_but_delivered_after_37weeks.csv", row.names = F)
# 
# # Remove fake_sPTB37_IDs from df
# clin1 <- clin_unique %>% dplyr::filter(!Participant.ID %in% fake_sPTB37_IDs)
# dim(clin1) #1984   73
# tapply(clin1$Gestation.at.delivery.wks.dec, clin1$sPTB37, summary)
# #########################


######################### Look at possible iPTBs
# Identify iPTBs
iPTBs_df <- clin_unique %>% dplyr::filter(sPTB37=="No" & Gestation.at.delivery.wks.dec<37)
dim(iPTBs_df) #90
tapply(iPTBs_df$Gestation.at.delivery.wks.dec, iPTBs_df$sPTB37, summary)
iPTB_IDs <- iPTBs_df$Participant.ID
#########################

######################### Make new column for term/preterm & for sPTB/iPTB/TB
clin2 <- clin_unique %>% 
  dplyr::mutate(PTB37=case_when(Gestation.at.delivery.wks.dec <37 ~ "Preterm",
                                Gestation.at.delivery.wks.dec >=37 ~ "Term",
                                TRUE ~ "NA"),
                PTB_type=case_when(Gestation.at.delivery.wks.dec <37 & sPTB37=="Yes" ~ "sPTB",
                                   Gestation.at.delivery.wks.dec <37 & sPTB37=="No" ~ "iPTB",
                                   Gestation.at.delivery.wks.dec >=37 ~ "TB",
                                   TRUE ~ "NA")) %>% 
  dplyr::select(Participant.ID, PTB37, PTB_type, everything())

# Checking
tapply(clin2$Gestation.at.delivery.wks.dec, clin2$PTB37, summary)
table(clin2$PTB37, useNA="always")
table(clin2$PTB_type, useNA="always")
table(clin2$PTB37, clin2$PTB_type, useNA="always") 
table(clin2$PTB37, clin2$sPTB37, useNA="always")
#########################






###########################################################################################################################
##                                                                                                                       ##
##                                                 Clean data columns                                                    ##
##                                                                                                                       ##
###########################################################################################################################

######################### Risk category
# Make clearer labels
clin2$Risk <- clin2$Low.risk.at.enrolment
table(clin2$Risk)
clin2$Risk[clin2$Risk=="Yes"] <- "Low"
clin2$Risk[clin2$Risk=="No"] <- "High"
table(clin2$Risk, clin2$Low.risk.at.enrolment)
clin2$Low.risk.at.enrolment <- NULL
#########################


######################### Unknowns to NAs in different columns
# Gender
table(clin2$Gender, useNA = "always")
clin2$Gender[clin2$Gender=="Indeterminate"] <- NA
clin2$Gender[clin2$Gender=="Unknown"] <- NA
table(clin2$Gender, useNA = "always")

# Domestic_violence
table(clin2$Domestic_violence, useNA = "always")
clin2$Domestic_violence[clin2$Domestic_violence=="Unknown"] <- NA
table(clin2$Domestic_violence, useNA = "always")

# Mode.of.delivery
table(clin2$Mode.of.delivery, useNA = "always")
clin2$Mode.of.delivery[clin2$Mode.of.delivery=="Unknown"] <- NA
table(clin2$Mode.of.delivery, useNA = "always")

# Pregnancy.Outcome
table(clin2$Pregnancy.Outcome, useNA = "always")
clin2$Pregnancy.Outcome[clin2$Pregnancy.Outcome=="Unknown"] <- NA
table(clin2$Pregnancy.Outcome, useNA = "always")

# Onset.of.labour
table(clin2$Onset.of.labour, useNA = "always")
clin2$Onset.of.labour[clin2$Onset.of.labour=="Unknown"] <- NA
table(clin2$Onset.of.labour, useNA = "always")
#########################



######################### BV
# Looking at nas
table(clin2$Diagnosed.with.BV, useNA="always") 
table(clin2$History.of.BV, useNA="always") # Note that "Past.or.present.history.of.BV" was renamed "History.of.BV" 

# Looking at overlap
table(clin2$History.of.BV, clin2$Diagnosed.with.BV, useNA="always") 
######################### 


######################### Primigravida
# Primigravida
table(clin2$Primigravida, useNA = "always")
clin2$Primigravida[clin2$Primigravida=="Yes"] <- "Primigravida"
clin2$Primigravida[clin2$Primigravida=="No"] <- "Multigravida"
table(clin2$Primigravida, useNA = "always") 
######################### 






###########################################################################################################################
##                                                                                                                       ##
##                                Add Socio-economic deprivation (IMD 2019) data                                         ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(clin2)
dim(clin2) #2056   66
clin2 %>% dplyr::select(Participant.ID, LSOA) %>% head
clin2 %>% dplyr::select(Participant.ID, LSOA) %>% tail
clin2$LSOA
#They all begin with "E01..." apart from "S99999999" & NAs


################################## Make IMD df for merge
# LSOA to IMD data - Flavia told me about this website
# Website that hosts the spreadsheets: https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019
# Direct link to spreadsheet https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/833970/File_1_-_IMD2019_Index_of_Multiple_Deprivation.xlsx
# Then I saved the tab in the xlsx of IMD scores as IMD2019.csv. I only got rid of the commas in the column "Index of Multiple Deprivation (IMD) Rank", and made the data type numeric in excel

# Import into R
IMD2019_gov <- read.csv("data/IMD/IMD2019.csv")

# Look
names(IMD2019_gov)
summary(IMD2019_gov)
str(IMD2019_gov)

# Keep only info we need & rename cols
IMD2019 <- IMD2019_gov %>% 
  dplyr::select(LSOA.code..2011., Index.of.Multiple.Deprivation..IMD..Rank, Index.of.Multiple.Deprivation..IMD..Decile) %>% 
  dplyr::rename("LSOA"="LSOA.code..2011.",
                "IMD_rank"="Index.of.Multiple.Deprivation..IMD..Rank",
                "IMD_decile"="Index.of.Multiple.Deprivation..IMD..Decile")
# IMD of rank 1 is the most deprived LSOA & 32844 is the least deprived LSOA
# IMD in decile 1 is the the most deprived 10% of LSOAs nationally & decile 10 is the least deprived 10% of LSOAs nationally
summary(IMD2019)
head(IMD2019)
##################################


################################## Merge
# Checks
clin2 %>% dplyr::select(Participant.ID, LSOA)
str(clin2$LSOA)
str(IMD2019)
dim(clin2) #2056   66
dim(IMD2019) #32844     3

# Make into df
clin2 <- as.data.frame(clin2)

# Merge
?base::merge
clin_all_info <- base::merge(IMD2019, clin2, by="LSOA",
                   all.x=F, all.y=T)

# Look
dim(clin_all_info) #2056   67
#View(clin_all_info)
names(clin_all_info)
str(clin_all_info)

# Looking at NAs
table(is.na(clin_all_info$LSOA)) #140 NAs
summary(clin_all_info$IMD_rank) #261 NAs
summary(clin_all_info$IMD_decile) #261 NAs
#There are more NAs as some LSOA areas have split into multiple areas. 


# Rearrange df
?relocate 
clin_all_info <- clin_all_info %>% 
  arrange(Participant.ID) %>% #order rows by ID
  relocate(Cerclage, LSOA, IMD_rank, IMD_decile, .after = Risk) #move to after risk
##################################






###########################################################################################################################
##                                                                                                                       ##
##                                       Make outcomes_clin (df of all outcome data)                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(clin_all_info)

######################### Make outcomes_clin
# dfs to keep for interest
outcomes_clin <- clin_all_info %>% 
  dplyr::select(Participant.ID:SGA) %>% 
  dplyr::select(Participant.ID, Gestation.at.delivery.wks.dec, PTB37, PTB_type, sPTB37, sPTB34, 
                Customised.birthweight.centiles, Apgar.Score.1.min, Apgar.Score.5.min, everything())
outcomes_clin$Gestational.Diabetes <- NULL
dim(outcomes_clin) #2056   21
names(outcomes_clin)
outcomes_clin <- as.data.frame(outcomes_clin)
#########################


######################### Factor levels
# Check
summary(outcomes_clin$Gestation.at.delivery.wks.dec) #looks reasonable
str(outcomes_clin)

# Character data to factor data
str(outcomes_clin)
cha_names_outcomes <- names(outcomes_clin %>% dplyr::select(PTB37:sPTB34, Premature.Prelabour.Rupture.of.Membranes:SGA))
outcomes_clin[cha_names_outcomes] <- lapply(outcomes_clin[cha_names_outcomes], factor)
str(outcomes_clin)
summary(outcomes_clin)

# PTB37
summary(outcomes_clin$PTB37)
outcomes_clin$PTB37 <- factor(outcomes_clin$PTB37, levels=c("Term", "Preterm"))
summary(outcomes_clin$PTB37)

# PTB_type
summary(outcomes_clin$PTB_type)
outcomes_clin$PTB_type <- factor(outcomes_clin$PTB_type, levels=c("TB", "sPTB", "iPTB"))
summary(outcomes_clin$PTB_type)

# Mode.of.delivery
summary(outcomes_clin$Mode.of.delivery)
outcomes_clin$Mode.of.delivery <- factor(outcomes_clin$Mode.of.delivery, levels=c("Spontaneous vaginal", "Instrumental vaginal", "Emergency caesarean in labour", "Emergency caesarean prior to onset of labour", "Elective caesarean"))
summary(outcomes_clin$Mode.of.delivery)

# Pregnancy.Outcome
summary(outcomes_clin$Pregnancy.Outcome)
outcomes_clin$Pregnancy.Outcome <- factor(outcomes_clin$Pregnancy.Outcome, levels=c("Live birth", "Stillbirth", "Late miscarriage", "Termination of Pregnancy"))
summary(outcomes_clin$Pregnancy.Outcome)

# Onset.of.labour
summary(outcomes_clin$Onset.of.labour)
outcomes_clin$Onset.of.labour <- factor(outcomes_clin$Onset.of.labour, levels=c("Spontaneous", "Induced", "Pre-labour caesarean", "Other surgical"))
summary(outcomes_clin$Onset.of.labour)
#########################





###########################################################################################################################
##                                                                                                                       ##
##                     Make clin_input (df to pass into splitting up for training/validation/test sets)                  ##
##                                                                                                                       ##
###########################################################################################################################

# Look
names(clin_all_info)


######################### Clinical input data for models (clin_input)
# dfs as model input - getting rid of extra outcome data
clin_input_int <- clin_all_info %>% 
  dplyr::select(!c(PTB_type:Onset.of.labour, Chorioamnionitis:SGA, #outcome data
                   LSOA, #as we now have IMD and LSOA won't go into model
                   Short_cervix_at_enrolment)) #only looking at short cervix after enrolment, not before
names(clin_input_int)
dim(clin_input_int) #2056   47

# Rearrange columns
clin_input <- clin_input_int %>% 
  dplyr::rename("pH"="pH_mean") %>% 
  dplyr::select(Participant.ID, #ID
                PTB37,  #Label
                Age, BMI, pH, CL_minimum, CRP, WCC, Early_miscarriages_number, Late_miscarriages_number, # Numeric/int data (including duplicate BMI data)
                IMD_rank, IMD_decile, #duplicates which I will reduce later - treating as integer data
                BMI_category, Short_cervix, #duplicates which I will reduce later
                Ethnicity1, Ethnicity2, #duplicates which I will reduce later 
                Risk, Previous.sPTB37, Previous.PPROM, Previous.cervical.surgery, Previous.late.miscarriage, #duplicates which I will reduce later
                Uterine.abnormality, Bicornuate.uterus, Double.cervix, Intra.uterine.Septum, Submucosal.fibroids, #duplicates which I will reduce later
                Cerclage, Progesterone, #intervention
                Preeclampsia,
                Gestational.Diabetes, Primigravida, Smoking, Diagnosed.with.BV, History.of.BV, History.of.UTIs.in.pregnancy, History.of.GBS, everything()) # Other factor data

# Look
dim(clin_input) #2056   48
names(clin_input)
clin_input <- as.data.frame(clin_input)
#########################


######################### Correcting data types in clin_input
str(clin_input)

# Integer to numeric data
clin_input$Age <- as.numeric(clin_input$Age)
clin_input$CRP <- as.numeric(clin_input$CRP)
clin_input$WCC <- as.numeric(clin_input$WCC)
str(clin_input)

# Character data to factor data
cha_names <- names(clin_input %>% dplyr::select(PTB37, Short_cervix, Risk:Chronic.viral.infection))
cha_names
clin_input[cha_names] <- lapply(clin_input[cha_names], factor)

# Look
str(clin_input)
summary(clin_input)
#########################

#########################  Changing factor levels
# PTB37
summary(clin_input$PTB37)
clin_input$PTB37 <- factor(clin_input$PTB37, levels=c("Term", "Preterm"))
summary(clin_input$PTB37)

# Risk
summary(clin_input$Risk)
clin_input$Risk <- factor(clin_input$Risk, levels=c("Low", "High"))
summary(clin_input$Risk)

# Smoking
summary(clin_input$Smoking)
clin_input$Smoking <- factor(clin_input$Smoking, levels=c("Never", "Ex - gave up before pregnancy", "Ex - gave up in pregnancy", "Current"))
summary(clin_input$Smoking)

# Gender
summary(clin_input$Gender)
clin_input$Gender <- factor(clin_input$Gender, levels=c("Male", "Female"))
summary(clin_input$Gender)

# Look
str(clin_input)
summary(clin_input)
#########################

######################### Finalising
# Making Participant.ID the row names
rownames(clin_input) <- paste("ID_", clin_input$Participant.ID, sep="")

# Duplicate cols note: I still have..
# BMI & BMI_category
# Risk & the criteria which made up risk
# Uterine.abnormality & the criteria which made up Uterine.abnormality
# Number of miscarriages & binary info on previous miscarriages/PTBs
# IMD_decile & IMD_rank
# CL_minimum & Short_cervix
#########################






###########################################################################################################################
##                                                                                                                       ##
##                                   Splitting clin_input into low & high risk women                                     ##
##                                                                                                                       ##
###########################################################################################################################

################################################################ Low risk women
# Make df of only low risk women
LR_clin_int <- clin_input %>% 
  filter(Risk=="Low")

# Looking
dim(LR_clin_int) #887  48
table(LR_clin_int$Risk, useNA="always")
table(LR_clin_int$PTB37, useNA="always")
54/887 #PTB rate of ~6%
summary(LR_clin_int)

# Deleting some columns
LR_clin <- LR_clin_int %>% 
  dplyr::select(!c(CRP, WCC, #data gained in labour 
                   Risk, #all low risk
                   CL_minimum, Short_cervix, # Not meant to be measured/present in low risk women 
                   Late_miscarriages_number, Previous.sPTB37, Previous.PPROM, Previous.cervical.surgery, Previous.late.miscarriage, # Not meant to be measured/present in low risk women 
                   Uterine.abnormality, Bicornuate.uterus, Double.cervix, Intra.uterine.Septum, Submucosal.fibroids, # Not meant to be measured/present in low risk women 
                   Chronic.renal.disease # only 1 ID is different from the rest
  ))
#note: I didn't delete data gained in labour from clin_input as we plot them in NAs_in_df.R

# Looking
summary(LR_clin)
dim(LR_clin) #887  32
table(LR_clin$PTB37, LR_clin$Ethnicity1, useNA="always")
################################################################



################################################################ High risk women
# Make df of only High risk women
HR_clin_int <- clin_input %>% 
  filter(Risk=="High")

# Looking
dim(HR_clin_int) #1169   48
table(HR_clin_int$Risk, useNA="always")
table(HR_clin_int$PTB37, useNA="always")
223/1169 #PTB rate of ~19%
summary(HR_clin_int)

# Deleting some columns
HR_clin <- HR_clin_int %>% 
  dplyr::select(!c(CRP, WCC, #data gained in labour 
                   Risk)) #all High risk
#note: I didn't delete data gained in labour from clin_input as we plot them in NAs_in_df.R

# Looking
summary(HR_clin)
dim(HR_clin) #1169   45
table(HR_clin$PTB37, HR_clin$Ethnicity1) 
################################################################





###########################################################################################################################
##                                                                                                                       ##
##                                                 Exporting data                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Saving data
save(clin_input, # input for models - high & low risk women
     LR_clin, HR_clin, # separate low risk df & high risk df
     outcomes_clin, # extra outcome data
     file="data/1Data_cleaning_output.RData")



# # Loading Rdata from the script "1Data_cleaning.R"
# lnames <- load(file="data/1Data_cleaning_output.RData")
# lnames






