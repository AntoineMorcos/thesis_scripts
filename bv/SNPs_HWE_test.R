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
library(HardyWeinberg) #install.packages("HardyWeinberg")


##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs") #pc
lnames = load(file="data_BV/main_dfs.RData") 
lnames



############################################################################################
##                                                                                        ##
##                      Test For Hardy-Weinberg Equilibrium                               ##
##                                                                                        ##
############################################################################################

#################### Using the HardyWeinberg package
#https://cran.r-project.org/web/packages/HardyWeinberg/HardyWeinberg.pdf
#https://www.rdocumentation.org/packages/HardyWeinberg/versions/1.6.3
####################



######################################################################## rs17333103
table(data_early$rs17333103, useNA="always")
rs17333103_counts <- c(CC=202, CT=63, TT=7)

# #Chi-squared test
# HWChisq(rs17333103_counts, x.linked=F, verbose=T, cc=F) # p-value =0.4390878 
# 
# #Exact Test for Hardy-Weinberg Equilibrium
# HWExact(rs17333103_counts, alternative="two.sided", pvaluetype="selome", x.linked=F, verbose=T) #p-value =  0.451552 

#p-values of all tests!!
HWAlltests(rs17333103_counts, x.linked=F, verbose=T, include.permutation.test=T)
#                                            Statistic   p-value
#Chi-square test:                            0.5986650 0.4390878
#Chi-square test with continuity correction: 0.3103869 0.5774429
#Likelihood-ratio test:                      0.5642175 0.4525659
#Exact test with selome p-value:                    NA 0.4515520
#Exact test with dost p-value:                      NA 0.5574784
#Exact test with mid p-value:                       NA 0.3847870
#Permutation test:                           0.5986650 0.4426471
########################################################################




######################################################################## rs1983649
table(data_early$rs1983649, useNA="always")
rs1983649_counts <- c(TT=58, AT=120, AA=95)

# #Chi-squared test
# HWChisq(rs1983649_counts, x.linked=F, verbose=T, cc=F) # p-value =0.0844474
# 
# #Exact Test for Hardy-Weinberg Equilibrium
# HWExact(rs1983649_counts, alternative="two.sided", pvaluetype="selome", x.linked=F, verbose=T) #p-value = 0.08513326

#p-values of all tests!!
HWAlltests(rs1983649_counts, x.linked=F, verbose=T, include.permutation.test=T)
#                                            Statistic    p-value
#Chi-square test:                             2.977159 0.08444740
#Chi-square test with continuity correction:  2.665664 0.10253500
#Likelihood-ratio test:                       2.975578 0.08452990
#Exact test with selome p-value:                    NA 0.08513326
#Exact test with dost p-value:                      NA 0.10231305
#Exact test with mid p-value:                       NA 0.07455957
#Permutation test:                            2.977159 0.10352941
########################################################################




########################################################################rs6032040
table(data_early$rs6032040, useNA="always")
rs6032040_counts <- c(TT=184, AT=87, AA=4)

# #Chi-squared test
# HWChisq(rs6032040_counts, x.linked=F, verbose=T, cc=F) # p-value = 0.07600271
# 
# #Exact Test for Hardy-Weinberg Equilibrium
# HWExact(rs6032040_counts, alternative="two.sided", pvaluetype="selome", x.linked=F, verbose=T) #p-value = 0.09180713

#p-values of all tests!!
HWAlltests(rs6032040_counts, x.linked=F, verbose=T, include.permutation.test=T)
#                                            Statistic    p-value
#Chi-square test:                             3.148373 0.07600271
#Chi-square test with continuity correction:  2.541549 0.11088594
#Likelihood-ratio test:                       3.625859 0.05688821
#Exact test with selome p-value:                    NA 0.09180713
#Exact test with dost p-value:                      NA 0.11074898
#Exact test with mid p-value:                       NA 0.07329002
#Permutation test:                            3.148373 0.09241176
########################################################################




########################################################################rs2664581
table(data_early$rs2664581, useNA="always")
rs2664581_counts <- c(CC=7, AC=68, AA=197)

# #Chi-squared test
# HWChisq(rs2664581_counts, x.linked=F, verbose=T, cc=F) # p-value =  0.6977768 
# 
# #Exact Test for Hardy-Weinberg Equilibrium
# HWExact(rs2664581_counts, alternative="two.sided", pvaluetype="selome", x.linked=F, verbose=T) #p-value = 0.639022 

#p-values of all tests!!
HWAlltests(rs2664581_counts, x.linked=F, verbose=T, include.permutation.test=T)
#                                             Statistic   p-value
#Chi-square test:                            0.15079501 0.6977768
#Chi-square test with continuity correction: 0.03572759 0.8500792
#Likelihood-ratio test:                      0.14658316 0.7018222
#Exact test with selome p-value:                     NA 0.6390220
#Exact test with dost p-value:                       NA 0.8255824
#Exact test with mid p-value:                        NA 0.5567163
#Permutation test:                           0.15079501 0.8113529
########################################################################




