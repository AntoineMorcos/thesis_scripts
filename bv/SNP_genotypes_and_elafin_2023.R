#clear environment
rm(list=ls())


library(car)
library(multcomp)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(stats) 
library(tidyr)
library(pracma) 
library(scales)
library(ggpubr) #devtools::install_github("kassambara/ggpubr")
library(rstatix)

##Loading .Rdata saved 
#setwd("~/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/BV_vs_OTUs/data_BV") #pc
lnames = load(file="main_dfs.RData") 
#setwd("~/OneDrive - King's College London/PhD/Projects/SNPs") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/SNPs") #pc
lnames


############################################################################################
##                                                                                        ##
##                                      df prep                                           ##
##                                                                                        ##
############################################################################################

# Cols we want 
snp_early <- data_early %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040,
         subst.elafin20.x, weeks.at.visit.group, 
         Ethnicity_short, Ethnicity_long, sPTB_34w, sPTB_37w,
         BV_categories) 

#check
str(snp_early)

##get rid of samples where we have NAs for ALL SNPs 
dim(snp_early) #302   
snp_early <- snp_early[!with(snp_early, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(snp_early) #275     - dropped 27 samples

#drop samples without elafin
snp_early <- snp_early %>% 
  drop_na(subst.elafin20.x)
dim(snp_early) #264  


##data shaping prep
str(snp_early)
snp_early$rs17333103 <- as.character(snp_early$rs17333103)
snp_early$rs1983649 <- as.character(snp_early$rs1983649)
snp_early$rs2664581 <- as.character(snp_early$rs2664581)
snp_early$rs6032040 <- as.character(snp_early$rs6032040)
str(snp_early)


# Gather df for plotting
?gather
?pivot_longer
snp_early_gather <- pivot_longer(snp_early,
                                 cols=rs17333103:rs6032040,
                                 names_to="SNP",
                                 values_to="Genotype")
head(snp_early_gather)


# Remove Genotype NAs in snp_early_gather to avoid them getting some of the %/proportion out of the SNP genotypes
table(snp_early_gather$Genotype, useNA="always") #8 NAs
snp_early_gather <- snp_early_gather[!is.na(snp_early_gather$Genotype),]
table(snp_early_gather$Genotype, useNA="always") #0 NAs



#Changing chr to factors
str(snp_early_gather)
snp_early_gather$Genotype <- factor(snp_early_gather$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA")) # these factor levels give the same outcome when doing the stat test as leaving it as character type (was wondering if excess levels would mess up p-values)
snp_early_gather$SNP <- factor(snp_early_gather$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))
str(snp_early_gather)
summary(snp_early_gather)




############################################################################################
##                                                                                        ##
##                                Genotypes vs elafin                                     ##
##                                                                                        ##
############################################################################################

# Looking
head(snp_early_gather)
?rstatix::pairwise_wilcox_test
?rstatix::adjust_pvalue
?rstatix::add_significance
?rstatix::add_xy_position
?stat_pvalue_manual
?facet_grid

# df of stats for plot - Genotype fdr
stats_SNP_elafin_early <- snp_early_gather %>%
  group_by(SNP) %>%
  pairwise_wilcox_test(subst.elafin20.x ~ Genotype, p.adjust.method = "fdr") %>%
  add_significance() %>% 
  add_xy_position(x = "Genotype", scales="free")

# Look
stats_SNP_elafin_early
table(stats_SNP_elafin_early$p.adj.signif) #an idea of how many are sig
str(stats_SNP_elafin_early)
tapply(stats_SNP_elafin_early$p.adj, stats_SNP_elafin_early$p.adj.signif, summary)


# Checks on p-values 
stats_SNP_elafin_early
#   SNP        .y.              group1 group2    n1    n2 statistic     p p.adj p.adj.signif y.position groups        xmin  xmax
#   <fct>      <chr>            <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>             <dbl> <named list> <dbl> <dbl>
# 1 rs17333103 subst.elafin20.x CC     CT       194    60     5114  0.156 0.234 ns             1251102. <chr [2]>        1     2
# 2 rs17333103 subst.elafin20.x CC     TT       194     7      388  0.055 0.164 ns             1534666. <chr [2]>        1     3
# 3 rs17333103 subst.elafin20.x CT     TT        60     7      162  0.33  0.33  ns             1818230. <chr [2]>        2     3
# 4 rs1983649  subst.elafin20.x TT     AT        55   114     3859  0.015 0.023 *              1062059. <chr [2]>        3     5
# 5 rs1983649  subst.elafin20.x TT     AA        55    93     3339  0.002 0.006 **             1062059. <chr [2]>        3     6
# 6 rs1983649  subst.elafin20.x AT     AA       114    93     5801  0.244 0.244 ns             1062059. <chr [2]>        5     6
# 7 rs6032040  subst.elafin20.x TT     AT       178    82     6381  0.104 0.312 ns             1062059. <chr [2]>        3     5
# 8 rs6032040  subst.elafin20.x TT     AA       178     4      258  0.349 0.523 ns             1062059. <chr [2]>        3     6
# 9 rs6032040  subst.elafin20.x AT     AA        82     4      136. 0.58  0.58  ns             1062059. <chr [2]>        5     6
# 10 rs2664581  subst.elafin20.x CC     AC         7    65      288  0.254 0.29  ns             1251102. <chr [2]>        1     4
# 11 rs2664581  subst.elafin20.x CC     AA         7   189      932  0.067 0.201 ns             1534666. <chr [2]>        1     6
# 12 rs2664581  subst.elafin20.x AC     AA        65   189     6684  0.29  0.29  ns             1818230. <chr [2]>        4     6
head(snp_early)
ggpubr::compare_means(subst.elafin20.x ~ rs17333103, data=snp_early, method = "wilcox.test", p.adjust.method = "fdr") # matches 
ggpubr::compare_means(subst.elafin20.x ~ rs1983649, data=snp_early, method = "wilcox.test", p.adjust.method = "fdr") # matches 





# Plot of elafin by Genotype 
SNP_elafin_concentration_early <- ggplot(snp_early_gather, aes(x=Genotype, y=subst.elafin20.x)) + 
  geom_violin(aes(fill=Genotype)) + 
  facet_grid(~SNP, scales = "free_y") + #scales = "free"  to not plot empty factor levels
  labs(y="Elafin concentration (pg/\u00B5L)") + 
  stat_pvalue_manual(stats_SNP_elafin_early, 
                     step.increase=0.15,# step.group.by=SNP,
                     hide.ns = T, label = "{p.adj.signif}", #fdr
                     size=4) + 
  theme(legend.position = "none", 
        axis.text.x=element_text(size=7),
        strip.text.x = element_text(size = 11)) +
  scale_y_continuous(limits=c(0, 1400000), expand = expansion(mult = c(0, 0.10))) # Add 10% spaces between the top of the y and the plot border (& 0% space to bottom of y)
ggsave("plots/elafin_early/SNP_elafin_concentration_early.pdf", SNP_elafin_concentration_early, height=3.5, width=8)


# Ns for plot - these are only including snp data which has elafin conc
dim(snp_early) #264
table(snp_early$rs17333103, useNA="always")
table(snp_early$rs1983649, useNA="always")
table(snp_early$rs6032040 , useNA="always")
table(snp_early$rs2664581, useNA="always")



















# Now exactly the same code but early -> late


############################################################################################
##                                                                                        ##
##                                      df prep                                           ##
##                                                                                        ##
############################################################################################

# Cols we want 
snp_late <- data_late %>% 
  select(participant_id, rs17333103, rs1983649, rs2664581, rs6032040,
         subst.elafin20.x, weeks.at.visit.group, 
         Ethnicity_short, Ethnicity_long, sPTB_34w, sPTB_37w,
         BV_categories) 

#check
str(snp_late)

##get rid of samples where we have NAs for ALL SNPs 
dim(snp_late) #302   
snp_late <- snp_late[!with(snp_late, is.na(rs17333103) & is.na(rs1983649) & is.na(rs2664581) & is.na(rs6032040) ),]
dim(snp_late) #275     - dropped 27 samples

#drop samples without elafin
snp_late <- snp_late %>% 
  drop_na(subst.elafin20.x)
dim(snp_late) #263 


##data shaping prep
str(snp_late)
snp_late$rs17333103 <- as.character(snp_late$rs17333103)
snp_late$rs1983649 <- as.character(snp_late$rs1983649)
snp_late$rs2664581 <- as.character(snp_late$rs2664581)
snp_late$rs6032040 <- as.character(snp_late$rs6032040)
str(snp_late)


# Gather df for plotting
?gather
?pivot_longer
snp_late_gather <- pivot_longer(snp_late,
                                cols=rs17333103:rs6032040,
                                names_to="SNP",
                                values_to="Genotype")
head(snp_late_gather)


# Remove Genotype NAs in snp_late_gather to avoid them getting some of the %/proportion out of the SNP genotypes
table(snp_late_gather$Genotype, useNA="always") #8 NAs
snp_late_gather <- snp_late_gather[!is.na(snp_late_gather$Genotype),]
table(snp_late_gather$Genotype, useNA="always") #0 NAs



#Changing chr to factors
str(snp_late_gather)
snp_late_gather$Genotype <- factor(snp_late_gather$Genotype, levels=c("CC", "CT", "TT", "AC", "AT", "AA")) # these factor levels give the same outcome when doing the stat test as leaving it as character type (was wondering if excess levels would mess up p-values)
snp_late_gather$SNP <- factor(snp_late_gather$SNP, levels=c("rs17333103", "rs1983649", "rs6032040", "rs2664581"))
str(snp_late_gather)
summary(snp_late_gather)




############################################################################################
##                                                                                        ##
##                                Genotypes vs elafin                                     ##
##                                                                                        ##
############################################################################################

# Looking
head(snp_late_gather)
?rstatix::pairwise_wilcox_test
?rstatix::adjust_pvalue
?rstatix::add_significance
?rstatix::add_xy_position
?stat_pvalue_manual
?facet_grid

# df of stats for plot - Genotype fdr
stats_SNP_elafin_late <- snp_late_gather %>%
  group_by(SNP) %>%
  pairwise_wilcox_test(subst.elafin20.x ~ Genotype, p.adjust.method = "fdr") %>%
  add_significance() %>% 
  add_xy_position(x = "Genotype", scales="free")

# Look
stats_SNP_elafin_late
#   SNP        .y.              group1 group2    n1    n2 statistic     p p.adj p.adj.signif y.position groups        xmin  xmax
#   <fct>      <chr>            <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>             <dbl> <named list> <dbl> <dbl>
# 1 rs17333103 subst.elafin20.x CC     CT       191    62     5032. 0.076 0.076 ns             1143117. <chr [2]>        1     2
# 2 rs17333103 subst.elafin20.x CC     TT       191     7      277  0.009 0.026 *              1264704. <chr [2]>        1     3
# 3 rs17333103 subst.elafin20.x CT     TT        62     7      120  0.055 0.076 ns             1386291. <chr [2]>        2     3
# 4 rs1983649  subst.elafin20.x TT     AT        57   115     3811  0.083 0.113 ns             1159115. <chr [2]>        3     5
# 5 rs1983649  subst.elafin20.x TT     AA        57    89     3280  0.003 0.009 **             1304699. <chr [2]>        3     6
# 6 rs1983649  subst.elafin20.x AT     AA       115    89     5780  0.113 0.113 ns             1450283. <chr [2]>        5     6
# 7 rs6032040  subst.elafin20.x TT     AT       175    84     7898. 0.332 0.895 ns             1146038. <chr [2]>        3     5
# 8 rs6032040  subst.elafin20.x TT     AA       175     4      336  0.895 0.895 ns             1272007. <chr [2]>        3     6
# 9 rs6032040  subst.elafin20.x AT     AA        84     4      149  0.711 0.895 ns             1397976. <chr [2]>        5     6
# 10 rs2664581  subst.elafin20.x CC     AC         7    67      335  0.065 0.065 ns             1143117. <chr [2]>        1     4
# 11 rs2664581  subst.elafin20.x CC     AA         7   186     1037  0.008 0.023 *              1264704. <chr [2]>        1     6
# 12 rs2664581  subst.elafin20.x AC     AA        67   186     7352. 0.029 0.043 *              1386291. <chr [2]>        4     6
table(stats_SNP_elafin_late$p.adj.signif) #an idea of how many are sig
str(stats_SNP_elafin_late)
tapply(stats_SNP_elafin_late$p.adj, stats_SNP_elafin_late$p.adj.signif, summary)


# Checks on p-values - adjusted p-values are v. similar but not identical
stats_SNP_elafin_late
# head(snp_late)
# ggpubr::compare_means(subst.elafin20.x ~ rs17333103, data=snp_late, method = "wilcox.test", p.adjust.method = "fdr") # matches 
# ggpubr::compare_means(subst.elafin20.x ~ rs1983649, data=snp_late, method = "wilcox.test", p.adjust.method = "fdr") 





# Plot of elafin by Genotype 
SNP_elafin_concentration_late <- ggplot(snp_late_gather, aes(x=Genotype, y=subst.elafin20.x)) + 
  geom_violin(aes(fill=Genotype)) + 
  facet_grid(~SNP, scales = "free_y") + #scales = "free"  to not plot empty factor levels
  labs(y="Elafin concentration (pg/\u00B5L)") + 
  stat_pvalue_manual(stats_SNP_elafin_late, 
                     # step.increase=0.005,# step.group.by=SNP,
                     hide.ns = T, label = "{p.adj.signif}", #fdr
                     size=4) + 
  theme(legend.position = "none", 
        axis.text.x=element_text(size=7),
        strip.text.x = element_text(size = 11)) +
  scale_y_continuous(limits=c(0, 1400000), expand = expansion(mult = c(0, 0.10))) # Add 10% spaces between the top of the y and the plot border (& 0% space to bottom of y)
ggsave("plots/elafin_late/SNP_elafin_concentration_late.pdf", SNP_elafin_concentration_late, height=3.5, width=8)
?scale_y_continuous



# Ns for plot - these are only including snp data which has elafin conc
dim(snp_late) #263
table(snp_late$rs17333103, useNA="always")
table(snp_late$rs1983649, useNA="always")
table(snp_late$rs6032040 , useNA="always")
table(snp_late$rs2664581, useNA="always")

