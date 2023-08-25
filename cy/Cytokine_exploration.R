# Clear environment
rm(list=ls())

# Load libraries
library(car)
library(multcomp)
library(reshape2)
library(tidyverse)
library(stats) 
library(tidyr)
library(BiocManager) #install.packages("BiocManager")
library(devtools)
library(ggbiplot) #install_github("vqv/ggbiplot")
library(ggrepel)
library(sva) #BiocManager::install("sva")
library(data.table) #filtering dfs
library(corrplot)#install.packages('corrplot')
library(ggpubr) #adj ps on ggplot 
library(rstatix) #adj ps on ggplot

#Versions of packages used
sessionInfo()


## Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects") #mac
setwd("C:/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines") #pc
getwd()
#list.files()







###########################################################################################################################
##                                                                                                                       ##
##                                             Load in data                                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Loading cytokine .Rdata created in the script "Immune__data_wrangling_cleaning_and_merging_with_metadata__CytokCompIgM.R"
lnames = load(file="data/immune_CytokinesComplementsIgM_2023.RData")
lnames





###########################################################################################################################
##                                                                                                                       ##
##                                        For checking N in violins/boxplots                                             ##
##                                                                                                                       ##
###########################################################################################################################
?base::merge
names(immu_MD)
MD_for_Ns <- immu_MD %>% dplyr::select(Participant.ID:Risk)
check_Ns <- base::merge(immu_small, MD_for_Ns, all.x=T, all.y=F)
# df of metadata where we have at least 1 immune marker




###########################################################################################################################
##                                                                                                                       ##
##                                             Heatmap of cytokines                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Info
?corrplot
citation("corrplot")
?psych::corr.test
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

# Use corr.test() from psych package to calculate the correlation matrix and corresponding p value matrix without adjustment.
ct_cytok <- psych::corr.test(immu_small[,-1], adjust = "fdr") #immu_small[,-1] is the df w/out the Participant.ID col
corr_cytok <- ct_cytok$r
p.mat_cytok <- ct_cytok$p
p.mat_cytok

# Finding n for captions
ct_cytok$n
min(ct_cytok$n)
max(ct_cytok$n)

####### Finalising FDR p-value matrix
# Duplicating p-value matrix - later this will be of only FDR p-values
p_adj_cytok <- p.mat_cytok

# Make nominal p-values (lower.tri) into NAs
p_adj_cytok[lower.tri(p_adj_cytok)] <- NA # don't need this step but I include to verify
p_adj_cytok 

# Make nominal p-values (lower.tri) into a symmetrical version of the FDR adjusted p-values 
p_adj_cytok[lower.tri(p_adj_cytok)] = t(p_adj_cytok)[lower.tri(p_adj_cytok)]
p_adj_cytok 

# Checking symmetry 
isSymmetric(p_adj_cytok) #TRUE
#######

# Plot corr_cytok with hclust - square, no p-values
dev.off()
pdf("plots/immune_basics/heatmap_immune_all.pdf", width=5, height=5)
corrplot(corr_cytok, #type = "upper", 
         order = "hclust", #diag = F,
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400)) 
dev.off()



# Plot corr_cytok with hclust - leave blank on non-significant coefficient & add text of significant correlation coefficients
pdf("plots/immune_basics/heatmap_immune_FDR_sig.pdf", width=5, height=5)
corrplot(corr_cytok, type = "upper", 
         order = "hclust", diag = F,
         p.mat = p_adj_cytok, sig.level = 0.05, insig ="blank",
         method = "square", tl.cex = 1, tl.col = 'black', #tl = axes text
         col= colorRampPalette(c("blue", "white", "red"))(400),
         addCoef.col = "grey10", number.cex = 0.7, number.digits = 2, number.font = 1) #number.font = 1 makes cor coefficients regular, not bold
dev.off()






###########################################################################################################################
##                                                                                                                       ##
##                                            Boxplot w stats prep                                                       ##
##                                                                                                                       ##
###########################################################################################################################

# Looking
?pivot_longer
names(immu_MD)
table(immu_MD$Ethnicity)
#View(immu_MD)

## Adding adjusted p-values to plots 
?rstatix::pairwise_wilcox_test
?rstatix::adjust_pvalue
?rstatix::add_significance
?rstatix::add_xy_position

# Color schemes for ggplots
Ethnicity3colors <- scale_color_manual(values = c("blue", "#42D5F9", "#C27BFA")) 
TB_PTB2colors <- scale_color_manual(values = c("#83EC88", "tomato"), na.value = "grey60") # no yes
Risk2colors <- scale_color_manual(values = c("#F99242", "#FFD44C"), na.value = "grey60")  # high low

# Color fill schemes for ggplots
Ethnicity3fills <- scale_fill_manual(values = c("blue", "#42D5F9", "#C27BFA")) 
TB_PTB2fills <- scale_fill_manual(values = c("#83EC88", "tomato"), na.value = "grey60") # no yes
Risk2fills <- scale_fill_manual(values = c("#47FFD8", "#FFD44C"), na.value = "grey60")  # high low

# Facet labels
Immune_marker_labs <- c("C5", "C5a", "G-CSF", "IFN-g", "IgM", "IL-10", "IL-13", "IL-17A", "IL-22", "IL-5", "IL-6", "IL-8", "MCP-1") #gamma symbol = \u03b3
names(Immune_marker_labs) <- c("C5", "C5a", "GCSF", "IFNg", "IgM", "IL10", "IL13", "IL17A", "IL22", "IL5", "IL6", "IL8", "MCP1")






###########################################################################################################################
##                                                                                                                       ##
##                                              Immune vs clinical metadata                                              ##
##                                                                                                                       ##
###########################################################################################################################

# Reorder risk factor
table(immu_MD$Risk)
immu_MD$Risk <- factor(immu_MD$Risk, levels=c("Low risk", "High risk"))
table(immu_MD$Risk)

# Gather df of immune markers
immune_gathered <- as.data.frame(pivot_longer(immu_MD, cols=C5:MCP1, names_to="Immune_marker", values_to="Quantity" ))
dim(immune_gathered) #[1] 1131  201
immune_gathered <- immune_gathered %>% dplyr::select(Participant.ID, Immune_marker, Quantity, Ethnicity, sPTB34, sPTB37, everything())
head(immune_gathered)
str(immune_gathered$Immune_marker)
str(immune_gathered$Quantity)


###################################### Ethnicity
# df of stats for plot - Ethnicity 
stats_Immune_marker_Quantity_Ethnicity <- immune_gathered %>%
  group_by(Immune_marker) %>%
  pairwise_wilcox_test(Quantity ~ Ethnicity) %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>% 
  add_xy_position(x = "Ethnicity", scales="free")
stats_Immune_marker_Quantity_Ethnicity
table(stats_Immune_marker_Quantity_Ethnicity$p.adj.signif) #1 sig
table(check_Ns$Ethnicity)

# Boxplot of immune markers by Ethnicity 
Immune_markers_by_Ethnicity <- ggplot(immune_gathered, aes(x=Ethnicity, y=Quantity))  + 
  geom_boxplot(outlier.shape=NA) + #outlier.shape=NA make no black outlier dots, so we just have 1 dot per ID cos of the jitter layer
  geom_jitter(aes(colour=Ethnicity), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "Ethnicity", color="Ethnicity") +
  Ethnicity3colors +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_Ethnicity, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% spaces between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_Ethnicity.pdf", Immune_markers_by_Ethnicity, height=5.5, width=10)

# Looking
tapply(immu_MD$IL5, immu_MD$Ethnicity, summary) #African-Caribbean women had higher IL-5 concentrations than African women

# Violin of immune markers by Ethnicity 
Immune_markers_by_Ethnicity_violin <- ggplot(immune_gathered, aes(x=Ethnicity, y=Quantity))  + 
  geom_violin(aes(fill=Ethnicity)) + #outlier.shape=NA make no black outlier dots, so we just have 1 dot per ID cos of the jitter layer
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "Ethnicity", color="Ethnicity") +
  Ethnicity3fills +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_Ethnicity, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% spaces between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_Ethnicity_violin.pdf", Immune_markers_by_Ethnicity_violin, height=5.5, width=10)

######################################



######################################sPTB34
# df of stats for plot - sPTB34 
stats_Immune_marker_Quantity_sPTB34 <- immune_gathered %>%
  group_by(Immune_marker) %>%
  pairwise_wilcox_test(Quantity ~ sPTB34) %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>% 
  add_xy_position(x = "sPTB34", scales="free")
stats_Immune_marker_Quantity_sPTB34
table(stats_Immune_marker_Quantity_sPTB34$p.adj.signif) #1 sig
table(check_Ns$sPTB34)

# Boxplot of immune markers by sPTB34 
Immune_markers_by_sPTB34 <- ggplot(immune_gathered, aes(x=sPTB34, y=Quantity))  + #label=Participant.ID
  geom_boxplot(outlier.shape=NA) + #outlier.shape=NA make no black outlier dots, so we just have 1 dot per ID cos of the jitter layer
  geom_jitter(aes(colour=sPTB34), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  #geom_text() +
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "sPTB <34 weeks' gestation", color="sPTB <34w") +
  TB_PTB2colors +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_sPTB34, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_sPTB34.pdf", Immune_markers_by_sPTB34, height=5.5, width=10)

# Violin of immune markers by sPTB34 
Immune_markers_by_sPTB34_violin <- ggplot(immune_gathered, aes(x=sPTB34, y=Quantity))  + #label=Participant.ID
  geom_violin(aes(fill=sPTB34)) + 
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "sPTB <34 weeks' gestation", color="sPTB <34w") +
  TB_PTB2fills +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_sPTB34, #bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_sPTB34_violin.pdf", Immune_markers_by_sPTB34_violin, height=5.5, width=10)

######################################


######################################sPTB37
# df of stats for plot - sPTB37 
stats_Immune_marker_Quantity_sPTB37 <- immune_gathered %>%
  group_by(Immune_marker) %>%
  pairwise_wilcox_test(Quantity ~ sPTB37) %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>% 
  add_xy_position(x = "sPTB37", scales="free")
stats_Immune_marker_Quantity_sPTB37
table(stats_Immune_marker_Quantity_sPTB37$p.adj.signif) #0 sig
table(check_Ns$sPTB37)

# Boxplot of immune markers by sPTB37 
Immune_markers_by_sPTB37 <- ggplot(immune_gathered, aes(x=sPTB37, y=Quantity))  + 
  geom_boxplot(outlier.shape=NA) + #outlier.shape=NA make no black outlier dots, so we just have 1 dot per ID cos of the jitter layer
  geom_jitter(aes(colour=sPTB37), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "sPTB <37 weeks' gestation", color="sPTB <37w") +
  TB_PTB2colors +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_sPTB37, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_sPTB37.pdf", Immune_markers_by_sPTB37, height=5.5, width=10)

# Violin of immune markers by sPTB37 
Immune_markers_by_sPTB37_violin <- ggplot(immune_gathered, aes(x=sPTB37, y=Quantity))  + #label=Participant.ID
  geom_violin(aes(fill=sPTB37)) + 
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "sPTB <37 weeks' gestation", color="sPTB <37w") +
  TB_PTB2fills +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_sPTB37, #bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_sPTB37_violin.pdf", Immune_markers_by_sPTB37_violin, height=5.5, width=10)

######################################


######################################Risk
# df of stats for plot - Risk 
stats_Immune_marker_Quantity_Risk <- immune_gathered %>%
  group_by(Immune_marker) %>%
  pairwise_wilcox_test(Quantity ~ Risk) %>%
  adjust_pvalue(method = "none") %>%
  add_significance() %>% 
  add_xy_position(x = "Risk", scales="free")
stats_Immune_marker_Quantity_Risk
table(stats_Immune_marker_Quantity_Risk$p.adj.signif) #0 sig
table(check_Ns$Risk)

# Boxplot of immune markers by Risk 
Immune_markers_by_Risk <- ggplot(immune_gathered, aes(x=Risk, y=Quantity))  + 
  geom_boxplot(outlier.shape=NA) + #outlier.shape=NA make no black outlier dots, so we just have 1 dot per ID cos of the jitter layer
  geom_jitter(aes(colour=Risk), alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "Risk category", color="Risk") +
  Risk2colors +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_Risk, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_Risk.pdf", Immune_markers_by_Risk, height=5.5, width=10)

# violin of immune markers by Risk 
Immune_markers_by_Risk_violin <- ggplot(immune_gathered, aes(x=Risk, y=Quantity))  + 
  geom_violin(aes(fill=Risk)) +
  facet_wrap(~Immune_marker, labeller=labeller(Immune_marker=Immune_marker_labs), strip.position="top", nrow=3, scales="free") + 
  labs(y="Concentration in cervicovaginal fluid", x = "Risk category", color="Risk") +
  Risk2fills +
  theme(axis.text.x=element_text(size=8), strip.text.x=element_text(angle = 0, size=12),
        legend.position = "none") +  
  stat_pvalue_manual(stats_Immune_marker_Quantity_Risk, bracket.nudge.y = -2, #Move down the bracket using `bracket.nudge.y`
                     hide.ns = T, label = "{p.adj}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Add 15% space between the p-value labels and the plot border
ggsave("plots/immune_basics/Immune_markers_by_Risk_violin.pdf", Immune_markers_by_Risk_violin, height=5.5, width=10)

######################################









###########################################################################################################################
##                                                                                                                       ##
##                                                     Immune vs pH                                                      ##
##                                                                                                                       ##
###########################################################################################################################

# Look
#View(immu_MD)
names(immu_MD)
immu_MD$medscinet_ph
tapply(immu_MD$medscinet_ph, immu_MD$Ethnicity, summary)

#Using Flavia's pH values, as some values in medscinet col are invalid


# How many women do we have their pH values
Immune_IDs <- immu_small$Participant.ID
length(Immune_IDs) #81
dim(immu_MD) #87 women
pH_n <- immu_MD %>% dplyr::filter(Participant.ID  %in% Immune_IDs & medscinet_ph != "NA") %>% 
  dplyr::select(Participant.ID, medscinet_ph, everything())
dim(pH_n) #65 women



################################################# Linear regression 
?lm 
names(immu_MD)
names(immu_small)

# C5
lm_C5 <- lm(C5 ~ medscinet_ph, data=immu_MD)
summary(lm_C5) #0.000732 

# C5a
lm_C5a <- lm(C5a ~ medscinet_ph, data=immu_MD)
summary(lm_C5a) #ns

# GCSF
lm_GCSF <- lm(GCSF ~ medscinet_ph, data=immu_MD)
summary(lm_GCSF) #ns

# IFNg
lm_IFNg <- lm(IFNg ~ medscinet_ph, data=immu_MD)
summary(lm_IFNg) #ns

# IgM
lm_IgM <- lm(IgM ~ medscinet_ph, data=immu_MD)
summary(lm_IgM) #ns

# IL10
lm_IL10 <- lm(IL10 ~ medscinet_ph, data=immu_MD)
summary(lm_IL10) #0.031

# IL13
lm_IL13 <- lm(IL13 ~ medscinet_ph, data=immu_MD)
summary(lm_IL13) #0.00462

# IL17A
lm_IL17A <- lm(IL17A ~ medscinet_ph, data=immu_MD)
summary(lm_IL17A) #0.0351

# IL22
lm_IL22 <- lm(IL22 ~ medscinet_ph, data=immu_MD)
summary(lm_IL22) #0.0046

# IL5
lm_IL5 <- lm(IL5 ~ medscinet_ph, data=immu_MD)
summary(lm_IL5) #ns

# IL6
lm_IL6 <- lm(IL6 ~ medscinet_ph, data=immu_MD)
summary(lm_IL6) #ns

# IL8
lm_IL8 <- lm(IL8 ~ medscinet_ph, data=immu_MD)
summary(lm_IL8) #ns

# MCP1
lm_MCP1 <- lm(MCP1 ~ medscinet_ph, data=immu_MD)
summary(lm_MCP1) #0.00353
#################################################



################################################# Plotting immune vs pH
#Creating facet labels 
vectorImmune
pH_labs <- c("C5 *", "C5a", "G-CSF", "IFN-g", "IgM", "IL-10 *", "IL-13 *", "IL-17A *", "IL-22 *", "IL-5", "IL-6", "IL-8", "MCP-1 *")
names(pH_labs) <- vectorImmune
pH_labs


# Scatter of immune by medscinet_ph - with linear regression line
Immune_markers_vs_pH <- ggplot(base::subset(immune_gathered, medscinet_ph != "NA"), aes(x=medscinet_ph, y=Quantity)) + 
  geom_point(alpha=0.3) +
  facet_wrap(~Immune_marker, strip.position="top", nrow=3, scales="free", labeller=labeller(Immune_marker=pH_labs)) +
  labs(y="Concentration in cervicovaginal fluid", x = "Cervicovaginal pH", color="Immune_marker") +
  theme(legend.position = "none", strip.text.x=element_text(angle = 0, size=10)) +
  scale_x_continuous(breaks=seq(3,6,1)) +
  geom_smooth(method='lm') # for some reason I cannot get it to be black instead of blue using aes(color=c(replicate(7, "black")))
ggsave("plots/immune_basics/Immune_markers_vs_pH.pdf", Immune_markers_vs_pH , height=5.5, width=10)
#################################################








###########################################################################################################################
##                                                                                                                       ##
##                                                   Immune vs BMI                                                       ##
##                                                                                                                       ##
###########################################################################################################################

# Looking more
dim(immu_MD) # 87 women
summary(immu_MD$BMI) #1 NA

# How many women do we have their BMI values
dim(immu_MD) #87 women
BMI_n <- immu_MD %>% dplyr::filter(Participant.ID  %in% Immune_IDs & BMI != "NA") %>% 
  dplyr::select(Participant.ID, BMI, everything())
dim(BMI_n) #80 women



################################################# Linear regression
?lm 
names(immu_MD)
names(immu_small)

# C5
lm_C5 <- lm(C5 ~ BMI, data=immu_MD)
summary(lm_C5) #0.000526

# Testing if C5 is still sig if I remove the outlier BMI ID
sort(immu_MD$BMI)
lm_C5 <- lm(C5 ~ BMI, data=base::subset(immu_MD, BMI<47))
summary(lm_C5) #ns 

# C5a
lm_C5a <- lm(C5a ~ BMI, data=immu_MD)
summary(lm_C5a) #ns

# GCSF
lm_GCSF <- lm(GCSF ~ BMI, data=immu_MD)
summary(lm_GCSF) #ns

# IFNg
lm_IFNg <- lm(IFNg ~ BMI, data=immu_MD)
summary(lm_IFNg) #ns

# IgM
lm_IgM <- lm(IgM ~ BMI, data=immu_MD)
summary(lm_IgM) #ns

# IL10
lm_IL10 <- lm(IL10 ~ BMI, data=immu_MD)
summary(lm_IL10) #ns

# IL13
lm_IL13 <- lm(IL13 ~ BMI, data=immu_MD)
summary(lm_IL13) #ns

# IL17A
lm_IL17A <- lm(IL17A ~ BMI, data=immu_MD)
summary(lm_IL17A) #ns

# IL22
lm_IL22 <- lm(IL22 ~ BMI, data=immu_MD)
summary(lm_IL22) #ns

# IL5
lm_IL5 <- lm(IL5 ~ BMI, data=immu_MD)
summary(lm_IL5) #ns

# IL6
lm_IL6 <- lm(IL6 ~ BMI, data=immu_MD)
summary(lm_IL6) #ns

# IL8
lm_IL8 <- lm(IL8 ~ BMI, data=immu_MD)
summary(lm_IL8) #ns

# MCP1
lm_MCP1 <- lm(MCP1 ~ BMI, data=immu_MD)
summary(lm_MCP1) #ns
#################################################


#Creating facet labels 
vectorImmune
BMI_labs <- c("C5", "C5a", "G-CSF", "IFN-g", "IgM", "IL-10", "IL-13", "IL-17A", "IL-22", "IL-5", "IL-6", "IL-8", "MCP-1")
names(BMI_labs) <- vectorImmune
BMI_labs

# Scatter of immune by BMI - with linear regression line
Immune_markers_vs_BMI <- ggplot(base::subset(immune_gathered, BMI != "NA"), aes(x=BMI, y=Quantity)) + 
  geom_point(alpha=0.4)  + 
  facet_wrap(~Immune_marker, strip.position="top", nrow=3, scales="free", labeller=labeller(Immune_marker=BMI_labs)) +
  labs(y="Concentration in cervicovaginal fluid", x = bquote("Body mass index "(kg/m^2)), color="Immune_marker") + 
  theme(legend.position = "none", strip.text.x=element_text(angle = 0, size=10)) +
  geom_smooth(method='lm') # for some reason I cannot get it to be black instead of blue using aes(color=c(replicate(7, "black")))
ggsave("plots/immune_basics/Immune_markers_vs_BMI.pdf", Immune_markers_vs_BMI , height=5, width=8)










###########################################################################################################################
##                                                                                                                       ##
##                                                     Save data                                                         ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data for easier re-running as this script takes a long time to run
save.image(file = "data/Cytokine_exploration_script_ALL_DATA.RData", compress=F) #resaved on 28/03/2023 to add extra violins and check_Ns df

# # Loading .Rdata saved
# lnames = load(file="data/Cytokine_exploration_script_ALL_DATA.RData")
# lnames



