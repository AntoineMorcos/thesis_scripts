
# Clear environment
rm(list=ls())

# Load libraries
library(reshape2) # heatmaps
library(BiocManager) #install.packages("BiocManager")
library(tidyverse) # ggplot2, dplyr, tidyr etc
library(ggrepel)
library(stats) # for Wilcoxon or other stats tests
library(ggpubr) # for stats test in the plots
library(lubridate) # for working with dates

library(car)
library(multcomp)
library(reshape2)
library(tidyverse)
library(stats) 
library(dplyr)
library(tidyr)

library(BiocManager) #install.packages("BiocManager")
library(DESeq2) #BiocManager::install("DESeq2")
library(ggrepel)
library(EnhancedVolcano) #BiocManager::install('EnhancedVolcano') #Nice volcano plots



## Set wd
# setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files("data/larger_cohort")







# Loading .Rdata saved
lnames = load(file="data/ALL_objects_in_Immune_cell_types_in_larger_cohort.RData")
lnames

# Loading .Rdata saved
lnames = load(file="data/ALL_objects_in_Plotting_MAs_volcanoes_heatmaps_and_PCAs_of_DESeq_data.RData")
lnames

# Loading .Rdata saved
lnames = load(file="data/ALL_objects_in_GO_and_GAGE_Pathview_Pathway_analysis_of_ethnicity_DEGs.RData")
lnames

# Loading .Rdata saved
lnames = load(file="data/ALL_objects_in_Metagenome_VIRGO_basics.RData")
lnames

# Loading .Rdata saved
lnames = load(file="data/ALL_objects_in_WGCNA_plotting_selected_GO_VIRGO.RData")
lnames

# Loading .Rdata saved
lnames = load(file="data/WGCNA/ALL_WGCNA_objects.RData")
lnames



dev.off()







###########################################################################################################################
##                                                                                                                       ##
##                                                 Figure 1 - done                                                       ##
##                                                                                                                       ##
###########################################################################################################################

# Figure 1A - flow cyto
# Make image in Adobe Illustrator (180w 33h)

# Figure 1B - cell proportions stacked bar
Figure_1B <- ggplot(summary_immune_AllTogether, aes(x=nominal, fill=Cell_type_pretty, y=Mean_percentage)) + 
  geom_bar(stat="identity", position="stack") + 
  labs(x=NULL, y="Mean percentage of CD45+ cells (%)", fill="CD45+ cells") + 
  fills_of_neutrophil +
  ylim(c(0,100)) +
  geom_text(aes(x=nominal, label=ifelse(Mean_percentage>4, paste0(round(Mean_percentage, 0), "%"),"")), #means label if it is above x% otherwise do not label
            colour = "white",  position=position_stack(vjust=0.5), size=2) + 
  theme(axis.ticks = element_blank(),
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
        legend.text=element_text(size=8), legend.title = element_text(size = 9), #legend text & title
        legend.key.size = unit(5, 'mm') #size of legend icons
  )
ggsave("plots/PAPER_FIGURES/Figure_1B.tiff", Figure_1B, height=55, width=85, units="mm", dpi=300)


# Figure 1C - PTB/TB cell proportions
Labs_1C <- c("Neutrophils", "NK cells", "B cells", "T cells", "Patrolling\nmonocytes", "Intermediate\nmonocytes", "Classical\nmonocytes")
names(Labs_1C) <- c("Neutrophils", "NKcells", "Bcells", "Tcells", "Patrolling_monocytes", "Intermediate_monocytes", "Classical_monocytes")

Figure_1C <- ggplot(immune_gather_uniq, aes(x=Outcome, y=Percentage))  + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(aes(colour=Outcome), size=1, alpha=0.6, width=0.25, height=0) + # Don't change height when jittering!
  facet_wrap(~Cell_type, strip.position="top", nrow=1, scale="free", labeller=labeller(Cell_type=Labs_1C)) +
  labs(y="Percentage of CD45+ cells (%)") +
  scale_color_manual(values = c("#00BFC4", "tomato", "grey50")) +
  theme(legend.position="none", 
        text = element_text(size = 8), 
        strip.text.x = element_text(size = 8), #facet labels text
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 8) #axes numbers text size
  )  
ggsave("plots/PAPER_FIGURES/Figure_1C.tiff", Figure_1C, height=55, width=180, units="mm", dpi=300)
ggsave("plots/PAPER_FIGURES/Figure_1C.pdf", Figure_1C, height=55, width=180, units="mm")


# Figure 1D - gestation lm cell proportions
Labs_1D <- c("Neutrophils *", "NK cells", "B cells", "T cells *", "Patrolling\nmonocytes", "Intermediate\nmonocytes", "Classical\nmonocytes")
names(Labs_1D) <- c("Neutrophils", "NKcells", "Bcells", "Tcells", "Patrolling_monocytes", "Intermediate_monocytes", "Classical_monocytes")

Figure_1D <- ggplot(immune_gather, aes(x=Gest.at.cyto.wks.dec, y=Percentage)) + 
  geom_point(aes(color=Cell_type), alpha=0.6, size=1)  + 
  facet_wrap(~Cell_type, strip.position="top", nrow=1, labeller=labeller(Cell_type=Labs_1D), scales="free") + #, scales="free"
  labs(x="Sampling time (weeks' gestation)", y="Percentage of CD45+ cells (%)", color="CD45+ cells") +
  colours_of_neutrophil +
  geom_smooth(method='lm') + # for some reason I cannot get it to be black instead of blue using aes(color=c(replicate(7, "black")))
  theme(legend.position="none", 
        text = element_text(size = 8), 
        strip.text.x = element_text(size = 8), #facet labels text
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8) #axes numbers text size
  )
ggsave("plots/PAPER_FIGURES/Figure_1D.tiff", Figure_1D, height=55, width=180, units="mm", dpi=300)
ggsave("plots/PAPER_FIGURES/Figure_1D.pdf", Figure_1D, height=55, width=180, units="mm")


# Combined heights (A4 has height 297mm)
33+(55*4)



###########################################################################################################################
##                                                                                                                       ##
##                                           Figure 2 - done but not ideal layout                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Figure_2A - PCA
Figure_2A <- ggplot(PCA_for_paper_data, aes(x=PC1, y=PC2, color=clusters_hc_k4, shape=Ethnicity)) +
  geom_point(size=2, stroke=1) +
  scale_color_manual(values = c("#FF226F", "#3D6A89", "#BEABFF", "#30B2FF")) +
  scale_shape_manual(values=c(3, 16), na.translate=F) +
  xlab(paste0("PC1: ", PCA_for_paper_percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", PCA_for_paper_percentVar[2], "% variance")) +
  guides(color=guide_legend(order=1), shape = guide_legend(order=2)) +
  labs(color="Cluster") +
  theme(text = element_text(size = 8), 
        strip.text.x = element_text(size = 8), #facet labels text
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size 
        legend.text=element_text(size=8), legend.title = element_text(size = 9), #legend text & title
        legend.key.size = unit(5, 'mm') #size of legend icons
        )
ggsave("plots/PAPER_FIGURES/Figure_2A.tiff", Figure_2A, height=70, width=85, units="mm", dpi=300)


# Figure_2B draft - cluster diagram  
dev.off()
pdf("plots/PAPER_FIGURES/drafts_DO_NOT_SUBMIT_THESE/Figure_2B_draft.pdf", height=6, width=10)
par(mfrow=c(1,1))
pheatmap(rld_corr_dds1, annotation_col=metadata_for_heatmap, border_color=NA,
         show_rownames=T, show_colnames=T, #Show IDs
         annotation_colors = my_colours,
         legend=T, legend_breaks=c(seq(0,1,0.05)),
         display_numbers=F #fontsize_number=6, #Prints the distance numbers over the cells
         )
dev.off()
# Final image has size 8 font everywhere apart from size 6 font for ID labels


# Figure_2C draft - volcano 
pdf("plots/PAPER_FIGURES/drafts_DO_NOT_SUBMIT_THESE/Figure_2C_draft.pdf", height=10, width=10)
EnhancedVolcano(res_sLFC_Ethnicity_sort, # final res from DESeq
                lab = rownames(res_sLFC_Ethnicity_sort),
                x = 'log2FoldChange',  # data to plot within res on the x-axis
                y = 'padj', # data to plot within res on the y-axis
                selectLab = c(top_10_up_and_down), #vector of the top 10 upregulated DEGs & top 10 downregulated DEGs
                labSize = 5, # font size of labels
                xlab = bquote(~Log[2]~ 'fold change'), # label x axis
                ylab = bquote("-" ~Log[10]~ '(FDR adjusted p-value)'), # label y axis
                pCutoff = 0.05, FCcutoff = 1.5, #lines
                pointSize = 2,
                colCustom = Colour_scheme_volcano, # the colour scheme we made above with our ifelse statements
                colAlpha = 0.6, # Alpha - transparency of points
                legendPosition = 'none', # no legend - makes it messy. Previously was "left"
                drawConnectors = TRUE, widthConnectors = 0.4, colConnectors = 'grey', arrowheads = F, # Connector line between gene text and point options
                #lengthConnectors = 5,
                gridlines.major = T, gridlines.minor = F,
                border = 'partial', borderWidth = 1, borderColour = 'black') # Border options
dev.off()
# Final image has size 8 font on axes, but size 5 font for top 10 DEG labels 


# Figure_2D - DEG GO
# Reverse log scales function - https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
# Creating facet labels 
GO.labs <- c("Biological process", "Molecular function", "Cellular component")
names(GO.labs) <- c("BP", "MF", "CC")
# Plotting point - up in red & down in blue
Figure_2D <- ggplot(GO_selected, aes(x=p.val, y=reorder(GO_term, -p.val)))  + #order GO_terms by p-values
  geom_point(aes(colour=Direction)) +
  facet_wrap(~ GO_type, scale="free", labeller=labeller(GO_type=GO.labs)) +
  labs(y=NULL, x="p-value") +
  scale_x_continuous(trans=reverselog_trans(10), limits=c(0.05, 0.03)) + #  xlim(0.05, 0.03) +breaks=c(0.05, 0.03), 
  theme(legend.position="none") +
  scale_colour_manual(values = c("red", "blue")) +
  theme(plot.title.position = "plot", text = element_text(size = 8), 
        strip.text.x = element_text(size = 7), #facet labels text
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
        legend.key.size = unit(5, 'mm') #size of legend icons
  )
ggsave("plots/PAPER_FIGURES/Figure_2D.pdf", Figure_2D, height=55, width=85, units="mm", dpi=300)



# Figure_2E - IPA network
# Double clicked plots/IPA/Immune_Network1_export.pdf to open it in AI. 
# C&P it to a 85mm(w)x80mm(h) artboard
# fonts<8


# Combined heights (A4 has height 297mm)
check2 <- 100 + # PCA (70mm) & dendro (100mm)
  80 + #volcano (80mm) & 2D DEG GO (55mm) 
  80 # & Network (80mm)
check2



###########################################################################################################################
##                                                                                                                       ##
##                                                Figure 3 - done                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# Figure_3 - WGCNA heatmap
# Made it by:
# Double clicked plots/WGCNA/VIRGO_in_datTraits/WGCNA_correlation_heatmap_p0.05_VIRGO_selectedMEs.pdf to open it in AI. 
# C&P it to a 180mmx180mm artboard
# Adjust font sizes (minimum of 8 apart from text in the heatmap coloured boxes (of cor&p-values))







###########################################################################################################################
##                                                                                                                       ##
##                                                Figure 4 - done                                                        ##
##                                                                                                                       ##
###########################################################################################################################

# The x-axis scales I want in all the below GO plots so that they are comparable
x_axis_for_GO_plots <- scale_x_continuous(trans=reverselog_trans(10), minor_breaks=lseq(1e-01,1e-18,18), limits=c(1e-01,1e-14))

theme_for_GO_plots <- theme(legend.position="none",
                            strip.text.y = element_text(angle = 0), #facet labels text
                            plot.title.position = "plot",
                            text = element_text(size = 8), 
                            axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
                            )

# Figure_4A - selected GO terms in magenta module [WAS FIGURE 3B]
Figure_4A <- ggplot(GO_magenta_selected, aes(x=Bonferroni, y=reorder(dataSetName, -Bonferroni)))  + #order dataSetNames by Bonferroni p-values
  geom_point(colour="magenta") +
  labs(y=NULL, x="Bonferroni adjusted p-value") +
  facet_grid(rows=vars(GO_type), space = "free_y", scales = "free") + 
  scale_x_continuous(trans=reverselog_trans(10), minor_breaks=lseq(1e-01,1e-18,18), limits=c(1e-01,1e-14)) + 
  theme(legend.position="none",
        strip.text.y = element_text(angle = 0), #facet labels text
        plot.title.position = "plot",
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
  )
ggsave("plots/PAPER_FIGURES/drafts_DO_NOT_SUBMIT_THESE/Figure_4A.pdf", Figure_4A, height=6.5, width=6)
# was: height=7, width=6


# Figure_4B - ALL sig p-values in ivory module [WAS FIGURE 3C]
Figure_4B <- ggplot(GO_ivory_all, aes(x=Bonferroni, y=reorder(dataSetName, -Bonferroni)))  + #order dataSetNames by Bonferroni p-values
  geom_point(colour="black") + #ivory points are too difficult to see
  labs(y=NULL, x="Bonferroni adjusted p-value")  +
  facet_grid(rows=vars(GO_type), space = "free_y", scales = "free") + 
  scale_x_continuous(trans=reverselog_trans(10), minor_breaks=lseq(1e-01,1e-18,18), limits=c(1e-01,1e-14)) +  
  theme(legend.position="none",
        strip.text.y = element_text(angle = 0), #facet labels text
        plot.title.position = "plot",
        text = element_text(size = 8), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
  )
ggsave("plots/PAPER_FIGURES/drafts_DO_NOT_SUBMIT_THESE/Figure_4B.pdf", Figure_4B, height=1.2, width=5.5) 
# was: height=1.7, width=5.5







###########################################################################################################################
##                                                                                                                       ##
##                                              Supplementary Figure 1                                                   ##
##                                                                                                                       ##
###########################################################################################################################

# Supplementary_Figure_1A - alpha diversity
Supplementary_Figure_1A <- ggplot(Shannon_gathered, aes(x=Ethnicity, y=diversity_score))  + 
  geom_boxplot() +
  scale_y_continuous(breaks=seq(0, 4, 1), minor_breaks=seq(0, 4, 0.5), limits=c(0,3.5)) +
  geom_jitter(aes(colour=Ethnicity), alpha=0.6, width=0.3, height=0) + #hollow circles and small jitter width. Don't change height when jittering!
  scale_color_manual(values = c("#FF226F", "blue")) +
  labs(y="Shannon's diversity index", x = "Transcriptomics cluster group") +
  ggpubr::stat_compare_means(comparisons = list( c("White", "Black")),
                             method = "wilcox.test", hide.ns = T, paired=F, label="p.format", vjust=-0.7) +
  scale_x_discrete(labels=c("White" = "Cluster 1", "Black" = "Clusters 2-4")) +
  theme(legend.position="none",
        plot.title = element_text(size=9),
        text = element_text(size = 9), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
        ) 
ggsave("plots/PAPER_FIGURES/Supplement_DO_NOT_SUBMIT_THIS_FOLDER/Supplementary_Figure_1A_VIRGO.png", Supplementary_Figure_1A, height=3, width=2, dpi=300)


# Supplementary_Figure_1B - Relative abundance by CST
CST.labs <- c("CST I", "CST III", "CST IV", "NA")
names(CST.labs) <- c("CST I", "CST III", "CST IV-B", "Unclassified")
Supplementary_Figure_1B <- ggplot(VIRGO_f1_gather, aes(fill=Taxa, y=Abundance, x=IDs_1)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~Ravel, scales = "free", space = "free",
             labeller=labeller(Ravel=CST.labs)) +
  labs(x="IDs", y="Relative abundance (%)") +
  scale_fill_manual(values = c("#FF0000", "yellow", "#00CC99", #Ato, Dia, Gard
                               "#B2F5EB", "#02FEF4", "#99CCFF", #crisp, gas, iners
                               "#0099FF", "#000EA1", "#BF64FD", "black")) + #jen, vag Prev b,  Other
  theme(strip.text.x = element_text(size = 8), #facet labels text
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), #axes numbers text size
        legend.key.size = unit(5, 'mm'), #size of legend icons
        legend.text = element_text(face="italic")
  )
ggsave("plots/PAPER_FIGURES/Supplement_DO_NOT_SUBMIT_THIS_FOLDER/Supplementary_Figure_1B_VIRGO.pdf", Supplementary_Figure_1B, height=3, width=5, dpi=300)








###########################################################################################################################
##                                                                                                                       ##
##                                              Supplementary Table 3                                                    ##
##                                                                                                                       ##
###########################################################################################################################

# Supplementary Table 3 - Which genes to which module 
table(moduleColors)
# moduleColors
# bisque4           black            blue           brown          brown4            cyan       darkgreen        darkgrey     darkmagenta 
#      69             662             870             954             113             418             232             283             303 
# darkolivegreen      darkorange     darkorange2         darkred   darkturquoise     floralwhite           green     greenyellow            grey 
#            134             192             182             242             273              80             723             484               6 
# grey60           ivory       lightcyan      lightcyan1 lightsteelblue1     lightyellow         magenta   mediumpurple3    midnightblue 
#    318            1034             321             272              93             247             582              97             323 
# navajowhite2          orange  palevioletred3            pink           plum2          purple       royalblue          salmon         sienna3 
#           35             200              36             626              55             552             246             448             132 
# skyblue        skyblue3       steelblue             tan        thistle1       turquoise          violet           white     yellowgreen 
#     183             109             170             482              44            1972             152             923             110
dim(table(moduleColors)) # 45 

# Make the df
genes_in_modules <- as.data.frame(cbind(names(datExpr), moduleColors)) %>%
  dplyr::rename(Gene=V1,
                Module = moduleColors) %>% 
  arrange(Gene)



# Nicer df format 
genes_in_modules_nice <- genes_in_modules %>% group_by(Module) %>% 
  dplyr::mutate(Gene = paste(Gene, collapse=", ")) %>%
  dplyr::distinct(Module, Gene) %>% 
  dplyr::select(Module, Gene)%>% 
  arrange(Module)

# Checking
genes_in_modules %>% filter(Module=="grey") 
# Gene Module
# 1     ABAT   grey
# 2   CORO1B   grey
# 3      MVP   grey
# 4    PPM1G   grey
# 5   TCF7L2   grey
# 6 TRNAU1AP   grey
genes_in_modules_nice %>% filter(Module=="grey")
#grey   ABAT, CORO1B, MVP, PPM1G, TCF7L2, TRNAU1AP #matching 


# Export to csv
write_csv(genes_in_modules_nice, "data/WGCNA/Gene_module_assignment__Supplementary_Table_3.csv")
