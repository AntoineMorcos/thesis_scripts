# Clear environment
rm(list=ls())


# Set wd
#setwd("~/OneDrive - King's College London/PhD/Projects/RNA-seq") #mac
setwd("/Users/alici/OneDrive - King's College London/PhD/Projects/RNA-seq") #pc
getwd()
#list.files()

# Libraries
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(GenomicAlignments)
library(tibble)
library(mygene) #BiocManager::install("mygene")
library(gage) #BiocManager::install("gage")
library(gageData) #BiocManager::install("gageData")
library(pathview) #BiocManager::install("pathview")
library("AnnotationDbi") 
library("org.Hs.eg.db")
library(tidyverse) 
library(emdbook) #lseq function

# Versions of packages used
sessionInfo()





############################################################################################
##                                                                                        ##
##                          General prep for GAGE/Pathview/GO etc                         ##
##                                                                                        ##
############################################################################################

# This RData is from "Downstream_RNA-seq_DESeq2.R"
# Loading .Rdata saved
lnames = load(file="data/DESeq2_essentials.RData")
lnames


#https://rpubs.com/barryus/class15

# Main info:
# dds <- DESeqDataSetFromMatrix(countData=counts_CVF_adjusted, colData=metadata_CVF, design= ~ Outcome + Ethnicity) 
# dds1 <- DESeq(dds_no0s)
# res_Ethnicity <- results(dds1, contrast=c("Ethnicity", "White", "Black"), alpha=0.05, lfcThreshold = 0.58) 
# res_sLFC_Ethnicity <- lfcShrink(dds1, res=res_Ethnicity, contrast=c("Ethnicity", "White", "Black"), type="ashr")


# look
head(res_sLFC_Ethnicity_sort)
# log2 fold change (MMSE): Ethnicity White vs Black 
# Wald test p-value: Ethnicity White vs Black 
# DataFrame with 6 rows and 5 columns
#               baseMean log2FoldChange     lfcSE      pvalue        padj
#               <numeric>      <numeric> <numeric>   <numeric>   <numeric>
# SPRR2E       6713.3328       -15.9505   1.36948 8.01802e-31 7.15047e-27
# UCKL1-AS1     647.5065       -12.1325   1.40288 1.84526e-17 7.06467e-14
# CD300LD       680.7159       -12.9662   1.51744 2.37654e-17 7.06467e-14
# LOC105372803  326.7821       -11.4491   1.53128 1.70841e-13 3.80891e-10
# KLK10         204.3586        -9.5961   1.29628 9.38485e-13 1.67388e-09
# ELOVL7         80.7799        21.5165   3.19902 1.35293e-12 2.01091e-09


# Since we mapped and counted against the gene SYMBOL annotation, our res related to this. 
# However, our pathway analysis downstream will use KEGG pathways, and genes in KEGG pathways are annotated with Entrez gene IDs. 
# So we need to convert our gene symbols to Entrez.

# The list of gene IDs there are available
columns(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
# [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIPROT"


# Adding the entrez ID which is needed for GAGE
res_sLFC_Ethnicity_sort$entrez <-  mapIds(org.Hs.eg.db,
                                     keys=row.names(res_sLFC_Ethnicity_sort), 
                                     column="ENTREZID", #the new column of the gene type ID we want
                                     keytype="SYMBOL", #the ID in res already
                                     multiVals="first")

# Adding the symbol ID which is needed for Reactome Pathway Analysis - this was a bit weird as we could have just made it equal to row.names
res_sLFC_Ethnicity_sort$symbol <-  mapIds(org.Hs.eg.db,
                                          keys=row.names(res_sLFC_Ethnicity_sort), 
                                          column="SYMBOL", #the new column of the gene type ID we want
                                          keytype="SYMBOL", #the ID in res already
                                          multiVals="first")
#Check
table(res_sLFC_Ethnicity_sort$symbol==row.names(res_sLFC_Ethnicity_sort), useNA="always")

# Adding the GENENAME ID 
res_sLFC_Ethnicity_sort$name <-  mapIds(org.Hs.eg.db,
                                          keys=row.names(res_sLFC_Ethnicity_sort), 
                                          column="GENENAME", #the new column of the gene type ID we want
                                          keytype="SYMBOL", #the ID in res already
                                          multiVals="first")


# I checked to see if GO or ONTOLOGYALL gives useful GO info, but it does not



# Looking
head(res_sLFC_Ethnicity_sort)
# log2 fold change (MMSE): Ethnicity White vs Black 
# Wald test p-value: Ethnicity White vs Black 
# DataFrame with 6 rows and 8 columns
#               baseMean log2FoldChange     lfcSE      pvalue        padj      entrez       symbol                   name
#              <numeric>      <numeric> <numeric>   <numeric>   <numeric> <character>  <character>            <character>
# SPRR2E       6713.3328       -15.9505   1.36948 8.01802e-31 7.15047e-27        6704       SPRR2E small proline rich p..
# UCKL1-AS1     647.5065       -12.1325   1.40288 1.84526e-17 7.06467e-14   100113386    UCKL1-AS1  UCKL1 antisense RNA 1
# CD300LD       680.7159       -12.9662   1.51744 2.37654e-17 7.06467e-14   100131439      CD300LD CD300 molecule like ..
# LOC105372803  326.7821       -11.4491   1.53128 1.70841e-13 3.80891e-10   105372803 LOC105372803 uncharacterized LOC1..
# KLK10         204.3586        -9.5961   1.29628 9.38485e-13 1.67388e-09        5655        KLK10 kallikrein related p..
# ELOVL7         80.7799        21.5165   3.19902 1.35293e-12 2.01091e-09       79993       ELOVL7 ELOVL fatty acid elo..



#Renaming for ease
res <- res_sLFC_Ethnicity_sort

# Looking
table(is.na(res$entrez)) #300 NAs (gene symbols with no matching entrez)






############################################################################################
##                                                                                        ##
##                                   GO - Prep                                            ##
##                                                                                        ##
############################################################################################

# go.sets.hs has all GO terms. 
# go.subs.hs is a named list containing indexes for the BP, CC, and MF ontologies. 

#Creating foldchanges vector for GO
foldchanges <-  res$log2FoldChange
names(foldchanges) <-  res$entrez
head(foldchanges)
length(foldchanges)

#load in data
data(go.sets.hs)
data(go.subs.hs)

# GO meaning of BP, MF, CC http://geneontology.org/docs/ontology-documentation/ 

# Looking at function
?gage



############################################################################################
##                                                                                        ##
##                        GO -  biological process (BP)                                   ##
##                                                                                        ##
############################################################################################

#Subsetting GO by BP
go.subs.hs_BP <-  go.sets.hs[go.subs.hs$BP]

#Making GO results for BP
GO_BP_res <-  gage(foldchanges, gsets=go.subs.hs_BP, same.dir=TRUE, FDR.adj=T)
#FDR.adj	  Boolean, whether to do adjust for multiple testing as to control FDR (False dicovery rate). Default to be TRUE.
# I added the FDR.adj=T argument to check and the code is using the FDR adjustment 

# Looking at head of results -(Note on 22/12/21: I've seen that some of these values have changed slightly since I last ran this, so I am running it again) 
lapply(GO_BP_res, head)
# $greater
#                                                                 p.geomean stat.mean       p.val     q.val set.size        exp1
# GO:0034622 cellular macromolecular complex assembly           0.004470647  2.622455 0.004470647 0.7951297      430 0.004470647
# GO:0072594 establishment of protein localization to organelle 0.011588556  2.279997 0.011588556 0.7951297      198 0.011588556
# GO:0030099 myeloid cell differentiation                       0.018644541  2.091205 0.018644541 0.7951297      199 0.018644541
# GO:0006612 protein targeting to membrane                      0.024048405  1.985594 0.024048405 0.7951297      137 0.024048405
# GO:0065004 protein-DNA complex assembly                       0.024645275  1.982245 0.024645275 0.7951297       93 0.024645275
# GO:0071824 protein-DNA complex subunit organization           0.026154046  1.954273 0.026154046 0.7951297      111 0.026154046
# 
# $less
#                                                     p.geomean stat.mean      p.val     q.val set.size       exp1
# GO:0031424 keratinization                          0.02009586 -2.125674 0.02009586 0.9163428       25 0.02009586
# GO:0060429 epithelium development                  0.02283611 -2.002318 0.02283611 0.9163428      355 0.02283611
# GO:0007606 sensory perception of chemical stimulus 0.03870743 -1.819806 0.03870743 0.9163428       29 0.03870743
# GO:0007009 plasma membrane organization            0.04275650 -1.736261 0.04275650 0.9163428       67 0.04275650
# GO:0030216 keratinocyte differentiation            0.04749580 -1.684359 0.04749580 0.9163428       64 0.04749580
# GO:0030855 epithelial cell differentiation         0.05047650 -1.645212 0.05047650 0.9163428      168 0.05047650
# 
# $stats
#                                                               stat.mean     exp1
# GO:0034622 cellular macromolecular complex assembly            2.622455 2.622455
# GO:0072594 establishment of protein localization to organelle  2.279997 2.279997
# GO:0030099 myeloid cell differentiation                        2.091205 2.091205
# GO:0006612 protein targeting to membrane                       1.985594 1.985594
# GO:0065004 protein-DNA complex assembly                        1.982245 1.982245
# GO:0071824 protein-DNA complex subunit organization            1.954273 1.954273



############################################################### upregulated

# Making the upregulated results into a df
GO_BP_up <- base::as.data.frame(GO_BP_res$greater)
dim(GO_BP_up) #11971     6
head(GO_BP_up) #matching above list of upregulated

# subset by significant p-values 
GO_BP_up_p.val <-  subset(GO_BP_up, p.val < 0.05) 
dim(GO_BP_up_p.val)# [1] 23  6 

###############################################################



############################################################### downregulated

# Making the downregulated results into a df
GO_BP_down <- base::as.data.frame(GO_BP_res$less)
dim(GO_BP_down) #11971     6
head(GO_BP_down) #matching above list of downregulated

# subset by significant p-values 
GO_BP_down_p.val <-  subset(GO_BP_down, p.val < 0.05) 
dim(GO_BP_down_p.val)# [1] 5 6

###############################################################






############################################################################################
##                                                                                        ##
##                            GO - molecular function (MF)                                ##
##                                                                                        ##
############################################################################################

#Subsetting GO by MF
go.subs.hs_MF <-  go.sets.hs[go.subs.hs$MF]

#Making GO results for MF
GO_MF_res <-  gage(foldchanges, gsets=go.subs.hs_MF, same.dir=TRUE, FDR.adj=T)

# Looking at head of results
lapply(GO_MF_res, head)
# $greater
#                                                    p.geomean stat.mean      p.val     q.val set.size       exp1
# GO:0042277 peptide binding                        0.02023628  2.064035 0.02023628 0.7339902       91 0.02023628
# GO:0033218 amide binding                          0.02562756  1.962346 0.02562756 0.7339902       93 0.02562756
# GO:0008134 transcription factor binding           0.04154998  1.736548 0.04154998 0.7339902      321 0.04154998
# GO:0000975 regulatory region DNA binding          0.04951379  1.654799 0.04951379 0.7339902      216 0.04951379
# GO:0001067 regulatory region nucleic acid binding 0.04951379  1.654799 0.04951379 0.7339902      216 0.04951379
# GO:0047485 protein N-terminus binding             0.04976414  1.658094 0.04976414 0.7339902       73 0.04976414
# 
# $less
#                                                          p.geomean stat.mean      p.val     q.val set.size       exp1
# GO:0004867 serine-type endopeptidase inhibitor activity 0.02124104 -2.058398 0.02124104 0.9029011       53 0.02124104
# GO:0051015 actin filament binding                       0.03128320 -1.888139 0.03128320 0.9029011       54 0.03128320
# GO:0004252 serine-type endopeptidase activity           0.03152597 -1.872940 0.03152597 0.9029011       86 0.03152597
# GO:0097110 scaffold protein binding                     0.03358571 -1.957305 0.03358571 0.9029011       14 0.03358571
# GO:0019842 vitamin binding                              0.04434695 -1.734880 0.04434695 0.9029011       39 0.04434695
# GO:0008236 serine-type peptidase activity               0.04595517 -1.694637 0.04595517 0.9029011       99 0.04595517
# 
# $stats
#                                                   stat.mean     exp1
# GO:0042277 peptide binding                         2.064035 2.064035
# GO:0033218 amide binding                           1.962346 1.962346
# GO:0008134 transcription factor binding            1.736548 1.736548
# GO:0000975 regulatory region DNA binding           1.654799 1.654799
# GO:0001067 regulatory region nucleic acid binding  1.654799 1.654799
# GO:0047485 protein N-terminus binding              1.658094 1.658094


############################################################### upregulated

# Making the upregulated results into a df
GO_MF_up <- base::as.data.frame(GO_MF_res$greater)
dim(GO_MF_up) #3903    6
head(GO_MF_up) #matching above list of upregulated

# subset by significant p-values 
GO_MF_up_p.val <-  subset(GO_MF_up, p.val < 0.05) 
dim(GO_MF_up_p.val)# [1] 6 6

###############################################################



############################################################### downregulated

# Making the downregulated results into a df
GO_MF_down <- base::as.data.frame(GO_MF_res$less)
dim(GO_MF_down) #3903    6
head(GO_MF_down) #matching above list of downregulated

# subset by significant p-values 
GO_MF_down_p.val <-  subset(GO_MF_down, p.val < 0.05) 
dim(GO_MF_down_p.val)# [1] 7 6

###############################################################





############################################################################################
##                                                                                        ##
##                            GO - cellluar components (CC)                               ##
##                                                                                        ##
############################################################################################


#Subsetting GO by CC
go.subs.hs_CC <-  go.sets.hs[go.subs.hs$CC]

#Making GO results for CC
GO_CC_res <-  gage(foldchanges, gsets=go.subs.hs_CC, same.dir=TRUE, FDR.adj=T)

# Looking at head of results
lapply(GO_CC_res, head)
# $greater
#                                     p.geomean stat.mean      p.val     q.val set.size       exp1
# GO:0044427 chromosomal part        0.00643020  2.493960 0.00643020 0.8161889      439 0.00643020
# GO:0000793 condensed chromosome    0.03415314  1.834621 0.03415314 0.8161889      138 0.03415314
# GO:0005921 gap junction            0.03649837  1.876824 0.03649837 0.8161889       13 0.03649837
# GO:0000228 nuclear chromosome      0.03759182  1.784095 0.03759182 0.8161889      232 0.03759182
# GO:0044454 nuclear chromosome part 0.04285416  1.723227 0.04285416 0.8161889      206 0.04285416
# GO:0005778 peroxisomal membrane    0.04558301  1.712031 0.04558301 0.8161889       43 0.04558301
# 
# $less
#                                                 p.geomean stat.mean       p.val     q.val set.size        exp1
# GO:0005882 intermediate filament              0.002336653 -2.939124 0.002336653 0.9429398       41 0.002336653
# GO:0045111 intermediate filament cytoskeleton 0.010313834 -2.344068 0.010313834 0.9429398       75 0.010313834
# GO:0045095 keratin filament                   0.015999281 -2.329625 0.015999281 0.9429398       14 0.015999281
# GO:0001533 cornified envelope                 0.025525694 -2.045613 0.025525694 0.9429398       19 0.025525694
# GO:0031225 anchored to membrane               0.043435389 -1.728234 0.043435389 0.9429398       73 0.043435389
# GO:0016459 myosin complex                     0.088040829 -1.373602 0.088040829 0.9429398       33 0.088040829
# 
# $stats
#                                    stat.mean     exp1
# GO:0044427 chromosomal part         2.493960 2.493960
# GO:0000793 condensed chromosome     1.834621 1.834621
# GO:0005921 gap junction             1.876824 1.876824
# GO:0000228 nuclear chromosome       1.784095 1.784095
# GO:0044454 nuclear chromosome part  1.723227 1.723227
# GO:0005778 peroxisomal membrane     1.712031 1.712031


############################################################### upregulated

# Making the upregulated results into a df
GO_CC_up <- base::as.data.frame(GO_CC_res$greater)
dim(GO_CC_up) #1325    6
head(GO_CC_up) #matching above list of upregulated

# subset by significant p-values 
GO_CC_up_p.val <-  subset(GO_CC_up, p.val < 0.05) 
dim(GO_CC_up_p.val)# [1] 7 6

###############################################################



############################################################### downregulated

# Making the downregulated results into a df
GO_CC_down <- base::as.data.frame(GO_CC_res$less)
dim(GO_CC_down) #1325    6
head(GO_CC_down) #matching above list of downregulated

# subset by significant p-values 
GO_CC_down_p.val <-  subset(GO_CC_down, p.val < 0.05) 
dim(GO_CC_down_p.val)# [1] 5 6

###############################################################





############################################################################################
##                                                                                        ##
##                                GO - combine to one df                                  ##
##                                                                                        ##
############################################################################################


#GO_BP_up - ready for merge
GO_BP_up$GO_type <- "BP"
GO_BP_up$Direction <- "Upregulated"
GO_BP_up$GO_ID <- (substring(rownames(GO_BP_up), 1, 10)) #keep characters between 1 and 10 inclusive
GO_BP_up$GO_term <- (substring(rownames(GO_BP_up), 12)) #keep all characters from the 12th character

#GO_BP_down - ready for merge
GO_BP_down$GO_type <- "BP"
GO_BP_down$Direction <- "Downregulated"
GO_BP_down$GO_ID <- (substring(rownames(GO_BP_down), 1, 10)) #keep characters between 1 and 10 inclusive
GO_BP_down$GO_term <- (substring(rownames(GO_BP_down), 12)) #keep all characters from the 12th character

#GO_MF_up - ready for merge
GO_MF_up$GO_type <- "MF"
GO_MF_up$Direction <- "Upregulated"
GO_MF_up$GO_ID <- (substring(rownames(GO_MF_up), 1, 10)) #keep characters between 1 and 10 inclusive
GO_MF_up$GO_term <- (substring(rownames(GO_MF_up), 12)) #keep all characters from the 12th character

#GO_MF_down - ready for merge
GO_MF_down$GO_type <- "MF"
GO_MF_down$Direction <- "Downregulated"
GO_MF_down$GO_ID <- (substring(rownames(GO_MF_down), 1, 10)) #keep characters between 1 and 10 inclusive
GO_MF_down$GO_term <- (substring(rownames(GO_MF_down), 12)) #keep all characters from the 12th character

#GO_CC_up - ready for merge
GO_CC_up$GO_type <- "CC"
GO_CC_up$Direction <- "Upregulated"
GO_CC_up$GO_ID <- (substring(rownames(GO_CC_up), 1, 10)) #keep characters between 1 and 10 inclusive
GO_CC_up$GO_term <- (substring(rownames(GO_CC_up), 12)) #keep all characters from the 12th character

#GO_CC_down - ready for merge
GO_CC_down$GO_type <- "CC"
GO_CC_down$Direction <- "Downregulated"
GO_CC_down$GO_ID <- (substring(rownames(GO_CC_down), 1, 10)) #keep characters between 1 and 10 inclusive
GO_CC_down$GO_term <- (substring(rownames(GO_CC_down), 12)) #keep all characters from the 12th character

# Looking
head(GO_BP_up)
head(GO_BP_down)
head(GO_MF_up)
head(GO_MF_down)
head(GO_CC_up)
head(GO_CC_down)

# Merge
GO_ALL_sig_and_nonsig <- bind_rows(GO_BP_up,
                                   GO_BP_down,
                                   GO_MF_up,
                                   GO_MF_down,
                                   GO_CC_up,
                                   GO_CC_down)


# Data types
str(GO_ALL_sig_and_nonsig)
GO_ALL_sig_and_nonsig$Direction <- factor(GO_ALL_sig_and_nonsig$Direction, levels=c("Upregulated", "Downregulated"))
GO_ALL_sig_and_nonsig$GO_type <- factor(GO_ALL_sig_and_nonsig$GO_type, levels=c("BP","MF","CC"))

#Looking
#View(GO_ALL_sig_and_nonsig)
dim(GO_ALL_sig_and_nonsig)
table(GO_ALL_sig_and_nonsig$GO_type, GO_ALL_sig_and_nonsig$Direction)
#    Upregulated Downregulated
# BP       11971         11971
# MF        3903          3903
# CC        1325          1325
head(GO_ALL_sig_and_nonsig)
names(GO_ALL_sig_and_nonsig)


# export to CSV
write.csv(GO_ALL_sig_and_nonsig, 'data/GO/All_GO_significant_or_not.csv')





############################################################################################
##                                                                                        ##
##                            GO - ALL sig - prep for plotting                            ##
##                                                                                        ##
############################################################################################

#GO_BP_up_p.val - ready for merge
GO_BP_up_p.val$GO_type <- "BP"
GO_BP_up_p.val$Direction <- "Upregulated"
GO_BP_up_p.val$GO_ID <- (substring(rownames(GO_BP_up_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_BP_up_p.val$GO_term <- (substring(rownames(GO_BP_up_p.val), 12)) #keep all characters from the 12th character

#GO_BP_down_p.val - ready for merge
GO_BP_down_p.val$GO_type <- "BP"
GO_BP_down_p.val$Direction <- "Downregulated"
GO_BP_down_p.val$GO_ID <- (substring(rownames(GO_BP_down_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_BP_down_p.val$GO_term <- (substring(rownames(GO_BP_down_p.val), 12)) #keep all characters from the 12th character

#GO_MF_up_p.val - ready for merge
GO_MF_up_p.val$GO_type <- "MF"
GO_MF_up_p.val$Direction <- "Upregulated"
GO_MF_up_p.val$GO_ID <- (substring(rownames(GO_MF_up_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_MF_up_p.val$GO_term <- (substring(rownames(GO_MF_up_p.val), 12)) #keep all characters from the 12th character

#GO_MF_down_p.val - ready for merge
GO_MF_down_p.val$GO_type <- "MF"
GO_MF_down_p.val$Direction <- "Downregulated"
GO_MF_down_p.val$GO_ID <- (substring(rownames(GO_MF_down_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_MF_down_p.val$GO_term <- (substring(rownames(GO_MF_down_p.val), 12)) #keep all characters from the 12th character

#GO_CC_up_p.val - ready for merge
GO_CC_up_p.val$GO_type <- "CC"
GO_CC_up_p.val$Direction <- "Upregulated"
GO_CC_up_p.val$GO_ID <- (substring(rownames(GO_CC_up_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_CC_up_p.val$GO_term <- (substring(rownames(GO_CC_up_p.val), 12)) #keep all characters from the 12th character

#GO_CC_down_p.val - ready for merge
GO_CC_down_p.val$GO_type <- "CC"
GO_CC_down_p.val$Direction <- "Downregulated"
GO_CC_down_p.val$GO_ID <- (substring(rownames(GO_CC_down_p.val), 1, 10)) #keep characters between 1 and 10 inclusive
GO_CC_down_p.val$GO_term <- (substring(rownames(GO_CC_down_p.val), 12)) #keep all characters from the 12th character

# Looking
GO_BP_up_p.val
GO_BP_down_p.val
GO_MF_up_p.val
GO_MF_down_p.val
GO_CC_up_p.val
GO_CC_down_p.val
#                                                 p.geomean stat.mean       p.val     q.val set.size        exp1 GO_type     Direction      GO_ID                            GO_term
# GO:0005882 intermediate filament              0.002336653 -2.939124 0.002336653 0.9429398       41 0.002336653      CC Downregulated GO:0005882              intermediate filament
# GO:0045111 intermediate filament cytoskeleton 0.010313834 -2.344068 0.010313834 0.9429398       75 0.010313834      CC Downregulated GO:0045111 intermediate filament cytoskeleton
# GO:0045095 keratin filament                   0.015999281 -2.329625 0.015999281 0.9429398       14 0.015999281      CC Downregulated GO:0045095                   keratin filament
# GO:0001533 cornified envelope                 0.025525694 -2.045613 0.025525694 0.9429398       19 0.025525694      CC Downregulated GO:0001533                 cornified envelope
# GO:0031225 anchored to membrane               0.043435389 -1.728234 0.043435389 0.9429398       73 0.043435389      CC Downregulated GO:0031225               anchored to membrane

# Merge
GO_ALL <- bind_rows(GO_BP_up_p.val,
                    GO_BP_down_p.val,
                    GO_MF_up_p.val,
                    GO_MF_down_p.val,
                    GO_CC_up_p.val,
                    GO_CC_down_p.val)


# Data types
str(GO_ALL)
GO_ALL$Direction <- factor(GO_ALL$Direction, levels=c("Upregulated", "Downregulated"))
GO_ALL$GO_type <- factor(GO_ALL$GO_type, levels=c("BP","MF","CC"))

#Looking
#View(GO_ALL)
dim(GO_ALL) #53 10
table(GO_ALL$GO_type, GO_ALL$Direction)
#    Upregulated Downregulated
# BP          23             5
# MF           6             7
# CC           7             5
head(GO_ALL)
#                                                                 p.geomean stat.mean       p.val     q.val set.size        exp1 GO_type   Direction      GO_ID
# GO:0034622 cellular macromolecular complex assembly           0.004470647  2.622455 0.004470647 0.7951297      430 0.004470647      BP Upregulated GO:0034622
# GO:0072594 establishment of protein localization to organelle 0.011588556  2.279997 0.011588556 0.7951297      198 0.011588556      BP Upregulated GO:0072594
# GO:0030099 myeloid cell differentiation                       0.018644541  2.091205 0.018644541 0.7951297      199 0.018644541      BP Upregulated GO:0030099
# GO:0006612 protein targeting to membrane                      0.024048405  1.985594 0.024048405 0.7951297      137 0.024048405      BP Upregulated GO:0006612
# GO:0065004 protein-DNA complex assembly                       0.024645275  1.982245 0.024645275 0.7951297       93 0.024645275      BP Upregulated GO:0065004
# GO:0071824 protein-DNA complex subunit organization           0.026154046  1.954273 0.026154046 0.7951297      111 0.026154046      BP Upregulated GO:0071824
#                                                                                                          GO_term
# GO:0034622 cellular macromolecular complex assembly                     cellular macromolecular complex assembly
# GO:0072594 establishment of protein localization to organelle establishment of protein localization to organelle
# GO:0030099 myeloid cell differentiation                                             myeloid cell differentiation
# GO:0006612 protein targeting to membrane                                           protein targeting to membrane
# GO:0065004 protein-DNA complex assembly                                             protein-DNA complex assembly
# GO:0071824 protein-DNA complex subunit organization                     protein-DNA complex subunit organization
names(GO_ALL)






############################################################################################
##                                                                                        ##
##                              GO - ALL sig - plotting                                   ##
##                                                                                        ##
############################################################################################

# Plotting points 
GO.labs <- c("Biological process", "Molecular function", "Cellular component")
names(GO.labs) <- c("BP", "MF", "CC")
GO_DEGs_for_thesis <- ggplot(GO_ALL, aes(x=p.val, y=reorder(GO_term, -p.val)))  + #order GO_terms by p-values
  geom_point(aes(colour=Direction)) +
  facet_grid(GO_type~., scale="free", space="free",
             labeller=labeller(GO_type=GO.labs)) +
  labs(y="Gene ontology term", x="p-value") +
  xlim(0.05, 0) +
  scale_colour_manual(values = c("red", "blue"))
ggsave("plots/GO/GO_DEGs.pdf", GO_DEGs_for_thesis, height=8, width=8)

#Actually I decided not to do log scales as it just makes it hard to work out what p-value it is as there's not much of a range:
summary(GO_ALL$p.val)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002337 0.025526 0.035288 0.032903 0.042757 0.049764 





############################################################################################
##                                                                                        ##
##                                 GO - SELECTED sig                                      ##
##                                                                                        ##
############################################################################################

# Looking at all the sig GO terms
GO_ALL$GO_term

# df of just GO terms we want to show in the paper
GO_selected <- GO_ALL %>% 
  filter(GO_term=="immune system development" | GO_term=="cytokine production" | 
           GO_term=="circadian rhythm" | GO_term=="regulation of binding" | 
           GO_term=="cell cycle checkpoint" | GO_term=="phospholipid metabolic process" | 
           GO_term=="hematopoietic or lymphoid organ development" | GO_term=="DNA repair" | 
           GO_term=="hemopoiesis")


# Look
GO_selected$GO_term
summary(GO_selected)








###########################################################################################################################
##                                                                                                                       ##
##                                                Exporting all data                                                     ##
##                                                                                                                       ##
###########################################################################################################################

# Saving all data
save.image(file = "data/ALL_objects_in_GO_and_GAGE_Pathview_Pathway_analysis_of_ethnicity_DEGs.RData", compress=F)


# # Loading .Rdata saved
# lnames = load(file="data/ALL_objects_in_GO_and_GAGE_Pathview_Pathway_analysis_of_ethnicity_DEGs.RData")
# lnames




