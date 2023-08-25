############################ Info
# In this script I wanted to check that the WGCNA data from July is the same as June to see if I need to redo cytoscape network. Specifically, I look at the files of which taxa are assigned to each module, and what is the weight b/ each taxa for input to cytoscape.
############################ 


############################ Taxa module assignment
cd "/c/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines/data/WGCNA"
diff taxa_module_assignment_Hel_sft10_mm5_June2022.csv taxa_module_assignment_Hel_sft10_mm5_July2022.csv #no output so no diff 
cksum taxa_module_assignment_Hel_sft10_mm5_Ju*
# 3663046562 6645 taxa_module_assignment_Hel_sft10_mm5_July2022.csv
# 3663046562 6645 taxa_module_assignment_Hel_sft10_mm5_June2022.csv

diff taxa_module_assignment_Hel_sft10_mm5_Dec2022.csv taxa_module_assignment_Hel_sft10_mm5_July2022.csv
#1c1
#< Taxa,Module
#---
#> Taxa,Module_when_5minModuleSize 
# So the only difference is the heading name 
############################ 



############################ Cytoscape input files
cd "/c/Users/alici/OneDrive - King's College London/PhD/Projects/metagenome_cytokines/Cytoscape"
diff CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_June2022.txt CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_July2022.txt #no output so no diff
diff CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_June2022.txt CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_July2022.txt #no output so no diff
cksum *
# 2348760675 179938 CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_July2022.txt
# 2348760675 179938 CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_June2022.txt
# 3252775961 6650 CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_July2022.txt
# 3252775961 6650 CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_June2022.txt
# 1479052959 54681 Cytoscape_June2022.pdf
# 2869612407 387340 Cytoscape_annotated_June2022.pdf
# 3493057282 135983 cyto_WGCNA_on_metagenome_n87_June2022.cys
# cksum: old: Is a directory
md5sum *
# 7fd899e136e7c666457af68d5206c1b9 *CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_July2022.txt
# 7fd899e136e7c666457af68d5206c1b9 *CytoInput_all_mods_thresh0.05_edges__Hel_sft12_mm5_June2022.txt
# c2e879a0770266ee23823726241fcaf1 *CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_July2022.txt
# c2e879a0770266ee23823726241fcaf1 *CytoInput_all_mods_thresh0.05_nodes__Hel_sft12_mm5_June2022.txt
# 84ba6c7c88c69115217dd1730cb9be9b *Cytoscape_June2022.pdf
# 60d954134595c9e1df47766a7fa1bc0a *Cytoscape_annotated_June2022.pdf
# d866a45f1fe545ca3b2a52d36849371f *cyto_WGCNA_on_metagenome_n87_June2022.cys
# md5sum: old: Is a directory
############################ 
