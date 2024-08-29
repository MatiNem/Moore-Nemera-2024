library(data.table)
library(dplyr)
library(ggplot2)
install.packages("reshape")
library(reshape)
#function for annotating significance
sig_function <- function(x){
  if(x < 0.05){
    sigx <- "*"
  }
  else{
    sigx <- "ns"
  }
  if(x < 0.01){
    sigx <- "**"
  }
  if(x < 0.001){
    sigx <- "***"
  }
  if(x < 0.0001){
    sigx <- "****"
  }
  return(sigx)
}



all_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_flat.txt")

#gene body and flank mCA/CA of all genes in different cell types
#mPv_snmcseq_gene_and_flank = fread("HG_lab/Mati/GabelLab/genesets/pv_Russell/mPv_snmcseq_gene_and_flank_mCA.txt")
#mSst_snmcseq_gene_and_flank = fread("HG_lab/Mati/GabelLab/genesets/Sst/mSst_all_snmcseq_gene_and_flank_mCA.txt")
#mL4_snmcseq_gene_and_flank = fread("HG_lab/Mati/GabelLab/genesets/L4/mL4_snmcseq_gene_and_flank_mCA.txt")
#mL5_snmcseq_gene_and_flank = fread("HG_lab/Mati/GabelLab/genesets/L5/mL5_all_snmcseq_gene_and_flank_mCA.txt")
#gene body and flank mCG/CG of all genes in different cell types
#mPv_snmcseq_gene_and_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/pv_Russell/mPv_snmcseq_gene_and_flank_mCG.txt")
#mSst_snmcseq_gene_and_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/Sst/mSst_all_snmcseq_gene_and_flank_mCG.txt")
#mL4_snmcseq_gene_and_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/L4/mL4_snmcseq_gene_and_flank_mCG.txt")
#mL5_snmcseq_gene_and_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/L5/mL5_all_snmcseq_gene_and_flank_mCG.txt")
#gene lists
pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")


#cCREs
mousebrain_union_cCREs_PV_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_PV_mCA[, cCRE_methylation := V5/V6]


mousebrain_union_cCREs_PV_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_PV_mCG[, cCRE_methylation := V5/V6]

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

#subtracting nonconversion rate
mousebrain_union_cCREs_PV_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCA$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_PV_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCG$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]

#turn negative methylation values into zeros
mousebrain_union_cCREs_PV_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_PV_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]

#non-promoter cCREs
nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")
#cell-type cCREs
PVGA_nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/PVGA/PVGA_nonPromoter_cCREs_mm9.bed")

#union cCREs linked to genes
mousebrain_union_nonPromoter_cCREs_linked_genes = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_cCRE_Cicero_linked_genes_mm9.txt")
#cell-type cCREs linked to genes
PVGA_nonPromoter_cCREs_linked_genes = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_all_cCRE_Cicero_linked_genes_mm9.txt")[cCRE_label %in% nonPromoter_cCREs_mm9[,V4],]

#intragenic non-promoter cCREs
intragenic_nonPromoter_cCREs_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_ensgene_mm9.txt")



#Isolate the connection between cCREs and the genes in which they reside
cCRE_intragenic_connections = intragenic_nonPromoter_cCREs_mm9[, paste0(cCRE_label,"|",Gene)]
intragenic_nonPromoter_cCREs_mm9[, Conns:=cCRE_intragenic_connections]
#Identify whether putative enhancers are intragenic or not

mousebrain_union_nonPromoter_cCREs_linked_genes[!(cCRE_label %in% intragenic_nonPromoter_cCREs_mm9[, cCRE_label]), Intragenic := 0]
mousebrain_union_nonPromoter_cCREs_linked_genes[cCRE_label %in% intragenic_nonPromoter_cCREs_mm9[, cCRE_label], Intragenic := 1]
#Identify whether putative enhancers are inside gene to which they are linked
mousebrain_union_nonPromoter_cCREs_linked_genes[!(Conns %in% cCRE_intragenic_connections), Intragenic_to_linked_gene := 0]
mousebrain_union_nonPromoter_cCREs_linked_genes[Conns %in% cCRE_intragenic_connections, Intragenic_to_linked_gene := 1]


intragenic_nonPromoter_cCREs_mm9[Conns %in% mousebrain_union_nonPromoter_cCREs_linked_genes[, Conns], cCRE_label]



#all cCREs linked to MeCP2-regulated genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")


#intragenic cCREs cognate-linked to MeCP2-regulated genes
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

##intragenic cCREs linked to their cognate genes
#PVGA
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")


#nonPVGA
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

##extragenic cCREs linked to their cognate genes
#PVGA
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")


#nonPVGA
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

#INTACT PV mCA of PV and non-PV intragenic-linked cCREs
PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked = data.table(rbind(cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="PV cCREs,\n unchanged genes"),
                                                                                                                 cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Non-PV cCREs,\n unchanged genes"),
                                                                                                                 cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="PV cCREs,\n MR genes"),
                                                                                                                 cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Non-PV cCREs,\n MR genes")))

PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked, aes(x = cCRE_group, y = as.numeric(cCRE_methylation_corrected), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "cCREs linked to:", values = c("PV cCREs,\n unchanged genes"="gray", "Non-PV cCREs,\n unchanged genes"="gray", "PV cCREs,\n MR genes"="red", "Non-PV cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(0,0.165))+
  ylab("PV INTACT mCA/CA") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x = element_text(size=14, angle=90), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intragenicLinked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCAperCA_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intragenicLinked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCAperCA_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked[cCRE_group=="PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked[cCRE_group=="Non-PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)])$p.value #p=4.684482e-282, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked[cCRE_group=="PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_intragenicLinked[cCRE_group=="Non-PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)])$p.value #p=7.331446e-69, ****


#INTACT PV mCA of PV and non-PV extragenic-linked cCREs
PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked = data.table(rbind(cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="PV cCREs,\n unchanged genes"),
                                                                                                                        cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Non-PV cCREs,\n unchanged genes"),
                                                                                                                        cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="PV cCREs,\n MR genes"),
                                                                                                                        cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCA[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Non-PV cCREs,\n MR genes")))

PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked, aes(x = cCRE_group, y = as.numeric(cCRE_methylation_corrected), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "cCREs linked to:", values = c("PV cCREs,\n unchanged genes"="gray", "Non-PV cCREs,\n unchanged genes"="gray", "PV cCREs,\n MR genes"="red", "Non-PV cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(0,0.165))+
  ylab("PV INTACT mCA/CA") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x = element_text(size=14, angle=90), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_extragenicLinked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCAperCA_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_extragenicLinked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCAperCA_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked[cCRE_group=="PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked[cCRE_group=="Non-PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)])$p.value #p-value < 2.2e-16, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked[cCRE_group=="PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCA_extragenicLinked[cCRE_group=="Non-PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)])$p.value #p=1.057513e-41, ****



#plot of mCG of intragenic and extragenic linked cCREs of PV MR genes
PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra = data.table(rbind(cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Intragenic PV cCREs,\n unchanged genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Intragenic non-PV cCREs,\n unchanged genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Intragenic PV cCREs,\n MR genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Intragenic non-PV cCREs,\n MR genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Extragenic PV cCREs,\n unchanged genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Extragenic non-PV cCREs,\n unchanged genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Extragenic PV cCREs,\n MR genes"),
                                                                                                                cbind(cCRE_methylation_corrected=mousebrain_union_cCREs_PV_mCG[(V4 %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), cCRE_methylation_corrected], cCRE_group="Extragenic non-PV cCREs,\n MR genes")))

PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra = PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Intragenic PV cCREs,\n unchanged genes", "Intragenic non-PV cCREs,\n unchanged genes", "Intragenic PV cCREs,\n MR genes", "Intragenic non-PV cCREs,\n MR genes",
                                                                                                                                                                                                                                                              "Extragenic PV cCREs,\n unchanged genes", "Extragenic non-PV cCREs,\n unchanged genes", "Extragenic PV cCREs,\n MR genes", "Extragenic non-PV cCREs,\n MR genes")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra, aes(x = cCRE_group, y = as.numeric(cCRE_methylation_corrected), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "cCREs:", values = c("Intragenic PV cCREs,\n unchanged genes"="gray", "Intragenic non-PV cCREs,\n unchanged genes"="gray", "Intragenic PV cCREs,\n MR genes"="red", "Intragenic non-PV cCREs,\n MR genes"="red",
                                                "Extragenic PV cCREs,\n unchanged genes"="gray", "Extragenic non-PV cCREs,\n unchanged genes"="gray", "Extragenic PV cCREs,\n MR genes"="red", "Extragenic non-PV cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(0,1))+
  ylab("PV INTACT mCG/CG") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x = element_text(size=14, angle=90), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intra_and_extra_linked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCGperCG_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intra_and_extra_linked_to_PV_MR_genes_coding_prefilt5_nondedup_PV_WT_KO_deep_INTACT_mCGperCG_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

options(scipen = 0)

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Intragenic PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Intragenic non-PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)])$p.value #p-value < 2.2e-16, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Intragenic PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Intragenic non-PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)])$p.value #p = 1.272077e-106, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Extragenic PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Extragenic non-PV cCREs,\n unchanged genes", as.numeric(cCRE_methylation_corrected)])$p.value #0, p-value < 2.2e-16, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Extragenic PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)], PVGA_and_nonPVGA_nonPromoter_cCRE_linked_to_Pv_MR_genes_noOther_nondedup_INTACT_mCG_intra_and_extra[cCRE_group=="Extragenic non-PV cCREs,\n MR genes", as.numeric(cCRE_methylation_corrected)])$p.value #p = 6.503104e-56, ****


