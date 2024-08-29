library(data.table)
library(dplyr)
install.packages("plotrix")
library(plotrix)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
union_mr_cCREs = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
union_ma_cCREs = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
union_all_cCREs = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

union_nonMR_cCREs = union_all_cCREs[!(V4 %in% union_mr_cCREs[,V4]), ]
union_nonMA_cCREs = union_all_cCREs[!(V4 %in% union_ma_cCREs[,V4]), ]

MR_genes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
MA_genes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")

mousebrain_union_cCREs_linked_genes = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_cCRE_Cicero_linked_genes_mm9.txt")
mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(cCRE_label %in% nonPromoter_cCREs_mm9[,V4]) & (Intragenic==1) & (Intragenic_to_linked_gene==1), ]
mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(cCRE_label %in% nonPromoter_cCREs_mm9[,V4]) & (Intragenic==0) & (Intragenic_to_linked_gene==0), ]

intragenic_nonPromoter_cCREs_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_ensgene_mm9.txt")
cCRE_gene_same_PvTAD = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_ensgene_mm9_jrmCE_PvTADs_Arrowhead40kb_KR_sameTADs.txt")


#intragenic
p_values_MR_cCRE_MR_genes_intragenic = rep(0,1000)
log2fc_MR_cCRE_MR_genes_intragenic = rep(0,1000)
p_values_MR_cCRE_MA_genes_intragenic = rep(0,1000)
log2fc_MR_cCRE_MA_genes_intragenic = rep(0,1000)
p_values_MA_cCRE_MR_genes_intragenic = rep(0,1000)
log2fc_MA_cCRE_MR_genes_intragenic = rep(0,1000)
p_values_MA_cCRE_MA_genes_intragenic = rep(0,1000)
log2fc_MA_cCRE_MA_genes_intragenic = rep(0,1000)
#Cicero
p_values_MR_cCRE_MR_genes_Cicero = rep(0,1000)
log2fc_MR_cCRE_MR_genes_Cicero = rep(0,1000)
p_values_MR_cCRE_MA_genes_Cicero = rep(0,1000)
log2fc_MR_cCRE_MA_genes_Cicero = rep(0,1000)
p_values_MA_cCRE_MR_genes_Cicero = rep(0,1000)
log2fc_MA_cCRE_MR_genes_Cicero = rep(0,1000)
p_values_MA_cCRE_MA_genes_Cicero = rep(0,1000)
log2fc_MA_cCRE_MA_genes_Cicero = rep(0,1000)
#same TAD
p_values_MR_cCRE_MR_genes_sameTAD = rep(0,1000)
log2fc_MR_cCRE_MR_genes_sameTAD = rep(0,1000)
p_values_MR_cCRE_MA_genes_sameTAD = rep(0,1000)
log2fc_MR_cCRE_MA_genes_sameTAD = rep(0,1000)
p_values_MA_cCRE_MR_genes_sameTAD = rep(0,1000)
log2fc_MA_cCRE_MR_genes_sameTAD = rep(0,1000)
p_values_MA_cCRE_MA_genes_sameTAD = rep(0,1000)
log2fc_MA_cCRE_MA_genes_sameTAD = rep(0,1000)
#Intragenic, Cicero cognate-linked
p_values_MR_cCRE_MR_genes_intragenicLinked = rep(0,1000)
log2fc_MR_cCRE_MR_genes_intragenicLinked = rep(0,1000)
p_values_MR_cCRE_MA_genes_intragenicLinked = rep(0,1000)
log2fc_MR_cCRE_MA_genes_intragenicLinked = rep(0,1000)
p_values_MA_cCRE_MR_genes_intragenicLinked = rep(0,1000)
log2fc_MA_cCRE_MR_genes_intragenicLinked = rep(0,1000)
p_values_MA_cCRE_MA_genes_intragenicLinked = rep(0,1000)
log2fc_MA_cCRE_MA_genes_intragenicLinked = rep(0,1000)
#extragenic, Cicero-linked
p_values_MR_cCRE_MR_genes_extragenicLinked = rep(0,1000)
log2fc_MR_cCRE_MR_genes_extragenicLinked = rep(0,1000)
p_values_MR_cCRE_MA_genes_extragenicLinked = rep(0,1000)
log2fc_MR_cCRE_MA_genes_extragenicLinked = rep(0,1000)
p_values_MA_cCRE_MR_genes_extragenicLinked = rep(0,1000)
log2fc_MA_cCRE_MR_genes_extragenicLinked = rep(0,1000)
p_values_MA_cCRE_MA_genes_extragenicLinked = rep(0,1000)
log2fc_MA_cCRE_MA_genes_extragenicLinked = rep(0,1000)

set.seed(seed=1)
for(i in 1:1000){
  Pv_mr_genes_PvTPM_MeCP2KO_df_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_df_resamp/resamp", i,".txt"), header=F)$V4
  Pv_ma_genes_PvTPM_MeCP2KO_df_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_ma_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_df_resamp/resamp", i,".txt"), header=F)$V4
  #MR cCREs in MR genes
  MR_cCRE_MR_genes_intragenic = length(intersect(union_mr_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MR_cCRE_MR_genes_intragenic_resamp = length(intersect(union_mr_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MR_genes_intragenic = length(intersect(union_nonMR_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MR_genes_intragenic_resamp = length(intersect(union_nonMR_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MR_genes_intragenic = data.frame(c(MR_cCRE_MR_genes_intragenic, nonMR_cCRE_MR_genes_intragenic),
                                               c(MR_cCRE_MR_genes_intragenic_resamp, nonMR_cCRE_MR_genes_intragenic_resamp), 
                                               row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MR_genes_intragenic) = c("MR_gene", "Resampled_gene")
  test_MR_cCRE_MR_genes_intragenic = fisher.test(dat_MR_cCRE_MR_genes_intragenic)
  p_values_MR_cCRE_MR_genes_intragenic[i] =  test_MR_cCRE_MR_genes_intragenic$p.value * sign(log2(test_MR_cCRE_MR_genes_intragenic$estimate))
  log2fc_MR_cCRE_MR_genes_intragenic[i] = log2(test_MR_cCRE_MR_genes_intragenic$estimate)
  
  #MR cCREs in MA genes
  MR_cCRE_MA_genes_intragenic = length(intersect(union_mr_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MR_cCRE_MA_genes_intragenic_resamp = length(intersect(union_mr_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MA_genes_intragenic = length(intersect(union_nonMR_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MA_genes_intragenic_resamp = length(intersect(union_nonMR_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MA_genes_intragenic = data.frame(c(MR_cCRE_MA_genes_intragenic, nonMR_cCRE_MA_genes_intragenic),
                                               c(MR_cCRE_MA_genes_intragenic_resamp, nonMR_cCRE_MA_genes_intragenic_resamp), 
                                               row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MA_genes_intragenic) = c("MA_gene", "Resampled_gene")
  test_MR_cCRE_MA_genes_intragenic = fisher.test(dat_MR_cCRE_MA_genes_intragenic)
  p_values_MR_cCRE_MA_genes_intragenic[i] =  test_MR_cCRE_MA_genes_intragenic$p.value * sign(log2(test_MR_cCRE_MA_genes_intragenic$estimate))
  log2fc_MR_cCRE_MA_genes_intragenic[i] = log2(test_MR_cCRE_MA_genes_intragenic$estimate)
  
  #MA cCREs in MR genes
  MA_cCRE_MR_genes_intragenic = length(intersect(union_ma_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MA_cCRE_MR_genes_intragenic_resamp = length(intersect(union_ma_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MR_genes_intragenic = length(intersect(union_nonMA_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MR_genes_intragenic_resamp = length(intersect(union_nonMA_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MR_genes_intragenic = data.frame(c(MA_cCRE_MR_genes_intragenic, nonMA_cCRE_MR_genes_intragenic),
                                               c(MA_cCRE_MR_genes_intragenic_resamp, nonMA_cCRE_MR_genes_intragenic_resamp), 
                                               row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MR_genes_intragenic) = c("MR_gene", "Resampled_gene")
  test_MA_cCRE_MR_genes_intragenic = fisher.test(dat_MA_cCRE_MR_genes_intragenic)
  p_values_MA_cCRE_MR_genes_intragenic[i] =  test_MA_cCRE_MR_genes_intragenic$p.value * sign(log2(test_MA_cCRE_MR_genes_intragenic$estimate))
  log2fc_MA_cCRE_MR_genes_intragenic[i] = log2(test_MA_cCRE_MR_genes_intragenic$estimate)
  #MA cCREs in MA genes
  MA_cCRE_MA_genes_intragenic = length(intersect(union_ma_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MA_cCRE_MA_genes_intragenic_resamp = length(intersect(union_ma_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MA_genes_intragenic = length(intersect(union_nonMA_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MA_genes_intragenic_resamp = length(intersect(union_nonMA_cCREs[, V4], intragenic_nonPromoter_cCREs_mm9[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MA_genes_intragenic = data.frame(c(MA_cCRE_MA_genes_intragenic, nonMA_cCRE_MA_genes_intragenic),
                                               c(MA_cCRE_MA_genes_intragenic_resamp, nonMA_cCRE_MA_genes_intragenic_resamp), 
                                               row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MA_genes_intragenic) = c("MA_gene", "Resampled_gene")
  test_MA_cCRE_MA_genes_intragenic = fisher.test(dat_MA_cCRE_MA_genes_intragenic)
  p_values_MA_cCRE_MA_genes_intragenic[i] =  test_MA_cCRE_MA_genes_intragenic$p.value * sign(log2(test_MA_cCRE_MA_genes_intragenic$estimate))
  log2fc_MA_cCRE_MA_genes_intragenic[i] = log2(test_MA_cCRE_MA_genes_intragenic$estimate)
  
  #MR cCREs linked to MR genes
  MR_cCRE_MR_genes_Cicero = length(intersect(union_mr_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MR_cCRE_MR_genes_Cicero_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MR_genes_Cicero = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MR_genes_Cicero_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MR_genes_Cicero = data.frame(c(MR_cCRE_MR_genes_Cicero, nonMR_cCRE_MR_genes_Cicero),
                                           c(MR_cCRE_MR_genes_Cicero_resamp, nonMR_cCRE_MR_genes_Cicero_resamp), 
                                           row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MR_genes_Cicero) = c("MR_gene", "Resampled_gene")
  test_MR_cCRE_MR_genes_Cicero = fisher.test(dat_MR_cCRE_MR_genes_Cicero)
  p_values_MR_cCRE_MR_genes_Cicero[i] =  test_MR_cCRE_MR_genes_Cicero$p.value * sign(log2(test_MR_cCRE_MR_genes_Cicero$estimate))
  log2fc_MR_cCRE_MR_genes_Cicero[i] = log2(test_MR_cCRE_MR_genes_Cicero$estimate)
  
  #MR cCREs linked to MA genes
  MR_cCRE_MA_genes_Cicero = length(intersect(union_mr_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MR_cCRE_MA_genes_Cicero_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MA_genes_Cicero = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MA_genes_Cicero_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MA_genes_Cicero = data.frame(c(MR_cCRE_MA_genes_Cicero, nonMR_cCRE_MA_genes_Cicero),
                                           c(MR_cCRE_MA_genes_Cicero_resamp, nonMR_cCRE_MA_genes_Cicero_resamp), 
                                           row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MA_genes_Cicero) = c("MA_gene", "Resampled_gene")
  test_MR_cCRE_MA_genes_Cicero = fisher.test(dat_MR_cCRE_MA_genes_Cicero)
  p_values_MR_cCRE_MA_genes_Cicero[i] =  test_MR_cCRE_MA_genes_Cicero$p.value * sign(log2(test_MR_cCRE_MA_genes_Cicero$estimate))
  log2fc_MR_cCRE_MA_genes_Cicero[i] = log2(test_MR_cCRE_MA_genes_Cicero$estimate)
  
  #MA cCREs linked to MR genes
  MA_cCRE_MR_genes_Cicero = length(intersect(union_ma_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MA_cCRE_MR_genes_Cicero_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MR_genes_Cicero = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MR_genes_Cicero_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MR_genes_Cicero = data.frame(c(MA_cCRE_MR_genes_Cicero, nonMA_cCRE_MR_genes_Cicero),
                                           c(MA_cCRE_MR_genes_Cicero_resamp, nonMA_cCRE_MR_genes_Cicero_resamp), 
                                           row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MR_genes_Cicero) = c("MR_gene", "Resampled_gene")
  test_MA_cCRE_MR_genes_Cicero = fisher.test(dat_MA_cCRE_MR_genes_Cicero)
  p_values_MA_cCRE_MR_genes_Cicero[i] =  test_MA_cCRE_MR_genes_Cicero$p.value * sign(log2(test_MA_cCRE_MR_genes_Cicero$estimate))
  log2fc_MA_cCRE_MR_genes_Cicero[i] = log2(test_MA_cCRE_MR_genes_Cicero$estimate)
  #MA cCREs linked to MA genes
  MA_cCRE_MA_genes_Cicero = length(intersect(union_ma_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MA_cCRE_MA_genes_Cicero_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MA_genes_Cicero = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MA_genes_Cicero_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_cCREs_linked_genes[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MA_genes_Cicero = data.frame(c(MA_cCRE_MA_genes_Cicero, nonMA_cCRE_MA_genes_Cicero),
                                           c(MA_cCRE_MA_genes_Cicero_resamp, nonMA_cCRE_MA_genes_Cicero_resamp), 
                                           row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MA_genes_Cicero) = c("MA_gene", "Resampled_gene")
  test_MA_cCRE_MA_genes_Cicero = fisher.test(dat_MA_cCRE_MA_genes_Cicero)
  p_values_MA_cCRE_MA_genes_Cicero[i] =  test_MA_cCRE_MA_genes_Cicero$p.value * sign(log2(test_MA_cCRE_MA_genes_Cicero$estimate))
  log2fc_MA_cCRE_MA_genes_Cicero[i] = log2(test_MA_cCRE_MA_genes_Cicero$estimate)
  
  #MR cCREs in same TAD as MR genes
  MR_cCRE_MR_genes_sameTAD = length(intersect(union_mr_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MR_cCRE_MR_genes_sameTAD_resamp = length(intersect(union_mr_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MR_genes_sameTAD = length(intersect(union_nonMR_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MR_genes_sameTAD_resamp = length(intersect(union_nonMR_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MR_genes_sameTAD = data.frame(c(MR_cCRE_MR_genes_sameTAD, nonMR_cCRE_MR_genes_sameTAD),
                                            c(MR_cCRE_MR_genes_sameTAD_resamp, nonMR_cCRE_MR_genes_sameTAD_resamp), 
                                            row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MR_genes_sameTAD) = c("MR_gene", "Resampled_gene")
  test_MR_cCRE_MR_genes_sameTAD = fisher.test(dat_MR_cCRE_MR_genes_sameTAD)
  p_values_MR_cCRE_MR_genes_sameTAD[i] =  test_MR_cCRE_MR_genes_sameTAD$p.value * sign(log2(test_MR_cCRE_MR_genes_sameTAD$estimate))
  log2fc_MR_cCRE_MR_genes_sameTAD[i] = log2(test_MR_cCRE_MR_genes_sameTAD$estimate)
  
  #MR cCREs in same TAD as MA genes
  MR_cCRE_MA_genes_sameTAD = length(intersect(union_mr_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MR_cCRE_MA_genes_sameTAD_resamp = length(intersect(union_mr_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MA_genes_sameTAD = length(intersect(union_nonMR_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MA_genes_sameTAD_resamp = length(intersect(union_nonMR_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MA_genes_sameTAD = data.frame(c(MR_cCRE_MA_genes_sameTAD, nonMR_cCRE_MA_genes_sameTAD),
                                            c(MR_cCRE_MA_genes_sameTAD_resamp, nonMR_cCRE_MA_genes_sameTAD_resamp), 
                                            row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MA_genes_sameTAD) = c("MA_gene", "Resampled_gene")
  test_MR_cCRE_MA_genes_sameTAD = fisher.test(dat_MR_cCRE_MA_genes_sameTAD)
  p_values_MR_cCRE_MA_genes_sameTAD[i] =  test_MR_cCRE_MA_genes_sameTAD$p.value * sign(log2(test_MR_cCRE_MA_genes_sameTAD$estimate))
  log2fc_MR_cCRE_MA_genes_sameTAD[i] = log2(test_MR_cCRE_MA_genes_sameTAD$estimate)
  
  #MA cCREs in same TAD as MR genes
  MA_cCRE_MR_genes_sameTAD = length(intersect(union_ma_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MA_cCRE_MR_genes_sameTAD_resamp = length(intersect(union_ma_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MR_genes_sameTAD = length(intersect(union_nonMA_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MR_genes_sameTAD_resamp = length(intersect(union_nonMA_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MR_genes_sameTAD = data.frame(c(MA_cCRE_MR_genes_sameTAD, nonMA_cCRE_MR_genes_sameTAD),
                                            c(MA_cCRE_MR_genes_sameTAD_resamp, nonMA_cCRE_MR_genes_sameTAD_resamp), 
                                            row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MR_genes_sameTAD) = c("MR_gene", "Resampled_gene")
  test_MA_cCRE_MR_genes_sameTAD = fisher.test(dat_MA_cCRE_MR_genes_sameTAD)
  p_values_MA_cCRE_MR_genes_sameTAD[i] =  test_MA_cCRE_MR_genes_sameTAD$p.value * sign(log2(test_MA_cCRE_MR_genes_sameTAD$estimate))
  log2fc_MA_cCRE_MR_genes_sameTAD[i] = log2(test_MA_cCRE_MR_genes_sameTAD$estimate)
  #MA cCREs in same TAD as MA genes
  MA_cCRE_MA_genes_sameTAD = length(intersect(union_ma_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MA_cCRE_MA_genes_sameTAD_resamp = length(intersect(union_ma_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MA_genes_sameTAD = length(intersect(union_nonMA_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MA_genes_sameTAD_resamp = length(intersect(union_nonMA_cCREs[, V4], cCRE_gene_same_PvTAD[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MA_genes_sameTAD = data.frame(c(MA_cCRE_MA_genes_sameTAD, nonMA_cCRE_MA_genes_sameTAD),
                                            c(MA_cCRE_MA_genes_sameTAD_resamp, nonMA_cCRE_MA_genes_sameTAD_resamp), 
                                            row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MA_genes_sameTAD) = c("MA_gene", "Resampled_gene")
  test_MA_cCRE_MA_genes_sameTAD = fisher.test(dat_MA_cCRE_MA_genes_sameTAD)
  p_values_MA_cCRE_MA_genes_sameTAD[i] =  test_MA_cCRE_MA_genes_sameTAD$p.value * sign(log2(test_MA_cCRE_MA_genes_sameTAD$estimate))
  log2fc_MA_cCRE_MA_genes_sameTAD[i] = log2(test_MA_cCRE_MA_genes_sameTAD$estimate)
  
  #MR cCREs cognate-linked to MR genes
  MR_cCRE_MR_genes_intragenicLinked = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MR_cCRE_MR_genes_intragenicLinked_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MR_genes_intragenicLinked = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MR_genes_intragenicLinked_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MR_genes_intragenicLinked = data.frame(c(MR_cCRE_MR_genes_intragenicLinked, nonMR_cCRE_MR_genes_intragenicLinked),
                                           c(MR_cCRE_MR_genes_intragenicLinked_resamp, nonMR_cCRE_MR_genes_intragenicLinked_resamp), 
                                           row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MR_genes_intragenicLinked) = c("MR_gene", "Resampled_gene")
  test_MR_cCRE_MR_genes_intragenicLinked = fisher.test(dat_MR_cCRE_MR_genes_intragenicLinked)
  p_values_MR_cCRE_MR_genes_intragenicLinked[i] =  test_MR_cCRE_MR_genes_intragenicLinked$p.value * sign(log2(test_MR_cCRE_MR_genes_intragenicLinked$estimate))
  log2fc_MR_cCRE_MR_genes_intragenicLinked[i] = log2(test_MR_cCRE_MR_genes_intragenicLinked$estimate)
  
  #MR cCREs cognate-linked to MA genes
  MR_cCRE_MA_genes_intragenicLinked = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MR_cCRE_MA_genes_intragenicLinked_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MA_genes_intragenicLinked = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MA_genes_intragenicLinked_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MA_genes_intragenicLinked = data.frame(c(MR_cCRE_MA_genes_intragenicLinked, nonMR_cCRE_MA_genes_intragenicLinked),
                                           c(MR_cCRE_MA_genes_intragenicLinked_resamp, nonMR_cCRE_MA_genes_intragenicLinked_resamp), 
                                           row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MA_genes_intragenicLinked) = c("MA_gene", "Resampled_gene")
  test_MR_cCRE_MA_genes_intragenicLinked = fisher.test(dat_MR_cCRE_MA_genes_intragenicLinked)
  p_values_MR_cCRE_MA_genes_intragenicLinked[i] =  test_MR_cCRE_MA_genes_intragenicLinked$p.value * sign(log2(test_MR_cCRE_MA_genes_intragenicLinked$estimate))
  log2fc_MR_cCRE_MA_genes_intragenicLinked[i] = log2(test_MR_cCRE_MA_genes_intragenicLinked$estimate)
  
  #MA cCREs cognate-linked to MR genes
  MA_cCRE_MR_genes_intragenicLinked = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MA_cCRE_MR_genes_intragenicLinked_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MR_genes_intragenicLinked = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MR_genes_intragenicLinked_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MR_genes_intragenicLinked = data.frame(c(MA_cCRE_MR_genes_intragenicLinked, nonMA_cCRE_MR_genes_intragenicLinked),
                                           c(MA_cCRE_MR_genes_intragenicLinked_resamp, nonMA_cCRE_MR_genes_intragenicLinked_resamp), 
                                           row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MR_genes_intragenicLinked) = c("MR_gene", "Resampled_gene")
  test_MA_cCRE_MR_genes_intragenicLinked = fisher.test(dat_MA_cCRE_MR_genes_intragenicLinked)
  p_values_MA_cCRE_MR_genes_intragenicLinked[i] =  test_MA_cCRE_MR_genes_intragenicLinked$p.value * sign(log2(test_MA_cCRE_MR_genes_intragenicLinked$estimate))
  log2fc_MA_cCRE_MR_genes_intragenicLinked[i] = log2(test_MA_cCRE_MR_genes_intragenicLinked$estimate)
  #MA cCREs cognate-linked to MA genes
  MA_cCRE_MA_genes_intragenicLinked = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MA_cCRE_MA_genes_intragenicLinked_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MA_genes_intragenicLinked = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MA_genes_intragenicLinked_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MA_genes_intragenicLinked = data.frame(c(MA_cCRE_MA_genes_intragenicLinked, nonMA_cCRE_MA_genes_intragenicLinked),
                                           c(MA_cCRE_MA_genes_intragenicLinked_resamp, nonMA_cCRE_MA_genes_intragenicLinked_resamp), 
                                           row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MA_genes_intragenicLinked) = c("MA_gene", "Resampled_gene")
  test_MA_cCRE_MA_genes_intragenicLinked = fisher.test(dat_MA_cCRE_MA_genes_intragenicLinked)
  p_values_MA_cCRE_MA_genes_intragenicLinked[i] =  test_MA_cCRE_MA_genes_intragenicLinked$p.value * sign(log2(test_MA_cCRE_MA_genes_intragenicLinked$estimate))
  log2fc_MA_cCRE_MA_genes_intragenicLinked[i] = log2(test_MA_cCRE_MA_genes_intragenicLinked$estimate)
  
  ##extragenic
  #MR cCREs extragenic and Cicero-linked to MR genes
  MR_cCRE_MR_genes_extragenicLinked = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MR_cCRE_MR_genes_extragenicLinked_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MR_genes_extragenicLinked = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MR_genes_extragenicLinked_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MR_genes_extragenicLinked = data.frame(c(MR_cCRE_MR_genes_extragenicLinked, nonMR_cCRE_MR_genes_extragenicLinked),
                                                     c(MR_cCRE_MR_genes_extragenicLinked_resamp, nonMR_cCRE_MR_genes_extragenicLinked_resamp), 
                                                     row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MR_genes_extragenicLinked) = c("MR_gene", "Resampled_gene")
  test_MR_cCRE_MR_genes_extragenicLinked = fisher.test(dat_MR_cCRE_MR_genes_extragenicLinked)
  p_values_MR_cCRE_MR_genes_extragenicLinked[i] =  test_MR_cCRE_MR_genes_extragenicLinked$p.value * sign(log2(test_MR_cCRE_MR_genes_extragenicLinked$estimate))
  log2fc_MR_cCRE_MR_genes_extragenicLinked[i] = log2(test_MR_cCRE_MR_genes_extragenicLinked$estimate)
  
  #MR cCREs extragenic and Cicero-linked to MA genes
  MR_cCRE_MA_genes_extragenicLinked = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MR_cCRE_MA_genes_extragenicLinked_resamp = length(intersect(union_mr_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMR_cCRE_MA_genes_extragenicLinked = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMR_cCRE_MA_genes_extragenicLinked_resamp = length(intersect(union_nonMR_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MR_cCRE_MA_genes_extragenicLinked = data.frame(c(MR_cCRE_MA_genes_extragenicLinked, nonMR_cCRE_MA_genes_extragenicLinked),
                                                     c(MR_cCRE_MA_genes_extragenicLinked_resamp, nonMR_cCRE_MA_genes_extragenicLinked_resamp), 
                                                     row.names=c("MR_cCRE", "NonMR_cCRE"))
  names(dat_MR_cCRE_MA_genes_extragenicLinked) = c("MA_gene", "Resampled_gene")
  test_MR_cCRE_MA_genes_extragenicLinked = fisher.test(dat_MR_cCRE_MA_genes_extragenicLinked)
  p_values_MR_cCRE_MA_genes_extragenicLinked[i] =  test_MR_cCRE_MA_genes_extragenicLinked$p.value * sign(log2(test_MR_cCRE_MA_genes_extragenicLinked$estimate))
  log2fc_MR_cCRE_MA_genes_extragenicLinked[i] = log2(test_MR_cCRE_MA_genes_extragenicLinked$estimate)
  
  #MA cCREs extragenic and Cicero-linked to MR genes
  MA_cCRE_MR_genes_extragenicLinked = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  MA_cCRE_MR_genes_extragenicLinked_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MR_genes_extragenicLinked = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MR_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MR_genes_extragenicLinked_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_mr_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MR_genes_extragenicLinked = data.frame(c(MA_cCRE_MR_genes_extragenicLinked, nonMA_cCRE_MR_genes_extragenicLinked),
                                                     c(MA_cCRE_MR_genes_extragenicLinked_resamp, nonMA_cCRE_MR_genes_extragenicLinked_resamp), 
                                                     row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MR_genes_extragenicLinked) = c("MR_gene", "Resampled_gene")
  test_MA_cCRE_MR_genes_extragenicLinked = fisher.test(dat_MA_cCRE_MR_genes_extragenicLinked)
  p_values_MA_cCRE_MR_genes_extragenicLinked[i] =  test_MA_cCRE_MR_genes_extragenicLinked$p.value * sign(log2(test_MA_cCRE_MR_genes_extragenicLinked$estimate))
  log2fc_MA_cCRE_MR_genes_extragenicLinked[i] = log2(test_MA_cCRE_MR_genes_extragenicLinked$estimate)
  #MA cCREs extragenic and Cicero-linked to MA genes
  MA_cCRE_MA_genes_extragenicLinked = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  MA_cCRE_MA_genes_extragenicLinked_resamp = length(intersect(union_ma_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  nonMA_cCRE_MA_genes_extragenicLinked = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% MA_genes[,V4]), cCRE_label]))
  nonMA_cCRE_MA_genes_extragenicLinked_resamp = length(intersect(union_nonMA_cCREs[, V4], mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[(Gene %in% Pv_ma_genes_PvTPM_MeCP2KO_df_resamp), cCRE_label]))
  dat_MA_cCRE_MA_genes_extragenicLinked = data.frame(c(MA_cCRE_MA_genes_extragenicLinked, nonMA_cCRE_MA_genes_extragenicLinked),
                                                     c(MA_cCRE_MA_genes_extragenicLinked_resamp, nonMA_cCRE_MA_genes_extragenicLinked_resamp), 
                                                     row.names=c("MA_cCRE", "NonMA_cCRE"))
  names(dat_MA_cCRE_MA_genes_extragenicLinked) = c("MA_gene", "Resampled_gene")
  test_MA_cCRE_MA_genes_extragenicLinked = fisher.test(dat_MA_cCRE_MA_genes_extragenicLinked)
  p_values_MA_cCRE_MA_genes_extragenicLinked[i] =  test_MA_cCRE_MA_genes_extragenicLinked$p.value * sign(log2(test_MA_cCRE_MA_genes_extragenicLinked$estimate))
  log2fc_MA_cCRE_MA_genes_extragenicLinked[i] = log2(test_MA_cCRE_MA_genes_extragenicLinked$estimate)
  
}

all_pvals_intragenic_union = rbind(p_values_MR_cCRE_MR_genes_intragenic,p_values_MR_cCRE_MA_genes_intragenic,p_values_MA_cCRE_MR_genes_intragenic,p_values_MA_cCRE_MA_genes_intragenic)
median_pvals_intragenic_union = apply(all_pvals_intragenic_union,1,median,na.rm=TRUE)
plot_pvals_intragenic_union = -matrix(log10(abs(median_pvals_intragenic_union)) * sign(median_pvals_intragenic_union),ncol=2,byrow=TRUE)
all_fcs_intragenic_union = rbind(log2fc_MR_cCRE_MR_genes_intragenic,log2fc_MR_cCRE_MA_genes_intragenic,log2fc_MA_cCRE_MR_genes_intragenic,log2fc_MA_cCRE_MA_genes_intragenic)
median_fcs_intragenic_union = apply(all_fcs_intragenic_union,1,median,na.rm=TRUE)
plot_fcs_intragenic_union = matrix(median_fcs_intragenic_union, ncol=2, byrow=2)


all_pvals_Cicero_union = rbind(p_values_MR_cCRE_MR_genes_Cicero,p_values_MR_cCRE_MA_genes_Cicero,p_values_MA_cCRE_MR_genes_Cicero,p_values_MA_cCRE_MA_genes_Cicero)
median_pvals_Cicero_union = apply(all_pvals_Cicero_union,1,median,na.rm=TRUE)
plot_pvals_Cicero_union = -matrix(log10(abs(median_pvals_Cicero_union)) * sign(median_pvals_Cicero_union),ncol=2,byrow=TRUE)
all_fcs_Cicero_union = rbind(log2fc_MR_cCRE_MR_genes_Cicero,log2fc_MR_cCRE_MA_genes_Cicero,log2fc_MA_cCRE_MR_genes_Cicero,log2fc_MA_cCRE_MA_genes_Cicero)
median_fcs_Cicero_union = apply(all_fcs_Cicero_union,1,median,na.rm=TRUE)
plot_fcs_Cicero_union = matrix(median_fcs_Cicero_union, ncol=2, byrow=2)

all_pvals_tad_union = rbind(p_values_MR_cCRE_MR_genes_sameTAD,p_values_MR_cCRE_MA_genes_sameTAD,p_values_MA_cCRE_MR_genes_sameTAD,p_values_MA_cCRE_MA_genes_sameTAD)
median_pvals_tad_union = apply(all_pvals_tad_union,1,median,na.rm=TRUE)
plot_pvals_tad_union = -matrix(log10(abs(median_pvals_tad_union)) * sign(median_pvals_tad_union),ncol=2,byrow=TRUE)
all_fcs_tad_union = rbind(log2fc_MR_cCRE_MR_genes_sameTAD,log2fc_MR_cCRE_MA_genes_sameTAD,log2fc_MA_cCRE_MR_genes_sameTAD,log2fc_MA_cCRE_MA_genes_sameTAD)
median_fcs_tad_union = apply(all_fcs_tad_union,1,median,na.rm=TRUE)
plot_fcs_tad_union = matrix(median_fcs_tad_union, ncol=2, byrow=2)


all_pvals_intragenicLinked_union = rbind(p_values_MR_cCRE_MR_genes_intragenicLinked,p_values_MR_cCRE_MA_genes_intragenicLinked,p_values_MA_cCRE_MR_genes_intragenicLinked,p_values_MA_cCRE_MA_genes_intragenicLinked)
median_pvals_intragenicLinked_union = apply(all_pvals_intragenicLinked_union,1,median,na.rm=TRUE)
plot_pvals_intragenicLinked_union = -matrix(log10(abs(median_pvals_intragenicLinked_union)) * sign(median_pvals_intragenicLinked_union),ncol=2,byrow=TRUE)
all_fcs_intragenicLinked_union = rbind(log2fc_MR_cCRE_MR_genes_intragenicLinked,log2fc_MR_cCRE_MA_genes_intragenicLinked,log2fc_MA_cCRE_MR_genes_intragenicLinked,log2fc_MA_cCRE_MA_genes_intragenicLinked)
median_fcs_intragenicLinked_union = apply(all_fcs_intragenicLinked_union,1,median,na.rm=TRUE)
plot_fcs_intragenicLinked_union = matrix(median_fcs_intragenicLinked_union, ncol=2, byrow=2)

all_pvals_extragenicLinked_union = rbind(p_values_MR_cCRE_MR_genes_extragenicLinked,p_values_MR_cCRE_MA_genes_extragenicLinked,p_values_MA_cCRE_MR_genes_extragenicLinked,p_values_MA_cCRE_MA_genes_extragenicLinked)
median_pvals_extragenicLinked_union = apply(all_pvals_extragenicLinked_union,1,median,na.rm=TRUE)
plot_pvals_extragenicLinked_union = -matrix(log10(abs(median_pvals_extragenicLinked_union)) * sign(median_pvals_extragenicLinked_union),ncol=2,byrow=TRUE)
all_fcs_extragenicLinked_union = rbind(log2fc_MR_cCRE_MR_genes_extragenicLinked,log2fc_MR_cCRE_MA_genes_extragenicLinked,log2fc_MA_cCRE_MR_genes_extragenicLinked,log2fc_MA_cCRE_MA_genes_extragenicLinked)
median_fcs_extragenicLinked_union = apply(all_fcs_extragenicLinked_union,1,median,na.rm=TRUE)
plot_fcs_extragenicLinked_union = matrix(median_fcs_extragenicLinked_union, ncol=2, byrow=2)

total_plot_pvals_union = cbind(plot_pvals_intragenic_union, plot_pvals_Cicero_union, plot_pvals_tad_union, plot_pvals_intragenicLinked_union)
total_plot_fcs_union = cbind(plot_fcs_intragenic_union, plot_fcs_Cicero_union, plot_fcs_tad_union, plot_fcs_intragenicLinked_union)

union_cCRE_palette <- colorRampPalette(c("orchid", "white", "goldenrod2"))(n = 299)
union_cCRE_col_breaks = c(seq(-34,-3.1,length=100),  # for orchid
                       seq(-3,3,length=100),           # for white
                       seq(3.1,34,length=100))             # for goldenrod2
setEPS()
postscript("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_cicero_sameTAD_cognateLinked_heatmap.eps")
heatmap.2(total_plot_pvals_union,
          notecol="black",      # change font color of cell labels to black
          breaks=union_cCRE_col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=union_cCRE_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA", key=TRUE)            # turn off column clustering
dev.off()


total_plot_pvals_union_noTAD = cbind(plot_pvals_intragenic_union, plot_pvals_Cicero_union, plot_pvals_intragenicLinked_union)
total_plot_fcs_union_noTAD = cbind(plot_fcs_intragenic_union, plot_fcs_Cicero_union, plot_fcs_intragenicLinked_union)

png("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_cicero_cognateLinked_heatmap.png", width=6, height=4.3, units="in", res=72)
heatmap.2(total_plot_pvals_union_noTAD,
          notecol="black",      # change font color of cell labels to black
          breaks=union_cCRE_col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=union_cCRE_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA", key=TRUE) 
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_cicero_cognateLinked_heatmap.eps", width=6, height=4.3)
heatmap.2(total_plot_pvals_union_noTAD,
          notecol="black",      # change font color of cell labels to black
          breaks=union_cCRE_col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=union_cCRE_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA", key=TRUE) 
dev.off()


#Include extragenic linked
total_plot_pvals_union_noTAD_withExtra = cbind(plot_pvals_intragenic_union, plot_pvals_Cicero_union, plot_pvals_intragenicLinked_union, plot_pvals_extragenicLinked_union)
total_plot_fcs_union_noTAD_withExtra = cbind(plot_fcs_intragenic_union, plot_fcs_Cicero_union, plot_fcs_intragenicLinked_union, plot_pvals_extragenicLinked_union)


png("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_extragenic_cicero_Linked_heatmap.png", width=6, height=4.3, units="in", res=72)
heatmap.2(total_plot_pvals_union_noTAD_withExtra,
          notecol="black",      # change font color of cell labels to black
          breaks=union_cCRE_col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=union_cCRE_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA", key=TRUE) 
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_extragenic_cicero_Linked_heatmap.eps", width=6, height=4.3)
heatmap.2(total_plot_pvals_union_noTAD_withExtra,
          notecol="black",      # change font color of cell labels to black
          breaks=union_cCRE_col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=union_cCRE_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA", key=TRUE) 
dev.off()

round(total_plot_fcs_union_noTAD_withExtra, 1)

write.table(total_plot_pvals_union_noTAD_withExtra, "HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_extragenic_cicero_Linked_heatmap_log10pvals.txt", quote=F, row.names=F, col.names=F)
write.table(total_plot_fcs_union_noTAD_withExtra, "HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/cCRE_gene_association_enrichment_values/cell_confusion_paper_images/union_cCRE_MeCP2dys_intragenic_extragenic_cicero_Linked_heatmap_log2OddsRatios.txt", quote=F, row.names=F, col.names=F)
