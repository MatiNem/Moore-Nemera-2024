library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
coding_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")
pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")

mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")

mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(cCRE_label %in% nonPromoter_cCREs_mm9[,V4]) & (Intragenic==1) & (Intragenic_to_linked_gene==1), ]
PVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(cCRE_label %in% nonPromoter_cCREs_mm9[,V4]) & (Intragenic==1) & (Intragenic_to_linked_gene==1), ]
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(cCRE_label %in% nonPromoter_cCREs_mm9[,V4]) & (Intragenic==1) & (Intragenic_to_linked_gene==1), ]


PVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9 = unique(PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene)==1, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)])
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9 = unique(nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene)==1, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)])

write.table(PVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9, "HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9, "HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")

#significantly dysregulated union nonpromoter cCREs
union_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
union_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
union_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#significantly dysregulated PVGA nonpromoter cCREs
Pv_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
Pv_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
Pv_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#significantly dysregulated nonPVGA nonpromoter cCREs
nonPv_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
nonPv_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
nonPv_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#intragenic non-promoter cCREs
intragenic_nonPromoter_cCREs_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_ensgene_mm9.txt")

#Isolate the connection between cCREs and the genes in which they reside
cCRE_intragenic_connections = intragenic_nonPromoter_cCREs_mm9[, paste0(cCRE_label,"|",Gene)]
intragenic_nonPromoter_cCREs_mm9[, Conns:=cCRE_intragenic_connections]

nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_nonPromoter_cCREs_1500bpWindows_Pv_MeCP2KO_H3K27ac_ChIP_edgeR.txt")
nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_nonPromoter_cCREs_1500bpWindows_Pv_MeCP2TG_H3K27ac_ChIP_edgeR.txt")


PvMR_metaMR_genes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
PvUnchanged_metaMR_genes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")

coding_genes_mm9_with_UnioncCRE_numbers = coding_genes_mm9[, .(V1, V2, V3, V4)]
names(coding_genes_mm9_with_UnioncCRE_numbers) = c("gene_chrom", "gene_start", "gene_end", "Gene")

for(i in coding_genes_mm9_with_UnioncCRE_numbers[, Gene]){
  union_intragenic_cCREs = unique(intragenic_nonPromoter_cCREs_mm9[(Gene == i), cCRE_label])
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cCREs := length(union_intragenic_cCREs)]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cCREs_sigMR := length(intersect(union_intragenic_cCREs, union_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cCREs_sigMA := length(intersect(union_intragenic_cCREs, union_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_intragenic := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
            nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_union_intragenic_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_intragenic := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                        nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_union_intragenic_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_intragenic := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_intragenic := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cCREs, as.numeric(logFC)], na.rm=TRUE)]
}

mousebrain_union_nonPromoter_cCREs_genes_genicBooleans <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
for(i in coding_genes_mm9_with_UnioncCRE_numbers[, Gene]){
  #cognate-linked
  union_intragenic_cognateLinked_cCREs = unique(mousebrain_union_nonPromoter_cCREs_genes_genicBooleans[(Gene == i) & (Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label])
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs := length(union_intragenic_cognateLinked_cCREs)]         
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs_sigMR := length(intersect(union_intragenic_cognateLinked_cCREs, union_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs_sigMA := length(intersect(union_intragenic_cognateLinked_cCREs, union_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
                                             nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_union_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                             nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_union_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  
  #extragenic
  union_extragenic_cCREs = unique(mousebrain_union_nonPromoter_cCREs_genes_genicBooleans[(Gene == i) & (Intragenic==0), cCRE_label])
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_extragenic_cCREs := length(union_extragenic_cCREs)]         
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_extragenic_cCREs_sigMR := length(intersect(union_extragenic_cCREs, union_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_extragenic_cCREs_sigMA := length(intersect(union_extragenic_cCREs, union_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_extragenic := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_extragenic_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
                                                                      nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_extragenic_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_union_extragenic_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_extragenic := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_extragenic_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                                                      nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_extragenic_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_union_extragenic_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_extragenic := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% union_extragenic_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_extragenic := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% union_extragenic_cCREs, as.numeric(logFC)], na.rm=TRUE)]
}

#write.table(coding_genes_mm9_with_UnioncCRE_numbers, file="HG_lab/Mati/GabelLab/genesets/coding_genes_union_nonpromoter_cCRE_numbers_H3K27ac_data_table.txt", quote=F, row.names=F, sep="\t")
coding_genes_mm9_with_UnioncCRE_numbers = fread("HG_lab/Mati/GabelLab/genesets/coding_genes_union_nonpromoter_cCRE_numbers_H3K27ac_data_table.txt")


for(i in coding_genes_mm9_with_UnioncCRE_numbers[, Gene]){
  #cognate-linked
  union_intragenic_cognateLinked_cCREs = unique(mousebrain_union_nonPromoter_cCREs_genes_genicBooleans[(Gene == i) & (Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label])
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs := length(union_intragenic_cognateLinked_cCREs)]         
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs_sigMR := length(intersect(union_intragenic_cognateLinked_cCREs, union_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), num_union_intragenic_cognateLinked_cCREs_sigMA := length(intersect(union_intragenic_cognateLinked_cCREs, union_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
                                                                                                                            nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_union_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                                                                                                            nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% union_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_union_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_UnioncCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% union_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
}  

#gene lists
pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]) & (num_union_intragenic_cognateLinked_cCREs_sigMR >= 4) & (mean_H3K27ac_KO_fc_intragenic_cognateLinked > 0) & (mean_H3K27ac_TG_fc_intragenic_cognateLinked < 0), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (num_union_intragenic_cognateLinked_cCREs_sigMR >= 1) & (mean_H3K27ac_KO_fc_intragenic_cognateLinked > 0) & (mean_H3K27ac_TG_fc_intragenic_cognateLinked < 0), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (mean_H3K27ac_KO_fc_intragenic_cognateLinked > 0) & (mean_H3K27ac_TG_fc_intragenic_cognateLinked < 0), ])

View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (Gene %in% PvUnchanged_metaMR_genes[,V4]), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (Gene %in% PvUnchanged_metaMR_genes[,V4]) & (mean_H3K27ac_KO_fc_intragenic_cognateLinked > 0) & (mean_H3K27ac_TG_fc_intragenic_cognateLinked < 0), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (Gene %in% PvUnchanged_metaMR_genes[,V4]) & (abs(mean_H3K27ac_KO_fc_intragenic_cognateLinked) < 0.15) & (abs(mean_H3K27ac_TG_fc_intragenic_cognateLinked) < 0.15), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (Gene %in% PvUnchanged_metaMR_genes[,V4]), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (abs(mean_H3K27ac_KO_fc_intragenic_cognateLinked) <= 0.1) & (abs(mean_H3K27ac_TG_fc_intragenic_cognateLinked) <= 0.1), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (sign(mean_H3K27ac_KO_fc_intragenic_cognateLinked) == sign(mean_H3K27ac_TG_fc_intragenic_cognateLinked)), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]) & (abs(mean_H3K27ac_KO_fc_intragenic_cognateLinked) <= 0.15) & (abs(mean_H3K27ac_TG_fc_intragenic_cognateLinked) <= 0.15), ])
View(coding_genes_mm9_with_UnioncCRE_numbers[(Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4]), ])



gene_of_interest = "Msi2"
cCREs_of_interest = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]

nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest), as.numeric(logFC)]


ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="indianred4")+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity")+
  #coord_cartesian(ylim=c(0,0.18))+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

H3K27ac_cCRE_in_gene_func <- function(gene_of_interest, cCRE_gene_connections, title, directory){
  cCREs_of_interest = cCRE_gene_connections[(Gene == gene_of_interest), cCRE_label]
  ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
    ggtitle(title)+
    geom_bar(stat="identity", fill="maroon")+
    ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
    coord_cartesian(ylim=c(-2, 2))+
    theme_bw()+
    theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
  ggsave(filename = paste0("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/", directory, "/", gene_of_interest, "_PvMeCP2KO_H3K27ac_log2fc_barplot.png"),   width = 4, height = 5, dpi = 300, units = "in", device='png')
  ggsave(filename = paste0("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/", directory, "/", gene_of_interest, "_PvMeCP2KO_H3K27ac_log2fc_barplot.eps"),   width = 4, height = 5, dpi = 300, units = "in", device='eps')
  
  ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
    ggtitle(title)+
    geom_bar(stat="identity", fill="maroon")+
    ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
    coord_cartesian(ylim=c(-2, 2))+
    theme_bw()+
    theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
  ggsave(filename = paste0("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/", directory, "/", gene_of_interest, "_PvMeCP2TG_H3K27ac_log2fc_barplot.png"),   width = 4, height = 5, dpi = 300, units = "in", device='png')
  ggsave(filename = paste0("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/", directory, "/", gene_of_interest, "_PvMeCP2TG_H3K27ac_log2fc_barplot.eps"),   width = 4, height = 5, dpi = 300, units = "in", device='eps')
}

H3K27ac_cCRE_in_gene_func("Ptprg", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Ptprg, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Rassf8", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Rassf8, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Dscam", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Dscam, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Galnt14", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Galnt14, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Negr1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Negr1, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Cdh13", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Cdh13, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Esrrg", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Esrrg, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Fat1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Fat1, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Myo1b", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Myo1b, cognate-linked cCREs", directory="pv_mr_genes")
H3K27ac_cCRE_in_gene_func("Srgap1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Srgap1, cognate-linked cCREs", directory="pv_mr_genes")

H3K27ac_cCRE_in_gene_func("Kctd1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Kctd1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Asap1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Asap1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Zbtb20", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Zbtb20, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Foxo1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Foxo1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Rora", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Rora, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Dnm3", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Dnm3, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Ust", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Ust, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Fmn1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Fmn1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Lin7a", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Lin7a, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Plekha7", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Plekha7, cognate-linked cCREs", directory="otherCellType_mr_genes")

H3K27ac_cCRE_in_gene_func("Slc24a4", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Slc24a4, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Arhgef3", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Arhgef3, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Atp1a1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Atp1a1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Phactr1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Phactr1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Pde7b", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Pde7b, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Camk2a", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Camk2a, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Syt13", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Syt13, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Frmd5", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Frmd5, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Man1c1", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Man1c1, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Cntn5", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Cntn5, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Dgkg", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Dgkg, cognate-linked cCREs", directory="otherCellType_mr_genes")
H3K27ac_cCRE_in_gene_func("Slc8a3", cCRE_gene_connections=mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans, title="Slc8a3, cognate-linked cCREs", directory="otherCellType_mr_genes")



nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table =  data.table(rbind(cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[, .(cCRE_label, logFC)], mutant="KO/WT"),
                                                                                cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[, .(cCRE_label, logFC)], mutant="OE/WT")))
nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[, c("cCRE_label_number") := tstrsplit(cCRE_label, "cCREs", fixed=TRUE, keep=2L)]

nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[, cCRE_label_with_mutant := paste(cCRE_label, mutant)]

nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table = nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[order(as.numeric(cCRE_label_number))]
cCRE_label_with_mutant_order = nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[, cCRE_label_with_mutant]
cCRE_label_order = unique(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[, cCRE_label])

nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table = nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table %>% mutate(cCRE_label_with_mutant = factor(cCRE_label_with_mutant, levels=cCRE_label_with_mutant_order))
nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table = nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table %>% mutate(cCRE_label = factor(cCRE_label, levels=cCRE_label_order))
nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table = nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))


gene_of_interest = "Srgap1"
cCREs_of_interest = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label_with_mutant, y=as.numeric(logFC), fill=mutant))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("KO/WT"="seagreen3", "OE/WT"="seagreen1"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 mutant/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label_with_mutant, y=as.numeric(logFC), fill=mutant))+
  geom_col(width = 0.5)+
  ggtitle("Srgap1, cognate-linked cCREs")+
  scale_fill_manual(values=c("KO/WT"="black", "OE/WT"="gray"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 mutant/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  facet_grid(.~cCRE_label_number,
             switch="x") + 
  theme_bw()+
  theme(panel.spacing = unit(0, "in"), plot.title=element_text(hjust=0.5), strip.text.x = element_text(angle=90), legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Srgap1_PvMeCP2mut_H3K27ac_log2fc_barplot.png", width = 30, height = 5, dpi = 300, units = "in", device='png')

gene_of_interest = "Lin7a"
cCREs_of_interest = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label_with_mutant, y=as.numeric(logFC), fill=mutant))+
  geom_bar(stat="identity")+
  ggtitle("Lin7a, cognate-linked cCREs")+
  scale_fill_manual(values=c("KO/WT"="seagreen3", "OE/WT"="seagreen1"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 mutant/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  facet_grid(.~cCRE_label) + 
  theme_bw()+
  theme(panel.spacing = unit(0.001, "cm"), plot.title=element_text(hjust=0.5), strip.text.x = element_blank(), legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

mousebrain_union_cCREs_PV_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_PV_mCA[, cCRE_methylation := V5/V6]

avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

mousebrain_union_cCREs_PV_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCA$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]

#turning negative methylation values into zeros
mousebrain_union_cCREs_PV_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]






gene_of_interest = "Srgap1"
cCREs_of_interest = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(mousebrain_union_cCREs_mPv_mCA[(V4 %in% cCREs_of_interest)], aes(x=V4, y=as.numeric(cCRE_methylation)))+
  geom_col(width = 0.5)+
  ggtitle(paste0(gene_of_interest, ", cognate-linked cCREs"))+
  ylab("mPv mCA/CA") + xlab("")+
  coord_cartesian(ylim=c(0, 0.15))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/enhancer_meth_plots/Srgap1_cognateLinked_nonPromoter_cCREs_mPv_mCAperCA_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_meth_plots/Srgap1_cognateLinked_nonPromoter_cCREs_mPv_mCAperCA_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

gene_of_interest = "Syt13"
cCREs_of_interest = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(mousebrain_union_cCREs_mPv_mCA[(V4 %in% cCREs_of_interest)], aes(x=V4, y=as.numeric(cCRE_methylation)))+
  geom_col(width = 0.5)+
  ggtitle(paste0(gene_of_interest, ", cognate-linked cCREs"))+
  ylab("mPv mCA/CA") + xlab("")+
  coord_cartesian(ylim=c(0, 0.15))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/enhancer_meth_plots/Syt13_cognateLinked_nonPromoter_cCREs_mPv_mCAperCA_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_meth_plots/Syt13_cognateLinked_nonPromoter_cCREs_mPv_mCAperCA_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

#PVGA and nonPVGA cCREs
coding_genes_mm9_with_PVGAcCRE_numbers = coding_genes_mm9[, .(V1, V2, V3, V4)]
names(coding_genes_mm9_with_PVGAcCRE_numbers) = c("gene_chrom", "gene_start", "gene_end", "Gene")

coding_genes_mm9_with_nonPVGAcCRE_numbers = coding_genes_mm9[, .(V1, V2, V3, V4)]
names(coding_genes_mm9_with_nonPVGAcCRE_numbers) = c("gene_chrom", "gene_start", "gene_end", "Gene")


for(i in coding_genes_mm9[, V4]){
  #cognate-linked PVGA cCREs
  Pv_intragenic_cognateLinked_cCREs = unique(PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(Gene == i) & (Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label])
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), num_Pv_intragenic_cognateLinked_cCREs := length(Pv_intragenic_cognateLinked_cCREs)]         
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), num_Pv_intragenic_cognateLinked_cCREs_sigMR := length(intersect(Pv_intragenic_cognateLinked_cCREs, Pv_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), num_Pv_intragenic_cognateLinked_cCREs_sigMA := length(intersect(Pv_intragenic_cognateLinked_cCREs, Pv_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
                                                                                                                            nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_Pv_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                                                                                                            nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_Pv_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% Pv_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_PVGAcCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% Pv_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  
  #cognate-linked nonPVGA cCREs
  nonPv_intragenic_cognateLinked_cCREs = unique(nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[(Gene == i) & (Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label])
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), num_nonPv_intragenic_cognateLinked_cCREs := length(nonPv_intragenic_cognateLinked_cCREs)]         
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), num_nonPv_intragenic_cognateLinked_cCREs_sigMR := length(intersect(nonPv_intragenic_cognateLinked_cCREs, nonPv_mr_cCREs_sig[,V4]))]
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), num_nonPv_intragenic_cognateLinked_cCREs_sigMA := length(intersect(nonPv_intragenic_cognateLinked_cCREs, nonPv_ma_cCREs_sig[,V4]))]
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), H3K27ac_KOandTG_MR_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label], 
                                                                                                                           nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label]))/num_nonPv_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), H3K27ac_KOandTG_MA_prop_intragenic_cognateLinked := length(intersect(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)<0), cCRE_label], 
                                                                                                                           nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs) & (as.numeric(logFC)>0), cCRE_label]))/num_nonPv_intragenic_cognateLinked_cCREs]
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), mean_H3K27ac_KO_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
  coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene == i), mean_H3K27ac_TG_fc_intragenic_cognateLinked := mean(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% nonPv_intragenic_cognateLinked_cCREs, as.numeric(logFC)], na.rm=TRUE)]
}  

#write.table(coding_genes_mm9_with_PVGAcCRE_numbers, file="HG_lab/Mati/GabelLab/genesets/coding_genes_mm9_with_PVGAcCRE_numbers_H3K27ac_data_table.txt", quote=F, row.names=F, sep="\t")
#write.table(coding_genes_mm9_with_nonPVGAcCRE_numbers, file="HG_lab/Mati/GabelLab/genesets/coding_genes_mm9_with_nonPVGAcCRE_numbers_H3K27ac_data_table.txt", quote=F, row.names=F, sep="\t")

View(coding_genes_mm9_with_PVGAcCRE_numbers[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]), ])

View(coding_genes_mm9_with_nonPVGAcCRE_numbers[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]) & (num_nonPv_intragenic_cognateLinked_cCREs_sigMR >= 3) & (mean_H3K27ac_KO_fc_intragenic_cognateLinked >= 0.4), ])


gene_of_interest = "Slit1"
cCREs_of_interest = nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", non-Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", non-Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

gene_of_interest = "Galnt14"
cCREs_of_interest = nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", non-Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

gene_of_interest = "Galnt14"
cCREs_of_interest = PVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())


gene_of_interest = "Ptprg"
cCREs_of_interest = nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", non-Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_nonPVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_nonPVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

gene_of_interest = "Ptprg"
cCREs_of_interest = nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", non-Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_nonPVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_nonPVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

gene_of_interest = "Ptprg"
cCREs_of_interest = PVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

gene_of_interest = "Ptprg"
cCREs_of_interest = PVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene %in% gene_of_interest), cCRE_label]
ggplot(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% cCREs_of_interest)], aes(x=cCRE_label, y=as.numeric(logFC)))+
  geom_bar(stat="identity", fill="maroon")+
  ggtitle(paste0(gene_of_interest, ", Pv cognate-linked cCREs"))+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

PVGA_cCREs_Ptprg = PVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene == "Ptprg"),cCRE_label]
nonPVGA_cCREs_Ptprg = nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene == "Ptprg"),cCRE_label]

PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt = data.table(rbind(cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% PVGA_cCREs_Ptprg], cCRE_group="PV"),
                                                                            cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[cCRE_label %in% nonPVGA_cCREs_Ptprg], cCRE_group="Non-PV")))
PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt = PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt %>% mutate(cCRE_label = factor(cCRE_label, levels=c(c(PVGA_cCREs_Ptprg, nonPVGA_cCREs_Ptprg))))
PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt = PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV","Non-PV")))

ggplot(PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2KO_dt, aes(x=cCRE_label, y=as.numeric(logFC), fill=cCRE_group))+
  scale_fill_manual(name="Intragenic linked cCREs:", values=c("PV"="forestgreen", "Non-PV"="gray"))+
  geom_bar(stat="identity")+
  ggtitle("Ptprg")+
  ylab("Log2 fold change H3K27ac\n MeCP2 KO/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_and_nonPVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_and_nonPVGA_cognateLinked_cCREs_PvMeCP2KO_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt = data.table(rbind(cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% PVGA_cCREs_Ptprg], cCRE_group="PV"),
                                                                                      cbind(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[cCRE_label %in% nonPVGA_cCREs_Ptprg], cCRE_group="Non-PV")))
PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt = PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt %>% mutate(cCRE_label = factor(cCRE_label, levels=c(c(PVGA_cCREs_Ptprg, nonPVGA_cCREs_Ptprg))))
PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt = PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV","Non-PV")))

ggplot(PVGA_and_nonPVGA_cCREs_intragenicLinked_Ptprg_H3K27ac_PvMeCP2TG_dt, aes(x=cCRE_label, y=as.numeric(logFC), fill=cCRE_group))+
  scale_fill_manual(name="Intragenic linked cCREs:", values=c("PV"="forestgreen", "Non-PV"="gray"))+
  geom_bar(stat="identity")+
  ggtitle("Ptprg")+
  ylab("Log2 fold change H3K27ac\n MeCP2 OE/WT") + xlab("")+
  coord_cartesian(ylim=c(-2, 2))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_and_nonPVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cCRE_in_gene_H3K27ac_fc_barplots/pv_mr_genes/Ptprg_PVGA_and_nonPVGA_cognateLinked_cCREs_PvMeCP2TG_H3K27ac_log2fc_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##re-doing barplots as heatmaps
#adding genotype info
nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table2 <- copy(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table)
nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table2$condition <- "KO/WT"

nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table2 <- copy(nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table)
nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table2$condition <- "OE/WT"

nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table2 <- rbind(nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table2[, .(cCRE_label, logFC, condition)], nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table2[, .(cCRE_label, logFC, condition)])

Srgap1_cCREs = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene == "Srgap1"), cCRE_label]
Syt13_cCREs = mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[(Gene == "Syt13"), cCRE_label]

#reshape the data
Srgap1_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide <- dcast(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table2[cCRE_label %in% Srgap1_cCREs],
                   condition ~ cCRE_label,
                   value.var = "logFC",
                   drop = FALSE) %>% data.frame(row.names="condition") %>% as.matrix

Syt13_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide <- dcast(nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table2[cCRE_label %in% Syt13_cCREs],
                                                                                 condition ~ cCRE_label,
                                                                                 value.var = "logFC",
                                                                                 drop = FALSE) %>% data.frame(row.names="condition") %>% as.matrix



H3K27ac_col_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
H3K27ac_col_breaks = c(seq(-2.2, -0.5,length=100),  # for blue
               seq(-0.49,0.49,length=100), #for white
               seq(0.5, 2.2, length=100)) #for red

setEPS()
postscript("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Srgap1_PvMeCP2mut_H3K27ac_log2fc_heatmap.eps")
heatmap.2(Srgap1_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide,
          main = "Srgap1, cognate-linked cCREs", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=H3K27ac_col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=H3K27ac_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.8)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Syt13_PvMeCP2mut_H3K27ac_log2fc_heatmap.eps")
heatmap.2(Syt13_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide,
          main = "Syt13, cognate-linked cCREs", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=H3K27ac_col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=H3K27ac_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.8)
dev.off()

Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA <- rbind(
  cbind(id="row1", mousebrain_union_cCREs_PV_mCA[V4 %in% Srgap1_cCREs], cCRE_list="Srgap1"),
  cbind(id="row1", mousebrain_union_cCREs_PV_mCA[V4 %in% Syt13_cCREs], cCRE_list="Syt13"),
  cbind(id="row2", mousebrain_union_cCREs_PV_mCA[V4 %in% Srgap1_cCREs], cCRE_list="Srgap1"),
  cbind(id="row2", mousebrain_union_cCREs_PV_mCA[V4 %in% Syt13_cCREs], cCRE_list="Syt13")
)

Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA$cCRE_label_new <- Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA[, paste0(cCRE_list,"_", V4)]

Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA[(id=="row1") & (cCRE_list=="Syt13"), cCRE_methylation := 0]
Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA[(id=="row2") & (cCRE_list=="Srgap1"), cCRE_methylation := 0]

Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA_dcast <- dcast(Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA, 
                                               formula = id ~ cCRE_label_new, 
                                               value.var = "cCRE_methylation")

Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA_matrix <- as.matrix(Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA_dcast, rownames="id")


#methylation heatmaps
PV_mCA_col_palette <- colorRampPalette(c("white", "seagreen3"))(n = 199)
#mPv_mCA_col_breaks = c(seq(0, ,length=100),  # for white
#                       seq(-0.49,0.49,length=100), #for white
#                       seq(0.5, 2.2, length=100)) #for red

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/enhancer_heatmaps/Srgap1_Syt13_cognateLinked_nonPromoter_cCREs_PV_WT_KO_mCAperCA_heatmap.eps")
heatmap.2(Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA_matrix,
          main = "Srgap1 and Syt13", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=H3K27ac_col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=PV_mCA_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.4)
dev.off()

##accentuating differences in heatmap more
##H3K27ac

H3K27ac_col_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
H3K27ac_col_breaks2 = c(seq(-2.2, -1,length=100),  # for blue
                       seq(-0.99,0.99,length=100), #for white
                       seq(1, 2.2, length=100)) #for red

H3K27ac_col_breaks3 = c(seq(-2.2, -1.8,length=100),  # for blue
                        seq(-1.79,1.79,length=100), #for white
                        seq(1.8, 2.2, length=100)) #for red

setEPS()
postscript("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Srgap1_PvMeCP2mut_H3K27ac_log2fc_heatmap_lighter.eps")
heatmap.2(Srgap1_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide,
          main = "Srgap1, cognate-linked cCREs", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=H3K27ac_col_breaks2,
          #breaks=H3K27ac_col_breaks3,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=H3K27ac_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.8)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Syt13_PvMeCP2mut_H3K27ac_log2fc_heatmap_lighter.eps")
heatmap.2(Syt13_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_new_table_wide,
          main = "Syt13, cognate-linked cCREs", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=H3K27ac_col_breaks2,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=H3K27ac_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.8)
dev.off()

#
PV_mCA_col_palette <- colorRampPalette(c("white", "seagreen3"))(n = 199)
PV_mCA_col_breaks = c(seq(0, 0.10,length=100),  # for white
                        seq(0.11,0.13,length=100)) #for seagreen3


setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/enhancer_heatmaps/Srgap1_Syt13_cognateLinked_nonPromoter_cCREs_PV_WT_KO_mCAperCA_heatmap_lighter.eps")
heatmap.2(Srgap1_Syt13_mousebrain_union_cCREs_PV_mCA_matrix,
          main = "Srgap1 and Syt13", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=PV_mCA_col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=PV_mCA_col_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE,
          cexRow = 0.8,             # decrease size of row labels
          cexCol = 0.4)
dev.off()

