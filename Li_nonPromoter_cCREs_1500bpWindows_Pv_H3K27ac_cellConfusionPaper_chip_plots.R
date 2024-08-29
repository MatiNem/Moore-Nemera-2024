library(edgeR)
library(data.table)
library(dplyr)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("MADAM")
library(MADAM)

nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")

#All PVGA cCREs
PVGA_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/PVGA/PVGA_cCRE_mm9.bed")
#All PVGA non-promoter cCREs
PVGA_nonPromoter_cCREs_mm9 = nonPromoter_cCREs_mm9[V4 %in% PVGA_cCREs_mm9[,V4], ]
#All nonPVGA nonpromoter cCREs
nonPVGA_nonPromoter_cCREs_mm9 = nonPromoter_cCREs_mm9[!(V4 %in% PVGA_cCREs_mm9[,V4]), ]

nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_nonPromoter_cCREs_1500bpWindows_Pv_MeCP2KO_H3K27ac_ChIP_edgeR.txt")
nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_nonPromoter_cCREs_1500bpWindows_Pv_MeCP2TG_H3K27ac_ChIP_edgeR.txt")

#union cCREs linked to non-deduplicated Pv MeCP2-repressed genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")

#PVGA cCREs linked to non-deduplicated Pv MeCP2-repressed genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")


#Non-PVGA cCREs linked to non-deduplicated Pv MeCP2-repressed genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

##intragenic cCREs linked to their cognate genes
#union
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

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
#union
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

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

##TG and KO together, intragenic-linked cCREs linked to unchanged and MR genes only 
#Make a data table containing union cCREs linked to MeCP2-regulated genes separately
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv unchanged genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv MeCP2-repressed genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv unchanged genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Pv unchanged genes", "Pv MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv unchanged genes"="gray", "Pv MeCP2-repressed genes"="red")) +
  facet_grid(.~mutant,
             labeller = label_value,
             switch = "x"
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value

#extragenic linked only
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv unchanged genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv MeCP2-repressed genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv unchanged genes"),
                                                                                                              cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Pv unchanged genes", "Pv MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv unchanged genes"="gray", "Pv MeCP2-repressed genes"="red")) +
  facet_grid(.~mutant,
             labeller = label_value,
             switch = "x"
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.851781e-07, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=8.801089e-06, ****


##including both intragenic and extragenic
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv unchanged genes", genomic_loc="Intragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv unchanged genes", genomic_loc="Intragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv unchanged genes", genomic_loc="Extragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Pv MeCP2-repressed genes", genomic_loc="Extragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv unchanged genes", genomic_loc="Extragenic"),
                                                                                                             cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Pv MeCP2-repressed genes", genomic_loc="Extragenic")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra %>% mutate(geneset = factor(geneset, levels=c("Pv unchanged genes", "Pv MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra %>% mutate(genomic_loc = factor(genomic_loc, levels=c("Intragenic", "Extragenic")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv unchanged genes"="gray", "Pv MeCP2-repressed genes"="red")) +
  facet_grid(.~mutant+genomic_loc,
             labeller = label_value,
             switch = "x"
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Li2021_union_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.778946e-163, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.851781e-07, ****


wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=2.686386e-156, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Pv unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Pv MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=8.801089e-06, ****

##TG
#Make a data table containing the Pv and non-Pv cCREs linked to MeCP2-regulated genes separately
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Other-cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Other-cell-type MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MeCP2-repressed genes", "Other-cell-type MeCP2-repressed genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Cell-type unchanged genes"="gray", "Cell-type MeCP2-repressed genes"="red", "Other-cell-type MeCP2-repressed genes"="lightpink" ,"Cell-type MeCP2-activated genes"="blue", "Other-cell-type MeCP2-activated genes"="cyan3")) +
  facet_grid(.~cCRE_group,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (OE/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2TG_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2TG_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value


##KO
#Make a data table containing the Pv and non-Pv cCREs linked to MeCP2-regulated genes separately
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Pv cCREs", geneset="Other-cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], cCRE_group="Non-Pv cCREs", geneset="Other-cell-type MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MeCP2-repressed genes", "Other-cell-type MeCP2-repressed genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Cell-type unchanged genes"="gray", "Cell-type MeCP2-repressed genes"="red", "Other-cell-type MeCP2-repressed genes"="lightpink" ,"Cell-type MeCP2-activated genes"="blue", "Other-cell-type MeCP2-activated genes"="cyan3")) +
  facet_grid(.~cCRE_group,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (KO/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2KO_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2KO_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(cCRE_group=="Non-Pv cCREs") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value


##
##union nonPromoter cognate-linked cCRE H3K27ac ChIP boxplot, mutant
#Make a data table containing the cCREs linked to PV MR and other-cell-type MR genes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Other-cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type unchanged genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                                  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Other-cell-type MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MeCP2-repressed genes", "Other-cell-type MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Cell-type unchanged genes"="gray", "Cell-type MeCP2-repressed genes"="red", "Other-cell-type MeCP2-repressed genes"="lightpink" ,"Cell-type MeCP2-activated genes"="blue", "Other-cell-type MeCP2-activated genes"="cyan3")) +
  facet_grid(.~mutant,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value

#extragenic linked only
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type unchanged genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Other-cell-type MeCP2-repressed genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type unchanged genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type MeCP2-repressed genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Other-cell-type MeCP2-repressed genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MeCP2-repressed genes", "Other-cell-type MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Cell-type unchanged genes"="gray", "Cell-type MeCP2-repressed genes"="red", "Other-cell-type MeCP2-repressed genes"="lightpink" ,"Cell-type MeCP2-activated genes"="blue", "Other-cell-type MeCP2-activated genes"="cyan3")) +
  facet_grid(.~mutant,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.851781e-07, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.5937957, ns
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.01653431, *

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=8.801089e-06, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.1397587, ns
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.3523627, ns

##include intragenic and extragenic
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type unchanged genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Other-cell-type MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type unchanged genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Other-cell-type MeCP2-repressed genes", genomic_loc="Intragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type unchanged genes", genomic_loc="Extragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Cell-type MeCP2-repressed genes", genomic_loc="Extragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", geneset="Other-cell-type MeCP2-repressed genes", genomic_loc="Extragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type unchanged genes", genomic_loc="Extragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Cell-type MeCP2-repressed genes", genomic_loc="Extragenic"),
                                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", geneset="Other-cell-type MeCP2-repressed genes", genomic_loc="Extragenic")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MeCP2-repressed genes", "Other-cell-type MeCP2-repressed genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra %>% mutate(genomic_loc = factor(genomic_loc, levels=c("Intragenic", "Extragenic")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Cell-type unchanged genes"="gray", "Cell-type MeCP2-repressed genes"="red", "Other-cell-type MeCP2-repressed genes"="lightpink" ,"Cell-type MeCP2-activated genes"="blue", "Other-cell-type MeCP2-activated genes"="cyan3")) +
  facet_grid(.~mutant+genomic_loc,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_cellvsOther_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.778946e-163, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=4.985179e-39, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.04061951, *

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=7.851781e-07, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.5937957, ns
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="KO/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.01653431, *


wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=2.686386e-156, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=8.2771e-24, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Intragenic") & (mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=4.028235e-08, ****

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=8.801089e-06, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Cell-type MeCP2-repressed genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.1397587, ns
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Cell-type unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_otherCellType_MR_genes_nondedup_intra_and_extra[(genomic_loc=="Extragenic") & (mutant=="OE/WT") & (geneset=="Other-cell-type MeCP2-repressed genes"), as.numeric(logFC)])$p.value #p=0.3523627, ns


##changing plotting order around
#KO
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Pv cCREs,\n unchanged genes"),
                                                                                                                         cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Non-Pv cCREs,\n unchanged genes"),
                                                                                                                         cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Pv cCREs,\n MR genes"),
                                                                                                                         cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Non-Pv cCREs,\n MR genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))


ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked, aes(x = cCRE_group, y = as.numeric(logFC), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (KO/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle=90))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2fcH3K27ac_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2fcH3K27ac_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

#TG
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", cCRE_group="Pv cCREs,\n unchanged genes"),
                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", cCRE_group="Non-Pv cCREs,\n unchanged genes"),
                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", cCRE_group="Pv cCREs,\n MR genes"),
                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="OE/WT", cCRE_group="Non-Pv cCREs,\n MR genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))


ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked, aes(x = cCRE_group, y = as.numeric(logFC), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (OE/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle=90))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2fcH3K27ac_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2fcH3K27ac_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value


wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intragenicLinked[cCRE_group=="Non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value

#intragenic-linked cCRE table
mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(Intragenic==1) & (Intragenic_to_linked_gene==1), ]
mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(Intragenic==0) & (Intragenic_to_linked_gene==0), ]
##PV meta-MR and non-PV meta-MR genes
PvMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
PvUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")

#Make a data table containing the cCREs linked to PV meta-MR and non-PV meta-MR genes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked = data.table(rbind(
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="Unchanged genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="PV meta-MR genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="Non-PV meta-MR genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="Unchanged genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="PV meta-MR genes"),
                                                                                                                                        cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="Non-PV meta-MR genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV meta-MR genes", "Non-PV meta-MR genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV meta-MR genes"="palevioletred4", "Non-PV meta-MR genes"="palevioletred2")) +
  facet_grid(.~mutant,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.7,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_metaMR_genes_coding_prefilt5_nondedup_cellType_vs_nonCellType_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_intragenicLinked_to_metaMR_genes_coding_prefilt5_nondedup_cellType_vs_nonCellType_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="KO/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_intragenicLinked[(mutant=="OE/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value

data.table(val=mousebrain_union_cCREs_mPv_mCA[is.finite(cCRE_methylation) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes")


#extragenic
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked = data.table(rbind(
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="Unchanged genes"),
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="PV meta-MR genes"),
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[,V4], cCRE_label]), logFC], mutant="KO/WT", geneset="Non-PV meta-MR genes"),
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="Unchanged genes"),
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="PV meta-MR genes"),
  cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[,V4], cCRE_label]), logFC], mutant="OE/WT", geneset="Non-PV meta-MR genes")))
#Define the ordering of the factors for plotting purposes
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV meta-MR genes", "Non-PV meta-MR genes")))
union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked = union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked %>% mutate(mutant = factor(mutant, levels=c("KO/WT", "OE/WT")))

ggplot(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked, aes(x = geneset, y = as.numeric(logFC), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV meta-MR genes"="palevioletred4", "Non-PV meta-MR genes"="palevioletred2")) +
  facet_grid(.~mutant,
             switch = "x",
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-1.7,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (mutant/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_metaMR_genes_coding_prefilt5_nondedup_cellType_vs_nonCellType_PvMeCP2mut_log2fcH3K27ac_boxplot.png", width = 4.1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/union_nonPromoter_1500bpWindow_cCREs_extragenicLinked_to_metaMR_genes_coding_prefilt5_nondedup_cellType_vs_nonCellType_PvMeCP2mut_log2fcH3K27ac_boxplot.eps", width = 4.1, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)])$p.value #p=0.004731562, **
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value #p=0.03257501, *
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="KO/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value #p=0.646045, ns

wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)])$p.value #p=4.899978e-08, ****
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="PV meta-MR genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value #p=0.002874359, **"
wilcox.test(union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Unchanged genes"), as.numeric(logFC)], union_nonPromoter_cCREs_1500bp_Pv_MeCP2mut_H3K27ac_ChIP_cellType_vs_nonCellType_metaMR_genes_nondedup_extragenicLinked[(mutant=="OE/WT") & (geneset=="Non-PV meta-MR genes"), as.numeric(logFC)])$p.value #P=0.0006141346, ***






###adding extragenic cCREs
#KO
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Intragenic Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Intragenic non-Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Intragenic Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Intragenic non-Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Extragenic Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Extragenic non-Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Extragenic Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="KO/WT", cCRE_group="Extragenic non-Pv cCREs,\n MR genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Intragenic Pv cCREs,\n unchanged genes", "Intragenic non-Pv cCREs,\n unchanged genes", "Intragenic Pv cCREs,\n MR genes", "Intragenic non-Pv cCREs,\n MR genes",
                                                                                                                                                                                                                                                              "Extragenic Pv cCREs,\n unchanged genes", "Extragenic non-Pv cCREs,\n unchanged genes", "Extragenic Pv cCREs,\n MR genes", "Extragenic non-Pv cCREs,\n MR genes")))


ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra, aes(x = cCRE_group, y = as.numeric(logFC), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Intragenic Pv cCREs,\n unchanged genes"="gray", "Intragenic non-Pv cCREs,\n unchanged genes"="gray", "Intragenic Pv cCREs,\n MR genes"="red", "Intragenic non-Pv cCREs,\n MR genes"="red",
                                                          "Extragenic Pv cCREs,\n unchanged genes"="gray", "Extragenic non-Pv cCREs,\n unchanged genes"="gray", "Extragenic Pv cCREs,\n MR genes"="red", "Extragenic non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (KO/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle=90))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2fcH3K27ac_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2fcH3K27ac_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=1.106535e-10, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=1.132884e-12, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=6.500396e-12, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=0.01834795, *

#TG
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = data.table(rbind(cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Intragenic Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Intragenic non-Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Intragenic Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Intragenic non-Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Extragenic Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Extragenic non-Pv cCREs,\n unchanged genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Extragenic Pv cCREs,\n MR genes"),
                                                                                                                       cbind(logFC=nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_new_table[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4]), logFC], mutant="TG/WT", cCRE_group="Extragenic non-Pv cCREs,\n MR genes")))
PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra = PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Intragenic Pv cCREs,\n unchanged genes", "Intragenic non-Pv cCREs,\n unchanged genes", "Intragenic Pv cCREs,\n MR genes", "Intragenic non-Pv cCREs,\n MR genes",
                                                                                                                                                                                                                                                              "Extragenic Pv cCREs,\n unchanged genes", "Extragenic non-Pv cCREs,\n unchanged genes", "Extragenic Pv cCREs,\n MR genes", "Extragenic non-Pv cCREs,\n MR genes")))


ggplot(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra, aes(x = cCRE_group, y = as.numeric(logFC), fill=cCRE_group))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Intragenic Pv cCREs,\n unchanged genes"="gray", "Intragenic non-Pv cCREs,\n unchanged genes"="gray", "Intragenic Pv cCREs,\n MR genes"="red", "Intragenic non-Pv cCREs,\n MR genes"="red",
                                                          "Extragenic Pv cCREs,\n unchanged genes"="gray", "Extragenic non-Pv cCREs,\n unchanged genes"="gray", "Extragenic Pv cCREs,\n MR genes"="red", "Extragenic non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-1.6,1.55))+
  ylab("Log2 H3K27ac ChIPseq fold change (TG/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle=90))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2fcH3K27ac_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2fcH3K27ac_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=4.211322e-08, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=0.0001977753, ***
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=1.189417e-34, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2TG_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=0.0001125783, ***
