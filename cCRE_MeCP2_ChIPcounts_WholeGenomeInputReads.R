library(data.table)
library(dplyr)
library(ggplot2)

#read counts in PV MeCP2 ChIP and input
PV_MeCP2_ChIP_total_read_numbers <- fread("HG_lab/Mati/GabelLab/ChIPseq/Pv_MeCP2/Pv_MeCP2_ChIP_total_read_numbers.csv")

#non-promoter cCREs
nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")

#PV MeCP2 ChIP counts
mousebrain_union_cCREs_mecp2_pv_stroud_8wk_chipCounts_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_MeCP2_chipCounts/mousebrain_union_cCREs_mecp2_pv_stroud_8wk_chipCounts_mm9.bed")
#PV MeCP2 Input counts
mousebrain_union_cCREs_input_pv_stroud_8wk_chipCounts_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_MeCP2_chipCounts/mousebrain_union_cCREs_input_pv_stroud_8wk_chipCounts_mm9.bed")
#cpm calculation using total number of reads mapped to the genome as denominator
mousebrain_union_cCREs_mecp2_pv_stroud_8wk_chipCounts_mm9[, cpm := (10**6) * (V5/PV_MeCP2_ChIP_total_read_numbers[filename=="GSM2803629_MECP2_ChIP_Pv_8wk_Cortex.txt", read_count])]
mousebrain_union_cCREs_input_pv_stroud_8wk_chipCounts_mm9[, cpm := (10**6) * (V5/PV_MeCP2_ChIP_total_read_numbers[filename=="GSM2803631_Input_Pv_8wk_Cortex.txt", read_count])]


mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9 = data.table(inner_join(x=mousebrain_union_cCREs_mecp2_pv_stroud_8wk_chipCounts_mm9, y=mousebrain_union_cCREs_input_pv_stroud_8wk_chipCounts_mm9, by=c("V1", "V2", "V3", "V4")))
names(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "MeCP2_count", "MeCP2_cpm", "Input_count", "Input_cpm")
mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[, log2_ChIP_Input_ratio_pseud := log2((MeCP2_cpm + 1)/(Input_cpm + 1))]

union_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
union_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
union_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% union_unchanged_cCREs_sig[,V4]], cCRE_set = "Unchanged cCREs"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% union_mr_cCREs_sig[,V4]], cCRE_set = "MR cCREs"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% union_ma_cCREs_sig[,V4]], cCRE_set = "MA cCREs")))

union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input = union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input %>% mutate(cCRE_set = factor(cCRE_set, levels=c("Unchanged cCREs", "MR cCREs", "MA cCREs")))

ggplot(union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input, aes(x = cCRE_set, y = log2_ChIP_Input_ratio_pseud, fill=cCRE_set))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("Unchanged cCREs"="gray", "MR cCREs"="red", "MA cCREs"="blue")) +
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
#ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercomb_MeCP2dys_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercomb_MeCP2dys_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="Unchanged cCREs", log2_ChIP_Input_ratio_pseud], union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="MR cCREs", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="MR cCREs", log2_ChIP_Input_ratio_pseud], union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="MA cCREs", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="Unchanged cCREs", log2_ChIP_Input_ratio_pseud], union_MeCP2reg_cCREs_MeCP2_ChIP_over_Input[cCRE_set=="MA cCREs", log2_ChIP_Input_ratio_pseud])$p.value



#all cCREs linked to MeCP2-regulated genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")



cCREs_linked_to_genes_MeCP2_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "MR genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "MA genes")))

cCREs_linked_to_genes_MeCP2_ChIP_over_Input = cCREs_linked_to_genes_MeCP2_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "MR genes", "MA genes")))

ggplot(cCREs_linked_to_genes_MeCP2_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "MR genes"="red", "MA genes"="blue")) +
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_Pv_MeCP2reg_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_Pv_MeCP2reg_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="MR genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="MR genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="MA genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_genes_MeCP2_ChIP_over_Input[geneset=="MA genes", log2_ChIP_Input_ratio_pseud])$p.value


intragenic_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain.union.cCRE.ensgene.mm9.bed")
names(intragenic_cCREs_mm9) = c("cCRE_chrom","cCRE_start","cCRE_end","cCRE_label","gene_chrom","gene_start","gene_end","Gene","gene_num_transcripts","gene_strand")



union_mr_cCREs_sig[, genic_location := "Extragenic"]
union_mr_cCREs_sig[V4 %in% intragenic_cCREs_mm9[, cCRE_label], genic_location := "Intragenic"]
union_ma_cCREs_sig[, genic_location := "Extragenic"]
union_ma_cCREs_sig[V4 %in% intragenic_cCREs_mm9[, cCRE_label], genic_location := "Intragenic"]
union_unchanged_cCREs_sig[, genic_location := "Extragenic"]
union_unchanged_cCREs_sig[V4 %in% intragenic_cCREs_mm9[, cCRE_label], genic_location := "Intragenic"]


union_unchanged_intragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_table[genic_location=="Intragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_table), 3)
union_unchanged_extragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_table[genic_location=="Extragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_table), 3)
union_MR_intragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_table[genic_location=="Intragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_table), 3)
union_MR_extragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_table[genic_location=="Extragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_table), 3)
union_MA_intragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_table[genic_location=="Intragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_table), 3)
union_MA_extragenicProp = round(nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_table[genic_location=="Extragenic", ])/nrow(PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_table), 3)

#cCRE ChIP table with cCREs linked to Pv MR and other-cell-type MR genes
cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Other-cell-type MR genes")))

cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV MR genes", "Other-cell-type MR genes")))

ggplot(cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV MR genes"="red", "Other-cell-type MR genes"="lightpink")) +
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="PV MR genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="PV MR genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Other-cell-type MR genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_linked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Other-cell-type MR genes", log2_ChIP_Input_ratio_pseud])$p.value


#intragenic cognate-linked cCREs
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

#extragenic linked cCREs
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")


#cCRE ChIP table with cCREs linked to Pv MR and other-cell-type MR genes
cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Other-cell-type MR genes")))

cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV MR genes", "Other-cell-type MR genes")))

ggplot(cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV MR genes"="red", "Other-cell-type MR genes"="lightpink")) +
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="PV MR genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="PV MR genes", log2_ChIP_Input_ratio_pseud], cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Other-cell-type MR genes", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Unchanged genes", log2_ChIP_Input_ratio_pseud], cCREs_intragenicLinked_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[geneset=="Other-cell-type MR genes", log2_ChIP_Input_ratio_pseud])$p.value


Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

#extragenic linked
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, PV cCREs"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, Non-PV cCREs"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, PV cCREs"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, Non-PV cCREs")))


PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input = PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes, PV cCREs", "Unchanged genes, Non-PV cCREs", "PV MR genes, PV cCREs", "PV MR genes, Non-PV cCREs")))

ggplot(PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes, PV cCREs"="gray", "Unchanged genes, Non-PV cCREs"="gray", "PV MR genes, PV cCREs"="red", "PV MR genes, Non-PV cCREs"="red")) +
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intragenicLinked_to_PvMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intragenicLinked_to_PvMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input[geneset=="Unchanged genes, PV cCREs", log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input[geneset=="Unchanged genes, Non-PV cCREs", log2_ChIP_Input_ratio_pseud])$p.value
wilcox.test(PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input[geneset=="PV MR genes, PV cCREs", log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intragenicLinked_to_PvMR_genes_ChIP_over_Input[geneset=="PV MR genes, Non-PV cCREs", log2_ChIP_Input_ratio_pseud])$p.value


###cell-type meta-MR genes vs non-cell-type meta MR genes
##PV meta-MR and non-PV meta-MR genes
PvMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
PvUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")

mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(Intragenic==1) & (Intragenic_to_linked_gene==1), ]

cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(log2_ChIP_Input_ratio_pseud=mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[,V4], cCRE_label], log2_ChIP_Input_ratio_pseud], geneset = "Unchanged genes"),
  cbind(log2_ChIP_Input_ratio_pseud=mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[,V4], cCRE_label], log2_ChIP_Input_ratio_pseud], geneset = "PV meta-MR genes"),
  cbind(log2_ChIP_Input_ratio_pseud=mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[,V4], cCRE_label], log2_ChIP_Input_ratio_pseud], geneset = "Non-PV meta-MR genes")))

cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input = cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV meta-MR genes", "Non-PV meta-MR genes")))

ggplot(cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input, aes(x = geneset, y = as.numeric(log2_ChIP_Input_ratio_pseud), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV meta-MR genes"="palevioletred4", "Non-PV meta-MR genes"="palevioletred2")) +
  #coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="Unchanged genes", as.numeric(log2_ChIP_Input_ratio_pseud)], cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="PV meta-MR genes", as.numeric(log2_ChIP_Input_ratio_pseud)])$p.value
wilcox.test(cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="PV meta-MR genes", as.numeric(log2_ChIP_Input_ratio_pseud)], cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="Non-PV meta-MR genes", as.numeric(log2_ChIP_Input_ratio_pseud)])$p.value
wilcox.test(cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="Unchanged genes", as.numeric(log2_ChIP_Input_ratio_pseud)], cCREs_linked_to_PvMetaMR_nonPvMetaMR_genes_ChIP_over_Input[geneset=="Non-PV meta-MR genes", as.numeric(log2_ChIP_Input_ratio_pseud)])$p.value




#####
cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Other-cell-type MR genes", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes", genic_loc="Extragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes", genic_loc="Extragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Other-cell-type MR genes", genic_loc="Extragenic")))

cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes", "PV MR genes", "Other-cell-type MR genes")))
cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input = cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))

ggplot(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes"="gray", "PV MR genes"="red", "Other-cell-type MR genes"="lightpink")) +
  coord_cartesian(ylim=c(-0.55,0.55))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  facet_grid(.~genic_loc,
             labeller = label_value,
             switch = "x"
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_intra_and_extra_linked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_WholeGenome_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/mousebrain_union_nonPromoter_cCREs_intra_and_extra_linked_to_PvMR_otherCellTypeMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_WholeGenome_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Unchanged genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="PV MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=2.75219e-71, ****
wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="PV MR genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Other-cell-type MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=5.50507e-13, ****
wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Unchanged genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Other-cell-type MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=0.006846553, **

wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Unchanged genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="PV MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=0.5707518, ns
wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="PV MR genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Other-cell-type MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=0.2021423, ns
wilcox.test(cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Unchanged genes"), log2_ChIP_Input_ratio_pseud], cCREs_intra_and_extra_to_PvMR_otherCellTypeMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Other-cell-type MR genes"), log2_ChIP_Input_ratio_pseud])$p.value #p=0.07092498, ns


##
PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input = data.table(rbind(
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, PV cCREs", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, Non-PV cCREs", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, PV cCREs", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, Non-PV cCREs", genic_loc="Intragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, PV cCREs", genic_loc="Extragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "Unchanged genes, Non-PV cCREs", genic_loc="Extragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, PV cCREs", genic_loc="Extragenic"),
  cbind(mousebrain_union_cCREs_mecp2_input_pv_stroud_8wk_chipCounts_mm9[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[,V4]], geneset = "PV MR genes, Non-PV cCREs", genic_loc="Extragenic")))


PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input = PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input %>% mutate(geneset = factor(geneset, levels=c("Unchanged genes, PV cCREs", "Unchanged genes, Non-PV cCREs", "PV MR genes, PV cCREs", "PV MR genes, Non-PV cCREs")))
PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input = PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))

ggplot(PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input, aes(x = geneset, y = log2_ChIP_Input_ratio_pseud, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Unchanged genes, PV cCREs"="gray", "Unchanged genes, Non-PV cCREs"="gray", "PV MR genes, PV cCREs"="red", "PV MR genes, Non-PV cCREs"="red")) +
  coord_cartesian(ylim=c(-0.55,0.55))+
  ylab("Log2 MeCP2 ChIP/Input cpm ratio") + xlab("")+
  facet_grid(.~genic_loc,
             labeller = label_value,
             switch = "x"
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intra_and_extra_linked_to_PvMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_wholeGenomeReads_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_MeCP2/cell_confusion_paper_chip_plots/PVGA_and_nonPVGA_nonPromoter_cCREs_intra_and_extra_linked_to_PvMR_genes_H3K27ac_ChIP_log2MeCP2_ChIP_over_Input_pv_stroud_8wk_wholeGenomeReads_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Unchanged genes, PV cCREs"), log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="Unchanged genes, Non-PV cCREs"), log2_ChIP_Input_ratio_pseud])$p.value #p=1.285976e-189, ****
wilcox.test(PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="PV MR genes, PV cCREs"), log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Intragenic") & (geneset=="PV MR genes, Non-PV cCREs"), log2_ChIP_Input_ratio_pseud])$p.value #p=8.14783e-38, ****

wilcox.test(PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Unchanged genes, PV cCREs"), log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="Unchanged genes, Non-PV cCREs"), log2_ChIP_Input_ratio_pseud])$p.value #p=1.545262e-203, ****
wilcox.test(PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="PV MR genes, PV cCREs"), log2_ChIP_Input_ratio_pseud], PVGA_and_nonPVGA_cCREs_intra_and_extra_to_PvMR_genes_ChIP_over_Input[(genic_loc=="Extragenic") & (geneset=="PV MR genes, Non-PV cCREs"), log2_ChIP_Input_ratio_pseud])$p.value #p=1.304967e-18, ****
