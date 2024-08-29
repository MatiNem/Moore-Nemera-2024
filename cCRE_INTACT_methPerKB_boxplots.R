 #This code produces boxplots of PV mCA per KB and mCG per KB for cCREs
library(data.table)
library(dplyr)
library(ggplot2)

cCRE_INTACT_PV_mCA_withCAcounts = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCA_mm9_withCAcounts.bed_s")
cCRE_INTACT_PV_mCG_withCGcounts = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCG_mm9_withCGcounts.bed_s")


cCRE_INTACT_PV_mCA_withCAcounts[, cCRE_methylation := as.integer(V5)/as.integer(V6)]
cCRE_INTACT_PV_mCG_withCGcounts[, cCRE_methylation := as.integer(V5)/as.integer(V6)]

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

#subtracting nonconversion rate
cCRE_INTACT_PV_mCA_withCAcounts$cCRE_methylation_corrected <- cCRE_INTACT_PV_mCA_withCAcounts$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
cCRE_INTACT_PV_mCG_withCGcounts$cCRE_methylation_corrected <- cCRE_INTACT_PV_mCG_withCGcounts$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]


#turn negative methylation values into zeros
cCRE_INTACT_PV_mCA_withCAcounts[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
cCRE_INTACT_PV_mCG_withCGcounts[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]

#calculating CA and CG methylation density across cCREs
cCRE_INTACT_PV_mCA_withCAcounts[, meth_per_kb:=(1000*as.integer(V7)/(V3-V2)*cCRE_methylation_corrected)]
cCRE_INTACT_PV_mCG_withCGcounts[, meth_per_kb:=(1000*as.integer(V7)/(V3-V2)*cCRE_methylation_corrected)]

#Fisher-combined MeCP2-regulated cCREs
union_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
union_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
union_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#data table of mCA/kb of MeCP2-regulated cCREs
MeCP2reg_cCREs_INTACT_mCAperKB = rbind(
           cbind(cCRE_INTACT_PV_mCA_withCAcounts[V4 %in% union_unchanged_cCREs_sig[, V4],], cCRE_set = "All other cCREs"),
           cbind(cCRE_INTACT_PV_mCA_withCAcounts[V4 %in% union_mr_cCREs_sig[, V4],], cCRE_set = "MR cCREs"),
           cbind(cCRE_INTACT_PV_mCA_withCAcounts[V4 %in% union_ma_cCREs_sig[, V4],], cCRE_set = "MA cCREs"))
MeCP2reg_cCREs_INTACT_mCAperKB = MeCP2reg_cCREs_INTACT_mCAperKB %>% mutate(cCRE_set = factor(cCRE_set, levels=c("All other cCREs", "MR cCREs", "MA cCREs")))
#boxplot of mCA/kb
ggplot(MeCP2reg_cCREs_INTACT_mCAperKB ,aes(x = cCRE_set, y = meth_per_kb, fill = cCRE_set)) +
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("All other cCREs"="gray", "MR cCREs"="red", "MA cCREs"="blue")) +
  coord_cartesian(ylim = c(0, 25)) +
  ylab("PV INTACT mCA/kb") +
  theme_bw() +
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/MeCP2dys_nonPromoter_cCREs_Fishercomb_PV_WT_KO_deep_INTACT_mCAperKB_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device="png")
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/MeCP2dys_nonPromoter_cCREs_Fishercomb_PV_WT_KO_deep_INTACT_mCAperKB_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device="eps")

wilcox.test(MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="All other cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="MR cCREs", as.numeric(meth_per_kb)])$p.value #p= 2.007492e-21, ****
wilcox.test(MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="MR cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="MA cCREs", as.numeric(meth_per_kb)])$p.value #p=<2.2e-16, ****
wilcox.test(MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="All other cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCAperKB[cCRE_set=="MA cCREs", as.numeric(meth_per_kb)])$p.value #p=2.504629e-303, ****


#data table of mCG/kb of MeCP2-regulated cCREs
MeCP2reg_cCREs_INTACT_mCGperKB = rbind(
  cbind(cCRE_INTACT_PV_mCG_withCGcounts[V4 %in% union_unchanged_cCREs_sig[, V4],], cCRE_set = "All other cCREs"),
  cbind(cCRE_INTACT_PV_mCG_withCGcounts[V4 %in% union_mr_cCREs_sig[, V4],], cCRE_set = "MR cCREs"),
  cbind(cCRE_INTACT_PV_mCG_withCGcounts[V4 %in% union_ma_cCREs_sig[, V4],], cCRE_set = "MA cCREs"))
MeCP2reg_cCREs_INTACT_mCGperKB = MeCP2reg_cCREs_INTACT_mCGperKB %>% mutate(cCRE_set = factor(cCRE_set, levels=c("All other cCREs", "MR cCREs", "MA cCREs")))
#boxplot of mCG/kb
ggplot(MeCP2reg_cCREs_INTACT_mCGperKB ,aes(x = cCRE_set, y = meth_per_kb, fill = cCRE_set)) +
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("All other cCREs"="gray", "MR cCREs"="red", "MA cCREs"="blue")) +
  coord_cartesian(ylim = c(0, 40)) +
  ylab("PV INTACT mCG/kb") +
  theme_bw() +
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/MeCP2dys_nonPromoter_cCREs_Fishercomb_PV_WT_KO_deep_INTACT_mCGperKB_boxplot.png", width = 3.7, height = 5, dpi = 300, units = "in", device="png")
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/enhancer_meth_plots/MeCP2dys_nonPromoter_cCREs_Fishercomb_PV_WT_KO_deep_INTACT_mCGperKB_boxplot.eps", width = 3.7, height = 5, dpi = 300, units = "in", device="eps")

wilcox.test(MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="All other cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="MR cCREs", as.numeric(meth_per_kb)])$p.value #p= 1.380723e-16, ****
wilcox.test(MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="MR cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="MA cCREs", as.numeric(meth_per_kb)])$p.value #p=1.357855e-230, ****
wilcox.test(MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="All other cCREs", as.numeric(meth_per_kb)], MeCP2reg_cCREs_INTACT_mCGperKB[cCRE_set=="MA cCREs", as.numeric(meth_per_kb)])$p.value #p=1.514743e-263, ****
