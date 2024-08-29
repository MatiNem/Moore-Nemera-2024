library(data.table)
library(dplyr)
library(ggplot2)

#all lambda nonconversion rates of INTACT bisulfite experiments
lambda_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/intact_lambda_nonconversion_table.csv")


#cCREs with mCA
union_cCRE_PV_WT_rep1_mCA <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA_mm9.bed")
union_cCRE_PV_WT_rep2_mCA <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA_mm9.bed")
union_cCRE_PV_KO_rep1_mCA <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA_mm9.bed")
union_cCRE_PV_KO_rep2_mCA <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA_mm9.bed")



#cCREs with mCG
union_cCRE_PV_WT_rep1_mCG <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCG_mm9.bed")
union_cCRE_PV_WT_rep2_mCG <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCG_mm9.bed")
union_cCRE_PV_KO_rep1_mCG <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCG_mm9.bed")
union_cCRE_PV_KO_rep2_mCG <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCG_mm9.bed")


#naming methylation table columns and calculating methylation levels
cCRE_meth_calc_func <- function(meth_table, nonconv_table, label_column){
  meth_table$cCRE_methylation <- as.integer(meth_table[[5]])/as.integer(meth_table[[6]])
  meth_table$cCRE_methylation_corrected <- meth_table$cCRE_methylation - nonconv_table[label==label_column, nonconversion_rate]
  names(meth_table) = c("chrom", "start", "end", "cCRE_label", "meth_reads", "cyto_reads", "cCRE_methylation", "cCRE_methylation_corrected")
  meth_table[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
  return(meth_table)
}

#methylation calculation with background nonconversion rate subtraction
union_cCRE_PV_WT_rep1_mCA <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_WT_rep1_mCA, nonconv_table = lambda_nonconv, label_column="PV_WT_rep1")
union_cCRE_PV_WT_rep2_mCA <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_WT_rep2_mCA, nonconv_table = lambda_nonconv, label_column="PV_WT_rep2")
union_cCRE_PV_KO_rep1_mCA <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_KO_rep1_mCA, nonconv_table = lambda_nonconv, label_column="PV_KO_rep1")
union_cCRE_PV_KO_rep2_mCA <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_KO_rep2_mCA, nonconv_table = lambda_nonconv, label_column="PV_KO_rep2")

union_cCRE_PV_WT_rep1_mCG <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_WT_rep1_mCG, nonconv_table = lambda_nonconv, label_column="PV_WT_rep1")
union_cCRE_PV_WT_rep2_mCG <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_WT_rep2_mCG, nonconv_table = lambda_nonconv, label_column="PV_WT_rep2")
union_cCRE_PV_KO_rep1_mCG <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_KO_rep1_mCG, nonconv_table = lambda_nonconv, label_column="PV_KO_rep1")
union_cCRE_PV_KO_rep2_mCG <- cCRE_meth_calc_func(meth_table=union_cCRE_PV_KO_rep2_mCG, nonconv_table = lambda_nonconv, label_column="PV_KO_rep2")


#
#Fisher-combined MeCP2-regulated cCREs
union_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
union_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
union_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#data table of mCA/kb of MeCP2-regulated cCREs
#mCA
reps_MeCP2reg_cCREs_INTACT_mCAperCA <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MA cCREs")
) %>% data.table


reps_MeCP2reg_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_cCREs_INTACT_mCAperCA %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_MeCP2reg_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_cCREs_INTACT_mCAperCA %>% mutate(subclass = factor(cCRE_group, levels=c("All other cCREs", "MR cCREs", "MA cCREs")))
reps_MeCP2reg_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_cCREs_INTACT_mCAperCA %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))


#plotting mean INTACT PV methylation (background subtracted) of unchanged, MR, and MA cCREs
ggplot(reps_MeCP2reg_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.07)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="All other cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="All other cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MR cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MR cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333
wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MA cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MA cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1

#
#mCG
reps_MeCP2reg_cCREs_INTACT_mCGperCG <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% union_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="All other cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% union_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MR cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MA cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% union_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MA cCREs")
) %>% data.table


reps_MeCP2reg_cCREs_INTACT_mCGperCG <- reps_MeCP2reg_cCREs_INTACT_mCGperCG %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_MeCP2reg_cCREs_INTACT_mCGperCG <- reps_MeCP2reg_cCREs_INTACT_mCGperCG %>% mutate(subclass = factor(cCRE_group, levels=c("All other cCREs", "MR cCREs", "MA cCREs")))
reps_MeCP2reg_cCREs_INTACT_mCGperCG <- reps_MeCP2reg_cCREs_INTACT_mCGperCG %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))


#plotting mean INTACET PV methylation (background subtracted) of unchanged, MR, and MA cCREs
ggplot(reps_MeCP2reg_cCREs_INTACT_mCGperCG, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("PV INTACT mCG/CG") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_cCRE_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_cCRE_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="WT") & (cCRE_group=="All other cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="KO") & (cCRE_group=="All other cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="WT") & (cCRE_group=="MR cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="KO") & (cCRE_group=="MR cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333
wilcox.test(reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="WT") & (cCRE_group=="MA cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_cCREs_INTACT_mCGperCG[(genotype=="KO") & (cCRE_group=="MA cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667


###
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
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n MR genes")
) %>% data.table

reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))

#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs intragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333



#INTACT PV mCA of PV and non-PV extragenic-linked cCREs
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n MR genes")
) %>% data.table

reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))

#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs extragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667

#mCG
#INTACT PV mCG of PV and non-PV intragenic-linked cCREs
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n MR genes")
) %>% data.table

reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))

#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs intragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("PV INTACT mCG/CG") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667



#INTACT PV mCG of PV and non-PV extragenic-linked cCREs
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n unchanged genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs,\n MR genes"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCG[cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs,\n MR genes")
) %>% data.table

reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs,\n unchanged genes", "Non-PV cCREs,\n unchanged genes", "PV cCREs,\n MR genes", "Non-PV cCREs,\n MR genes")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra = reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))

#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs extragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("PV INTACT mCG/CG") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCGperCG_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCGperCG_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1

##
#separated Mecp2-regulated cCREs by PV and non-PV cCREs
#PV
PVGA_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
PVGA_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
PVGA_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")
#non-PV
nonPVGA_mr_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")
nonPVGA_ma_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")
nonPVGA_unchanged_cCREs_sig = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")

#
#mCA
reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% PVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="All other PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% PVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="All other PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% PVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="All other PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% PVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="All other PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% PVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MR PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% PVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MR PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% PVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MR PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% PVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MR PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% PVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MA PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% PVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MA PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% PVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MA PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% PVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MA PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% nonPVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="All other non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% nonPVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="All other non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% nonPVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="All other non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% nonPVGA_unchanged_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="All other non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% nonPVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MR non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% nonPVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MR non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% nonPVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MR non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% nonPVGA_mr_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MR non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% nonPVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="MA non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% nonPVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="MA non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% nonPVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="MA non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% nonPVGA_ma_cCREs_sig$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="MA non-PV cCREs")
) %>% data.table


reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(cCRE_group = factor(cCRE_group, levels=c("All other PV cCREs", "All other non-PV cCREs", "MR PV cCREs", "MR non-PV cCREs", "MA PV cCREs", "MA non-PV cCREs")))
reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))


#plotting mean INTACT PV methylation (background subtracted) of unchanged, MR, and MA cCREs
ggplot(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/MeCP2reg_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="All other PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="All other PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MR PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MR PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MA PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MA PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="All other non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="All other non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MR non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MR non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MA non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MA non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1



####
#all PVGA and non-PVGA nonpromoter cCREs
PVGA_nonpromoter_cCREs <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/PVGA/PVGA_nonPromoter_cCREs_mm9.bed")
nonPVGA_nonpromoter_cCREs <- fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nonPVGA_nonPromoter_cCREs_mm9.bed")
#PVGA and non-PVGA nonpromoter cCREs that are linked to genes
PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
#

#all PVGA and non-PVGA nonpromoter cCREs
reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% PVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% PVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% PVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% PVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% nonPVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% nonPVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% nonPVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% nonPVGA_nonpromoter_cCREs$V4, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Non-PV cCREs")
) %>% data.table

reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(cCRE_group = factor(cCRE_group, levels=c("PV cCREs", "Non-PV cCREs")))
reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))


#plotting mean INTACT PV methylation (background subtracted) of all PVGA and non-PVGA non-promoter cCREs
ggplot(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/all_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/all_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="PV cCREs"), as.numeric(mean_cCRE_meth)], reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667

###
reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- rbind(
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Linked PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Linked PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Linked PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Linked PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep1_mCA[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041642", cCRE_group="Linked non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_WT_rep2_mCA[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="WT", gtac_id="LIB041644", cCRE_group="Linked non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep1_mCA[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041641", cCRE_group="Linked non-PV cCREs"),
  cbind(mean_cCRE_meth=mean(union_cCRE_PV_KO_rep2_mCA[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans$cCRE_label, cCRE_methylation_corrected], na.rm=TRUE), genotype="KO", gtac_id="LIB041643", cCRE_group="Linked non-PV cCREs")
) %>% data.table

reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Linked PV cCREs", "Linked non-PV cCREs")))
reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA <- reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA %>% mutate(gtac_id = factor(gtac_id, levels=c("LIB041642", "LIB041644", "LIB041641", "LIB041643")))


#plotting mean INTACT PV methylation (background subtracted) of linked PVGA and non-PVGA non-promoter cCREs
ggplot(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=15), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/linked_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/linked_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Linked PV cCREs"), as.numeric(mean_cCRE_meth)], reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Linked PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Linked non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Linked non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667



###bigger versions
#plotting mean INTACT PV methylation (background subtracted) of unchanged, MR, and MA cCREs
ggplot(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_MeCP2reg_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_MeCP2reg_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="All other PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="All other PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MR PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MR PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MA PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MA PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="All other non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="All other non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MR non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MR non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333
wilcox.test(reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="MA non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_MeCP2reg_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="MA non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1

#plotting mean INTACT PV methylation (background subtracted) of all PVGA and non-PVGA non-promoter cCREs
ggplot(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_all_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 2, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_all_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 2, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="PV cCREs"), as.numeric(mean_cCRE_meth)], reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=1
wilcox.test(reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_all_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667


ggplot(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=15), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_linked_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 2, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_linked_PV_and_nonPV_cCRE_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 2, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Linked PV cCREs"), as.numeric(mean_cCRE_meth)], reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Linked PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="WT") & (cCRE_group=="Linked non-PV cCREs"), as.numeric(mean_cCRE_meth)], reps_linked_PVGA_and_nonPVGA_cCREs_INTACT_mCAperCA[(genotype=="KO") & (cCRE_group=="Linked non-PV cCREs"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667




#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs intragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 2.5, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_PVGA_and_nonPVGA_cCRE_intragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 2.5, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_intra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.3333333



#plotting mean INTACT PV methylation (background subtracted) of PV and non-PV cCREs extragenic linked to PV MR genes or unchanged genes
ggplot(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra, aes(x = genotype, y = as.numeric(mean_cCRE_meth), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.08)) +
  ylab("PV INTACT mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype\n", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~cCRE_group, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.png", width = 2.5, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/bigger_PVGA_and_nonPVGA_cCRE_extragenicLinked_PV_MR_genes_PV_INTACT_mCAperCA_backgroundSub_reps_stripplot.eps", width = 2.5, height = 6, dpi = 300, units = "in", device=cairo_ps)


wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n unchanged genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
wilcox.test(reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="WT") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)], reps_PVGA_nonPVGA_nonPromoter_cCREs_linked_to_PV_MR_genes_INTACT_mCAperCA_extra[(genotype=="KO") & (cCRE_group=="Non-PV cCREs,\n MR genes"), as.numeric(mean_cCRE_meth)])$p.value #p.value=0.6666667
