


library(tidyr)
library(dplyr)


mL4_snmcseq_gene_and_flank_mC = read.table("~/Documents/Gabel_Lab/Data/Meth Analysis/Files/MN/mL4_snmcseq_gene_and_flank_mCA.txt", header=TRUE)
mL4_snmcseq_gene_and_flank_mC$gene_methylation = as.numeric(as.character(mL4_snmcseq_gene_and_flank_mC$gene_methylation))
mean(mL4_snmcseq_gene_and_flank_mC$gene_methylation, na.rm = TRUE)

mL5_snmcseq_gene_and_flank_mC = read.table("~/Documents/Gabel_Lab/Data/Meth Analysis/Files/MN/mL5_all_snmcseq_gene_and_flank_mCA.txt", header=TRUE)
mL5_snmcseq_gene_and_flank_mC$gene_methylation = as.numeric(as.character(mL5_snmcseq_gene_and_flank_mC$gene_methylation))
mean(mL5_snmcseq_gene_and_flank_mC$gene_methylation, na.rm = TRUE)

mSst_all_snmcseq_gene_and_flank_mC = read.table("~/Documents/Gabel_Lab/Data/Meth Analysis/Files/MN/mSst_all_snmcseq_gene_and_flank_mCA.txt", header=TRUE)
mSst_all_snmcseq_gene_and_flank_mC$gene_methylation = as.numeric(as.character(mSst_all_snmcseq_gene_and_flank_mC$gene_methylation))
mean(mSst_all_snmcseq_gene_and_flank_mC$gene_methylation, na.rm = TRUE)

mPv_snmcseq_gene_and_flank_mC = read.table("~/Documents/Gabel_Lab/Data/Meth Analysis/Files/MN/mPv_snmcseq_gene_and_flank_mCA.txt", header=TRUE)
mPv_snmcseq_gene_and_flank_mC$gene_methylation = as.numeric(as.character(mPv_snmcseq_gene_and_flank_mC$gene_methylation))
mean(mPv_snmcseq_gene_and_flank_mC$gene_methylation, na.rm = TRUE)


setwd("~/Documents/Gabel_Lab/Analysis/Mecp2_cs_plots/deseq_outputs/pv_ko/nondedup/exon/deseq_files/")
pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921 = read.table("pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv", header=TRUE)
pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921$Gene = rownames(pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921)

setwd("~/Documents/Gabel_Lab/Data/RNA seq/GDT/")
#intact_GDT_113020 = read.table("GDT_INTACT_RNA_113020.tsv", sep="\t", header=TRUE)
intact_GDT_081421 = read.table("GDT_INTACT_RNA_081421.tsv", sep="\t", header=TRUE)

head(intact_GDT_081421$MR_meta)

intact_GDT_081421_pv = as.data.frame(cbind(Gene = intact_GDT_081421$Gene,
                                           Chr = intact_GDT_081421$chr,
                                           Start = intact_GDT_081421$start,
                                           Stop = intact_GDT_081421$stop,
                                           luo_pv_snmc_mCA = intact_GDT_081421$luo_pv_snmc_mCA,
                                           luo_pv_n_mCA = intact_GDT_081421$luo_pv_n_mCA,
                                           luo_2_sst_mCA = intact_GDT_081421$luo_2_sst_mCA,
                                           mCA = intact_GDT_081421$luo_pv_n_mCA,
                                           stroud_8wk_mCA = intact_GDT_081421$stroud_8wk_mCA,
                                           pv_kotg_enh_count = intact_GDT_081421$pv_kotg_enh_count,
                                           pv_total_enh_count = intact_GDT_081421$pv_total_enh_count,
                                           li_atac_PV_cCRE_count = intact_GDT_081421$li_atac_PV_ccre_count,
                                           li_atac_union_ccre_nonTSS_count = intact_GDT_081421$li_atac_union_ccre_nonTSS_count,
                                           li_atac_union_ccre_all_count = intact_GDT_081421$li_atac_union_ccre_all_count,
                                           enh_count =intact_GDT_081421$enh_count,
                                           Length = intact_GDT_081421$Length,
                                           MR_meta = intact_GDT_081421$MR_meta))
intact_GDT_081421_pv$mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_n_mCA))
intact_GDT_081421_pv$luo_pv_snmc_mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_snmc_mCA))
intact_GDT_081421_pv$luo_pv_n_mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_n_mCA))
intact_GDT_081421_pv$luo_2_sst_mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_2_sst_mCA))
intact_GDT_081421_pv$stroud_8wk_mCA = as.numeric(as.character(intact_GDT_081421_pv$stroud_8wk_mCA))
intact_GDT_081421_pv$Length = as.numeric(as.character(intact_GDT_081421_pv$Length))

intact_GDT_081421_pv_m = merge(intact_GDT_081421_pv, pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921, by="Gene")
intact_GDT_081421_pv_m = merge(intact_GDT_081421_pv_m, mPv_snmcseq_gene_and_flank_mC, by.x="Gene", by.y="gene")








setwd("~/Documents/Gabel_Lab/Analysis/Mecp2_CS_plots/deseq_outputs/rbp4_ko/nondedup/deseq_tables")

rbp4_ko_exon_nondedup_coding_d3_dseq_rep_prefilt5_res_df_all_110921 = read.table("rbp4_ko_exon_nondedeup_coding_d3_dseq_rep_prefilt5_res_df_all_110921.tsv", header=TRUE)
rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921 = read.table("rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

mL5_ensgene_mm9_sum = read.table("~/Documents/Gabel_Lab/Data/Meth Analysis/Files/MN/mL5_ensgene_mm9_sum.bed", sep="\t")
head(mL5_ensgene_mm9_sum)

mL5_ensgene_mm9_sum$V7 = as.numeric(as.character(mL5_ensgene_mm9_sum$V7))
mL5_ensgene_mm9_sum$V8 = as.numeric(as.character(mL5_ensgene_mm9_sum$V8))
mL5_ensgene_mm9_sum$mCA = mL5_ensgene_mm9_sum$V7/mL5_ensgene_mm9_sum$V8
colnames(mL5_ensgene_mm9_sum)[4] = "Gene"


intact_GDT_081421_rbp4 = as.data.frame(cbind(Gene = intact_GDT_081421$Gene,
                                             Chr = intact_GDT_081421$chr,
                                             Start = intact_GDT_081421$start,
                                             Stop = intact_GDT_081421$stop,
                                             stroud_8wk_mCA = intact_GDT_081421$stroud_8wk_mCA,                                             
                                             luo_L4_mCA = intact_GDT_081421$luo_L4_mCA,
                                             pv_kotg_enh_count = intact_GDT_081421$pv_kotg_enh_count,
                                             pv_total_enh_count = intact_GDT_081421$pv_total_enh_count,
                                             li_atac_ITL4_cCRE_count = intact_GDT_081421$li_atac_ITL4_cCRE_count,
                                             li_atac_union_ccre_nonTSS_count = intact_GDT_081421$li_atac_union_ccre_nonTSS_count,
                                             li_atac_union_ccre_all_count = intact_GDT_081421$li_atac_union_ccre_all_count,
                                             enh_count =intact_GDT_081421$enh_count,
                                             Length = intact_GDT_081421$Length,
                                             MR_meta = intact_GDT_081421$MR_meta))
intact_GDT_081421_rbp4$Length = as.numeric(as.character(intact_GDT_081421_rbp4$Length))
intact_GDT_081421_rbp4$log2_length = log2(intact_GDT_081421_rbp4$Length)
intact_GDT_081421_rbp4$stroud_8wk_mCA = as.numeric(as.character(intact_GDT_081421_rbp4$stroud_8wk_mCA))

intact_GDT_081421_rbp4_m = merge(intact_GDT_081421_rbp4, mL5_ensgene_mm9_sum, by="Gene")
intact_GDT_081421_rbp4_m$luo_L5_mCA = as.numeric(as.character(intact_GDT_081421_rbp4_m$mCA))
head(intact_GDT_081421_rbp4_m$luo_L5_mCA)
rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921$Gene = rownames(rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921)
rbp4_ko_exon_nondedup_coding_d3_dseq_rep_prefilt5_res_df_all_110921$Gene = rownames(rbp4_ko_exon_nondedup_coding_d3_dseq_rep_prefilt5_res_df_all_110921)
intact_GDT_081421_rbp4_m1 = merge(intact_GDT_081421_rbp4_m, rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921, by="Gene")

intact_GDT_081421_rbp4_m = merge(intact_GDT_081421_rbp4_m1, mL5_snmcseq_gene_and_flank_mC, by.x="Gene", by.y="gene")




setwd("~/Documents/Gabel_Lab/Analysis/Mecp2_cs_plots/deseq_outputs/sst_ko/nondedup/exon/deseq_tables/")

sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921 = read.table("sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv", header=TRUE)


intact_GDT_081421_sst = as.data.frame(cbind(Gene = intact_GDT_081421$Gene,
                                            Chr = intact_GDT_081421$chr,
                                            Start = intact_GDT_081421$start,
                                            Stop = intact_GDT_081421$stop,
                                            luo_n_sst_mCA = intact_GDT_081421$luo_n_sst_mCA,
                                            luo_2_sst_mCA = intact_GDT_081421$luo_2_sst_mCA,
                                            stroud_8wk_mCA = intact_GDT_081421$stroud_8wk_mCA, 
                                            pv_kotg_enh_count = intact_GDT_081421$pv_kotg_enh_count,
                                            pv_total_enh_count = intact_GDT_081421$pv_total_enh_count,
                                            li_atac_SST_ccre_count = intact_GDT_081421$li_atac_SST_ccre_count,
                                            li_atac_union_ccre_nonTSS_count = intact_GDT_081421$li_atac_union_ccre_nonTSS_count,
                                            li_atac_union_ccre_all_count = intact_GDT_081421$li_atac_union_ccre_all_count,
                                            enh_count =intact_GDT_081421$enh_count,
                                            Length = intact_GDT_081421$Length,
                                            MR_meta = intact_GDT_081421$MR_meta))
intact_GDT_081421_sst$Length = as.numeric(as.character(intact_GDT_081421_sst$Length))
intact_GDT_081421_sst$log2_length = log2(intact_GDT_081421_sst$Length)
intact_GDT_081421_sst$mCA = as.numeric(as.character(intact_GDT_081421_sst$luo_n_sst_mCA))
intact_GDT_081421_sst$stroud_8wk_mCA = as.numeric(as.character(intact_GDT_081421_sst$stroud_8wk_mCA))

sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921$Gene = row.names(sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921)
intact_GDT_081421_sst_m = merge(intact_GDT_081421_sst,sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921, by="Gene")

intact_GDT_081421_sst_m = merge(intact_GDT_081421_sst, mSst_all_snmcseq_gene_and_flank_mC, by.x="Gene", by.y="gene")





setwd("~/Documents/Gabel_Lab/Analysis/Mecp2_CS_plots/deseq_outputs/nr5a1_ko/nondedup/deseq_tables/")
nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921= read.table("nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv", header=TRUE)


intact_GDT_081421_nr5a1 = as.data.frame(cbind(Gene = intact_GDT_081421$Gene,
                                              Chr = intact_GDT_081421$chr,
                                              Start = intact_GDT_081421$start,
                                              Stop = intact_GDT_081421$stop,
                                              luo_L4_mCA = intact_GDT_081421$luo_L4_mCA,
                                              stroud_8wk_mCA = intact_GDT_081421$stroud_8wk_mCA, 
                                              pv_kotg_enh_count = intact_GDT_081421$pv_kotg_enh_count,
                                              pv_total_enh_count = intact_GDT_081421$pv_total_enh_count,
                                              li_atac_ITL4_cCRE_count = intact_GDT_081421$li_atac_ITL4_cCRE_count,
                                              li_atac_union_ccre_nonTSS_count = intact_GDT_081421$li_atac_union_ccre_nonTSS_count,
                                              li_atac_union_ccre_all_count = intact_GDT_081421$li_atac_union_ccre_all_count,
                                              enh_count =intact_GDT_081421$enh_count,
                                              Length = intact_GDT_081421$Length,
                                              MR_meta = intact_GDT_081421$MR_meta))
intact_GDT_081421_nr5a1$Length = as.numeric(as.character(intact_GDT_081421_nr5a1$Length))
intact_GDT_081421_nr5a1$log2_length = log2(intact_GDT_081421_nr5a1$Length)
intact_GDT_081421_nr5a1$mCA = as.numeric(as.character(intact_GDT_081421_nr5a1$luo_L4_mCA))


nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921$Gene = rownames(nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921)
intact_GDT_081421_nr5a1_m = merge(intact_GDT_081421_nr5a1, nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921, by="Gene")
intact_GDT_081421_nr5a1_m = merge(intact_GDT_081421_nr5a1_m, mL4_snmcseq_gene_and_flank_mC, by="Gene", by.y="gene")





intact_GDT_081421_nr5a1_m_g100 = intact_GDT_081421_nr5a1_m[which(intact_GDT_081421_nr5a1_m$Length > 100000),]
intact_GDT_081421_sst_m_g100 = intact_GDT_081421_sst_m[which(intact_GDT_081421_sst_m$Length > 100000),]
intact_GDT_081421_pv_m_g100 = intact_GDT_081421_pv_m[which(intact_GDT_081421_pv_m$Length > 100000),]
intact_GDT_081421_rbp4_m_g100 = intact_GDT_081421_rbp4_m1[which(intact_GDT_081421_rbp4_m1$Length > 100000),]

lm_nr5a1 = lm(intact_GDT_081421_nr5a1_m_g100$luo_L4_mCA ~ intact_GDT_081421_nr5a1_m_g100$ashr_log2FoldChange)
summary(lm_nr5a1)

lm(intact_GDT_081421_sst_m_g100$luo_n_sst_mCA ~ intact_GDT_081421_sst_m_g100$ashr_log2FoldChange)
lm_sst = lm(intact_GDT_081421_sst_m_g100$luo_2_sst_mCA ~ intact_GDT_081421_sst_m_g100$ashr_log2FoldChange)
summary(lm_sst)

lm_pv = lm(intact_GDT_081421_pv_m_g100$luo_pv_snmc_mCA ~ intact_GDT_081421_pv_m_g100$ashr_log2FoldChange)
lm(intact_GDT_081421_pv_m_g100$luo_pv_n_mCA ~ intact_GDT_081421_pv_m_g100$ashr_log2FoldChange)
summary(lm_pv)

lm_rbp4 = lm(intact_GDT_081421_rbp4_m_g100$luo_L5_mCA ~ intact_GDT_081421_rbp4_m_g100$ashr_log2FoldChange)
summary(lm_rbp4)

long_lm_coef = c(.01452, .03814, 0.05230, .07010)
long_lm_stderr = c(0.0052805, 0.0043776, 0.0037086, 0.0048672)

type = c("L4", "L5.", "PV", "SST")
type_1 = c("Layer 4", "Layer 5", "PV", "SST")
type_2 = c("Nr5a1", "Rbp4","PV", "SST")

df_lm = as.data.frame(long_lm_coef)
df_lm_t = as.data.frame(t(df_lm))
colnames(df_lm_t) = type
head(df_lm_t)

col_og_1 = c("skyblue", "orchid", "forestgreen", "darkorange")
col_og_2 = c("cornflowerblue", "mediumpurple3", "forestgreen", "darkorange")
col_og_2 = c("cornflowerblue", "#8F00FF", "seagreen", "coral1")
col_npg = c("#4DBBD5FF", "#8491B4FF","#00A087FF", "#F39B7FFF")
par(mar=c(10,9,5,4), mgp = c(5, 1, 0), lty=0)

slope_ylab= bquote(atop("Slope linear model mCA","by FC for Genes >100kb"))
dev.off()
long_lm_coef_bplot = barplot(long_lm_coef,
        names=type,
        col=col_npg,
        ylab= slope_ylab,
        cex.names=1.5,las=2,ylim=c(0,.08),
        cex.lab=1.5, cex.axis = 1.5)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

error.bar(long_lm_coef_bplot, long_lm_coef, long_lm_stderr)

intact_GDT_081421_pv_m$decile_luo_pv_snmc_mCA= ntile(intact_GDT_081421_pv_m$luo_pv_snmc_mCA, 10)
intact_GDT_081421_sst_m$decile_luo_2_sst_mCA= ntile(intact_GDT_081421_sst_m$luo_2_sst_mCA, 10)
intact_GDT_081421_nr5a1_m$decile_luo_L4_mCA= ntile(intact_GDT_081421_nr5a1_m$luo_L4_mCA, 10)
intact_GDT_081421_rbp4_m1$decile_luo_L5_mCA= ntile(intact_GDT_081421_rbp4_m1$luo_L5_mCA, 10)



intact_GDT_081421_pv_m_g100_d10_mca = intact_GDT_081421_pv_m[which(intact_GDT_081421_pv_m$Length > 100000 & intact_GDT_081421_pv_m$decile_luo_pv_snmc_mCA==10),]
intact_GDT_081421_sst_m_g100_d10_mca = intact_GDT_081421_sst_m[which(intact_GDT_081421_sst_m$Length > 100000 & intact_GDT_081421_sst_m$decile_luo_2_sst_mCA==10),]
intact_GDT_081421_nr5a1_m_g100_d10_mca = intact_GDT_081421_nr5a1_m[which(intact_GDT_081421_nr5a1_m$Length > 100000 & intact_GDT_081421_nr5a1_m$decile_luo_L4_mCA==10),]
intact_GDT_081421_rbp4_m_g100_d10_mca = intact_GDT_081421_rbp4_m1[which(intact_GDT_081421_rbp4_m1$Length > 100000 & intact_GDT_081421_rbp4_m1$decile_luo_L5_mCA==10),]


mean(intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange)
mean(intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange)
mean(intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)
mean(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange)

std <- function(x) sd(x)/sqrt(length(x))

std(intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange)
std(intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange)
std(intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)
std(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange)

t.test(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)
t.test(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange)
t.test(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange)
t.test(intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)
t.test(intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange)
t.test(intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange, intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)


mean_fc_long_high_meth = c(0.02181681, 0.04942076, 0.06938747, 0.1094242)

stderr_fc_l_high_meth = c(.004669986, .007273569, 0.008892548, .009342636)

slope_ylab= bquote(atop("Mean log2FC genes>100kb","and top decile mCA"))

par(mar=c(12,10,5,4), mgp = c(5, 1, 0), lty=0)
barplot(mean_fc_long_high_meth,
        names=type,border="white",
        col=col_og_2,
        ylab= slope_ylab,
        cex.names=2,las=2,ylim=c(0,.12),
        cex.lab=2, cex.axis = 2)

bplot_fc_lm = barplot(mean_fc_long_high_meth,
        names=type,border="white",
        col=col_npg,
        ylab= slope_ylab,
        cex.names=1.5,las=2,ylim=c(0,.12),
        cex.lab=1.5, cex.axis = 1.5)

error.bar(bplot_fc_lm, mean_fc_long_high_meth, stderr_fc_l_high_meth)

intact_GDT_081421_nr5a1$luo_L4_mCA = as.numeric(as.character(intact_GDT_081421_nr5a1$luo_L4_mCA))
intact_GDT_081421_pv$luo_pv_snmc_mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_snmc_mCA))
intact_GDT_081421_pv$luo_pv_n_mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_n_mCA))
intact_GDT_081421_sst$luo_2_sst_mCA = as.numeric(as.character(intact_GDT_081421_sst$luo_2_sst_mCA))
intact_GDT_081421_sst$luo_n_sst_mCA = as.numeric(as.character(intact_GDT_081421_sst$luo_n_sst_mCA))

mL5_ensgene_mm9_sum$mCA = as.numeric(as.character(intact_GDT_081421_pv$luo_pv_snmc_mCA))


mean(intact_GDT_081421_nr5a1$luo_L4_mCA, na.rm = T)
mean(intact_GDT_081421_pv$luo_pv_snmc_mCA, na.rm = T)
mean(intact_GDT_081421_pv$luo_pv_n_mCA, na.rm = T)
mean(intact_GDT_081421_sst$luo_n_sst_mCA, na.rm = T)
mean(intact_GDT_081421_sst$luo_2_sst_mCA, na.rm = T)
mean(mL5_ensgene_mm9_sum$mCA, na.rm=T)

ylab_meth = "mean genomic mCA/CA"
sc_meth_cell_types = c(0.02802428, 0.04024396, 0.0446, .05311)
sc_meth_cell_types = c(0.0280238, 0.0402341, 0.04374827, 0.04604574)
par(mar=c(10,10,5,4), mgp = c(5, 1, 0))
barplot(sc_meth_cell_types,
                      names=type,border="white",lwd=2,
                      col=col_og_2,
                      ylab= ylab_meth,
                      cex.names=1.5,las=2,ylim=c(0,.05),
                      cex.lab=1.5, cex.axis = 1.5)

bplot_fc_lm = barplot(sc_meth_cell_types, lwd=2,
                      names=type,border=c("skyblue", "orchid", "seagreen", "darkorange"),
                      col="gray",
                      ylab= ylab_meth,
                      cex.names=1.5,las=2,ylim=c(0,.06),
                      cex.lab=1.5, cex.axis = 1.5)

mean(intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange)
mean(intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange)
mean(intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange)



library(ggplot2)

head(intact_GDT_081421_pv_m_g100_d10_mca)
intact_GDT_081421_rbp4_m_g100_d10_mca
intact_GDT_rbp4_long_d10_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_rbp4_m_g100_d10_mca$normal_contrast_log2FoldChange,
                                           "apeglm_gWK_log2FoldChange" = intact_GDT_081421_rbp4_m_g100_d10_mca$apeglm_gWK_log2FoldChange, 
                                        "ashr_log2FoldChange" = intact_GDT_081421_rbp4_m_g100_d10_mca$ashr_log2FoldChange,
                                        "mean_cs_mCA" = 0.0402341,
                                        "no_sig_genes" =249,
                                        "no_MR_genes" = 142,
                                        "Cell_type" = "L5"))
intact_GDT_nr5a1_long_d10_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_nr5a1_m_g100_d10_mca$normal_contrast_log2FoldChange,
                                           "apeglm_gWK_log2FoldChange" = intact_GDT_081421_nr5a1_m_g100_d10_mca$apeglm_gWK_log2FoldChange, 
                                           "ashr_log2FoldChange" = intact_GDT_081421_nr5a1_m_g100_d10_mca$ashr_log2FoldChange,
                                           "mean_cs_mCA" = 0.0280238,
                                           "no_sig_genes" =95,
                                           "no_MR_genes" = 53,
                                           "Cell_type" = "L4"))
intact_GDT_pv_long_d10_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_pv_m_g100_d10_mca$normal_contrast_log2FoldChange,
                                           "apeglm_gWK_log2FoldChange" = intact_GDT_081421_pv_m_g100_d10_mca$apeglm_gWK_log2FoldChange, 
                                           "ashr_log2FoldChange" = intact_GDT_081421_pv_m_g100_d10_mca$ashr_log2FoldChange,
                                             "mean_cs_mCA" = 0.04374827,
                                           "no_sig_genes" = 431,
                                           "no_MR_genes" = 234,
                                           "Cell_type" = "PV"))
intact_GDT_sst_long_d10_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_sst_m_g100_d10_mca$normal_contrast_log2FoldChange,
                                           "apeglm_gWK_log2FoldChange" = intact_GDT_081421_sst_m_g100_d10_mca$apeglm_gWK_log2FoldChange, 
                                           "ashr_log2FoldChange" = intact_GDT_081421_sst_m_g100_d10_mca$ashr_log2FoldChange,
                                           "mean_cs_mCA" = 0.04604574,
                                           "no_sig_genes" = 568,
                                           "no_MR_genes" = 377,
                                           "Cell_type" = "SST"))
  
fc_df_full = as.data.frame(rbind(intact_GDT_nr5a1_long_d10_fc_s,
                   intact_GDT_rbp4_long_d10_fc_s,
                   intact_GDT_pv_long_d10_fc_s,
                   intact_GDT_sst_long_d10_fc_s))
head(fc_df_full)
table(fc_df_full$Cell_type)
fc_df_full$Cell_type = as.factor(fc_df_full$Cell_type)

fc_df_full$mean_cs_mCA = as.numeric(as.character(fc_df_full$mean_cs_mCA))
fc_df_full$no_sig_genes = as.numeric(as.character(fc_df_full$no_sig_genes))
fc_df_full$ashr_log2FoldChange = as.numeric(as.character(fc_df_full$ashr_log2FoldChange))
fc_df_full$apeglm_gWK_log2FoldChange = as.numeric(as.character(fc_df_full$apeglm_gWK_log2FoldChange))
fc_df_full$normal_contrast_log2FoldChange = as.numeric(as.character(fc_df_full$normal_contrast_log2FoldChange))
fc_df_full$ashr_log2FoldChange = as.numeric(as.character(fc_df_full$ashr_log2FoldChange))




setwd("~/Documents/Documents_New/INTACT_RNA_Analysis/")
p = ggplot(fc_df_full, aes(x=Cell_type, y=apeglm_gWK_log2FoldChange))+
  geom_violin(trim=FALSE)
p = p + stat_summary(fun.data=mean_sdl, mult=1, 
                 geom="pointrange", color="red")+ 
  stat_summary(fun.data = mean_se, geom = "errorbar")
p

p = ggplot(fc_df_full, aes(x=Cell_type, y=normal_contrast_log2FoldChange))+
  geom_point(position = position_jitter(width = 0.1), 
             alpha = .1)
p+stat_summary(fun = mean, geom = "point", shape="diamond", color=col_og_2,size=5) +
  stat_summary(fun.data = mean_se, color=col_og_2, geom = "errorbar")+theme_classic()

p = ggplot(fc_df_full, aes(x=mean_cs_mCA, y=ashr_log2FoldChange, color=Cell_type))+
  geom_point(stat="summary", fun="mean", shape="diamond", size=12)+
  scale_color_manual(values=col_og_2)
p  =p+stat_summary(fun.data = mean_se, geom = "errorbar")+
  theme_classic()
p =p+coord_cartesian(ylim=c(0,0.15), xlim=c(0,0.050))

ggsave(filename = "cs_mCA_ashr_fc_l100_d10_corr_0_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")

p = p + theme_classic()+ coord_cartesian(ylim=c(0,0.15), xlim=c(0.025,0.050))

ggsave(filename = "cs_mCA_ashr_fc_corr_l100_d10_unscal_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")


p = ggplot(fc_df_full, aes(x=mean_cs_mCA, y=no_sig_genes, color=Cell_type))+
  geom_point(stat="summary", fun="mean", shape="diamond", size=12)+
  scale_color_manual(values=col_og_2)
p + theme_classic()+ coord_cartesian(ylim=c(0,600), xlim=c(0,0.050))
p = p + theme_classic()+ coord_cartesian(ylim=c(0,600), xlim=c(0,0.050))


ggsave(filename = "cs_mCA_sig_gene_count_corr_0_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")

p = p + theme_classic()+ coord_cartesian(ylim=c(0,600), xlim=c(0.025,0.050))

ggsave(filename = "cs_mCA_sig_gene_count_corr_unscal_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")







intact_GDT_081421_nr5a1_mr_meta = intact_GDT_081421_nr5a1_m[which(intact_GDT_081421_nr5a1_m$MR_meta == TRUE),]
intact_GDT_081421_sst_mr_meta = intact_GDT_081421_sst_m[which(intact_GDT_081421_sst_m$MR_meta == TRUE),]
intact_GDT_081421_pv_mr_meta = intact_GDT_081421_pv_m[which(intact_GDT_081421_pv_m$MR_meta == TRUE),]
intact_GDT_081421_rbp4_mr_meta = intact_GDT_081421_rbp4_m1[which(intact_GDT_081421_rbp4_m1$MR_meta == TRUE),]


mean(intact_GDT_081421_nr5a1_mr_meta$ashr_log2FoldChange)
mean(intact_GDT_081421_rbp4_mr_meta$ashr_log2FoldChange)
mean(intact_GDT_081421_pv_mr_meta$ashr_log2FoldChange)
mean(intact_GDT_081421_sst_mr_meta$ashr_log2FoldChange)


intact_GDT_rbp4_meta_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_rbp4_mr_meta$normal_contrast_log2FoldChange,
                                                    "apeglm_gWK_log2FoldChange" = intact_GDT_081421_rbp4_mr_meta$apeglm_gWK_log2FoldChange, 
                                                    "ashr_log2FoldChange" = intact_GDT_081421_rbp4_mr_meta$ashr_log2FoldChange,
                                                    "mean_cs_mCA" = 0.0402341,
                                                    "Cell_type" = "L5"))
intact_GDT_nr5a1_meta_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_nr5a1_mr_meta$normal_contrast_log2FoldChange,
                                                     "apeglm_gWK_log2FoldChange" = intact_GDT_081421_nr5a1_mr_meta$apeglm_gWK_log2FoldChange, 
                                                     "ashr_log2FoldChange" = intact_GDT_081421_nr5a1_mr_meta$ashr_log2FoldChange,
                                                     "mean_cs_mCA" = 0.0280238,
                                                     "Cell_type" = "L4"))
intact_GDT_pv_meta_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_pv_mr_meta$normal_contrast_log2FoldChange,
                                                  "apeglm_gWK_log2FoldChange" = intact_GDT_081421_pv_mr_meta$apeglm_gWK_log2FoldChange, 
                                                  "ashr_log2FoldChange" = intact_GDT_081421_pv_mr_meta$ashr_log2FoldChange,
                                                  "mean_cs_mCA" = 0.04374827,
                                                  "Cell_type" = "PV"))
intact_GDT_sst_meta_fc_s = as.data.frame(cbind("normal_contrast_log2FoldChange" = intact_GDT_081421_sst_mr_meta$normal_contrast_log2FoldChange,
                                                   "apeglm_gWK_log2FoldChange" = intact_GDT_081421_sst_mr_meta$apeglm_gWK_log2FoldChange, 
                                                   "ashr_log2FoldChange" = intact_GDT_081421_sst_mr_meta$ashr_log2FoldChange,
                                                   "mean_cs_mCA" = 0.04604574,
                                                   "Cell_type" = "SST"))

fc_df_meta = as.data.frame(rbind(intact_GDT_nr5a1_meta_fc_s,
                                 intact_GDT_rbp4_meta_fc_s,
                                 intact_GDT_pv_meta_fc_s,
                                 intact_GDT_sst_meta_fc_s))
head(fc_df_full)
table(fc_df_meta$Cell_type)
fc_df_meta$Cell_type = as.factor(fc_df_meta$Cell_type)

fc_df_meta$mean_cs_mCA = as.numeric(as.character(fc_df_meta$mean_cs_mCA))
fc_df_meta$no_sig_genes = as.numeric(as.character(fc_df_meta$no_sig_genes))
fc_df_meta$ashr_log2FoldChange = as.numeric(as.character(fc_df_meta$ashr_log2FoldChange))
fc_df_meta$apeglm_gWK_log2FoldChange = as.numeric(as.character(fc_df_meta$apeglm_gWK_log2FoldChange))
fc_df_meta$normal_contrast_log2FoldChange = as.numeric(as.character(fc_df_meta$normal_contrast_log2FoldChange))
fc_df_meta$ashr_log2FoldChange = as.numeric(as.character(fc_df_meta$ashr_log2FoldChange))




setwd("~/Documents/Documents_New/INTACT_RNA_Analysis/")
p = ggplot(fc_df_meta, aes(x=Cell_type, y=apeglm_gWK_log2FoldChange))+
  geom_violin(trim=FALSE)
p = p + stat_summary(fun.data=mean_sdl, mult=1, 
                     geom="pointrange", color="red")+ 
  stat_summary(fun.data = mean_se, geom = "errorbar")
p

p = ggplot(fc_df_meta, aes(x=Cell_type, y=normal_contrast_log2FoldChange))+
  geom_point(position = position_jitter(width = 0.1), 
             alpha = .1)
p+stat_summary(fun = mean, geom = "point", shape="diamond", color=col_og_2,size=5) +
  stat_summary(fun.data = mean_se, color=col_og_2, geom = "errorbar")+theme_classic()

p = ggplot(fc_df_meta, aes(x=mean_cs_mCA, y=ashr_log2FoldChange, color=Cell_type))+
  geom_point(stat="summary", fun="mean", shape="diamond", size=12)+
  scale_color_manual(values=col_og_2)
p=p+stat_summary(fun.data = mean_se, geom = "errorbar")+
  theme_classic()
p= p+coord_cartesian(ylim=c(0,0.18), xlim=c(0,0.050))

ggsave(filename = "cs_mCA_ashr_fc_meta_corr_0_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")

p = p + coord_cartesian(ylim=c(0,0.18), xlim=c(0.025,0.050))

ggsave(filename = "cs_mCA_ashr_fc_meta_corr_unscal_thin.eps",
       plot= p, device= "eps", dpi=1200, width=15, height=20, units="cm")


