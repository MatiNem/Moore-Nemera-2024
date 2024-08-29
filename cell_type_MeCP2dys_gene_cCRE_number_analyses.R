library(data.table)
library(dplyr)
library(ggplot2)

#gene lists
pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")


sst_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
sst_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

L4_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L4_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

L5_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L5_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")


#intragenic non-promoter cCREs
intragenic_nonPromoter_cCREs_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_ensgene_mm9.txt")

#all cCREs linked to MeCP2-regulated genes
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")

Sst_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Sst_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Sst_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")

L4_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L4_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L4_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L4_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")

L5_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L5_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L5_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L5_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9.bed")

#table of numbers of cCREs associated with genes
coding_genes_mm9_with_UnioncCRE_numbers = fread("HG_lab/Mati/GabelLab/genesets/coding_genes_union_nonpromoter_cCRE_numbers_H3K27ac_data_table.txt")

mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")

#cCREs
union_cCREs_mPv_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mPv_Luo2017_snmcseq_CA_merged_mm9.bed")
union_cCREs_mPv_mCA[, cCRE_methylation := V5/V6]
union_cCREs_mSst_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mSst_all_Luo2017_snmcseq_CA_merged_mm9.bed")
union_cCREs_mSst_mCA[, cCRE_methylation := V5/V6]
union_cCREs_mL4_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mL4_Luo2017_snmcseq_CA_merged_mm9.bed")
union_cCREs_mL4_mCA[, cCRE_methylation := V5/V6]
union_cCREs_mL5_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mL5_all_Luo2017_snmcseq_CA_merged_mm9.bed")
union_cCREs_mL5_mCA[, cCRE_methylation := V5/V6]

names(union_cCREs_mPv_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "Pv_methylation")
names(union_cCREs_mSst_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "Sst_methylation")
names(union_cCREs_mL4_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "L4_methylation")
names(union_cCREs_mL5_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "L5_methylation")

mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mCA =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans, y=union_cCREs_mPv_mCA[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, Pv_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label"), nomatch=NULL))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mCA =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mCA, y=union_cCREs_mSst_mCA[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, Sst_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label"), nomatch=NULL))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mCA =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mCA, y=union_cCREs_mL4_mCA[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, L4_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label"), nomatch=NULL))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mL5_mCA =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mCA, y=union_cCREs_mL5_mCA[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, L5_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label"), nomatch=NULL))
#write.table(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mL5_mCA, file="HG_lab/Mati/GabelLab/For_Russell/cell_type_cCRE_gene_linkages/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans_mPv_mSst_mL4_mL5_mCA.txt", quote=F, row.names=F, sep="\t")

#cCRE mCG
union_cCREs_mPv_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mPv_Luo2017_snmcseq_CG_merged_mm9.bed")
union_cCREs_mPv_mCG[, cCRE_methylation := V5/V6]
union_cCREs_mSst_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mSst_all_Luo2017_snmcseq_CG_merged_mm9.bed")
union_cCREs_mSst_mCG[, cCRE_methylation := V5/V6]
union_cCREs_mL4_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mL4_Luo2017_snmcseq_CG_merged_mm9.bed")
union_cCREs_mL4_mCG[, cCRE_methylation := V5/V6]
union_cCREs_mL5_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mL5_all_Luo2017_snmcseq_CG_merged_mm9.bed")
union_cCREs_mL5_mCG[, cCRE_methylation := V5/V6]

names(union_cCREs_mPv_mCG) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "Pv_methylation")
names(union_cCREs_mSst_mCG) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "Sst_methylation")
names(union_cCREs_mL4_mCG) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "L4_methylation")
names(union_cCREs_mL5_mCG) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "meth_reads", "cyto_reads", "L5_methylation")

mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mCG =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans, y=union_cCREs_mPv_mCG[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, Pv_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mCG =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mCG, y=union_cCREs_mSst_mCG[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, Sst_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mCG =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mCG, y=union_cCREs_mL4_mCG[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, L4_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")))
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mL5_mCG =  data.table(inner_join(x=mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mCG, y=union_cCREs_mL5_mCG[, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label, L5_methylation)], by=c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")))
#write.table(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans_mPv_mSst_mL4_mL5_mCG, file="HG_lab/Mati/GabelLab/For_Russell/cell_type_cCRE_gene_linkages/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans_mPv_mSst_mL4_mL5_mCG.txt", quote=F, row.names=F, sep="\t")

gene_allLinkedcCRE_counts = mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, .(count = .N), by = Gene]
names(gene_allLinkedcCRE_counts)[2] = "num_union_linked_cCREs"

gene_intragenic_noncognateLinked_cCRE_counts = mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene==0), .(count = .N), by = Gene]
names(gene_intragenic_noncognateLinked_cCRE_counts)[2] = "num_union_intragenic_noncognateLinked_cCREs"

coding_genes_mm9_with_UnioncCRE_numbers_total <- unique(data.table(left_join(x=coding_genes_mm9_with_UnioncCRE_numbers , y=gene_allLinkedcCRE_counts, by=c("Gene"), nomatch=NA)))
coding_genes_mm9_with_UnioncCRE_numbers_total[is.na(num_union_linked_cCREs), num_union_linked_cCREs := 0]

coding_genes_mm9_with_UnioncCRE_numbers_total <- unique(data.table(left_join(x=coding_genes_mm9_with_UnioncCRE_numbers_total, y=gene_intragenic_noncognateLinked_cCRE_counts, by=c("Gene"), nomatch=NA)))
coding_genes_mm9_with_UnioncCRE_numbers_total[is.na(num_union_intragenic_noncognateLinked_cCREs), num_union_intragenic_noncognateLinked_cCREs := 0]

cellTypeMR_genes_with_UnioncCRE_numbers = rbind(
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="PV MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="SST MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L4 MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L5 MR genes"))

cellTypeMR_genes_with_UnioncCRE_numbers = cellTypeMR_genes_with_UnioncCRE_numbers %>% mutate(geneset = factor(geneset, levels=c("L4 MR genes", "L5 MR genes", "PV MR genes", "SST MR genes")))


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers, aes(x = geneset, y = num_union_intragenic_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L5 MR genes"="orchid", "PV MR genes"="forestgreen", "SST MR genes"="darkorange")) +
  coord_cartesian(ylim=c(0,250))+
  ylab("Number of intragenic cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(cellTypeMR_genes_with_UnioncCRE_numbers, aes(x = geneset, y = num_union_intragenic_cognateLinked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L5 MR genes"="orchid", "PV MR genes"="forestgreen", "SST MR genes"="darkorange")) +
  coord_cartesian(ylim=c(0,45))+
  ylab("Number of intragenic cognate-linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(cellTypeMR_genes_with_UnioncCRE_numbers, aes(x = geneset, y = num_union_extragenic_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L5 MR genes"="orchid", "PV MR genes"="forestgreen", "SST MR genes"="darkorange")) +
  coord_cartesian(ylim=c(0,17))+
  ylab("Number of extragenic linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_extragenicLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_extragenicLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(cellTypeMR_genes_with_UnioncCRE_numbers, aes(x = geneset, y = num_union_linked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L5 MR genes"="orchid", "PV MR genes"="forestgreen", "SST MR genes"="darkorange")) +
  coord_cartesian(ylim=c(0,75))+
  ylab("Number of linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_allLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_allLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers, aes(x = geneset, y = num_union_intragenic_noncognateLinked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=FALSE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L5 MR genes"="orchid", "PV MR genes"="forestgreen", "SST MR genes"="darkorange")) +
  coord_cartesian(ylim=c(0,11))+
  ylab("Number of intragenic\n noncognate-linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_noncognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_noncognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


#
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_linked_union_nonPromoter_cCREcoords_mm9


pv_mr_genes_resamp = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
sst_mr_genes_resamp = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Sst_mr_genes_coding_prefilt5_nondedup_SstTPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
L4_mr_genes_resamp = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L4_mr_genes_coding_prefilt5_nondedup_L4TPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
L5_mr_genes_resamp = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L5_mr_genes_coding_prefilt5_nondedup_L5TPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1

pv_mr_genes_resamp2 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_nondedup/resamp2.txt", header=FALSE)$V1
sst_mr_genes_resamp2 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Sst_mr_genes_coding_prefilt5_nondedup_SstTPM_MeCP2KO_nondedup/resamp2.txt", header=FALSE)$V1
L4_mr_genes_resamp2 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L4_mr_genes_coding_prefilt5_nondedup_L4TPM_MeCP2KO_nondedup/resamp2.txt", header=FALSE)$V1
L5_mr_genes_resamp2 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L5_mr_genes_coding_prefilt5_nondedup_L5TPM_MeCP2KO_nondedup/resamp2.txt", header=FALSE)$V1


cellTypeMR_genes_with_UnioncCRE_numbers_withResamp = rbind(
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="PV MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="SST MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L4 MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L5 MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% pv_mr_genes_resamp,], geneset="PV resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% sst_mr_genes_resamp,], geneset="SST resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L4_mr_genes_resamp,], geneset="L4 resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L5_mr_genes_resamp,], geneset="L5 resamp"))

cellTypeMR_genes_with_UnioncCRE_numbers_withResamp = cellTypeMR_genes_with_UnioncCRE_numbers_withResamp %>% mutate(geneset = factor(geneset, levels=c("L4 MR genes", "L4 resamp", "L5 MR genes", "L5 resamp", "PV MR genes", "PV resamp", "SST MR genes", "SST resamp")))

ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = num_union_intragenic_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,250))+
  ylab("Number of intragenic cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


cellTypeMR_genes_with_UnioncCRE_numbers_withResamp2 = rbind(
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="PV MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="SST MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L4 MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4],], geneset="L5 MR genes"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% pv_mr_genes_resamp2,], geneset="PV resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% sst_mr_genes_resamp2,], geneset="SST resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L4_mr_genes_resamp2,], geneset="L4 resamp"),
  cbind(coding_genes_mm9_with_UnioncCRE_numbers_total[Gene %in% L5_mr_genes_resamp2,], geneset="L5 resamp"))

cellTypeMR_genes_with_UnioncCRE_numbers_withResamp2 = cellTypeMR_genes_with_UnioncCRE_numbers_withResamp2 %>% mutate(geneset = factor(geneset, levels=c("L4 MR genes", "L4 resamp", "L5 MR genes", "L5 resamp", "PV MR genes", "PV resamp", "SST MR genes", "SST resamp")))

ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp2, aes(x = geneset, y = num_union_intragenic_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,250))+
  ylab("Number of intragenic cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))



ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = num_union_linked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,75))+
  ylab("Number of linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_allLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_allLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = num_union_intragenic_cognateLinked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,45))+
  ylab("Number of intragenic cognate-linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_cognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = num_union_extragenic_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,17))+
  ylab("Number of extragenic linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_extragenicLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_extragenicLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = num_union_intragenic_noncognateLinked_cCREs, fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(0,13))+
  ylab("Number of intragenic\n noncogate-linked cCREs") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_noncognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/enhancer_number_plots/cell_confusion_paper_number_plots/number_intragenic_noncognateLinked_cCREs_in_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp, aes(x = geneset, y = log10(gene_end - gene_start), fill=geneset))+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape=NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("L4 MR genes"="cyan", "L4 resamp"="gray", "L5 MR genes"="orchid", "L5 resamp"="gray", "PV MR genes"="forestgreen", "PV resamp"="gray", "SST MR genes"="darkorange", "SST resamp"="gray")) +
  coord_cartesian(ylim=c(2.5,6.5))+
  ylab("Log10 gene length") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/gene_body_plots/Log10GeneLength_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_body_plots/Log10GeneLength_cellTypeMR_genes_q0.1_coding_prefilt5_nondedup_boxplot_withResamp.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L4 MR genes", num_union_intragenic_cCREs], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L4 resamp", num_union_intragenic_cCREs])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L5 MR genes", num_union_intragenic_cCREs], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L5 resamp", num_union_intragenic_cCREs])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="PV MR genes", num_union_intragenic_cCREs], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="PV resamp", num_union_intragenic_cCREs])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="SST MR genes", num_union_intragenic_cCREs], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="SST resamp", num_union_intragenic_cCREs])$p.value

wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L4 MR genes", log10(gene_end - gene_start)], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L4 resamp", log10(gene_end - gene_start)])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L5 MR genes", log10(gene_end - gene_start)], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="L5 resamp", log10(gene_end - gene_start)])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="PV MR genes", log10(gene_end - gene_start)], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="PV resamp", log10(gene_end - gene_start)])$p.value
wilcox.test(cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="SST MR genes", log10(gene_end - gene_start)], cellTypeMR_genes_with_UnioncCRE_numbers_withResamp[geneset=="SST resamp", log10(gene_end - gene_start)])$p.value
