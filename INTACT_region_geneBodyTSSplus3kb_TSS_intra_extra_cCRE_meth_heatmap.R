library(data.table)
library(dplyr)
library(ggplot2)

pval_fc_heatmap <- function(pval_matrix, fc_matrix, column_names, color, pval_max, title, x_size=18, y_size=18, x_angle=0, y_angle=0, tile_text_size=6){
  test <- as.data.frame(pval_matrix)
  test_fc <- as.data.frame(fc_matrix)
  colnames(test) = column_names
  colnames(test_fc) = column_names
  test$cellType = row.names(test)
  test_fc$cellType = row.names(test_fc)
  #Melting data so we can plot it with GGplot
  test.m <- suppressWarnings(melt(test,id.vars = c("cellType")))
  test_fc.m <- suppressWarnings(melt(test_fc,id.vars = c("cellType")))
  #Resetting factors
  test.m = test.m %>% mutate(cellType = factor(cellType, levels=c("L5", "L4", "Sst", "Pv")))
  test_fc.m = test_fc.m %>% mutate(cellType = factor(cellType, levels=c("L5", "L4", "Sst", "Pv")))
  test_pval_and_fc = data.table(cbind(test.m, fc = test_fc.m$value))
  test_pval_and_fc$pvalMax = test_pval_and_fc$value
  test_pval_and_fc[value > pval_max, pvalMax := pval_max]
  
  #Creating the plot itself
  ggplot(test_pval_and_fc,aes(variable,cellType)) + geom_tile(aes(fill=pvalMax),color = "white") +
    ggtitle(title)+
    #Creating legend
    guides(fill=guide_colorbar(title="-log10 p-value", title.position="top")) +
    #Creating color range
    scale_fill_gradientn(limits = c(0,pval_max), colors=c("white",color),guide="colorbar") +
    geom_text(aes(label = fc), color = "black", size = tile_text_size)+
    theme_bw()+
    #Rotating labels
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
          axis.text.x=element_text(size=x_size, angle=x_angle), axis.text.y=element_text(size=y_size, angle=y_angle), axis.ticks = element_blank(), axis.title=element_blank())
}

cellType_elem_fc_cCREtypes_TSS <- function(data_table, geneset_numerator, geneset_denominator){
  Pv = c(
    data_table[(element == "Region") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "Region") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Gene body") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "Gene body") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val],
    data_table[(element == "TSS") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "TSS") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val],
    data_table[(element == "All cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "All cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Intragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "Intragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Extragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), median.val] / data_table[(element == "Extragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), median.val])
  Sst = c(
    data_table[(element == "Region") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "Region") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Gene body") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "Gene body") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val],
    data_table[(element == "TSS") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "TSS") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val],
    data_table[(element == "All cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "All cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Intragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "Intragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Extragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), median.val] / data_table[(element == "Extragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), median.val])
  L4 = c(
    data_table[(element == "Region") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "Region") & (cell_type=="L4") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Gene body") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "Gene body") & (cell_type=="L4") & (geneset==geneset_denominator), median.val],
    data_table[(element == "TSS") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "TSS") & (cell_type=="L4") & (geneset==geneset_denominator), median.val],
    data_table[(element == "All cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "All cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Intragenic cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "Intragenic cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Extragenic cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), median.val] / data_table[(element == "Extragenic cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), median.val])
  L5 = c(
    data_table[(element == "Region") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "Region") & (cell_type=="L5") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Gene body") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "Gene body") & (cell_type=="L5") & (geneset==geneset_denominator), median.val],
    data_table[(element == "TSS") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "TSS") & (cell_type=="L5") & (geneset==geneset_denominator), median.val],
    data_table[(element == "All cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "All cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Intragenic cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "Intragenic cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), median.val],
    data_table[(element == "Extragenic cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), median.val] / data_table[(element == "Extragenic cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), median.val])
  return(rbind(Pv, Sst, L4, L5))
}

cellType_elem_pvals_cCREtypes_TSS <- function(data_table, geneset_numerator, geneset_denominator){
  Pv = c(
    wilcox.test(data_table[(element == "Region") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "Region") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Gene body") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "Gene body") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "TSS") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "TSS") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "All cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "All cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Intragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "Intragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Extragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_numerator), val], data_table[(element == "Extragenic cCREs") & (cell_type=="Pv") & (geneset==geneset_denominator), val])$p.value)
  Sst = c(
    wilcox.test(data_table[(element == "Region") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "Region") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Gene body") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "Gene body") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "TSS") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "TSS") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "All cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "All cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Intragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "Intragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Extragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_numerator), val], data_table[(element == "Extragenic cCREs") & (cell_type=="Sst") & (geneset==geneset_denominator), val])$p.value)
  L4 = c(
    wilcox.test(data_table[(element == "Region") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "Region") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Gene body") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "Gene body") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "TSS") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "TSS") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "All cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "All cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Intragenic cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "Intragenic cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Extragenic cCREs") & (cell_type=="L4") & (geneset==geneset_numerator), val], data_table[(element == "Extragenic cCREs") & (cell_type=="L4") & (geneset==geneset_denominator), val])$p.value)
  L5 = c(
    wilcox.test(data_table[(element == "Region") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "Region") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Gene body") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "Gene body") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "TSS") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "TSS") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "All cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "All cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Intragenic cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "Intragenic cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value,
    wilcox.test(data_table[(element == "Extragenic cCREs") & (cell_type=="L5") & (geneset==geneset_numerator), val], data_table[(element == "Extragenic cCREs") & (cell_type=="L5") & (geneset==geneset_denominator), val])$p.value)
  return(rbind(Pv, Sst, L4, L5))
}
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

#meta-MeCP2-regulated genes
meta_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
meta_MA_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MA_genes_geneColumn4_mm9.bed")
meta_unchanged_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/MeCP2_unchanged_meta_genes_mm9.bed")

coding_genes = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')

non_metaMR_genes <- fread("HG_lab/Mati/GabelLab/genesets/meta_genes/nonMetaMR_genes_geneColumn4_mm9.bed")

#cCREs
mousebrain_union_cCREs_PV_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_PV_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_SST_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_SST_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L4_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_L4_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L5_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L5_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_L5_mCA[, cCRE_methylation := V5/V6]

mousebrain_union_cCREs_PV_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_PV_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_SST_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_SST_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L4_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_L4_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L5_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L5_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_L5_mCG[, cCRE_methylation := V5/V6]

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

mousebrain_union_cCREs_PV_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCA$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_SST_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_SST_mCA$cCRE_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L4_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_L4_mCA$cCRE_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L5_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_L5_mCA$cCRE_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

mousebrain_union_cCREs_PV_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCG$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_SST_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_SST_mCG$cCRE_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L4_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_L4_mCG$cCRE_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L5_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_L5_mCG$cCRE_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

#turning negative methylation values into zeros
mousebrain_union_cCREs_PV_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_SST_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L4_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L5_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]

mousebrain_union_cCREs_PV_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_SST_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L4_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L5_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]

#non-promoter cCREs
nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(Intragenic==1) & (Intragenic_to_linked_gene==1), ]
mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")[(Intragenic==0), ]



#union non-promoter cCREs linked to MR genes
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


#intragenic cCREs cognate-linked to MeCP2-regulated genes
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Sst_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Sst_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L4_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L4_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L5_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L5_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/intragenicLinked/L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

#extragenic cCREs
Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Pv_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Sst_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Sst_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Sst_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/Sst_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L4_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L4_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L4_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L4_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L5_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L5_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L5_otherCellType_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")
L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nondedup/extragenicLinked/L5_otherCellType_MA_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_union_nonPromoter_cCREcoords_mm9.bed")

#gene body and flank methylation tables, using genic methylation from TSS+3kb to TES
PV_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L4_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
SST_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L5_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")

PV_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
L4_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
SST_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
L5_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")

#replacing negative methylation values with zeros
PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
PV_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]

PV_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
PV_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]

#methylation of 1kb promoter windows centered at TSS
promoterWindows_PV_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
promoterWindows_PV_mCA[, promoter_methylation := (V7/V8)]
promoterWindows_SST_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
promoterWindows_SST_mCA[, promoter_methylation := (V7/V8)]
promoterWindows_L4_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
promoterWindows_L4_mCA[, promoter_methylation := (V7/V8)]
promoterWindows_L5_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_L5_WT_KO_deep_INTACT_mCA_mm9.bed")
promoterWindows_L5_mCA[, promoter_methylation := (V7/V8)]

promoterWindows_PV_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
promoterWindows_PV_mCG[, promoter_methylation := (V7/V8)]
promoterWindows_SST_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
promoterWindows_SST_mCG[, promoter_methylation := (V7/V8)]
promoterWindows_L4_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
promoterWindows_L4_mCG[, promoter_methylation := (V7/V8)]
promoterWindows_L5_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_1kb_promoterWindows_L5_WT_KO_deep_INTACT_mCG_mm9.bed")
promoterWindows_L5_mCG[, promoter_methylation := (V7/V8)]

#nonconversion rate correction
promoterWindows_PV_mCA$promoter_methylation_corrected <- promoterWindows_PV_mCA$promoter_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
promoterWindows_SST_mCA$promoter_methylation_corrected <- promoterWindows_SST_mCA$promoter_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
promoterWindows_L4_mCA$promoter_methylation_corrected <- promoterWindows_L4_mCA$promoter_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
promoterWindows_L5_mCA$promoter_methylation_corrected <- promoterWindows_L5_mCA$promoter_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

promoterWindows_PV_mCG$promoter_methylation_corrected <- promoterWindows_PV_mCG$promoter_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
promoterWindows_SST_mCG$promoter_methylation_corrected <- promoterWindows_SST_mCG$promoter_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
promoterWindows_L4_mCG$promoter_methylation_corrected <- promoterWindows_L4_mCG$promoter_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
promoterWindows_L5_mCG$promoter_methylation_corrected <- promoterWindows_L5_mCG$promoter_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

#turning negative methylation values into zeros
promoterWindows_PV_mCA[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_SST_mCA[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_L4_mCA[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_L5_mCA[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]

promoterWindows_PV_mCG[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_SST_mCG[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_L4_mCG[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]
promoterWindows_L5_mCG[promoter_methylation_corrected < 0, promoter_methylation_corrected := 0]

all_cells_metaMR_cCREtypes_TSS_mCA = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Meta-MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="All other genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Meta-MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="All other genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Meta-MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="All other genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Meta-MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="All other genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="All other genes")
)

all_cells_metaMR_cCREtypes_TSS_mCA_summary = all_cells_metaMR_cCREtypes_TSS_mCA %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_metaMR_cCREtypes_TSS_mCA_summary = all_cells_metaMR_cCREtypes_TSS_mCA_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_metaMR_cCREtypes_TSS_mCA_summary = all_cells_metaMR_cCREtypes_TSS_mCA_summary %>% mutate(geneset = factor(geneset, levels=c("Meta-MR genes", "All other genes")))
all_cells_metaMR_cCREtypes_TSS_mCA_summary = all_cells_metaMR_cCREtypes_TSS_mCA_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_metaMR_cCREtypes_TSS_mCA_summary_dt = data.table(all_cells_metaMR_cCREtypes_TSS_mCA_summary)


cellType_elem_metaMR_allOther_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_metaMR_cCREtypes_TSS_mCA_summary_dt, "Meta-MR genes", "All other genes"),1)
cellType_elem_metaMR_allOther_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_metaMR_cCREtypes_TSS_mCA, "Meta-MR genes", "All other genes"))

pval_fc_heatmap(cellType_elem_metaMR_allOther_cCREtypes_TSS_mCA_pvals, cellType_elem_metaMR_allOther_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 40, "Ratio mCA (Meta-MR)/(All other)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/metaMR_over_allOther_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/metaMR_over_allOther_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

all_cells_metaMR_cCREtypes_TSS_mCG = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Meta-MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="All other genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Meta-MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="All other genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Meta-MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="All other genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Meta-MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="All other genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% meta_MR_genes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Meta-MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% non_metaMR_genes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="All other genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% meta_MR_genes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Meta-MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% non_metaMR_genes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="All other genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="All other genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% meta_MR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Meta-MR genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% non_metaMR_genes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="All other genes")
)

all_cells_metaMR_cCREtypes_TSS_mCG_summary = all_cells_metaMR_cCREtypes_TSS_mCG %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_metaMR_cCREtypes_TSS_mCG_summary = all_cells_metaMR_cCREtypes_TSS_mCG_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_metaMR_cCREtypes_TSS_mCG_summary = all_cells_metaMR_cCREtypes_TSS_mCG_summary %>% mutate(geneset = factor(geneset, levels=c("Meta-MR genes", "All other genes")))
all_cells_metaMR_cCREtypes_TSS_mCG_summary = all_cells_metaMR_cCREtypes_TSS_mCG_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_metaMR_cCREtypes_TSS_mCG_summary_dt = data.table(all_cells_metaMR_cCREtypes_TSS_mCG_summary)


cellType_elem_metaMR_allOther_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_metaMR_cCREtypes_TSS_mCG_summary_dt, "Meta-MR genes", "All other genes"),1)
cellType_elem_metaMR_allOther_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_metaMR_cCREtypes_TSS_mCG, "Meta-MR genes", "All other genes"))

pval_fc_heatmap(cellType_elem_metaMR_allOther_cCREtypes_TSS_mCG_pvals, cellType_elem_metaMR_allOther_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 40, "Ratio mCG (Meta-MR)/(All other)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/metaMR_over_allOther_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/metaMR_over_allOther_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')



PvMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
PvUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
SstMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/SstMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
SstUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/SstUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L4MR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4MR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4Unchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4Unchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L5MR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5MR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5Unchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5Unchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")

all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% PvUnchanged_metaMRgenes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Non-Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% PvUnchanged_metaMRgenes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% PvUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% SstMR_metaMRgenes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% SstUnchanged_metaMRgenes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Non-Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% SstMR_metaMRgenes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% SstUnchanged_metaMRgenes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% SstMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% SstUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% SstMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% SstUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% SstMR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% SstUnchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4MR_metaMRgenes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4Unchanged_metaMRgenes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Non-Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4MR_metaMRgenes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4Unchanged_metaMRgenes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5MR_metaMRgenes[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5Unchanged_metaMRgenes[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Non-Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5MR_metaMRgenes[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5Unchanged_metaMRgenes[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5MR_metaMRgenes[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5Unchanged_metaMRgenes[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Non-Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% PvMR_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% PvUnchanged_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Non-Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% SstMR_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% SstUnchanged_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Non-Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4MR_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4Unchanged_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Non-Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5MR_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5Unchanged_metaMRgenes[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Non-Cell-type MR genes")
)

all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type MR genes", "Non-Cell-type MR genes")))
all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary_dt = data.table(all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary)

#element mCA/CA fold difference of non-cell-type meta MR over cell-type meta MR genes
cellType_elem_cellTypevsNonCellType_metaMR_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA_summary_dt, "Non-Cell-type MR genes", "Cell-type MR genes"),1)
cellType_elem_cellTypevsNonCellType_metaMR_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_cellTypeVsNonCellType_MetaMR_cCREtypes_TSS_mCA, "Non-Cell-type MR genes", "Cell-type MR genes"))

pval_fc_heatmap(cellType_elem_cellTypevsNonCellType_metaMR_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypevsNonCellType_metaMR_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "orchid", 20, "Ratio mCA (Non-cell-type meta-MR)/(Cell-type meta-MR)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/nonCellType_over_cellType_metaMR_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/nonCellType_over_cellType_metaMR_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#element mCA/CA of cell-type vs other-cell-type MR genes
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Other-cell-type MR genes")
)

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type MR genes", "Other-cell-type MR genes")))
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary_dt = data.table(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary)

#element mCA/CA fold difference of non-cell-type meta MR over cell-type meta MR genes
cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA_summary_dt, "Other-cell-type MR genes", "Cell-type MR genes"),1)
cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCA, "Other-cell-type MR genes", "Cell-type MR genes"))

pval_fc_heatmap(cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "orchid", 20, "Ratio mCA (Other-cell-type MR)/(Cell-type MR)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/otherCellType_over_cellType_MR_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/otherCellType_over_cellType_MR_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')


#element mCG/CG of cell-type vs other-cell-type MR genes
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Other-cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Other-cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_otherCellType_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Other-cell-type MR genes")
)

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type MR genes", "Other-cell-type MR genes")))
all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary = all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary_dt = data.table(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary)

#element mCG/CG fold difference of non-cell-type meta MR over cell-type meta MR genes
cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG_summary_dt, "Other-cell-type MR genes", "Cell-type MR genes"),1)
cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_cellTypeVsOtherCellType_MR_cCREtypes_TSS_mCG, "Other-cell-type MR genes", "Cell-type MR genes"))

pval_fc_heatmap(cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCG_pvals, cellType_elem_cellTypevsOtherCellType_MR_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "orchid", 20, "Ratio mCG (Other-cell-type MR)/(Cell-type MR)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/otherCellType_over_cellType_MR_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/otherCellType_over_cellType_MR_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')


#element mCA/CA of MeCP2-regulated genes
all_cells_MeCP2dys_cCREtypes_TSS_mCA = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
)

all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MR genes", "Cell-type MA genes")))
all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary)

#element mCA/CA fold difference of cell-type meta MR over unchanged genes
cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes"),1)
cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA, "Cell-type MR genes", "Cell-type unchanged genes"))

pval_fc_heatmap(cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MR)/(Unchanged)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchanged_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchanged_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#element mCA/CA fold difference of cell-type meta MA over unchanged genes
cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes"),1)
cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA, "Cell-type MA genes", "Cell-type unchanged genes"))

pval_fc_heatmap(cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MA)/(Unchanged)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchanged_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchanged_genes_INTACT_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#mCG
#element mCG/CG of MeCP2-regulated genes
all_cells_MeCP2dys_cCREtypes_TSS_mCG = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
)

all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes", "Cell-type MR genes", "Cell-type MA genes")))
all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary)

#element mCG/CG fold difference of cell-type meta MR over unchanged genes
cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes"),1)
cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG, "Cell-type MR genes", "Cell-type unchanged genes"))

pval_fc_heatmap(cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCG_pvals, cellType_elem_cellTypeMRvsUnchanged_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 20, "Ratio mCG (Cell-type MR)/(Unchanged)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchanged_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchanged_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#element mCG/CG fold difference of cell-type meta MA over unchanged genes
cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes"),1)
cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG, "Cell-type MA genes", "Cell-type unchanged genes"))

pval_fc_heatmap(cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCG_pvals, cellType_elem_cellTypeMAvsUnchanged_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 20, "Ratio mCG (Cell-type MA)/(Unchanged)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchanged_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchanged_genes_INTACT_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')


#barplot of cell-type vs other-cell-type MR gene element mCA/CA

all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE = data.table(rbind(
  cbind(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L4"),
  cbind(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L5"),
  cbind(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="Pv"),
  cbind(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="Sst"),
  cbind(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L4"),
  cbind(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L5"),
  cbind(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="Pv"),
  cbind(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="Sst"),
  cbind(val=median(mousebrain_union_cCREs_L4_mCA[(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L4"),
  cbind(val=median(mousebrain_union_cCREs_L5_mCA[(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L5"),
  cbind(val=median(mousebrain_union_cCREs_PV_mCA[(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE),  element="Intragenic cCREs", cell_type="Pv"),
  cbind(val=median(mousebrain_union_cCREs_SST_mCA[(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="Sst")
))


all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE = all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(element = factor(element, levels=c("Region", "Gene body", "Intragenic cCREs")))
all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE = all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(cell_type = factor(cell_type, levels=c("L4", "L5", "Pv", "Sst")))

ggplot(all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE[cell_type=="L5"] ,aes(x = element, y = as.numeric(val))) +
  geom_point(shape = 18, size = 10) +  # Diamond shape
  geom_text(stat='identity', aes(label=round(as.numeric(val), 3)), vjust=-0.5)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_mCAperCA_diamondplot.png', width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_mCAperCA_diamondplot.eps', width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE[cell_type=="Pv"] ,aes(x = element, y = as.numeric(val))) +
  geom_point(shape = 18, size = 10) +  # Diamond shape
  geom_text(stat='identity', aes(label=round(as.numeric(val), 3)), vjust=-0.5)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("PV INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_mCAperCA_diamondplot.png', width = 3.7, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_mCAperCA_diamondplot.eps', width = 3.7, height = 5, dpi = 300, units = "in", device='eps')


###resample
#Performs resampling and returns randomized sample of rownames with equivalent expression values - full_data - data to sample from, samp = sample to resample against (determines expression values to sample from and resample size). Samp format is a matrix or dataframe of expression values in the same format as Full data.. Full data is expected to be a matrix or dataframe of expression values, with rownames of gene IDs. Duplicate option tells script to look for duplicate rownames in the resampling, and repeats resampling to avoid duplicates (off by default, not recommended to use if data to sample is >20% of the "full data")
resample <- function(full_data, samp, dup=FALSE) {
  expr_range = 10;
  if(ncol(samp) == 1){
    samp_vals = samp
  }
  else{
    samp_vals = rowMeans(samp)
  }
  full_df = as.data.frame(full_data)
  if(ncol(full_df) == 1){
    full_vals = full_df
  }
  else{
    full_vals = rowMeans(full_df)
  }
  samp_names = rownames(samp)
  full_df$means = full_vals
  order_data = full_df[order(full_df$means),]
  order_data$index = 1:nrow(order_data)
  order_data = order_data[!(rownames(order_data) %in% samp_names),]
  #ind = match(samp,names(ranks))
  resamp = samp_names
  ind = unlist(lapply(lapply(samp_vals,">",order_data$means),sum))
  resamps = sample(c(-expr_range:expr_range,1),length(ind),replace = TRUE)
  #resamp_ind = lapply(ind, FUN = function(x) x+sample(-expr_range:-1,1:expr_range,1))
  resamp_ranks = ind + resamps
  resamp_ranks[resamp_ranks<1] = sample(1:10, sum(resamp_ranks<1))
  resamp_ranks[resamp_ranks>nrow(order_data)] = sample(nrow(order_data)-0:9, sum(resamp_ranks>nrow(order_data)))
  resamp = rownames(order_data)[resamp_ranks]
  run_tries = 0
  while(dup & sum(duplicated(resamp))>0 & run_tries < 100){
    resamps = sample(c(-expr_range:expr_range,1),length(ind),replace = TRUE)
    resamp_ranks = ind + resamps
    resamp_ranks[resamp_ranks<1] = sample(1:expr_range, sum(resamp_ranks<1))
    resamp_ranks[resamp_ranks>nrow(order_data)] = sample(nrow(order_data)-1:expr_range+1, sum(resamp_ranks>nrow(order_data)))
    resamp = rownames(order_data)[resamp_ranks]
    run_tries = run_tries+1
  }
  return(resamp)
}

#TPMs
Pv_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")
Sst_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_Mecp2KO_gene_TPMs_nondedup.txt")
L4_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_Mecp2KO_gene_TPMs_nondedup.txt")
L5_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_Mecp2KO_gene_TPMs_nondedup.txt")

Pv_TPM_df = data.frame(Pv_TPM[, .(Gene, Pv_WT_TPM_avg)], row.names="Gene")
Sst_TPM_df = data.frame(Sst_TPM[, .(Gene, Sst_WT_TPM_avg)], row.names="Gene")
L4_TPM_df = data.frame(L4_TPM[, .(Gene, L4_WT_TPM_avg)], row.names="Gene")
L5_TPM_df = data.frame(L5_TPM[, .(Gene, L5_WT_TPM_avg)], row.names="Gene")

Pv_TPM_df$Pv_WT_TPM_avg2 <- Pv_TPM_df$Pv_WT_TPM_avg
Sst_TPM_df$Sst_WT_TPM_avg2 <- Sst_TPM_df$Sst_WT_TPM_avg
L4_TPM_df$L4_WT_TPM_avg2 <- L4_TPM_df$L4_WT_TPM_avg
L5_TPM_df$L5_WT_TPM_avg2 <- L5_TPM_df$L5_WT_TPM_avg

#unchanged and Mecp2-regulated gene TPMs
Pv_TPM_df_5kb_unchanged <- Pv_TPM_df[pv_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],]
Pv_TPM_df_5kb_MR <- Pv_TPM_df[pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]
Pv_TPM_df_5kb_MA <- Pv_TPM_df[pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]

Sst_TPM_df_5kb_unchanged <- Sst_TPM_df[sst_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],]
Sst_TPM_df_5kb_MR <- Sst_TPM_df[sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]
Sst_TPM_df_5kb_MA <- Sst_TPM_df[sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]

L4_TPM_df_5kb_unchanged <- L4_TPM_df[L4_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],]
L4_TPM_df_5kb_MR <- L4_TPM_df[L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]
L4_TPM_df_5kb_MA <- L4_TPM_df[L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]

L5_TPM_df_5kb_unchanged <- L5_TPM_df[L5_unchanged_genes_p0.5_nondedup_mm9[V3-V2 >=5000,V4],]
L5_TPM_df_5kb_MR <- L5_TPM_df[L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]
L5_TPM_df_5kb_MA <- L5_TPM_df[L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4],]

pv_mr_resamp1 <- resample(Pv_TPM_df_5kb_unchanged, Pv_TPM_df_5kb_MR )

multi_resamp_func <- function(full_df, sample_df, num_resamp=100, save_resamp=FALSE, output_resamp=NULL){
  set.seed(1)
  #set number of resamplings
  #call matrix that you will populate with resamplings
  resamplings = matrix(nrow=nrow(sample_df),ncol=num_resamp)
  #for loop to populate the resampling matrix
  for(n in 1:num_resamp){
    resamplings[,n] = resample(full_df, sample_df)
  }
  if(save_resamp){
    write.table(resamplings, output_resamp, quote=F, row.names=F, col.names=F, sep="\t")
  }
  return(resamplings)
}

#PV
pv_mr_5kb_resamps <- multi_resamp_func(full_df=Pv_TPM_df_5kb_unchanged, sample_df=Pv_TPM_df_5kb_MR, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/Pv_mr_genes_coding_prefilt5_nondedup_min5kb_genes_PvTPM_100resamplings.txt")



pv_ma_5kb_resamps <- multi_resamp_func(full_df=Pv_TPM_df_5kb_unchanged, sample_df=Pv_TPM_df_5kb_MA, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/Pv_ma_genes_coding_prefilt5_nondedup_min5kb_genes_PvTPM_100resamplings.txt")

#SST
sst_mr_5kb_resamps <- multi_resamp_func(full_df=Sst_TPM_df_5kb_unchanged, sample_df=Sst_TPM_df_5kb_MR, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/Sst_mr_genes_coding_prefilt5_nondedup_min5kb_genes_SstTPM_100resamplings.txt")



sst_ma_5kb_resamps <- multi_resamp_func(full_df=Sst_TPM_df_5kb_unchanged, sample_df=Sst_TPM_df_5kb_MA, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/Sst_ma_genes_coding_prefilt5_nondedup_min5kb_genes_SstTPM_100resamplings.txt")

#L4
L4_mr_5kb_resamps <- multi_resamp_func(full_df=L4_TPM_df_5kb_unchanged, sample_df=L4_TPM_df_5kb_MR, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L4_mr_genes_coding_prefilt5_nondedup_min5kb_genes_L4TPM_100resamplings.txt")



L4_ma_5kb_resamps <- multi_resamp_func(full_df=L4_TPM_df_5kb_unchanged, sample_df=L4_TPM_df_5kb_MA, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L4_ma_genes_coding_prefilt5_nondedup_min5kb_genes_L4TPM_100resamplings.txt")

#L5
L5_mr_5kb_resamps <- multi_resamp_func(full_df=L5_TPM_df_5kb_unchanged, sample_df=L5_TPM_df_5kb_MR, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L5_mr_genes_coding_prefilt5_nondedup_min5kb_genes_L5TPM_100resamplings.txt")



L5_ma_5kb_resamps <- multi_resamp_func(full_df=L5_TPM_df_5kb_unchanged, sample_df=L5_TPM_df_5kb_MA, num_resamp=100, save_resamp=TRUE, 
                                       output_resamp="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L5_ma_genes_coding_prefilt5_nondedup_min5kb_genes_L5TPM_100resamplings.txt")

#read in the resampled genes
pv_mr_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/pv_mr_genes_coding_prefilt5_nondedup_min5kb_genes_PvTPM_100resamplings.txt", header=FALSE)
pv_ma_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/pv_ma_genes_coding_prefilt5_nondedup_min5kb_genes_PvTPM_100resamplings.txt", header=FALSE)

sst_mr_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/sst_mr_genes_coding_prefilt5_nondedup_min5kb_genes_SstTPM_100resamplings.txt", header=FALSE)
sst_ma_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/sst_ma_genes_coding_prefilt5_nondedup_min5kb_genes_SstTPM_100resamplings.txt", header=FALSE)

L4_mr_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L4_mr_genes_coding_prefilt5_nondedup_min5kb_genes_L4TPM_100resamplings.txt", header=FALSE)
L4_ma_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L4_ma_genes_coding_prefilt5_nondedup_min5kb_genes_L4TPM_100resamplings.txt", header=FALSE)

L5_mr_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L5_mr_genes_coding_prefilt5_nondedup_min5kb_genes_L5TPM_100resamplings.txt", header=FALSE)
L5_ma_5kb_resamps <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/5kb_genes/L5_ma_genes_coding_prefilt5_nondedup_min5kb_genes_L5TPM_100resamplings.txt", header=FALSE)

####subclass Mecp2-reg genes over unchanged resampled genes
#mCA
all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
)

all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes, MR resampled", "Cell-type unchanged genes, MA resampled", "Cell-type MR genes", "Cell-type MA genes")))
all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary)

#element mCA/CA fold difference of cell-type MR over unchanged genes
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"),1)
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"))

pval_fc_heatmap(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MR)/(MR-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchangedMRresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchangedMRresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#element mCA/CA fold difference of cell-type MA over unchanged genes
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"),1)
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"))

pval_fc_heatmap(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals, cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MA)/(MA-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchangedMAresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchangedMAresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')


##
#mCG
all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp = rbind(
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[1]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[1]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
  data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_5kb_resamps[[1]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
  data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[1]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
  data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
)

all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp %>% 
  group_by(element, cell_type, geneset, ) %>%
  summarise(median.val = median(val, na.rm = TRUE),
            sd.val = sd(val, na.rm = TRUE),
            n.val = n()) %>%
  mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))

all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes, MR resampled", "Cell-type unchanged genes, MA resampled", "Cell-type MR genes", "Cell-type MA genes")))
all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))

all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary)

#element mCG/CG fold difference of cell-type MR over unchanged genes
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"),1)
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"))

pval_fc_heatmap(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals, cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCG (Cell-type MR)/(MR-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchangedMRresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMR_over_unchangedMRresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#element mCG/CG fold difference of cell-type MA over unchanged genes
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc = round(cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"),1)
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"))

pval_fc_heatmap(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals, cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCG (Cell-type MA)/(MA-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchangedMAresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/CellTypeMA_over_unchangedMAresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')



####trying to use medians of resamplings
#mCA
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc_list <- list()
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals_list <- list()

cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc_list <- list()
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals_list <- list() 


for(i in 1:100){
  #mCA
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp = rbind(
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_PV_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_SST_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_L4_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCA[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_L5_mCA[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCA[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
  )
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp %>% 
    group_by(element, cell_type, geneset, ) %>%
    summarise(median.val = median(val, na.rm = TRUE),
              sd.val = sd(val, na.rm = TRUE),
              n.val = n()) %>%
    mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes, MR resampled", "Cell-type unchanged genes, MA resampled", "Cell-type MR genes", "Cell-type MA genes")))
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary)
  
  #element mCA/CA fold difference of cell-type MR over unchanged genes
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc_i = cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled")
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals_i = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"))
  
  
  #element mCA/CA fold difference of cell-type MA over unchanged genes
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc_i = cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled")
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals_i = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCA_resamp, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"))
  
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc_list[[i]] <- cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc_i
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals_list[[i]] <- cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals_i
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc_list[[i]] <- cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc_i
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals_list[[i]] <- cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals_i
}

# Convert the list of matrices to an array
array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc <- array(unlist(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals <- array(unlist(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc <- array(unlist(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals <- array(unlist(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals_list), dim = c(4, 6, 100))

# Compute the median along the third dimension (which corresponds to the different matrices)
median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc <- round(apply(array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc, c(1, 2), median),1)
median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals <- apply(array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals, c(1, 2), median)
median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc <- round(apply(array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc, c(1, 2), median),1)
median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals <- apply(array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals, c(1, 2), median)

rownames(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals) <- c("Pv", "Sst", "L4", "L5")


pval_fc_heatmap(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_pvals, median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MR)/(MR-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMR_over_unchangedMRresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMR_over_unchangedMRresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

pval_fc_heatmap(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_pvals, median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCA_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 80, "Ratio mCA (Cell-type MA)/(MA-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMA_over_unchangedMAresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMA_over_unchangedMAresamp_genes_mCAperCA_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

#mCG
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc_list <- list()
cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals_list <- list()

cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc_list <- list()
cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals_list <- list() 


for(i in 1:100){
  #mCG
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp = rbind(
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MR genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Pv", element="Region", geneset="Cell-type MA genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=PV_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Pv", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_PV_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Pv", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_PV_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% pv_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Pv", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MR genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="Sst", element="Region", geneset="Cell-type MA genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=SST_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="Sst", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_SST_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="Sst", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_SST_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% sst_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="Sst", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MR genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L4", element="Region", geneset="Cell-type MA genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=L4_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L4", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_L4_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L4", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L4_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L4_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L4", element="Extragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[i]]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MR genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(flank_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), flank_methylation_corrected], cell_type="L5", element="Region", geneset="Cell-type MA genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_5kb_resamps[[i]]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MR genes"),
    data.table(val=L5_INTACT_gene_body_TSSplus3kb_flank_mCG[is.finite(gene_methylation_corrected) & (gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), gene_methylation_corrected], cell_type="L5", element="Gene body", geneset="Cell-type MA genes"),
    data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_5kb_resamps[[i]]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MR genes"),
    data.table(val=promoterWindows_L5_mCG[is.finite(promoter_methylation_corrected) & (V4 %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >=5000,V4]), promoter_methylation_corrected], cell_type="L5", element="TSS", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="All cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Intragenic cCREs", geneset="Cell-type MA genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MR resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_5kb_resamps[[i]],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type unchanged genes, MA resampled"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_mr_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MR genes"),
    data.table(val=mousebrain_union_cCREs_L5_mCG[is.finite(cCRE_methylation_corrected) &(V4 %in% mousebrain_union_nonPromoter_cCREs_extragenicLinked_genes_genicBooleans[Gene %in% L5_ma_genes_q0.1_nondedup_mm9[V3-V2 >= 5000,V4],cCRE_label]), cCRE_methylation_corrected], cell_type="L5", element="Extragenic cCREs", geneset="Cell-type MA genes")
  )
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp %>% 
    group_by(element, cell_type, geneset, ) %>%
    summarise(median.val = median(val, na.rm = TRUE),
              sd.val = sd(val, na.rm = TRUE),
              n.val = n()) %>%
    mutate(error.val = qnorm(0.95)*sd.val / sqrt(n.val))
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(element = factor(element, levels=c("Region", "Gene body", "TSS", "All cCREs", "Intragenic cCREs", "Extragenic cCREs")))
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(geneset = factor(geneset, levels=c("Cell-type unchanged genes, MR resampled", "Cell-type unchanged genes, MA resampled", "Cell-type MR genes", "Cell-type MA genes")))
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary = all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary %>% mutate(cell_type = factor(cell_type, levels=c("Pv", "Sst", "L4", "L5")))
  
  all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt = data.table(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary)
  
  #element mCG/CG fold difference of cell-type MR over unchanged genes
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc_i = cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled")
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals_i = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp, "Cell-type MR genes", "Cell-type unchanged genes, MR resampled"))
  
  
  #element mCG/CG fold difference of cell-type MA over unchanged genes
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc_i = cellType_elem_fc_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp_summary_dt, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled")
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals_i = -log10(cellType_elem_pvals_cCREtypes_TSS(all_cells_MeCP2dys_cCREtypes_TSS_mCG_resamp, "Cell-type MA genes", "Cell-type unchanged genes, MA resampled"))
  
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc_list[[i]] <- cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc_i
  cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals_list[[i]] <- cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals_i
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc_list[[i]] <- cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc_i
  cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals_list[[i]] <- cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals_i
}

# Convert the list of matrices to an array
array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc <- array(unlist(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals <- array(unlist(cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc <- array(unlist(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc_list), dim = c(4, 6, 100))
array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals <- array(unlist(cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals_list), dim = c(4, 6, 100))

# Compute the median along the third dimension (which corresponds to the different matrices)
median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc <- round(apply(array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc, c(1, 2), median),1)
median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals <- apply(array_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals, c(1, 2), median)
median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc <- round(apply(array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc, c(1, 2), median),1)
median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals <- apply(array_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals, c(1, 2), median)

rownames(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc) <- c("Pv", "Sst", "L4", "L5")
rownames(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals) <- c("Pv", "Sst", "L4", "L5")


pval_fc_heatmap(median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_pvals, median_cellType_elem_cellTypeMRvsresampMR_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 20, "Ratio mCG (Cell-type MR)/(MR-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMR_over_unchangedMRresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMR_over_unchangedMRresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

pval_fc_heatmap(median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_pvals, median_cellType_elem_cellTypeMAvsresampMA_cCREtypes_TSS_mCG_fc, c("Region", "Gene body", "TSS", "All cCREs", "Intragenic CREs", "Extragenic cCREs"), "goldenrod2", 20, "Ratio mCG (Cell-type MA)/(MA-resampled)", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMA_over_unchangedMAresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.png", width = 5, height = 5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/median_CellTypeMA_over_unchangedMAresamp_genes_mCGperCG_region_geneBodyTSSplus3kb_TSS_allLinked_intragenicLinked_extragenicLinked_cCREs_wilcoxPvals_heatmap.eps", width = 5, height = 5.5, dpi = 300, units = "in", device='eps')

####Cacna1i diamond plots using replicates of PV and L5

#gene body and flank mCA
PV_WT_LIB041642_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_WT_LIB041642_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
PV_WT_LIB041644_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_WT_LIB041644_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
PV_KO_LIB041641_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_KO_LIB041641_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
PV_KO_LIB041643_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_KO_LIB041643_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")

L5_WT_LIB041645_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_WT_LIB041645_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L5_WT_LIB041648_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_WT_LIB041648_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L5_KO_LIB041646_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_KO_LIB041646_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L5_KO_LIB041647_gene_and_flank_mCA=fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_KO_LIB041647_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")


#cCRE methylation calculation
cCRE_meth_calc_func_ID <- function(meth_table, nonconv_table, label_column){
  meth_table$cCRE_methylation <- as.integer(meth_table[[5]])/as.integer(meth_table[[6]])
  meth_table$cCRE_methylation_corrected <- meth_table$cCRE_methylation - nonconv_table[GTAC_ESP_ID==label_column, nonconversion_rate]
  names(meth_table) = c("chrom", "start", "end", "cCRE_label", "meth_reads", "cyto_reads", "cCRE_methylation", "cCRE_methylation_corrected")
  meth_table[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
  return(meth_table)
}

#cCREs
mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/PV_reps/mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA_mm9.bed")

mousebrain_union_cCRE_L5_WT_LIB041645_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/L5_reps/mousebrain_union_cCRE_L5_WT_LIB041645_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_L5_WT_LIB041648_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/L5_reps/mousebrain_union_cCRE_L5_WT_LIB041648_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_L5_KO_LIB041646_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/L5_reps/mousebrain_union_cCRE_L5_KO_LIB041646_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCRE_L5_KO_LIB041647_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/L5_reps/mousebrain_union_cCRE_L5_KO_LIB041647_deep_INTACT_mCA_mm9.bed")

#nonconversion rates
lambda_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/intact_lambda_nonconversion_table.csv")
#calculating cCRE methylation with background nonconversion rate subtraction, column labeling
mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041642")
mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041644")
mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041641")
mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041643")

mousebrain_union_cCRE_L5_WT_LIB041645_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_L5_WT_LIB041645_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041645")
mousebrain_union_cCRE_L5_WT_LIB041648_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_L5_WT_LIB041648_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041648")
mousebrain_union_cCRE_L5_KO_LIB041646_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_L5_KO_LIB041646_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041646")
mousebrain_union_cCRE_L5_KO_LIB041647_deep_INTACT_mCA = cCRE_meth_calc_func_ID(meth_table=mousebrain_union_cCRE_L5_KO_LIB041647_deep_INTACT_mCA, nonconv_table=lambda_nonconv, label_column="LIB041647")

##
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE = data.table(rbind(
  cbind(val=L5_WT_LIB041645_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041645"),
  cbind(val=L5_WT_LIB041648_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041648"),
  cbind(val=L5_KO_LIB041646_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041646"),
  cbind(val=L5_KO_LIB041647_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041647"),
  cbind(val=L5_WT_LIB041645_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041645"),
  cbind(val=L5_WT_LIB041648_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041648"),
  cbind(val=L5_KO_LIB041646_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041646"),
  cbind(val=L5_KO_LIB041647_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041647"),
  cbind(val=median(mousebrain_union_cCRE_L5_WT_LIB041645_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041645"),
  cbind(val=median(mousebrain_union_cCRE_L5_WT_LIB041648_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L5", genotype="WT", GTAC_ESP_ID="LIB041648"),
  cbind(val=median(mousebrain_union_cCRE_L5_KO_LIB041646_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041646"),
  cbind(val=median(mousebrain_union_cCRE_L5_KO_LIB041647_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="L5", genotype="KO", GTAC_ESP_ID="LIB041647"),
  cbind(val=PV_WT_LIB041642_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041642"),
  cbind(val=PV_WT_LIB041644_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041644"),
  cbind(val=PV_KO_LIB041641_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041641"),
  cbind(val=PV_KO_LIB041643_gene_and_flank_mCA[gene=="Cacna1i", flank_methylation_corrected], element="Region", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041643"),
  cbind(val=PV_WT_LIB041642_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041642"),
  cbind(val=PV_WT_LIB041644_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041644"),
  cbind(val=PV_KO_LIB041641_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041641"),
  cbind(val=PV_KO_LIB041643_gene_and_flank_mCA[gene=="Cacna1i", gene_methylation_corrected], element="Gene body", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041643"),
  cbind(val=median(mousebrain_union_cCRE_PV_WT_LIB041642_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041642"),
  cbind(val=median(mousebrain_union_cCRE_PV_WT_LIB041644_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="PV", genotype="WT", GTAC_ESP_ID="LIB041644"),
  cbind(val=median(mousebrain_union_cCRE_PV_KO_LIB041641_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041641"),
  cbind(val=median(mousebrain_union_cCRE_PV_KO_LIB041643_deep_INTACT_mCA[(cCRE_label %in% mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_genicBooleans[Gene=="Cacna1i",cCRE_label]), cCRE_methylation_corrected], na.rm=TRUE), element="Intragenic cCREs", cell_type="PV", genotype="KO", GTAC_ESP_ID="LIB041643")
))


L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE = L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(element = factor(element, levels=c("Region", "Gene body", "Intragenic cCREs")))
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE = L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(cell_type = factor(cell_type, levels=c("L5", "PV")))
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE = L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))



L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE$val <- as.numeric(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE$val)


ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE[cell_type=="L5"] ,aes(x = element, y = as.numeric(val))) +
  #geom_point(color="orchid") +
  #geom_jitter(color='orchid', size=2.5, width = 0.1)+
  stat_summary(fun = mean, 
               geom = "point", size=4, color="orchid") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="orchid", size=0.5, width=0.25)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_reps_mean_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_reps_mean_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')



ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE[cell_type=="PV"] ,aes(x = element, y = as.numeric(val))) +
  #geom_point(color="orchid") +
  #geom_jitter(color='forestgreen', size=2.5, width = 0.1)+
  stat_summary(fun = mean, 
               geom = "point", size=4, color="forestgreen") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="forestgreen", size=0.5, width=0.25)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("PV INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_reps_mean_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_reps_mean_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')




###


#mean of WT and KO for each element
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_WT_KO_mean <- group_by(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE, element, cell_type, genotype) %>%
  summarise(
    mean.val = mean(val),
    sd.val = sd(val),
    n.val = n(),
    se.val = sd(val)/ sqrt(n.val)
  ) %>% data.table




ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_WT_KO_mean[cell_type == "L5"], 
       aes(x = element, y = as.numeric(mean.val), color = genotype)) +
  #geom_errorbar(aes(ymin = mean.val - se.val, ymax = mean.val + se.val), 
   #             width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "Genotype", values = c("WT" = "purple", "KO" = "orange")) +
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_separate_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_separate_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')

ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_WT_KO_mean[cell_type == "PV"], 
       aes(x = element, y = as.numeric(mean.val), color = genotype)) +
  #geom_errorbar(aes(ymin = mean.val - se.val, ymax = mean.val + se.val), 
   #             width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "Genotype", values = c("WT" = "purple", "KO" = "orange")) +
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("PV INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_separate_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_separate_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')


#pairs PV LIB041641-LIB041642, LIB041643-LIB041644
#pairs L5 LIB041645-LIB041646, LIB041647-LIB041648
##
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE$pair_group <- "NA"
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE[GTAC_ESP_ID %in% c("LIB041641", "LIB041642", "LIB041645", "LIB041646") , pair_group := "Pair group 1"]
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE[GTAC_ESP_ID %in% c("LIB041643", "LIB041644", "LIB041647", "LIB041648") , pair_group := "Pair group 2"]

L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE <- L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE %>% mutate(pair_group = factor(pair_group, levels=c("Pair group 1", "Pair group 2")))
##
#calculate mean and standard error of the mean of regional, gene body, and median intragenic cCRE methylation of Cacna1i for each pair
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_pair_summ <- group_by(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE, element, cell_type, pair_group) %>%
  summarise(
    mean.val = mean(val),
    sd.val = sd(val),
    n.val = n(),
    se.val = sd(val)/ sqrt(n.val)
  ) %>% data.table


ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_pair_summ[cell_type == "L5"], 
       aes(x = element, y = as.numeric(mean.val), color=pair_group)) +
  #geom_errorbar(aes(ymin = mean.val - se.val, ymax = mean.val + se.val), 
  #             width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.6), color="orchid") +
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_rep_pairs_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_rep_pairs_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')




ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_pair_summ[cell_type == "L5"], 
       aes(x = element, y = as.numeric(mean.val), color = pair_group)) +
  #geom_errorbar(aes(ymin = mean.val - se.val, ymax = mean.val + se.val), 
  #             width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "", values = c("Pair group 1" = "orchid", "Pair group 2" = "orchid")) +
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_separate_mCAperCA_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_separate_mCAperCA_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')


ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_pair_summ[cell_type == "PV"], 
       aes(x = element, y = as.numeric(mean.val), color = pair_group)) +
  #geom_errorbar(aes(ymin = mean.val - se.val, ymax = mean.val + se.val), 
  #             width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(name = "", values = c("Pair group 1" = "forestgreen", "Pair group 2" = "forestgreen")) +
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("PV INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))






##


#calculate mean and standard error of the mean of regional, gene body, and median intragenic cCRE methylation of Cacna1i across reps
L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_summ <- group_by(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE, element, cell_type, ) %>%
  summarise(
    mean.val = mean(val),
    sd.val = sd(val),
    n.val = n(),
    se.val = sd(val)/ sqrt(n.val)
  ) %>% data.table




ggplot(L5_PV_reps_Cacna1i_region_geneBody_intragenicLinkedcCRE_summ[cell_type=="L5"] ,aes(x = element, y = as.numeric(val))) +
  geom_point(shape = 16, size = 10) +  # Diamond shape
  geom_errorbar(aes(ymin = long_highmCA_gene_logfc_mean - long_highmCA_gene_logfc_se, ymax = long_highmCA_gene_logfc_mean + long_highmCA_gene_logfc_se), 
                width = 0.0025) +  # Error bars for standard error
  geom_text(stat='identity', aes(label=round(as.numeric(val), 3)), vjust=-0.5)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("L5 INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_mCAperCA_diamondplot.png', width = 3.7, height = 5, dpi = 300, units = "in", device='png')
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_L5_INTACT_WT_KO_mCAperCA_diamondplot.eps', width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_cells_Cacna1i_region_geneBody_intragenicLinkedcCRE[cell_type=="Pv"] ,aes(x = element, y = as.numeric(val))) +
  geom_point(shape = 18, size = 10) +  # Diamond shape
  geom_text(stat='identity', aes(label=round(as.numeric(val), 3)), vjust=-0.5)+
  coord_cartesian(ylim = c(0, 0.07)) +
  ylab("PV INTACT mCA/CA") +
  theme_bw() +
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x=element_text(size=15, angle=90))
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_mCAperCA_diamondplot.png', width = 3.7, height = 5, dpi = 300, units = "in", device='png')
#ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/Cacna1i_region_geneBody_intragenicLinkedcCRE_PV_INTACT_WT_KO_mCAperCA_diamondplot.eps', width = 3.7, height = 5, dpi = 300, units = "in", device='eps')

###
#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")
SST_deseq2 <- data.table(SST_deseq, keep.rownames="Gene")
L4_deseq2 <- data.table(L4_deseq, keep.rownames="Gene")
L5_deseq2 <- data.table(L5_deseq, keep.rownames="Gene")

#
PV_L5_deseq2 = rbind(
  cbind(PV_deseq2, subclass="PV"),
  cbind(L5_deseq2, subclass="L5")
) %>% data.table

PV_L5_deseq2 <- PV_L5_deseq2 %>% mutate(subclass = factor(subclass, levels=c("PV", "L5")))

ggplot(PV_L5_deseq2[Gene == "Cacna1i"], 
       aes(x = subclass, y = as.numeric(ashr_log2FoldChange), color=subclass)) +
  geom_errorbar(aes(ymin = ashr_log2FoldChange - ashr_lfcSE, ymax = ashr_log2FoldChange + ashr_lfcSE), 
               width = 0.25) +
  geom_point(size = 4) +
  coord_cartesian(ylim = c(-0.05, 0.3)) +
  geom_hline(yintercept = 0, color = "black") +
  scale_color_manual(name = "", values = c("PV" = "forestgreen", "L5" = "orchid")) +
  ylab("Cacna1i log2 fold change (KO/WT)") +
  theme_bw() +
  theme(legend.position = "None", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(color = "black"), 
        axis.ticks.x = element_blank(), axis.ticks.y=element_line(color="black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90))
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_PV_Cacna1i_ashr_log2FoldChange_lfcSE_dotplot.png', width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave('HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_PV_Cacna1i_ashr_log2FoldChange_lfcSE_dotplot.eps', width = 3, height = 5, dpi = 300, units = "in", device='eps')
