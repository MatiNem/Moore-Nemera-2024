library(data.table)
library(dplyr)
library(ggplot2)
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
chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

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

PvMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
PvUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/PvUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
SstMR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/SstMR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
SstUnchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/SstUnchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L4MR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4MR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4Unchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4Unchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L5MR_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5MR_metaMR_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5Unchanged_metaMRgenes = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5Unchanged_metaMR_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")

meta_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
meta_MA_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MA_genes_geneColumn4_mm9.bed")
meta_unchanged_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/MeCP2_unchanged_meta_genes_mm9.bed")

nonmeta_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/nonMetaMR_genes_geneColumn4_mm9.bed")

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

#using only chromosomes 1 to 19 and X
PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX =  PV_INTACT_gene_body_TSSplus3kb_flank_mCA[chrom %in% chrom_list, ]
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX =  SST_INTACT_gene_body_TSSplus3kb_flank_mCA[chrom %in% chrom_list, ]
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX =  L4_INTACT_gene_body_TSSplus3kb_flank_mCA[chrom %in% chrom_list, ]
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX =  L5_INTACT_gene_body_TSSplus3kb_flank_mCA[chrom %in% chrom_list, ]

PV_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX =  PV_INTACT_gene_body_TSSplus3kb_flank_mCG[chrom %in% chrom_list, ]
SST_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX =  SST_INTACT_gene_body_TSSplus3kb_flank_mCG[chrom %in% chrom_list, ]
L4_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX =  L4_INTACT_gene_body_TSSplus3kb_flank_mCG[chrom %in% chrom_list, ]
L5_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX =  L5_INTACT_gene_body_TSSplus3kb_flank_mCG[chrom %in% chrom_list, ]

#deciles of flank methylation
PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_meth_decile <- as.character(ntile(PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_methylation_corrected, 10))
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_meth_decile <- as.character(ntile(SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_methylation_corrected, 10))
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_meth_decile <- as.character(ntile(L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_methylation_corrected, 10))
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_meth_decile <- as.character(ntile(L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX$flank_methylation_corrected, 10))

PV_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_meth_decile <- as.character(ntile(PV_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_methylation_corrected, 10))
SST_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_meth_decile <- as.character(ntile(SST_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_methylation_corrected, 10))
L4_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_meth_decile <- as.character(ntile(L4_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_methylation_corrected, 10))
L5_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_meth_decile <- as.character(ntile(L5_INTACT_gene_body_TSSplus3kb_flank_mCG_chr1toX$flank_methylation_corrected, 10))

PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX = PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX = SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX = L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX = L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))




meth_decile_oddsRatio_pval = function(flank_meth_data_table, cellTypeMR_genes, cellTypeUnchanged_genes){
  p_values_cellTypeMR_metaMR = rep(0,10)
  log2fc_cellTypeMR_metaMR = rep(0,10)
  for(i in 1:10){
  cellTypeMR_metaMR = length(intersect(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% cellTypeMR_genes[,V4]), gene], meta_MR_genes[,V4]))
  cellTypeMR_nonmetaMR = length(intersect(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% cellTypeMR_genes[,V4]), gene], nonmeta_MR_genes[,V4]))
  cellTypeUnchanged_metaMR = length(intersect(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% cellTypeUnchanged_genes[,V4]), gene], meta_MR_genes[,V4]))
  cellTypeUnchanged_nonmetaMR = length(intersect(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% cellTypeUnchanged_genes[,V4]), gene], nonmeta_MR_genes[,V4]))
  dat_cellTypeMR_metaMR = data.frame(c(cellTypeMR_metaMR, cellTypeMR_nonmetaMR),
                   c(cellTypeUnchanged_metaMR, cellTypeUnchanged_nonmetaMR))
  test_cellTypeMR_metaMR = fisher.test(dat_cellTypeMR_metaMR)
  p_values_cellTypeMR_metaMR[i] =  test_cellTypeMR_metaMR$p.value
  log2fc_cellTypeMR_metaMR[i] = log2(test_cellTypeMR_metaMR$estimate)
  }
  pval_fc_table = data.table(cbind(p_values=p_values_cellTypeMR_metaMR, log2fc=log2fc_cellTypeMR_metaMR))
  return(pval_fc_table)
}


decile_pval_fc_heatmap <- function(pval_matrix, fc_matrix, column_names, color, pval_max, title, x_size=18, y_size=18, x_angle=0, y_angle=0, tile_text_size=6){
  test <- as.data.frame(pval_matrix)
  test_fc <- as.data.frame(fc_matrix)
  colnames(test) = column_names
  colnames(test_fc) = column_names
  test$decile = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  test_fc$decile = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  #Melting data so we can plot it with GGplot
  test.m <- suppressWarnings(melt(test,id.vars = c("decile")))
  test_fc.m <- suppressWarnings(melt(test_fc,id.vars = c("decile")))
  #Resetting factors
  test.m = test.m %>% mutate(decile = factor(decile, levels=c("10", "9", "8", "7", "6", "5", "4", "3", "2", "1")))
  test_fc.m = test_fc.m %>% mutate(decile = factor(decile, levels=c("10", "9", "8", "7", "6", "5", "4", "3", "2", "1")))
  test_pval_and_fc = data.table(cbind(test.m, fc = test_fc.m$value))
  test_pval_and_fc$pvalMax = test_pval_and_fc$value
  test_pval_and_fc[value > pval_max, pvalMax := pval_max]

  #Creating the plot itself
  ggplot(test_pval_and_fc,aes(variable,decile)) + geom_tile(aes(fill=pvalMax),color = "white") +
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

decile_pval_fc_heatmap(cellType_metaMR_flank_mCA_decile_log10pvals, cellType_metaMR_flank_mCA_decile_log2fc, c("L5", "L4", "PV", "SST"), "red4", 20, "", x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/cellTypeMR_metaMR_overlap_in_flank_mL5_mCAperCA_deciles_wilcoxPvals_heatmap.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/cellTypeMR_metaMR_overlap_in_flank_mL5_mCAperCA_deciles_wilcoxPvals_heatmap.eps", width = 4, height = 6, dpi = 300, units = "in", device='eps')

PV_decileCounts <- PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX[is.finite(flank_meth_decile) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]),][, .N, by = flank_meth_decile] %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
SST_decileCounts <- SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX[is.finite(flank_meth_decile) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4]),][, .N, by = flank_meth_decile] %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L4_decileCounts <- L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX[is.finite(flank_meth_decile) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4]),][, .N, by = flank_meth_decile] %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L5_decileCounts <- L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX[is.finite(flank_meth_decile) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4]),][, .N, by = flank_meth_decile] %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))

#cell type MR genes in same cell-type regional mCA/CA deciles
ggplot(data=PV_decileCounts, aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 10, stat = "identity") +  # Diamond shape
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank PV mCA/CA decile") + ylab("Number of PV MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_PVMR_genes_in_PV_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_PVMR_genes_in_PV_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=SST_decileCounts, aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 10, stat = "identity") +  # Diamond shape
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank SST mCA/CA decile") + ylab("Number of SST MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_SSTMR_genes_in_SST_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_SSTMR_genes_in_SST_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=L4_decileCounts, aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 10, stat = "identity") +  # Diamond shape
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank L4 mCA/CA decile") + ylab("Number of L4 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L4MR_genes_in_L4_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L4MR_genes_in_L4_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=L5_decileCounts, aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 10, stat = "identity") +  # Diamond shape
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank L5 mCA/CA decile") + ylab("Number of L5 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L5MR_genes_in_L5_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L5MR_genes_in_L5_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')



pv_mr_genes_resamp1 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
sst_mr_genes_resamp1 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Sst_mr_genes_coding_prefilt5_nondedup_SstTPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
L4_mr_genes_resamp1 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L4_mr_genes_coding_prefilt5_nondedup_L4TPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1
L5_mr_genes_resamp1 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L5_mr_genes_coding_prefilt5_nondedup_L5TPM_MeCP2KO_nondedup/resamp1.txt", header=FALSE)$V1


#function for enrichment of MR genes being in deciles over expression-resampled unchanged genes
meth_decile_pval_resamp = function(flank_meth_data_table, MR_genes, resamp_genes){
  p_values_cellTypeMR = rep(0,10)
  for(i in 1:10){
    decile_MR = length(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% MR_genes), gene])
    nondecile_MR = length(flank_meth_data_table[(flank_meth_decile!=i) & (gene %in% MR_genes), gene])
    decile_resamp = length(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% resamp_genes), gene])
    nondecile_resamp = length(flank_meth_data_table[(flank_meth_decile!=i) & (gene %in% resamp_genes), gene])
    dat_cellTypeMR = data.frame(c(decile_MR, decile_resamp),
                                       c(nondecile_MR, nondecile_resamp))
    test_cellTypeMR = fisher.test(dat_cellTypeMR)
    p_values_cellTypeMR[i] = test_cellTypeMR$p.value
  }
  return(p_values_cellTypeMR) 
}

meth_decile_oddsRatio_resamp = function(flank_meth_data_table, MR_genes, resamp_genes){
  p_values_cellTypeMR = rep(0,10)
  log2fc_cellTypeMR = rep(0,10)
  for(i in 1:10){
    decile_MR = length(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% MR_genes), gene])
    nondecile_MR = length(flank_meth_data_table[(flank_meth_decile!=i) & (gene %in% MR_genes), gene])
    decile_resamp = length(flank_meth_data_table[(flank_meth_decile==i) & (gene %in% resamp_genes), gene])
    nondecile_resamp = length(flank_meth_data_table[(flank_meth_decile!=i) & (gene %in% resamp_genes), gene])
    dat_cellTypeMR = data.frame(c(decile_MR, decile_resamp),
                                c(nondecile_MR, nondecile_resamp))
    test_cellTypeMR = fisher.test(dat_cellTypeMR)
    #p_values_cellTypeMR[i] = test_cellTypeMR$p.value
    log2fc_cellTypeMR[i] = log2(test_cellTypeMR$estimate)
  }
  #pval_fc_table = data.table(cbind(p_values=p_values_cellTypeMR, log2fc=log2fc_cellTypeMR))
  #return(pval_fc_table)
  return(log2fc_cellTypeMR) 
}


pv_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
sst_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
L4_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
L5_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))

for (i in 1:1000){
  #Pv
  pv_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  pv_df_new = meth_decile_pval_resamp(mPv_snmcseq_gene_and_flank_mCA_chr1toX, pv_mr_genes_q0.1_nondedup_mm9[,V4],  pv_mr_genes_resamp)
  pv_df = rbind(pv_df, as.list(pv_df_new))
  #Sst
  sst_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Sst_mr_genes_coding_prefilt5_nondedup_SstTPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  sst_df_new = meth_decile_pval_resamp(mSst_snmcseq_gene_and_flank_mCA_chr1toX, sst_mr_genes_q0.1_nondedup_mm9[,V4],  sst_mr_genes_resamp)
  sst_df = rbind(sst_df, as.list(sst_df_new))
  #L4
  L4_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L4_mr_genes_coding_prefilt5_nondedup_L4TPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  L4_df_new = meth_decile_pval_resamp(mL4_snmcseq_gene_and_flank_mCA_chr1toX, L4_mr_genes_q0.1_nondedup_mm9[,V4],  L4_mr_genes_resamp)
  L4_df = rbind(L4_df, as.list(L4_df_new))
  #L5
  L5_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L5_mr_genes_coding_prefilt5_nondedup_L5TPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  L5_df_new = meth_decile_pval_resamp(mL5_snmcseq_gene_and_flank_mCA_chr1toX, L5_mr_genes_q0.1_nondedup_mm9[,V4],  L5_mr_genes_resamp)
  L5_df = rbind(L5_df, as.list(L5_df_new))
}

pv_df_medians = sapply(pv_df, median)
sst_df_medians = sapply(sst_df, median)
L4_df_medians = sapply(L4_df, median)
L5_df_medians = sapply(L5_df, median)

write.table(pv_df_medians, file="HG_lab/Mati/GabelLab/gene_flank_plots/pv_mr_genes_in_mPv_flank_mCAperCA_deciles_median_pvals.txt", quote=F, row.names=F, sep="\t")
write.table(sst_df_medians, file="HG_lab/Mati/GabelLab/gene_flank_plots/sst_mr_genes_in_mSst_flank_mCAperCA_deciles_median_pvals.txt", quote=F, row.names=F, sep="\t")
write.table(L4_df_medians, file="HG_lab/Mati/GabelLab/gene_flank_plots/L4_mr_genes_in_mL4_flank_mCAperCA_deciles_median_pvals.txt", quote=F, row.names=F, sep="\t")
write.table(L5_df_medians, file="HG_lab/Mati/GabelLab/gene_flank_plots/L5_mr_genes_in_mL5_flank_mCAperCA_deciles_median_pvals.txt", quote=F, row.names=F, sep="\t")

pv_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), pvals = pv_df_medians, sigs = sapply(pv_df_medians, sig_function)))
sst_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), pvals = sst_df_medians, sigs = sapply(sst_df_medians, sig_function)))
L4_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), pvals = L4_df_medians, sigs = sapply(L4_df_medians, sig_function)))
L5_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), pvals = L5_df_medians, sigs = sapply(L5_df_medians, sig_function)))

mPv_snmcseq_gene_and_flank_mCA_chr1toX = data.table(left_join(mPv_snmcseq_gene_and_flank_mCA_chr1toX, pv_df_medians_flankLabels, by=c("flank_meth_decile")))
mSst_snmcseq_gene_and_flank_mCA_chr1toX = data.table(left_join(mSst_snmcseq_gene_and_flank_mCA_chr1toX, sst_df_medians_flankLabels, by=c("flank_meth_decile")))
mL4_snmcseq_gene_and_flank_mCA_chr1toX = data.table(left_join(mL4_snmcseq_gene_and_flank_mCA_chr1toX, L4_df_medians_flankLabels, by=c("flank_meth_decile")))
mL5_snmcseq_gene_and_flank_mCA_chr1toX = data.table(left_join(mL5_snmcseq_gene_and_flank_mCA_chr1toX, L5_df_medians_flankLabels, by=c("flank_meth_decile")))

mPv_snmcseq_gene_and_flank_mCA_chr1toX = mPv_snmcseq_gene_and_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mSst_snmcseq_gene_and_flank_mCA_chr1toX = mSst_snmcseq_gene_and_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mL4_snmcseq_gene_and_flank_mCA_chr1toX = mL4_snmcseq_gene_and_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mL5_snmcseq_gene_and_flank_mCA_chr1toX = mL5_snmcseq_gene_and_flank_mCA_chr1toX %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))

decile_list = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
ggplot(data=mPv_snmcseq_gene_and_flank_mCA_chr1toX[(flank_meth_decile %in% decile_list) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=sigs), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mPv mCA/CA decile") + ylab("Number of Pv MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_PvMR_genes_in_mPv_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_PvMR_genes_in_mPv_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=mSst_snmcseq_gene_and_flank_mCA_chr1toX[(flank_meth_decile %in% decile_list) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=sigs), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mSst mCA/CA decile") + ylab("Number of Sst MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_SstMR_genes_in_mSst_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_SstMR_genes_in_mSst_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=mL4_snmcseq_gene_and_flank_mCA_chr1toX[(flank_meth_decile %in% decile_list) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=sigs), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL4 mCA/CA decile") + ylab("Number of L4 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L4MR_genes_in_mL4_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L4MR_genes_in_mL4_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=mL5_snmcseq_gene_and_flank_mCA_chr1toX[(flank_meth_decile %in% decile_list) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=sigs), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL5 mCA/CA decile") + ylab("Number of L5 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L5MR_genes_in_mL5_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L5MR_genes_in_mL5_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


pv_fc_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
sst_fc_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
L4_fc_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
L5_fc_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))

for (i in 1:1000){
  #Pv
  pv_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Pv_mr_genes_coding_prefilt5_nondedup_PvTPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  pv_fc_df_new = meth_decile_oddsRatio_resamp(PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, pv_mr_genes_q0.1_nondedup_mm9[,V4],  pv_mr_genes_resamp)
  pv_fc_df = rbind(pv_fc_df, as.list(pv_fc_df_new))
  #Sst
  sst_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/Sst_mr_genes_coding_prefilt5_nondedup_SstTPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  sst_fc_df_new = meth_decile_oddsRatio_resamp(SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, sst_mr_genes_q0.1_nondedup_mm9[,V4],  sst_mr_genes_resamp)
  sst_fc_df = rbind(sst_fc_df, as.list(sst_fc_df_new))
  #L4
  L4_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L4_mr_genes_coding_prefilt5_nondedup_L4TPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  L4_fc_df_new = meth_decile_oddsRatio_resamp(L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, L4_mr_genes_q0.1_nondedup_mm9[,V4],  L4_mr_genes_resamp)
  L4_fc_df = rbind(L4_fc_df, as.list(L4_fc_df_new))
  #L5
  L5_mr_genes_resamp = fread(paste0("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/L5_mr_genes_coding_prefilt5_nondedup_L5TPM_MeCP2KO_nondedup/resamp",i,".txt"), header=FALSE)$V1
  L5_fc_df_new = meth_decile_oddsRatio_resamp(L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, L5_mr_genes_q0.1_nondedup_mm9[,V4],  L5_mr_genes_resamp)
  L5_fc_df = rbind(L5_fc_df, as.list(L5_fc_df_new))
}

pv_fc_df_medians = sapply(pv_fc_df, median)
sst_fc_df_medians = sapply(sst_fc_df, median)
L4_fc_df_medians = sapply(L4_fc_df, median)
L5_fc_df_medians = sapply(L5_fc_df, median)

write.table(pv_fc_df_medians, file="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/pv_mr_genes_in_PV_flank_INTACT_WT_KO_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(sst_fc_df_medians, file="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/sst_mr_genes_in_SST_flank_INTACT_WT_KO_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(L4_fc_df_medians, file="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L4_mr_genes_in_L4_flank_INTACT_WT_KO_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(L5_fc_df_medians, file="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_mr_genes_in_L5_flank_INTACT_WT_KO_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")

pv_fc_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = pv_fc_df_medians))
sst_fc_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = sst_fc_df_medians))
L4_fc_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = L4_fc_df_medians))
L5_fc_df_medians_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = L5_fc_df_medians))

PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = data.table(left_join(PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, pv_fc_df_medians_flankLabels, by=c("flank_meth_decile")))
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = data.table(left_join(SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, sst_fc_df_medians_flankLabels, by=c("flank_meth_decile")))
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = data.table(left_join(L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, L4_fc_df_medians_flankLabels, by=c("flank_meth_decile")))
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = data.table(left_join(L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX, L5_fc_df_medians_flankLabels, by=c("flank_meth_decile")))

PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels = L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))


PV_decileCounts_withLabels <- PV_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels[is.finite(flank_meth_decile) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]),][, .(N = .N, oddsRatio = unique(oddsRatio)), by = flank_meth_decile]
SST_decileCounts_withLabels <- SST_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels[is.finite(flank_meth_decile) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4]),][, .(N = .N, oddsRatio = unique(oddsRatio)), by = flank_meth_decile]
L4_decileCounts_withLabels <- L4_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels[is.finite(flank_meth_decile) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4]),][, .(N = .N, oddsRatio = unique(oddsRatio)), by = flank_meth_decile]
L5_decileCounts_withLabels <- L5_INTACT_gene_body_TSSplus3kb_flank_mCA_chr1toX_withLabels[is.finite(flank_meth_decile) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4]),][, .(N = .N, oddsRatio = unique(oddsRatio)), by = flank_meth_decile]




#cell type MR genes in same cell-type regional mCA/CA deciles
ggplot(data=PV_decileCounts_withLabels , aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 8, stat = "identity") +  # Diamond shape
  geom_text(aes(label=round(as.numeric(oddsRatio), 1)), vjust=-2, size=2) +  # Text labels
  coord_cartesian(ylim = c(0, 50)) +
  xlab("Gene flank PV mCA/CA decile") + ylab("Number of PV MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_PVMR_genes_in_PV_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_PVMR_genes_in_PV_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=SST_decileCounts_withLabels , aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 8, stat = "identity") +  # Diamond shape
  geom_text(aes(label=round(as.numeric(oddsRatio), 1)), vjust=-2, size=2) +  # Text labels
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Gene flank SST mCA/CA decile") + ylab("Number of SST MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_SSTMR_genes_in_SST_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_SSTMR_genes_in_SST_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=L4_decileCounts_withLabels , aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 8, stat = "identity") +  # Diamond shape
  geom_text(aes(label=round(as.numeric(oddsRatio), 1)), vjust=-2, size=2) +  # Text labels
  coord_cartesian(ylim = c(0, 20)) +
  xlab("Gene flank L4 mCA/CA decile") + ylab("Number of L4 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L4MR_genes_in_L4_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L4MR_genes_in_L4_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=L5_decileCounts_withLabels , aes(x=flank_meth_decile, y=N)) +
  geom_point(shape = 18, size = 8, stat = "identity") +  # Diamond shape
  geom_text(aes(label=round(as.numeric(oddsRatio), 1)), vjust=-2, size=2) +  # Text labels
  coord_cartesian(ylim = c(0, 35)) +
  xlab("Gene flank L5 mCA/CA decile") + ylab("Number of L5 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L5MR_genes_in_L5_flank_INTACT_WT_KO_mCA_deciles_diamondplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/number_L5MR_genes_in_L5_flank_INTACT_WT_KO_mCA_deciles_diamondplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')





ggplot(data=mPv_snmcseq_gene_and_flank_mCA_chr1toX_withLabels[(flank_meth_decile %in% decile_list) & (gene %in% pv_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2, size=4)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mPv mCA/CA decile") + ylab("Number of Pv MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_PvMR_genes_in_mPv_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_PvMR_genes_in_mPv_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=mSst_snmcseq_gene_and_flank_mCA_chr1toX_withLabels[(flank_meth_decile %in% decile_list) & (gene %in% sst_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2, size=4)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mSst mCA/CA decile") + ylab("Number of Sst MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_SstMR_genes_in_mSst_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_SstMR_genes_in_mSst_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=mL4_snmcseq_gene_and_flank_mCA_chr1toX_withLabels[(flank_meth_decile %in% decile_list) & (gene %in% L4_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2, size=4)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL4 mCA/CA decile") + ylab("Number of L4 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L4MR_genes_in_mL4_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L4MR_genes_in_mL4_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=mL5_snmcseq_gene_and_flank_mCA_chr1toX_withLabels[(flank_meth_decile %in% decile_list) & (gene %in% L5_mr_genes_q0.1_nondedup_mm9[,V4]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2, size=4)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL5 mCA/CA decile") + ylab("Number of L5 MR genes")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L5MR_genes_in_mL5_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/number_L5MR_genes_in_mL5_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##
m1 = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_1_bed_na.bed")
m2 = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_2_bed_na.bed")
m3 = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_3_bed_na.bed")
g3 = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_bed_na.bed")

oddsRatio_median_func = function(number_resamp, resamp_file_directory, meth_file, genes_of_interest){
  fc_df = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"))
  for (i in 1:number_resamp){
    #Pv
    resamp = fread(paste0(resamp_file_directory, "resamp",i,".txt"), header=FALSE)$V1
    fc_df_new = meth_decile_oddsRatio_resamp(meth_file, genes_of_interest,  resamp)
    fc_df = rbind(fc_df, as.list(fc_df_new))
  }
  fc_df_median = sapply(fc_df, median)
}

pv_fc_df_medians_g3 = oddsRatio_median_func(1000, "HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/resampled_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_PvTPM_MeCP2KO_nondedup/", mPv_snmcseq_gene_and_flank_mCA_chr1toX, g3[,V6])
sst_fc_df_medians_g3 = oddsRatio_median_func(1000, "HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/resampled_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_SstTPM_MeCP2KO_nondedup/", mSst_snmcseq_gene_and_flank_mCA_chr1toX, g3[,V6])
L4_fc_df_medians_g3 = oddsRatio_median_func(1000, "HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/resampled_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_L4TPM_MeCP2KO_nondedup/", mL4_snmcseq_gene_and_flank_mCA_chr1toX, g3[,V6])
L5_fc_df_medians_g3 = oddsRatio_median_func(1000, "HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/resampled_genes/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_L5TPM_MeCP2KO_nondedup/", mL5_snmcseq_gene_and_flank_mCA_chr1toX, g3[,V6])
  

write.table(pv_fc_df_medians_g3, file="HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mPv_flank_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(sst_fc_df_medians_g3, file="HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mSst_flank_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(L4_fc_df_medians_g3, file="HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL4_flank_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")
write.table(L5_fc_df_medians_g3, file="HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL5_flank_mCAperCA_deciles_median_log2oddsRatio.txt", quote=F, row.names=F, sep="\t")

pv_fc_df_medians_g3_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = pv_fc_df_medians_g3))
sst_fc_df_medians_g3_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = sst_fc_df_medians_g3))
L4_fc_df_medians_g3_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = L4_fc_df_medians_g3))
L5_fc_df_medians_g3_flankLabels = data.table(cbind(flank_meth_decile=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), oddsRatio = L5_fc_df_medians_g3))

mPv_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = data.table(left_join(mPv_snmcseq_gene_and_flank_mCA_chr1toX, pv_fc_df_medians_g3_flankLabels, by=c("flank_meth_decile")))
mSst_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = data.table(left_join(mSst_snmcseq_gene_and_flank_mCA_chr1toX, sst_fc_df_medians_g3_flankLabels, by=c("flank_meth_decile")))
mL4_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = data.table(left_join(mL4_snmcseq_gene_and_flank_mCA_chr1toX, L4_fc_df_medians_g3_flankLabels, by=c("flank_meth_decile")))
mL5_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = data.table(left_join(mL5_snmcseq_gene_and_flank_mCA_chr1toX, L5_fc_df_medians_g3_flankLabels, by=c("flank_meth_decile")))

mPv_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = mPv_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mSst_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = mSst_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mL4_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = mL4_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))
mL5_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels = mL5_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels %>% mutate(flank_meth_decile = factor(flank_meth_decile, levels=c("1","2","3","4","5","6","7","8","9","10")))


ggplot(data=mPv_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels[(flank_meth_decile %in% decile_list) & (gene %in% g3[,V6]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mPv mCA/CA decile") + ylab("Number of genes used >3 times")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mPv_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mPv_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data=mSst_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels[(flank_meth_decile %in% decile_list) & (gene %in% g3[,V6]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mSst mCA/CA decile") + ylab("Number of genes used >3 times")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mSst_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mSst_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=mL4_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels[(flank_meth_decile %in% decile_list) & (gene %in% g3[,V6]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL4 mCA/CA decile") + ylab("Number of genes used >3 times")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL4_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL4_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=mL5_snmcseq_gene_and_flank_mCA_chr1toX_g3_flankLabels[(flank_meth_decile %in% decile_list) & (gene %in% g3[,V6]),], aes(x=flank_meth_decile)) +
  geom_bar()+
  geom_text(stat='count', aes(label=round(as.numeric(oddsRatio),1)), vjust=-0.2)+
  scale_x_discrete(drop=FALSE, na.translate = FALSE)+
  xlab("Gene flank mL5 mCA/CA decile") + ylab("Number of genes used >3 times")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.text=element_text(size=15), axis.line=element_line(color="black"), axis.title=element_text(size=15))
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL5_flank_mCA_deciles_barplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/gene_flank_plots/allen18_L5_sst_pv_cluster_n123_g_tabled_m_g3_in_mL5_flank_mCA_deciles_barplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')
