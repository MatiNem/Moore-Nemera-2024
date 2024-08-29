library(data.table)
library(dplyr)
library(edgeR)
library(vioplot)
#install.packages("vioplot")

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

cpm_function <- function(data_table, read_column){
  data_table[, cpm := (10**6) * data_table[, get(read_column)]/sum(data_table[, get(read_column)])]
}

nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")


#promoters of genes enriched in Exc over Pv cells or Pv over exc cells
Exc_over_Pv_genes_top200_promoterWindows = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top200_promoterWindows_mm9.bed")
Pv_over_Exc_genes_top200_promoterWindows = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top200_promoterWindows_mm9.bed")

Exc_over_Pv_genes_top200_mm9 = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top200_mm9.bed")
Pv_over_Exc_genes_top200_mm9 = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top200_mm9.bed")

Exc_over_Pv_genes_top100_q0.01_mm9 = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top100_q0.01_mm9.bed")
Pv_over_Exc_genes_top100_q0.01_mm9 = fread("HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top100_q0.01_mm9.bed")


#cCREs linked to these genes
Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9.bed")
Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")
Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed")

Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9.bed")
Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9.bed")


#non-promoter cCREs linked to genes
mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")


prom_H3K27ac_PV_WTKO_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WTko_B1_AAGGAGTC_chipCounts.bed")
prom_H3K27ac_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WTKO_B2_AAGGAGTC_chipCounts.bed")
prom_H3K27ac_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WTKO_B3_CTTCCCAT_chipCounts.bed")
prom_H3K27ac_PV_KO_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_KO_B1_GCTGAGTC_chipCounts.bed")
prom_H3K27ac_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_KO_B2_GCTGAGTC_chipCounts.bed")
prom_H3K27ac_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_KO_B3_TCGCGAGG_chipCounts.bed")


#counts per million manually
prom_H3K27ac_PV_WTKO_B1[, cpm := (10**6) * (prom_H3K27ac_PV_WTKO_B1[,V7]/sum(prom_H3K27ac_PV_WTKO_B1[,V7]))]
prom_H3K27ac_PV_WTKO_B2[, cpm := (10**6) * (prom_H3K27ac_PV_WTKO_B2[,V7]/sum(prom_H3K27ac_PV_WTKO_B2[,V7]))]
prom_H3K27ac_PV_WTKO_B3[, cpm := (10**6) * (prom_H3K27ac_PV_WTKO_B3[,V7]/sum(prom_H3K27ac_PV_WTKO_B3[,V7]))]
prom_H3K27ac_PV_KO_B1[, cpm := (10**6) * (prom_H3K27ac_PV_KO_B1[,V7]/sum(prom_H3K27ac_PV_KO_B1[,V7]))]
prom_H3K27ac_PV_KO_B2[, cpm := (10**6) * (prom_H3K27ac_PV_KO_B2[,V7]/sum(prom_H3K27ac_PV_KO_B2[,V7]))]
prom_H3K27ac_PV_KO_B3[, cpm := (10**6) * (prom_H3K27ac_PV_KO_B3[,V7]/sum(prom_H3K27ac_PV_KO_B3[,V7]))]


prom_Input_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_Input_PV_WTKO_B2_CGACCTAA_chipCounts.bed")
prom_Input_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_Input_PV_WTKO_B3_CTCGAACA_chipCounts.bed")
prom_Input_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_Input_PV_KO_B2_TACATCGG_chipCounts.bed")
prom_Input_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_KO/ensgene_mm9_1kb_promoterWindows_Input_PV_KO_B3_TCTAGGAG_chipCounts.bed")

prom_Input_PV_WTKO_B2[, cpm := (10**6) * (prom_Input_PV_WTKO_B2[,V7]/sum(prom_Input_PV_WTKO_B2[,V7]))]
prom_Input_PV_WTKO_B3[, cpm := (10**6) * (prom_Input_PV_WTKO_B3[,V7]/sum(prom_Input_PV_WTKO_B3[,V7]))]
prom_Input_PV_KO_B2[, cpm := (10**6) * (prom_Input_PV_KO_B2[,V7]/sum(prom_Input_PV_KO_B2[,V7]))]
prom_Input_PV_KO_B3[, cpm := (10**6) * (prom_Input_PV_KO_B3[,V7]/sum(prom_Input_PV_KO_B3[,V7]))]


prom_H3K27ac_Total_WTKO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_WTKO_B2_CTCAGGGC_chipCounts.bed")
prom_H3K27ac_Total_WTKO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_WTKO_B3_CTATCTTG_chipCounts.bed")
prom_H3K27ac_Total_KO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_KO_B2_TGCATGTA_chipCounts.bed")
prom_H3K27ac_Total_KO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_KO_B3_ATGTCGGC_chipCounts.bed")

prom_H3K27ac_Total_WTKO_B2[, cpm := (10**6) * (prom_H3K27ac_Total_WTKO_B2[,V7]/sum(prom_H3K27ac_Total_WTKO_B2[,V7]))]
prom_H3K27ac_Total_WTKO_B3[, cpm := (10**6) * (prom_H3K27ac_Total_WTKO_B3[,V7]/sum(prom_H3K27ac_Total_WTKO_B3[,V7]))]
prom_H3K27ac_Total_KO_B2[, cpm := (10**6) * (prom_H3K27ac_Total_KO_B2[,V7]/sum(prom_H3K27ac_Total_KO_B2[,V7]))]
prom_H3K27ac_Total_KO_B3[, cpm := (10**6) * (prom_H3K27ac_Total_KO_B3[,V7]/sum(prom_H3K27ac_Total_KO_B3[,V7]))]

prom_Input_Total_WTKO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_Input_Total_WTKO_B2_ATCGTCTC_chipCounts.bed")
prom_Input_Total_WTKO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_Input_Total_WTKO_B3_CCTTAGGT_chipCounts.bed")
prom_Input_Total_KO_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_Input_Total_KO_B2_CCAACACT_chipCounts.bed")
prom_Input_Total_KO_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_KO/ensgene_mm9_1kb_promoterWindows_Input_Total_KO_B3_AACCGAAC_chipCounts.bed")

prom_Input_Total_WTKO_B2[, cpm := (10**6) * (prom_Input_Total_WTKO_B2[,V7]/sum(prom_Input_Total_WTKO_B2[,V7]))]
prom_Input_Total_WTKO_B3[, cpm := (10**6) * (prom_Input_Total_WTKO_B3[,V7]/sum(prom_Input_Total_WTKO_B3[,V7]))]
prom_Input_Total_KO_B2[, cpm := (10**6) * (prom_Input_Total_KO_B2[,V7]/sum(prom_Input_Total_KO_B2[,V7]))]
prom_Input_Total_KO_B3[, cpm := (10**6) * (prom_Input_Total_KO_B3[,V7]/sum(prom_Input_Total_KO_B3[,V7]))]

prom_H3K27ac_Input_Total_WTKO = cbind(prom_H3K27ac_Total_WTKO_B2[, .(V4, cpm)], 
                                      prom_H3K27ac_Total_WTKO_B3[, .(cpm)],
                                      prom_Input_Total_WTKO_B2[, .(cpm)],
                                      prom_Input_Total_WTKO_B3[, .(cpm)])
names(prom_H3K27ac_Input_Total_WTKO) = c("Gene", "WT1_cpm", "WT2_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_Total_WTKO[, WT_avg_cpm := (WT1_cpm + WT2_cpm)/2]
prom_H3K27ac_Input_Total_WTKO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_Input_Total_KO = cbind(prom_H3K27ac_Total_KO_B2[, .(V4, cpm)], 
                                      prom_H3K27ac_Total_KO_B3[, .(cpm)],
                                      prom_Input_Total_KO_B2[, .(cpm)],
                                      prom_Input_Total_KO_B3[, .(cpm)])
names(prom_H3K27ac_Input_Total_KO) = c("Gene", "KO1_cpm", "KO2_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_Total_KO[, KO_avg_cpm := (KO1_cpm + KO2_cpm)/2]
prom_H3K27ac_Input_Total_KO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_Total_KO[, log2_cpm_ratio_pseud := log2((KO_avg_cpm + 1)/(Input_avg_cpm+1))]


prom_H3K27ac_Input_PV_WTKO = cbind(prom_H3K27ac_PV_WTKO_B1[, .(V4, cpm)],
                                   prom_H3K27ac_PV_WTKO_B2[, .(cpm)], 
                                   prom_H3K27ac_PV_WTKO_B3[, .(cpm)],
                                   prom_Input_PV_WTKO_B2[, .(cpm)],
                                   prom_Input_PV_WTKO_B3[, .(cpm)])
names(prom_H3K27ac_Input_PV_WTKO) = c("Gene", "WT1_cpm", "WT2_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_PV_WTKO[, WT_avg_cpm := (WT1_cpm + WT2_cpm + WT3_cpm)/3]
prom_H3K27ac_Input_PV_WTKO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_Input_PV_KO = cbind(prom_H3K27ac_PV_KO_B1[, .(V4, cpm)],
                                   prom_H3K27ac_PV_KO_B2[, .(cpm)], 
                                   prom_H3K27ac_PV_KO_B3[, .(cpm)],
                                   prom_Input_PV_KO_B2[, .(cpm)],
                                   prom_Input_PV_KO_B3[, .(cpm)])
names(prom_H3K27ac_Input_PV_KO) = c("Gene", "KO1_cpm", "KO2_cpm", "KO3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_PV_KO[, KO_avg_cpm := (KO1_cpm + KO2_cpm + KO3_cpm)/3]
prom_H3K27ac_Input_PV_KO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_PV_KO[, log2_cpm_ratio_pseud := log2((KO_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WT_B1_AAGGAGTC_chipCounts.bed")
prom_H3K27ac_PV_WTtg_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WTtg_B2_CGTCTAAC_chipCounts.bed")
prom_H3K27ac_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_WTtg_B3_TCGCGAGG_chipCounts.bed")
prom_H3K27ac_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_Tg_B1_GCTGAGTC_chipCounts.bed")
prom_H3K27ac_PV_TG_B2 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_Tg_B2_GACTACGA_chipCounts.bed")
prom_H3K27ac_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_PV_Tg_B3_TGCATGTA_chipCounts.bed")

prom_H3K27ac_PV_WTtg_B1[, cpm := (10**6) * (prom_H3K27ac_PV_WTtg_B1[,V7]/sum(prom_H3K27ac_PV_WTtg_B1[,V7]))]
prom_H3K27ac_PV_WTtg_B2[, cpm := (10**6) * (prom_H3K27ac_PV_WTtg_B2[,V7]/sum(prom_H3K27ac_PV_WTtg_B2[,V7]))]
prom_H3K27ac_PV_WTtg_B3[, cpm := (10**6) * (prom_H3K27ac_PV_WTtg_B3[,V7]/sum(prom_H3K27ac_PV_WTtg_B3[,V7]))]
prom_H3K27ac_PV_TG_B1[, cpm := (10**6) * (prom_H3K27ac_PV_TG_B1[,V7]/sum(prom_H3K27ac_PV_TG_B1[,V7]))]
prom_H3K27ac_PV_TG_B2[, cpm := (10**6) * (prom_H3K27ac_PV_TG_B2[,V7]/sum(prom_H3K27ac_PV_TG_B2[,V7]))]
prom_H3K27ac_PV_TG_B3[, cpm := (10**6) * (prom_H3K27ac_PV_TG_B3[,V7]/sum(prom_H3K27ac_PV_TG_B3[,V7]))]

prom_Input_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_Input_PV_WT_B1_TCGCGAGG_chipCounts.bed")
prom_Input_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_Input_PV_WTtg_B3_AACTCGGA_chipCounts.bed")
prom_Input_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_Input_PV_Tg_B1_CTTCCCAT_chipCounts.bed")
prom_Input_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/PV_TG/ensgene_mm9_1kb_promoterWindows_Input_PV_Tg_B3_ACTCCTAC_chipCounts.bed")

prom_Input_PV_WTtg_B1[, cpm := (10**6) * (prom_Input_PV_WTtg_B1[,V7]/sum(prom_Input_PV_WTtg_B1[,V7]))]
prom_Input_PV_WTtg_B3[, cpm := (10**6) * (prom_Input_PV_WTtg_B3[,V7]/sum(prom_Input_PV_WTtg_B3[,V7]))]
prom_Input_PV_TG_B1[, cpm := (10**6) * (prom_Input_PV_TG_B1[,V7]/sum(prom_Input_PV_TG_B1[,V7]))]
prom_Input_PV_TG_B3[, cpm := (10**6) * (prom_Input_PV_TG_B3[,V7]/sum(prom_Input_PV_TG_B3[,V7]))]


prom_H3K27ac_Input_PV_WTtg = cbind(prom_H3K27ac_PV_WTtg_B1[, .(V4, cpm)],
                                   prom_H3K27ac_PV_WTtg_B2[, .(cpm)],
                                   prom_H3K27ac_PV_WTtg_B3[, .(cpm)],
                                   prom_Input_PV_WTtg_B1[, .(cpm)],
                                   prom_Input_PV_WTtg_B3[, .(cpm)])
names(prom_H3K27ac_Input_PV_WTtg) = c("Gene", "WT1_cpm", "WT2_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_PV_WTtg[, WT_avg_cpm := (WT1_cpm + WT2_cpm + WT3_cpm)/3]
prom_H3K27ac_Input_PV_WTtg[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_PV_WTtg[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_Input_PV_TG = cbind(prom_H3K27ac_PV_TG_B1[, .(V4, cpm)],
                                   prom_H3K27ac_PV_TG_B2[, .(cpm)],
                                   prom_H3K27ac_PV_TG_B3[, .(cpm)],
                                   prom_Input_PV_TG_B1[, .(cpm)],
                                   prom_Input_PV_TG_B3[, .(cpm)])
names(prom_H3K27ac_Input_PV_TG) = c("Gene", "TG1_cpm", "TG2_cpm", "TG3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_PV_TG[, TG_avg_cpm := (TG1_cpm + TG2_cpm + TG3_cpm)/3]
prom_H3K27ac_Input_PV_TG[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_PV_TG[, log2_cpm_ratio_pseud := log2((TG_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_Total_WTtg_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_WT_B1_CTCAGGGC_chipCounts.bed")
prom_H3K27ac_Total_WTtg_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_WTtg_B3_ATGTCGGC_chipCounts.bed")
prom_H3K27ac_Total_TG_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_Tg_B1_TGCATGTA_chipCounts.bed")
prom_H3K27ac_Total_TG_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_H3K27ac_Total_Tg_B3_CTATCTTG_chipCounts.bed")

prom_H3K27ac_Total_WTtg_B1[, cpm := (10**6) * (prom_H3K27ac_Total_WTtg_B1[,V7]/sum(prom_H3K27ac_Total_WTtg_B1[,V7]))]
prom_H3K27ac_Total_WTtg_B3[, cpm := (10**6) * (prom_H3K27ac_Total_WTtg_B3[,V7]/sum(prom_H3K27ac_Total_WTtg_B3[,V7]))]
prom_H3K27ac_Total_TG_B1[, cpm := (10**6) * (prom_H3K27ac_Total_TG_B1[,V7]/sum(prom_H3K27ac_Total_TG_B1[,V7]))]
prom_H3K27ac_Total_TG_B3[, cpm := (10**6) * (prom_H3K27ac_Total_TG_B3[,V7]/sum(prom_H3K27ac_Total_TG_B3[,V7]))]

prom_Input_Total_WTtg_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_Input_Total_WT_B1_ATGTCGGC_chipCounts.bed")
prom_Input_Total_WTtg_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_Input_Total_WTtg_B3_ACAACAGC_chipCounts.bed")
prom_Input_Total_TG_B1 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_Input_Total_Tg_B1_CTATCTTG_chipCounts.bed")
prom_Input_Total_TG_B3 = fread("HG_lab/Mati/GabelLab/genesets/Pv_H3K27ac_chipCounts/Total_TG/ensgene_mm9_1kb_promoterWindows_Input_Total_Tg_B3_ACCATCCT_chipCounts.bed")

prom_Input_Total_WTtg_B1[, cpm := (10**6) * (prom_Input_Total_WTtg_B1[,V7]/sum(prom_Input_Total_WTtg_B1[,V7]))]
prom_Input_Total_WTtg_B3[, cpm := (10**6) * (prom_Input_Total_WTtg_B3[,V7]/sum(prom_Input_Total_WTtg_B3[,V7]))]
prom_Input_Total_TG_B1[, cpm := (10**6) * (prom_Input_Total_TG_B1[,V7]/sum(prom_Input_Total_TG_B1[,V7]))]
prom_Input_Total_TG_B3[, cpm := (10**6) * (prom_Input_Total_TG_B3[,V7]/sum(prom_Input_Total_TG_B3[,V7]))]

prom_H3K27ac_Input_Total_WTtg = cbind(prom_H3K27ac_Total_WTtg_B1[, .(V4, cpm)], 
                                    prom_H3K27ac_Total_WTtg_B3[, .(cpm)],
                                    prom_Input_Total_WTtg_B1[, .(cpm)],
                                    prom_Input_Total_WTtg_B3[, .(cpm)])
names(prom_H3K27ac_Input_Total_WTtg) = c("Gene", "WT1_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_Total_WTtg[, WT_avg_cpm := (WT1_cpm + WT3_cpm)/2]
prom_H3K27ac_Input_Total_WTtg[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_Total_WTtg[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

prom_H3K27ac_Input_Total_TG = cbind(prom_H3K27ac_Total_TG_B1[, .(V4, cpm)], 
                                    prom_H3K27ac_Total_TG_B3[, .(cpm)],
                                    prom_Input_Total_TG_B1[, .(cpm)],
                                    prom_Input_Total_TG_B3[, .(cpm)])
names(prom_H3K27ac_Input_Total_TG) = c("Gene", "TG1_cpm", "TG3_cpm", "Input1_cpm", "Input2_cpm")
prom_H3K27ac_Input_Total_TG[, TG_avg_cpm := (TG1_cpm + TG3_cpm)/2]
prom_H3K27ac_Input_Total_TG[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
prom_H3K27ac_Input_Total_TG[, log2_cpm_ratio_pseud := log2((TG_avg_cpm + 1)/(Input_avg_cpm+1))]


png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top200_DEG_1kb_promoterWindows_PVMeCP2WT_vs_TotalMeCP2WT_log2_H3K27ac_ChIP_over_input.png")
smoothScatter(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="gold", pch=16, cex=1)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="forestgreen", pch=16, cex=1)
abline(coef = c(0,1))
#abline(lm(prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud]~prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud]))
legend("topleft", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), cex=1,
       col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top200_DEG_1kb_promoterWindows_PVMeCP2WT_vs_TotalMeCP2WT_log2_H3K27ac_ChIP_over_input.eps")
smoothScatter(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="gold", pch=16, cex=1)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top200_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="forestgreen", pch=16, cex=1)
abline(coef = c(0,1))
#abline(lm(prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud]~prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud]))
legend("topleft", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), cex=1,
       col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top100_DEG_q0.01_1kb_promoterWindows_PVMeCP2WTKO_vs_TotalMeCP2WTKO_log2_H3K27ac_ChIP_over_input.png")
smoothScatter(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WTKO log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="gold", pch=16, cex=1)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="forestgreen", pch=16, cex=1)
abline(coef = c(0,1))
#abline(lm(prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud]~prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud]))
legend("topleft", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), cex=1,
       col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top100_DEG_q0.01_1kb_promoterWindows_PVMeCP2WTKO_vs_TotalMeCP2WTKO_log2_H3K27ac_ChIP_over_input.eps")
smoothScatter(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WTKO log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="gold", pch=16, cex=1)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="forestgreen", pch=16, cex=1)
abline(coef = c(0,1))
#abline(lm(prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud]~prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud]))
legend("topleft", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), cex=1,
       col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

#Read in files with Pv MeCP2 KO H3K27ac ChIP counts in cCREs
cCREs_1500bp_H3K27ac_PV_WTKO_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTko_B1_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTKO_B2_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTKO_B3_CTTCCCAT_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B1_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B2_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B3_TCGCGAGG_chipCounts_mm9.bed")

cCREs_1500bp_H3K27ac_PV_WTKO_B1[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B1[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B1[,V5]))]
cCREs_1500bp_H3K27ac_PV_WTKO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B2[,V5]))]
cCREs_1500bp_H3K27ac_PV_WTKO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B3[,V5]))]
cCREs_1500bp_H3K27ac_PV_KO_B1[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B1[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B1[,V5]))]
cCREs_1500bp_H3K27ac_PV_KO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B2[,V5]))]
cCREs_1500bp_H3K27ac_PV_KO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B3[,V5]))]

cCREs_1500bp_Input_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_WTKO_B2_CGACCTAA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_WTKO_B3_CTCGAACA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_KO_B2_TACATCGG_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_KO_B3_TCTAGGAG_chipCounts_mm9.bed")

cCREs_1500bp_Input_PV_WTKO_B2[, cpm := (10**6) * (cCREs_1500bp_Input_PV_WTKO_B2[,V5]/sum(cCREs_1500bp_Input_PV_WTKO_B2[,V5]))]
cCREs_1500bp_Input_PV_WTKO_B3[, cpm := (10**6) * (cCREs_1500bp_Input_PV_WTKO_B3[,V5]/sum(cCREs_1500bp_Input_PV_WTKO_B3[,V5]))]
cCREs_1500bp_Input_PV_KO_B2[, cpm := (10**6) * (cCREs_1500bp_Input_PV_KO_B2[,V5]/sum(cCREs_1500bp_Input_PV_KO_B2[,V5]))]
cCREs_1500bp_Input_PV_KO_B3[, cpm := (10**6) * (cCREs_1500bp_Input_PV_KO_B3[,V5]/sum(cCREs_1500bp_Input_PV_KO_B3[,V5]))]

cCREs_1500bp_H3K27ac_Input_PV_WTKO = cbind(cCREs_1500bp_H3K27ac_PV_WTKO_B1[, .(V4, cpm)],
                                           cCREs_1500bp_H3K27ac_PV_WTKO_B2[, .(cpm)], 
                                           cCREs_1500bp_H3K27ac_PV_WTKO_B3[, .(cpm)],
                                           cCREs_1500bp_Input_PV_WTKO_B2[, .(cpm)],
                                           cCREs_1500bp_Input_PV_WTKO_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_PV_WTKO) = c("cCRE_label", "WT1_cpm", "WT2_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_PV_WTKO[, WT_avg_cpm := (WT1_cpm + WT2_cpm + WT3_cpm)/3]
cCREs_1500bp_H3K27ac_Input_PV_WTKO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

cCREs_1500bp_H3K27ac_Input_PV_KO = cbind(cCREs_1500bp_H3K27ac_PV_KO_B1[, .(V4, cpm)],
                                           cCREs_1500bp_H3K27ac_PV_KO_B2[, .(cpm)], 
                                           cCREs_1500bp_H3K27ac_PV_KO_B3[, .(cpm)],
                                           cCREs_1500bp_Input_PV_KO_B2[, .(cpm)],
                                           cCREs_1500bp_Input_PV_KO_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_PV_KO) = c("cCRE_label", "KO1_cpm", "KO2_cpm", "KO3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_PV_KO[, KO_avg_cpm := (KO1_cpm + KO2_cpm + KO3_cpm)/3]
cCREs_1500bp_H3K27ac_Input_PV_KO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_PV_KO[, log2_cpm_ratio_pseud := log2((KO_avg_cpm + 1)/(Input_avg_cpm+1))]


cCREs_1500bp_H3K27ac_Total_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_WTKO_B2_CTCAGGGC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_WTKO_B3_CTATCTTG_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_KO_B2_TGCATGTA_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_KO_B3_ATGTCGGC_chipCounts_mm9.bed")

cCREs_1500bp_H3K27ac_Total_WTKO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_Total_WTKO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_Total_WTKO_B2[,V5]))]
cCREs_1500bp_H3K27ac_Total_WTKO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_Total_WTKO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_Total_WTKO_B3[,V5]))]
cCREs_1500bp_H3K27ac_Total_KO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_Total_KO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_Total_KO_B2[,V5]))]
cCREs_1500bp_H3K27ac_Total_KO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_Total_KO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_Total_KO_B3[,V5]))]

cCREs_1500bp_Input_Total_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_Input_Total_WTKO_B2_ATCGTCTC_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_Input_Total_WTKO_B3_CCTTAGGT_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_Input_Total_KO_B2_CCAACACT_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_KO/mousebrain_union_cCREs_1500bpWindows_Input_Total_KO_B3_AACCGAAC_chipCounts_mm9.bed")

cCREs_1500bp_Input_Total_WTKO_B2[, cpm := (10**6) * (cCREs_1500bp_Input_Total_WTKO_B2[,V5]/sum(cCREs_1500bp_Input_Total_WTKO_B2[,V5]))]
cCREs_1500bp_Input_Total_WTKO_B3[, cpm := (10**6) * (cCREs_1500bp_Input_Total_WTKO_B3[,V5]/sum(cCREs_1500bp_Input_Total_WTKO_B3[,V5]))]
cCREs_1500bp_Input_Total_KO_B2[, cpm := (10**6) * (cCREs_1500bp_Input_Total_KO_B2[,V5]/sum(cCREs_1500bp_Input_Total_KO_B2[,V5]))]
cCREs_1500bp_Input_Total_KO_B3[, cpm := (10**6) * (cCREs_1500bp_Input_Total_KO_B3[,V5]/sum(cCREs_1500bp_Input_Total_KO_B3[,V5]))]

cCREs_1500bp_H3K27ac_Input_Total_WTKO = cbind(cCREs_1500bp_H3K27ac_Total_WTKO_B2[, .(V4,cpm)], 
                                           cCREs_1500bp_H3K27ac_Total_WTKO_B3[, .(cpm)],
                                           cCREs_1500bp_Input_Total_WTKO_B2[, .(cpm)],
                                           cCREs_1500bp_Input_Total_WTKO_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_Total_WTKO) = c("cCRE_label", "WT1_cpm", "WT2_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_Total_WTKO[, WT_avg_cpm := (WT1_cpm + WT2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_WTKO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

cCREs_1500bp_H3K27ac_Input_Total_KO = cbind(cCREs_1500bp_H3K27ac_Total_KO_B2[, .(V4,cpm)], 
                                              cCREs_1500bp_H3K27ac_Total_KO_B3[, .(cpm)],
                                              cCREs_1500bp_Input_Total_KO_B2[, .(cpm)],
                                              cCREs_1500bp_Input_Total_KO_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_Total_KO) = c("cCRE_label", "KO1_cpm", "KO2_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_Total_KO[, KO_avg_cpm := (KO1_cpm + KO2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_KO[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_KO[, log2_cpm_ratio_pseud := log2((KO_avg_cpm + 1)/(Input_avg_cpm+1))]

cCREs_1500bp_H3K27ac_Total_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_WT_B1_CTCAGGGC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_WTtg_B3_ATGTCGGC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_Tg_B1_TGCATGTA_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_Total_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_Total_Tg_B3_CTATCTTG_chipCounts_mm9.bed")

cpm_function(cCREs_1500bp_H3K27ac_Total_WTtg_B1, "V5")
cpm_function(cCREs_1500bp_H3K27ac_Total_WTtg_B3, "V5")
cpm_function(cCREs_1500bp_H3K27ac_Total_TG_B1, "V5")
cpm_function(cCREs_1500bp_H3K27ac_Total_TG_B3, "V5")

cCREs_1500bp_Input_Total_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_Input_Total_WT_B1_ATGTCGGC_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_Input_Total_WTtg_B3_ACAACAGC_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_Input_Total_Tg_B1_CTATCTTG_chipCounts_mm9.bed")
cCREs_1500bp_Input_Total_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/Total_TG/mousebrain_union_cCREs_1500bpWindows_Input_Total_Tg_B3_ACCATCCT_chipCounts_mm9.bed")

cpm_function(cCREs_1500bp_Input_Total_WTtg_B1, "V5")
cpm_function(cCREs_1500bp_Input_Total_WTtg_B3, "V5")
cpm_function(cCREs_1500bp_Input_Total_TG_B1, "V5")
cpm_function(cCREs_1500bp_Input_Total_TG_B3, "V5")

cCREs_1500bp_H3K27ac_Input_Total_WTtg = cbind(cCREs_1500bp_H3K27ac_Total_WTtg_B1[, .(V4,cpm)], 
                                              cCREs_1500bp_H3K27ac_Total_WTtg_B3[, .(cpm)],
                                              cCREs_1500bp_Input_Total_WTtg_B1[, .(cpm)],
                                              cCREs_1500bp_Input_Total_WTtg_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_Total_WTtg) = c("cCRE_label", "WT1_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_Total_WTtg[, WT_avg_cpm := (WT1_cpm + WT3_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_WTtg[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_WTtg[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

cCREs_1500bp_H3K27ac_Input_Total_TG = cbind(cCREs_1500bp_H3K27ac_Total_TG_B1[, .(V4,cpm)], 
                                              cCREs_1500bp_H3K27ac_Total_TG_B3[, .(cpm)],
                                              cCREs_1500bp_Input_Total_TG_B1[, .(cpm)],
                                              cCREs_1500bp_Input_Total_TG_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_Total_TG) = c("cCRE_label", "TG1_cpm", "TG3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_Total_TG[, TG_avg_cpm := (TG1_cpm + TG3_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_TG[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_Total_TG[, log2_cpm_ratio_pseud := log2((TG_avg_cpm + 1)/(Input_avg_cpm+1))]


cCREs_1500bp_H3K27ac_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WT_B1_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTtg_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTtg_B2_CGTCTAAC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTtg_B3_TCGCGAGG_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B1_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B2_GACTACGA_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B3_TGCATGTA_chipCounts_mm9.bed")

cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B1, "V5")
cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B2, "V5")
cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B3, "V5")
cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B1, "V5")
cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B2, "V5")
cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B3, "V5")

cCREs_1500bp_Input_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_WT_B1_TCGCGAGG_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_WTtg_B3_AACTCGGA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_Tg_B1_CTTCCCAT_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_Tg_B3_ACTCCTAC_chipCounts_mm9.bed")

cpm_function(cCREs_1500bp_Input_PV_WTtg_B1, "V5")
cpm_function(cCREs_1500bp_Input_PV_WTtg_B3, "V5")
cpm_function(cCREs_1500bp_Input_PV_TG_B1, "V5")
cpm_function(cCREs_1500bp_Input_PV_TG_B3, "V5")

cCREs_1500bp_H3K27ac_Input_PV_WTtg = cbind(cCREs_1500bp_H3K27ac_PV_WTtg_B1[, .(V4,cpm)], 
                                           cCREs_1500bp_H3K27ac_PV_WTtg_B2[, .(cpm)],
                                           cCREs_1500bp_H3K27ac_PV_WTtg_B3[, .(cpm)],
                                           cCREs_1500bp_Input_PV_WTtg_B1[, .(cpm)],
                                           cCREs_1500bp_Input_PV_WTtg_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_PV_WTtg) = c("cCRE_label", "WT1_cpm", "WT2_cpm", "WT3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_PV_WTtg[, WT_avg_cpm := (WT1_cpm + WT2_cpm + WT3_cpm)/3]
cCREs_1500bp_H3K27ac_Input_PV_WTtg[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_PV_WTtg[, log2_cpm_ratio_pseud := log2((WT_avg_cpm + 1)/(Input_avg_cpm+1))]

cCREs_1500bp_H3K27ac_Input_PV_TG = cbind(cCREs_1500bp_H3K27ac_PV_TG_B1[, .(V4,cpm)],
                                         cCREs_1500bp_H3K27ac_PV_TG_B2[, .(cpm)],
                                         cCREs_1500bp_H3K27ac_PV_TG_B3[, .(cpm)],
                                         cCREs_1500bp_Input_PV_TG_B1[, .(cpm)],
                                         cCREs_1500bp_Input_PV_TG_B3[, .(cpm)])
names(cCREs_1500bp_H3K27ac_Input_PV_TG) = c("cCRE_label", "TG1_cpm", "TG2_cpm", "TG3_cpm", "Input1_cpm", "Input2_cpm")
cCREs_1500bp_H3K27ac_Input_PV_TG[, TG_avg_cpm := (TG1_cpm + TG2_cpm + TG3_cpm)/3]
cCREs_1500bp_H3K27ac_Input_PV_TG[, Input_avg_cpm := (Input1_cpm + Input2_cpm)/2]
cCREs_1500bp_H3K27ac_Input_PV_TG[, log2_cpm_ratio_pseud := log2((TG_avg_cpm + 1)/(Input_avg_cpm+1))]

png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top200_DEG_linked_PVGA_cCREs_1500bpWindows_PVMeCP2WT_vs_TotalMeCP2WT_log2_H3K27ac_ChIP_over_input.png")
smoothScatter(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), pch=16, cex=0.3, col="lightgray")
par(new=TRUE)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="gold", pch=16, cex=0.5)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="forestgreen", pch=16, cex=0.5)
abline(coef = c(0,1))
legend("topleft", title="PV cCREs linked to:", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), 
                               cex=1, col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top200_DEG_linked_nonPVGA_cCREs_1500bpWindows_PVMeCP2WT_vs_TotalMeCP2WT_log2_H3K27ac_ChIP_over_input.png")
smoothScatter(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), pch=16, cex=0.3, col="lightgray")
par(new=TRUE)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="gold", pch=16, cex=0.5)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WT\n log2 H3K27ac ChIP/Input", ylab="Total WT log2 H3K27ac ChIP/Input", xlim=c(0, 8.5), ylim=c(0, 8.5), col="forestgreen", pch=16, cex=0.5)
abline(coef = c(0,1))
legend("topleft", title="Non-PV cCREs linked to:", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), 
       cex=1, col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top100_DEG_q0.01_linked_union_cCREs_1500bpWindows_PVMeCP2WTKO_vs_TotalMeCP2WTKO_log2_H3K27ac_ChIP_over_input.png")
smoothScatter(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WTKO log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2.5, 4), pch=16, cex=0.3, col="lightgray")
par(new=TRUE)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2, 4), col="gold", pch=16, cex=0.5)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2, 4), col="forestgreen", pch=16, cex=0.5)
abline(coef = c(0,1))
legend("topleft", title="cCREs linked to:", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), 
       cex=1, col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/Huntley2020_Exc_Pv_top100_DEG_q0.01_linked_union_cCREs_1500bpWindows_PVMeCP2WTKO_vs_TotalMeCP2WTKO_log2_H3K27ac_ChIP_over_input.eps")
smoothScatter(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WTKO log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2.5, 4), pch=16, cex=0.3, col="lightgray")
par(new=TRUE)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2, 4), col="gold", pch=16, cex=0.5)
points(x=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], y=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-2.5, 4), ylim=c(-2, 4), col="forestgreen", pch=16, cex=0.5)
abline(coef = c(0,1))
legend("topleft", title="cCREs linked to:", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), 
       cex=1, col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[,cCRE_label], log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")
dev.off()
#boxplots of H3K27ac enrichment
enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))
enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt[, val := as.numeric(val)]
ggplot(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (PV WTKO / Total WTKO)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_PV_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_PV_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))
enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt[, val := as.numeric(val)]

ggplot(enrich_PV_KO_Total_WTKO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (PV KO / Total WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2WT_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2WT_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_Total_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_PV_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_Total_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_PV_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_Total_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_Total_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt = enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt = enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))

ggplot(enrich_Total_PV_KO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (Total/PV KO)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())


enrich_Total_WTKO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="Exc"),
  cbind(val=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="PV"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="PV")))
enrich_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_Total_WTKO_prom_cCRE_H3K27ac_dt = enrich_Total_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))

ggplot(enrich_Total_WTKO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-4,4))+
  ylab("Total WT Log2 H3K27ac ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())

enrich_PV_WTKO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="Exc"),
  cbind(val=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="PV"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="PV")))
enrich_PV_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_WTKO_prom_cCRE_H3K27ac_dt = enrich_PV_WTKO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))

ggplot(enrich_PV_WTKO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-4,4))+
  ylab("PV WT Log2 H3K27ac ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())


enrich_PV_KO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_PV_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="Exc"),
  cbind(val=prom_H3K27ac_Input_PV_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], element="Promoters", enriched_cell="PV"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud], element="Linked cCREs", enriched_cell="PV")))
enrich_PV_KO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_KO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))

ggplot(enrich_PV_KO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-4,4))+
  ylab("PV KO Log2 H3K27ac ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())


enrich_PV_WTKO_prom_cCRE_H3K27ac_dt_all  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], element="Promoters"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], element="Linked cCREs")))


ggplot(enrich_PV_WTKO_prom_cCRE_H3K27ac_dt_all, aes(x = element, y = as.numeric(val), fill=element))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-6,6))+
  ylab("PV WT Log2 H3K27ac ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())

enrich_Total_WTKO_prom_cCRE_H3K27ac_dt_all  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], element="Promoters"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], element="Linked cCREs")))


ggplot(enrich_Total_WTKO_prom_cCRE_H3K27ac_dt_all, aes(x = element, y = as.numeric(val), fill=element))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-6,6))+
  ylab("Total WT Log2 H3K27ac ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())

enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], WT_avg_cpm], element="Promoters", enriched_cell="Exc", genotype="Total WT"),
  cbind(val=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], WT_avg_cpm], element="Promoters", enriched_cell="PV", genotype="Total WT"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], WT_avg_cpm], element="Linked cCREs", enriched_cell="Exc", genotype="Total WT"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_Total_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], WT_avg_cpm], element="Linked cCREs", enriched_cell="PV", genotype="Total WT")))
enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))


enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], WT_avg_cpm], element="Promoters", enriched_cell="Exc", genotype="PV WT"),
  cbind(val=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], WT_avg_cpm], element="Promoters", enriched_cell="PV", genotype="PV WT"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], WT_avg_cpm], element="Linked cCREs", enriched_cell="Exc", genotype="PV WT"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_WTKO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], WT_avg_cpm], element="Linked cCREs", enriched_cell="PV", genotype="PV WT")))
enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))

enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt  = data.table(rbind(
  cbind(val=prom_H3K27ac_Input_PV_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], KO_avg_cpm], element="Promoters", enriched_cell="Exc",  genotype="PV KO"),
  cbind(val=prom_H3K27ac_Input_PV_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], KO_avg_cpm], element="Promoters", enriched_cell="PV",  genotype="PV KO"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], KO_avg_cpm], element="Linked cCREs", enriched_cell="Exc",  genotype="PV KO"),
  cbind(val=cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], KO_avg_cpm], element="Linked cCREs", enriched_cell="PV",  genotype="PV KO")))
enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt = enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))


ggplot(enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt[element=="Promoters", ], aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  #coord_cartesian(ylim=c(-4,4))+
  ylab("PV KO Log2 H3K27ac ChIP") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())

boxplot(as.numeric(val)~enriched_cell, data=enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt[element=="Promoters"], col=(c("gold","forestgreen")), main="Promoters", xlab="Enriched cell type", ylab="PV WT Log2 H3K27ac ChIP", notch=TRUE)
boxplot(as.numeric(val)~enriched_cell, data=enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt[element=="Linked cCREs"], col=(c("gold","forestgreen")), main="Linked cCREs", xlab="Enriched cell type", ylab="PV WT Log2 H3K27ac ChIP", notch=TRUE)


enrich_Total_PV_rawCPM_dt = rbind(enrich_Total_WTKO_prom_cCRE_H3K27ac_rawCPM_dt, enrich_PV_WTKO_prom_cCRE_H3K27ac_rawCPM_dt, enrich_PV_KO_prom_cCRE_H3K27ac_rawCPM_dt)
enrich_Total_PV_rawCPM_dt = enrich_Total_PV_rawCPM_dt %>% mutate(genotype = factor(genotype, levels=c("Total WT", "PV WT", "PV KO")))

ggplot(enrich_Total_PV_rawCPM_dt[element=="Promoters"], aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  ggtitle("Promoters")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~genotype,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(0,230))+
  ylab("Log2 H3K27ac ChIP cpm") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())

ggplot(enrich_Total_PV_rawCPM_dt[element=="Linked cCREs"], aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  ggtitle("Linked cCREs")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~genotype,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(0,32))+
  ylab("Log2 H3K27ac ChIP cpm") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())


##PV TG vs total TG
enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_PV_TG[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_TG[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_PV_TG[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_TG[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_TG[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_TG[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_TG[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_TG[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt = enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt = enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))
enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[, val := as.numeric(val)]

ggplot(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (PV TG / Total TG)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##PV KO vs total KO
enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_PV_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_KO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_PV_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_KO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_KO[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_KO[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt = enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))
enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt[, val := as.numeric(val)]
ggplot(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (PV KO / Total KO)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##PV WTtg vs total WTtg
enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt  = data.table(rbind(
  cbind(val=(prom_H3K27ac_Input_PV_WTtg[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTtg[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="Exc"),
  cbind(val=(prom_H3K27ac_Input_PV_WTtg[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud] - prom_H3K27ac_Input_Total_WTtg[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud]), element="Promoters", enriched_cell="PV"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_WTtg[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTtg[cCRE_label %in% Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="Exc"),
  cbind(val=(cCREs_1500bp_H3K27ac_Input_PV_WTtg[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud] - cCREs_1500bp_H3K27ac_Input_Total_WTtg[cCRE_label %in% Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9[, V4], log2_cpm_ratio_pseud]), element="Linked cCREs", enriched_cell="PV")))
enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt = enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt %>% mutate(element = factor(element, levels=c("Promoters", "Linked cCREs")))
enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt = enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt %>% mutate(enriched_cell = factor(enriched_cell, levels=c("Exc", "PV")))
enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt[, val := as.numeric(val)]

ggplot(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt, aes(x = enriched_cell, y = as.numeric(val), fill=enriched_cell))+ 
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  facet_grid(.~element,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels to the cell-type and element variables
  ) + 
  coord_cartesian(ylim=c(-5,5))+
  ylab("Log2 H3K27ac FD (PV WTtg / Total WTtg)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=15), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_linkedcCREs_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Promoters"), as.numeric(val)], enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Promoters"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Promoters"), as.numeric(val)], enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Promoters"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Linked cCREs"), as.numeric(val)], enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Linked cCREs"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Linked cCREs"), as.numeric(val)], enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Linked cCREs"), as.numeric(val)])$p.value


wilcox.test(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Promoters"), as.numeric(val)], enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Promoters"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Promoters"), as.numeric(val)], enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Promoters"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Linked cCREs"), as.numeric(val)], enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[(enriched_cell=="Exc") & (element=="Linked cCREs"), as.numeric(val)])$p.value
wilcox.test(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Linked cCREs"), as.numeric(val)], enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[(enriched_cell=="PV") & (element=="Linked cCREs"), as.numeric(val)])$p.value

smoothScatter(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WTKO log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), pch=16, cex=0.6, col="lightgray", bty="n")
par(new=TRUE)
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Exc_over_Pv_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="gold", pch=16, cex=1, bty="n")
points(x=prom_H3K27ac_Input_PV_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], y=prom_H3K27ac_Input_Total_WTKO[Gene %in% Pv_over_Exc_genes_top100_q0.01_mm9[, V4], log2_cpm_ratio_pseud], xlab="PV WTKO\n log2 H3K27ac ChIP/Input", ylab="Total WTKO log2 H3K27ac ChIP/Input", xlim=c(-4, 4), ylim=c(-4, 4), col="forestgreen", pch=16, cex=1, bty="n")
abline(coef = c(0,1))
#abline(lm(prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud]~prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud]))
legend("topleft", legend=c('EXC-enriched genes', 'PV-enriched genes', 'All genes'), cex=1,
       col=c('gold', 'forestgreen', 'lightgray'), pch=c(16,16), bty="n")
legend("bottomright", legend = paste0("rho=", round(cor(prom_H3K27ac_Input_PV_WTKO[, log2_cpm_ratio_pseud], prom_H3K27ac_Input_Total_WTKO[, log2_cpm_ratio_pseud], method="spearman", use="complete.obs"), 3)), bty = "n")

elem_boxplot_func <- function(data_table, elem, title, x_label, x_min=-5, x_max=5){
  ggplot(data_table[element==elem], aes(y = enriched_cell, x = as.numeric(val), fill=enriched_cell))+ 
    ggtitle(elem)+
    stat_boxplot(geom='errorbar', width=0.25)+
    geom_boxplot(outlier.shape = NA, notch=TRUE)+
    scale_x_continuous(position = "top") +
    scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
    coord_cartesian(xlim=c(as.numeric(x_min),as.numeric(x_max)))+
    xlab(x_label) + ylab("")+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.x=element_line(color="black"), axis.ticks.y = element_blank(), axis.text.x=element_text(size=15), axis.text.y = element_blank())
}


elem_boxplot_func(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV WTKO / Total WTKO)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV KO / Total KO)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV WTtg / Total WTtg)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV TG / Total TG)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV WTKO / Total WTKO)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV KO / Total KO)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV WTtg / Total WTtg)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_boxplot_func(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV TG / Total TG)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_boxplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')


with(mPv_mL4_snmcseq_gene_flank_cCRE_mCA , vioplot( 
  val[element=="Regional"] , val[element=="Gene body"], val[element=="cCRE"],  
  names=c("Regional","Gene body","cCRE"),
  ylab="(mPv mCA/CA)/(mL4 mCA/CA)"
))

median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75), na.rm=TRUE)
  names(out) <- c("ymin","y","ymax")
  return(out) 
}
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
elem_viol_func <- function(data_table, elem, title, x_label, x_min=-5, x_max=5){
  data = data_table[element==elem]
  IQR = quantile(data[, val], na.rm=TRUE)[4] - quantile(data[, val], na.rm=TRUE)[2]
  outlier_min = quantile(data[, val], na.rm=TRUE)[2] - 1.5*IQR
  outlier_max = quantile(data[, val], na.rm=TRUE)[4] + 1.5*IQR
  data = data[(val >= outlier_min) & (val <= outlier_max)]
  pval <- wilcox.test(data[(enriched_cell=="PV"), as.numeric(val)], data[(enriched_cell=="Exc"), as.numeric(val)])$p.value
  sig_pval <- sig_function(pval)
  print(c(pval, sig_pval))
  p <- ggplot(data, aes(y = enriched_cell, x = as.numeric(val), fill=enriched_cell))+ 
    scale_x_continuous(position = "top")+ 
    ggtitle(elem)+
    geom_violin()+
    stat_summary(fun=median, geom="point")+
    scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
    coord_cartesian(xlim=c(as.numeric(x_min),as.numeric(x_max)))+
    xlab(x_label) + ylab("")+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.x=element_line(color="black"), axis.ticks.y = element_blank(), axis.text.x=element_text(size=15), axis.text.y = element_blank())
  return(p)
}


elem_viol_func(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV WTKO / Total WTKO)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV KO / Total KO)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV WTtg / Total WTtg)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt, elem="Linked cCREs", x_label="Log2 H3K27ac FD (PV TG / Total TG)", x_min=-4, x_max=4)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_linkedcCREs_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV WTKO / Total WTKO)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTKO_TotalMeCP2WTKO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV KO / Total KO)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2KO_TotalMeCP2KO_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')


elem_viol_func(enrich_PV_WTtg_Total_WTtg_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV WTtg / Total WTtg)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2WTtg_TotalMeCP2WTtg_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

elem_viol_func(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt, elem="Promoters", x_label="Log2 H3K27ac FD (PV TG / Total TG)", x_min=-5, x_max=5)
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/chip_plots/Pv_H3K27ac/cell_confusion_paper_chip_plots/Log2_H3K27ac_PVMeCP2TG_TotalMeCP2TG_foldDiff_Huntley2020_Exc_Pv_enriched_promoters_vioplot.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')



with(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt , vioplot( 
  val[enriched_cell=="Exc"] , val[enriched_cell=="PV"],
  xlab="",
  names=c("Exc","PV"),
  col=c("gold", "forestgreen"),
  horizontal = TRUE,
))
xaxt = "n"
axis(3)

with(enrich_PV_KO_Total_KO_prom_cCRE_H3K27ac_dt , vioplot( 
  val[enriched_cell=="Exc"] , val[enriched_cell=="PV"],
  xlab="",
  names=c("Exc","PV"),
  col=c("gold", "forestgreen"),
  xaxt="n",
  yaxt="n",
  horizontal=TRUE
))
axis(3)
mtext("Log2 H3K27ac FD (PV KO / Total KO)", side=3, line=3)

with(enrich_PV_Total_WTKO_prom_cCRE_H3K27ac_dt , vioplot( 
  val[enriched_cell=="Exc"] , val[enriched_cell=="PV"],
  xlab="",
  names=c("Exc","PV"),
  col=c("gold", "forestgreen"),
  xaxt="n",
  yaxt="n",
  horizontal=TRUE
))
axis(3)
mtext("Log2 H3K27ac FD (PV WTKO / Total WTKO)", side=3, line=3)

ggplot(enrich_PV_TG_Total_TG_prom_cCRE_H3K27ac_dt[element=="Promoters", ], aes(y = enriched_cell, x = as.numeric(val), fill=enriched_cell))+ 
  ggtitle("Promoters")+
  geom_violin()+
  geom_boxplot(width=0.1, outlier.shape=NA)+
  scale_x_continuous(position="top")+
  stat_summary(fun=median, geom="point", color="black")+
  scale_fill_manual(name = "Genes enriched in:", values = c("Exc"="gold", "PV"="forestgreen")) +
  coord_cartesian(xlim=c(-6,6))+
  xlab("Log2 H3K27ac FD (PV TG / Total TG)") + ylab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x=element_line(color="black"), axis.ticks.y = element_blank(), axis.text.x=element_text(size=15), axis.text.y = element_blank())
