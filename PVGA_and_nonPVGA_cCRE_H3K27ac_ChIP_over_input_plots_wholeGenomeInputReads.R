library(data.table)
library(dplyr)
library(ggplot2)


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

#total read count tables, PV MeCP2 KO and MeCP2 TG H3K27ac ChIP-seq experiments
PV_KO_H3K27ac_ChIP_total_read_numbers <- fread("HG_lab/Mati/GabelLab/ChIPseq/Pv_H3K27ac/PV_KO//PV_KO_H3K27ac_ChIP_total_read_numbers.csv")
PV_TG_H3K27ac_ChIP_total_read_numbers <- fread("HG_lab/Mati/GabelLab/ChIPseq/Pv_H3K27ac/PV_TG//PV_TG_H3K27ac_ChIP_total_read_numbers.csv")



cpm_function <- function(data_table, read_column){
  data_table[, cpm := (10**6) * data_table[, get(read_column)]/sum(data_table[, get(read_column)])]
}

#Read in files with Pv MeCP2 KO H3K27ac ChIP counts in cCREs
cCREs_1500bp_H3K27ac_PV_WTKO_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTko_B1_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTKO_B2_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTKO_B3_CTTCCCAT_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B1_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B2_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_KO_B3_TCGCGAGG_chipCounts_mm9.bed")


#cCREs_1500bp_H3K27ac_PV_WTKO_B1[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B1[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B1[,V5]))]
#cCREs_1500bp_H3K27ac_PV_WTKO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B2[,V5]))]
#cCREs_1500bp_H3K27ac_PV_WTKO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_WTKO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_PV_WTKO_B3[,V5]))]
#cCREs_1500bp_H3K27ac_PV_KO_B1[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B1[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B1[,V5]))]
#cCREs_1500bp_H3K27ac_PV_KO_B2[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B2[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B2[,V5]))]
#cCREs_1500bp_H3K27ac_PV_KO_B3[, cpm := (10**6) * (cCREs_1500bp_H3K27ac_PV_KO_B3[,V5]/sum(cCREs_1500bp_H3K27ac_PV_KO_B3[,V5]))]

#cpm calculation using total number of reads mapped to the genome as denominator
cCREs_1500bp_H3K27ac_PV_WTKO_B1[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WTko_B1_AAGGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_WTKO_B2[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WTKO_B2_AAGGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_WTKO_B3[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WTKO_B3_CTTCCCAT.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_KO_B1[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_KO_B1_GCTGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_KO_B2[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_KO_B2_GCTGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_KO_B3[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_KO_B3_TCGCGAGG.U.bam", read_count])]


#cpm_function(cCREs_1500bp_H3K27ac_PV_WTKO_B1, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_WTKO_B2, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_WTKO_B3, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_KO_B1, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_KO_B2, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_KO_B3, "V5")


cCREs_1500bp_Input_PV_WTKO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_WTKO_B2_CGACCTAA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_WTKO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_WTKO_B3_CTCGAACA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_KO_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_KO_B2_TACATCGG_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_KO_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_KO/mousebrain_union_cCREs_Input_PV_KO_B3_TCTAGGAG_chipCounts_mm9.bed")


cCREs_1500bp_Input_PV_WTKO_B2[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_WTKO_B2_CGACCTAA.U.bam", read_count])]
cCREs_1500bp_Input_PV_WTKO_B3[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_WTKO_B3_CTCGAACA.U.bam", read_count])]
cCREs_1500bp_Input_PV_KO_B2[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_KO_B2_TACATCGG.U.bam", read_count])]
cCREs_1500bp_Input_PV_KO_B3[, cpm := (10**6) * (V5/PV_KO_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_KO_B3_TCTAGGAG.U.bam", read_count])]

#cpm_function(cCREs_1500bp_Input_PV_WTKO_B2, "V5")
#cpm_function(cCREs_1500bp_Input_PV_WTKO_B3, "V5")
#cpm_function(cCREs_1500bp_Input_PV_KO_B2, "V5")
#cpm_function(cCREs_1500bp_Input_PV_KO_B3, "V5")


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


#H3K27ac signal in 1500bp windows centered on cCRE centers for PV MeCP2 OE (TG) and its wild-type littermates
cCREs_1500bp_H3K27ac_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WT_B1_AAGGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTtg_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTtg_B2_CGTCTAAC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_WTtg_B3_TCGCGAGG_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B1_GCTGAGTC_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B2 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B2_GACTACGA_chipCounts_mm9.bed")
cCREs_1500bp_H3K27ac_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_H3K27ac_PV_Tg_B3_TGCATGTA_chipCounts_mm9.bed")


#cpm calculation using total number of reads mapped to the genome as denominator
cCREs_1500bp_H3K27ac_PV_WTtg_B1[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WT_B1_AAGGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_WTtg_B2[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WTtg_B2_CGTCTAAC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_WTtg_B3[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_WTtg_B3_TCGCGAGG.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_TG_B1[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_Tg_B1_GCTGAGTC.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_TG_B2[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_Tg_B2_GACTACGA.U.bam", read_count])]
cCREs_1500bp_H3K27ac_PV_TG_B3[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="H3K27ac_PV_Tg_B3_TGCATGTA.U.bam", read_count])]


#cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B1, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B2, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_WTtg_B3, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B1, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B2, "V5")
#cpm_function(cCREs_1500bp_H3K27ac_PV_TG_B3, "V5")


cCREs_1500bp_Input_PV_WTtg_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_WT_B1_TCGCGAGG_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_WTtg_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_WTtg_B3_AACTCGGA_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_TG_B1 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_Tg_B1_CTTCCCAT_chipCounts_mm9.bed")
cCREs_1500bp_Input_PV_TG_B3 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chipCounts/PV_TG/mousebrain_union_cCREs_1500bpWindows_Input_PV_Tg_B3_ACTCCTAC_chipCounts_mm9.bed")


cCREs_1500bp_Input_PV_WTtg_B1[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_WT_B1_TCGCGAGG.U.bam", read_count])]
cCREs_1500bp_Input_PV_WTtg_B3[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_WTtg_B3_AACTCGGA.U.bam", read_count])]
cCREs_1500bp_Input_PV_TG_B1[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_Tg_B1_CTTCCCAT.U.bam", read_count])]
cCREs_1500bp_Input_PV_TG_B3[, cpm := (10**6) * (V5/PV_TG_H3K27ac_ChIP_total_read_numbers[filename=="Input_PV_Tg_B3_ACTCCTAC.U.bam", read_count])]


#cpm_function(cCREs_1500bp_Input_PV_WTtg_B1, "V5")
#cpm_function(cCREs_1500bp_Input_PV_WTtg_B3, "V5")
#cpm_function(cCREs_1500bp_Input_PV_TG_B1, "V5")
#cpm_function(cCREs_1500bp_Input_PV_TG_B3, "V5")

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

##
cCREs_1500bp_H3K27ac_Input_PV_WTKO
cCREs_1500bp_H3K27ac_Input_PV_KO
cCREs_1500bp_H3K27ac_Input_PV_WTtg
cCREs_1500bp_H3K27ac_Input_PV_TG


#cCRE lists
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



#intragenic
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

##extragenic cCREs linked to genes
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


###PV WTKO H3K27ac signal in cCREs linked to PV MR and unchanged genes
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Extragenic"),
                     cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTKO", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Extragenic")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_intra_and_extra_linked_to_PV_MR_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("WT littermate of PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-2.4,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  facet_grid(.~genic_loc,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

#wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=1.106535e-10, ****
#wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Intragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=1.132884e-12, ****
#wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n unchanged genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n unchanged genes", as.numeric(logFC)])$p.value #p=6.500396e-12, ****
#wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic Pv cCREs,\n MR genes", as.numeric(logFC)], PVGA_and_nonPVGA_nonPromoter_cCREs_1500bp_Pv_MeCP2KO_H3K27ac_ChIP_MR_genes_nondedup_intra_and_extra[cCRE_group=="Extragenic non-Pv cCREs,\n MR genes", as.numeric(logFC)])$p.value #p=0.01834795, *


##PV KO
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="KO", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Extragenic")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_intra_and_extra_linked_to_PV_MR_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-2.4,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  facet_grid(.~genic_loc,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2_H3K27ac_CPM_ChIP_over_input_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2_H3K27ac_CPM_ChIP_over_input_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


###
###PV WTtg
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Extragenic"),
                                                                                                      cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="WTtg", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Extragenic")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_intra_and_extra_linked_to_PV_MR_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("WT littermate of PV MeCP2 OE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-2.4,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  facet_grid(.~genic_loc,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2WTtg_log2_H3K27ac_CPM_ChIP_over_input_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2WTtg_log2_H3K27ac_CPM_ChIP_over_input_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

###PV TG
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_intragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Intragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Non-Pv cCREs,\n unchanged genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_PVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Pv cCREs,\n MR genes", genic_loc="Extragenic"),
                                                                                                        cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% Pv_MR_genes_q0.1_coding_prefilt5_nondedup_extragenicLinked_nonPVGA_nonPromoter_cCREcoords_mm9[, V4])], condition="TG", cCRE_group="Non-Pv cCREs,\n MR genes", genic_loc="Extragenic")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs,\n unchanged genes", "Non-Pv cCREs,\n unchanged genes", "Pv cCREs,\n MR genes", "Non-Pv cCREs,\n MR genes")))
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes %>% mutate(genic_loc = factor(genic_loc, levels=c("Intragenic", "Extragenic")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_intra_and_extra_linked_to_PV_MR_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("PV MeCP2 OE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs,\n unchanged genes"="gray", "Non-Pv cCREs,\n unchanged genes"="gray", "Pv cCREs,\n MR genes"="red", "Non-Pv cCREs,\n MR genes"="red")) +
  coord_cartesian(ylim=c(-2.4,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  facet_grid(.~genic_loc,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2_H3K27ac_CPM_ChIP_over_input_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_intra_and_extra_linked_to_MR_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2_H3K27ac_CPM_ChIP_over_input_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


##PVGA and nonPVGA nonPromoter cCREs
PVGA_nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/PVGA/PVGA_nonPromoter_cCREs_mm9.bed")
nonPVGA_nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/nonPVGA_nonPromoter_cCREs_mm9.bed")

names(PVGA_nonPromoter_cCREs_mm9) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")
names(nonPVGA_nonPromoter_cCREs_mm9) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label")

#linked cCREs
PVGA_nonPromoter_cCREs_linked_genes=fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_linked_genes=fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
#split by intragenic and extragenic
PVGA_nonPromoter_cCREs_intragenicLinked_genes = PVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==1)]
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes = nonPVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==1)]

PVGA_nonPromoter_cCREs_extragenicLinked_genes = PVGA_nonPromoter_cCREs_linked_genes[(Intragenic==0) & (Intragenic_to_linked_gene==0)]
nonPVGA_nonPromoter_cCREs_extragenicLinked_genes = nonPVGA_nonPromoter_cCREs_linked_genes[(Intragenic==0) & (Intragenic_to_linked_gene==0)]



###don't split up intragenic and extragenic linked cCREs

###PV WTKO H3K27ac signal in PVGA and nonPVGA cCREs 
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% PVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="WTKO", cCRE_group="Pv cCREs"),
                                                                                  cbind(cCREs_1500bp_H3K27ac_Input_PV_WTKO[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="WTKO", cCRE_group="Non-Pv cCREs")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("WT littermate of PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("Pv cCREs"="red", "Non-Pv cCREs"="red")) +
  coord_cartesian(ylim=c(-2,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("WT littermate of PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("Pv cCREs"="aquamarine4", "Non-Pv cCREs"="aquamarine4")) +
  coord_cartesian(ylim=c(-2,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot_newColor.png", width = 2.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTKO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot_newColor.eps", width = 2.5, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes[(cCRE_group=="Pv cCREs"), log2_cpm_ratio_pseud], nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTKO_linked_to_genes[(cCRE_group=="Non-Pv cCREs"), log2_cpm_ratio_pseud])$p.value #p < 2.2e-16, ****

#PV KO H3K27ac signal in PVGA and nonPVGA cCREs 
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% PVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="KO", cCRE_group="Pv cCREs"),
                                                                                cbind(cCREs_1500bp_H3K27ac_Input_PV_KO[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="KO", cCRE_group="Non-Pv cCREs")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "", values = c("Pv cCREs"="red", "Non-Pv cCREs"="red")) +
  coord_cartesian(ylim=c(-2,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2KO_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes[(cCRE_group=="Pv cCREs"), log2_cpm_ratio_pseud], nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_KO_linked_to_genes[(cCRE_group=="Non-Pv cCREs"), log2_cpm_ratio_pseud])$p.value #p < 2.2e-16, ****


###PV WTtg H3K27ac signal in PVGA and nonPVGA cCREs 
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% PVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="WTtg", cCRE_group="Pv cCREs"),
                                                                                  cbind(cCREs_1500bp_H3K27ac_Input_PV_WTtg[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="WTtg", cCRE_group="Non-Pv cCREs")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("WT littermate of PV MeCP2 OE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs"="red", "Non-Pv cCREs"="red")) +
  coord_cartesian(ylim=c(-2,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTtg_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2WTtg_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes[(cCRE_group=="Pv cCREs"), log2_cpm_ratio_pseud], nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_WTtg_linked_to_genes[(cCRE_group=="Non-Pv cCREs"), log2_cpm_ratio_pseud])$p.value #p < 2.2e-16, ****

###PV TG H3K27ac signal in PVGA and nonPVGA cCREs 
nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes = data.table(rbind(cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% PVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="TG", cCRE_group="Pv cCREs"),
                                                                                cbind(cCREs_1500bp_H3K27ac_Input_PV_TG[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_linked_genes[, cCRE_label])], condition="TG", cCRE_group="Non-Pv cCREs")))

nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes = nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes %>% mutate(cCRE_group = factor(cCRE_group, levels=c("Pv cCREs", "Non-Pv cCREs")))


ggplot(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes, aes(x = cCRE_group, y = as.numeric(log2_cpm_ratio_pseud), fill=cCRE_group))+
  ggtitle("PV MeCP2 OE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  scale_fill_manual(name = "cCREs linked to:", values = c("Pv cCREs"="red", "Non-Pv cCREs"="red")) +
  coord_cartesian(ylim=c(-2,4))+
  ylab("Log2 H3K27ac CPM ChIP/Input") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=10, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PVGA_and_nonPVGA_nonPromoter_1500bpWindow_cCREs_linked_to_genes_coding_prefilt5_nondedup_PvMeCP2TG_log2_H3K27ac_CPM_ChIP_over_input_wholeGenomeReads_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes[(cCRE_group=="Pv cCREs"), log2_cpm_ratio_pseud], nonPromoter_cCREs_1500bp_H3K27ac_Input_PV_TG_linked_to_genes[(cCRE_group=="Non-Pv cCREs"), log2_cpm_ratio_pseud])$p.value #p < 2.2e-16, ****
