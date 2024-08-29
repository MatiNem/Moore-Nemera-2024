library(data.table)
library(dplyr)

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")


chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

nonPromoter_cCREs_mm9 = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_using_ensgene_mm9_1kb_promoterWindows.bed")
mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
mousebrain_union_nonPromoter_cCREs_Cicero_intragenicLinked_genes_mm9 = mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene==1)]

PV_TAD_mCA = fread("HG_lab/Mati/GabelLab/hic_data/cortex_domains/HiC_ncx_CN_ArrowheadContacts_5kb_KR/Bonev_cortex_ArrowheadContacts_5kb_KR_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
SST_TAD_mCA = fread("HG_lab/Mati/GabelLab/hic_data/cortex_domains/HiC_ncx_CN_ArrowheadContacts_5kb_KR/Bonev_cortex_ArrowheadContacts_5kb_KR_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
L4_TAD_mCA = fread("HG_lab/Mati/GabelLab/hic_data/cortex_domains/HiC_ncx_CN_ArrowheadContacts_5kb_KR/Bonev_cortex_ArrowheadContacts_5kb_KR_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
L5_TAD_mCA = fread("HG_lab/Mati/GabelLab/hic_data/cortex_domains/HiC_ncx_CN_ArrowheadContacts_5kb_KR/Bonev_cortex_ArrowheadContacts_5kb_KR_L5_WT_KO_deep_INTACT_mCA_mm9.bed")

names(PV_TAD_mCA)= c("TAD_chrom", "TAD_start", "TAD_end", "TAD_coords", "TAD_meth", "TAD_cov")
names(SST_TAD_mCA)= c("TAD_chrom", "TAD_start", "TAD_end", "TAD_coords", "TAD_meth", "TAD_cov")
names(L4_TAD_mCA)= c("TAD_chrom", "TAD_start", "TAD_end", "TAD_coords", "TAD_meth", "TAD_cov")
names(L5_TAD_mCA)= c("TAD_chrom", "TAD_start", "TAD_end", "TAD_coords", "TAD_meth", "TAD_cov")

PV_TAD_mCA[, TAD_ID := paste0(TAD_chrom,":",TAD_start,"-",TAD_end)]
SST_TAD_mCA[, TAD_ID := paste0(TAD_chrom,":",TAD_start,"-",TAD_end)]
L4_TAD_mCA[, TAD_ID := paste0(TAD_chrom,":",TAD_start,"-",TAD_end)]
L5_TAD_mCA[, TAD_ID := paste0(TAD_chrom,":",TAD_start,"-",TAD_end)]

union_cCREs_TADs_overlapRegions_only = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR_overlapRegionsOnly.bed")
union_cCREs_TADs_fullCoords = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR.bed")                                             

 

union_cCREs_TADs_overlapRegions_only_PV_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR_overlapRegionsOnly_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
union_cCREs_TADs_overlapRegions_only_SST_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR_overlapRegionsOnly_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
union_cCREs_TADs_overlapRegions_only_L4_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR_overlapRegionsOnly_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
union_cCREs_TADs_overlapRegions_only_L5_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_mm9_Bonev_cortexTADs_Arrowhead5kb_KR_overlapRegionsOnly_L5_WT_KO_deep_INTACT_mCA_mm9.bed")

names(union_cCREs_TADs_overlapRegions_only_PV_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "cCRE_meth", "cCRE_cov")
names(union_cCREs_TADs_overlapRegions_only_SST_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "cCRE_meth", "cCRE_cov")
names(union_cCREs_TADs_overlapRegions_only_L4_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "cCRE_meth", "cCRE_cov")
names(union_cCREs_TADs_overlapRegions_only_L5_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "cCRE_meth", "cCRE_cov")


#union_cCREs_TADs_overlapRegions_only_TADcoords = cbind(union_cCREs_TADs_overlapRegions_only, union_cCREs_TADs_fullCoords)

union_cCREs_TADs_overlapRegions_only_PV_mCA = cbind(union_cCREs_TADs_overlapRegions_only_PV_mCA, TAD_ID=union_cCREs_TADs_fullCoords[, V5])
union_cCREs_TADs_overlapRegions_only_PV_mCA = unique(data.table(inner_join(x=union_cCREs_TADs_overlapRegions_only_PV_mCA, y=PV_TAD_mCA[, .(TAD_ID, TAD_meth, TAD_cov)], by=c("TAD_ID"))))
union_cCREs_TADs_overlapRegions_only_PV_mCA[, cCRE_methylation := cCRE_meth/cCRE_cov]
union_cCREs_TADs_overlapRegions_only_PV_mCA[, TAD_methylation := (TAD_meth - cCRE_meth)/(TAD_cov - cCRE_cov)]
union_cCREs_TADs_overlapRegions_only_PV_mCA$cCRE_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_PV_mCA$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
union_cCREs_TADs_overlapRegions_only_PV_mCA$TAD_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_PV_mCA$TAD_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_PV_mCA[TAD_methylation_corrected < 0, TAD_methylation_corrected := 0]

union_cCREs_TADs_overlapRegions_only_SST_mCA = cbind(union_cCREs_TADs_overlapRegions_only_SST_mCA, TAD_ID=union_cCREs_TADs_fullCoords[, V5])
union_cCREs_TADs_overlapRegions_only_SST_mCA = unique(data.table(inner_join(x=union_cCREs_TADs_overlapRegions_only_SST_mCA, y=SST_TAD_mCA[, .(TAD_ID, TAD_meth, TAD_cov)], by=c("TAD_ID"))))
union_cCREs_TADs_overlapRegions_only_SST_mCA[, cCRE_methylation := cCRE_meth/cCRE_cov]
union_cCREs_TADs_overlapRegions_only_SST_mCA[, TAD_methylation := (TAD_meth - cCRE_meth)/(TAD_cov - cCRE_cov)]
union_cCREs_TADs_overlapRegions_only_SST_mCA$cCRE_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_SST_mCA$cCRE_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
union_cCREs_TADs_overlapRegions_only_SST_mCA$TAD_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_SST_mCA$TAD_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_SST_mCA[TAD_methylation_corrected < 0, TAD_methylation_corrected := 0]


union_cCREs_TADs_overlapRegions_only_L4_mCA = cbind(union_cCREs_TADs_overlapRegions_only_L4_mCA, TAD_ID=union_cCREs_TADs_fullCoords[, V5])
union_cCREs_TADs_overlapRegions_only_L4_mCA = unique(data.table(inner_join(x=union_cCREs_TADs_overlapRegions_only_L4_mCA, y=L4_TAD_mCA[, .(TAD_ID, TAD_meth, TAD_cov)], by=c("TAD_ID"))))
union_cCREs_TADs_overlapRegions_only_L4_mCA[, cCRE_methylation := cCRE_meth/cCRE_cov]
union_cCREs_TADs_overlapRegions_only_L4_mCA[, TAD_methylation := (TAD_meth - cCRE_meth)/(TAD_cov - cCRE_cov)]
union_cCREs_TADs_overlapRegions_only_L4_mCA$cCRE_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_L4_mCA$cCRE_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
union_cCREs_TADs_overlapRegions_only_L4_mCA$TAD_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_L4_mCA$TAD_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_L4_mCA[TAD_methylation_corrected < 0, TAD_methylation_corrected := 0]



union_cCREs_TADs_overlapRegions_only_L5_mCA = cbind(union_cCREs_TADs_overlapRegions_only_L5_mCA, TAD_ID=union_cCREs_TADs_fullCoords[, V5])
union_cCREs_TADs_overlapRegions_only_L5_mCA = unique(data.table(inner_join(x=union_cCREs_TADs_overlapRegions_only_L5_mCA, y=L5_TAD_mCA[, .(TAD_ID, TAD_meth, TAD_cov)], by=c("TAD_ID"))))
union_cCREs_TADs_overlapRegions_only_L5_mCA[, cCRE_methylation := cCRE_meth/cCRE_cov]
union_cCREs_TADs_overlapRegions_only_L5_mCA[, TAD_methylation := (TAD_meth - cCRE_meth)/(TAD_cov - cCRE_cov)]
union_cCREs_TADs_overlapRegions_only_L5_mCA$cCRE_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_L5_mCA$cCRE_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
union_cCREs_TADs_overlapRegions_only_L5_mCA$TAD_methylation_corrected <- union_cCREs_TADs_overlapRegions_only_L5_mCA$TAD_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
union_cCREs_TADs_overlapRegions_only_L5_mCA[TAD_methylation_corrected < 0, TAD_methylation_corrected := 0]


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_PV_WT_KO_deep_INTACT_mCAperCA_scatterplot.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="PV TAD mCA/CA", ylab="PV cCRE mCA/CA", main="", asp=1, xlim=c(0,0.2), ylim=c(0, 0.2))
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_SST_WT_KO_deep_INTACT_mCAperCA_scatterplot.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="SST TAD mCA/CA", ylab="SST cCRE mCA/CA", main="", asp=1, xlim=c(0,0.2), ylim=c(0, 0.2))
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_L4_WT_KO_deep_INTACT_mCAperCA_scatterplot.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="L4 TAD mCA/CA", ylab="L4 cCRE mCA/CA", main="", asp=1, xlim=c(0,0.2), ylim=c(0, 0.2))
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_L5_WT_KO_deep_INTACT_mCAperCA_scatterplot.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="L5 TAD mCA/CA", ylab="L5 cCRE mCA/CA", main="", asp=1, xlim=c(0,0.2), ylim=c(0, 0.2))
dev.off()



#correlations
round(cor(union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.495
round(cor(union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.445
round(cor(union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.436
round(cor(union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.437


###0 to 0.15
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_PV_WT_KO_deep_INTACT_mCAperCA_scatterplot_0to0.15.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_PV_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="PV TAD mCA/CA", ylab="PV cCRE mCA/CA", main="", asp=1, xlim=c(0,0.15), ylim=c(0, 0.15))
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_SST_WT_KO_deep_INTACT_mCAperCA_scatterplot_0to0.15.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_SST_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="SST TAD mCA/CA", ylab="SST cCRE mCA/CA", main="", asp=1, xlim=c(0,0.15), ylim=c(0, 0.15))
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_L4_WT_KO_deep_INTACT_mCAperCA_scatterplot_0to0.15.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L4_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="L4 TAD mCA/CA", ylab="L4 cCRE mCA/CA", main="", asp=1, xlim=c(0,0.15), ylim=c(0, 0.15))
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/mousebrain_union_nonPromoter_cCRE_vs_TAD_L5_WT_KO_deep_INTACT_mCAperCA_scatterplot_0to0.15.png", width=2000, height=2000, res=300)
par(bg="blue")
smoothScatter(union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], TAD_methylation_corrected], union_cCREs_TADs_overlapRegions_only_L5_mCA[cCRE_label %in% nonPromoter_cCREs_mm9[,V4], cCRE_methylation_corrected], colramp = colorRampPalette(c("blue","yellow","red")), transformation = function(x) x^.5, xlab="L5 TAD mCA/CA", ylab="L5 cCRE mCA/CA", main="", asp=1, xlim=c(0,0.15), ylim=c(0, 0.15))
dev.off()
