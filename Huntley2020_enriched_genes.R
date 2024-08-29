library(data.table)
library(dplyr)
#coding genes
coding_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")
#Load data table containing expression of genes in excitatory, Pv, Sst, and Vip cell types
cellType_gene_exp = fread("HG_lab/Mati/GabelLab/RNAseq/Huntley2020_mouse/Huntley2020_EXC_PV_SST_VIP_expr.txt")
#genes upregulated in excitatory neurons relative to Pv neurons
Exc_over_Pv_genes = cellType_gene_exp[(adj.p.EXC.vs.PV <= 0.05) & (lfc.EXC.vs.PV > 0) & (rpkm.EXC >= 10), SYMBOL]
#genes upregulated in Pv neurons relative to excitatory neurons
Pv_over_Exc_genes = cellType_gene_exp[(adj.p.EXC.vs.PV <= 0.05) & (lfc.EXC.vs.PV < 0) & (rpkm.PV >= 10), SYMBOL]

#write out genes enriched in excitatory neurons over Pv neurons
write.table(coding_genes_mm9[V4 %in% Exc_over_Pv_genes], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
#write out genes enriched in Pv neurons over excitatory over Pv neurons
write.table(coding_genes_mm9[V4 %in% Exc_over_Pv_genes], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")


promoterWindows = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_1kb_promoterWindows.bed")
promoterWindows_coding = promoterWindows[V4 %in% coding_genes_mm9[,V4]]

Exc_over_Pv_genes_promoterWindows = promoterWindows_coding[V4 %in% Exc_over_Pv_genes,]
Pv_over_Exc_genes_promoterWindows = promoterWindows_coding[V4 %in% Pv_over_Exc_genes,]

write.table(Exc_over_Pv_genes_promoterWindows, file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_promoterWindows_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Pv_over_Exc_genes_promoterWindows, file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_promoterWindows_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")


#top 200 differential genes
Exc_over_Pv_genes_top200 = cellType_gene_exp[!is.na(SYMBOL) & (adj.p.EXC.vs.PV <= 0.05) & (lfc.EXC.vs.PV > 0) & (rpkm.EXC >= 10) & (SYMBOL %in% coding_genes_mm9[,V4])][order(-lfc.EXC.vs.PV), SYMBOL][1:200]
Pv_over_Exc_genes_top200 = cellType_gene_exp[!is.na(SYMBOL) & (adj.p.EXC.vs.PV <= 0.05) & (lfc.EXC.vs.PV < 0) & (rpkm.PV >= 10) & (SYMBOL %in% coding_genes_mm9[,V4])][order(lfc.EXC.vs.PV), SYMBOL][1:200]
write.table(coding_genes_mm9[V4 %in% Exc_over_Pv_genes_top200], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top200_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(coding_genes_mm9[V4 %in% Pv_over_Exc_genes_top200], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top200_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")


Exc_over_Pv_genes_top200_promoterWindows = promoterWindows_coding[V4 %in% Exc_over_Pv_genes_top200,]
Pv_over_Exc_genes_top200_promoterWindows = promoterWindows_coding[V4 %in% Pv_over_Exc_genes_top200,]

write.table(Exc_over_Pv_genes_top200_promoterWindows, file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top200_promoterWindows_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Pv_over_Exc_genes_top200_promoterWindows, file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top200_promoterWindows_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")

#non-promoter cCREs linked to genes
mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")


Exc_over_Pv_genes_top200_linked_union_nonPromoter_cCREcoords_mm9 = mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Exc_over_Pv_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]
Pv_over_Exc_genes_top200_linked_union_nonPromoter_cCREcoords_mm9 = mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Pv_over_Exc_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]

Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9 = PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Exc_over_Pv_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]
Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9 = PVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Pv_over_Exc_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]

Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Exc_over_Pv_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]
Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9 = nonPVGA_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Pv_over_Exc_genes_top200, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]


write.table(Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Exc_over_Pv_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Pv_over_Exc_genes_top200_linked_PVGA_nonPromoter_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Exc_over_Pv_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Huntley2020_genes/Pv_over_Exc_genes_top200_linked_nonPVGA_nonPromoter_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")

#more stringent selection of genes
#top 200 differential genes
Exc_over_Pv_genes_top100_q0.01 = cellType_gene_exp[!is.na(SYMBOL) & (adj.p.EXC.vs.PV <= 0.01) & (lfc.EXC.vs.PV > 0) & (rpkm.EXC >= 10) & (SYMBOL %in% coding_genes_mm9[,V4])][order(-lfc.EXC.vs.PV), SYMBOL][1:100]
Pv_over_Exc_genes_top100_q0.01 = cellType_gene_exp[!is.na(SYMBOL) & (adj.p.EXC.vs.PV <= 0.01) & (lfc.EXC.vs.PV < 0) & (rpkm.PV >= 10) & (SYMBOL %in% coding_genes_mm9[,V4])][order(lfc.EXC.vs.PV), SYMBOL][1:100]
write.table(coding_genes_mm9[V4 %in% Exc_over_Pv_genes_top100_q0.01], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Exc_over_Pv_genes_top100_q0.01_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(coding_genes_mm9[V4 %in% Pv_over_Exc_genes_top100_q0.01], file="HG_lab/Mati/GabelLab/genesets/Huntley2020_mouse/Huntley2020_Pv_over_Exc_genes_top100_q0.01_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")


Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9 = mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Exc_over_Pv_genes_top100_q0.01, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]
Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9 = mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_genicBooleans[Gene %in% Pv_over_Exc_genes_top100_q0.01, .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)]
write.table(Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Exc_over_Pv_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_over_Exc_genes_top100_q0.01_linked_union_nonPromoter_cCREcoords_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

