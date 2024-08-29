library(data.table)
library(ggplot2)

chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
all_genes_mm9 =fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_flat.txt")
coding_genes_mm9 = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')
#should have same gene order
Pv_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")
Sst_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_Mecp2KO_gene_TPMs_nondedup.txt")
L4_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_Mecp2KO_gene_TPMs_nondedup.txt")
L5_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_Mecp2KO_gene_TPMs_nondedup.txt")

#
#gene body and flank methylation tables, using genic methylation from TSS+3kb to TES, should have the same gene order
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

#adding a gene body methylation column with the subclass's name for downstream analysis
PV_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol <- copy(PV_INTACT_gene_body_TSSplus3kb_flank_mCA)
PV_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, PV_gene_methylation_corrected := gene_methylation_corrected]

SST_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol <- copy(SST_INTACT_gene_body_TSSplus3kb_flank_mCA)
SST_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, SST_gene_methylation_corrected := gene_methylation_corrected]

L4_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol <- copy(L4_INTACT_gene_body_TSSplus3kb_flank_mCA)
L4_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, L4_gene_methylation_corrected := gene_methylation_corrected]

L5_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol <- copy(L5_INTACT_gene_body_TSSplus3kb_flank_mCA)
L5_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, L5_gene_methylation_corrected := gene_methylation_corrected]


names(coding_genes_mm9) = c("chrom", "start", "end", "Gene", "transcripts", "strand")
#coding genes with chromosomes in specific chromosome list
coding_genes_mm9_chr1toX = coding_genes_mm9[chrom %in% chrom_list]
#add subclass gene body methylation data to list of gene names, using genes at least 5kb long
coding_genes_mm9_TSSplus3kb_mCA_subclasses = data.table(inner_join(x=coding_genes_mm9_chr1toX[end-start >= 5000, .(Gene)], y=PV_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, .(gene, PV_gene_methylation_corrected)], by=c("Gene"="gene")))
coding_genes_mm9_TSSplus3kb_mCA_subclasses = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses, y=SST_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, .(gene, SST_gene_methylation_corrected)], by=c("Gene"="gene")))
coding_genes_mm9_TSSplus3kb_mCA_subclasses = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses, y=L4_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, .(gene, L4_gene_methylation_corrected)], by=c("Gene"="gene")))
coding_genes_mm9_TSSplus3kb_mCA_subclasses = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses, y=L5_INTACT_gene_body_TSSplus3kb_flank_mCA_newCol[, .(gene, L5_gene_methylation_corrected)], by=c("Gene"="gene")))


coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses, y=Pv_TPM[, .(Gene, Pv_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs, y=Sst_TPM[, .(Gene, Sst_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs, y=L4_TPM[, .(Gene, L4_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs, y=L5_TPM[, .(Gene, L5_WT_TPM_avg)], by="Gene"))


Pv_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/Sst_WT_TPM_avg > 5), Gene]
Pv_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]
Pv_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

Sst_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/Pv_WT_TPM_avg > 5), Gene]
Sst_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]
Sst_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

L4_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/Pv_WT_TPM_avg > 5), Gene]
L4_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/Sst_WT_TPM_avg > 5), Gene]
L4_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

L5_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/Pv_WT_TPM_avg> 5), Gene]
L5_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/Sst_WT_TPM_avg> 5), Gene]
L5_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]

xlim=c(0,0.14), ylim=c(0,0.14)

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_SST_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="SST INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_Sst_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_Sst_enriched_genes, SST_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="SST INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_Pv_enriched_genes, SST_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="SST INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.934


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L4_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Pv_enriched_genes, L4_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.793


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_L5_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="L5 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L5_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L5_enriched_genes, L5_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Pv_enriched_genes, L5_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.836

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L4_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Sst_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Sst_enriched_genes, L4_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.831


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_vs_L5_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="L5 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L5_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L5_enriched_genes, L5_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Sst_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Sst_enriched_genes, L5_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.868

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_L4_enriched_genes, L5_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="L5 INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_L5_enriched_genes, L5_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_L5_enriched_genes, L4_gene_methylation_corrected], xlab="L5 INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.935


###with limits
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_SST_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="SST INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_Sst_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_Sst_enriched_genes, SST_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="SST INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_Pv_enriched_genes, SST_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="SST INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.934


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
par(new=TRUE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L4_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Pv_enriched_genes, L4_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.793


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_vs_L5_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="L5 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L5_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Pv_over_L5_enriched_genes, L5_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Pv_enriched_genes, PV_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Pv_enriched_genes, L5_gene_methylation_corrected], xlab="PV INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.836

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L4_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Sst_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_Sst_enriched_genes, L4_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.831


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_vs_L5_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="L5 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L5_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% Sst_over_L5_enriched_genes, L5_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="darkorange", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Sst_enriched_genes, SST_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_Sst_enriched_genes, L5_gene_methylation_corrected], xlab="SST INTACT mCA/CA", ylab="L5 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) #0.868

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_vs_L4_WT_KO_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCA/mCA", ylab="L4 INTACT mCA/CA",  pch=16, cex=0.6, col="lightgray", asp=1, xlim=c(0,0.14), ylim=c(0,0.14),axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.14, by = 0.02))
axis(2, at = seq(0, 0.14, by = 0.02))
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_L4_enriched_genes, L5_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L5_over_L4_enriched_genes, L4_gene_methylation_corrected], xlab="L5 INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="orchid", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_L5_enriched_genes, L5_gene_methylation_corrected], y=coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[Gene %in% L4_over_L5_enriched_genes, L4_gene_methylation_corrected], xlab="L5 INTACT mCA/CA", ylab="L4 INTACT mCA/CA", xlim=c(0, 0.14), ylim=c(0, 0.14), col="skyblue", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected]~coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected]))
dev.off()

round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="pearson", use="complete.obs"), 3) # 0.935

#spearman corrs
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #0.952
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #0.856
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, PV_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #0.881
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #0.881
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, SST_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #0.907
round(cor(coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L5_gene_methylation_corrected], coding_genes_mm9_TSSplus3kb_mCA_subclasses_TPMs[, L4_gene_methylation_corrected], method="spearman", use="complete.obs"), 3) # 0.95


