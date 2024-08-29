library(data.table)
library(dplyr)
library(ggplot2)

chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
coding_genes_mm9 = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')
names(coding_genes_mm9) = c("chrom", "start", "end", "Gene", "transcripts", "strand")
#coding genes with chromosomes in specific chromosome list
coding_genes_mm9_chr1toX = coding_genes_mm9[chrom %in% chrom_list]

coding_genes_mm9_chr1toX_minLength5kb <- coding_genes_mm9_chr1toX[end-start >= 500, Gene]


PV_WT_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_WT_deep_INTACT_mCA_mm9.bed")
PV_KO_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_KO_deep_INTACT_mCA_mm9.bed")
L4_WT_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_WT_deep_INTACT_mCA_mm9.bed")
L4_KO_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_KO_deep_INTACT_mCA_mm9.bed")
L5_WT_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_WT_deep_INTACT_mCA_mm9.bed")
L5_KO_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_KO_deep_INTACT_mCA_mm9.bed")
SST_WT_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_WT_deep_INTACT_mCA_mm9.bed")
SST_KO_INTACT_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_KO_deep_INTACT_mCA_mm9.bed")



PV_WT_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_WT_deep_INTACT_mCG_mm9.bed")
PV_KO_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_KO_deep_INTACT_mCG_mm9.bed")
L4_WT_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_WT_deep_INTACT_mCG_mm9.bed")
L4_KO_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_KO_deep_INTACT_mCG_mm9.bed")
L5_WT_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_WT_deep_INTACT_mCG_mm9.bed")
L5_KO_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_KO_deep_INTACT_mCG_mm9.bed")
SST_WT_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_WT_deep_INTACT_mCG_mm9.bed")
SST_KO_INTACT_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_KO_deep_INTACT_mCG_mm9.bed")


PV_WT_INTACT_geneBody_TSSplus3kb_mCA$V9 <- PV_WT_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
PV_KO_INTACT_geneBody_TSSplus3kb_mCA$V9 <- PV_KO_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
L4_WT_INTACT_geneBody_TSSplus3kb_mCA$V9 <- L4_WT_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
L4_KO_INTACT_geneBody_TSSplus3kb_mCA$V9 <- L4_KO_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
L5_WT_INTACT_geneBody_TSSplus3kb_mCA$V9 <- L5_WT_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
L5_KO_INTACT_geneBody_TSSplus3kb_mCA$V9 <- L5_KO_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
SST_WT_INTACT_geneBody_TSSplus3kb_mCA$V9 <- SST_WT_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
SST_KO_INTACT_geneBody_TSSplus3kb_mCA$V9 <- SST_KO_INTACT_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]

PV_WT_INTACT_geneBody_TSSplus3kb_mCG$V9 <- PV_WT_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
PV_KO_INTACT_geneBody_TSSplus3kb_mCG$V9 <- PV_KO_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
L4_WT_INTACT_geneBody_TSSplus3kb_mCG$V9 <- L4_WT_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
L4_KO_INTACT_geneBody_TSSplus3kb_mCG$V9 <- L4_KO_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
L5_WT_INTACT_geneBody_TSSplus3kb_mCG$V9 <- L5_WT_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
L5_KO_INTACT_geneBody_TSSplus3kb_mCG$V9 <- L5_KO_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
SST_WT_INTACT_geneBody_TSSplus3kb_mCG$V9 <- SST_WT_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
SST_KO_INTACT_geneBody_TSSplus3kb_mCG$V9 <- SST_KO_INTACT_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]


#PV, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCA/mCA", ylab="PV KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#L4, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCA/mCA", ylab="L4 KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#L5, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCA/mCA", ylab="L5 KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#SST, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCA/mCA", ylab="SST KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#PV, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCG/mCG", ylab="PV KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#L4, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCG/mCG", ylab="L4 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#L5, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCG/mCG", ylab="L5 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#SST, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCG/mCG", ylab="SST KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
#abline(lm(L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]~L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9]))
legend("topleft", legend = paste0("r=", round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#versions of the above plots without correlation labeled or lines of unity
#PV, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCA/mCA", ylab="PV KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L4, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCA/mCA", ylab="L4 KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L5, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCA/mCA", ylab="L5 KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#SST, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCA/mCA", ylab="SST KO mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()


#PV, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCG/mCG", ylab="PV KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L4, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCG/mCG", ylab="L4 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L5, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCG/mCG", ylab="L5 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#SST, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCG/mCG", ylab="SST KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()


#correlations
round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.951
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.916
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.962
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.958

round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.881
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.892
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.93
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.864

#average nonconversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

##bisulfite nonconversion rate subtraction
meth_calc_func_avg_nonconv_sub <- function(meth_table, label_column, nonconv_table){
  meth_table$gene_methylation_corrected <- meth_table[[9]] - nonconv_table[label==label_column, nonconversion_rate]
  meth_table[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
  return(meth_table)
}

PV_WT_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=PV_WT_INTACT_geneBody_TSSplus3kb_mCA, label_column="PV_WT_lambda", nonconv_table=avg_nonconv)
PV_KO_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=PV_KO_INTACT_geneBody_TSSplus3kb_mCA, label_column="PV_KO_lambda", nonconv_table=avg_nonconv)
L4_WT_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=L4_WT_INTACT_geneBody_TSSplus3kb_mCA, label_column="L4_WT_lambda", nonconv_table=avg_nonconv)
L4_KO_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=L4_KO_INTACT_geneBody_TSSplus3kb_mCA, label_column="L4_KO_lambda", nonconv_table=avg_nonconv)
L5_WT_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=L5_WT_INTACT_geneBody_TSSplus3kb_mCA, label_column="L5_WT_lambda", nonconv_table=avg_nonconv)
L5_KO_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=L5_KO_INTACT_geneBody_TSSplus3kb_mCA, label_column="L5_KO_lambda", nonconv_table=avg_nonconv)
SST_WT_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=SST_WT_INTACT_geneBody_TSSplus3kb_mCA, label_column="SST_WT_lambda", nonconv_table=avg_nonconv)
SST_KO_INTACT_geneBody_TSSplus3kb_mCA <- meth_calc_func_avg_nonconv_sub(meth_table=SST_KO_INTACT_geneBody_TSSplus3kb_mCA, label_column="SST_KO_lambda", nonconv_table=avg_nonconv)

PV_WT_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=PV_WT_INTACT_geneBody_TSSplus3kb_mCG, label_column="PV_WT_lambda", nonconv_table=avg_nonconv)
PV_KO_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=PV_KO_INTACT_geneBody_TSSplus3kb_mCG, label_column="PV_KO_lambda", nonconv_table=avg_nonconv)
L4_WT_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=L4_WT_INTACT_geneBody_TSSplus3kb_mCG, label_column="L4_WT_lambda", nonconv_table=avg_nonconv)
L4_KO_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=L4_KO_INTACT_geneBody_TSSplus3kb_mCG, label_column="L4_KO_lambda", nonconv_table=avg_nonconv)
L5_WT_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=L5_WT_INTACT_geneBody_TSSplus3kb_mCG, label_column="L5_WT_lambda", nonconv_table=avg_nonconv)
L5_KO_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=L5_KO_INTACT_geneBody_TSSplus3kb_mCG, label_column="L5_KO_lambda", nonconv_table=avg_nonconv)
SST_WT_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=SST_WT_INTACT_geneBody_TSSplus3kb_mCG, label_column="SST_WT_lambda", nonconv_table=avg_nonconv)
SST_KO_INTACT_geneBody_TSSplus3kb_mCG <- meth_calc_func_avg_nonconv_sub(meth_table=SST_KO_INTACT_geneBody_TSSplus3kb_mCG, label_column="SST_KO_lambda", nonconv_table=avg_nonconv)



#versions of the above plots using background nonconversion rate subtracted methylation

#PV, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCA/mCA", ylab="PV KO mCA/CA", xlim=c(0,0.12), ylim=c(0,0.12), pch=16, cex=0.6, col="lightgray", asp=1, axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.12, by = 0.02))
axis(2, at = seq(0, 0.12, by = 0.02))
dev.off()

#L4, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCA/mCA", ylab="L4 KO mCA/CA", xlim=c(0,0.12), ylim=c(0,0.12), pch=16, cex=0.6, col="lightgray", asp=1, axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.12, by = 0.02))
axis(2, at = seq(0, 0.12, by = 0.02))
dev.off()

#L5, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCA/mCA", ylab="L5 KO mCA/CA", xlim=c(0,0.12), ylim=c(0,0.12), pch=16, cex=0.6, col="lightgray", asp=1, axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.12, by = 0.02))
axis(2, at = seq(0, 0.12, by = 0.02))
dev.off()

#SST, WT vs KO, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCA/mCA", ylab="SST KO mCA/CA", xlim=c(0,0.12), ylim=c(0,0.12), pch=16, cex=0.6, col="lightgray", asp=1, axes=FALSE)
# Add custom axes
axis(1, at = seq(0, 0.12, by = 0.02))
axis(2, at = seq(0, 0.12, by = 0.02))
dev.off()

##mCG
#PV, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_PV_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV WT mCG/mCG", ylab="PV KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L4, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_L4_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 WT mCG/mCG", ylab="L4 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#L5, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_L5_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 WT mCG/mCG", ylab="L5 KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()

#SST, WT vs KO, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/backgroundSub_SST_WT_vs_KO_deep_INTACT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST WT mCG/mCG", ylab="SST KO mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray", asp=1)
dev.off()


#correlations
round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], PV_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.951
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L4_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.917
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L5_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.962
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], SST_KO_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.958

round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], PV_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.881
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L4_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.892
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], L5_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.93
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], SST_KO_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, gene_methylation_corrected], method="spearman", use="complete.obs"), 3) #rho=0.864



###comparing Luo 2017 pseudobulk methylomes to INTACT methylomes
#pseudobulk methylomes
mPv_Luo2017_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mPv_Luo2017_snmcseq_CA_mm9_halfOpen.bed")
mL4_Luo2017_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL4_Luo2017_snmcseq_CA_mm9_halfOpen.bed")
mL5_all_Luo2017_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL5_all_Luo2017_snmcseq_CA_mm9_halfOpen.bed")
mSst_all_Luo2017_geneBody_TSSplus3kb_mCA <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mSst_all_Luo2017_snmcseq_CA_mm9_halfOpen.bed")

mPv_Luo2017_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mPv_Luo2017_snmcseq_CG_mm9_halfOpen.bed")
mL4_Luo2017_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL4_Luo2017_snmcseq_CG_mm9_halfOpen.bed")
mL5_all_Luo2017_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL5_all_Luo2017_snmcseq_CG_mm9_halfOpen.bed")
mSst_all_Luo2017_geneBody_TSSplus3kb_mCG <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mSst_all_Luo2017_snmcseq_CG_mm9_halfOpen.bed")

mPv_Luo2017_geneBody_TSSplus3kb_mCA$V9 <- mPv_Luo2017_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
mL4_Luo2017_geneBody_TSSplus3kb_mCA$V9 <- mL4_Luo2017_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
mL5_all_Luo2017_geneBody_TSSplus3kb_mCA$V9 <- mL5_all_Luo2017_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]
mSst_all_Luo2017_geneBody_TSSplus3kb_mCA$V9 <- mSst_all_Luo2017_geneBody_TSSplus3kb_mCA[, as.numeric(V7)/as.numeric(V8)]

mPv_Luo2017_geneBody_TSSplus3kb_mCG$V9 <- mPv_Luo2017_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
mL4_Luo2017_geneBody_TSSplus3kb_mCG$V9 <- mL4_Luo2017_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
mL5_all_Luo2017_geneBody_TSSplus3kb_mCG$V9 <- mL5_all_Luo2017_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]
mSst_all_Luo2017_geneBody_TSSplus3kb_mCG$V9 <- mSst_all_Luo2017_geneBody_TSSplus3kb_mCG[, as.numeric(V7)/as.numeric(V8)]

#PV, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="PV snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#L4, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 INTACT mCA/mCA", ylab="L4 snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#L5, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCA/mCA", ylab="L5 snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#SST, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="SST snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

##mCG
#PV, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCG/mCG", ylab="PV snmcseq mCG/CA", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#L4, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 INTACT mCG/mCG", ylab="L4 snmcseq mCG/CA", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()


#L5, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCG/mCG", ylab="L5 snmcseq mCG/CA", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#SST, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCG/mCG", ylab="SST snmcseq mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
abline(a=0, b=1)
legend("topleft", legend = paste0("r=", round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#versions of the above plots without correlation labeled or lines of unity
#PV, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCA/mCA", ylab="PV snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
dev.off()


#L4, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 INTACT mCA/mCA", ylab="L4 snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
dev.off()


#L5, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCA/mCA", ylab="L5 snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
dev.off()

#SST, INTACT vs pseudobulk, mCA
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCA_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCA/mCA", ylab="SST snmcseq mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
dev.off()

##mCG
#PV, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/PV_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="PV INTACT mCG/mCG", ylab="PV snmcseq mCG/CA", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
dev.off()


#L4, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L4_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L4 INTACT mCG/mCG", ylab="L4 snmcseq mCG/CA", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
dev.off()


#L5, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/L5_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="L5 INTACT mCG/mCG", ylab="L5 snmcseq mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
dev.off()

#SST, INTACT vs pseudobulk, mCG
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/scatterplots/SST_INTACT_vs_snmcseq_WT_geneBody_TSSplus3kb_mCG_scatterplot_noLegend.png", width=2000, height=2000, res=300)
smoothScatter(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], 
              colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="SST INTACT mCG/mCG", ylab="SST snmcseq mCG/CG", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.6, col="lightgray")
dev.off()


#Spearman correlations
round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.952
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.945
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho= 0.965
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCA[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.959

round(cor(PV_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mPv_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.886
round(cor(L4_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL4_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.924
round(cor(L5_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mL5_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.925
round(cor(SST_WT_INTACT_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], mSst_all_Luo2017_geneBody_TSSplus3kb_mCG[V4 %in% coding_genes_mm9_chr1toX_minLength5kb, V9], method="spearman", use="complete.obs"), 3) #rho=0.865
