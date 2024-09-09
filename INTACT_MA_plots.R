library(data.table)
library(dplyr)
library(DESeq2)
library(ggplot2)
BiocManager::install("geneplotter", version = "3.18")
library(geneplotter)
library(gplots)


Pv_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_ko_nondedup_exon_all_rawcounts.tsv', header=TRUE)
Sst_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_ko_nondedup_exon_rawcounts_011722.tsv', header=TRUE)
L4_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/Nr5a1_ko_nondedup_exon_rawcounts.tsv', header=TRUE)
L5_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/rbp4_all_exon_nondedup_counts.tsv')

Pv_raw <- data.frame(Pv_raw, row.names="Gene")




#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")


Pv_raw_keep <- Pv_raw[rownames(PV_deseq), c(1:4, 7:10)]
Sst_raw_keep <- Sst_raw[rownames(SST_deseq),]
L4_raw_keep <- L4_raw[rownames(L4_deseq),]
L5_raw_keep <- L5_raw[rownames(L5_deseq), c(1,2,4,5, 6, 7, 9, 10)]

Pv_coldata = data.frame(names(Pv_raw_keep))
row.names(Pv_coldata) <- names(Pv_raw_keep)

condition = c("KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT")
Pv_coldata = cbind(Pv_coldata, condition)
Pv_dds <- DESeqDataSetFromMatrix(countData = Pv_raw_keep,
                                 colData = Pv_coldata,
                                 design = ~ condition)
Pv_dds <- DESeq(Pv_dds)
Pv_res <- results(Pv_dds)
resultsNames(Pv_dds)
Pv_resLFC <- lfcShrink(Pv_dds, contrast=c("condition", "KO", "WT"), type="ashr")

Pv_resLFC = cbind(data.frame(Pv_resLFC, row.names=NULL), Gene=rownames(Pv_resLFC))



####
dds_func <- function(raw_counts){
  coldata = data.frame(names(raw_counts))
  row.names(coldata) <- names(raw_counts)
  
  condition = c("KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT")
  coldata = cbind(coldata, condition)
  dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                   colData = coldata,
                                   design = ~ condition)
  dds <- DESeq(dds)
  resLFC <- lfcShrink(dds, contrast=c("condition", "KO", "WT"), type="ashr")
  return(resLFC)
}


#MA plot
Pv_dds2 <- dds_func(Pv_raw_keep)
Sst_dds <- dds_func(Sst_raw_keep)
L4_dds <- dds_func(L4_raw_keep)
L5_dds <- dds_func(L5_raw_keep)

plotMA(Pv_dds2, main="PV MA plot, ashr", ylim=c(-2,2))

plotMA(Sst_dds, main="SST MA plot, ashr", ylim=c(-2,2))

plotMA(L4_dds, main="L4 MA plot, ashr", ylim=c(-2,2))

plotMA(L5_dds, main="L5 MA plot, ashr", ylim=c(-2,2))




#pairwise MA plots
Pv_raw_WT <- data.table(Pv_raw[,7:10], keep.rownames="Gene")
Sst_raw_WT <- data.table(Sst_raw[,5:8], keep.rownames="Gene")
L4_raw_WT <- data.table(L4_raw[,5:8], keep.rownames="Gene")
L5_raw_WT <- data.table(L5_raw[,c(6, 7, 9, 10)], keep.rownames="Gene")

Pv_Sst_merge <- data.frame(left_join(Pv_raw_WT, Sst_raw_WT, by="Gene"), row.names="Gene")

Pv_Sst_merge_coldata = data.frame(names(Pv_Sst_merge))
row.names(Pv_Sst_merge_coldata) <- names(Pv_Sst_merge)

condition_merge = c("PV", "PV", "PV", "PV", "SST", "SST", "SST", "SST")
Pv_Sst_merge_coldata = cbind(Pv_Sst_merge_coldata, condition_merge)
Pv_Sst_dds <- DESeqDataSetFromMatrix(countData = Pv_Sst_merge,
                                 colData = Pv_Sst_merge_coldata,
                                 design = ~ condition_merge)
Pv_Sst_dds  <- DESeq(Pv_Sst_dds)
Pv_Sst_res <- results(Pv_Sst_dds)
resultsNames(Pv_Sst_dds )
Pv_Sst_resLFC <- lfcShrink(Pv_Sst_dds, contrast=c("condition_merge", "PV", "SST"), type="ashr")

Pv_Sst_resLFC = cbind(data.frame(Pv_Sst_resLFC, row.names=NULL), Gene=rownames(Pv_Sst_resLFC)) %>% data.table

plotMA(Pv_Sst_dds, main="PV over Sst MA plot, ashr")

subclass_pair_DE_func <- function(subclass1, subclass2, subclass1_gene_counts, subclass2_gene_counts, output_plot, output_table, ymin=-4, ymax=4, save_table=FALSE, save_plot=FALSE){
  count_merge <- data.frame(left_join(subclass1_gene_counts, subclass2_gene_counts, by="Gene"), row.names="Gene")
  
  count_merge_coldata = data.frame(names(count_merge))
  row.names(count_merge_coldata) <- names(count_merge)
  
  condition = c(rep(subclass1, 4), rep(subclass2, 4))
  count_merge_coldata = cbind(count_merge_coldata, condition)
  dds <- DESeqDataSetFromMatrix(countData = count_merge,
                                       colData = count_merge_coldata,
                                       design = ~ condition)
  dds  <- DESeq(dds)
  resLFC <- lfcShrink(dds, contrast=c("condition", subclass1, subclass2), type="ashr")
  
  resLFC = cbind(data.frame(resLFC, row.names=NULL), Gene=rownames(resLFC)) %>% data.table
  if(save_table){
    write.csv(resLFC, file=output_table, row.names=F)
  }
  
  plotMA(dds, main=paste(subclass1, "over", subclass2, "WT MA plot, ashr"), ylim=c(ymin, ymax))
  if(save_plot){
    png(paste0(output_plot, ".png"), width=2000, height=2000, res=300)
    plotMA(dds, main=paste(subclass1, "over", subclass2, "WT MA plot, ashr"), ylim=c(ymin, ymax))
    dev.off()
    
    postscript(paste0(output_plot, ".eps"), width=2000, height=2000)
    plotMA(dds, main=paste(subclass1, "over", subclass2, "WT MA plot, ashr"), ylim=c(ymin, ymax))
    dev.off()
  }
  
}

subclass_pair_DE_func(subclass1="PV", subclass2="SST", subclass1_gene_counts=Pv_raw_WT, subclass2_gene_counts=Sst_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/PV_over_SST_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_SST_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)


subclass_pair_DE_func(subclass1="PV", subclass2="L4", subclass1_gene_counts=Pv_raw_WT, subclass2_gene_counts=L4_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/PV_over_L4_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L4_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)

subclass_pair_DE_func(subclass1="PV", subclass2="L5", subclass1_gene_counts=Pv_raw_WT, subclass2_gene_counts=L5_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/PV_over_L5_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L5_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)

subclass_pair_DE_func(subclass1="SST", subclass2="L4", subclass1_gene_counts=Sst_raw_WT, subclass2_gene_counts=L4_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/SST_over_L4_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L4_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)

subclass_pair_DE_func(subclass1="SST", subclass2="L5", subclass1_gene_counts=Sst_raw_WT, subclass2_gene_counts=L5_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/SST_over_L5_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L5_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)

subclass_pair_DE_func(subclass1="L5", subclass2="L4", subclass1_gene_counts=L5_raw_WT, subclass2_gene_counts=L4_raw_WT,
                      output_table="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/L5_over_L4_WT_deseq2_output.csv",
                      output_plot="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/L5_over_L4_WT_MA_ashr_plot",
                      ymin=-5, ymax=5,
                      save_plot=TRUE, save_table=TRUE)


###
chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
coding_genes_mm9 = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')
##using TSS+3kb to TES for gene body methylation
mPv_geneBody_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mPv_Luo2017_snmcseq_CA_merged_mm9.bed")
mSst_geneBody_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mSst_all_Luo2017_snmcseq_CA_merged_mm9.bed")
mL4_geneBody_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL4_Luo2017_snmcseq_CA_merged_mm9.bed")
mL5_geneBody_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_TSS_plus_3kb_mL5_all_Luo2017_snmcseq_CA_merged_mm9.bed")

names(mPv_geneBody_TSSplus3kb_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mSst_geneBody_TSSplus3kb_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mL4_geneBody_TSSplus3kb_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mL5_geneBody_TSSplus3kb_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")

mPv_geneBody_TSSplus3kb_mCA[, Pv_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mL4_geneBody_TSSplus3kb_mCA[, L4_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mSst_geneBody_TSSplus3kb_mCA[, Sst_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mL5_geneBody_TSSplus3kb_mCA[, L5_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]


names(coding_genes_mm9) = c("chrom", "start", "end", "Gene", "transcripts", "strand")
#coding genes with chromosomes in specific chromosome list
coding_genes_mm9_chr1toX = coding_genes_mm9[chrom %in% chrom_list]
#add cell-type gene body methylation data to list of gene names, using genes at least 5kb long
coding_genes_mm9_TSSplus3kb_mCA = data.table(inner_join(x=coding_genes_mm9_chr1toX[end-start >= 5000, .(Gene)], y=mPv_geneBody_TSSplus3kb_mCA[, .(Gene, Pv_methylation)], by=c("Gene")))
coding_genes_mm9_TSSplus3kb_mCA = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA, y=mSst_geneBody_TSSplus3kb_mCA[, .(Gene, Sst_methylation)], by=c("Gene")))
coding_genes_mm9_TSSplus3kb_mCA = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA, y=mL4_geneBody_TSSplus3kb_mCA[, .(Gene, L4_methylation)], by=c("Gene")))
coding_genes_mm9_TSSplus3kb_mCA = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA, y=mL5_geneBody_TSSplus3kb_mCA[, .(Gene, L5_methylation)], by=c("Gene")))

coding_genes_mm9_TSSplus3kb_mCA_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA, y=Pv_TPM[, .(Gene, Pv_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_TPMs, y=Sst_TPM[, .(Gene, Sst_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_TPMs, y=L4_TPM[, .(Gene, L4_WT_TPM_avg)], by="Gene"))
coding_genes_mm9_TSSplus3kb_mCA_TPMs = data.table(inner_join(x=coding_genes_mm9_TSSplus3kb_mCA_TPMs, y=L5_TPM[, .(Gene, L5_WT_TPM_avg)], by="Gene"))

Pv_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/Sst_WT_TPM_avg > 5), Gene]
Pv_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]
Pv_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Pv_WT_TPM_avg >= 1) & (Pv_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

Sst_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/Pv_WT_TPM_avg > 5), Gene]
Sst_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]
Sst_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(Sst_WT_TPM_avg >= 1) & (Sst_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

L4_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/Pv_WT_TPM_avg > 5), Gene]
L4_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/Sst_WT_TPM_avg > 5), Gene]
L4_over_L5_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L4_WT_TPM_avg >= 1) & (L4_WT_TPM_avg/L5_WT_TPM_avg > 5), Gene]

L5_over_Pv_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/Pv_WT_TPM_avg> 5), Gene]
L5_over_Sst_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/Sst_WT_TPM_avg> 5), Gene]
L5_over_L4_enriched_genes = coding_genes_mm9_TSSplus3kb_mCA_TPMs[(L5_WT_TPM_avg >= 1) & (L5_WT_TPM_avg/L4_WT_TPM_avg > 5), Gene]


#png("HG_lab/Mati/GabelLab/cell_confusion_corrPlots/mPv_snmcseq_geneBody_mCA_vs_mSst_snmcseq_geneBody_TSSplus3kb_mCA_scatterplot_gene5foldDE_coding_nondedup_noLegend.png")
smoothScatter(coding_genes_mm9_TSSplus3kb_mCA_TPMs[, Pv_methylation], coding_genes_mm9_TSSplus3kb_mCA_TPMs[, Sst_methylation], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, xlab="mPv mCA/mCA", ylab="mSst mCA/CA", xlim=c(0,0.11), ylim=c(0,0.11), pch=16, cex=0.6, col="lightgray")
par(new=TRUE)
points(x=coding_genes_mm9_TSSplus3kb_mCA_TPMs[Gene %in% Pv_over_Sst_enriched_genes, Pv_methylation], y=coding_genes_mm9_TSSplus3kb_mCA_TPMs[Gene %in% Pv_over_Sst_enriched_genes, Sst_methylation], xlab="mPv mCA/CA", ylab="mSst mCA/CA", xlim=c(0, 0.11), ylim=c(0, 0.11), col="forestgreen", pch=16, cex=1)
points(x=coding_genes_mm9_TSSplus3kb_mCA_TPMs[Gene %in% Sst_over_Pv_enriched_genes, Pv_methylation], y=coding_genes_mm9_TSSplus3kb_mCA_TPMs[Gene %in% Sst_over_Pv_enriched_genes, Sst_methylation], xlab="mPv mCA/CA", ylab="mSst mCA/CA", xlim=c(0, 0.11), ylim=c(0, 0.11), col="darkorange", pch=16, cex=1)
abline(lm(coding_genes_mm9_TSSplus3kb_mCA_TPMs[, Sst_methylation]~coding_genes_mm9_TSSplus3kb_mCA_TPMs[, Pv_methylation]))
#dev.off()

#function for outputing dds table only
subclass_pair_dds_func <- function(subclass1, subclass2, subclass1_gene_counts, subclass2_gene_counts){
  count_merge <- data.frame(left_join(subclass1_gene_counts, subclass2_gene_counts, by="Gene"), row.names="Gene")
  
  count_merge_coldata = data.frame(names(count_merge))
  row.names(count_merge_coldata) <- names(count_merge)
  
  condition = c(rep(subclass1, 4), rep(subclass2, 4))
  count_merge_coldata = cbind(count_merge_coldata, condition)
  dds <- DESeqDataSetFromMatrix(countData = count_merge,
                                colData = count_merge_coldata,
                                design = ~ condition)
  dds  <- DESeq(dds)
  resLFC <- lfcShrink(dds, contrast=c("condition", subclass1, subclass2), type="ashr")
  return(resLFC)
}

#function for outputing dds table only, no shrink
noShrink_subclass_pair_dds_func <- function(subclass1, subclass2, subclass1_gene_counts, subclass2_gene_counts){
  count_merge <- data.frame(left_join(subclass1_gene_counts, subclass2_gene_counts, by="Gene"), row.names="Gene")
  
  count_merge_coldata = data.frame(names(count_merge))
  row.names(count_merge_coldata) <- names(count_merge)
  
  condition = c(rep(subclass1, 4), rep(subclass2, 4))
  count_merge_coldata = cbind(count_merge_coldata, condition)
  dds <- DESeqDataSetFromMatrix(countData = count_merge,
                                colData = count_merge_coldata,
                                design = ~ condition)
  dds  <- DESeq(dds)
  res <- results(dds, contrast=c("condition", subclass1, subclass2))
  return(res)
}

#"skyblue", "orchid", "darkorange", "forestgreen"

Pv_raw_WT_coding = Pv_raw_WT[Gene %in% coding_genes_mm9$Gene]
Sst_raw_WT_coding = Sst_raw_WT[Gene %in% coding_genes_mm9$Gene]
L4_raw_WT_coding = L4_raw_WT[Gene %in% coding_genes_mm9$Gene]
L5_raw_WT_coding = L5_raw_WT[Gene %in% coding_genes_mm9$Gene]

Pv_Sst_dds <- subclass_pair_dds_func(subclass1="PV", subclass2="SST", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=Sst_raw_WT_coding)
Pv_L4_dds <- subclass_pair_dds_func(subclass1="PV", subclass2="L4", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)
Pv_L5_dds <- subclass_pair_dds_func(subclass1="PV", subclass2="L5", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=L5_raw_WT_coding)
Sst_L4_dds <- subclass_pair_dds_func(subclass1="SST", subclass2="L4", subclass1_gene_counts=Sst_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)
Sst_L5_dds <- subclass_pair_dds_func(subclass1="SST", subclass2="L5", subclass1_gene_counts=Sst_raw_WT_coding, subclass2_gene_counts=L5_raw_WT_coding)
L5_L4_dds <- subclass_pair_dds_func(subclass1="L5", subclass2="L4", subclass1_gene_counts=L5_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)

Pv_Sst_MA_table <- plotMA(Pv_Sst_dds, returnData=TRUE)

plotMA(Pv_Sst_dds, ylim=c(-6, 6))

geneplotter::plotMA(Pv_Sst_dds, colSig="gray32", colNonSig="gray32")

Pv_Sst_resLFC = cbind(data.frame(Pv_Sst_dds, row.names=NULL), Gene=rownames(Pv_Sst_dds)) %>% data.table
Pv_L4_resLFC = cbind(data.frame(Pv_L4_dds, row.names=NULL), Gene=rownames(Pv_L4_dds)) %>% data.table
Pv_L5_resLFC = cbind(data.frame(Pv_L5_dds, row.names=NULL), Gene=rownames(Pv_L5_dds)) %>% data.table
Sst_L4_resLFC = cbind(data.frame(Sst_L4_dds, row.names=NULL), Gene=rownames(Sst_L4_dds)) %>% data.table
Sst_L5_resLFC = cbind(data.frame(Sst_L5_dds, row.names=NULL), Gene=rownames(Sst_L5_dds)) %>% data.table
L5_L4_resLFC = cbind(data.frame(L5_L4_dds, row.names=NULL), Gene=rownames(L5_L4_dds)) %>% data.table


plotMA(Pv_Sst_dds, xlim=c(0, 1e05))

plot(x=log2(Pv_Sst_resLFC$baseMean), y=Pv_Sst_resLFC$log2FoldChange, main="PV over SST WT MA plot", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4)
plot(x=Pv_Sst_MA_table$mean, y=Pv_Sst_MA_table$lfc, main="PV over SST WT MA plot2", xlab="Mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, xlim=c(0, 1e5))
plot(x=log2(Pv_Sst_MA_table$mean), y=Pv_Sst_MA_table$lfc, main="PV over SST WT MA plot2", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, ylim=c(-6, 6))


plotMA(Pv_Sst_dds, col=ifelse(rownames(Pv_Sst_dds) %in% Pv_over_Sst_genes, "forestgreen", 
                                             ifelse(rownames(Pv_Sst_dds) %in% Sst_over_Pv_genes, "darkorange", "black")))



plot(x=Pv_Sst_resLFC[, baseMean], y=Pv_Sst_resLFC[,log2FoldChange], main="PV over SST WT MA plot", xlab="Mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, log="x")

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_SST_WT_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=Pv_Sst_resLFC[, log2(baseMean)], y=Pv_Sst_resLFC[,log2FoldChange], main="PV over SST WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=Pv_Sst_resLFC[Gene %in% Pv_over_Sst_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC[Gene %in% Pv_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="forestgreen")
points(x=Pv_Sst_resLFC[Gene %in% Sst_over_Pv_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC[Gene %in% Sst_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="darkorange")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=Pv_L4_resLFC[, log2(baseMean)], y=Pv_L4_resLFC[,log2FoldChange], main="PV over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=Pv_L4_resLFC[Gene %in% Pv_over_L4_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC[Gene %in% Pv_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="forestgreen")
points(x=Pv_L4_resLFC[Gene %in% L4_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC[Gene %in% L4_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L5_WT_MA_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=Pv_L5_resLFC[, log2(baseMean)], y=Pv_L5_resLFC[,log2FoldChange], main="PV over L5 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=Pv_L5_resLFC[Gene %in% Pv_over_L5_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC[Gene %in% Pv_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="forestgreen")
points(x=Pv_L5_resLFC[Gene %in% L5_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC[Gene %in% L5_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=Sst_L4_resLFC[, log2(baseMean)], y=Sst_L4_resLFC[,log2FoldChange], main="SST over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=Sst_L4_resLFC[Gene %in% Sst_over_L4_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC[Gene %in% Sst_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="darkorange")
points(x=Sst_L4_resLFC[Gene %in% L4_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC[Gene %in% L4_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L5_WT_MA_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=Sst_L5_resLFC[, log2(baseMean)], y=Sst_L5_resLFC[,log2FoldChange], main="SST over L5 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=Sst_L5_resLFC[Gene %in% Sst_over_L5_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC[Gene %in% Sst_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="darkorange")
points(x=Sst_L5_resLFC[Gene %in% L5_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC[Gene %in% L5_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/L5_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_plot.png", width=2000, height=2000, res=300)
plot(x=L5_L4_resLFC[, log2(baseMean)], y=L5_L4_resLFC[,log2FoldChange], main="L5 over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="gray", ylim=c(-8, 8))
par(new=TRUE)
points(x=L5_L4_resLFC[Gene %in% L5_over_L4_enriched_genes, log2(baseMean)], y=L5_L4_resLFC[Gene %in% L5_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="orchid")
points(x=L5_L4_resLFC[Gene %in% L4_over_L5_enriched_genes, log2(baseMean)], y=L5_L4_resLFC[Gene %in% L4_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

plot(x=log2(Pv_L4_resLFC$baseMean), y=Pv_L4_resLFC$log2FoldChange, main="PV over L4 WT MA plot", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.4)

df <- results(Pv_Sst_dds)

coding_genes_mm9_TSSplus3kb_mCA_TPMs[]
L5_over_L4_enriched_genes


Pv_Sst_MAplot_data <- plotMA(Pv_Sst_dds, plot=FALSE)

length(Pv_over_Sst_enriched_genes)
length(Sst_over_Pv_enriched_genes)

DEG_lengths <- rbind(
  c(0, length(Pv_over_Sst_enriched_genes), length(Pv_over_L4_enriched_genes), length(Pv_over_L5_enriched_genes)),
  c(length(Sst_over_Pv_enriched_genes), 0, length(Sst_over_L4_enriched_genes), length(Sst_over_L5_enriched_genes)),
  c(length(L4_over_Pv_enriched_genes), length(L4_over_Sst_enriched_genes), 0, length(L4_over_L5_enriched_genes)),
  c(length(L5_over_Pv_enriched_genes), length(L5_over_Sst_enriched_genes), length(L5_over_L4_enriched_genes), 0)
)

rownames(DEG_lengths) = c("PV", "SST", "L4", "L5")
colnames(DEG_lengths) = c("PV", "SST", "L4", "L5")
DEG_lengths <- as.matrix(DEG_lengths)


col_palette <- colorRampPalette(c("white", "red"))(n = 299)

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/heatmaps/INTACT_betweenSubclass_DEG_numbers_heatmap.eps")
heatmap.2(DEG_lengths,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=col_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE) 
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_SST_WT_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_Sst_resLFC[, log2(baseMean)], y=Pv_Sst_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over SST WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_Sst_resLFC[Gene %in% Pv_over_Sst_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC[Gene %in% Pv_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_Sst_resLFC[Gene %in% Sst_over_Pv_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC[Gene %in% Sst_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_L4_resLFC[, log2(baseMean)], y=Pv_L4_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_L4_resLFC[Gene %in% Pv_over_L4_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC[Gene %in% Pv_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_L4_resLFC[Gene %in% L4_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC[Gene %in% L4_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L5_WT_MA_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_L5_resLFC[, log2(baseMean)], y=Pv_L5_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over L5 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_L5_resLFC[Gene %in% Pv_over_L5_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC[Gene %in% Pv_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_L5_resLFC[Gene %in% L5_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC[Gene %in% L5_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Sst_L4_resLFC[, log2(baseMean)], y=Sst_L4_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="SST over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Sst_L4_resLFC[Gene %in% Sst_over_L4_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC[Gene %in% Sst_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
points(x=Sst_L4_resLFC[Gene %in% L4_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC[Gene %in% L4_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L5_WT_MA_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Sst_L5_resLFC[, log2(baseMean)], y=Sst_L5_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="SST over L5 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Sst_L5_resLFC[Gene %in% Sst_over_L5_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC[Gene %in% Sst_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
points(x=Sst_L5_resLFC[Gene %in% L5_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC[Gene %in% L5_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/L5_over_L4_WT_MA_ashr_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=L5_L4_resLFC[, log2(baseMean)], y=L5_L4_resLFC[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="L5 over L4 WT MA plot, ashr", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=L5_L4_resLFC[Gene %in% L5_over_L4_enriched_genes, log2(baseMean)], y=L5_L4_resLFC[Gene %in% L5_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
points(x=L5_L4_resLFC[Gene %in% L4_over_L5_enriched_genes, log2(baseMean)], y=L5_L4_resLFC[Gene %in% L4_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

###no-shrink deseq2 tables
Pv_Sst_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="PV", subclass2="SST", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=Sst_raw_WT_coding)
Pv_L4_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="PV", subclass2="L4", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)
Pv_L5_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="PV", subclass2="L5", subclass1_gene_counts=Pv_raw_WT_coding, subclass2_gene_counts=L5_raw_WT_coding)
Sst_L4_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="SST", subclass2="L4", subclass1_gene_counts=Sst_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)
Sst_L5_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="SST", subclass2="L5", subclass1_gene_counts=Sst_raw_WT_coding, subclass2_gene_counts=L5_raw_WT_coding)
L5_L4_dds_noShrink <- noShrink_subclass_pair_dds_func(subclass1="L5", subclass2="L4", subclass1_gene_counts=L5_raw_WT_coding, subclass2_gene_counts=L4_raw_WT_coding)

Pv_Sst_resLFC_noShrink = cbind(data.frame(Pv_Sst_dds_noShrink, row.names=NULL), Gene=rownames(Pv_Sst_dds_noShrink)) %>% data.table
Pv_L4_resLFC_noShrink = cbind(data.frame(Pv_L4_dds_noShrink, row.names=NULL), Gene=rownames(Pv_L4_dds_noShrink)) %>% data.table
Pv_L5_resLFC_noShrink = cbind(data.frame(Pv_L5_dds_noShrink, row.names=NULL), Gene=rownames(Pv_L5_dds_noShrink)) %>% data.table
Sst_L4_resLFC_noShrink = cbind(data.frame(Sst_L4_dds_noShrink, row.names=NULL), Gene=rownames(Sst_L4_dds_noShrink)) %>% data.table
Sst_L5_resLFC_noShrink = cbind(data.frame(Sst_L5_dds_noShrink, row.names=NULL), Gene=rownames(Sst_L5_dds_noShrink)) %>% data.table
L5_L4_resLFC_noShrink = cbind(data.frame(L5_L4_dds_noShrink, row.names=NULL), Gene=rownames(L5_L4_dds_noShrink)) %>% data.table


png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_SST_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_Sst_resLFC_noShrink[, log2(baseMean)], y=Pv_Sst_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over SST WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_Sst_resLFC_noShrink[Gene %in% Pv_over_Sst_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC_noShrink[Gene %in% Pv_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_Sst_resLFC_noShrink[Gene %in% Sst_over_Pv_enriched_genes, log2(baseMean)], y=Pv_Sst_resLFC_noShrink[Gene %in% Sst_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L4_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_L4_resLFC_noShrink[, log2(baseMean)], y=Pv_L4_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over L4 WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_L4_resLFC_noShrink[Gene %in% Pv_over_L4_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC_noShrink[Gene %in% Pv_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_L4_resLFC_noShrink[Gene %in% L4_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L4_resLFC_noShrink[Gene %in% L4_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_over_L5_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_L5_resLFC_noShrink[, log2(baseMean)], y=Pv_L5_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV over L5 WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Pv_L5_resLFC_noShrink[Gene %in% Pv_over_L5_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC_noShrink[Gene %in% Pv_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="forestgreen")
points(x=Pv_L5_resLFC_noShrink[Gene %in% L5_over_Pv_enriched_genes, log2(baseMean)], y=Pv_L5_resLFC_noShrink[Gene %in% L5_over_Pv_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L4_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Sst_L4_resLFC_noShrink[, log2(baseMean)], y=Sst_L4_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="SST over L4 WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Sst_L4_resLFC_noShrink[Gene %in% Sst_over_L4_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC_noShrink[Gene %in% Sst_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
points(x=Sst_L4_resLFC_noShrink[Gene %in% L4_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L4_resLFC_noShrink[Gene %in% L4_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/SST_over_L5_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Sst_L5_resLFC_noShrink[, log2(baseMean)], y=Sst_L5_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="SST over L5 WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=Sst_L5_resLFC_noShrink[Gene %in% Sst_over_L5_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC_noShrink[Gene %in% Sst_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="darkorange")
points(x=Sst_L5_resLFC_noShrink[Gene %in% L5_over_Sst_enriched_genes, log2(baseMean)], y=Sst_L5_resLFC_noShrink[Gene %in% L5_over_Sst_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/L5_over_L4_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=L5_L4_resLFC_noShrink[, log2(baseMean)], y=L5_L4_resLFC_noShrink[,log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="L5 over L4 WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(-4,15), ylim=c(-8, 8), asp=1)
par(new=TRUE)
points(x=L5_L4_resLFC_noShrink[Gene %in% L5_over_L4_enriched_genes, log2(baseMean)], y=L5_L4_resLFC_noShrink[Gene %in% L5_over_L4_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orchid")
points(x=L5_L4_resLFC_noShrink[Gene %in% L4_over_L5_enriched_genes, log2(baseMean)], y=L5_L4_resLFC_noShrink[Gene %in% L4_over_L5_enriched_genes,log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="skyblue")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()


L5_L4_resLFC_noShrink_more <- left_join(x=L5_L4_resLFC_noShrink, y=coding_genes_mm9_TSSplus3kb_mCA_TPMs, by="Gene")

L5_L4_resLFC_noShrink_more[Gene %in% L5_over_L4_enriched_genes]