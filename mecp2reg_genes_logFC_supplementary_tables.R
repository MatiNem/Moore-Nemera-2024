library(data.table)
library(dplyr)
library(ggplot2)


#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")
SST_deseq2 <- data.table(SST_deseq, keep.rownames="Gene")
L4_deseq2 <- data.table(L4_deseq, keep.rownames="Gene")
L5_deseq2 <- data.table(L5_deseq, keep.rownames="Gene")

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

#column names from deseq tables to include  in supplementary table
minimial_cols <- c("Gene", "ashr_log2FoldChange", "ashr_lfcSE", "ashr_pvalue", "ashr_padj")
#
PV_Mecp2Reg_genes_deseq <- rbind(
  cbind(PV_deseq2[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="PV MR gene"),
  cbind(PV_deseq2[Gene %in% pv_ma_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="PV MA gene"),
  cbind(PV_deseq2[Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4, ..minimial_cols], gene_class="PV unchanged gene")
) %>% data.table


SST_Mecp2Reg_genes_deseq <- rbind(
  cbind(SST_deseq2[Gene %in% sst_mr_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="SST MR gene"),
  cbind(SST_deseq2[Gene %in% sst_ma_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="SST MA gene"),
  cbind(SST_deseq2[Gene %in% sst_unchanged_genes_p0.5_nondedup_mm9$V4, ..minimial_cols], gene_class="SST unchanged gene")
) %>% data.table

L4_Mecp2Reg_genes_deseq <- rbind(
  cbind(L4_deseq2[Gene %in% L4_mr_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="L4 MR gene"),
  cbind(L4_deseq2[Gene %in% L4_ma_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="L4 MA gene"),
  cbind(L4_deseq2[Gene %in% L4_unchanged_genes_p0.5_nondedup_mm9$V4, ..minimial_cols], gene_class="L4 unchanged gene")
) %>% data.table

L5_Mecp2Reg_genes_deseq <- rbind(
  cbind(L5_deseq2[Gene %in% L5_mr_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="L5 MR gene"),
  cbind(L5_deseq2[Gene %in% L5_ma_genes_q0.1_nondedup_mm9$V4, ..minimial_cols], gene_class="L5 MA gene"),
  cbind(L5_deseq2[Gene %in% L5_unchanged_genes_p0.5_nondedup_mm9$V4, ..minimial_cols], gene_class="L5 unchanged gene")
) %>% data.table

write.csv(PV_Mecp2Reg_genes_deseq, "HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/PV_mecp2Reg_genes_coding_prefilt5_nondedup_deseq_table.csv", row.names=F)
write.csv(SST_Mecp2Reg_genes_deseq, "HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/SST_mecp2Reg_genes_coding_prefilt5_nondedup_deseq_table.csv", row.names=F)
write.csv(L4_Mecp2Reg_genes_deseq, "HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/L4_mecp2Reg_genes_coding_prefilt5_nondedup_deseq_table.csv", row.names=F)
write.csv(L5_Mecp2Reg_genes_deseq, "HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/L5_mecp2Reg_genes_coding_prefilt5_nondedup_deseq_table.csv", row.names=F)

#adding core MeCP2-regulated genes
core_MR_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
core_MA_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MA_genes_geneColumn4_mm9.bed")
#core_unchanged_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/MeCP2_unchanged_meta_genes_mm9.bed")

names(core_MR_genes_mm9) <- c("chrom", "start", "end", "gene", "gene_length", "strand")
names(core_MA_genes_mm9) <- c("chrom", "start", "end", "gene", "gene_length", "strand")
names(core_unchanged_genes_mm9) <- c("chrom", "start", "end", "gene", "gene_length", "strand")

master_gene_list <- unique(c(PV_Mecp2Reg_genes_deseq$Gene, SST_Mecp2Reg_genes_deseq$Gene, 
                          L4_Mecp2Reg_genes_deseq$Gene, L5_Mecp2Reg_genes_deseq$Gene,
                          core_MR_genes_mm9$gene, core_MA_genes_mm9$gene, core_unchanged_genes_mm9$gene))


all_genes_mm9 <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_flat.txt")
coding_genes_mm9 <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")

length(setdiff(all_genes_mm9$V4, c(core_MR_genes_mm9$gene, core_MA_genes_mm9$gene)))

#defining unchanged genes that aren't significant
L4_nonsig_genes <- coding_genes_mm9[!(V4 %in% c(L4_mr_genes_q0.1_nondedup_mm9$V4, L4_ma_genes_q0.1_nondedup_mm9$V4, L4_unchanged_genes_p0.5_nondedup_mm9$V4))]
L5_nonsig_genes <- coding_genes_mm9[!(V4 %in% c(L5_mr_genes_q0.1_nondedup_mm9$V4, L5_ma_genes_q0.1_nondedup_mm9$V4, L5_unchanged_genes_p0.5_nondedup_mm9$V4))]
pv_nonsig_genes <- coding_genes_mm9[!(V4 %in% c(pv_mr_genes_q0.1_nondedup_mm9$V4, pv_ma_genes_q0.1_nondedup_mm9$V4, pv_unchanged_genes_p0.5_nondedup_mm9$V4))]
sst_nonsig_genes <- coding_genes_mm9[!(V4 %in% c(sst_mr_genes_q0.1_nondedup_mm9$V4, sst_ma_genes_q0.1_nondedup_mm9$V4, sst_unchanged_genes_p0.5_nondedup_mm9$V4))]

#coding gene that is not core MR or core MA
nonCore_coding_genes_mm9 <- coding_genes_mm9[!(V4 %in% c(core_MR_genes_mm9$gene, core_MA_genes_mm9$gene))]

INTACT_core_gene_labels <- rbind(
  cbind(gene=coding_genes_mm9[V4 %in% L4_mr_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MR", category="L4"),
  cbind(gene=coding_genes_mm9[V4 %in% L4_ma_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MA", category="L4"),
  cbind(gene=coding_genes_mm9[V4 %in% L4_unchanged_genes_p0.5_nondedup_mm9$V4, V4], gene_class = "Unchanged", category="L4"),
  cbind(gene=coding_genes_mm9[V4 %in% L4_nonsig_genes$V4, V4], gene_class = "Nonsignificant", category="L4"),
  cbind(gene=coding_genes_mm9[V4 %in% L5_mr_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MR", category="L5"),
  cbind(gene=coding_genes_mm9[V4 %in% L5_ma_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MA", category="L5"),
  cbind(gene=coding_genes_mm9[V4 %in% L5_unchanged_genes_p0.5_nondedup_mm9$V4, V4], gene_class = "Unchanged", category="L5"),
  cbind(gene=coding_genes_mm9[V4 %in% L5_nonsig_genes$V4, V4], gene_class = "Nonsignificant", category="L5"),
  cbind(gene=coding_genes_mm9[V4 %in% pv_mr_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MR", category="PV"),
  cbind(gene=coding_genes_mm9[V4 %in% pv_ma_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MA", category="PV"),
  cbind(gene=coding_genes_mm9[V4 %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4, V4], gene_class = "Unchanged", category="PV"),
  cbind(gene=coding_genes_mm9[V4 %in% pv_nonsig_genes$V4, V4], gene_class = "Nonsignificant", category="PV"),
  cbind(gene=coding_genes_mm9[V4 %in% sst_mr_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MR", category="SST"),
  cbind(gene=coding_genes_mm9[V4 %in% sst_ma_genes_q0.1_nondedup_mm9$V4, V4], gene_class = "MA", category="SST"),
  cbind(gene=coding_genes_mm9[V4 %in% sst_unchanged_genes_p0.5_nondedup_mm9$V4, V4], gene_class = "Unchanged", category="SST"),
  cbind(gene=coding_genes_mm9[V4 %in% sst_nonsig_genes$V4, V4], gene_class = "Nonsignificant", category="SST"),
  cbind(gene=coding_genes_mm9[V4 %in% core_MR_genes_mm9$gene, V4], gene_class = "MR", category="Core"),
  cbind(gene=coding_genes_mm9[V4 %in% core_MA_genes_mm9$gene, V4], gene_class = "MA", category="Core"),
  cbind(gene=coding_genes_mm9[V4 %in% nonCore_coding_genes_mm9$V4, V4], gene_class = "Unchanged", category="Core")
) %>% data.table


INTACT_core_gene_labels_dcast <- dcast(INTACT_core_gene_labels, gene ~ category, value.var="gene_class") %>% data.table
INTACT_core_gene_labels_dcast <- INTACT_core_gene_labels_dcast[, .(gene, L4, L5, PV, SST, Core)]

View(INTACT_core_gene_labels_dcast[gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4])
View(INTACT_core_gene_labels_dcast[gene %in% L4_ma_genes_q0.1_nondedup_mm9$V4])
View(INTACT_core_gene_labels_dcast[gene %in% L5_unchanged_genes_p0.5_nondedup_mm9$V4])

write.csv(INTACT_core_gene_labels_dcast, "HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/INTACT_coding_sig_and_nonsig_genes_with_core_table.csv", quote=FALSE, row.names=FALSE)
