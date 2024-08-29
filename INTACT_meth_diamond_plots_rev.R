library(data.table)
library(dplyr)
library(ggplot2)

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")


#whole gene, mCA
PV_INTACT_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
SST_INTACT_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
L4_INTACT_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
L5_INTACT_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_L5_WT_KO_deep_INTACT_mCA_mm9.bed")

names(PV_INTACT_gene_body_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(SST_INTACT_gene_body_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L4_INTACT_gene_body_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L5_INTACT_gene_body_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")

PV_INTACT_gene_body_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
SST_INTACT_gene_body_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L4_INTACT_gene_body_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L5_INTACT_gene_body_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

#subtracting background non-conversion rate
PV_INTACT_gene_body_mCA$gene_methylation_corrected <- PV_INTACT_gene_body_mCA$gene_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
SST_INTACT_gene_body_mCA$gene_methylation_corrected <- SST_INTACT_gene_body_mCA$gene_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
L4_INTACT_gene_body_mCA$gene_methylation_corrected <- L4_INTACT_gene_body_mCA$gene_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L5_INTACT_gene_body_mCA$gene_methylation_corrected <- L5_INTACT_gene_body_mCA$gene_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

#making negative methylations zero for biological relevance
PV_INTACT_gene_body_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]

all_INTACT_geneBody_mCA <- rbind(cbind(L4_INTACT_gene_body_mCA, id.var="L4"),
                                  cbind(L5_INTACT_gene_body_mCA, id.var="L5"),
                                  cbind(PV_INTACT_gene_body_mCA, id.var="PV"),
                                  cbind(SST_INTACT_gene_body_mCA, id.var="SST"))

all_INTACT_geneBody_mCA = all_INTACT_geneBody_mCA %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))



#whole gene, mCG
PV_INTACT_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
SST_INTACT_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
L4_INTACT_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
L5_INTACT_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_L5_WT_KO_deep_INTACT_mCG_mm9.bed")

names(PV_INTACT_gene_body_mCG) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(SST_INTACT_gene_body_mCG) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L4_INTACT_gene_body_mCG) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L5_INTACT_gene_body_mCG) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")

PV_INTACT_gene_body_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
SST_INTACT_gene_body_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L4_INTACT_gene_body_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L5_INTACT_gene_body_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

#subtracting background non-conversion rate
PV_INTACT_gene_body_mCG$gene_methylation_corrected <- PV_INTACT_gene_body_mCG$gene_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
SST_INTACT_gene_body_mCG$gene_methylation_corrected <- SST_INTACT_gene_body_mCG$gene_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
L4_INTACT_gene_body_mCG$gene_methylation_corrected <- L4_INTACT_gene_body_mCG$gene_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L5_INTACT_gene_body_mCG$gene_methylation_corrected <- L5_INTACT_gene_body_mCG$gene_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

#making negative methylations zero for biological relevance
PV_INTACT_gene_body_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]



all_INTACT_geneBody_mCG <- rbind(cbind(L4_INTACT_gene_body_mCG, id.var="L4"),
                                 cbind(L5_INTACT_gene_body_mCG, id.var="L5"),
                                 cbind(PV_INTACT_gene_body_mCG, id.var="PV"),
                                 cbind(SST_INTACT_gene_body_mCG, id.var="SST"))

all_INTACT_geneBody_mCG = all_INTACT_geneBody_mCG %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))

#mean and standard deviation of genic mCA levels
mean_data_full_mCA <- aggregate(as.numeric(gene_methylation_corrected) ~ id.var, data = all_INTACT_geneBody_mCA, FUN = mean)
se_data_full_mCA <- aggregate(as.numeric(gene_methylation_corrected) ~ id.var, data = all_INTACT_geneBody_mCA, FUN = function(x) sd(x) / sqrt(length(x)))

data_full_mCA <- merge(mean_data_full_mCA, se_data_full_mCA, by = "id.var", suffixes = c("_mean", "_se"))
names(data_full_mCA) <- c("id.var", "gene_methylation_corrected_mean", "gene_methylation_corrected_se")
# Plot with ggplot
ggplot(data_full_mCA, aes(x = id.var, y = gene_methylation_corrected_mean, color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = gene_methylation_corrected_mean - gene_methylation_corrected_se, ymax = gene_methylation_corrected_mean + gene_methylation_corrected_se), 
                width = 0.6) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean genic mCA/CA")+xlab("")+
  coord_cartesian(ylim = c(0, 0.05)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


#mean and standard deviation of global mCG levels
mean_data_full_mCG <- aggregate(as.numeric(gene_methylation_corrected) ~ id.var, data = all_INTACT_geneBody_mCG, FUN = mean)
se_data_full_mCG <- aggregate(as.numeric(gene_methylation_corrected) ~ id.var, data = all_INTACT_geneBody_mCG, FUN = function(x) sd(x) / sqrt(length(x)))

data_full_mCG <- merge(mean_data_full_mCG, se_data_full_mCG, by = "id.var", suffixes = c("_mean", "_se"))
names(data_full_mCG) <- c("id.var", "gene_methylation_corrected_mean", "gene_methylation_corrected_se")
# Plot with ggplot
ggplot(data_full_mCG, aes(x = id.var, y = gene_methylation_corrected_mean, color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = gene_methylation_corrected_mean - gene_methylation_corrected_se, ymax = gene_methylation_corrected_mean + gene_methylation_corrected_se), 
                width = 0.6) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean genic mCG/CG")+xlab("")+
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')



PV_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
SST_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L4_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L5_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4

PV_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
SST_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L4_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L5_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4


#number of dysregulated genes
num_dys_subclass <- rbind(
  cbind(id.var="L4", num_dys_gene=length(c(L4_MR_genes, L4_MA_genes))),
  cbind(id.var="L5", num_dys_gene=length(c(L5_MR_genes, L5_MA_genes))),
  cbind(id.var="PV", num_dys_gene=length(c(PV_MR_genes, PV_MA_genes))),
  cbind(id.var="SST", num_dys_gene=length(c(SST_MR_genes, SST_MA_genes)))
) %>% data.table

data_full_mCA_numDys <- left_join(x=data_full_mCA, y=num_dys_subclass, by=c("id.var"))
data_full_mCG_numDys <- left_join(x=data_full_mCG, y=num_dys_subclass, by=c("id.var"))

#mCA
ggplot(data_full_mCA_numDys, aes(x = as.numeric(gene_methylation_corrected_mean), y = as.numeric(num_dys_gene), color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Number of dysregulated genes")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 600)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_numDys_genes_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_numDys_genes_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

#mCG
ggplot(data_full_mCG_numDys, aes(x = as.numeric(gene_methylation_corrected_mean), y = as.numeric(num_dys_gene), color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Number of dysregulated genes")+xlab("Mean genic mCG/CG")+
  coord_cartesian(xlim=c(0.0,1), ylim = c(0, 600)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_vs_numDys_genes_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_vs_numDys_genes_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")
SST_deseq2 <- data.table(SST_deseq, keep.rownames="Gene")
L4_deseq2 <- data.table(L4_deseq, keep.rownames="Gene")
L5_deseq2 <- data.table(L5_deseq, keep.rownames="Gene")



#core MR genes
core_MR_genes_mm9 <- fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
names(core_MR_genes_mm9) <- c("chrom", "start", "end", "gene", "gene_length", "strand")



#longest, highest mCA genes in genome for each subclass 
#mPv_snmcseq_geneBody_mCA_filt <- mPv_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mSst_all_snmcseq_geneBody_mCA_filt <- mSst_all_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mL4_snmcseq_geneBody_mCA_filt <- mL4_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mL5_all_snmcseq_geneBody_mCA_filt <- mL5_all_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]


###
#standard error function
std <- function(x) sd(x)/sqrt(length(x))


core_MR_deseq <- rbind(
  cbind(id.var="L4", core_logfc_mean=mean(L4_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange]), core_logfc_se=std(L4_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange])),
  cbind(id.var="L5",core_logfc_mean=mean(L5_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange]), core_logfc_se=std(L5_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange])),
  cbind(id.var="PV",core_logfc_mean=mean(PV_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange]), core_logfc_se=std(PV_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange])),
  cbind(id.var="SST",core_logfc_mean=mean(SST_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange]), core_logfc_se=std(SST_deseq2[Gene %in% core_MR_genes_mm9$gene, ashr_log2FoldChange]))
) %>% data.table

data_full_mCA_numDys_core <- left_join(x=data_full_mCA_numDys, y=core_MR_deseq, by=c("id.var"))
data_full_mCA_numDys_core$core_logfc_mean <- as.numeric(data_full_mCA_numDys_core$core_logfc_mean)
data_full_mCA_numDys_core$core_logfc_se <- as.numeric(data_full_mCA_numDys_core$core_logfc_se)

data_full_mCG_numDys_core <- left_join(x=data_full_mCG_numDys, y=core_MR_deseq, by=c("id.var"))
data_full_mCG_numDys_core$core_logfc_mean <- as.numeric(data_full_mCG_numDys_core$core_logfc_mean)
data_full_mCG_numDys_core$core_logfc_se <- as.numeric(data_full_mCG_numDys_core$core_logfc_se)

#mean genic mCA/CA vs core MR gene expression log2FC
ggplot(data_full_mCA_numDys_core, aes(x = as.numeric(gene_methylation_corrected_mean), y = as.numeric(core_logfc_mean), color = id.var)) +
  ggtitle("Core MeCP2-repressed genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = core_logfc_mean - core_logfc_se, ymax = core_logfc_mean + core_logfc_se), 
                width = 0.0025) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 0.18)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_core_MR_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_core_MR_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data_full_mCG_numDys_core, aes(x = as.numeric(gene_methylation_corrected_mean), y = as.numeric(core_logfc_mean), color = id.var)) +
  ggtitle("Core MeCP2-repressed genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = core_logfc_mean - core_logfc_se, ymax = core_logfc_mean + core_logfc_se), 
                width = 0.05) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCG/CG")+
  coord_cartesian(xlim=c(0.0,1.0), ylim = c(0, 0.18)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_vs_core_MR_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCGperCG_vs_core_MR_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

##determining longest, most high mCA/CA genes in the genome
PV_INTACT_gene_body_mCA$gene_length <- PV_INTACT_gene_body_mCA[, end-start]
SST_INTACT_gene_body_mCA$gene_length <- SST_INTACT_gene_body_mCA[, end-start]
L4_INTACT_gene_body_mCA$gene_length <- L4_INTACT_gene_body_mCA[, end-start]
L5_INTACT_gene_body_mCA$gene_length <- L5_INTACT_gene_body_mCA[, end-start]

PV_INTACT_gene_body_mCA$meth_decile <- as.character(ntile(PV_INTACT_gene_body_mCA$gene_methylation_corrected, 10))
SST_INTACT_gene_body_mCA$meth_decile <- as.character(ntile(SST_INTACT_gene_body_mCA$gene_methylation_corrected, 10))
L4_INTACT_gene_body_mCA$meth_decile <- as.character(ntile(L4_INTACT_gene_body_mCA$gene_methylation_corrected, 10))
L5_INTACT_gene_body_mCA$meth_decile <- as.character(ntile(L5_INTACT_gene_body_mCA$gene_methylation_corrected, 10))

PV_INTACT_gene_body_mCA <- PV_INTACT_gene_body_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
SST_INTACT_gene_body_mCA <- SST_INTACT_gene_body_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
L4_INTACT_gene_body_mCA <- L4_INTACT_gene_body_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
L5_INTACT_gene_body_mCA <- L5_INTACT_gene_body_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))

#filtering to longest, most highly methylated genes
PV_INTACT_gene_body_mCA_filt <- PV_INTACT_gene_body_mCA[(gene_length > 1e5) & (meth_decile ==10)]
SST_INTACT_gene_body_mCA_filt <- SST_INTACT_gene_body_mCA[(gene_length > 1e5) & (meth_decile ==10)]
L4_INTACT_gene_body_mCA_filt <- L4_INTACT_gene_body_mCA[(gene_length > 1e5) & (meth_decile ==10)]
L5_INTACT_gene_body_mCA_filt <- L5_INTACT_gene_body_mCA[(gene_length > 1e5) & (meth_decile ==10)]

#just grabbing the genes
L4_long_highmCA_genes <- L4_INTACT_gene_body_mCA_filt[, Gene]
L5_long_highmCA_genes <- L5_INTACT_gene_body_mCA_filt[, Gene]
PV_long_highmCA_genes <- PV_INTACT_gene_body_mCA_filt[, Gene]
SST_long_highmCA_genes <- SST_INTACT_gene_body_mCA_filt[, Gene]

long_highmCA_genes_deseq <- rbind(
  cbind(id.var="L4", long_highmCA_gene_logfc_mean=mean(L4_deseq2[Gene %in% L4_long_highmCA_genes, ashr_log2FoldChange]), long_highmCA_gene_logfc_se=std(L4_deseq2[Gene %in% L4_long_highmCA_genes, ashr_log2FoldChange])),
  cbind(id.var="L5", long_highmCA_gene_logfc_mean=mean(L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange]), long_highmCA_gene_logfc_se=std(L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])),
  cbind(id.var="PV", long_highmCA_gene_logfc_mean=mean(PV_deseq2[Gene %in% PV_long_highmCA_genes, ashr_log2FoldChange]), long_highmCA_gene_logfc_se=std(PV_deseq2[Gene %in% PV_long_highmCA_genes, ashr_log2FoldChange])),
  cbind(id.var="SST", long_highmCA_gene_logfc_mean=mean(SST_deseq2[Gene %in% SST_long_highmCA_genes, ashr_log2FoldChange]), long_highmCA_gene_logfc_se=std(SST_deseq2[Gene %in% SST_long_highmCA_genes, ashr_log2FoldChange]))
) %>% data.table

data_full_mCA_numDys_long_highmCA <- left_join(x=data_full_mCA_numDys, y=long_highmCA_genes_deseq, by=c("id.var"))

data_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean <- as.numeric(data_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean)
data_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_se <- as.numeric(data_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_se)

data_full_mCA_numDys_long_highmCA = data_full_mCA_numDys_long_highmCA %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))


ggplot(data_full_mCA_numDys_long_highmCA, aes(x = as.numeric(gene_methylation_corrected_mean), y = as.numeric(long_highmCA_gene_logfc_mean), color = id.var)) +
  ggtitle("Long, highly methylated genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = long_highmCA_gene_logfc_mean - long_highmCA_gene_logfc_se, ymax = long_highmCA_gene_logfc_mean + long_highmCA_gene_logfc_se), 
                width = 0.0025) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 0.15)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_long_top10percentHighmCA_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT_corrected_full_genic_mCAperCA_vs_long_top10percentHighmCA_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(L4_deseq2[Gene %in% L4_long_highmCA_genes, ashr_log2FoldChange], L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=0.000243138
wilcox.test(PV_deseq2[Gene %in% PV_long_highmCA_genes, ashr_log2FoldChange], L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=0.7181604
wilcox.test(SST_deseq2[Gene %in% SST_long_highmCA_genes, ashr_log2FoldChange], PV_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=5.258381e-10


##re-doing plots but using replicate means as for the standard error of the mean
#mCA
PV_WT_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_WT_LIB041642_deep_INTACT_mCA_mm9.bed")
PV_WT_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_WT_LIB041644_deep_INTACT_mCA_mm9.bed")
PV_KO_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_KO_LIB041641_deep_INTACT_mCA_mm9.bed")
PV_KO_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_KO_LIB041643_deep_INTACT_mCA_mm9.bed")

SST_WT_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB042979_deep_INTACT_mCA_mm9.bed")
SST_WT_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB042980_deep_INTACT_mCA_mm9.bed")
SST_WT_rep3_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB045574_deep_INTACT_mCA_mm9.bed")
SST_KO_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_KO_LIB042978_deep_INTACT_mCA_mm9.bed")
SST_KO_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_KO_LIB045575_deep_INTACT_mCA_mm9.bed")

L4_WT_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_WT_LIB041634_deep_INTACT_mCA_mm9.bed")
L4_WT_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_WT_LIB042981_deep_INTACT_mCA_mm9.bed")
L4_KO_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_KO_LIB041633_deep_INTACT_mCA_mm9.bed")
L4_KO_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_KO_LIB042982_deep_INTACT_mCA_mm9.bed")


L5_WT_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_WT_LIB041645_deep_INTACT_mCA_mm9.bed")
L5_WT_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_WT_LIB041648_deep_INTACT_mCA_mm9.bed")
L5_KO_rep1_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_KO_LIB041646_deep_INTACT_mCA_mm9.bed")
L5_KO_rep2_gene_body_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_KO_LIB041647_deep_INTACT_mCA_mm9.bed")

#mCG
PV_WT_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_WT_LIB041642_deep_INTACT_mCG_mm9.bed")
PV_WT_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_WT_LIB041644_deep_INTACT_mCG_mm9.bed")
PV_KO_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_KO_LIB041641_deep_INTACT_mCG_mm9.bed")
PV_KO_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_PV_KO_LIB041643_deep_INTACT_mCG_mm9.bed")

SST_WT_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB042979_deep_INTACT_mCG_mm9.bed")
SST_WT_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB042980_deep_INTACT_mCG_mm9.bed")
SST_WT_rep3_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_WT_LIB045574_deep_INTACT_mCG_mm9.bed")
SST_KO_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_KO_LIB042978_deep_INTACT_mCG_mm9.bed")
SST_KO_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_reps/ensgene_mm9_SST_KO_LIB045575_deep_INTACT_mCG_mm9.bed")

L4_WT_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_WT_LIB041634_deep_INTACT_mCG_mm9.bed")
L4_WT_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_WT_LIB042981_deep_INTACT_mCG_mm9.bed")
L4_KO_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_KO_LIB041633_deep_INTACT_mCG_mm9.bed")
L4_KO_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_reps/ensgene_mm9_L4_KO_LIB042982_deep_INTACT_mCG_mm9.bed")


L5_WT_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_WT_LIB041645_deep_INTACT_mCG_mm9.bed")
L5_WT_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_WT_LIB041648_deep_INTACT_mCG_mm9.bed")
L5_KO_rep1_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_KO_LIB041646_deep_INTACT_mCG_mm9.bed")
L5_KO_rep2_gene_body_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_L5_KO_LIB041647_deep_INTACT_mCG_mm9.bed")


#all nonconversion rates
lambda_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/intact_lambda_nonconversion_table.csv")


#naming methylation table columns and calculating methylation levels
meth_calc_func <- function(meth_table, label_column){
  meth_table$gene_methylation <- as.integer(meth_table[[7]])/as.integer(meth_table[[8]])
  meth_table$gene_methylation_corrected <- meth_table$gene_methylation - lambda_nonconv[label==label_column, nonconversion_rate]
  names(meth_table) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads", "gene_methylation", "gene_methylation_corrected")
  meth_table[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
  return(meth_table)
}

#mCA
PV_WT_rep1_gene_body_mCA <- meth_calc_func(meth_table=PV_WT_rep1_gene_body_mCA, label_column="PV_WT_rep1")
PV_WT_rep2_gene_body_mCA <- meth_calc_func(meth_table=PV_WT_rep2_gene_body_mCA, label_column="PV_WT_rep2")
PV_KO_rep1_gene_body_mCA <- meth_calc_func(meth_table=PV_KO_rep1_gene_body_mCA, label_column="PV_KO_rep1")
PV_KO_rep2_gene_body_mCA <- meth_calc_func(meth_table=PV_KO_rep2_gene_body_mCA, label_column="PV_KO_rep2")

SST_WT_rep1_gene_body_mCA <- meth_calc_func(meth_table=SST_WT_rep1_gene_body_mCA, label_column="SST_WT_rep1")
SST_WT_rep2_gene_body_mCA <- meth_calc_func(meth_table=SST_WT_rep2_gene_body_mCA, label_column="SST_WT_rep2")
SST_WT_rep3_gene_body_mCA <- meth_calc_func(meth_table=SST_WT_rep3_gene_body_mCA, label_column="SST_WT_rep3")
SST_KO_rep1_gene_body_mCA <- meth_calc_func(meth_table=SST_KO_rep1_gene_body_mCA, label_column="SST_KO_rep1")
SST_KO_rep2_gene_body_mCA <- meth_calc_func(meth_table=SST_KO_rep2_gene_body_mCA, label_column="SST_KO_rep2")

L4_WT_rep1_gene_body_mCA <- meth_calc_func(meth_table=L4_WT_rep1_gene_body_mCA, label_column="L4_WT_rep1")
L4_WT_rep2_gene_body_mCA <- meth_calc_func(meth_table=L4_WT_rep2_gene_body_mCA, label_column="L4_WT_rep2")
L4_KO_rep1_gene_body_mCA <- meth_calc_func(meth_table=L4_KO_rep1_gene_body_mCA, label_column="L4_KO_rep1")
L4_KO_rep2_gene_body_mCA <- meth_calc_func(meth_table=L4_KO_rep2_gene_body_mCA, label_column="L4_KO_rep2")

L5_WT_rep1_gene_body_mCA <- meth_calc_func(meth_table=L5_WT_rep1_gene_body_mCA, label_column="L5_WT_rep1")
L5_WT_rep2_gene_body_mCA <- meth_calc_func(meth_table=L5_WT_rep2_gene_body_mCA, label_column="L5_WT_rep2")
L5_KO_rep1_gene_body_mCA <- meth_calc_func(meth_table=L5_KO_rep1_gene_body_mCA, label_column="L5_KO_rep1")
L5_KO_rep2_gene_body_mCA <- meth_calc_func(meth_table=L5_KO_rep2_gene_body_mCA, label_column="L5_KO_rep2")

#mCG
PV_WT_rep1_gene_body_mCG <- meth_calc_func(meth_table=PV_WT_rep1_gene_body_mCG, label_column="PV_WT_rep1")
PV_WT_rep2_gene_body_mCG <- meth_calc_func(meth_table=PV_WT_rep2_gene_body_mCG, label_column="PV_WT_rep2")
PV_KO_rep1_gene_body_mCG <- meth_calc_func(meth_table=PV_KO_rep1_gene_body_mCG, label_column="PV_KO_rep1")
PV_KO_rep2_gene_body_mCG <- meth_calc_func(meth_table=PV_KO_rep2_gene_body_mCG, label_column="PV_KO_rep2")

SST_WT_rep1_gene_body_mCG <- meth_calc_func(meth_table=SST_WT_rep1_gene_body_mCG, label_column="SST_WT_rep1")
SST_WT_rep2_gene_body_mCG <- meth_calc_func(meth_table=SST_WT_rep2_gene_body_mCG, label_column="SST_WT_rep2")
SST_WT_rep3_gene_body_mCG <- meth_calc_func(meth_table=SST_WT_rep3_gene_body_mCG, label_column="SST_WT_rep3")
SST_KO_rep1_gene_body_mCG <- meth_calc_func(meth_table=SST_KO_rep1_gene_body_mCG, label_column="SST_KO_rep1")
SST_KO_rep2_gene_body_mCG <- meth_calc_func(meth_table=SST_KO_rep2_gene_body_mCG, label_column="SST_KO_rep2")

L4_WT_rep1_gene_body_mCG <- meth_calc_func(meth_table=L4_WT_rep1_gene_body_mCG, label_column="L4_WT_rep1")
L4_WT_rep2_gene_body_mCG <- meth_calc_func(meth_table=L4_WT_rep2_gene_body_mCG, label_column="L4_WT_rep2")
L4_KO_rep1_gene_body_mCG <- meth_calc_func(meth_table=L4_KO_rep1_gene_body_mCG, label_column="L4_KO_rep1")
L4_KO_rep2_gene_body_mCG <- meth_calc_func(meth_table=L4_KO_rep2_gene_body_mCG, label_column="L4_KO_rep2")

L5_WT_rep1_gene_body_mCG <- meth_calc_func(meth_table=L5_WT_rep1_gene_body_mCG, label_column="L5_WT_rep1")
L5_WT_rep2_gene_body_mCG <- meth_calc_func(meth_table=L5_WT_rep2_gene_body_mCG, label_column="L5_WT_rep2")
L5_KO_rep1_gene_body_mCG <- meth_calc_func(meth_table=L5_KO_rep1_gene_body_mCG, label_column="L5_KO_rep1")
L5_KO_rep2_gene_body_mCG <- meth_calc_func(meth_table=L5_KO_rep2_gene_body_mCG, label_column="L5_KO_rep2")

reps_INTACT_geneBody_mCA <- rbind(
  cbind(mean_gene_meth=PV_WT_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_WT_rep1", id.var="PV"),
  cbind(mean_gene_meth=PV_WT_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_WT_rep2", id.var="PV"),
  cbind(mean_gene_meth=PV_KO_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_KO_rep1", id.var="PV"),
  cbind(mean_gene_meth=PV_KO_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_KO_rep2", id.var="PV"),
  cbind(mean_gene_meth=SST_WT_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep1", id.var="SST"),
  cbind(mean_gene_meth=SST_WT_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep2", id.var="SST"),
  cbind(mean_gene_meth=SST_WT_rep3_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep3", id.var="SST"),
  cbind(mean_gene_meth=SST_KO_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_KO_rep1", id.var="SST"),
  cbind(mean_gene_meth=SST_KO_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_KO_rep2", id.var="SST"),
  cbind(mean_gene_meth=L4_WT_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_WT_rep1", id.var="L4"),
  cbind(mean_gene_meth=L4_WT_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_WT_rep2", id.var="L4"),
  cbind(mean_gene_meth=L4_KO_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_KO_rep1", id.var="L4"),
  cbind(mean_gene_meth=L4_KO_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_KO_rep2", id.var="L4"),
  cbind(mean_gene_meth=L5_WT_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_WT_rep1", id.var="L5"),
  cbind(mean_gene_meth=L5_WT_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_WT_rep2", id.var="L5"),
  cbind(mean_gene_meth=L5_KO_rep1_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_KO_rep1", id.var="L5"),
  cbind(mean_gene_meth=L5_KO_rep2_gene_body_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_KO_rep2", id.var="L5")
) %>% data.table

reps_INTACT_geneBody_mCG <- rbind(
  cbind(mean_gene_meth=PV_WT_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_WT_rep1", id.var="PV"),
  cbind(mean_gene_meth=PV_WT_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_WT_rep2", id.var="PV"),
  cbind(mean_gene_meth=PV_KO_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_KO_rep1", id.var="PV"),
  cbind(mean_gene_meth=PV_KO_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="PV_KO_rep2", id.var="PV"),
  cbind(mean_gene_meth=SST_WT_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep1", id.var="SST"),
  cbind(mean_gene_meth=SST_WT_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep2", id.var="SST"),
  cbind(mean_gene_meth=SST_WT_rep3_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_WT_rep3", id.var="SST"),
  cbind(mean_gene_meth=SST_KO_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_KO_rep1", id.var="SST"),
  cbind(mean_gene_meth=SST_KO_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="SST_KO_rep2", id.var="SST"),
  cbind(mean_gene_meth=L4_WT_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_WT_rep1", id.var="L4"),
  cbind(mean_gene_meth=L4_WT_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_WT_rep2", id.var="L4"),
  cbind(mean_gene_meth=L4_KO_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_KO_rep1", id.var="L4"),
  cbind(mean_gene_meth=L4_KO_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L4_KO_rep2", id.var="L4"),
  cbind(mean_gene_meth=L5_WT_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_WT_rep1", id.var="L5"),
  cbind(mean_gene_meth=L5_WT_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_WT_rep2", id.var="L5"),
  cbind(mean_gene_meth=L5_KO_rep1_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_KO_rep1", id.var="L5"),
  cbind(mean_gene_meth=L5_KO_rep2_gene_body_mCG[, mean(gene_methylation_corrected, na.rm=TRUE)], label="L5_KO_rep2", id.var="L5")
) %>% data.table


reps_INTACT_geneBody_mCA = reps_INTACT_geneBody_mCA %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))
reps_INTACT_geneBody_mCG = reps_INTACT_geneBody_mCG %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))

#grand mean and standard error of genic mCA levels across reps
mean_reps_data_full_mCA <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_INTACT_geneBody_mCA, FUN = mean)
se_reps_data_full_mCA <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_INTACT_geneBody_mCA, FUN = function(x) sd(x) / sqrt(length(x)))

data_reps_full_mCA <- merge(mean_reps_data_full_mCA, se_reps_data_full_mCA, by = "id.var", suffixes = c("_mean", "_se"))
names(data_reps_full_mCA) <- c("id.var", "gene_meth_grand_mean", "gene_meth_grand_se")
# Plot with ggplot
ggplot(data_reps_full_mCA, aes(x = id.var, y = gene_meth_grand_mean, color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = gene_meth_grand_mean - gene_meth_grand_se, ymax = gene_meth_grand_mean + gene_meth_grand_se), 
                width = 0.4) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean genic mCA/CA")+xlab("")+
  coord_cartesian(ylim = c(0, 0.05)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


#grand mean and standard error of genic mCG levels across reps
mean_reps_data_full_mCG <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_INTACT_geneBody_mCG, FUN = mean)
se_reps_data_full_mCG <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_INTACT_geneBody_mCG, FUN = function(x) sd(x) / sqrt(length(x)))

data_reps_full_mCG <- merge(mean_reps_data_full_mCG, se_reps_data_full_mCG, by = "id.var", suffixes = c("_mean", "_se"))
names(data_reps_full_mCG) <- c("id.var", "gene_meth_grand_mean", "gene_meth_grand_se")
# Plot with ggplot
ggplot(data_reps_full_mCG, aes(x = id.var, y = gene_meth_grand_mean, color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = gene_meth_grand_mean - gene_meth_grand_se, ymax = gene_meth_grand_mean + gene_meth_grand_se), 
                width = 0.4) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean genic mCG/CG")+xlab("")+
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

##
data_reps_full_mCA_numDys <- left_join(x=data_reps_full_mCA, y=num_dys_subclass, by=c("id.var"))
data_reps_full_mCG_numDys <- left_join(x=data_reps_full_mCG, y=num_dys_subclass, by=c("id.var"))

#mCA
ggplot(data_reps_full_mCA_numDys, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(num_dys_gene), color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Number of dysregulated genes")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 600)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_numDys_genes_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_numDys_genes_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

#mCG
ggplot(data_reps_full_mCG_numDys, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(num_dys_gene), color = id.var)) +
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Number of dysregulated genes")+xlab("Mean genic mCG/CG")+
  coord_cartesian(xlim=c(0.0,1), ylim = c(0, 600)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_vs_numDys_genes_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_vs_numDys_genes_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

#
data_reps_full_mCA_numDys_core <- left_join(x=data_reps_full_mCA_numDys, y=core_MR_deseq, by=c("id.var"))
data_reps_full_mCA_numDys_core$core_logfc_mean <- as.numeric(data_reps_full_mCA_numDys_core$core_logfc_mean)
data_reps_full_mCA_numDys_core$core_logfc_se <- as.numeric(data_reps_full_mCA_numDys_core$core_logfc_se)

data_reps_full_mCG_numDys_core <- left_join(x=data_reps_full_mCG_numDys, y=core_MR_deseq, by=c("id.var"))
data_reps_full_mCG_numDys_core$core_logfc_mean <- as.numeric(data_reps_full_mCG_numDys_core$core_logfc_mean)
data_reps_full_mCG_numDys_core$core_logfc_se <- as.numeric(data_reps_full_mCG_numDys_core$core_logfc_se)

#mean genic mCA/CA vs core MR gene expression log2FC
ggplot(data_reps_full_mCA_numDys_core, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(core_logfc_mean), color = id.var)) +
  ggtitle("Core MeCP2-repressed genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = core_logfc_mean - core_logfc_se, ymax = core_logfc_mean + core_logfc_se), 
                width = 0.0025) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 0.16)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_core_MR_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_core_MR_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(data_reps_full_mCG_numDys_core, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(core_logfc_mean), color = id.var)) +
  ggtitle("Core MeCP2-repressed genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = core_logfc_mean - core_logfc_se, ymax = core_logfc_mean + core_logfc_se), 
                width = 0.05) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCG/CG")+
  coord_cartesian(xlim=c(0.0,1.0), ylim = c(0, 0.16)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_vs_core_MR_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCGperCG_vs_core_MR_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


#long, highly methylated genes 
data_reps_full_mCA_numDys_long_highmCA <- left_join(x=data_reps_full_mCA_numDys, y=long_highmCA_genes_deseq, by=c("id.var"))

data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean <- as.numeric(data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean)
data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_se <- as.numeric(data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_se)

data_reps_full_mCA_numDys_long_highmCA = data_reps_full_mCA_numDys_long_highmCA %>% mutate(id.var = factor(id.var, levels=c("L4", "L5", "PV", "SST")))


ggplot(data_reps_full_mCA_numDys_long_highmCA, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(long_highmCA_gene_logfc_mean), color = id.var)) +
  ggtitle("Long, highly methylated genes")+
  geom_point(shape = 18, size = 10) +  # Diamond shape for mean
  geom_errorbar(aes(ymin = long_highmCA_gene_logfc_mean - long_highmCA_gene_logfc_se, ymax = long_highmCA_gene_logfc_mean + long_highmCA_gene_logfc_se), 
                width = 0.0025) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("L4"="skyblue", "L5"="purple3", "PV"="forestgreen", "SST"="darkorange")) +
  ylab("Mean log2 FC MeCP2 KO/WT")+xlab("Mean genic mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.05), ylim = c(0, 0.16)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_long_top10percentHighmCA_gene_Log2FC_diamond_plot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/grandMean_by_reps_INTACT_corrected_full_genic_mCAperCA_vs_long_top10percentHighmCA_gene_Log2FC_diamond_plot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(L4_deseq2[Gene %in% L4_long_highmCA_genes, ashr_log2FoldChange], L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=0.000243138
wilcox.test(PV_deseq2[Gene %in% PV_long_highmCA_genes, ashr_log2FoldChange], L5_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=0.7181604
wilcox.test(SST_deseq2[Gene %in% SST_long_highmCA_genes, ashr_log2FoldChange], PV_deseq2[Gene %in% L5_long_highmCA_genes, ashr_log2FoldChange])$p.value #p=5.258381e-10

###pearson and spearman correlations for above plots
##mCA
round(cor(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="pearson", use="complete.obs"), 3) #0.837
round(cor(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="spearman", use="complete.obs"), 3) #0.8

round(cor(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="pearson", use="complete.obs"), 3) #0.679
round(cor(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="spearman", use="complete.obs"), 3) #0.8

round(cor(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="pearson", use="complete.obs"), 3) #0.729
round(cor(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="spearman", use="complete.obs"), 3) #0.8

gene_mCA_corr <- rbind(
  cbind(comparison="Mean_gene_mCA_num_dys_genes", 
        pearson.corr=cor(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="pearson", use="complete.obs"),
        pearson.pval=cor.test(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="pearson", use="complete.obs")$p.value,
        spearman.corr=cor(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="spearman", use="complete.obs"), 
        spearman.pval=cor.test(x=data_reps_full_mCA_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCA_numDys$num_dys_gene), method="spearman", use="complete.obs")$p.value),
  cbind(comparison="Mean_gene_mCA_long_highmCA_gene_logfc_mean", 
        pearson.corr=cor(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="pearson", use="complete.obs"),
        pearson.pval=cor.test(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="pearson", use="complete.obs")$p.value,
        spearman.corr=cor(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="spearman", use="complete.obs"),
        spearman.corr=cor.test(x=data_reps_full_mCA_numDys_long_highmCA$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_long_highmCA$long_highmCA_gene_logfc_mean, method="spearman", use="complete.obs")$p.value),
  cbind(comparison="Mean_gene_mCA_core_logfc_mean", 
        pearson.corr=cor(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="pearson", use="complete.obs"),
        pearson.pval=cor.test(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="pearson", use="complete.obs")$p.value,
        spearman.corr=cor(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="spearman", use="complete.obs"),
        spearman.corr=cor.test(x=data_reps_full_mCA_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCA_numDys_core$core_logfc_mean, method="spearman", use="complete.obs")$p.value)
) %>% data.table

##mCG

gene_mCG_corr <- rbind(
  cbind(comparison="Mean_gene_mCG_num_dys_genes", 
        pearson.corr=cor(x=data_reps_full_mCG_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCG_numDys$num_dys_gene), method="pearson", use="complete.obs"),
        pearson.pval=cor.test(x=data_reps_full_mCG_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCG_numDys$num_dys_gene), method="pearson", use="complete.obs")$p.value,
        spearman.corr=cor(x=data_reps_full_mCG_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCG_numDys$num_dys_gene), method="spearman", use="complete.obs"), 
        spearman.pval=cor.test(x=data_reps_full_mCG_numDys$gene_meth_grand_mean, y=as.numeric(data_reps_full_mCG_numDys$num_dys_gene), method="spearman", use="complete.obs")$p.value),
  cbind(comparison="Mean_gene_mCG_core_logfc_mean", 
        pearson.corr=cor(x=data_reps_full_mCG_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCG_numDys_core$core_logfc_mean, method="pearson", use="complete.obs"),
        pearson.pval=cor.test(x=data_reps_full_mCG_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCG_numDys_core$core_logfc_mean, method="pearson", use="complete.obs")$p.value,
        spearman.corr=cor(x=data_reps_full_mCG_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCG_numDys_core$core_logfc_mean, method="spearman", use="complete.obs"),
        spearman.corr=cor.test(x=data_reps_full_mCG_numDys_core$gene_meth_grand_mean, y=data_reps_full_mCG_numDys_core$core_logfc_mean, method="spearman", use="complete.obs")$p.value)
) %>% data.table