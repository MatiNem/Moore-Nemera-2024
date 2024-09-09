library(data.table)
library(dplyr)
library(ggplot2)

Pv_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")
Sst_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_Mecp2KO_gene_TPMs_nondedup.txt")
L4_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_Mecp2KO_gene_TPMs_nondedup.txt")
L5_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_Mecp2KO_gene_TPMs_nondedup.txt")

#whole gene 
mPv_snmcseq_geneBody_mCA = fread("HG_lab/Mati/GabelLab/genesets/pv_Russell/all_genes_mPv_Luo2017_snmcseq_CA_merged_mm9.bed")
mSst_all_snmcseq_geneBody_mCA = fread("HG_lab/Mati/GabelLab/genesets/Sst/all_genes_mSst_all_Luo2017_snmcseq_CA_merged_mm9.bed")
mL4_snmcseq_geneBody_mCA = fread("HG_lab/Mati/GabelLab/genesets/L4/all_genes_mL4_Luo2017_snmcseq_CA_merged_mm9.bed")
mL5_all_snmcseq_geneBody_mCA = fread("HG_lab/Mati/GabelLab/genesets/L5/all_genes_mL5_all_Luo2017_snmcseq_CA_merged_mm9.bed")

names(mPv_snmcseq_geneBody_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mSst_all_snmcseq_geneBody_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mL4_snmcseq_geneBody_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(mL5_all_snmcseq_geneBody_mCA) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")

mPv_snmcseq_geneBody_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mSst_all_snmcseq_geneBody_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mL4_snmcseq_geneBody_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
mL5_all_snmcseq_geneBody_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]


mPv_snmcseq_geneBody_mCA$gene_length <- mPv_snmcseq_geneBody_mCA[, end-start]
mSst_all_snmcseq_geneBody_mCA$gene_length <- mSst_all_snmcseq_geneBody_mCA[, end-start]
mL4_snmcseq_geneBody_mCA$gene_length <- mL4_snmcseq_geneBody_mCA[, end-start]
mL5_all_snmcseq_geneBody_mCA$gene_length <- mL5_all_snmcseq_geneBody_mCA[, end-start]

mPv_snmcseq_geneBody_mCA$meth_decile <- as.character(ntile(mPv_snmcseq_geneBody_mCA$gene_methylation, 10))
mSst_all_snmcseq_geneBody_mCA$meth_decile <- as.character(ntile(mSst_all_snmcseq_geneBody_mCA$gene_methylation, 10))
mL4_snmcseq_geneBody_mCA$meth_decile <- as.character(ntile(mL4_snmcseq_geneBody_mCA$gene_methylation, 10))
mL5_all_snmcseq_geneBody_mCA$meth_decile <- as.character(ntile(mL5_all_snmcseq_geneBody_mCA$gene_methylation, 10))

#longest, highest mCA genes in genome for each subclass 
#mPv_snmcseq_geneBody_mCA_filt <- mPv_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mSst_all_snmcseq_geneBody_mCA_filt <- mSst_all_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mL4_snmcseq_geneBody_mCA_filt <- mL4_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]
#mL5_all_snmcseq_geneBody_mCA_filt <- mL5_all_snmcseq_geneBody_mCA[(meth_decile == 10) & (gene_length > 1e5)]


###
mPv_snmcseq_geneBody_mCA_filt0 <- mPv_snmcseq_geneBody_mCA[(gene_length > 1e5) & (meth_decile ==10)]
mSst_all_snmcseq_geneBody_mCA_filt0 <- mSst_all_snmcseq_geneBody_mCA[(gene_length > 1e5)& (meth_decile ==10)]
mL4_snmcseq_geneBody_mCA_filt0 <- mL4_snmcseq_geneBody_mCA[(gene_length > 1e5)& (meth_decile ==10)]
mL5_all_snmcseq_geneBody_mCA_filt0 <- mL5_all_snmcseq_geneBody_mCA[(gene_length > 1e5)& (meth_decile ==10)]

#longest, highest mCA genes in genome for each subclass 
mPv_snmcseq_geneBody_mCA_filt <- mPv_snmcseq_geneBody_mCA[(gene_length > 1e5)]
mSst_all_snmcseq_geneBody_mCA_filt <- mSst_all_snmcseq_geneBody_mCA[(gene_length > 1e5)]
mL4_snmcseq_geneBody_mCA_filt <- mL4_snmcseq_geneBody_mCA[(gene_length > 1e5)]
mL5_all_snmcseq_geneBody_mCA_filt <- mL5_all_snmcseq_geneBody_mCA[(gene_length > 1e5)]

mPv_snmcseq_geneBody_mCA_filt$meth_decile <- as.character(ntile(mPv_snmcseq_geneBody_mCA_filt$gene_methylation, 10))
mSst_all_snmcseq_geneBody_mCA_filt$meth_decile <- as.character(ntile(mSst_all_snmcseq_geneBody_mCA_filt$gene_methylation, 10))
mL4_snmcseq_geneBody_mCA_filt$meth_decile <- as.character(ntile(mL4_snmcseq_geneBody_mCA_filt$gene_methylation, 10))
mL5_all_snmcseq_geneBody_mCA_filt$meth_decile <- as.character(ntile(mL5_all_snmcseq_geneBody_mCA_filt$gene_methylation, 10))

mPv_snmcseq_geneBody_mCA_filt2 <- mPv_snmcseq_geneBody_mCA_filt[meth_decile == 10]
mSst_all_snmcseq_geneBody_mCA_filt2 <- mSst_all_snmcseq_geneBody_mCA_filt[meth_decile == 10]
mL4_snmcseq_geneBody_mCA_filt2 <- mL4_snmcseq_geneBody_mCA_filt[meth_decile == 10]
mL5_all_snmcseq_geneBody_mCA_filt2 <- mL5_all_snmcseq_geneBody_mCA_filt[meth_decile == 10]
#data tables with subclass WT and KO TPM 
Pv_TPM_melt <- melt(Pv_TPM, id.vars = "Gene", measure.vars = c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"), 
                      variable.name = "Condition", value.name = "TPM_Value") %>% data.table

Sst_TPM_melt <- melt(Sst_TPM, id.vars = "Gene", measure.vars = c("Sst_WT_TPM_avg", "Sst_KO_TPM_avg"), 
                    variable.name = "Condition", value.name = "TPM_Value") %>% data.table

L4_TPM_melt <- melt(L4_TPM, id.vars = "Gene", measure.vars = c("L4_WT_TPM_avg", "L4_KO_TPM_avg"), 
                     variable.name = "Condition", value.name = "TPM_Value") %>% data.table

L5_TPM_melt <- melt(L5_TPM, id.vars = "Gene", measure.vars = c("L5_WT_TPM_avg", "L5_KO_TPM_avg"), 
                    variable.name = "Condition", value.name = "TPM_Value") %>% data.table


# Add the Genotype column
Pv_TPM_melt$Genotype = "NA"
Pv_TPM_melt[Condition == "Pv_WT_TPM_avg", Genotype := "WT"]
Pv_TPM_melt[Condition == "Pv_KO_TPM_avg", Genotype := "KO"]

Sst_TPM_melt$Genotype = "NA"
Sst_TPM_melt[Condition == "Sst_WT_TPM_avg", Genotype := "WT"]
Sst_TPM_melt[Condition == "Sst_KO_TPM_avg", Genotype := "KO"]

L4_TPM_melt$Genotype = "NA"
L4_TPM_melt[Condition == "L4_WT_TPM_avg", Genotype := "WT"]
L4_TPM_melt[Condition == "L4_KO_TPM_avg", Genotype := "KO"]

L5_TPM_melt$Genotype = "NA"
L5_TPM_melt[Condition == "L5_WT_TPM_avg", Genotype := "WT"]
L5_TPM_melt[Condition == "L5_KO_TPM_avg", Genotype := "KO"]

INTACT_TPM_long_mCAhigh_dt <- rbind(
  cbind(Pv_TPM_melt[Gene %in% mPv_snmcseq_geneBody_mCA_filt$Gene], subclass="PV"),
  cbind(Sst_TPM_melt[Gene %in% mSst_all_snmcseq_geneBody_mCA_filt$Gene], subclass="SST"),
  cbind(L4_TPM_melt[Gene %in% mL4_snmcseq_geneBody_mCA_filt$Gene], subclass="L4"),
  cbind(L5_TPM_melt[Gene %in% mL5_all_snmcseq_geneBody_mCA_filt$Gene], subclass="L5"))

INTACT_TPM_long_mCAhigh_dt = INTACT_TPM_long_mCAhigh_dt %>% mutate(Genotype = factor(Genotype, levels=c("WT","KO")))

ggplot(INTACT_TPM_long_mCAhigh_dt, aes(x = subclass, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(0,1.0))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.long.mCAhigh.genes.log2TPMplus1.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.long.mCAhigh.genes.log2TPMplus1.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

###
#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")
SST_deseq2 <- data.table(SST_deseq, keep.rownames="Gene")
L4_deseq2 <- data.table(L4_deseq, keep.rownames="Gene")
L5_deseq2 <- data.table(L5_deseq, keep.rownames="Gene")

INTACT_deseq_long_mCAhigh <- rbind(
  cbind(L4_deseq2[Gene %in% mL4_snmcseq_geneBody_mCA_filt2$Gene], subclass="L4"),
  cbind(L5_deseq2[Gene %in% mL5_all_snmcseq_geneBody_mCA_filt2$Gene], subclass="L5"),
  cbind(PV_deseq2[Gene %in% mPv_snmcseq_geneBody_mCA_filt2$Gene], subclass="PV"),
  cbind(SST_deseq2[Gene %in% mSst_all_snmcseq_geneBody_mCA_filt2$Gene], subclass="SST"))

INTACT_deseq_long_mCAhigh_summ <- group_by(INTACT_deseq_long_mCAhigh, subclass) %>%
  summarise(
    mean_ashr_log2FoldChange = mean(ashr_log2FoldChange)
  ) %>% data.table


###
INTACT_deseq_long_mCAhigh_new <- rbind(
  cbind(L4_deseq2[Gene %in% mL4_snmcseq_geneBody_mCA_filt0$Gene], subclass="L4"),
  cbind(L5_deseq2[Gene %in% mL5_all_snmcseq_geneBody_mCA_filt0$Gene], subclass="L5"),
  cbind(PV_deseq2[Gene %in% mPv_snmcseq_geneBody_mCA_filt0$Gene], subclass="PV"),
  cbind(SST_deseq2[Gene %in% mSst_all_snmcseq_geneBody_mCA_filt0$Gene], subclass="SST"))

INTACT_deseq_long_mCAhigh_new_summ <- group_by(INTACT_deseq_long_mCAhigh_new, subclass) %>%
  summarise(
    mean_ashr_log2FoldChange = mean(ashr_log2FoldChange)
  ) %>% data.table


INTACT_TPM_long_mCAhigh_new_dt <- rbind(
  cbind(Pv_TPM_melt[Gene %in% mPv_snmcseq_geneBody_mCA_filt0$Gene], subclass="PV"),
  cbind(Sst_TPM_melt[Gene %in% mSst_all_snmcseq_geneBody_mCA_filt0$Gene], subclass="SST"),
  cbind(L4_TPM_melt[Gene %in% mL4_snmcseq_geneBody_mCA_filt0$Gene], subclass="L4"),
  cbind(L5_TPM_melt[Gene %in% mL5_all_snmcseq_geneBody_mCA_filt0$Gene], subclass="L5"))

INTACT_TPM_long_mCAhigh_new_dt = INTACT_TPM_long_mCAhigh_new_dt %>% mutate(Genotype = factor(Genotype, levels=c("WT","KO")))

ggplot(INTACT_TPM_long_mCAhigh_new_dt, aes(x = subclass, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("Long (>100kb) high-mCA (top 10%) genes\nof subclasses")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(0,1.0))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.long.mCAhigh.genes.log2TPMplus1.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.long.mCAhigh.genes.log2TPMplus1.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##MR genes
PV_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
SST_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L4_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L5_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4

INTACT_TPM_MR_genes <- rbind(
  cbind(Pv_TPM_melt[Gene %in% PV_MR_genes], subclass="PV"),
  cbind(Sst_TPM_melt[Gene %in% SST_MR_genes], subclass="SST"),
  cbind(L4_TPM_melt[Gene %in% L4_MR_genes], subclass="L4"),
  cbind(L5_TPM_melt[Gene %in% L5_MR_genes], subclass="L5"))

INTACT_TPM_MR_genes = INTACT_TPM_MR_genes %>% mutate(Genotype = factor(Genotype, levels=c("WT","KO")))

ggplot(INTACT_TPM_MR_genes, aes(x = subclass, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("MR genes of subclasses")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(0,1.0))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.MR.genes.log2TPMplus1.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_boxplots/INTACT.subclass.WT.KO.MR.genes.log2TPMplus1.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


INTACT_TPM_MR_genes[Genotype=="WT"]
wilcox.test(INTACT_TPM_MR_genes[(subclass=="L4") & (Genotype=="WT"), TPM_Value], INTACT_TPM_MR_genes[(subclass=="L4") & (Genotype=="KO"), TPM_Value])$p.value
wilcox.test(INTACT_TPM_MR_genes[(subclass=="L5") & (Genotype=="WT"), TPM_Value], INTACT_TPM_MR_genes[(subclass=="L5") & (Genotype=="KO"), TPM_Value])$p.value
wilcox.test(INTACT_TPM_MR_genes[(subclass=="PV") & (Genotype=="WT"), TPM_Value], INTACT_TPM_MR_genes[(subclass=="PV") & (Genotype=="KO"), TPM_Value])$p.value
wilcox.test(INTACT_TPM_MR_genes[(subclass=="SST") & (Genotype=="WT"), TPM_Value], INTACT_TPM_MR_genes[(subclass=="SST") & (Genotype=="KO"), TPM_Value])$p.value


##MA genes
PV_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
SST_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L4_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4
L5_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4

INTACT_TPM_MA_genes <- rbind(
  cbind(Pv_TPM_melt[Gene %in% PV_MA_genes], subclass="PV"),
  cbind(Sst_TPM_melt[Gene %in% SST_MA_genes], subclass="SST"),
  cbind(L4_TPM_melt[Gene %in% L4_MA_genes], subclass="L4"),
  cbind(L5_TPM_melt[Gene %in% L5_MA_genes], subclass="L5"))

INTACT_TPM_MA_genes = INTACT_TPM_MA_genes %>% mutate(Genotype = factor(Genotype, levels=c("WT","KO")))


##try doing violin plot

#violin plot of MR gene TPMs
ggplot(INTACT_TPM_MR_genes, aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("MR genes of subclasses")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.genes.log2TPMplus1.violinplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.genes.log2TPMplus1.violinplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

#violin plot of MA gene TPMs
ggplot(INTACT_TPM_MA_genes, aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("MA genes of subclasses")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MA.genes.log2TPMplus1.violinplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MA.genes.log2TPMplus1.violinplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


##all other genes
coding_genes_mm9 <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")
names(coding_genes_mm9) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")

all_other_genes_PV <- setdiff(coding_genes_mm9$gene,  c(PV_MR_genes, PV_MA_genes))
all_other_genes_SST <- setdiff(coding_genes_mm9$gene,  c(SST_MR_genes, SST_MA_genes))
all_other_genes_L4 <- setdiff(coding_genes_mm9$gene,  c(L4_MR_genes, L4_MA_genes))
all_other_genes_L5 <- setdiff(coding_genes_mm9$gene,  c(L5_MR_genes, L5_MA_genes))

#all genes
INTACT_TPM_all_genes <- rbind(
  cbind(Pv_TPM_melt[Gene %in% all_other_genes_PV], subclass="PV", gene_class="All other genes"),
  cbind(Sst_TPM_melt[Gene %in% all_other_genes_SST], subclass="SST", gene_class="All other genes"),
  cbind(L4_TPM_melt[Gene %in% all_other_genes_L4], subclass="L4", gene_class="All other genes"),
  cbind(L5_TPM_melt[Gene %in% all_other_genes_L5 ], subclass="L5", gene_class="All other genes"),
  cbind(Pv_TPM_melt[Gene %in% PV_MR_genes], subclass="PV", gene_class="MR genes"),
  cbind(Sst_TPM_melt[Gene %in% SST_MR_genes], subclass="SST", gene_class="MR genes"),
  cbind(L4_TPM_melt[Gene %in% L4_MR_genes], subclass="L4", gene_class="MR genes"),
  cbind(L5_TPM_melt[Gene %in% L5_MR_genes], subclass="L5", gene_class="MR genes"),
  cbind(Pv_TPM_melt[Gene %in% PV_MA_genes], subclass="PV", gene_class="MA genes"),
  cbind(Sst_TPM_melt[Gene %in% SST_MA_genes], subclass="SST", gene_class="MA genes"),
  cbind(L4_TPM_melt[Gene %in% L4_MA_genes], subclass="L4", gene_class="MA genes"),
  cbind(L5_TPM_melt[Gene %in% L5_MA_genes], subclass="L5", gene_class="MA genes"))

#mutate to set plotting order
INTACT_TPM_all_genes = INTACT_TPM_all_genes %>% mutate(Genotype = factor(Genotype, levels=c("WT","KO")))
INTACT_TPM_all_genes = INTACT_TPM_all_genes %>% mutate(subclass = factor(subclass, levels=c("L4","L5", "PV", "SST")))
INTACT_TPM_all_genes = INTACT_TPM_all_genes %>% mutate(gene_class = factor(gene_class, levels=c("All other genes","MR genes", "MA genes")))

##try doing violin plot

#violin plot of MR gene TPMs
ggplot(INTACT_TPM_all_genes, aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass + gene_class,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = median, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.all.genes.log2TPMplus1.violinplot.png", width = 10, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.all.genes.log2TPMplus1.violinplot.eps", width = 10, height = 5, dpi = 300, units = "in", device='eps')


#no all other
ggplot(INTACT_TPM_all_genes[gene_class %in% c("MR genes", "MA genes")], aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("MR and MA only")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass + gene_class,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = median, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.MA.genes.log2TPMplus1.violinplot.png", width = 10, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.MA.genes.log2TPMplus1.violinplot.eps", width = 10, height = 5, dpi = 300, units = "in", device='eps')


ggplot(INTACT_TPM_all_genes[gene_class %in% c("MR genes", "MA genes")], aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("MR and MA only")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass + gene_class,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = median, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.MA.genes.log2TPMplus1.violinplot.thinner.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.MR.MA.genes.log2TPMplus1.violinplot.thinner.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


#minimum TPM 1 
ggplot(INTACT_TPM_all_genes[TPM_Value >=1], aes(x = Genotype, y = log2(as.numeric(TPM_Value) + 1), fill=Genotype))+
  ggtitle("TPM >= 1")+
  geom_violin(trim=TRUE)+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Gene expression, log2(TPM + 1)") + xlab("")+
  #separate by subclass
  facet_grid(.~subclass + gene_class,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  #trim=TRUE means the data is trimmed to fit the range of observation
  #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
  scale_fill_manual(name = "Genotype:", values = c("WT"="purple", "KO"="orange")) +
  #geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = median, 
               geom = "point", size=0.4, color="black") + 
  #stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.all.genes.TPMmin1.log2TPMplus1.violinplot.png", width = 10, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename = "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/violinplots/INTACT.subclass.WT.KO.all.genes.TPMmin1.log2TPMplus1.violinplot.eps", width = 10, height = 5, dpi = 300, units = "in", device='eps')
