library(data.table)
library(dplyr)
library(ggplot2)

Pv_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")
Sst_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_Mecp2KO_gene_TPMs_nondedup.txt")
L4_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_Mecp2KO_gene_TPMs_nondedup.txt")
L5_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_Mecp2KO_gene_TPMs_nondedup.txt")

#PV TPMs melted
Pv_TPM_melt <- unique(melt(Pv_TPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table

#coding genes
coding_genes_mm9 <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")
names(coding_genes_mm9) <- c("chrom", "start", "end", "gene", "num_transcripts", "strand")


PVGA_nonPromoter_cCREs_linked_genes=fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/PVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
nonPVGA_nonPromoter_cCREs_linked_genes=fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/nonPVGA_nonPromoter_cCRE_Cicero_linked_genes_mm9_genicBooleans.txt")
#split by intragenic and extragenic
PVGA_nonPromoter_cCREs_intragenicLinked_genes = PVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==1)]
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes = nonPVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==1)]

PVGA_nonPromoter_cCREs_extragenicLinked_genes = PVGA_nonPromoter_cCREs_linked_genes[(Intragenic==0) & (Intragenic_to_linked_gene==0)]
nonPVGA_nonPromoter_cCREs_extragenicLinked_genes = nonPVGA_nonPromoter_cCREs_linked_genes[(Intragenic==0) & (Intragenic_to_linked_gene==0)]


PVGA_nonPromoter_cCREs_intragenic_noncognateLinked_genes = PVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==0)]
nonPVGA_nonPromoter_cCREs_intragenic_noncognateLinked_genes = nonPVGA_nonPromoter_cCREs_linked_genes[(Intragenic==1) & (Intragenic_to_linked_gene==0)]


#join cCRE-gene linkage tables with PV TPM table
PVGA_nonPromoter_cCREs_linked_genes_PVTPM <- left_join(x=PVGA_nonPromoter_cCREs_linked_genes, y=Pv_TPM[, .(Gene, Pv_WT_TPM_avg, Pv_KO_TPM_avg)], by="Gene")
nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM <- left_join(x=nonPVGA_nonPromoter_cCREs_linked_genes, y=Pv_TPM[, .(Gene, Pv_WT_TPM_avg, Pv_KO_TPM_avg)], by="Gene")

#split by intragenic and extragenic
PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM = PVGA_nonPromoter_cCREs_linked_genes_PVTPM[(Intragenic==1) & (Intragenic_to_linked_gene==1)]
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM = nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM[(Intragenic==1) & (Intragenic_to_linked_gene==1)]

PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM = PVGA_nonPromoter_cCREs_linked_genes_PVTPM[(Intragenic==0) & (Intragenic_to_linked_gene==0)]
nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM = nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM[(Intragenic==0) & (Intragenic_to_linked_gene==0)]


PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt <- unique(melt(PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt <- unique(melt(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt <- unique(melt(PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt <- unique(melt(nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table

#gene lists
pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")


##intragenic cCREs linked to their cognate genes
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

##extragenic cCREs linked to their cognate genes
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


###intragenic linked
PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = rbind(
  cbind(PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV unchanged genes\nof PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV unchanged genes\nof PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV unchanged genes\nof non-PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV unchanged genes\nof non-PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV MR genes\nof PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(PVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV MR genes\nof PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV MR genes\nof non-PV ccREs", cCRE_genic_loc="Intragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_intragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV MR genes\nof non-PV ccREs", cCRE_genic_loc="Intragenic")
)

names(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO)[3] = "TPM"

PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO %>% mutate(gene_class = factor(gene_class, levels=c("PV unchanged genes\nof PV ccREs", "PV unchanged genes\nof non-PV ccREs", "PV MR genes\nof PV ccREs", "PV MR genes\nof non-PV ccREs")))
PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO, aes(x = gene_class, y = log2(as.numeric(TPM)+1), fill=gene_class))+
  ggtitle("Genes of intragenic linked cCREs")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "Genotype", values=c("WT"="purple", "KO"="orange"))+
  scale_fill_manual(name = "", values = c("PV unchanged genes\nof PV ccREs"="gray", "PV unchanged genes\nof non-PV ccREs"="gray", "PV MR genes\nof PV ccREs"="red", "PV MR genes\nof non-PV ccREs"="red"))+
  coord_cartesian(ylim=c(0,10.5))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  facet_grid(.~genotype,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MR_and_unchanged_genes_intragenicLinked_to_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MR_and_unchanged_genes_intragenicLinked_to_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


options(scipen = 0)

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV unchanged genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV unchanged genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 9.658813e-30, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV MR genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV MR genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.2592805, ns
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV unchanged genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV unchanged genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 9.274982e-30, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV MR genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_intragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV MR genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.2508031, ****

##
PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = rbind(
  cbind(PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV unchanged genes\nof PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV unchanged genes\nof PV ccREs", cCRE_genic_loc="Exragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV unchanged genes\nof non-PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_unchanged_genes_p0.5_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV unchanged genes\nof non-PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV MR genes\nof PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(PVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV MR genes\nof PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_class="PV MR genes\nof non-PV ccREs", cCRE_genic_loc="Extragenic"),
  cbind(nonPVGA_nonPromoter_cCREs_extragenicLinked_genes_PVTPM_melt[(Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_class="PV MR genes\nof non-PV ccREs", cCRE_genic_loc="Extragenic")
)

names(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO)[3] = "TPM"

PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO %>% mutate(gene_class = factor(gene_class, levels=c("PV unchanged genes\nof PV ccREs", "PV unchanged genes\nof non-PV ccREs", "PV MR genes\nof PV ccREs", "PV MR genes\nof non-PV ccREs")))
PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO = PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))

ggplot(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO, aes(x = gene_class, y = log2(as.numeric(TPM)+1), fill=gene_class))+
  ggtitle("Genes of extragenic linked cCREs")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "Genotype", values=c("WT"="purple", "KO"="orange"))+
  scale_fill_manual(name = "", values = c("PV unchanged genes\nof PV ccREs"="gray", "PV unchanged genes\nof non-PV ccREs"="gray", "PV MR genes\nof PV ccREs"="red", "PV MR genes\nof non-PV ccREs"="red"))+
  coord_cartesian(ylim=c(0,10.5))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  facet_grid(.~genotype,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MR_and_unchanged_genes_extragenicLinked_to_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MR_and_unchanged_genes_extragenicLinked_to_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


options(scipen = 0)

wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV unchanged genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV unchanged genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 2.537803e-37, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV MR genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="WT") & (gene_class=="PV MR genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.04235916, *
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV unchanged genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV unchanged genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 3.764909e-37, ****
wilcox.test(PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV MR genes\nof PV ccREs"), log2(as.numeric(TPM)+1)], PVGA_and_nonPVGA_nonPromoter_cCRE_extragenicLinked_to_Pv_MR_genes_noOther_nondedup_PVTPM_WT_KO[(genotype=="KO") & (gene_class=="PV MR genes\nof non-PV ccREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.04215707, *

###
unique(PVGA_nonPromoter_cCREs_linked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene])
unique(nonPVGA_nonPromoter_cCREs_linked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene])


length(intersect(PVGA_nonPromoter_cCREs_intragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene], nonPVGA_nonPromoter_cCREs_intragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene]))/length(unique(PVGA_nonPromoter_cCREs_intragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene]))

length(intersect(PVGA_nonPromoter_cCREs_extragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene], nonPVGA_nonPromoter_cCREs_extragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene]))/length(unique(PVGA_nonPromoter_cCREs_extragenicLinked_genes[Gene %in% pv_mr_genes_q0.1_nondedup_mm9$V4, Gene]))
##

##
mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")

mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")

mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/mousebrain_union_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")


PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")

PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")

PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")


nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_mm9.bed")

nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMA_mm9.bed")

nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/Pv_H3K27ac_chip_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged_mm9.bed")




#MR PV and non-PV cCREs linked to genes, PV TPMs
MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM = PVGA_nonPromoter_cCREs_linked_genes_PVTPM[(cCRE_label %in% PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR$V4), ]
MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM = nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR$V4), ]

unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM = PVGA_nonPromoter_cCREs_linked_genes_PVTPM[(cCRE_label %in% PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged$V4), ]
unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM = nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM[(cCRE_label %in% nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombUnchanged$V4), ]



#melt
MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt <- unique(melt(MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt <- unique(melt(MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt <- unique(melt(unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table
unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt <- unique(melt(unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM, id.vars="Gene", measure.vars=c("Pv_WT_TPM_avg", "Pv_KO_TPM_avg"))) %>% data.table

genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO = rbind(
  cbind(unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes of unchanged\nPV cCREs"),
  cbind(unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes of unchanged\nnon-PV cCREs"),
  cbind(MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes of MR\nPV cCREs"),
  cbind(MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes of MR\nnon-PV cCREs"),
  cbind(unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes of unchanged\nPV cCREs"),
  cbind(unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes of unchanged\nnon-PV cCREs"),
  cbind(MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes of MR\nPV cCREs"),
  cbind(MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt[(variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes of MR\nnon-PV cCREs")
)

names(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO)[3] = "TPM"

genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO = genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("Genes of unchanged\nPV cCREs", "Genes of unchanged\nnon-PV cCREs", "Genes of MR\nPV cCREs", "Genes of MR\nnon-PV cCREs")))
genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO = genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))

ggplot(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO, aes(x = gene_cCRE_group, y = log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("Genes linked to all MR and unchanged cCREs")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "Genotype", values=c("WT"="purple", "KO"="orange"))+
  scale_fill_manual(name = "", values = c("Genes of unchanged\nPV cCREs"="gray", "Genes of unchanged\nnon-PV cCREs"="gray", "Genes of MR\nPV cCREs"="red", "Genes of MR\nnon-PV cCREs"="red"))+
  coord_cartesian(ylim=c(0,12))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  facet_grid(.~genotype,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/genes_linked_to_MR_and_unchanged_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/genes_linked_to_MR_and_unchanged_PVGA_and_nonPVGA_nonPromoter_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')



length(intersect(MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene, MR_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene))/length(unique(MR_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene))
length(intersect(unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene, unchanged_nonPVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene))/length(unique(unchanged_PVGA_nonPromoter_cCREs_linked_genes_PVTPM_melt$Gene))

wilcox.test(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="WT") & (gene_cCRE_group=="Genes of unchanged\nPV cCREs"), log2(as.numeric(TPM)+1)], genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="WT") & (gene_cCRE_group=="Genes of unchanged\nnon-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 3.247674e-164, ****
wilcox.test(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="WT") & (gene_cCRE_group=="Genes of MR\nPV cCREs"), log2(as.numeric(TPM)+1)], genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="WT") & (gene_cCRE_group=="Genes of MR\nnon-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.1487358, ns
wilcox.test(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="KO") & (gene_cCRE_group=="Genes of unchanged\nPV cCREs"), log2(as.numeric(TPM)+1)], genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="KO") & (gene_cCRE_group=="Genes of unchanged\nnon-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 1.458822e-165, ****
wilcox.test(genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="KO") & (gene_cCRE_group=="Genes of MR\nPV cCREs"), log2(as.numeric(TPM)+1)], genes_linked_to_PVGA_and_nonPVGA_nonPromoter_MR_cCREs_PVTPM_WT_KO[(genotype=="KO") & (gene_cCRE_group=="Genes of MR\nnon-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 0.1142921, ns

###
#genes not linked to PV cCREs or non-PV cCREs
genes_unlinked_to_nonPromoter_cCREs <- setdiff(coding_genes_mm9$gene, c(PVGA_nonPromoter_cCREs_linked_genes[, Gene], nonPVGA_nonPromoter_cCREs_linked_genes[, Gene]))

PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs <- rbind(
 cbind(Pv_TPM_melt[(Gene %in% genes_unlinked_to_nonPromoter_cCREs) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes unlinked\nto cCREs"),
 cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_intragenicLinked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes intragenic linked\nto PV cCREs"),
 cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_intragenicLinked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes intragenic linked\nto non-PV cCREs"),
 cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_extragenicLinked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes extragenic linked\nto PV cCREs"),
 cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_extragenicLinked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes extragenic linked\nto non-PV cCREs")
)

names(PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs)[3] = "TPM"

PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs <- PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("Genes unlinked\nto cCREs", "Genes intragenic linked\nto PV cCREs", "Genes intragenic linked\nto non-PV cCREs", "Genes extragenic linked\nto PV cCREs", "Genes extragenic linked\nto non-PV cCREs")))


ggplot(PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs, aes(x=gene_cCRE_group, y=log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("PV MeCP2 WT")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "Genotype", values=c("WT"="purple", "KO"="orange"))+
  scale_fill_manual(name = "", values = c("Genes unlinked\nto cCREs"="gray", "Genes intragenic linked\nto PV cCREs"="red", "Genes intragenic linked\nto non-PV cCREs"="red", "Genes extragenic linked\nto PV cCREs"="red", "Genes extragenic linked\nto non-PV cCREs"="red"))+
  coord_cartesian(ylim=c(0,10.5))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_intra_and_extra_linked_to_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

###KO
PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs <- rbind(
  cbind(Pv_TPM_melt[(Gene %in% genes_unlinked_to_nonPromoter_cCREs) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes unlinked\nto cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_intragenicLinked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes intragenic linked\nto PV cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_intragenicLinked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes intragenic linked\nto non-PV cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_extragenicLinked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes extragenic linked\nto PV cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_extragenicLinked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes extragenic linked\nto non-PV cCREs")
)

names(PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs)[3] = "TPM"

PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs <- PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("Genes unlinked\nto cCREs", "Genes intragenic linked\nto PV cCREs", "Genes intragenic linked\nto non-PV cCREs", "Genes extragenic linked\nto PV cCREs", "Genes extragenic linked\nto non-PV cCREs")))


ggplot(PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs, aes(x=gene_cCRE_group, y=log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "Genotype", values=c("KO"="purple", "KO"="orange"))+
  scale_fill_manual(name = "", values = c("Genes unlinked\nto cCREs"="gray", "Genes intragenic linked\nto PV cCREs"="red", "Genes intragenic linked\nto non-PV cCREs"="red", "Genes extragenic linked\nto PV cCREs"="red", "Genes extragenic linked\nto non-PV cCREs"="red"))+
  coord_cartesian(ylim=c(0,10.5))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))

ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_KO_TPM_of_genes_intra_and_extra_linked_to_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

###expression of genes linked to all PV and non-PV cCREs, not split up by genic location

#WT
PV_WT_TPM_of_genes_linked_to_cCREs <- rbind(
  cbind(Pv_TPM_melt[(Gene %in% coding_genes_mm9$gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="All genes"),
  cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_linked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes linked\nto PV cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_linked_genes$Gene) & (variable=="Pv_WT_TPM_avg")], genotype="WT", gene_cCRE_group="Genes linked\nto non-PV cCREs")
)

names(PV_WT_TPM_of_genes_linked_to_cCREs)[3] = "TPM"

PV_WT_TPM_of_genes_linked_to_cCREs <- PV_WT_TPM_of_genes_linked_to_cCREs %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("All genes", "Genes linked\nto PV cCREs", "Genes linked\nto non-PV cCREs")))


ggplot(PV_WT_TPM_of_genes_linked_to_cCREs, aes(x=gene_cCRE_group, y=log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("PV MeCP2 WT")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("All genes"="gray", "Genes linked\nto PV cCREs"="red", "Genes linked\nto non-PV cCREs"="red"))+
  coord_cartesian(ylim=c(0,13))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="All genes"), log2(as.numeric(TPM)+1)], PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="Genes linked\nto PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p < 2.2e-16, ****
wilcox.test(PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="Genes linked\nto PV cCREs"), log2(as.numeric(TPM)+1)], PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="Genes linked\nto non-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 2.405899e-164, ****
wilcox.test(PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="All genes"), log2(as.numeric(TPM)+1)], PV_WT_TPM_of_genes_linked_to_cCREs[(genotype=="WT") & (gene_cCRE_group=="Genes linked\nto non-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p < 2.2e-16, ****
#
#KO
PV_KO_TPM_of_genes_linked_to_cCREs <- rbind(
  cbind(Pv_TPM_melt[(Gene %in% coding_genes_mm9$gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="All genes"),
  cbind(Pv_TPM_melt[(Gene %in% PVGA_nonPromoter_cCREs_linked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes linked\nto PV cCREs"),
  cbind(Pv_TPM_melt[(Gene %in% nonPVGA_nonPromoter_cCREs_linked_genes$Gene) & (variable=="Pv_KO_TPM_avg")], genotype="KO", gene_cCRE_group="Genes linked\nto non-PV cCREs")
)

names(PV_KO_TPM_of_genes_linked_to_cCREs)[3] = "TPM"

PV_KO_TPM_of_genes_linked_to_cCREs <- PV_KO_TPM_of_genes_linked_to_cCREs %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("All genes", "Genes linked\nto PV cCREs", "Genes linked\nto non-PV cCREs")))


ggplot(PV_KO_TPM_of_genes_linked_to_cCREs, aes(x=gene_cCRE_group, y=log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("PV MeCP2 KO")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("All genes"="gray", "Genes linked\nto PV cCREs"="red", "Genes linked\nto non-PV cCREs"="red"))+
  coord_cartesian(ylim=c(0,13))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_KO_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_KO_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="All genes"), log2(as.numeric(TPM)+1)], PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="Genes linked\nto PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p < 2.2e-16, ****
wilcox.test(PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="Genes linked\nto PV cCREs"), log2(as.numeric(TPM)+1)], PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="Genes linked\nto non-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p = 5.667932e-166, ****
wilcox.test(PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="All genes"), log2(as.numeric(TPM)+1)], PV_KO_TPM_of_genes_linked_to_cCREs[(genotype=="KO") & (gene_cCRE_group=="Genes linked\nto non-PV cCREs"), log2(as.numeric(TPM)+1)])$p.value #p < 2.2e-16, ****

##
#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")

PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs <- rbind(
  cbind(PV_deseq2[(Gene %in% coding_genes_mm9$gene)], gene_cCRE_group="All genes"),
  cbind(PV_deseq2[(Gene %in% PVGA_nonPromoter_cCREs_linked_genes$Gene)], gene_cCRE_group="Genes linked\nto PV cCREs"),
  cbind(PV_deseq2[(Gene %in% nonPVGA_nonPromoter_cCREs_linked_genes$Gene)], gene_cCRE_group="Genes linked\nto non-PV cCREs")
)

names(PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs)[3] = "TPM"

PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs <- PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs %>% mutate(gene_cCRE_group = factor(gene_cCRE_group, levels=c("All genes", "Genes linked\nto PV cCREs", "Genes linked\nto non-PV cCREs")))


ggplot(PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs, aes(x=gene_cCRE_group, y=ashr_log2FoldChange, fill=gene_cCRE_group))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("All genes"="gray", "Genes linked\nto PV cCREs"="red", "Genes linked\nto non-PV cCREs"="red"))+
  coord_cartesian(ylim=c(-0.15,0.15))+
  ylab("Gene expression log2 fold change (KO/WT)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


wilcox.test(PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="All genes"), ashr_log2FoldChange], PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="Genes linked\nto PV cCREs"), ashr_log2FoldChange])$p.value #p = 0.002934, **
wilcox.test(PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="Genes linked\nto PV cCREs"), ashr_log2FoldChange], PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="Genes linked\nto non-PV cCREs"), ashr_log2FoldChange])$p.value #p = 0.05765979, ns
wilcox.test(PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="All genes"), ashr_log2FoldChange], PV_WT_KO_logFC_exp_of_genes_linked_to_cCREs[(gene_cCRE_group=="Genes linked\nto non-PV cCREs"), ashr_log2FoldChange])$p.value #p = 0.3835446, ns
#

ggplot(PV_WT_TPM_of_genes_linked_to_cCREs[gene_cCRE_group != "All genes"], aes(x=gene_cCRE_group, y=log2(as.numeric(TPM)+1), fill=gene_cCRE_group))+
  ggtitle("PV MeCP2 WT")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("Genes linked\nto PV cCREs"="aquamarine", "Genes linked\nto non-PV cCREs"="aquamarine"))+
  coord_cartesian(ylim=c(0,13))+
  ylab("Gene expression, log2(TPM+1)") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_text(size=12, angle=90), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot_newColor.png", width = 2.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_WT_TPM_of_genes_linked_to_cCREs_log2TPMplus1_boxplot_newColor.eps", width = 2.5, height = 5, dpi = 300, units = "in", device='eps')
