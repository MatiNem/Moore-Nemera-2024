library(data.table)
library(dplyr)
library(ggplot2)

tissue_meths_all=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/brain_tissue_summary_global_meth.tsv")

transpose(tissue_meths_all)

tissue_meths_all[nucleotides=="CA",] - tissue_meths_all[nucleotides=="Lambda",]

tissue_meths_all2 <- melt(tissue_meths_all, id.var="nucleotides") %>% data.table


tissue_meths_all3 <- left_join(x=tissue_meths_all2[nucleotides=="CA", .(variable, value)], y=tissue_meths_all2[nucleotides=="Lambda", .(variable, value)], by="variable", suffix=c(".CA", ".Lambda"))


tissue_meths_all3$mCA_corrected <- tissue_meths_all3[, value.CA - value.Lambda]

tissue_meths_all3$tissue <- "NA"
tissue_meths_all3$genotype <- "NA"

tissue_meths_all3[grep("HYPO", variable), tissue:="HYPO"]
tissue_meths_all3[grep("STR", variable), tissue:="STR"]
tissue_meths_all3[grep("CB", variable), tissue:="CB"]

tissue_meths_all3[grep("WT", variable), genotype := "WT"]
tissue_meths_all3[grep("KO", variable), genotype := "KO"]


tissue_meths_all3 = tissue_meths_all3 %>% mutate(tissue = factor(tissue, levels=c("CB", "HYPO", "STR")))
tissue_meths_all3 = tissue_meths_all3 %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))


#bigger
ggplot(tissue_meths_all3, aes(x = genotype, y = as.numeric(mCA_corrected), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.47, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.03)) +
  ylab("Global mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~tissue, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_global_mCA_backgroundSub_WT_KO_stripplot_bigger.png", width = 3.5, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_global_mCA_backgroundSub_WT_KO_stripplot_bigger.eps", width = 3.5, height = 6, dpi = 300, units = "in", device=cairo_ps)



wilcox.test(tissue_meths_all3[(tissue=="CB") & (genotype=="WT"), mCA_corrected], tissue_meths_all3[(tissue=="CB") & (genotype=="KO"), mCA_corrected])$p.value #p-value=1
wilcox.test(tissue_meths_all3[(tissue=="HYPO") & (genotype=="WT"), mCA_corrected], tissue_meths_all3[(tissue=="HYPO") & (genotype=="KO"), mCA_corrected])$p.value #p-value=1
wilcox.test(tissue_meths_all3[(tissue=="STR") & (genotype=="WT"), mCA_corrected], tissue_meths_all3[(tissue=="STR") & (genotype=="KO"), mCA_corrected])$p.value #p-value=0.3333

#without background subtracting, testing
ggplot(tissue_meths_all3, aes(x = genotype, y = as.numeric(value.CA), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.04)) +
  ylab("Global mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~tissue, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_global_mCA_backgroundSub_WT_KO_stripplot_bigger.png", width = 3, height = 6, dpi = 300, units = "in", device='png')
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_global_mCA_backgroundSub_WT_KO_stripplot_bigger.eps", width = 3, height = 6, dpi = 300, units = "in", device=cairo_ps)

#####
#table with gene expression fold changes of brain tissues
tissue_fc_table <- fread("HG_lab/Mati/GabelLab/genesets/Large_Data_table_Gabel_et_al_1.txt")
tissue_fc_table[1, CB_KO_FC]
tissue_fc_table[1, HPTH_KO_FC]
tissue_fc_table[1, STR_KO_FC]
###

coding_genes = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')

##methylation of whole genes, brain tissues
CB_WT_KO_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_CB_WT_KO_mCA_mm9.bed")
HYPO_WT_KO_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_HYPO_WT_KO_mCA_mm9.bed")
STR_WT_KO_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_STR_WT_KO_mCA_mm9.bed")

#names of columns with each tissue
cols_CB = names(tissue_meths_all)[grep("CB", names(tissue_meths_all))]
cols_HYPO = names(tissue_meths_all)[grep("HYPO", names(tissue_meths_all))]
cols_STR = names(tissue_meths_all)[grep("STR", names(tissue_meths_all))]

#average lambda values
CB_WT_KO_avg_lambda <- mean(tissue_meths_all3[variable %in% cols_CB, value.Lambda])
HYPO_WT_KO_avg_lambda <- mean(tissue_meths_all3[variable %in% cols_HYPO, value.Lambda])
STR_WT_KO_avg_lambda <- mean(tissue_meths_all3[variable %in% cols_STR, value.Lambda])

#average lambda-subtracted mCA
CB_WT_KO_avg_mCA_corrected <- mean(tissue_meths_all3[variable %in% cols_CB, mCA_corrected])
HYPO_WT_KO_avg_mCA_corrected <- mean(tissue_meths_all3[variable %in% cols_HYPO, mCA_corrected])
STR_WT_KO_avg_mCA_corrected <- mean(tissue_meths_all3[variable %in% cols_STR, mCA_corrected])


#standard error function
std <- function(x) sd(x)/sqrt(length(x))

CB_WT_KO_se_mCA_corrected <- std(tissue_meths_all3[variable %in% cols_CB, mCA_corrected])
HYPO_WT_KO_se_mCA_corrected <- std(tissue_meths_all3[variable %in% cols_HYPO, mCA_corrected])
STR_WT_KO_se_mCA_corrected <- std(tissue_meths_all3[variable %in% cols_STR, mCA_corrected])

#naming columns
names(CB_WT_KO_fullGene_mCA) <- c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(HYPO_WT_KO_fullGene_mCA) <- c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(STR_WT_KO_fullGene_mCA) <- c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads")

CB_WT_KO_fullGene_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
HYPO_WT_KO_fullGene_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
STR_WT_KO_fullGene_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

#subtracting background non-conversion rate
CB_WT_KO_fullGene_mCA$gene_methylation_corrected <- CB_WT_KO_fullGene_mCA$gene_methylation - CB_WT_KO_avg_lambda
HYPO_WT_KO_fullGene_mCA$gene_methylation_corrected <- HYPO_WT_KO_fullGene_mCA$gene_methylation - HYPO_WT_KO_avg_lambda
STR_WT_KO_fullGene_mCA$gene_methylation_corrected <- STR_WT_KO_fullGene_mCA$gene_methylation - STR_WT_KO_avg_lambda

#making negative methylations zero for biological relevance
CB_WT_KO_fullGene_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
HYPO_WT_KO_fullGene_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
STR_WT_KO_fullGene_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]



##determining longest, most high mCA/CA genes in the genome
CB_WT_KO_fullGene_mCA$gene_length <- CB_WT_KO_fullGene_mCA[, end-start]
HYPO_WT_KO_fullGene_mCA$gene_length <- HYPO_WT_KO_fullGene_mCA[, end-start]
STR_WT_KO_fullGene_mCA$gene_length <- STR_WT_KO_fullGene_mCA[, end-start]


CB_WT_KO_fullGene_mCA$meth_decile <- as.character(ntile(CB_WT_KO_fullGene_mCA$gene_methylation_corrected, 10))
HYPO_WT_KO_fullGene_mCA$meth_decile <- as.character(ntile(HYPO_WT_KO_fullGene_mCA$gene_methylation_corrected, 10))
STR_WT_KO_fullGene_mCA$meth_decile <- as.character(ntile(STR_WT_KO_fullGene_mCA$gene_methylation_corrected, 10))


CB_WT_KO_fullGene_mCA <- CB_WT_KO_fullGene_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
HYPO_WT_KO_fullGene_mCA <- HYPO_WT_KO_fullGene_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))
STR_WT_KO_fullGene_mCA <- STR_WT_KO_fullGene_mCA %>%  mutate(meth_decile = factor(meth_decile, levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)))



#filtering to longest, most highly methylated genes
CB_WT_KO_fullGene_mCA_filt <- CB_WT_KO_fullGene_mCA[(gene_length > 1e5) & (meth_decile ==10)]
HYPO_WT_KO_fullGene_mCA_filt <- HYPO_WT_KO_fullGene_mCA[(gene_length > 1e5) & (meth_decile ==10)]
STR_WT_KO_fullGene_mCA_filt <- STR_WT_KO_fullGene_mCA[(gene_length > 1e5) & (meth_decile ==10)]

#just grabbing the genes
CB_long_highmCA_genes <- CB_WT_KO_fullGene_mCA_filt[, Gene]
HYPO_long_highmCA_genes <- HYPO_WT_KO_fullGene_mCA_filt[, Gene]
STR_long_highmCA_genes <- STR_WT_KO_fullGene_mCA_filt[, Gene]



tissue_long_highmCA_genes_FC <- rbind(
  cbind(id.var="CB", long_highmCA_gene_logfc_mean=mean(tissue_fc_table[Gene %in% CB_long_highmCA_genes , CB_KO_FC]), long_highmCA_gene_logfc_se=std(tissue_fc_table[Gene %in% CB_long_highmCA_genes, CB_KO_FC]), mean_mCA_corrected=CB_WT_KO_avg_mCA_corrected, se_mCA_corrected=CB_WT_KO_se_mCA_corrected),
  cbind(id.var="HYPO", long_highmCA_gene_logfc_mean=mean(tissue_fc_table[Gene %in% HYPO_long_highmCA_genes , HPTH_KO_FC]), long_highmCA_gene_logfc_se=std(tissue_fc_table[Gene %in% HYPO_long_highmCA_genes, HPTH_KO_FC]), mean_mCA_corrected=HYPO_WT_KO_avg_mCA_corrected, se_mCA_corrected=HYPO_WT_KO_se_mCA_corrected),
  cbind(id.var="STR", long_highmCA_gene_logfc_mean=mean(tissue_fc_table[Gene %in% STR_long_highmCA_genes , STR_KO_FC]), long_highmCA_gene_logfc_se=std(tissue_fc_table[Gene %in% STR_long_highmCA_genes, STR_KO_FC]), mean_mCA_corrected=STR_WT_KO_avg_mCA_corrected, se_mCA_corrected=STR_WT_KO_se_mCA_corrected)
) %>% data.table


ggplot(tissue_long_highmCA_genes_FC, aes(x = as.numeric(mean_mCA_corrected), y = as.numeric(long_highmCA_gene_logfc_mean), color = id.var)) +
  ggtitle("Long, highly methylated genes")+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = as.numeric(long_highmCA_gene_logfc_mean) - as.numeric(long_highmCA_gene_logfc_se), ymax = as.numeric(long_highmCA_gene_logfc_mean) + as.numeric(long_highmCA_gene_logfc_se)), 
                width = 0.0025) +  # Error bars for standard error
  geom_errorbar(aes(xmin = as.numeric(mean_mCA_corrected) - as.numeric(se_mCA_corrected), xmax = as.numeric(mean_mCA_corrected) + as.numeric(se_mCA_corrected)), 
                width = 0.0025) +  # Error bars for standard error
  scale_color_manual(name = "", values = c("CB"="lightgoldenrod3", "HYPO"="gray50", "STR"="black")) +
  ylab("Fold-change MeCP2 KO/WT")+xlab("Mean region mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.03), ylim = c(0, 0.1)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_region_corrected_mCA_long_top10percentHighmCA_gene_MeCP2KOWT_foldChange_dotplot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/CB_HYPO_STR_region_corrected_mCA_long_top10percentHighmCA_gene_MeCP2KOWT_foldChange_dotplot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')


#methylation of whole genes by rep
#CB
CB_WT_506_3_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/CB_reps/ensgene_mm9_CB_WT_Gabel_6_506-3_CB_D706_D506GAATTCGTAT_TAAGATTA_S59_CA.bed")
CB_WT_598_1_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/CB_reps/ensgene_mm9_CB_WT_Gabel_19_598-1_CB_D705_D505_ATTCAGAAAT_CTTCGCCT_S101_CA.bed")
CB_KO_506_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/CB_reps/ensgene_mm9_CB_KO_Gabel_5_506-2_CB_D705_D505ATTCAGAAAT_CTTCGCCT_S58_CA.bed")
CB_KO_598_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/CB_reps/ensgene_mm9_CB_KO_Gabel_20_598-2_CB_D706_D506_GAATTCGTAT_TAAGATTA_S102_CA.bed")

#HYPO
HYPO_WT_506_3_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/HYPO_reps/ensgene_mm9_HYPO_WT_Gabel_4_506-3_HYPO_D704_D504GAGATTCCAT_TCAGAGCC_S57_CA.bed")
HYPO_WT_598_1_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/HYPO_reps/ensgene_mm9_HYPO_WT_Gabel_17_598-1_HYPO_D703_D503_CGCTCATTAT_AGGATAGG_S99_CA.bed")
HYPO_KO_506_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/HYPO_reps/ensgene_mm9_HYPO_KO_Gabel_3_506-2_HYPO_D703_D503CGCTCATTAT_AGGATAGG_S56_CA.bed")
HYPO_KO_598_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/HYPO_reps/ensgene_mm9_HYPO_KO_Gabel_18_598-2_HYPO_D704_D504_GAGATTCCAT_TCAGAGCC_S100_CA.bed")

#STR
STR_WT_506_3_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/STR_reps/ensgene_mm9_STR_WT_Gabel_2_506-3_STR_D702_D502TCCGGAGAAT_GCCTCTAT_S55_CA.bed")
STR_WT_598_1_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/STR_reps/ensgene_mm9_STR_WT_Gabel_15_598-1_STR_D701_D501_ATTACTCGAT_AGGCTATA_S97_CA.bed")
STR_KO_506_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/STR_reps/ensgene_mm9_STR_KO_Gabel_1_506-2_STR_D701_D501ATTACTCGAT_AGGCTATA_S54_CA.bed")
STR_KO_598_2_fullGene_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/STR_reps/ensgene_mm9_STR_KO_Gabel_16_598-2_STR_D702_D502_TCCGGAGAAT_GCCTCTAT_S98_CA.bed")




#naming methylation table columns and calculating methylation levels, with decile calculation
tissue_meth_calc_func <- function(meth_table, tissue_global_meth_table, sample_name){
  meth_table$gene_methylation <- as.integer(meth_table[[7]])/as.integer(meth_table[[8]])
  meth_table$gene_methylation_corrected <- meth_table$gene_methylation - tissue_global_meth_table[variable==sample_name, value.Lambda]
  names(meth_table) = c("chrom", "start", "end", "Gene", "transcripts", "strand", "meth_reads", "cyto_reads", "gene_methylation", "gene_methylation_corrected")
  meth_table[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
  meth_table$meth_decile <- as.character(ntile(meth_table$gene_methylation_corrected, 10))
  return(meth_table)
}




#mCA
CB_WT_506_3_fullGene_mCA <- tissue_meth_calc_func(meth_table=CB_WT_506_3_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="CB_WT_Gabel_6_506-3_CB_D706_D506GAATTCGTAT_TAAGATTA_S59_summary")
CB_WT_598_1_fullGene_mCA <- tissue_meth_calc_func(meth_table=CB_WT_598_1_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="CB_WT_Gabel_19_598-1_CB_D705_D505_ATTCAGAAAT_CTTCGCCT_S101_summary")
CB_KO_506_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=CB_KO_506_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="CB_KO_Gabel_5_506-2_CB_D705_D505ATTCAGAAAT_CTTCGCCT_S58_summary")
CB_KO_598_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=CB_KO_598_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="CB_KO_Gabel_20_598-2_CB_D706_D506_GAATTCGTAT_TAAGATTA_S102_summary")

HYPO_WT_506_3_fullGene_mCA <- tissue_meth_calc_func(meth_table=HYPO_WT_506_3_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="HYPO_WT_Gabel_4_506-3_HYPO_D704_D504GAGATTCCAT_TCAGAGCC_S57_summary")
HYPO_WT_598_1_fullGene_mCA <- tissue_meth_calc_func(meth_table=HYPO_WT_598_1_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="HYPO_WT_Gabel_17_598-1_HYPO_D703_D503_CGCTCATTAT_AGGATAGG_S99_summary")
HYPO_KO_506_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=HYPO_KO_506_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="HYPO_KO_Gabel_3_506-2_HYPO_D703_D503CGCTCATTAT_AGGATAGG_S56_summary")
HYPO_KO_598_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=HYPO_KO_598_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="HYPO_KO_Gabel_18_598-2_HYPO_D704_D504_GAGATTCCAT_TCAGAGCC_S100_summary")

STR_WT_506_3_fullGene_mCA <- tissue_meth_calc_func(meth_table=STR_WT_506_3_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="STR_WT_Gabel_2_506-3_STR_D702_D502TCCGGAGAAT_GCCTCTAT_S55_summary")
STR_WT_598_1_fullGene_mCA <- tissue_meth_calc_func(meth_table=STR_WT_598_1_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="STR_WT_Gabel_15_598-1_STR_D701_D501_ATTACTCGAT_AGGCTATA_S97_summary")
STR_KO_506_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=STR_KO_506_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="STR_KO_Gabel_1_506-2_STR_D701_D501ATTACTCGAT_AGGCTATA_S54_summary")
STR_KO_598_2_fullGene_mCA <- tissue_meth_calc_func(meth_table=STR_KO_598_2_fullGene_mCA, tissue_global_meth_table=tissue_meths_all3, sample_name="STR_KO_Gabel_16_598-2_STR_D702_D502_TCCGGAGAAT_GCCTCTAT_S98_summary")




reps_tissues_fullGene_mCA <- rbind(
  cbind(mean_gene_meth=CB_WT_506_3_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="CB_WT_rep1", id.var="CB"),
  cbind(mean_gene_meth=CB_WT_598_1_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="CB_WT_rep1", id.var="CB"),
  cbind(mean_gene_meth=CB_KO_506_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="CB_KO_rep1", id.var="CB"),
  cbind(mean_gene_meth=CB_KO_598_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="CB_KO_rep2", id.var="CB"),
  cbind(mean_gene_meth=HYPO_WT_506_3_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="HYPO_WT_rep1", id.var="HYPO"),
  cbind(mean_gene_meth=HYPO_WT_598_1_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="HYPO_WT_rep1", id.var="HYPO"),
  cbind(mean_gene_meth=HYPO_KO_506_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="HYPO_KO_rep1", id.var="HYPO"),
  cbind(mean_gene_meth=HYPO_KO_598_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="HYPO_KO_rep2", id.var="HYPO"),
  cbind(mean_gene_meth=STR_WT_506_3_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="STR_WT_rep1", id.var="STR"),
  cbind(mean_gene_meth=STR_WT_598_1_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="STR_WT_rep1", id.var="STR"),
  cbind(mean_gene_meth=STR_KO_506_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="STR_KO_rep1", id.var="STR"),
  cbind(mean_gene_meth=STR_KO_598_2_fullGene_mCA[, mean(gene_methylation_corrected, na.rm=TRUE)], label="STR_KO_rep2", id.var="STR")
) %>% data.table
  
#grand mean and standard error of genic mCG levels across reps
tissues_mean_reps_data_full_mCA <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_tissues_fullGene_mCA, FUN = mean)
tissues_se_reps_data_full_mCA <- aggregate(as.numeric(mean_gene_meth) ~ id.var, data = reps_tissues_fullGene_mCA, FUN = function(x) sd(x) / sqrt(length(x)))

tissues_reps_full_mCA <- merge(tissues_mean_reps_data_full_mCA, tissues_se_reps_data_full_mCA, by = "id.var", suffixes = c("_mean", "_se"))
names(tissues_reps_full_mCA) <- c("id.var", "gene_meth_grand_mean", "gene_meth_grand_se")


reps_tissue_long_highmCA_genes_FC <- left_join(x=tissue_long_highmCA_genes_FC, y=tissues_reps_full_mCA, by=c("id.var"))



ggplot(reps_tissue_long_highmCA_genes_FC, aes(x = as.numeric(gene_meth_grand_mean), y = as.numeric(long_highmCA_gene_logfc_mean), color = id.var)) +
  ggtitle("Long, highly methylated genes")+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = as.numeric(long_highmCA_gene_logfc_mean) - as.numeric(long_highmCA_gene_logfc_se), ymax = as.numeric(long_highmCA_gene_logfc_mean) + as.numeric(long_highmCA_gene_logfc_se),
                width = 0.0025)) +  # Error bars for standard error
  geom_errorbar(aes(xmin = as.numeric(gene_meth_grand_mean) - as.numeric(gene_meth_grand_se), xmax = as.numeric(gene_meth_grand_mean) + as.numeric(gene_meth_grand_se))) +
  scale_color_manual(name = "", values = c("CB"="lightgoldenrod3", "HYPO"="gray50", "STR"="black")) +
  ylab("Fold-change MeCP2 KO/WT")+xlab("Mean region mCA/CA")+
  coord_cartesian(xlim=c(0.0,0.03), ylim = c(0, 0.1)) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/repErrors_CB_HYPO_STR_region_corrected_mCA_long_top10percentHighmCA_gene_MeCP2KOWT_foldChange_dotplot.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/repErrors_CB_HYPO_STR_region_corrected_mCA_long_top10percentHighmCA_gene_MeCP2KOWT_foldChange_dotplot.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')
