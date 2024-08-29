library(data.table)
library(dplyr)
library(ggplot2)
library("biomaRt")
#coding genes
coding_genes_mm9 = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')
names(coding_genes_mm9) = c("chrom", "start", "end", "Gene", "transcripts", "strand")
#DEseq outputs of each INTACT-isolated subclass
PV_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
SST_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L4_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv")
L5_deseq <- read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv")

PV_deseq2 <- data.table(PV_deseq, keep.rownames="Gene")
SST_deseq2 <- data.table(SST_deseq, keep.rownames="Gene")
L4_deseq2 <- data.table(L4_deseq, keep.rownames="Gene")
L5_deseq2 <- data.table(L5_deseq, keep.rownames="Gene")

#transcription factor gene expression
#HOMER motifs enriched in PV MR cCREs relative to unchanged cCREs
PV_MR_cCRE_motifs = unique(fread("HG_lab/Mati/GabelLab/HOMER_motif_analyses/Li_Ren_2021_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_bgFishercombUnchanged_mm10/knownResults.txt"))
names(PV_MR_cCRE_motifs) = c("HomerID", "Consensus", "P.value", "Log.P.value", "q.value", "num_target_motif", "percent_target_motif", "num_bg_motif", "percent_bg_motif")
PV_MR_cCRE_motifs_enriched = PV_MR_cCRE_motifs[P.value <= 1e-2, ]
#HOMER motifs enriched in PV MR cCREs relative to unchanged cCREs
nonPV_MR_cCRE_motifs = unique(fread("HG_lab/Mati/GabelLab/HOMER_motif_analyses/Li_Ren_2021_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_bgFishercombUnchanged_mm10/knownResults.txt"))
names(nonPV_MR_cCRE_motifs) = c("HomerID", "Consensus", "P.value", "Log.P.value", "q.value", "num_target_motif", "percent_target_motif", "num_bg_motif", "percent_bg_motif")
nonPV_MR_cCRE_motifs_enriched = nonPV_MR_cCRE_motifs[P.value <= 1e-2, ]


#homer_MR_PVGA_cCREs = unique(fread("HG_lab/Mati/GabelLab/HOMER_motif_analyses/Li_Ren_2021_cCREs/PVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_bgFishercombUnchanged_mm10/knownResults.txt", sep="\t"))
#homer_MR_nonPVGA_cCREs = unique(fread("HG_lab/Mati/GabelLab/HOMER_motif_analyses/Li_Ren_2021_cCREs/nonPVGA_nonPromoter_cCREs_H3K27ac_ChIP_fishercombMR_bgFishercombUnchanged_mm10/knownResults.txt", sep="\t"))


#HOMER matrix IDs mapped to TFs in CIS-BP database
homer_ref = fread("HG_lab/Mati/GabelLab/HOMER_motif_analyses/Homer_DB_motifs_TF_With_CISBP_info_merged_FinalRefTable.csv")
homer_ref$Gene = homer_ref$TF_Name


#homer_MR_PVGA_cCREs_expand <- inner_join(x=homer_MR_PVGA_cCREs, y=homer_ref, by=c("Motif Name"="HomerID"))
#homer_MR_nonPVGA_cCREs_expand <- inner_join(x=homer_MR_nonPVGA_cCREs, y=homer_ref, by=c("Motif Name"="HomerID"))




length(unique(PV_MR_cCRE_motifs_enriched_more$Gene))
length(unique(nonPV_MR_cCRE_motifs_enriched_more$Gene))

#low intersection of gene names between mm9 and this TF list
intersect(PV_deseq2$Gene, PV_MR_cCRE_motifs_enriched_more$Gene)

#
mouse_mm9 = useMart("ensembl", host="https://may2009.archive.ensembl.org/", dataset = "mmusculus_gene_ensembl")
#mm9 Ensembl IDs
#mm9_ensIDs = getBM(attributes=c("ensembl_gene_id", "mgi_symbol"), mart=mouse_mm9)
attributes_vec <- c("mgi_symbol", "external_gene_id", "chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id")
#attributes_vec <- c("mgi_symbol")
mm9_ensIDs = getBM(attributes=attributes_vec, mart=mouse_mm9) %>% data.table

chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X")
mm9_ensIDs_chr1toX <- mm9_ensIDs[chromosome_name %in% chroms]  



#http://jul2019.archive.ensembl.org/


#mm9_ensIDs = fread("HG_lab/Mati/GabelLab/genesets/mm9_may2009_ensemblIDs.txt")

#mm9 Ensembl IDs

mouse_mart_2019 = useMart("ensembl", host="https://jul2019.archive.ensembl.org/", dataset = "mmusculus_gene_ensembl")

mouse_mart_2019_ensIDs = getBM(attributes=c("mgi_symbol", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id"), mart=mouse_mart_2019) %>% data.table
mouse_mart_2019_ensIDs_chr1toX <- mouse_mart_2019_ensIDs[chromosome_name %in% chroms]  

mm9_mouse2019_join <-  inner_join(x=mm9_ensIDs_chr1toX, y=mouse_mart_2019_ensIDs_chr1toX, by=c("external_gene_id" = "external_gene_name"))



human_mart_hg38 = useMart("ensembl", host="http://useast.ensembl.org/", dataset = "hsapiens_gene_ensembl")
mm39_ensIDs = data.table(getBM(attributes=c("ensembl_gene_id"), mart=mouse_mart_mm39))




#jan19 marts for mouse and human
mouse_jan19 = useMart("ensembl", host="https://jul2019.archive.ensembl.org/", dataset = "mmusculus_gene_ensembl")
human_jan19 = useMart("ensembl", host="https://jul2019.archive.ensembl.org/", dataset = "hsapiens_gene_ensembl")
mouse_jan19_ensIDs = data.table(getBM(attributes=c("ensembl_gene_id"), mart=mouse_jan19))
human_jan19_ensIDs = data.table(getBM(attributes=c("ensembl_gene_id"), mart=human_jan19))


mouse_jan19_ensIDs_mgi = data.table(getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "external_gene_name"), mart=mouse_jan19))

convertHumanIDs <- function(x){
  
  human = useMart("ensembl", host="https://jul2019.archive.ensembl.org/", dataset="hsapiens_gene_ensembl")
  mouse = useMart("ensembl", host="https://jul2019.archive.ensembl.org/", dataset = "mmusculus_gene_ensembl")
  
  IDs = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x, mart = human, attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows=T)
  mousex <- unique(IDs[, 2])
  
  return(mousex)
}

#table of human ensembl IDs (column 1) converted to mouse ensembl IDs (column 2) and mgi symbols (column 3)
human_to_mouse_jan19_ensIDs = data.table(getLDS(attributes = c("ensembl_gene_id"), filters="ensembl_gene_id", values=human_jan19_ensIDs$ensembl_gene_id, mart = human_jan19, attributesL = c("ensembl_gene_id", "mgi_symbol", "external_gene_name"), martL = mouse_jan19, uniqueRows=T))
#data table of mouse and human ensembl IDs (column 1) and mouse ensembl ID equivalents (column 2) for purpose of merging with homer reference table
human_mouse_dt =  data.table(rbind(
  cbind(DBID=mouse_jan19_ensIDs_mgi$ensembl_gene_id, mouse_ensID=mouse_jan19_ensIDs_mgi$ensembl_gene_id, mgi_symbol=mouse_jan19_ensIDs_mgi$mgi_symbol, mouse_external_gene_name=mouse_jan19_ensIDs_mgi$external_gene_name),
  cbind(DBID=human_to_mouse_jan19_ensIDs$Gene.stable.ID, mouse_ensID=human_to_mouse_jan19_ensIDs$Gene.stable.ID.1, mgi_symbol=human_to_mouse_jan19_ensIDs$MGI.symbol, mouse_external_gene_name=human_to_mouse_jan19_ensIDs$Gene.name)
))
#homer table with mouse ensembl ID equivalent column added through merging
homer_ref2 = data.table(inner_join(x=homer_ref, y=human_mouse_dt, by=c("DBID")))

setdiff(homer_ref2$mgi_symbol, coding_genes_mm9$Gene)

setdiff(homer_ref2$mgi_symbol, coding_genes_mm9$Gene)


PV_MR_cCRE_motifs_enriched_more <- left_join(x=PV_MR_cCRE_motifs_enriched, y=homer_ref2, by="HomerID") %>% unique
nonPV_MR_cCRE_motifs_enriched_more <- left_join(x=nonPV_MR_cCRE_motifs_enriched, y=homer_ref2, by="HomerID") %>% unique

PV_MR_cCRE_motifs_enriched_TFs <- PV_MR_cCRE_motifs_enriched_more[mgi_symbol != "NA", mgi_symbol]
nonPV_MR_cCRE_motifs_enriched_TFs <- nonPV_MR_cCRE_motifs_enriched_more[mgi_symbol != "NA", mgi_symbol]
non_enriched_cCRE_TFs <- setdiff(coding_genes_mm9$Gene, c(PV_MR_cCRE_motifs_enriched_more$mgi_symbol, nonPV_MR_cCRE_motifs_enriched_more$mgi_symbol))


setdiff(PV_MR_cCRE_motifs_enriched_more$mgi_symbol, coding_genes_mm9$Gene)
setdiff(nonPV_MR_cCRE_motifs_enriched_more$mgi_symbol, coding_genes_mm9$Gene)

PV_deseq2_cCRE_TF_fc <- rbind(
  cbind(PV_deseq2[Gene %in% non_enriched_cCRE_TFs], gene_class = "All other genes"),
  cbind(PV_deseq2[Gene %in% PV_MR_cCRE_motifs_enriched_TFs], gene_class = "PV MR cCRE TFs"),
  cbind(PV_deseq2[Gene %in% nonPV_MR_cCRE_motifs_enriched_TFs], gene_class = "Non-PV MR cCRE TFs")
)

PV_deseq2_cCRE_TF_fc = PV_deseq2_cCRE_TF_fc %>% mutate(gene_class = factor(gene_class, levels=c("All other genes", "PV MR cCRE TFs", "Non-PV MR cCRE TFs")))



ggplot(PV_deseq2_cCRE_TF_fc, aes(x = gene_class, y = ashr_log2FoldChange, fill=gene_class))+
  ggtitle("")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.1,0.1))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "gene_class:", values = c("All other genes"="gray", "PV MR cCRE TFs"="red", "Non-PV MR cCRE TFs"="lightpink1")) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/MR_cCRE_enriched_TFs_PV_log2FoldChange.png", width = 3, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/MR_cCRE_enriched_TFs_PV_log2FoldChange.eps", width = 3, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PV_deseq2_cCRE_TF_fc[gene_class=="All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc[gene_class=="PV MR cCRE TFs", ashr_log2FoldChange])$p.value
wilcox.test(PV_deseq2_cCRE_TF_fc[gene_class=="All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc[gene_class=="Non-PV MR cCRE TFs", ashr_log2FoldChange])$p.value
wilcox.test(PV_deseq2_cCRE_TF_fc[gene_class=="PV MR cCRE TFs", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc[gene_class=="Non-PV MR cCRE TFs", ashr_log2FoldChange])$p.value


#PV TPMs
Pv_TPM = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")

#PV MR and MA genes
PV_MR_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_mr_genes_q0.1_nondedup_mm9_promoters.bed")$V4
PV_MA_genes <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_ma_genes_q0.1_nondedup_mm9_promoters.bed")$V4



#png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/L5_over_L4_WT_noShrink_gene5foldDE_coding_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=Pv_TPM[, log2(Pv_WT_TPM_avg + 1)], y=Pv_TPM[, log2(Pv_KO_TPM_avg + 1)], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="", xlab="WT log2(TPM + 1)", ylab="KO log2(TPM + 1)", pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=Pv_TPM[Gene %in% PV_MR_genes, log2(Pv_WT_TPM_avg + 1)], y=Pv_TPM[Gene %in% PV_MR_genes, log2(Pv_KO_TPM_avg + 1)], pch=16, cex=1, col="purple")
points(x=Pv_TPM[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, log2(Pv_WT_TPM_avg + 1)], y=Pv_TPM[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, log2(Pv_KO_TPM_avg + 1)], pch=16, cex=1, col="orange")
abline(a = 0, b = 1, col = "red", lty = 2)  # 'h' specifies the y-value for the horizontal line
#dev.off()


smoothScatter(x=Pv_TPM[, Pv_WT_TPM_avg], y=Pv_TPM[, Pv_KO_TPM_avg], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="", xlab="WT TPM", ylab="KO TPM", pch=16, cex=0.6, col="lightgray", asp=1)
par(new=TRUE)
points(x=Pv_TPM[Gene %in% PV_MR_genes, Pv_WT_TPM_avg], y=Pv_TPM[Gene %in% PV_MR_genes, Pv_KO_TPM_avg], pch=16, cex=1, col="purple")
points(x=Pv_TPM[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, Pv_WT_TPM_avg], y=Pv_TPM[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, Pv_KO_TPM_avg], pch=16, cex=1, col="orange")
abline(a = 0, b = 1, col = "red", lty = 2)



#edgeR pseudocounts?
#
#MA plot with PV MR genes and PV MR cCRE TFs highlighted
png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/MA_plots/PV_MeCP2KO_over_WT_noShrink_PV_MR_genes_MR_cCRE_TFs_MA_smoothscatterplot.png", width=2000, height=2000, res=300)
smoothScatter(x=PV_deseq2[, log2(noshrink_contrast_baseMean)], y=PV_deseq2[, noshrink_contrast_log2FoldChange], colramp = colorRampPalette(c("white","gray","black")), transformation = function(x) x^.5, main="PV, MeCP2 KO over WT MA plot, no shrink", xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=0.6, col="lightgray", xlim=c(0,15), ylim=c(-4, 5.6), asp=1)
par(new=TRUE)
points(x=PV_deseq2[Gene %in% PV_MR_genes, log2(noshrink_contrast_baseMean)], y=PV_deseq2[Gene %in% PV_MR_genes, noshrink_contrast_log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="#4daf4a")
points(x=PV_deseq2[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, log2(noshrink_contrast_baseMean)], y=PV_deseq2[Gene %in% PV_MR_cCRE_motifs_enriched_TFs, noshrink_contrast_log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="orange")
points(x=PV_deseq2[Gene %in% nonPV_MR_cCRE_motifs_enriched_TFs, log2(noshrink_contrast_baseMean)], y=PV_deseq2[Gene %in% nonPV_MR_cCRE_motifs_enriched_TFs, noshrink_contrast_log2FoldChange], xlab="Log2 mean of normalized counts", ylab="Log2 fold change", pch=16, cex=1, col="purple")
abline(h = 0, col = "gray32", lty = 2)  # 'h' specifies the y-value for the horizontal line
dev.off()


#genes that aren't PV MeCP2-regulated or TFs associated with MR cCREs
non_PV_MeCP2reg_or_TF_genes <- setdiff(coding_genes_mm9$Gene, c(PV_MR_genes, PV_MA_genes, PV_MR_cCRE_motifs_enriched_more$mgi_symbol, nonPV_MR_cCRE_motifs_enriched_more$mgi_symbol))


PV_deseq2_cCRE_TF_fc2 <- rbind(
  cbind(PV_deseq2[Gene %in% non_PV_MeCP2reg_or_TF_genes], gene_class = "All other genes"),
  cbind(PV_deseq2[Gene %in% PV_MR_genes], gene_class = "PV MR genes"),
  cbind(PV_deseq2[Gene %in% PV_MA_genes], gene_class = "PV MA genes"),
  cbind(PV_deseq2[Gene %in% PV_MR_cCRE_motifs_enriched_TFs], gene_class = "PV MR cCRE TFs"),
  cbind(PV_deseq2[Gene %in% nonPV_MR_cCRE_motifs_enriched_TFs], gene_class = "Non-PV MR cCRE TFs")
)

PV_deseq2_cCRE_TF_fc2 = PV_deseq2_cCRE_TF_fc2 %>% mutate(gene_class = factor(gene_class, levels=c("All other genes", "PV MR genes", "PV MA genes", "PV MR cCRE TFs", "Non-PV MR cCRE TFs")))


ggplot(PV_deseq2_cCRE_TF_fc2, aes(x = gene_class, y = ashr_log2FoldChange, fill=gene_class))+
  ggtitle("")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.5,0.8))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "gene_class:", values = c("All other genes"="gray", "PV MR genes"="red", "PV MA genes"="blue", "PV MR cCRE TFs"="red4", "Non-PV MR cCRE TFs"="pink")) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MeCP2reg_genes_MR_cCRE_enriched_TFs_PV_log2FoldChange.png", width = 3.5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/PV_MeCP2reg_genes_MR_cCRE_enriched_TFs_PV_log2FoldChange.eps", width = 3.5, height = 5, dpi = 300, units = "in", device='eps')

wilcox.test(PV_deseq2_cCRE_TF_fc2[gene_class == "All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc2[gene_class == "PV MR genes", ashr_log2FoldChange])$p.value #p=2.014057e-151, ****
wilcox.test(PV_deseq2_cCRE_TF_fc2[gene_class == "PV MR genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc2[gene_class == "PV MA genes", ashr_log2FoldChange])$p.value #p=1.290525e-70, ****
wilcox.test(PV_deseq2_cCRE_TF_fc2[gene_class == "All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc2[gene_class == "PV MA genes", ashr_log2FoldChange])$p.value #p=5.650609e-125, ****
wilcox.test(PV_deseq2_cCRE_TF_fc2[gene_class == "All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc2[gene_class == "PV MR cCRE TFs", ashr_log2FoldChange])$p.value #p=0.7943756, ns
wilcox.test(PV_deseq2_cCRE_TF_fc2[gene_class == "All other genes", ashr_log2FoldChange], PV_deseq2_cCRE_TF_fc2[gene_class == "Non-PV MR cCRE TFs", ashr_log2FoldChange])$p.value #p=0.3656825, ns


ggplot(PV_deseq2_cCRE_TF_fc2, aes(x = gene_class, y = ashr_log2FoldChange, fill=gene_class))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.5,0.8))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "gene_class:", values = c("All other genes"="gray", "PV MR genes"="red", "PV MA genes"="blue", "PV MR cCRE TFs"="red4", "Non-PV MR cCRE TFs"="pink")) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))