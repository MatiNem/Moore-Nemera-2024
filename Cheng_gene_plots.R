library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
options(scipy=999)

gene_panel_classes = fread("HG_lab/Vizgen/Gene_lists_for_analysis/MERFISH_gene_panel_classes.txt")
sig_function <- function(x){
  if(x < 0.05){
    sigx <- "*"
  }
  else{
    sigx <- "ns"
  }
  if(x < 0.01){
    sigx <- "**"
  }
  if(x < 0.001){
    sigx <- "***"
  }
  if(x < 0.0001){
    sigx <- "****"
  }
  return(sigx)
}
deseq2_WT <- data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq2_nsd1_wt_cko_8wk_normalized_counts.txt"), keep.rownames="Gene")


deseq2_WT$wt_avg <- deseq2_WT[, (wt1+wt2+wt3+wt4+wt5+wt6)/6]

ensgene_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed", header=FALSE)$V4
L23_genes <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.all.genes.txt", header=FALSE)$V1
L23_genes_mm9 <- intersect(ensgene_mm9, L23_genes)
non_L23_genes <- setdiff(ensgene_mm9, L23_genes_mm9)

#Performs resampling and returns randomized sample of rownames with equivalent expression values - full_data - data to sample from, samp = sample to resample against (determines expression values to sample from and resample size). Samp format is a matrix or dataframe of expression values in the same format as Full data.. Full data is expected to be a matrix or dataframe of expression values, with rownames of gene IDs. Duplicate option tells script to look for duplicate rownames in the resampling, and repeats resampling to avoid duplicates (off by default, not recommended to use if data to sample is >20% of the "full data")
resample <- function(full_data, samp, dup=FALSE) {
  expr_range = 10;
  if(ncol(samp) == 1){
    samp_vals = samp
  }
  else{
    samp_vals = rowMeans(samp)
  }
  full_df = as.data.frame(full_data)
  if(ncol(full_df) == 1){
    full_vals = full_df
  }
  else{
    full_vals = rowMeans(full_df)
  }
  samp_names = rownames(samp)
  full_df$means = full_vals
  order_data = full_df[order(full_df$means),]
  order_data$index = 1:nrow(order_data)
  order_data = order_data[!(rownames(order_data) %in% samp_names),]
  #ind = match(samp,names(ranks))
  resamp = samp_names
  ind = unlist(lapply(lapply(samp_vals,">",order_data$means),sum))
  resamps = sample(c(-expr_range:expr_range,1),length(ind),replace = TRUE)
  #resamp_ind = lapply(ind, FUN = function(x) x+sample(-expr_range:-1,1:expr_range,1))
  resamp_ranks = ind + resamps
  resamp_ranks[resamp_ranks<1] = sample(1:10, sum(resamp_ranks<1))
  resamp_ranks[resamp_ranks>nrow(order_data)] = sample(nrow(order_data)-0:9, sum(resamp_ranks>nrow(order_data)))
  resamp = rownames(order_data)[resamp_ranks]
  run_tries = 0
  while(dup & sum(duplicated(resamp))>0 & run_tries < 100){
    resamps = sample(c(-expr_range:expr_range,1),length(ind),replace = TRUE)
    resamp_ranks = ind + resamps
    resamp_ranks[resamp_ranks<1] = sample(1:expr_range, sum(resamp_ranks<1))
    resamp_ranks[resamp_ranks>nrow(order_data)] = sample(nrow(order_data)-1:expr_range+1, sum(resamp_ranks>nrow(order_data)))
    resamp = rownames(order_data)[resamp_ranks]
    run_tries = run_tries+1
  }
  return(resamp)
}

deseq2_WT_L23_df <- data.frame(deseq2_WT[Gene %in% L23_genes_mm9, .(Gene, wt_avg)], row.names="Gene")
deseq2_WT_nonL23_df <- data.frame(deseq2_WT[Gene %in% non_L23_genes, .(Gene, wt_avg)], row.names="Gene")

deseq2_WT_L23_df$wt_avg2 <- deseq2_WT_L23_df$wt_avg
deseq2_WT_nonL23_df$wt_avg2 <- deseq2_WT_nonL23_df$wt_avg

#example of one resampling
resamp1 <- resample(deseq2_WT_nonL23_df, deseq2_WT_L23_df)

set.seed(1)
#set number of resamplings
num_resamp = 100
#call matrix that you will populate with resamplings
resamplings = matrix(nrow=nrow(deseq2_WT_L23_df),ncol=num_resamp)
#for loop to populate the resampling matrix
for(n in 1:num_resamp){
  resamplings[,n] = resample(deseq2_WT_nonL23_df, deseq2_WT_L23_df)
}
write.table(resamplings, file="HG_lab/Mati/GabelLab/genesets/Cheng2022/L23.resampled.genes.txt", quote=F, row.names=F, col.names=F, sep="\t")

resamplings=fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L23.resampled.genes.txt", header=FALSE)

L23_and_resamp_genes <- rbind(
  cbind(gene=L23_genes_mm9, id.var="L2/3"),
  cbind(gene=resamplings$V1, id.var="Resampled")
) %>% data.table

genes_cCRE_number <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_genes_number_of_nonPromoter_cCREs.bed")[, 1:7]
names(genes_cCRE_number) = c("chrom", "start", "end", "gene", "num_transcripts", "strand", "num_cCREs")

L23_and_resamp_genes <- left_join(L23_and_resamp_genes, genes_cCRE_number, by=c("gene"))
L23_and_resamp_genes$gene_length <- L23_and_resamp_genes[, end-start]


genes_mm9_mca <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2_GSM2757472_Cortex_BS_8wk_ca.bed")
names(genes_mm9_mca) = c("chrom", "start", "end", "gene", "num_transcripts", "strand", "meth", "cov")
genes_mm9_mca$meth_level_mca <- genes_mm9_mca[, meth/cov]

genes_mm9_mcg <- fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2_GSM2757472_Cortex_BS_8wk_cg.bed")
names(genes_mm9_mcg) = c("chrom", "start", "end", "gene", "num_transcripts", "strand", "meth", "cov")
genes_mm9_mcg$meth_level_mcg <- genes_mm9_mcg[, meth/cov]


L23_and_resamp_genes <- left_join(L23_and_resamp_genes, genes_mm9_mca[, .(chrom, start, end, gene, num_transcripts, strand, meth_level_mca)], by=c("chrom", "start", "end", "gene", "num_transcripts", "strand"))
L23_and_resamp_genes <- left_join(L23_and_resamp_genes, genes_mm9_mcg[, .(chrom, start, end, gene, num_transcripts, strand, meth_level_mcg)], by=c("chrom", "start", "end", "gene", "num_transcripts", "strand"))

L23_and_resamp_genes = L23_and_resamp_genes %>% mutate(id.var = factor(id.var, levels=c("Resampled", "L2/3")))
#gene length plot
ggplot(L23_and_resamp_genes, aes(x = id.var, y = as.numeric(log10(gene_length)), fill=id.var))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Log10 gene length") + xlab("")+
  scale_fill_manual(name = "Genes:", values = c("Resampled"="gray", "L2/3"="salmon"))+ 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_log10GeneLength_boxplot.png", width = 2.5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_log10GeneLength_boxplot.eps", width = 2.5, height = 5 , dpi = 300, units = "in", device='eps')

#enhancer number plot
ggplot(L23_and_resamp_genes, aes(x = id.var, y = as.numeric(num_cCREs), fill=id.var))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,210))+
  ylab("Number of intragenic cCREs") + xlab("")+
  scale_fill_manual(name = "Genes:", values = c("Resampled"="gray", "L2/3"="salmon"))+ 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_number_cCREs_boxplot.png", width = 2, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_number_cCREs_boxplot.eps", width = 2, height = 5 , dpi = 300, units = "in", device='eps')

#mCA/CA plot
ggplot(L23_and_resamp_genes, aes(x = id.var, y = as.numeric(meth_level_mca), fill=id.var))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,0.06))+
  ylab("mCA/CA") + xlab("")+
  scale_fill_manual(name = "Genes:", values = c("Resampled"="gray", "L2/3"="salmon"))+ 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_mCA_boxplot.png", width = 2, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_mCA_boxplot.eps", width = 2, height = 5 , dpi = 300, units = "in", device='eps')

#mCG/CG plot
ggplot(L23_and_resamp_genes, aes(x = id.var, y = as.numeric(meth_level_mcg), fill=id.var))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0.25,1.0))+
  ylab("mCG/CG") + xlab("")+
  scale_fill_manual(name = "Genes:", values = c("Resampled"="gray", "L2/3"="salmon"))+ 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_mCG_boxplot.png", width = 3.5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23_Cheng2022_mCG_boxplot.eps", width = 3.5, height = 5 , dpi = 300, units = "in", device='eps')

wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(log10(gene_length))], L23_and_resamp_genes[id.var=="L2/3", as.numeric(log10(gene_length))])$p.value
wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(num_cCREs)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(num_cCREs)])$p.value
wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(meth_level_mca)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(meth_level_mca)])$p.value
wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(meth_level_mcg)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(meth_level_mcg)])$p.value

sig_function(wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(log10(gene_length))], L23_and_resamp_genes[id.var=="L2/3", as.numeric(log10(gene_length))])$p.value)
sig_function(wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(num_cCREs)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(num_cCREs)])$p.value)
sig_function(wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(meth_level_mca)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(meth_level_mca)])$p.value)
sig_function(wilcox.test(L23_and_resamp_genes[id.var=="Resampled", as.numeric(meth_level_mcg)], L23_and_resamp_genes[id.var=="L2/3", as.numeric(meth_level_mcg)])$p.value)

##overlap analysis
#L2/3 A, B, and C genes
L23_genes_A <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.A.genes.txt", header=FALSE)$V1
L23_genes_B <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.B.genes.txt", header=FALSE)$V1
L23_genes_C <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.C.genes.txt", header=FALSE)$V1


L23_genes_mm9_A <- intersect(ensgene_mm9, L23_genes_A)
L23_genes_mm9_B <- intersect(ensgene_mm9, L23_genes_B)
L23_genes_mm9_C <- intersect(ensgene_mm9, L23_genes_C)

non_L23_genes_mm9_A <- setdiff(ensgene_mm9, L23_genes_mm9_A)
non_L23_genes_mm9_B <- setdiff(ensgene_mm9, L23_genes_mm9_B)
non_L23_genes_mm9_C <- setdiff(ensgene_mm9, L23_genes_mm9_C)

#MERFISH gene panel A/B/C genes
L23_all_merfish <- intersect(gene_panel_classes$Gene, L23_genes)
L23_A_merfish <- intersect(gene_panel_classes$Gene, L23_genes_A)
L23_B_merfish <- intersect(gene_panel_classes$Gene, L23_genes_B)
L23_C_merfish <- intersect(gene_panel_classes$Gene, L23_genes_C)

nonL23_A_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_A)
nonL23_B_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_B)
nonL23_C_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_C)

#meta MR genes
meta_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
names(meta_MR_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")
meta_MA_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MA_genes_geneColumn4_mm9.bed")
names(meta_MA_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")

#cortex MR genes
cortex_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/mr_fishercomb_mm9.1.bed")
names(cortex_MR_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")
cortex_MA_genes = fread("HG_lab/Mati/GabelLab/genesets/ma_fishercomb_mm9.1.bed")
names(cortex_MA_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")


#non-MR or non-MA gene lists
non_meta_MR_gene_list = setdiff(ensgene_mm9, meta_MR_genes$gene)
non_cortex_MR_gene_list = setdiff(ensgene_mm9, cortex_MR_genes$gene)
non_meta_MA_gene_list = setdiff(ensgene_mm9, meta_MA_genes$gene)
non_cortex_MA_gene_list = setdiff(ensgene_mm9, cortex_MA_genes$gene)

#calculating length of intersections
meta_MR_cortex_MR = length(intersect(meta_MR_genes$Gene, cortex_MR_genes$Gene))
meta_MR_non_cortex_MR = length(intersect(meta_MR_genes$Gene, non_cortex_MR_gene_list))
non_meta_MR_cortex_MR = length(intersect(non_meta_MR_gene_list, cortex_MR_genes$Gene))
non_meta_MR_non_cortex_MR = length(intersect(non_meta_MR_gene_list, non_cortex_MR_gene_list))

#building the contingency table of overlaps
contingency_table = data.frame(c(meta_MR_cortex_MR, meta_MR_non_cortex_MR),
                               c(non_meta_MR_cortex_MR, non_meta_MR_non_cortex_MR),
                               row.names=c("cortex_MR", "non_cortex_MR"))
names(contingency_table) = c("meta_MR", "non_meta_MR")
test_contingency_table = fisher.test(contingency_table)
#p-value
p_values_contingency_table =  test_contingency_table$p.value
#log2 Fisher odds ratio
log2_OR = log2(test_contingency_table$estimate)


gene_overlap_func <- function(full_geneList, geneList1, geneList2){
  non_geneList1 <- setdiff(full_geneList, geneList1)
  non_geneList2 <- setdiff(full_geneList, geneList2)
  
  #calculating length of intersections
  geneList1_geneList2 = length(intersect(geneList1, geneList2))
  geneList1_nongeneList2 = length(intersect(geneList1, non_geneList2))
  nongeneList1_geneList2 = length(intersect(non_geneList1, geneList2))
  nongeneList1_nongeneList2 = length(intersect(non_geneList1, non_geneList2))
  
  #building the contingency table of overlaps
  contingency_table = data.frame(c(geneList1_geneList2, geneList1_nongeneList2),
                                 c(nongeneList1_geneList2, nongeneList1_nongeneList2),
                                 row.names=c("geneList2", "nongeneList2"))
  names(contingency_table) = c("geneList1", "nongeneList1")
  test_contingency_table = fisher.test(contingency_table)
  #p-value
  p_values_contingency_table =  test_contingency_table$p.value
  #log2 Fisher odds ratio
  log2_OR = log2(test_contingency_table$estimate)
  return(list(p_values_contingency_table, log2_OR, contingency_table))
}

gene_overlap_func(ensgene_mm9, L23_genes_mm9_A, cortex_MR_genes$gene)
gene_overlap_func(ensgene_mm9, L23_genes_mm9_A, cortex_MA_genes$gene)
gene_overlap_func(ensgene_mm9, L23_genes_mm9_B, cortex_MR_genes$gene)
gene_overlap_func(ensgene_mm9, L23_genes_mm9_B, cortex_MA_genes$gene)
gene_overlap_func(ensgene_mm9, L23_genes_mm9_C, cortex_MR_genes$gene)
gene_overlap_func(ensgene_mm9, L23_genes_mm9_C, cortex_MA_genes$gene)

A_metaMR_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_A, meta_MR_genes$gene)
A_metaMA_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_A, meta_MA_genes$gene)
B_metaMR_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_B, meta_MR_genes$gene)
B_metaMA_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_B, meta_MA_genes$gene)
C_metaMR_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_C, meta_MR_genes$gene)
C_metaMA_overlap <- gene_overlap_func(ensgene_mm9, L23_genes_mm9_C, meta_MA_genes$gene)

overlap_pvals_matrix <- matrix(c(unlist(A_metaMR_overlap[1]), 
                                 unlist(A_metaMA_overlap[1]),
                                 unlist(B_metaMR_overlap[1]),
                                 unlist(B_metaMA_overlap[1]),
                                 unlist(C_metaMR_overlap[1]),
                                 unlist(C_metaMA_overlap[1])), 
                               nrow=2, ncol=3)
rownames(overlap_pvals_matrix) = c("MR", "MA")
colnames(overlap_pvals_matrix) = c("A", "B", "C")

overlap_log10pvals_matrix <- -log10(overlap_pvals_matrix)

overlap_log2OR_matrix <- matrix(c(unlist(A_metaMR_overlap[2]), 
                                 unlist(A_metaMA_overlap[2]),
                                 unlist(B_metaMR_overlap[2]),
                                 unlist(B_metaMA_overlap[2]),
                                 unlist(C_metaMR_overlap[2]),
                                 unlist(C_metaMA_overlap[2])), 
                               nrow=2, ncol=3)
rownames(overlap_log2OR_matrix) = c("MR", "MA")
colnames(overlap_log2OR_matrix) = c("A", "B", "C")


col_palette <- colorRampPalette(c("white", "red", "darkred"))(n = 299)
#col_breaks = c(seq(-3, -1,length=100),  # for blue
#               seq(-0.99,0.99,length=100), #for white
#               seq(1, 3, length=100)) #for red

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_genes_Cheng2022_metaMR_metaMA_fisher_overlap_pvals_heatmap.eps")
heatmap.2(overlap_log10pvals_matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
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


###
Vis_Liu2023 <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/Liu2023_regions/VIS.2.Liu2023.subtypeAgg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

types_of_interest = c("166_L2/3 IT CTX", "171_L2/3 IT CTX", "168_L2/3 IT CTX", "179_L4 IT CTX")
Vis_Liu2023_TOI <- Vis_Liu2023[subtype %in% types_of_interest]

Vis_Liu2023_TOI = Vis_Liu2023_TOI %>% mutate(subtype = factor(subtype, levels=types_of_interest))

CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
CTX_HIP_annot[, cl := as.character(cl)]
#colors subtypes, subclasses, neighborhoods for ggplots
name_vec_subtype = unique(CTX_HIP_annot$cluster_label)
col_vec_subtype = unique(CTX_HIP_annot$cluster_color)
colors_subtype = setNames(col_vec_subtype, name_vec_subtype)

ggplot(Vis_Liu2023_TOI[gene %in% L23_genes_mm9_A], aes(x = subtype, y = as.numeric(logFC), fill=subtype))+
  ggtitle("L23_A genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1, 1))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Type:", values = colors_subtype[types_of_interest]) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(angle=90), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.A.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 4, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.A.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 4, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(Vis_Liu2023_TOI[gene %in% L23_genes_mm9_B], aes(x = subtype, y = as.numeric(logFC), fill=subtype))+
  ggtitle("L23_B genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1, 1))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Type:", values = colors_subtype[types_of_interest]) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(angle=90), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.B.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 4, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.B.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 4, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(Vis_Liu2023_TOI[gene %in% L23_genes_mm9_C], aes(x = subtype, y = as.numeric(logFC), fill=subtype))+
  ggtitle("L23_C genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1, 1.8))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Type:", values = colors_subtype[types_of_interest]) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(angle=90), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.C.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 4, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.C.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 4, height = 5 , dpi = 300, units = "in", device='eps')

Vis_Liu2023_TOI_A <- rbind(
 cbind(Vis_Liu2023_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-A"),
 cbind(Vis_Liu2023_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "A"), 
 cbind(Vis_Liu2023_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-A"),
 cbind(Vis_Liu2023_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "A"), 
 cbind(Vis_Liu2023_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-A"),
 cbind(Vis_Liu2023_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "A"), 
 cbind(Vis_Liu2023_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-A"),
 cbind(Vis_Liu2023_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "A") 
) %>% data.table

Vis_Liu2023_TOI_B <- rbind(
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "B"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "B"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "B"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-B"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "B") 
) %>% data.table

Vis_Liu2023_TOI_C <- rbind(
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "C"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "C"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "C"), 
  cbind(Vis_Liu2023_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-C"),
  cbind(Vis_Liu2023_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "C") 
) %>% data.table


Vis_Liu2023_TOI_A = Vis_Liu2023_TOI_A %>% mutate(subtype = factor(subtype, levels=types_of_interest))
Vis_Liu2023_TOI_A = Vis_Liu2023_TOI_A %>% mutate(gene_list = factor(gene_list, levels=c("non-A", "A")))
Vis_Liu2023_TOI_B = Vis_Liu2023_TOI_B %>% mutate(subtype = factor(subtype, levels=types_of_interest))
Vis_Liu2023_TOI_B = Vis_Liu2023_TOI_B %>% mutate(gene_list = factor(gene_list, levels=c("non-B", "B")))
Vis_Liu2023_TOI_C = Vis_Liu2023_TOI_C %>% mutate(subtype = factor(subtype, levels=types_of_interest))
Vis_Liu2023_TOI_C = Vis_Liu2023_TOI_C %>% mutate(gene_list = factor(gene_list, levels=c("non-C", "C")))

ggplot(Vis_Liu2023_TOI_A, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_A genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","purple3")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonA.vs.A.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonA.vs.A.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(Vis_Liu2023_TOI_B, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_B genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","red")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonB.vs.B.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonB.vs.B.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(Vis_Liu2023_TOI_C, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_C genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","darkgreen")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonC.vs.C.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonC.vs.C.genes.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


wilcox.test(Vis_Liu2023_TOI_C[(gene_list=="non-C") & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)], Vis_Liu2023_TOI_C[(gene_list=="C") & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)])$p.value

#function for performing Wilcoxon rank-sum test between expression fold change of L2/3 layer-specific genes and expression fold change of non-L2/3 layer-specific genes
L23_vs_nonL23_sigs_func <- function(logfc_table, geneList1, geneList2){
  sig1 <- pval1 <- NA
  sig2 <- pval2 <- NA
  sig3 <- pval3 <- NA
  sig4 <- pval4 <- NA
  
  #number of genes that go into Wilcoxon rank-sum tests
  n.geneList1.166 <- length(logfc_table[(gene_list==geneList1) & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)])
  n.geneList2.166 <- length(logfc_table[(gene_list==geneList2) & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)])
  
  n.geneList1.171 <- length(logfc_table[(gene_list==geneList1) & (subtype=="171_L2/3 IT CTX"), as.numeric(logFC)])
  n.geneList2.171 <- length(logfc_table[(gene_list==geneList2) & (subtype=="171_L2/3 IT CTX"), as.numeric(logFC)])
  
  n.geneList1.168 <- length(logfc_table[(gene_list==geneList1) & (subtype=="168_L2/3 IT CTX"), as.numeric(logFC)])
  n.geneList2.168 <- length(logfc_table[(gene_list==geneList2) & (subtype=="168_L2/3 IT CTX"), as.numeric(logFC)])
  
  n.geneList1.179 <- length(logfc_table[(gene_list==geneList1) & (subtype=="179_L4 IT CTX"), as.numeric(logFC)])
  n.geneList2.179 <- length(logfc_table[(gene_list==geneList2) & (subtype=="179_L4 IT CTX"), as.numeric(logFC)])
  
  tryCatch({
    pval1 <- wilcox.test(logfc_table[(gene_list==geneList1) & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)], logfc_table[(gene_list==geneList2) & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)])$p.value
    sig1 <- sig_function(pval1)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval2 <- wilcox.test(logfc_table[(gene_list==geneList1) & (subtype=="171_L2/3 IT CTX"), as.numeric(logFC)], logfc_table[(gene_list==geneList2) & (subtype=="171_L2/3 IT CTX"), as.numeric(logFC)])$p.value
    sig2 <- sig_function(pval2)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval3 <- wilcox.test(logfc_table[(gene_list==geneList1) & (subtype=="168_L2/3 IT CTX"), as.numeric(logFC)], logfc_table[(gene_list==geneList2) & (subtype=="168_L2/3 IT CTX"), as.numeric(logFC)])$p.value
    sig3 <- sig_function(pval3)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval4 <- wilcox.test(logfc_table[(gene_list==geneList1) & (subtype=="179_L4 IT CTX"), as.numeric(logFC)], logfc_table[(gene_list==geneList2) & (subtype=="179_L4 IT CTX"), as.numeric(logFC)])$p.value
    sig4 <- sig_function(pval4)
  }, error = function(e) {
    # Handle the error
  })
  
  sigs <- data.table(rbind(cbind(subtype="166_L2/3 IT CTX", gene_list1 = geneList1, gene_list2 = geneList2, n_gene_list1=n.geneList1.166, n_gene_list2=n.geneList2.166, sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(subtype="171_L2/3 IT CTX", gene_list1 = geneList1, gene_list2 = geneList2, n_gene_list1=n.geneList1.171, n_gene_list2=n.geneList2.171, sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(subtype="168_L2/3 IT CTX", gene_list1 = geneList1, gene_list2 = geneList2, n_gene_list1=n.geneList1.168, n_gene_list2=n.geneList2.168, sig_symbol=sig3, wilcox.pval=pval3),
                           cbind(subtype="179_L4 IT CTX", gene_list1 = geneList1, gene_list2 = geneList2, n_gene_list1=n.geneList1.179, n_gene_list2=n.geneList2.179, sig_symbol=sig4, wilcox.pval=pval4)))
  return(sigs)
}

sigs_Vis_Liu2023_TOI_A <- L23_vs_nonL23_sigs_func(logfc_table=Vis_Liu2023_TOI_A, geneList1="A", geneList2="non-A")
sigs_Vis_Liu2023_TOI_B <- L23_vs_nonL23_sigs_func(logfc_table=Vis_Liu2023_TOI_B, geneList1="B", geneList2="non-B")
sigs_Vis_Liu2023_TOI_C <- L23_vs_nonL23_sigs_func(logfc_table=Vis_Liu2023_TOI_C, geneList1="C", geneList2="non-C")

sigs_Vis_Liu2023_TOI_dt <- rbind(cbind(sigs_Vis_Liu2023_TOI_A, pseudobulkDGE_nCells_min = 10),
                                  cbind(sigs_Vis_Liu2023_TOI_B, pseudobulkDGE_nCells_min = 10),
                                  cbind(sigs_Vis_Liu2023_TOI_C, pseudobulkDGE_nCells_min = 10)) %>% data.table


###overlaps between statistically significant MR genes and A/B/C genes

#subtype significantly MR gene overlap with L2/3 A/B/C genes
MR_L23_166_A <-  gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish)
MR_L23_166_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish)
MR_L23_166_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish)

MR_L23_171_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish)
MR_L23_171_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish)
MR_L23_171_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish)

MR_L23_168_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish)
MR_L23_168_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish)
MR_L23_168_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish)

MR_L4_179_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish)
MR_L4_179_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish)
MR_L4_179_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish)

#subtype significantly MA gene overlap with L2/3 A/B/C genes
MA_L23_166_A <-  gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish)
MA_L23_166_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish)
MA_L23_166_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish)

MA_L23_171_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish)
MA_L23_171_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish)
MA_L23_171_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish)

MA_L23_168_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish)
MA_L23_168_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish)
MA_L23_168_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish)

MA_L4_179_A <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish)
MA_L4_179_B <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish)
MA_L4_179_C <- gene_overlap_func(gene_panel_classes[,Gene], Vis_Liu2023_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish)
#matrix of p-values of overlaps between statistically Mecp2-regulated genes and L2/3 A/B/C genes
sig_ABC_pvals_matrix <- matrix(c(unlist(MR_L23_166_A[1]), unlist(MR_L23_166_B[1]), unlist(MR_L23_166_C[1]),
                                 unlist(MR_L23_171_A[1]), unlist(MR_L23_171_B[1]), unlist(MR_L23_171_C[1]),
                                 unlist(MR_L23_168_A[1]), unlist(MR_L23_168_B[1]), unlist(MR_L23_168_C[1]),
                                 unlist(MR_L4_179_A[1]), unlist(MR_L4_179_B[1]), unlist(MR_L4_179_C[1]),
                                 unlist(MA_L23_166_A[1]), unlist(MA_L23_166_B[1]), unlist(MA_L23_166_C[1]),
                                 unlist(MA_L23_171_A[1]), unlist(MA_L23_171_B[1]), unlist(MA_L23_171_C[1]),
                                 unlist(MA_L23_168_A[1]), unlist(MA_L23_168_B[1]), unlist(MA_L23_168_C[1]),
                                 unlist(MA_L4_179_A[1]), unlist(MA_L4_179_B[1]), unlist(MA_L4_179_C[1])), 
                               nrow=3, ncol=8)
rownames(sig_ABC_pvals_matrix) = c("A", "B", "C")
colnames(sig_ABC_pvals_matrix) = c("MR_L23_166", "MR_L23_171", "MR_L23_168", "MR_L4_179",
                                   "MA_L23_166", "MA_L23_171", "MA_L23_168", "MA_L4_179")

write.table(sig_ABC_pvals_matrix, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmap/L23_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_merfish_pvals.txt", row.names=TRUE, col.names=TRUE, quote=F, sep="\t")
sig_ABC_log10pvals_matrix <- -log10(sig_ABC_pvals_matrix + 1 - 1)

#log2 odds ratio
sig_ABC_log2OR_matrix <- matrix(c(unlist(MR_L23_166_A[2]), unlist(MR_L23_166_B[2]), unlist(MR_L23_166_C[2]),
                                 unlist(MR_L23_171_A[2]), unlist(MR_L23_171_B[2]), unlist(MR_L23_171_C[2]),
                                 unlist(MR_L23_168_A[2]), unlist(MR_L23_168_B[2]), unlist(MR_L23_168_C[2]),
                                 unlist(MR_L4_179_A[2]), unlist(MR_L4_179_B[2]), unlist(MR_L4_179_C[2]),
                                 unlist(MA_L23_166_A[2]), unlist(MA_L23_166_B[2]), unlist(MA_L23_166_C[2]),
                                 unlist(MA_L23_171_A[2]), unlist(MA_L23_171_B[2]), unlist(MA_L23_171_C[2]),
                                 unlist(MA_L23_168_A[2]), unlist(MA_L23_168_B[2]), unlist(MA_L23_168_C[2]),
                                 unlist(MA_L4_179_A[2]), unlist(MA_L4_179_B[2]), unlist(MA_L4_179_C[2])), 
                               nrow=3, ncol=8)
rownames(sig_ABC_log2OR_matrix) = c("A", "B", "C")
colnames(sig_ABC_log2OR_matrix) = c("MR_L23_166", "MR_L23_171", "MR_L23_168", "MR_L4_179",
                                   "MA_L23_166", "MA_L23_171", "MA_L23_168", "MA_L4_179")

round(sig_ABC_log2OR_matrix,2)
col_palette2 <- colorRampPalette(c("white", "red", "darkred"))(n = 299)
#col_breaks = c(seq(-3, -1,length=100),  # for blue
#               seq(-0.99,0.99,length=100), #for white
#               seq(1, 3, length=100)) #for red



setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_merfish_heatmap.eps")
heatmap.2(sig_ABC_log10pvals_matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=col_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, key.title="Log10 p-value")
dev.off()

####heatmap of average expression
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


WT_log_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
KO_log_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

WT_and_KO_log_avgExp_matrix_TOI <- left_join(x=WT_log_avgExp_matrix[, c("Gene","166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX")],
                                         y=KO_log_avgExp_matrix[, c("Gene","166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX")],
                                         by=c("Gene"), suffix = c(".WT", ".KO"))

#formatting L23 gene names to fit data.frame standards
L23_genes_newNames <- make.names(L23_genes)
L23_A_merfish_newNames <- make.names(L23_A_merfish)
L23_B_merfish_newNames <- make.names(L23_B_merfish)
L23_C_merfish_newNames <- make.names(L23_C_merfish)

WT_and_KO_log_avgExp_matrix_TOI_ABC <- data.frame(WT_and_KO_log_avgExp_matrix_TOI[Gene %in% L23_genes_newNames], row.names="Gene")

WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore <- t(scale(t(as.matrix(WT_and_KO_log_avgExp_matrix_TOI_ABC))))


my_palette <- colorRampPalette(c("blue", "white",  "red" ))(n = 299)


type_order <- c("X166_L2.3.IT.CTX.WT", "X166_L2.3.IT.CTX.KO", 
                "X171_L2.3.IT.CTX.WT", "X171_L2.3.IT.CTX.KO",
                "X168_L2.3.IT.CTX.WT", "X168_L2.3.IT.CTX.KO",
                "X179_L4.IT.CTX.WT", "X179_L4.IT.CTX.KO")

#create heat map with key
heatmap.2(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[make.names(L23_A_merfish), ], main = "A",
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.7, 
          cexCol = 0.7, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5)) 

heatmap.2(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[make.names(L23_B_merfish), ], main = "B",
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.7, 
          cexCol = 0.7, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5)) 
heatmap.2(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[make.names(L23_C_merfish), ], main = "C",
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.7, 
          cexCol = 0.7, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5)) 

#setEPS()
#postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_ABC_logAverageExpression_zscores_merfish_heatmap.eps")
#heatmap.2(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[c(make.names(L23_A_merfish), make.names(L23_B_merfish), make.names(L23_C_merfish)),],
#          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
#          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.2, 
#          cexCol = 0.6, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
#          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
#          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
#dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/Vis_L23_ABC_AverageExpression_zscores_merfish_heatmap.eps")
heatmap.2(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[c(make.names(L23_A_merfish), make.names(L23_B_merfish), make.names(L23_C_merfish)), type_order],
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.2, 
          cexCol = 0.6, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()
#WT_log_avgExp_dt <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
#KO_log_avgExp_dt <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

#WT_log_avgExp_dt$subtype.t.type <- WT_log_avgExp_dt[, paste0(subtype, "_WT")]
#KO_log_avgExp_dt$subtype.t.type <- KO_log_avgExp_dt[, paste0(subtype, "_KO")]
###
exp1_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp2_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp6_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp6.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp7_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp7.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")

exp1_sct_counts2 <- left_join(x=exp1_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp2_sct_counts2 <- left_join(x=exp2_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp6_sct_counts2 <- left_join(x=exp6_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp7_sct_counts2 <- left_join(x=exp7_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))


all_sct_counts <- rbind(exp1_sct_counts, exp2_sct_counts, exp6_sct_counts, exp7_sct_counts)
#visual cortex cells
All_Vis_layerdepth <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/layer_depth/All_Vis_layerdepth_fixed.csv")

exp1_sct_counts_Vis <- exp1_sct_counts2[index %in% All_Vis_layerdepth$cell.id,]
exp2_sct_counts_Vis <- exp2_sct_counts2[index %in% All_Vis_layerdepth$cell.id,]
exp6_sct_counts_Vis <- exp6_sct_counts2[index %in% All_Vis_layerdepth$cell.id,]
exp7_sct_counts_Vis <- exp7_sct_counts2[index %in% All_Vis_layerdepth$cell.id,]

exp1_sct_counts_Vis_glut <- exp1_sct_counts_Vis[class_label=="Glutamatergic",]
exp2_sct_counts_Vis_glut <- exp2_sct_counts_Vis[class_label=="Glutamatergic",]
exp6_sct_counts_Vis_glut <- exp6_sct_counts_Vis[class_label=="Glutamatergic",]
exp7_sct_counts_Vis_glut <- exp7_sct_counts_Vis[class_label=="Glutamatergic",]

exp1_sct_counts_Vis_glut_genesOnly <- data.frame(exp1_sct_counts_Vis_glut[, 1:551], row.names="index")
exp1_sct_counts_Vis_glut_scaled <- scale(exp1_sct_counts_Vis_glut_genesOnly)

exp2_sct_counts_Vis_glut_genesOnly <- data.frame(exp2_sct_counts_Vis_glut[, 1:551], row.names="index")
exp2_sct_counts_Vis_glut_scaled <- scale(exp2_sct_counts_Vis_glut_genesOnly)

exp6_sct_counts_Vis_glut_genesOnly <- data.frame(exp6_sct_counts_Vis_glut[, 1:551], row.names="index")
exp6_sct_counts_Vis_glut_scaled <- scale(exp6_sct_counts_Vis_glut_genesOnly)

exp7_sct_counts_Vis_glut_genesOnly <- data.frame(exp7_sct_counts_Vis_glut[, 1:551], row.names="index")
exp7_sct_counts_Vis_glut_scaled <- scale(exp7_sct_counts_Vis_glut_genesOnly)

#exp2_sct_counts_Vis_WT <- exp2_sct_counts_Vis[t.type.Under1Over1=="WT",]
#exp2_sct_counts_Vis_KO <- exp2_sct_counts_Vis[t.type.Under1Over1=="KO",]

all_sct_counts_Vis <- all_sct_counts[index %in% All_Vis_layerdepth$cell.id,]

all_sct_counts_Vis_genesOnly <- data.frame(all_sct_counts_Vis[, 1:551], row.names="index")

#calculate z-scores of SCTransform-corrected count
all_sct_counts_Vis_scaled <- scale(all_sct_counts_Vis_genesOnly)



all_sct_counts_Vis_scaled_A <- rowMeans(all_sct_counts_Vis_scaled[, L23_A_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(all_sct_counts_Vis_scaled_A)[2] = "avg_zScore"
#exp2_sct_counts_Vis_scaled_A <- inner_join(x=all_sct_counts_Vis_scaled_A, y=exp2_sct_counts_Vis[, .(index, t.type.Under1Over1, center_x, center_y)], by="index")

all_sct_counts_Vis_scaled_B <- rowMeans(all_sct_counts_Vis_scaled[, L23_B_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(all_sct_counts_Vis_scaled_B)[2] = "avg_zScore"
#exp2_sct_counts_Vis_scaled_B <- inner_join(x=all_sct_counts_Vis_scaled_B, y=exp2_sct_counts_Vis[, .(index, t.type.Under1Over1, center_x, center_y)], by="index")

all_sct_counts_Vis_scaled_C <- rowMeans(all_sct_counts_Vis_scaled[, L23_C_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(all_sct_counts_Vis_scaled_C)[2] = "avg_zScore"
#exp2_sct_counts_Vis_scaled_C <- inner_join(x=all_sct_counts_Vis_scaled_C, y=exp2_sct_counts_Vis[, .(index, t.type.Under1Over1, center_x, center_y)], by="index")


limx<-c(0,10000)
limy<-c(0,10000)

#vis glut
exp2_sct_counts_Vis_glut_scaled_A <- rowMeans(exp2_sct_counts_Vis_glut_scaled[, L23_A_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(exp2_sct_counts_Vis_glut_scaled_A)[2] = "avg_zScore"
exp2_sct_counts_Vis_glut_scaled_A  <- inner_join(x=exp2_sct_counts_Vis_glut_scaled_A, y=exp2_sct_counts_Vis_glut, by="index")

exp2_sct_counts_Vis_glut_scaled_B <- rowMeans(exp2_sct_counts_Vis_glut_scaled[, L23_B_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(exp2_sct_counts_Vis_glut_scaled_B)[2] = "avg_zScore"
exp2_sct_counts_Vis_glut_scaled_B  <- inner_join(x=exp2_sct_counts_Vis_glut_scaled_B, y=exp2_sct_counts_Vis_glut, by="index")

exp2_sct_counts_Vis_glut_scaled_C <- rowMeans(exp2_sct_counts_Vis_glut_scaled[, L23_C_merfish_newNames]) %>% data.frame %>% data.table(keep.rownames="index")
names(exp2_sct_counts_Vis_glut_scaled_C)[2] = "avg_zScore"
exp2_sct_counts_Vis_glut_scaled_C  <- inner_join(x=exp2_sct_counts_Vis_glut_scaled_C, y=exp2_sct_counts_Vis_glut, by="index")

A_pal <- colorRampPalette(c('white','blue'))
B_pal <- colorRampPalette(c('white','goldenrod4'))
C_pal <- colorRampPalette(c('white','green'))

exp2_sct_counts_Vis_glut_scaled_A$plot_color <- A_pal(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_scaled_A$avg_zScore, breaks=20))]
exp2_sct_counts_Vis_glut_scaled_B$plot_color <- B_pal(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_scaled_B$avg_zScore, breaks=20))]
exp2_sct_counts_Vis_glut_scaled_C$plot_color <- C_pal(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_scaled_C$avg_zScore, breaks=20))]

#function for averaging scaled SCTransform-corrected counts
exp_counts_scale_func <- function(exp_scaled, exp_counts_summary, gene_list, color_palette){
  exp_scaled_gene_list <- rowMeans(exp_scaled[, gene_list]) %>% data.frame %>% data.table(keep.rownames="index")
  names(exp_scaled_gene_list)[2] = "avg_zScore"
  exp_scaled_gene_list <- inner_join(x=exp_scaled_gene_list, y=exp_counts_summary, by="index")
  
  exp_scaled_gene_list$plot_color <- get(color_palette)(20)[as.numeric(cut(exp_scaled_gene_list$avg_zScore, breaks=20))]
  return(exp_scaled_gene_list)
}


exp1_sct_counts_Vis_glut_scaled_A = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut, gene_list=L23_A_merfish_newNames, color_palette="A_pal")
exp1_sct_counts_Vis_glut_scaled_B = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut, gene_list=L23_B_merfish_newNames, color_palette="B_pal")
exp1_sct_counts_Vis_glut_scaled_C = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut, gene_list=L23_C_merfish_newNames, color_palette="C_pal")

exp2_sct_counts_Vis_glut_scaled_A = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut, gene_list=L23_A_merfish_newNames, color_palette="A_pal")
exp2_sct_counts_Vis_glut_scaled_B = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut, gene_list=L23_B_merfish_newNames, color_palette="B_pal")
exp2_sct_counts_Vis_glut_scaled_C = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut, gene_list=L23_C_merfish_newNames, color_palette="C_pal")

exp6_sct_counts_Vis_glut_scaled_A = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut, gene_list=L23_A_merfish_newNames, color_palette="A_pal")
exp6_sct_counts_Vis_glut_scaled_B = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut, gene_list=L23_B_merfish_newNames, color_palette="B_pal")
exp6_sct_counts_Vis_glut_scaled_C = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut, gene_list=L23_C_merfish_newNames, color_palette="C_pal")

exp7_sct_counts_Vis_glut_scaled_A = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut, gene_list=L23_A_merfish_newNames, color_palette="A_pal")
exp7_sct_counts_Vis_glut_scaled_B = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut, gene_list=L23_B_merfish_newNames, color_palette="B_pal")
exp7_sct_counts_Vis_glut_scaled_C = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut, gene_list=L23_C_merfish_newNames, color_palette="C_pal")


limx<-c(1800,4400)
limy<-c(7200,9900)

#exp1
exp1_L23_L <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_L_VIS_L23.csv")
exp1_L23_R <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_R_VIS_L23.csv")
names(exp1_L23_L)="index"
names(exp1_L23_R)="index"
exp1_L23_cells <- c(exp1_L23_L$index, exp1_L23_R$index)

#exp2
exp2_L23_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS_L23.csv")
exp2_L23_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS_L23.csv")
names(exp2_L23_L)="index"
names(exp2_L23_R)="index"
exp2_L23_cells <- c(exp2_L23_L$index, exp2_L23_R$index)


exp2_VIS_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS.csv")
exp2_VIS_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS.csv")
names(exp2_VIS_L)="index"
names(exp2_VIS_R)="index"
exp2_VIS_cells <- c(exp2_VIS_L$index, exp2_VIS_R$index)

#exp6
exp6_L23_L <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_L_VIS_L23.csv")
exp6_L23_R <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_R_VIS_L23.csv")
names(exp6_L23_L)="index"
names(exp6_L23_R)="index"
exp6_L23_cells <- c(exp6_L23_L$index, exp6_L23_R$index)

#exp7
exp7_L23_L <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_L_VIS_L23.csv")
exp7_L23_R <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_R_VIS_L23.csv")
names(exp7_L23_L)="index"
names(exp7_L23_R)="index"
exp7_L23_cells <- c(exp7_L23_L$index, exp7_L23_R$index)

max(exp2_sct_counts_Vis_glut[index %in% exp2_L23_L$index,center_x])
min(exp2_sct_counts_Vis_glut[index %in% exp2_L23_L$index,center_x])

limx_L23 <- c(2000,4200)
limy_L23 <- c(7800, 9200)
#A
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_WT_Cheng_A_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color])
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_KO_Cheng_A_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color])
dev.off()
#B
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_WT_Cheng_B_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_KO_Cheng_B_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()
#C
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_WT_Cheng_C_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_avgZscore_KO_Cheng_C_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()
#subtypes
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_subtypeColors_WT_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),cluster_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_subtypeColors_KO_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),cluster_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_subtypeColors_WT_scatterplot.png")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),cluster_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2_L23_subtypeColors_KO_scatterplot.png")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),cluster_color],xlim=limx_L23, ylim=limy_L23,main="")
dev.off()


plot(exp2_sct_counts_Vis_glut[,center_x],exp2_sct_counts_Vis_glut[,center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT")
points(exp2_sct_counts_Vis_glut[t.type.Under1Over1=="WT",center_x],exp2_sct_counts_Vis_glut[t.type.Under1Over1=="WT",center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[t.type.Under1Over1=="WT",cluster_color],xlim=limx_L23, ylim=limy_L23,main="")
points(exp2_sct_counts_Vis_glut[index %in% exp2_L23_cells,center_x],exp2_sct_counts_Vis_glut[index %in% exp2_L23_cells,center_y],pch = 20,cex=.5, col = "black",xlim=limx_L23, ylim=limy_L23,main="")


plot(exp2_sct_counts_Vis_glut[,center_x],exp2_sct_counts_Vis_glut[,center_y],pch = 20,cex=.3,col = "white",xlim=limx, ylim=limy,main="KO")
points(exp2_sct_counts_Vis_glut[t.type.Under1Over1=="KO",center_x],exp2_sct_counts_Vis_glut[t.type.Under1Over1=="KO",center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut[t.type.Under1Over1=="KO",cluster_color],xlim=limx, ylim=limy,main="")


C_pal2 <- colorRampPalette(c('white','green3'))
exp2_sct_counts_Vis_glut_scaled_C2 <- exp2_sct_counts_Vis_glut_scaled_C
exp2_sct_counts_Vis_glut_scaled_C2$plot_color <- C_pal2(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_scaled_C2$avg_zScore, breaks=20))]

plot(exp2_sct_counts_Vis_glut_scaled_C2[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C2[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_C2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1, height=5)
  plot(c(0,5), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,5,y+1/scale, col=lut[i], border=NA)
  }
}

#png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/colorBar_exp2_L23_avgZscore_Cheng_A_genes_scatterplot.png")
color.bar(A_pal(20), min=signif(min(exp2_sct_counts_Vis_glut_scaled_A$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_scaled_A$avg_zScore),2), nticks=2, title="A")

#color.bar(A_pal(20), min=signif(min(exp2_sct_counts_Vis_glut_scaled_A$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_scaled_A$avg_zScore),2), nticks=3, title="A" , output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_Cheng_A_zScore_colorKey.png")


color.bar(B_pal(20), min=signif(min(exp2_sct_counts_Vis_glut_scaled_B$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_scaled_B$avg_zScore),2), nticks=2, title="B")

color.bar(C_pal(20), min=signif(min(exp2_sct_counts_Vis_glut_scaled_C$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_scaled_C$avg_zScore),2), nticks=2, title="C")

#install.packages("fields")
library(fields)

exp2_sct_counts_Vis_glut_scaled_A[order(avg_zScore), avg_zScore]
#trying to plot color bar on scatterplot, doesn't look very good
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
colorbar.plot(x=3500, y=7800, strip=exp2_sct_counts_Vis_glut_scaled_A[order(avg_zScore), avg_zScore], col=exp2_sct_counts_Vis_glut_scaled_A[order(avg_zScore), plot_color])

plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
colorbar.plot(x=3500, y=7800, strip=exp2_sct_counts_Vis_glut_scaled_B[order(avg_zScore), avg_zScore], col=exp2_sct_counts_Vis_glut_scaled_B[order(avg_zScore), plot_color])

plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.5, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")
colorbar.plot(x=3500, y=7800, strip=exp2_sct_counts_Vis_glut_scaled_C[order(avg_zScore), avg_zScore], col=exp2_sct_counts_Vis_glut_scaled_C[order(avg_zScore), plot_color])





##exp1, left hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_L_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()



##exp1, right hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_A[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_B[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp1_L23_R_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp1, x", ylab="Exp1, y")
points(exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp1_sct_counts_Vis_glut_scaled_C[(index %in% exp1_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()


##exp2, left hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_A_genes_scatterplot.png")
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_A_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_A_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_B_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_B_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_Cheng_C_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_Cheng_C_genes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#subtypes
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_WT_subtypes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),cluster_color],main="")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscore_KO_subtypes_scatterplot.eps")
plot(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO", xlab="Exp2, x", ylab="Exp2, y", asp=1, xlim=limx_L23, ylim=limy_L23)
points(exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),cluster_color],main="")
dev.off()


##exp2, right hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_A[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_B[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_R_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()


##exp6, left hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_L_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()



##exp6, right hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_A[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_B[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp6_L23_R_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp6, x", ylab="Exp6, y")
points(exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp6_sct_counts_Vis_glut_scaled_C[(index %in% exp6_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()


##exp7, left hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_L_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()



##exp7, right hemisphere
#A
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_WT_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, A", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_KO_Cheng_A_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, A", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_A[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#B
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_WT_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, B", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_KO_Cheng_B_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, B", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_B[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()

#C
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_WT_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
#plot(exp7_sct_counts_Vis_glut_scaled_C[,center_x],exp7_sct_counts_Vis_glut_scaled_C[,center_y],pch = 20,cex=.3,col = "gray", main="WT, C", xlab="Exp7, x", ylab="Exp7, y")
plot(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="WT, C", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="WT"),plot_color],main="")
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp7_L23_R_avgZscore_KO_Cheng_C_genes_scatterplot.png", width=2000, height=2000, res=300)
plot(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index),center_y],pch = 20,cex=.3,col = "white", main="KO, C", xlab="Exp7, x", ylab="Exp7, y")
points(exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp7_L23_R$index) & (t.type.Under1Over1=="KO"),plot_color],main="")
dev.off()




plot(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_cells),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_cells),center_y],pch = 20,cex=.3,col = "gray", main="KO, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_R$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = "black",main="")
points(exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=.8, col = exp7_sct_counts_Vis_glut_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],main="")

rbPal <- colorRampPalette(c('yellow','darkgreen'))

exp2_sct_counts_Vis_glut_scaled_C3 <- exp2_sct_counts_Vis_glut_scaled_C
exp2_sct_counts_Vis_glut_scaled_C3$plot_color <- rbPal(50)[as.numeric(cut(exp2_sct_counts_Vis_glut_scaled_C3$avg_zScore, breaks=50))]

plot(exp2_sct_counts_Vis_glut_scaled_C3[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_scaled_C3[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_scaled_C3[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_scaled_C3[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_scaled_C3[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")

###Restricting z-score calculation to just L2/3
exp1_sct_counts_Vis_glut_L23 <- exp1_sct_counts_Vis_glut[index %in% exp1_L23_cells]
exp1_sct_counts_Vis_glut_L23_genesOnly <- data.frame(exp1_sct_counts_Vis_glut_L23[, 1:551], row.names="index")
exp1_sct_counts_Vis_glut_L23_scaled <- scale(exp1_sct_counts_Vis_glut_L23_genesOnly)

exp2_sct_counts_Vis_glut_L23 <- exp2_sct_counts_Vis_glut[index %in% exp2_L23_cells]
exp2_sct_counts_Vis_glut_L23_genesOnly <- data.frame(exp2_sct_counts_Vis_glut_L23[, 1:551], row.names="index")
exp2_sct_counts_Vis_glut_L23_scaled <- scale(exp2_sct_counts_Vis_glut_L23_genesOnly)

exp6_sct_counts_Vis_glut_L23 <- exp6_sct_counts_Vis_glut[index %in% exp6_L23_cells]
exp6_sct_counts_Vis_glut_L23_genesOnly <- data.frame(exp6_sct_counts_Vis_glut_L23[, 1:551], row.names="index")
exp6_sct_counts_Vis_glut_L23_scaled <- scale(exp6_sct_counts_Vis_glut_L23_genesOnly)

exp7_sct_counts_Vis_glut_L23 <- exp7_sct_counts_Vis_glut[index %in% exp7_L23_cells]
exp7_sct_counts_Vis_glut_L23_genesOnly <- data.frame(exp7_sct_counts_Vis_glut_L23[, 1:551], row.names="index")
exp7_sct_counts_Vis_glut_L23_scaled <- scale(exp7_sct_counts_Vis_glut_L23_genesOnly)


A_pal3 <- colorRampPalette(c('lightblue','purple3'))
#B_pal3 <- colorRampPalette(c('pink', 'red'))
B_pal3 <- colorRampPalette(c('yellow', 'orangered1'))
C_pal3 <- colorRampPalette(c('yellow','darkgreen'))

exp1_sct_counts_Vis_glut_L23_scaled_A = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut_L23, gene_list=L23_A_merfish_newNames, color_palette="A_pal3")
exp1_sct_counts_Vis_glut_L23_scaled_B = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut_L23, gene_list=L23_B_merfish_newNames, color_palette="B_pal3")
exp1_sct_counts_Vis_glut_L23_scaled_C = exp_counts_scale_func(exp_scaled=exp1_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp1_sct_counts_Vis_glut_L23, gene_list=L23_C_merfish_newNames, color_palette="C_pal3")

exp2_sct_counts_Vis_glut_L23_scaled_A = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut_L23, gene_list=L23_A_merfish_newNames, color_palette="A_pal3")
exp2_sct_counts_Vis_glut_L23_scaled_B = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut_L23, gene_list=L23_B_merfish_newNames, color_palette="B_pal3")
exp2_sct_counts_Vis_glut_L23_scaled_C = exp_counts_scale_func(exp_scaled=exp2_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp2_sct_counts_Vis_glut_L23, gene_list=L23_C_merfish_newNames, color_palette="C_pal3")

exp6_sct_counts_Vis_glut_L23_scaled_A = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut_L23, gene_list=L23_A_merfish_newNames, color_palette="A_pal3")
exp6_sct_counts_Vis_glut_L23_scaled_B = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut_L23, gene_list=L23_B_merfish_newNames, color_palette="B_pal3")
exp6_sct_counts_Vis_glut_L23_scaled_C = exp_counts_scale_func(exp_scaled=exp6_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp6_sct_counts_Vis_glut_L23, gene_list=L23_C_merfish_newNames, color_palette="C_pal3")

exp7_sct_counts_Vis_glut_L23_scaled_A = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut_L23, gene_list=L23_A_merfish_newNames, color_palette="A_pal3")
exp7_sct_counts_Vis_glut_L23_scaled_B = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut_L23, gene_list=L23_B_merfish_newNames, color_palette="B_pal3")
exp7_sct_counts_Vis_glut_L23_scaled_C = exp_counts_scale_func(exp_scaled=exp7_sct_counts_Vis_glut_L23_scaled, exp_counts_summary=exp7_sct_counts_Vis_glut_L23, gene_list=L23_C_merfish_newNames, color_palette="C_pal3")


plot(exp2_sct_counts_Vis_glut_L23_scaled_A[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_A[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, A", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_L23_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_A[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")


plot(exp2_sct_counts_Vis_glut_L23_scaled_B[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_B[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, B", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_L23_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_B[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")


plot(exp2_sct_counts_Vis_glut_L23_scaled_C[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_C[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, C", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_L23_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_C[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")

color.bar(A_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_A$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_A$avg_zScore),2), nticks=2, title="A")
color.bar(B_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_B$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_B$avg_zScore),2), nticks=2, title="B")
color.bar(C_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_C$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_C$avg_zScore),2), nticks=2, title="C")

L23_plot_func <- function(exp_data_table, cell.ids, t.type, color_column, output_file){
  setEPS()
  postscript(output_file)
  plot(exp_data_table[(index %in% cell.ids),center_x],exp_data_table[(index %in% cell.ids),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23, main="", xlab="x", ylab="y", asp=1)
  points(exp_data_table[(index %in% cell.ids) & (t.type.Under1Over1==t.type),center_x],exp_data_table[(index %in% cell.ids) & (t.type.Under1Over1==t.type),center_y],pch = 20,cex=1.5, col = exp_data_table[(index %in% cell.ids) & (t.type.Under1Over1==t.type),get(color_column)])
  dev.off()
}
#all
L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23, cell.ids=exp2_L23_L$index, t.type="WT", color_column="cluster_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_WT_subtypes_scatterplot.eps")
L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23, cell.ids=exp2_L23_L$index, t.type="KO", color_column="cluster_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_KO_subtypes_scatterplot.eps")

L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_A, cell.ids=exp2_L23_L$index, t.type="WT", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_WT_Cheng_A_genes_scatterplot.eps")
L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_A, cell.ids=exp2_L23_L$index, t.type="KO", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_KO_Cheng_A_genes_scatterplot.eps")
color.bar(A_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_A$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_A$avg_zScore),2), nticks=2, title="A")

L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_B, cell.ids=exp2_L23_L$index, t.type="WT", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_WT_Cheng_B_genes_scatterplot.eps")
L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_B, cell.ids=exp2_L23_L$index, t.type="KO", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_KO_Cheng_B_genes_scatterplot.eps")
color.bar(B_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_B$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_B$avg_zScore),2), nticks=2, title="B")


L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_C, cell.ids=exp2_L23_L$index, t.type="WT", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_WT_Cheng_C_genes_scatterplot.eps")
L23_plot_func(exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_C, cell.ids=exp2_L23_L$index, t.type="KO", color_column="plot_color", output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_avgZscoreAcrossL23_KO_Cheng_C_genes_scatterplot.eps")
color.bar(C_pal3(20), min=signif(min(exp2_sct_counts_Vis_glut_L23_scaled_C$avg_zScore),2), max=signif(max(exp2_sct_counts_Vis_glut_L23_scaled_C$avg_zScore),2), nticks=2, title="C")


#specific genes
View(Vis_Liu2023_TOI[gene %in% L23_A_merfish_newNames])
View(Vis_Liu2023_TOI[gene %in% L23_B_merfish_newNames])
View(Vis_Liu2023_TOI[gene %in% L23_C_merfish_newNames])

join_cols_exp2_sct_counts_Vis_glut_L23 <- c("index", names(exp2_sct_counts_Vis_glut_L23)[552:573])

exp2_sct_counts_Vis_glut_L23_scaled_dt <- data.frame(exp2_sct_counts_Vis_glut_L23_scaled) %>% data.table(keep.rownames="index")
exp2_sct_counts_Vis_glut_L23_scaled_dt <- inner_join(x=exp2_sct_counts_Vis_glut_L23_scaled_dt, y=exp2_sct_counts_Vis_glut_L23[, ..join_cols_exp2_sct_counts_Vis_glut_L23], by="index")

#A example
A_example = "Spon1"
exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example <- exp2_sct_counts_Vis_glut_L23_scaled_dt
exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example$plot_color <- A_pal3(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[, get(A_example)], breaks=20))]

plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, Spon1", xlab="Exp2, x", ylab="Exp2, y", asp=1)
points(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")

plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO, Spon1", xlab="Exp2, x", ylab="Exp2, y")
points(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")

example_gene_plot_png_func <- function(example_gene, exp_data_table, cell.id, color_pal, color_column, xlim=limx_L23, ylim=limy_L23, output_file){
  exp_data_table_example <- exp_data_table
  exp_data_table_example$plot_color <- get(color_pal)(20)[as.numeric(cut(exp_data_table_example[, get(example_gene)], breaks=20))]
  
  png(paste0(output_file, "_WT.png"), width=2000, height=2000, res=300)
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="WT", xlab="x", ylab="y", asp=1)
  points(exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),center_x],exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),get(color_column)],xlim=xlim, ylim=ylim)
  dev.off()
  
  png(paste0(output_file, "_KO.png"), width=2000, height=2000, res=300)
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="KO", xlab="x", ylab="y", asp=1)
  points(exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),center_x],exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=1.5, col = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),get(color_column)],xlim=xlim, ylim=ylim)
  dev.off()
  
  color.bar(get(color_pal)(20), min=signif(min(exp_data_table_example[, get(example_gene)]),2), max=signif(max(exp_data_table_example[, get(example_gene)]),2), nticks=2, title=example_gene)
}

example_gene_plot_eps_func <- function(example_gene, exp_data_table, cell.id, color_pal, color_column, xlim=limx_L23, ylim=limy_L23, output_file, dot_size=1.5){
  exp_data_table_example <- exp_data_table
  exp_data_table_example$plot_color <- get(color_pal)(20)[as.numeric(cut(exp_data_table_example[, get(example_gene)], breaks=20))]
  
  setEPS()
  postscript(paste0(output_file, "_WT.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="WT", xlab="x", ylab="y", asp=1)
  points(exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),center_x],exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=dot_size, col = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT"),get(color_column)],xlim=xlim, ylim=ylim)
  dev.off()
  
  setEPS()
  postscript(paste0(output_file, "_KO.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="KO", xlab="x", ylab="y", asp=1)
  points(exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),center_x],exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),center_y],pch = 20,cex=dot_size, col = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO"),get(color_column)],xlim=xlim, ylim=ylim)
  dev.off()
}

example_gene_plot_png_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1")
example_gene_plot_png_func(example_gene="Kcnn3", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Kcnn3")
example_gene_plot_png_func(example_gene="Cdh13", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Cdh13")
example_gene_plot_png_func(example_gene="Rfx3", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Rfx3")
example_gene_plot_png_func(example_gene="Egfem1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Egfem1")
example_gene_plot_png_func(example_gene="Syt10", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Syt10")
example_gene_plot_png_func(example_gene="Col23a1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Col23a1")


example_gene_plot_png_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6")
example_gene_plot_png_func(example_gene="Nell1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Nell1")
example_gene_plot_png_func(example_gene="Rgs8", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Rgs8")


example_gene_plot_png_func(example_gene="Bmper", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Bmper")
example_gene_plot_png_func(example_gene="Zmat4", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Zmat4")
example_gene_plot_png_func(example_gene="Tshz1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1")

#eps
example_gene_plot_eps_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1")
example_gene_plot_eps_func(example_gene="Kcnn3", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Kcnn3")
example_gene_plot_eps_func(example_gene="Cdh13", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Cdh13")

example_gene_plot_eps_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1_2", dot_size=2.5)


example_gene_plot_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6")
example_gene_plot_eps_func(example_gene="Nell1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Nell1")
example_gene_plot_eps_func(example_gene="Rgs8", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Rgs8", dot_size=2.5)

example_gene_plot_eps_func(example_gene="Bmper", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Bmper")
example_gene_plot_eps_func(example_gene="Zmat4", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Zmat4")
example_gene_plot_eps_func(example_gene="Tshz1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1", dot_size=2.5)


example_gene_plot_eps_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1", dot_size=2)
example_gene_plot_eps_func(example_gene="Rgs8", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Rgs8", dot_size=2)
example_gene_plot_eps_func(example_gene="Tshz1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1", dot_size=2)

example_gene_plot_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6", dot_size=2)
example_gene_plot_eps_func(example_gene="Rfx3", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Rfx3", dot_size=2)
##
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.4, xlab="x", ylab="y", col="gray")
points(exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=0.8, col = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),cluster_color])

plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT, Spon1", xlab="Exp2, x", ylab="Exp2, y", asp=1)
points(exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=1.5, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_A_example[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),plot_color],xlim=limx_L23, ylim=limy_L23,main="")


exp2_sct_counts_Vis_glut_L23_scaled_dt2 <- copy(exp2_sct_counts_Vis_glut_L23_scaled_dt)


library(scales)
hue_pal()(3)
type_colors_new <- cbind(predicted.id=unique(exp2_sct_counts_Vis_glut_L23_scaled_dt2[, predicted.id]),
      new_color=hue_pal()(length(unique(exp2_sct_counts_Vis_glut_L23_scaled_dt2[, predicted.id])))) %>% data.table



exp2_sct_counts_Vis_glut_L23_scaled_dt2 <- left_join(x=exp2_sct_counts_Vis_glut_L23_scaled_dt2, y=type_colors_new)

plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.4, xlab="x", ylab="y", col="gray")
points(exp2_sct_counts_Vis_glut_L23_scaled_dt2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),center_y],pch = 20,cex=0.5, col = exp2_sct_counts_Vis_glut_L23_scaled_dt2[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT"),new_color])
legend(x=6000,
       y=4000,
       legend=type_colors_new$predicted.id,
       col=type_colors_new$new_color, lty=1:2, cex=0.8)

exp2_sct_counts_Vis_glut_L23_scaled_dt3 <- copy(exp2_sct_counts_Vis_glut_L23_scaled_dt)

type_palette <- rainbow(length(unique(exp2_sct_counts_Vis_glut_L23_scaled_dt2[, predicted.id])))                                 # Apply rainbow function
type_colors_new2 <- cbind(predicted.id=unique(exp2_sct_counts_Vis_glut_L23_scaled_dt3[, predicted.id]),
                          new_color=type_palette) %>% data.table

exp2_sct_counts_Vis_glut_L23_scaled_dt3 <- left_join(x=exp2_sct_counts_Vis_glut_L23_scaled_dt3, y=type_colors_new2)

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_subtypes_newColors_scatterplot.eps")
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.25, xlab="x", ylab="y", col="gray", asp=1)
points(exp2_sct_counts_Vis_glut_L23_scaled_dt3[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt3[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=0.35, col = exp2_sct_counts_Vis_glut_L23_scaled_dt3[(index %in% exp2_L23_L$index),new_color])
legend(x=6000,
       y=4000,
       legend=type_colors_new2$predicted.id,
       col=type_colors_new2$new_color, lty=1:2, cex=0.4)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_VIS_L_colored_scatterplot.eps")
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.25, xlab="x", ylab="y", col="lightgray", asp=1)
points(exp2_sct_counts[(index %in% exp2_VIS_L$index),center_x],exp2_sct_counts[(index %in% exp2_VIS_L$index),center_y],pch = 20,cex=0.35, col = "darkgray")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_colored_scatterplot.eps")
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.25, xlab="x", ylab="y", col="lightgray", asp=1)
points(exp2_sct_counts[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=0.35, col = "darkgray")
dev.off()

plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.25, xlab="x", ylab="y", col="lightgray", asp=1)
points(exp2_sct_counts[(index %in% exp2_VIS_L$index),center_x],exp2_sct_counts[(index %in% exp2_VIS_L$index),center_y],pch = 20,cex=0.35, col = "gray1")
points(exp2_sct_counts[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=0.35, col = "black")


#plotting only relevant cell types
exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt <- exp2_sct_counts_Vis_glut_L23_scaled_dt[predicted.id %in% types_of_interest]

hue_pal()(4)
new_colors_TOI <- cbind(predicted.id=types_of_interest,
                        new_color = c("magenta", "cyan", "orange", "darkred")) %>% data.table

exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt <- left_join(x=exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt, y=new_colors_TOI)

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/exp2_L23_L_WTandKO_L23_scatterplot_newColors_166_171_168_179.eps")
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.3, xlab="x", ylab="y", col="white", xlim=limx_L23, ylim=limy_L23, asp=1, main="WT and KO")
points(exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=2.5, col = exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index),new_color])
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/cex2_exp2_L23_L_WTandKO_L23_scatterplot_newColors_166_171_168_179.eps")
plot(exp2_sct_counts[, center_x], exp2_sct_counts[, center_y], pch=20, cex=0.3, xlab="x", ylab="y", col="white", xlim=limx_L23, ylim=limy_L23, asp=1, main="WT and KO")
points(exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("168_L2/3 IT CTX", "179_L4 IT CTX")),center_x],exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("168_L2/3 IT CTX", "179_L4 IT CTX")),center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("168_L2/3 IT CTX", "179_L4 IT CTX")),new_color])
points(exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("166_L2/3 IT CTX", "171_L2/3 IT CTX")),center_x],exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("166_L2/3 IT CTX", "171_L2/3 IT CTX")),center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_TOI_dt[(index %in% exp2_L23_L$index) & (predicted.id %in% c("166_L2/3 IT CTX", "171_L2/3 IT CTX")),new_color])
dev.off()

#cell type identities
exp2_sct_counts_Vis_glut_L23_scaled_dt

exp2_sct_counts_Vis_glut_L23_scaled_dt_L_summary = exp2_sct_counts_Vis_glut_L23_scaled_dt[index %in% exp2_L23_L$index] %>% 
                                                    group_by(predicted.id) %>%
                                                    summarise(n.val = n(),
                                                              n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
                                                              n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
                                                              ) %>% data.table

exp1_sct_counts_Vis_glut_L23_summary = exp1_sct_counts_Vis_glut_L23[index %in% exp1_L23_cells] %>% 
  group_by(predicted.id) %>%
  summarise(n.total = n(),
            n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
            n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
  ) %>% data.table

exp2_sct_counts_Vis_glut_L23_summary = exp2_sct_counts_Vis_glut_L23[index %in% exp2_L23_cells] %>% 
  group_by(predicted.id) %>%
  summarise(n.total = n(),
            n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
            n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
  ) %>% data.table

exp6_sct_counts_Vis_glut_L23_summary = exp6_sct_counts_Vis_glut_L23[index %in% exp6_L23_cells] %>% 
  group_by(predicted.id) %>%
  summarise(n.total = n(),
            n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
            n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
  ) %>% data.table

exp7_sct_counts_Vis_glut_L23_summary = exp7_sct_counts_Vis_glut_L23[index %in% exp7_L23_cells] %>% 
  group_by(predicted.id) %>%
  summarise(n.total = n(),
            n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
            n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
  ) %>% data.table

exp1_sct_counts_Vis_glut_L23_summary$WT_KO_ratio <- exp1_sct_counts_Vis_glut_L23_summary[, n.WT/n.KO]
exp2_sct_counts_Vis_glut_L23_summary$WT_KO_ratio <- exp2_sct_counts_Vis_glut_L23_summary[, n.WT/n.KO]
exp6_sct_counts_Vis_glut_L23_summary$WT_KO_ratio <- exp6_sct_counts_Vis_glut_L23_summary[, n.WT/n.KO]
exp7_sct_counts_Vis_glut_L23_summary$WT_KO_ratio <- exp7_sct_counts_Vis_glut_L23_summary[, n.WT/n.KO]

all_exp_sct_counts_Vis_glut_L23 <- rbind(exp1_sct_counts_Vis_glut_L23, exp2_sct_counts_Vis_glut_L23, exp6_sct_counts_Vis_glut_L23, exp7_sct_counts_Vis_glut_L23)
all_exp_sct_counts_Vis_glut_L23_summary = all_exp_sct_counts_Vis_glut_L23 %>% 
  group_by(predicted.id) %>%
  summarise(n.total = n(),
            n.WT = sum((t.type.Under1Over1=="WT"), na.rm=TRUE),
            n.KO = sum((t.type.Under1Over1=="KO"), na.rm=TRUE)
  ) %>% data.table
all_exp_sct_counts_Vis_glut_L23_summary$WT_KO_ratio <- all_exp_sct_counts_Vis_glut_L23_summary[, n.WT/n.KO]

write.csv(all_exp_sct_counts_Vis_glut_L23_summary, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/type_counts/L23.Vis.glut.cell.type.counts.exp1.exp2.exp6.exp7.100vol.300counts.CTX.HC.pred0.2.csv", quote=F, row.names=F)
write.csv(exp1_sct_counts_Vis_glut_L23_summary, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/type_counts/L23.Vis.glut.cell.type.counts.exp1.100vol.300counts.CTX.HC.pred0.2.csv", quote=F, row.names=F)
write.csv(exp2_sct_counts_Vis_glut_L23_summary, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/type_counts/L23.Vis.glut.cell.type.counts.exp2.100vol.300counts.CTX.HC.pred0.2.csv", quote=F, row.names=F)
write.csv(exp6_sct_counts_Vis_glut_L23_summary, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/type_counts/L23.Vis.glut.cell.type.counts.exp6.100vol.300counts.CTX.HC.pred0.2.csv", quote=F, row.names=F)
write.csv(exp7_sct_counts_Vis_glut_L23_summary, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/type_counts/L23.Vis.glut.cell.type.counts.exp7.100vol.300counts.CTX.HC.pred0.2.csv", quote=F, row.names=F)
##


ggplot(Vis_Liu2023_TOI_A[gene=="Spon1"], aes(x = subtype, y = as.numeric(logFC)))+
  ggtitle("Spon1")+
  geom_col(fill="purple3")+
  coord_cartesian(ylim=c(-1.7, 1.0))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y=element_text(size=10))

ggplot(Vis_Liu2023_TOI_B[gene=="Rgs8"], aes(x = subtype, y = as.numeric(logFC)))+
  ggtitle("Rgs8")+
  geom_col(fill="red")+
  coord_cartesian(ylim=c(-1.7, 1.0))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y=element_text(size=10))


ggplot(Vis_Liu2023_TOI_C[gene=="Tshz1"], aes(x = subtype, y = as.numeric(logFC)))+
  ggtitle("Tshz1")+
  geom_col(fill="darkgreen")+
  coord_cartesian(ylim=c(-1.7, 1.0))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title.y=element_text(size=10))

library(reshape2)

Vis_Liu2023_TOI[gene %in% c("Spon1", "Rgs8", "Tshz1")]

example_genes_dcast <- dcast(Vis_Liu2023_TOI[gene %in% c("Spon1", "Rgs8", "Tshz1")], subtype ~ gene, value.var = "logFC") %>% data.table
example_genes_dcast = example_genes_dcast[, .(subtype, Spon1, Rgs8, Tshz1)]
example_genes_matrix = t(as.matrix(example_genes_dcast, rownames="subtype"))

example_genes_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
#example_genes_col_breaks = c(seq(-1., 1.7,length=100),  # for blue
#               seq(-0.99,0.99,length=100), #for white
#               seq(1, 3, length=100)) #for red
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_ABC_Spon1_Rgs8_Tshz1_logFC_merfish_heatmap.eps")
heatmap.2(example_genes_matrix,
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, col=example_genes_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()


#example genes with Trpc6 instead of Rgs8
Spon1_Trpc6_Tshz1_dcast <- dcast(Vis_Liu2023_TOI[gene %in% c("Spon1", "Trpc6", "Tshz1")], subtype ~ gene, value.var = "logFC") %>% data.table
Spon1_Trpc6_Tshz1_dcast = Spon1_Trpc6_Tshz1_dcast[, .(subtype, Spon1, Trpc6, Tshz1)]
Spon1_Trpc6_Tshz1_matrix = t(as.matrix(Spon1_Trpc6_Tshz1_dcast, rownames="subtype"))

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_ABC_Spon1_Trpc6_Tshz1_logFC_merfish_heatmap.eps")
heatmap.2(Spon1_Trpc6_Tshz1_matrix,
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, col=example_genes_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()
###
summary_table_cols <- c("index", "avg_zScore", names(exp1_sct_counts_Vis_glut_L23_scaled_A)[553:570])
summ_exp_sct_counts_Vis_glut_L23_scaled <- rbind(
  cbind(exp1_sct_counts_Vis_glut_L23_scaled_A[, ..summary_table_cols], gene_list="A_avg_zScore"),
  cbind(exp1_sct_counts_Vis_glut_L23_scaled_B[, ..summary_table_cols], gene_list="B_avg_zScore"),
  cbind(exp1_sct_counts_Vis_glut_L23_scaled_C[, ..summary_table_cols], gene_list="C_avg_zScore"),
  cbind(exp2_sct_counts_Vis_glut_L23_scaled_A[, ..summary_table_cols], gene_list="A_avg_zScore"),
  cbind(exp2_sct_counts_Vis_glut_L23_scaled_B[, ..summary_table_cols], gene_list="B_avg_zScore"),
  cbind(exp2_sct_counts_Vis_glut_L23_scaled_C[, ..summary_table_cols], gene_list="C_avg_zScore"),
  cbind(exp6_sct_counts_Vis_glut_L23_scaled_A[, ..summary_table_cols], gene_list="A_avg_zScore"),
  cbind(exp6_sct_counts_Vis_glut_L23_scaled_B[, ..summary_table_cols], gene_list="B_avg_zScore"),
  cbind(exp6_sct_counts_Vis_glut_L23_scaled_C[, ..summary_table_cols], gene_list="C_avg_zScore"),
  cbind(exp7_sct_counts_Vis_glut_L23_scaled_A[, ..summary_table_cols], gene_list="A_avg_zScore"),
  cbind(exp7_sct_counts_Vis_glut_L23_scaled_B[, ..summary_table_cols], gene_list="B_avg_zScore"),
  cbind(exp7_sct_counts_Vis_glut_L23_scaled_C[, ..summary_table_cols], gene_list="C_avg_zScore")) %>% data.table

summ_exp_sct_counts_Vis_glut_L23_scaled <- left_join(x=summ_exp_sct_counts_Vis_glut_L23_scaled, y=All_Vis_layerdepth[, .(cell.id, NDR)], by=c("index"="cell.id"))
summ_exp_sct_counts_Vis_glut_L23_scaled <- summ_exp_sct_counts_Vis_glut_L23_scaled[, .(index, gene_list, avg_zScore, rep,
                                            t.type.Under1Over1, t.type.Under1Over2, prediction.score.max,
                                            predicted.id, subclass_label, neighborhood_label, class_label, NDR,
                                            total_sct_counts,total_raw_counts, fov, volume, 
                                            center_x, center_y,
                                            min_x, max_x, min_y, max_y)]

summ_exp_sct_counts_Vis_glut_L23_scaled_dcast <- dcast(summ_exp_sct_counts_Vis_glut_L23_scaled, index + rep +
      t.type.Under1Over1 + t.type.Under1Over2 + prediction.score.max +
      predicted.id + subclass_label + neighborhood_label + class_label + NDR +
      total_sct_counts + total_raw_counts + fov + volume + 
      center_x + center_y +
      min_x + max_x + min_y + max_y 
        ~ gene_list, 
      value.var = "avg_zScore") %>% data.table

write.csv(summ_exp_sct_counts_Vis_glut_L23_scaled_dcast, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/mecp2Het_VIS_glut_L23_Cheng2022_genes_sct_avgZscores_layerDepths.csv", quote=F, row.names=F)

#
allExp_L23_cells = data.table(c(exp1_L23_cells, exp2_L23_cells, exp6_L23_cells, exp7_L23_cells))
names(allExp_L23_cells)="cell.id"
write.csv(allExp_L23_cells, file = "HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/exp1_exp2_exp6_exp7_L23_cells.csv", quote=F, row.names=F)

####L23 and VIS only for calculation of gene expression fold changes
L23_Vis_logFC <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/L23.Vis.subtypeAgg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

L23_Vis_logFC_TOI <- L23_Vis_logFC[subtype %in% types_of_interest]

L23_Vis_logFC_TOI = L23_Vis_logFC_TOI %>% mutate(subtype = factor(subtype, levels=types_of_interest))



L23_Vis_logFC_TOI_A <- rbind(
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-A"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "A"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-A"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "A"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-A"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "A"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_A_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-A"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_A_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "A") 
) %>% data.table

L23_Vis_logFC_TOI_B <- rbind(
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "B"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "B"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-B"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "B"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_B_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-B"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_B_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "B") 
) %>% data.table

L23_Vis_logFC_TOI_C <- rbind(
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="166_L2/3 IT CTX")], gene_list = "C"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="171_L2/3 IT CTX")], gene_list = "C"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "non-C"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="168_L2/3 IT CTX")], gene_list = "C"), 
  cbind(L23_Vis_logFC_TOI[(gene %in% nonL23_C_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "non-C"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_C_merfish_newNames) & (subtype=="179_L4 IT CTX")], gene_list = "C") 
) %>% data.table


L23_Vis_logFC_TOI_A = L23_Vis_logFC_TOI_A %>% mutate(subtype = factor(subtype, levels=types_of_interest))
L23_Vis_logFC_TOI_A = L23_Vis_logFC_TOI_A %>% mutate(gene_list = factor(gene_list, levels=c("non-A", "A")))
L23_Vis_logFC_TOI_B = L23_Vis_logFC_TOI_B %>% mutate(subtype = factor(subtype, levels=types_of_interest))
L23_Vis_logFC_TOI_B = L23_Vis_logFC_TOI_B %>% mutate(gene_list = factor(gene_list, levels=c("non-B", "B")))
L23_Vis_logFC_TOI_C = L23_Vis_logFC_TOI_C %>% mutate(subtype = factor(subtype, levels=types_of_interest))
L23_Vis_logFC_TOI_C = L23_Vis_logFC_TOI_C %>% mutate(gene_list = factor(gene_list, levels=c("non-C", "C")))

ggplot(L23_Vis_logFC_TOI_A, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_A genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","purple3")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonA.vs.A.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonA.vs.A.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


ggplot(L23_Vis_logFC_TOI_B, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_B genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","red")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonB.vs.B.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonB.vs.B.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(L23_Vis_logFC_TOI_C, aes(x = gene_list, y = as.numeric(logFC), fill=gene_list))+
  ggtitle("L23_C genes")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Log2 fold change (KO/WT) ") + xlab("")+
  scale_fill_manual(name="Genes:", values = c("gray","darkgreen")) +
  facet_grid(.~subtype,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonC.vs.C.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/L23.Cheng2022.nonC.vs.C.genes.L23.VIS.2.Liu2023.subtypeAgg.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


sigs_nCells10_L23_Vis_logFC_TOI_A <- L23_vs_nonL23_sigs_func(logfc_table=L23_Vis_logFC_TOI_A, geneList1="A", geneList2="non-A")
sigs_nCells10_L23_Vis_logFC_TOI_B <- L23_vs_nonL23_sigs_func(logfc_table=L23_Vis_logFC_TOI_B, geneList1="B", geneList2="non-B")
sigs_nCells10_L23_Vis_logFC_TOI_C <- L23_vs_nonL23_sigs_func(logfc_table=L23_Vis_logFC_TOI_C, geneList1="C", geneList2="non-C")


sigs_L23_Vis_logFC_TOI_dt <- rbind(cbind(sigs_nCells10_L23_Vis_logFC_TOI_A, pseudobulkDGE_nCells_min = 10),
                                  cbind(sigs_nCells10_L23_Vis_logFC_TOI_B, pseudobulkDGE_nCells_min = 10),
                                  cbind(sigs_nCells10_L23_Vis_logFC_TOI_C, pseudobulkDGE_nCells_min = 10)) %>% data.table

write.csv(sigs_L23_Vis_logFC_TOI_dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/nCells10_VIS_L23_vs_nonL23_Cheng_boxplot_pvals.csv", row.names=F, quote=F)

#checking the above
wilcox.test(L23_Vis_logFC_TOI_A[(gene_list=="A") & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)], L23_Vis_logFC_TOI_A[(gene_list=="non-A") & (subtype=="166_L2/3 IT CTX"), as.numeric(logFC)])


#
#example genes with Trpc6 instead of Rgs8, using L2/3 VIS pseudobulkDGE
L23_Vis_Spon1_Trpc6_Tshz1_dcast <- dcast(L23_Vis_logFC_TOI[gene %in% c("Spon1", "Trpc6", "Tshz1")], subtype ~ gene, value.var = "logFC") %>% data.table
L23_Vis_Spon1_Trpc6_Tshz1_dcast = L23_Vis_Spon1_Trpc6_Tshz1_dcast[, .(subtype, Spon1, Trpc6, Tshz1)]
L23_Vis_Spon1_Trpc6_Tshz1_matrix = t(as.matrix(L23_Vis_Spon1_Trpc6_Tshz1_dcast, rownames="subtype"))

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_Vis_ABC_Spon1_Trpc6_Tshz1_logFC_merfish_heatmap.eps")
heatmap.2(L23_Vis_Spon1_Trpc6_Tshz1_matrix,
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, col=example_genes_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

#fold changes for just A, B, and C genes in types of interest
L23_Vis_logFC_TOI_ABC <- rbind(
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_A_merfish_newNames)], gene_list = "Superficial"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_B_merfish_newNames)], gene_list = "Intermediate"),
  cbind(L23_Vis_logFC_TOI[(gene %in% L23_C_merfish_newNames)], gene_list = "Deep")
  ) %>% data.table

names(L23_Vis_logFC_TOI_ABC)[1] = "type"
write.csv(L23_Vis_logFC_TOI_ABC, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/L23_Vis_pseudobulkDGE_logFC_Cheng2022_genes.csv", quote=F, row.names=F)


####
L23_Vis_WT_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
L23_Vis_KO_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

L23_Vis_WT_and_KO_avgExp_matrix_TOI <- left_join(x=L23_Vis_WT_avgExp_matrix[, c("Gene","166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX")],
                                             y=L23_Vis_KO_avgExp_matrix[, c("Gene","166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX")],
                                             by=c("Gene"), suffix = c(".WT", ".KO"))
L23_Vis_WT_and_KO_avgExp_matrix_TOI_ABC <- data.frame(L23_Vis_WT_and_KO_avgExp_matrix_TOI[Gene %in% L23_genes_newNames], row.names="Gene")

L23_Vis_WT_and_KO_avgExp_matrix_TOI_ABC_zscore <- t(scale(t(as.matrix(L23_Vis_WT_and_KO_avgExp_matrix_TOI_ABC ))))


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_Vis_only_Cheng2022_ABC_AverageExpression_zscores_merfish_heatmap.eps")
heatmap.2(L23_Vis_WT_and_KO_avgExp_matrix_TOI_ABC_zscore[c(L23_A_merfish_newNames, L23_B_merfish_newNames, L23_C_merfish_newNames), type_order],
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.2, 
          cexCol = 0.6, labRow = NULL, labCol = NULL, col=my_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

L23_Vis_gene_order <- rownames(L23_Vis_WT_and_KO_avgExp_matrix_TOI_ABC_zscore[c(L23_A_merfish_newNames, L23_B_merfish_newNames, L23_C_merfish_newNames), type_order])
original_gene_order <- rownames(WT_and_KO_log_avgExp_matrix_TOI_ABC_zscore[c(make.names(L23_A_merfish), make.names(L23_B_merfish), make.names(L23_C_merfish)), type_order])

####re-do overlap matrix of Cheng 2022 genes and significantly changed genes in Mecp2 KO/+, this time using values from only L2/3 Vis
#subtype significantly MR gene overlap with L2/3 A/B/C genes
L23_Vis_MR_L23_166_A <-  gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
L23_Vis_MR_L23_166_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
L23_Vis_MR_L23_166_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

L23_Vis_MR_L23_171_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
L23_Vis_MR_L23_171_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
L23_Vis_MR_L23_171_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

L23_Vis_MR_L23_168_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
L23_Vis_MR_L23_168_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
L23_Vis_MR_L23_168_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

L23_Vis_MR_L4_179_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
L23_Vis_MR_L4_179_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
L23_Vis_MR_L4_179_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

#subtype significantly MA gene overlap with L2/3 A/B/C genes
L23_Vis_MA_L23_166_A <-  gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
L23_Vis_MA_L23_166_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
L23_Vis_MA_L23_166_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="166_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

L23_Vis_MA_L23_171_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
L23_Vis_MA_L23_171_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
L23_Vis_MA_L23_171_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="171_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

L23_Vis_MA_L23_168_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
L23_Vis_MA_L23_168_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
L23_Vis_MA_L23_168_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="168_L2/3 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

L23_Vis_MA_L4_179_A <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
L23_Vis_MA_L4_179_B <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
L23_Vis_MA_L4_179_C <- gene_overlap_func(gene_panel_classes[,Gene], L23_Vis_logFC_TOI[(subtype=="179_L4 IT CTX") & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)
#matrix of p-values of overlaps between statistically Mecp2-regulated genes and L2/3 A/B/C genes
L23_Vis_sig_ABC_pvals_matrix <- matrix(c(unlist(L23_Vis_MR_L23_166_A[1]), unlist(L23_Vis_MR_L23_166_B[1]), unlist(L23_Vis_MR_L23_166_C[1]),
                                 unlist(L23_Vis_MR_L23_171_A[1]), unlist(L23_Vis_MR_L23_171_B[1]), unlist(L23_Vis_MR_L23_171_C[1]),
                                 unlist(L23_Vis_MR_L23_168_A[1]), unlist(L23_Vis_MR_L23_168_B[1]), unlist(L23_Vis_MR_L23_168_C[1]),
                                 unlist(L23_Vis_MR_L4_179_A[1]), unlist(L23_Vis_MR_L4_179_B[1]), unlist(L23_Vis_MR_L4_179_C[1]),
                                 unlist(L23_Vis_MA_L23_166_A[1]), unlist(L23_Vis_MA_L23_166_B[1]), unlist(L23_Vis_MA_L23_166_C[1]),
                                 unlist(L23_Vis_MA_L23_171_A[1]), unlist(L23_Vis_MA_L23_171_B[1]), unlist(L23_Vis_MA_L23_171_C[1]),
                                 unlist(L23_Vis_MA_L23_168_A[1]), unlist(L23_Vis_MA_L23_168_B[1]), unlist(L23_Vis_MA_L23_168_C[1]),
                                 unlist(L23_Vis_MA_L4_179_A[1]), unlist(L23_Vis_MA_L4_179_B[1]), unlist(L23_Vis_MA_L4_179_C[1])), 
                               nrow=3, ncol=8)
rownames(L23_Vis_sig_ABC_pvals_matrix) = c("A", "B", "C")
colnames(L23_Vis_sig_ABC_pvals_matrix) = c("MR_L23_166", "MR_L23_171", "MR_L23_168", "MR_L4_179",
                                   "MA_L23_166", "MA_L23_171", "MA_L23_168", "MA_L4_179")

write.table(L23_Vis_sig_ABC_pvals_matrix, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_Vis_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_merfish_pvals.txt", row.names=TRUE, col.names=TRUE, quote=F, sep="\t")
L23_Vis_sig_ABC_log10pvals_matrix <- -log10(L23_Vis_sig_ABC_pvals_matrix + 1 - 1)

#log2 odds ratio
L23_Vis_sig_ABC_log2OR_matrix <- matrix(c(unlist(L23_Vis_MR_L23_166_A[2]), unlist(L23_Vis_MR_L23_166_B[2]), unlist(L23_Vis_MR_L23_166_C[2]),
                                  unlist(L23_Vis_MR_L23_171_A[2]), unlist(L23_Vis_MR_L23_171_B[2]), unlist(L23_Vis_MR_L23_171_C[2]),
                                  unlist(L23_Vis_MR_L23_168_A[2]), unlist(L23_Vis_MR_L23_168_B[2]), unlist(L23_Vis_MR_L23_168_C[2]),
                                  unlist(L23_Vis_MR_L4_179_A[2]), unlist(L23_Vis_MR_L4_179_B[2]), unlist(L23_Vis_MR_L4_179_C[2]),
                                  unlist(L23_Vis_MA_L23_166_A[2]), unlist(L23_Vis_MA_L23_166_B[2]), unlist(L23_Vis_MA_L23_166_C[2]),
                                  unlist(L23_Vis_MA_L23_171_A[2]), unlist(L23_Vis_MA_L23_171_B[2]), unlist(L23_Vis_MA_L23_171_C[2]),
                                  unlist(L23_Vis_MA_L23_168_A[2]), unlist(L23_Vis_MA_L23_168_B[2]), unlist(L23_Vis_MA_L23_168_C[2]),
                                  unlist(L23_Vis_MA_L4_179_A[2]), unlist(L23_Vis_MA_L4_179_B[2]), unlist(L23_Vis_MA_L4_179_C[2])), 
                                nrow=3, ncol=8)
rownames(L23_Vis_sig_ABC_log2OR_matrix) = c("A", "B", "C")
colnames(L23_Vis_sig_ABC_log2OR_matrix) = c("MR_L23_166", "MR_L23_171", "MR_L23_168", "MR_L4_179",
                                    "MA_L23_166", "MA_L23_171", "MA_L23_168", "MA_L4_179")

round(L23_Vis_sig_ABC_log2OR_matrix,2)
col_palette2 <- colorRampPalette(c("white", "red", "darkred"))(n = 299)
#col_breaks = c(seq(-3, -1,length=100),  # for blue
#               seq(-0.99,0.99,length=100), #for white
#               seq(1, 3, length=100)) #for red



setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/L23_Vis_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_merfish_heatmap.eps")
heatmap.2(L23_Vis_sig_ABC_log10pvals_matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=col_palette2,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, key.title="Log10 p-value")
dev.off()

###ordering cell plotting by expression
#example_gene_plot_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6", dot_size=2)

#exp2_zScoreOrder_Trpc6 = exp2_sct_counts_Vis_glut_L23_scaled_dt[order(Trpc6), index]
exp2_zScoreOrder_WT_Trpc6  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT")][(order(Trpc6)), index]
exp2_zScoreOrder_KO_Trpc6  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO")][(order(Trpc6)), index]

exp2_zScoreOrder_WT_Spon1  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT")][(order(Spon1)), index]
exp2_zScoreOrder_KO_Spon1  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO")][(order(Spon1)), index]

exp2_zScoreOrder_WT_Tshz1  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="WT")][(order(Tshz1)), index]
exp2_zScoreOrder_KO_Tshz1  = exp2_sct_counts_Vis_glut_L23_scaled_dt[(index %in% exp2_L23_L$index) & (t.type.Under1Over1=="KO")][(order(Tshz1)), index]


exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1 <- copy(exp2_sct_counts_Vis_glut_L23_scaled_dt)
exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6 <- copy(exp2_sct_counts_Vis_glut_L23_scaled_dt)
exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1 <- copy(exp2_sct_counts_Vis_glut_L23_scaled_dt)

exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1$plot_color <- get("A_pal3")(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[, get("Spon1")], breaks=20))]
exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6$plot_color <- get("B_pal3")(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[, get("Trpc6")], breaks=20))]
exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[, get("Tshz1")], breaks=20))]


#
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1_WT.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_WT_Spon1){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,plot_color])
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1_KO.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_KO_Spon1){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Spon1[index==i,plot_color])
}
dev.off()
#
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6_WT.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_WT_Trpc6){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,plot_color])
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6_KO.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_KO_Trpc6){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Trpc6[index==i,plot_color])
}
dev.off()

#
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1_WT.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="WT", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_WT_Tshz1){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,plot_color])
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1_KO.eps")
plot(exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[(index %in% exp2_L23_L$index),center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[(index %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="KO", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_KO_Tshz1){
  points(exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,center_x],exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_Vis_glut_L23_scaled_dt_Tshz1[index==i,plot_color])
}
dev.off()



example_gene_plot_zScoreOrder_png_func <- function(example_gene, exp_data_table, cell.id, color_pal, color_column, xlim=limx_L23, ylim=limy_L23, output_file, dot_size=1.5){
  exp_data_table_example <- exp_data_table
  exp_data_table_example$plot_color <- get(color_pal)(20)[as.numeric(cut(exp_data_table_example[, get(example_gene)], breaks=20))]
  
  zScoreOrder_WT = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT")][(order(example_gene)), index]
  zScoreOrder_KO = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO")][(order(example_gene)), index]
  
  setEPS()
  postscript(paste0(output_file, "_WT.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="WT", xlab="x", ylab="y", asp=1)
  for(i in zScoreOrder_WT){
    points(exp_data_table_example[index==i,center_x],exp_data_table_example[index==i,center_y],pch = 20,cex=dot_size, col = exp_data_table_example[index==i,get(color_column)])
  }
  dev.off()
  
  setEPS()
  postscript(paste0(output_file, "_KO.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="KO", xlab="x", ylab="y", asp=1)
  for(i in zScoreOrder_KO){
    points(exp_data_table_example[index==i,center_x],exp_data_table_example[index==i,center_y],pch = 20,cex=dot_size, col = exp_data_table_example[index==i,get(color_column)])
  }
  dev.off()
  
  color.bar(get(color_pal)(20), min=signif(min(exp_data_table_example[, get(example_gene)]),2), max=signif(max(exp_data_table_example[, get(example_gene)]),2), nticks=2, title=example_gene)
}

example_gene_plot_zScoreOrder_eps_func <- function(example_gene, exp_data_table, cell.id, color_pal, color_column, xlim=limx_L23, ylim=limy_L23, output_file, dot_size=1.5){
  exp_data_table_example <- exp_data_table
  exp_data_table_example$plot_color <- get(color_pal)(20)[as.numeric(cut(exp_data_table_example[, get(example_gene)], breaks=20))]
  
  zScoreOrder_WT = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="WT")][(order(example_gene)), index]
  zScoreOrder_KO = exp_data_table_example[(index %in% cell.id) & (t.type.Under1Over1=="KO")][(order(example_gene)), index]
  
  setEPS()
  postscript(paste0(output_file, "_WT.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="WT", xlab="x", ylab="y", asp=1)
  for(i in zScoreOrder_WT){
    points(exp_data_table_example[index==i,center_x],exp_data_table_example[index==i,center_y],pch = 20,cex=dot_size, col = exp_data_table_example[index==i,get(color_column)])
  }
  dev.off()
  
  setEPS()
  postscript(paste0(output_file, "_KO.eps"))
  plot(exp_data_table_example[(index %in% cell.id),center_x],exp_data_table_example[(index %in% cell.id),center_y],pch = 20,cex=.3,col = "white",xlim=xlim, ylim=ylim,main="KO", xlab="x", ylab="y", asp=1)
  for(i in zScoreOrder_KO){
    points(exp_data_table_example[index==i,center_x],exp_data_table_example[index==i,center_y],pch = 20,cex=dot_size, col = exp_data_table_example[index==i,get(color_column)])
  }
  dev.off()
  
  color.bar(get(color_pal)(20), min=signif(min(exp_data_table_example[, get(example_gene)]),2), max=signif(max(exp_data_table_example[, get(example_gene)]),2), nticks=2, title=example_gene)
}

example_gene_plot_zScoreOrder_eps_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1", dot_size=2)
example_gene_plot_zScoreOrder_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6", dot_size=2)
example_gene_plot_zScoreOrder_eps_func(example_gene="Tshz1", exp_data_table=exp2_sct_counts_Vis_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="C_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_C_gene_Tshz1", dot_size=2)
