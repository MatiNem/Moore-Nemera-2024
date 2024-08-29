library(data.table)
library(dplyr)
library(ggplot2)

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
#TPMs
Pv_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_Mecp2KO_gene_TPMs_nondedup.txt")
Sst_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_Mecp2KO_gene_TPMs_nondedup.txt")
L4_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_Mecp2KO_gene_TPMs_nondedup.txt")
L5_TPM <- fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_Mecp2KO_gene_TPMs_nondedup.txt")

##gene body mCA
#gene body and flank methylation tables, using genic methylation from TSS+3kb to TES
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

#core MR genes
meta_MR_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/meta_MR_genes_geneColumn4_mm9.bed")
names(meta_MR_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")

meta_unchanged_genes = fread("HG_lab/Mati/GabelLab/genesets/meta_genes/MeCP2_unchanged_meta_genes_mm9.bed")
names(meta_unchanged_genes) = c("chrom", "start", "end", "gene", "num_transcripts", "strand")

#resampling function
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

#PV
Pv_TPM_metaMR_df <- data.frame(Pv_TPM[Gene %in% meta_MR_genes$gene, .(Gene, Pv_WT_TPM_avg)], row.names="Gene")
Pv_TPM_metaMR_df$Pv_WT_TPM_avg2 <- Pv_TPM_metaMR_df$Pv_WT_TPM_avg

Pv_TPM_metaUnchanged_df <- data.frame(Pv_TPM[Gene %in% meta_unchanged_genes$gene, .(Gene, Pv_WT_TPM_avg)], row.names="Gene")
Pv_TPM_metaUnchanged_df$Pv_WT_TPM_avg2 <- Pv_TPM_metaUnchanged_df$Pv_WT_TPM_avg

#Sst
Sst_TPM_metaMR_df <- data.frame(Sst_TPM[Gene %in% meta_MR_genes$gene, .(Gene, Sst_WT_TPM_avg)], row.names="Gene")
Sst_TPM_metaMR_df$Sst_WT_TPM_avg2 <- Sst_TPM_metaMR_df$Sst_WT_TPM_avg

Sst_TPM_metaUnchanged_df <- data.frame(Sst_TPM[Gene %in% meta_unchanged_genes$gene, .(Gene, Sst_WT_TPM_avg)], row.names="Gene")
Sst_TPM_metaUnchanged_df$Sst_WT_TPM_avg2 <- Sst_TPM_metaUnchanged_df$Sst_WT_TPM_avg

#L4
L4_TPM_metaMR_df <- data.frame(L4_TPM[Gene %in% meta_MR_genes$gene, .(Gene, L4_WT_TPM_avg)], row.names="Gene")
L4_TPM_metaMR_df$L4_WT_TPM_avg2 <- L4_TPM_metaMR_df$L4_WT_TPM_avg

L4_TPM_metaUnchanged_df <- data.frame(L4_TPM[Gene %in% meta_unchanged_genes$gene, .(Gene, L4_WT_TPM_avg)], row.names="Gene")
L4_TPM_metaUnchanged_df$L4_WT_TPM_avg2 <- L4_TPM_metaUnchanged_df$L4_WT_TPM_avg

#L5
L5_TPM_metaMR_df <- data.frame(L5_TPM[Gene %in% meta_MR_genes$gene, .(Gene, L5_WT_TPM_avg)], row.names="Gene")
L5_TPM_metaMR_df$L5_WT_TPM_avg2 <- L5_TPM_metaMR_df$L5_WT_TPM_avg

L5_TPM_metaUnchanged_df <- data.frame(L5_TPM[Gene %in% meta_unchanged_genes$gene, .(Gene, L5_WT_TPM_avg)], row.names="Gene")
L5_TPM_metaUnchanged_df$L5_WT_TPM_avg2 <- L5_TPM_metaUnchanged_df$L5_WT_TPM_avg


#example of one resampling
resamp1 <- resample(Pv_TPM_metaUnchanged_df, Pv_TPM_metaMR_df)

#set number of resamplings
num_resamp = 100
#call matrix that you will populate with resamplings
Pv_resamplings = matrix(nrow=nrow(Pv_TPM_metaMR_df),ncol=num_resamp)
set.seed(1)
#for loop to populate the resampling matrix
for(n in 1:num_resamp){
  Pv_resamplings[,n] = resample(Pv_TPM_metaUnchanged_df, Pv_TPM_metaMR_df)
}

#call matrix that you will populate with resamplings
Sst_resamplings = matrix(nrow=nrow(Sst_TPM_metaMR_df),ncol=num_resamp)
set.seed(1)
#for loop to populate the resampling matrix
for(n in 1:num_resamp){
  Sst_resamplings[,n] = resample(Sst_TPM_metaUnchanged_df, Sst_TPM_metaMR_df)
}

#call matrix that you will populate with resamplings
L4_resamplings = matrix(nrow=nrow(L4_TPM_metaMR_df),ncol=num_resamp)
set.seed(1)
#for loop to populate the resampling matrix
for(n in 1:num_resamp){
  L4_resamplings[,n] = resample(L4_TPM_metaUnchanged_df, L4_TPM_metaMR_df)
}

#call matrix that you will populate with resamplings
L5_resamplings = matrix(nrow=nrow(L5_TPM_metaMR_df),ncol=num_resamp)
set.seed(1)
#for loop to populate the resampling matrix
for(n in 1:num_resamp){
  L5_resamplings[,n] = resample(L5_TPM_metaUnchanged_df, L5_TPM_metaMR_df)
}

write.table(Pv_resamplings, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_PvTPM_resamplings_100.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(Sst_resamplings, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_SstTPM_resamplings_100.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(L4_resamplings, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_L4TPM_resamplings_100.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(L5_resamplings, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_L5TPM_resamplings_100.txt", quote=F, row.names=F, col.names=F, sep="\t")

#reading in resamplings
Pv_resamplings=fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_PvTPM_resamplings_100.txt", header=FALSE)
Sst_resamplings=fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_SstTPM_resamplings_100.txt", header=FALSE)
L4_resamplings=fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_L4TPM_resamplings_100.txt", header=FALSE)
L5_resamplings=fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/resampled_genes/meta_genes/meta_MR_L5TPM_resamplings_100.txt", header=FALSE)


gene_body_mCA_INTACT_resamp <- rbind(
  cbind(L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% L4_resamplings[[1]]], gene_list="resamp", mCA_label="L4"),
  cbind(L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% meta_MR_genes[, gene]], gene_list="meta_MR", mCA_label="L4"),
  cbind(L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% L5_resamplings[[1]]], gene_list="resamp", mCA_label="L5"),
  cbind(L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% meta_MR_genes[, gene]], gene_list="meta_MR", mCA_label="L5"),
  cbind(PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% Pv_resamplings[[1]]], gene_list="resamp", mCA_label="PV"),
  cbind(PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% meta_MR_genes[, gene]], gene_list="meta_MR", mCA_label="PV"),
  cbind(SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% Sst_resamplings[[1]]], gene_list="resamp", mCA_label="SST"),
  cbind(SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene %in% meta_MR_genes[, gene]], gene_list="meta_MR", mCA_label="SST")
)
gene_body_mCA_INTACT_resamp = gene_body_mCA_INTACT_resamp %>% mutate(gene_list = factor(gene_list, levels=c("resamp", "meta_MR")))
gene_body_mCA_INTACT_resamp = gene_body_mCA_INTACT_resamp %>% mutate(mCA_label = factor(mCA_label, levels=c("L4", "L5", "PV", "SST")))

ggplot(gene_body_mCA_INTACT_resamp, aes(x = mCA_label, y = as.numeric(gene_methylation_corrected), fill=gene_list))+
  ggtitle("")+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,0.15))+
  ylab("Gene body mCA/CA") + xlab("")+
  scale_fill_manual(name = "", values = c("resamp"="gray", "meta_MR"="red")) +
  facet_grid(.~mCA_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/INTACT.core.MR.resamp.gene.body.INTACT.WT.KO.deep.mCA.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/INTACT.core.MR.resamp.gene.body.INTACT.WT.KO.deep.mCA.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

options(scipen=0)
wilcox.test(gene_body_mCA_INTACT_resamp[(gene_list=="resamp") & (mCA_label=="L4"), as.numeric(gene_methylation_corrected)], gene_body_mCA_INTACT_resamp[(gene_list=="meta_MR") & (mCA_label=="L4"), as.numeric(gene_methylation_corrected)])$p.value #p=1.909372e-42
wilcox.test(gene_body_mCA_INTACT_resamp[(gene_list=="resamp") & (mCA_label=="L5"), as.numeric(gene_methylation_corrected)], gene_body_mCA_INTACT_resamp[(gene_list=="meta_MR") & (mCA_label=="L5"), as.numeric(gene_methylation_corrected)])$p.value #p=2.682717e-42
wilcox.test(gene_body_mCA_INTACT_resamp[(gene_list=="resamp") & (mCA_label=="PV"), as.numeric(gene_methylation_corrected)], gene_body_mCA_INTACT_resamp[(gene_list=="meta_MR") & (mCA_label=="PV"), as.numeric(gene_methylation_corrected)])$p.value #p=1.15077e-36
wilcox.test(gene_body_mCA_INTACT_resamp[(gene_list=="resamp") & (mCA_label=="SST"), as.numeric(gene_methylation_corrected)], gene_body_mCA_INTACT_resamp[(gene_list=="meta_MR") & (mCA_label=="SST"), as.numeric(gene_methylation_corrected)])$p.value #p=5.052809e-41


