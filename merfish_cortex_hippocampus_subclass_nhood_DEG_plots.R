#install.packages("Seurat")
install.packages('remotes')
#install.packages('Seurat', repos = c('https://satijalab.r-universe.dev', 'https://cloud.r-project.org'))
remotes::install_version(package = 'Seurat', version = package_version('4.3.0.1'))
library(data.table)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(reshape2)
library(edgeR)
library(scuttle)
library(scran)
library(tibble)
#significance symbol function
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

#Gene lists for analyses
gene_panel_classes = fread("HG_lab/Vizgen/Gene_lists_for_analysis/MERFISH_gene_panel_classes.txt")
unchanged_genes = gene_panel_classes[(MR_Meta==FALSE) & (MA_meta==FALSE) & (Rec_cellspec_MR==FALSE) &
                                       (PV_MR==FALSE) & (PV_MA==FALSE) & (PV_OTHERCELLTYPE==FALSE) &
                                       (SST_MR==FALSE) & (SST_MA==FALSE) & (SST_OTHERCELLTYPE==FALSE) &
                                       (L5_MR==FALSE) & (L5_MA==FALSE) & (L5_OTHERCELLTYPE==FALSE) &
                                       (L4_MR==FALSE) & (L4_MA==FALSE) & (L4_OTHERCELLTYPE==FALSE), Gene]
nonMR_genes = gene_panel_classes[(MR_Meta==FALSE) & (Rec_cellspec_MR==FALSE) &
                                   (PV_MR==FALSE) & (PV_OTHERCELLTYPE==FALSE) &
                                   (SST_MR==FALSE) & (SST_OTHERCELLTYPE==FALSE) &
                                   (L5_MR==FALSE) & (L5_OTHERCELLTYPE==FALSE) &
                                   (L4_MR==FALSE) & (L4_OTHERCELLTYPE==FALSE), Gene]
unchanged_nonShort_genes = gene_panel_classes[(MR_Meta==FALSE) & (MA_meta==FALSE) & (Rec_cellspec_MR==FALSE) &
                                                (PV_MR==FALSE) & (PV_MA==FALSE) & (PV_OTHERCELLTYPE==FALSE) &
                                                (SST_MR==FALSE) & (SST_MA==FALSE) & (SST_OTHERCELLTYPE==FALSE) &
                                                (L5_MR==FALSE) & (L5_MA==FALSE) & (L5_OTHERCELLTYPE==FALSE) &
                                                (L4_MR==FALSE) & (L4_MA==FALSE) & (L4_OTHERCELLTYPE==FALSE) & Short_lowmC_Marker==FALSE, Gene]
short_lowmC_genes = gene_panel_classes[Short_lowmC_Marker==TRUE,Gene]
any_meta_MR_genes = gene_panel_classes[(MR_Meta==TRUE) | (Rec_cellspec_MR==TRUE), Gene]

#experiment 1 metadata 
cell_metadata <- fread("HG_lab/Vizgen/experiment1_2022-08-18/merfish_output/cell_metadata.csv")
names(cell_metadata)[1] = "X"
#experiment gene counts
gene_counts <- fread("HG_lab/Vizgen/experiment1_2022-08-18/merfish_output/cell_by_gene.csv")
names(gene_counts)[1] = "X"
gene_names = setdiff(names(gene_counts), "X")
cell_metadata_gene_counts = left_join(x=cell_metadata, y=gene_counts, by="X")%>% 
  mutate(total_counts = rowSums(.[,10:ncol(.)]))

#experiment 2 metadata 
cell_metadata2 <- fread("HG_lab/Vizgen/experiment2_2022-09-08/merfish_output/cell_metadata.csv")
names(cell_metadata2)[1] = "X"
#experiment gene counts
gene_counts2 <- fread("HG_lab/Vizgen/experiment2_2022-09-08/merfish_output/cell_by_gene.csv")
names(gene_counts2)[1] = "X"
gene_names2 = setdiff(names(gene_counts2), "X")
cell_metadata_gene_counts2 = left_join(x=cell_metadata2, y=gene_counts2, by="X")%>% 
  mutate(total_counts = rowSums(.[,10:ncol(.)]))


#experiment 6 metadata 
cell_metadata6 <- fread("HG_lab/Vizgen/experiment6_2023-01-23/merfish_output/cell_metadata.csv")
names(cell_metadata6)[1] = "X"
#experiment gene counts
gene_counts6 <- fread("HG_lab/Vizgen/experiment6_2023-01-23/merfish_output/cell_by_gene.csv")
names(gene_counts6)[1] = "X"
gene_names6 = setdiff(names(gene_counts6), "X")
cell_metadata_gene_counts6 = left_join(x=cell_metadata6, y=gene_counts6, by="X")%>% 
  mutate(total_counts = rowSums(.[,10:ncol(.)]))

#experiment 7 metadata 
cell_metadata7 <- fread("HG_lab/Vizgen/experiment7_2023-01-27/merfish_output/cell_metadata.csv")
names(cell_metadata7)[1] = "X"
#experiment gene counts
gene_counts7 <- fread("HG_lab/Vizgen/experiment7_2023-01-27/merfish_output/cell_by_gene.csv")
names(gene_counts7)[1] = "X"
gene_names7 = setdiff(names(gene_counts7), "X")
cell_metadata_gene_counts7 = left_join(x=cell_metadata7, y=gene_counts7, by="X")%>% 
  mutate(total_counts = rowSums(.[,10:ncol(.)]))

#experiment metadata 
WT_exp_cell_metadata <- fread("HG_lab/Vizgen/Experiment3_WT_coronal_experiment_2022-12-13/merfish_output/cell_metadata.csv")
names(WT_exp_cell_metadata)[1] = "X"

#experiment gene counts
WT_exp_gene_counts <- fread("HG_lab/Vizgen/Experiment3_WT_coronal_experiment_2022-12-13/merfish_output/cell_by_gene.csv")
names(WT_exp_gene_counts)[1] = "X"
WT_exp_gene_names = setdiff(names(WT_exp_gene_counts), "X")
WT_exp_cell_metadata_gene_counts = left_join(x=WT_exp_cell_metadata, y=WT_exp_gene_counts, by="X")%>% 
  mutate(total_counts = rowSums(.[,10:ncol(.)]))


#Yao 2021 metadata
yao_meta_dt = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/metadata.csv", na.strings=c("", "NA"))
yao_meta_dt = yao_meta_dt[!is.na(cluster_label), ]
#subclass assignments of cluster labels from Yao 2021
cluster_subclass = unique(yao_meta_dt[, .(cluster_label, subclass_label)])
yao2021_DE = readRDS("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/comb.de.genes.Zizhen.Yao.rda")


CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
CTX_HIP_annot[, cl := as.character(cl)]

#Seurat object of MeCP2 KO/+ MERFISH cells and their labels
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_pred0.2_Under1Over1_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_pred0.2_Under1Over2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

#average expressions in 100vol, 300 count filter
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

#average expressions in 300vol, 600 count filter
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT = subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT = subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

#INTACT targets
creLine_dist = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/Yao_S3_CreLine_distribution.csv")
creLine_dist[, cl := as.character(cl)]
creLine_dist_size = left_join(x=creLine_dist, y=CTX_HIP_annot[, .(cl, cluster_label, cluster_size, X10X_cells.size, Smartseq_cells.size, subclass_label)], by=c("cl"))
creLine_dist_size$cluster_size_Rbp4 = creLine_dist_size[, as.numeric(`Rbp4-Cre_KL100`) * cluster_size]

creLine_dist_size$cluster_size_Rbp4 = creLine_dist_size[, as.numeric(`Rbp4-Cre_KL100`) * cluster_size]
Rbp4_sum = sum(creLine_dist_size[, cluster_size_Rbp4], na.rm=TRUE)
creLine_dist_size[, cluster_size_Rbp4_prop := cluster_size_Rbp4/Rbp4_sum]

creLine_dist_size$cluster_size_Nr5a1 = creLine_dist_size[, as.numeric(`Nr5a1-Cre`) * cluster_size]
Nr5a1_sum = sum(creLine_dist_size[, cluster_size_Nr5a1], na.rm=TRUE)
creLine_dist_size[, cluster_size_Nr5a1_prop := cluster_size_Nr5a1/Nr5a1_sum]
Rbp4_L5_targets = creLine_dist_size[cluster_size_Rbp4_prop >= 0.03, cluster_label.x]

Nr5a1_L4_targets = creLine_dist_size[cluster_size_Nr5a1_prop >= 0.03, cluster_label.x]
#Nr5a1_L4_targets = c("178_L4 IT CTX","179_L4 IT CTX","180_L4 IT CTX", "181_L4 IT CTX")

#

zscore_exp_acrossSubclass_func <- function(gene.by.subtype.avg.counts, subclass_names){
  subclass_subtypes = cluster_subclass[subclass_label%in%subclass_names, cluster_label]
  subclass_subtypes_in_data = intersect(colnames(gene.by.subtype.avg.counts), subclass_subtypes)
  if(length(subclass_subtypes_in_data) > 1){
    gene.by.subtype.avg.counts.subclass = gene.by.subtype.avg.counts[, c(subclass_subtypes_in_data)]
    gene.by.subtype.avg.counts.subclass.scale = scale(t(gene.by.subtype.avg.counts.subclass))
    gene.by.subtype.avg.counts.subclass.scale = data.table(gene.by.subtype.avg.counts.subclass.scale, keep.rownames = "subtype")
    gene.by.subtype.avg.counts.subclass.scale  = data.table(melt(gene.by.subtype.avg.counts.subclass.scale, id.vars=c("subtype")))
    names(gene.by.subtype.avg.counts.subclass.scale) = c("subtype", "gene", "zScore")
  }
  else{
    gene.by.subtype.avg.counts.subclass.scale = logical(0)
  }
  return(gene.by.subtype.avg.counts.subclass.scale)
}
scaled.DEG.subclass.func <- function(subclass_label, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table){
  gene.by.subtype.avg.counts.scale = zscore_exp_acrossSubclass_func(gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, subclass_names=subclass_label)
  if(length(gene.by.subtype.avg.counts.scale)>0){
    gene.by.subtype.avg.counts.scale  = left_join(x=gene.by.subtype.avg.counts.scale , y=de.results.table, by=c("subtype", "gene"))
    cluster.DEG.dt = setNames(data.table(matrix(nrow = 0, ncol = 3)), c("cl", "gene", "dir"))
    for(i in unique(binary_comparison_table[(subclass_label.x%in%subclass_label), cl.x])){
      for(j in binary_comparison_table[(subclass_label.x%in%subclass_label) & (cl.x==i), comparison]){
        up.genes = names(yao2021_DE[[c(j, "up.genes")]])
        down.genes= names(yao2021_DE[[c(j, "down.genes")]])
        non_DEGs = setdiff(gene.by.subtype.avg.counts.scale[, gene], c(up.genes, down.genes))
        if(length(up.genes)>0){
          cluster.DEG.dt = data.table(rbind(
            cluster.DEG.dt,
            cbind(cl=i, gene=up.genes, dir="High")))
        }
        if(length(down.genes)>0){
          cluster.DEG.dt = data.table(rbind(
            cluster.DEG.dt,
            cbind(cl=i, gene=down.genes, dir="Low")))
        }
        if(length(non_DEGs)>0){
          cluster.DEG.dt = data.table(rbind(
            cluster.DEG.dt,
            cbind(cl=i, gene=non_DEGs, dir="Non-DEG")))
        }
      }
    }
    cluster.DEG.dt = left_join(unique(cluster.DEG.dt), CTX_HIP_annot[,.(cl, cluster_label, subclass_label, neighborhood_label)], by=c("cl"))
    names(cluster.DEG.dt)[4] = "subtype"
    
    scaled.DEG.table = left_join(x=cluster.DEG.dt, y=gene.by.subtype.avg.counts.scale, by=c("subtype", "gene"))
    scaled.DEG.table <- left_join(scaled.DEG.table, avg.exp.table, by=c("subtype", "gene"))
    scaled.DEG.table$avgExp.subclass.med <- median(scaled.DEG.table[, as.numeric(avgExp)], na.rm=TRUE)
  }
  else{
    scaled.DEG.table = logical(0)
  }
  return(scaled.DEG.table)
}
#binary comparisons
yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2 <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Yao2021_pairwise_DEGs/yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.comparisons.csv")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
scaled.DEG.Pvalb <- scaled.DEG.subclass.func(subclass_label="Pvalb", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2, avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)
#scaled.DEG.Pvalb.old <- scaled.DEG.subclass.func(subclass_label="Pvalb", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2, avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)


group.scaled.DEG.func2 <- function(grouping="neighborhood", nhood_label=NULL, subclass_list=NULL, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table){
  scaled.DEG.table.group = setNames(data.table(matrix(nrow = 0, ncol = 15)), c("cl", "gene", "dir", "subtype", "subclass_label", "neighborhood_label", "zScore", "logFC", "logCPM", "F", "PValue", "FDR", "p.adj.BH.all", "avgExp", "avgExp.subclass.med"))
  if(grouping=="neighborhood"){
    for(i in unique(CTX_HIP_annot[neighborhood_label==nhood_label, subclass_label])){
      scaled.DEG.table.subclass = scaled.DEG.subclass.func(subclass_label=i, gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
      if(length(scaled.DEG.table.subclass) > 0){
        scaled.DEG.table.group = rbind(scaled.DEG.table.group, scaled.DEG.table.subclass)
      }
      else{
        next
      }
    }
  }
  if(grouping=="subclass_list"){
    for(i in unique(subclass_list)){
      scaled.DEG.table.subclass = scaled.DEG.subclass.func(subclass_label=i, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table)
      if(length(scaled.DEG.table.subclass) > 0){
        scaled.DEG.table.group = rbind(scaled.DEG.table.group, scaled.DEG.table.subclass)
      }
      else{
        next
      }
    }
  }
  return(scaled.DEG.table.group)
}

#scaled.DEG.table.MGE.test  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "MGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2, avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)
#scaled.DEG.table.CGE.test  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "CGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2, avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)

acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table <- function(scaled.DEG.table){
  subtype_zscore_MR_dt = data.table(rbind(
    cbind(scaled.DEG.table[(dir=="Non-DEG") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[PV_MR==TRUE,Gene])], gene_class="PV MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[PV_MR==TRUE,Gene])], gene_class="PV MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[SST_MR==TRUE,Gene])], gene_class="SST MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[SST_MR==TRUE,Gene])], gene_class="SST MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[L4_MR==TRUE,Gene])], gene_class="L4 MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[L4_MR==TRUE,Gene])], gene_class="L4 MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[L5_MR==TRUE,Gene])], gene_class="L5 MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[L5_MR==TRUE,Gene])], gene_class="L5 MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% short_lowmC_genes)], gene_class="Short low mC"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% short_lowmC_genes)], gene_class="Short low mC")))
  subtype_zscore_MR_dt_nonDEG = subtype_zscore_MR_dt[(dir=="Non-DEG")]
  subtype_zscore_MR_dt_up = subtype_zscore_MR_dt[(dir=="High") & (avgExp > avgExp.subclass.med),]
  subtype_zscore_MR_dt_down = subtype_zscore_MR_dt[(dir=="Low") & (avgExp < avgExp.subclass.med),]
  subtype_zscore_MR_dt = rbind(subtype_zscore_MR_dt_nonDEG, subtype_zscore_MR_dt_up, subtype_zscore_MR_dt_down)
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "PV MR", "SST MR", "L4 MR", "L5 MR", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(subtype_zscore_MR_dt[!(is.na(gene_class))])
}

all.nhoods.func3 <- function(gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table){
  scaled.DEG.table.L2.3.IT  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "L2/3 IT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  L2.3.IT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.L2.3.IT)
  
  scaled.DEG.table.L4.5.6.IT.Car3  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "L4/5/6 IT Car3", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  L4.5.6.IT.Car3.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.L4.5.6.IT.Car3)
  
  scaled.DEG.table.MGE  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "MGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  MGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.MGE)
  
  scaled.DEG.table.CGE  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "CGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  CGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.CGE)
  
  scaled.DEG.table.PT  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "PT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  PT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.PT)
  
  scaled.DEG.table.NP.CT.L6b  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "NP/CT/L6b", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  NP.CT.L6b.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.NP.CT.L6b)
  
  scaled.DEG.table.DG.SUB.CA  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "DG/SUB/CA", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  DG.SUB.CA.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.DG.SUB.CA)
  
  all.nhoods.dt = data.table(rbind(
    cbind(L2.3.IT.dt, neighborhood="L2/3 IT"),
    cbind(L4.5.6.IT.Car3.dt, neighborhood="L4/5/6 IT Car3"),
    cbind(MGE.dt, neighborhood="MGE"),
    cbind(CGE.dt, neighborhood="CGE"),
    cbind(PT.dt, neighborhood="PT"),
    cbind(NP.CT.L6b.dt, neighborhood="NP/CT/L6b"),
    cbind(DG.SUB.CA.dt, neighborhood="DG/SUB/CA")))
  
  all.nhoods.dt = all.nhoods.dt %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE","MGE","L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "PV MR", "SST MR", "L4 MR", "L5 MR", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(all.nhoods.dt)
}

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table = all.nhoods.func3(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                                                                                               de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, 
                                                                                               binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                                                                                               avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)

#write.csv(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table, file="HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.table(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table, file="HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.txt", quote=F, row.names=F, sep="\t")


subclass_low_high_nonDEG_violPlot_func <- function(subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c(paste(subclass_abbr, "MR"))) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==5){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"), "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  count = subtype_logfc_data_table[!(is.na(gene_class))] %>% 
    group_by(gene_class, dir) %>% 
    summarise(count = n()) %>% data.table
  
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
    ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observations
    geom_violin(fill=plot_color, trim=FALSE)+
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = 6, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = 6, height = 6, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  return(vplot)
}

Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 = subclass_low_high_nonDEG_violPlot_func(subclass_list=c("Pvalb"),
                                                                                                                           subclass_abbr="PV",
                                                                                                                           low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                                                           plot_color=colors_subclass["Pvalb"],
                                                                                                                           plot_save=TRUE,
                                                                                                                           plot_title="Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                                           plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.PvMR.recurrent.MR.logFC.violPlot2",
                                                                                                                           return_table=TRUE)



Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 = subclass_low_high_nonDEG_violPlot_func(subclass_list=c("Sst"),
                                                                                                                             subclass_abbr="SST",
                                                                                                                             low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                                                             plot_color=colors_subclass["Sst"],
                                                                                                                             plot_save=TRUE,
                                                                                                                             plot_title="Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                                             plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.SstMR.recurrent.MR.logFC.violPlot2",
                                                                                                                             return_table=TRUE)

L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 = subclass_low_high_nonDEG_violPlot_func(subclass_list=Rbp4_L5_targets,
                                                                                                                        custom_subtype_list=TRUE,
                                                                                                                        subclass_abbr="L5",
                                                                                                                        low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                                                        ymin=-2.5, ymax=2.5,
                                                                                                                        plot_color="#50B2AD",
                                                                                                                        plot_save=TRUE,
                                                                                                                        plot_title="L5\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                                        plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Rbp4.L5.targets/L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.L5MR.recurrent.MR.logFC.violPlot2",
                                                                                                                        return_table=TRUE)

subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs <- function(subtype_table, subclass_abbr){
  pval1 <- wilcox.test(x=subtype_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_table[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig1 <- sig_function(pval1)
  
  pval2 <- wilcox.test(x=subtype_table[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")), logFC],
                       y=subtype_table[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig2 <- sig_function(pval2)
  
  pval3 <- wilcox.test(x=subtype_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_table[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig3 <- sig_function(pval3)
  
  pval4 <- wilcox.test(x=subtype_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig4 <- sig_function(pval4)
  
  pval5 <- wilcox.test(x=subtype_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
                       y=subtype_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig5 <- sig_function(pval5)
  
  pval6 <- wilcox.test(x=subtype_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig6 <- sig_function(pval6)
  
  sigs <- data.table(rbind(cbind(compar=paste("Non-MR, non-MA Non-DEG;", subclass_abbr,"MR Low"), sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar=paste(subclass_abbr,"MR Low;", subclass_abbr,"MR High"), sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", subclass_abbr,"MR High"), sig_symbol=sig3, wilcox.pval=pval3),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR Low"), sig_symbol=sig4, wilcox.pval=pval4),
                           cbind(compar=paste("Any recurrent MR Low;", "Any recurrent MR High"), sig_symbol=sig5, wilcox.pval=pval5),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR High"), sig_symbol=sig6, wilcox.pval=pval6)))
  return(sigs)
}

View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2, subclass_abbr="PV"))
View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2, subclass_abbr="SST"))
View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2, subclass_abbr="L5"))

Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 %>% 
  group_by(gene_class, dir) %>% 
  summarise(count = n()) %>% data.table

Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 %>% 
  group_by(gene_class, dir) %>% 
  summarise(count = n()) %>% data.table

L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_2 %>% 
  group_by(gene_class, dir) %>% 
  summarise(count = n()) %>% data.table

#
subclass.MRlow.MRhigh.summ <- all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table %>% group_by(subclass_label) %>% 
  summarise(
    mean.logfc.nonMRnonMA = mean(as.numeric(logFC)[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG")], na.rm=TRUE),
    mean.logfc.MRlow = mean(as.numeric(logFC)[(gene_class=="Any recurrent MR") & (dir=="Low")], na.rm=TRUE),
    mean.logfc.MRhigh = mean(as.numeric(logFC)[(gene_class=="Any recurrent MR") & (dir=="High")], na.rm=TRUE),
    Low = mean.logfc.MRlow-mean.logfc.nonMRnonMA,
    High = mean.logfc.MRhigh-mean.logfc.nonMRnonMA
    ) %>% data.table



subclass.MRlow.MRhigh.summ.melt <- data.table(melt(subclass.MRlow.MRhigh.summ, id.vars = "subclass_label", measure.vars=c("Low", "High")))

ggplot(subclass.MRlow.MRhigh.summ.melt[is.finite(value)],aes(x=subclass_label,y=variable)) + geom_tile(aes(fill=value),color = "white") +
  #Creating legend
  guides(fill=guide_colorbar(title="Fold difference of mean fold changes", title.position="top")) +
  #Creating color range
  scale_fill_gradientn(colors=c("white","red"),guide="colorbar") +
  #scale_fill_gradientn(limits = c(0,20), colors=c("white","red"),guide="colorbar") +
  #geom_text(aes(label = fc), color = "black", size = 6)+
  theme_bw()+
  #Rotating labels
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.text.x=element_text(size=10, angle=90), axis.text.y=element_text(size=10), axis.ticks = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/subclasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMRlow.anyMetaMRhigh.vs.nonDEG.nonMR.nonMA.foldDiffFC.heatmap.png", width = 10, height = 5, dpi = 300, units = "in", device='png')


#violin plot colorless, dots on it
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_func <- function(subclass_list, custom_subtype_list=FALSE, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  count = subtype_logfc_data_table[!(is.na(gene_class))] %>% 
    group_by(gene_class, dir) %>% 
    summarise(count = n()) %>% data.table
  
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
    ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observations
    geom_violin(trim=FALSE)+
    geom_jitter(size=0.6, alpha=0.9, color=plot_color) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = 4, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  return(vplot)
}

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_func <- function(subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c(paste(subclass_abbr, "MR"))) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"))))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  count = subtype_logfc_data_table[!(is.na(gene_class))] %>% 
    group_by(gene_class, dir) %>% 
    summarise(count = n()) %>% data.table
  
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
    ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observations
    geom_violin(trim=FALSE)+
    geom_jitter(size=0.6, alpha=0.9, color=plot_color) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = 4, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  return(vplot)
}

#colorless violins, colored dot plot for figures for manuscript
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_func(subclass_list=c("Pvalb"),
                                       low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                       plot_color=colors_subclass["Pvalb"],
                                       plot_save=TRUE,
                                       plot_title="Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                       plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot",
                                       return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_func(subclass_list=c("Pvalb"),
                                                     subclass_abbr="PV",
                                                     low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                     plot_color=colors_subclass["Pvalb"],
                                                     plot_save=TRUE,
                                                     plot_title="Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                     plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.PvMR.logFC.viol.dotPlot",
                                                     return_table=FALSE)
#Sst
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_func(subclass_list=c("Sst"),
                                                     low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                     ymin=-3, ymax=3,
                                                     plot_color=colors_subclass["Sst"],
                                                     plot_save=TRUE,
                                                     plot_title="Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                     plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot",
                                                     return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_func(subclass_list=c("Sst"),
                                                    subclass_abbr="SST",
                                                    low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                    ymin=-3, ymax=3,
                                                    plot_color=colors_subclass["Sst"],
                                                    plot_save=TRUE,
                                                    plot_title="Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                    plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.SstMR.logFC.viol.dotPlot",
                                                    return_table=FALSE)
#Rbp4 L5 targets
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_func(subclass_list=Rbp4_L5_targets,
                                       custom_subtype_list=TRUE,
                                       low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                       ymin=-2.5, ymax=2.5,
                                       plot_color="#50B2AD",
                                       plot_save=TRUE,
                                       plot_title="L5\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                       plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Rbp4.L5.targets/L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot",
                                       return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_func(subclass_list=Rbp4_L5_targets,
                                                     custom_subtype_list=TRUE,
                                                     subclass_abbr="L5",
                                                     low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                     ymin=-2.5, ymax=2.5,
                                                     plot_color="#50B2AD",
                                                     plot_save=TRUE,
                                                     plot_title="L5\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                     plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Rbp4.L5.targets/L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.L5MR.logFC.viol.dotPlot",
                                                     return_table=FALSE)


#colors subtypes, subclasses, neighborhoods for ggplots
name_vec_subtype = unique(CTX_HIP_annot$cluster_label)
col_vec_subtype = unique(CTX_HIP_annot$cluster_color)
colors_subtype = setNames(col_vec_subtype, name_vec_subtype)

name_vec_subclass = unique(CTX_HIP_annot$subclass_label)
col_vec_subclass = unique(CTX_HIP_annot$subclass_color)
colors_subclass = setNames(col_vec_subclass, name_vec_subclass)

name_vec_nhood = unique(CTX_HIP_annot$neighborhood_label)
col_vec_nhood = unique(CTX_HIP_annot$neighborhood_color)
colors_nhood = setNames(col_vec_nhood, name_vec_nhood)

#neighborhood agg
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table <- rbind(
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG")]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Any recurrent MR") & (dir=="Low")]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Any recurrent MR") & (dir=="High")])
  )

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE","MGE","L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA")))



all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table, neighborhood, dir, gene_class) %>%
  summarise(
    IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
    outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
    outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[2] +1.5*IQR
  ) %>% data.table

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR <-left_join(x=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table,y=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR, by=c("neighborhood", "dir", "gene_class"))

ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table, aes(x = neighborhood, y = as.numeric(logFC), color=subtype))+
  ggtitle("")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  #separate by gene class and direction of change
  facet_grid(.~neighborhood + gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  #trim=TRUE means the data is trimmed to fit the range of observations
  geom_violin(trim=TRUE, color="black", fill=NA)+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[, subtype], cluster_label]))]) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.viol.dotPlot.png", width = 15.5, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.viol.dotPlot.eps", width = 15.5, height = 6, dpi = 300, units = "in", device=cairo_ps, fallback_resolution = 300)


ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR[(as.numeric(logFC)>=as.numeric(outlier_min)) & (as.numeric(logFC)<=as.numeric(outlier_max))], aes(x = neighborhood, y = as.numeric(logFC), color=subtype))+
  ggtitle("")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  #separate by gene class and direction of change
  facet_grid(.~neighborhood + gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  #trim=TRUE means the data is trimmed to fit the range of observations
  geom_violin(trim=TRUE, color="black", fill=NA)+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[, subtype], cluster_label]))]) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.png", width = 15.5, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.eps", width = 15.5, height = 6, dpi = 300, units = "in", device=cairo_ps, fallback_resolution = 300)


###color dots in violin plot by subtype, any meta MR genes
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func <- function(subclass_list, custom_subtype_list=FALSE, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_width=3, plot_height=5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  
  subtype_logfc_data_IQR <- group_by(subtype_logfc_data_table, dir, gene_class) %>%
    summarise(
      IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
      outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
      outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[2] +1.5*IQR
    ) %>% data.table
  
  subtype_logfc_data_table <- left_join(x=subtype_logfc_data_table,y=subtype_logfc_data_IQR, by=c("dir", "gene_class"))
  subtype_logfc_data_table <- subtype_logfc_data_table[(as.numeric(logFC) >= outlier_min) & (as.numeric(logFC) <= outlier_max)]
  
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #print(unique(subtype_logfc_data_table[, subtype]))
  plot_color2 = plot_color[subtype_reduced_list]
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table, aes(x = dir, y = as.numeric(logFC), color=subtype))+
    #ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observation
    #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
    scale_color_manual(name="", values = plot_color2)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    geom_violin(trim=TRUE, color="black", fill=NA)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())+
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  
  pval1 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig1 <- sig_function(pval1)
  pval2 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig2 <- sig_function(pval2)
  pval3 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig3 <- sig_function(pval3)
  sigs <- data.table(rbind(cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR Low"), sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar=paste("Any recurrent MR Low;", "Any recurrent MR High"), sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR High"), sig_symbol=sig3, wilcox.pval=pval3)))
  return(sigs)
  return(vplot)
}


subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=c("Pvalb"),
                                                     low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                     ymin=-1.2, ymax=1.2,
                                                     plot_color=colors_subtype,
                                                     plot_width=3, plot_height=5,
                                                     plot_save=TRUE,
                                                     plot_title="",
                                                     plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors",
                                                     return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=c("Sst"),
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=1.2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=3, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors",
                                                                  return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=Rbp4_L5_targets,
                                                                  custom_subtype_list=TRUE,
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=1.2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=3, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Rbp4.L5.targets.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors",
                                                                  return_table=FALSE)
#smaller
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=c("Pvalb"),
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=1.2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=c("Sst"),
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=1.2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func(subclass_list=Rbp4_L5_targets,
                                                                  custom_subtype_list=TRUE,
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=1.2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Rbp4.L5.targets.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)
###color dots in violin plot by subtype, INTACT MR genes
subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func <- function(subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_width=3, plot_height=5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c(paste(subclass_abbr, "MR"))) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_logfc_data_IQR <- group_by(subtype_logfc_data_table, dir, gene_class) %>%
    summarise(
      IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
      outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
      outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[2] +1.5*IQR
    ) %>% data.table
  
  subtype_logfc_data_table <- left_join(x=subtype_logfc_data_table,y=subtype_logfc_data_IQR, by=c("dir", "gene_class"))
  subtype_logfc_data_table <- subtype_logfc_data_table[(as.numeric(logFC) >= outlier_min) & (as.numeric(logFC) <= outlier_max)]
  
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }

  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"))))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #print(unique(subtype_logfc_data_table[, subtype]))
  plot_color2 = plot_color[subtype_reduced_list]
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table, aes(x = dir, y = as.numeric(logFC), color=subtype))+
    #ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observation
    #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
    scale_color_manual(name="", values = plot_color2)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    geom_violin(trim=FALSE, color="black", fill=NA)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())+
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  pval1 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig1 <- sig_function(pval1)
  pval2 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig2 <- sig_function(pval2)
  pval3 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")), logFC])$p.value
  sig3 <- sig_function(pval3)
  sigs <- data.table(rbind(cbind(compar=paste("Non-MR, non-MA Non-DEG;", paste(subclass_abbr, "MR Low")), sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar=paste(paste(subclass_abbr, "MR Low;"), paste(subclass_abbr, "MR High")), sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", paste(subclass_abbr, "MR High")), sig_symbol=sig3, wilcox.pval=pval3)))
  return(sigs)
  return(vplot)
}

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=c("Pvalb"),
                                                                 subclass_abbr="PV",
                                                                 low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                 ymin=-1.5, ymax=1.5,
                                                                 plot_width=3, plot_height=5,
                                                                 plot_color=colors_subtype,
                                                                 plot_save=TRUE,
                                                                 plot_title="",
                                                                 plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors",
                                                                 return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=c("Sst"),
                                                                 subclass_abbr="SST",
                                                                 low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                 ymin=-1.5, ymax=1.5,
                                                                 plot_width=3, plot_height=5,
                                                                 plot_color=colors_subtype,
                                                                 plot_save=TRUE,
                                                                 plot_title="",
                                                                 plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors",
                                                                 return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=Rbp4_L5_targets,
                                                                  subclass_abbr="L5",
                                                                  custom_subtype_list=TRUE,
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.5, ymax=1.5,
                                                                  plot_width=3, plot_height=5,
                                                                  plot_color=colors_subtype,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Rbp4.L5.targets.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors",
                                                                  return_table=FALSE)

##smaller
subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=c("Pvalb"),
                                                                 subclass_abbr="PV",
                                                                 low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                 ymin=-1.5, ymax=1.5,
                                                                 plot_width=2.5, plot_height=5,
                                                                 plot_color=colors_subtype,
                                                                 plot_save=TRUE,
                                                                 plot_title="",
                                                                 plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                 return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=c("Sst"),
                                                                 subclass_abbr="SST",
                                                                 low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                 ymin=-1.5, ymax=1.5,
                                                                 plot_width=2.5, plot_height=5,
                                                                 plot_color=colors_subtype,
                                                                 plot_save=TRUE,
                                                                 plot_title="",
                                                                 plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                 return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_INTACTMR_subtypeColor_func(subclass_list=Rbp4_L5_targets,
                                                                 subclass_abbr="L5",
                                                                 custom_subtype_list=TRUE,
                                                                 low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                 ymin=-1.5, ymax=1.5,
                                                                 plot_width=2.5, plot_height=5,
                                                                 plot_color=colors_subtype,
                                                                 plot_save=TRUE,
                                                                 plot_title="",
                                                                 plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/Rbp4.L5.targets.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.INTACT.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                 return_table=FALSE)

#restricting jitter inside violin plot
#install.packages("ggforce")
library(ggforce)
ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR[(as.numeric(logFC)>=as.numeric(outlier_min)) & (as.numeric(logFC)<=as.numeric(outlier_max))], aes(x = neighborhood, y = as.numeric(logFC), color=subtype))+
  ggtitle("")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  #separate by gene class and direction of change
  facet_grid(.~neighborhood + gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  
  #geom_jitter(size=0.6, alpha=0.9) +
  geom_sina(size=0.6, alpha=0.9)+
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  #trim=TRUE means the data is trimmed to fit the range of observations
  geom_violin(trim=TRUE, color="black", fill=NA)+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[, subtype], cluster_label]))]) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/geom.sina.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.png", width = 15.5, height = 6, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.png", width = 15.5, height = 6, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.eps", width = 15.5, height = 6, dpi = 300, units = "in", device=cairo_ps, fallback_resolution = 300)



###looking for examples
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table=fread("HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.txt")


yao2021_DE[[c("114_116", "down.genes")]]
yao2021_DE[[c("114_116", "up.genes")]]

yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(subclass_label.x=="Pvalb") & (cluster_label.x=="114_Pvalb")]

head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(subtype=="114_Pvalb") & (dir=="Low") & (p.adj.BH.all < 0.1)])

head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(subtype=="116_Pvalb") & (dir=="Low") & (p.adj.BH.all < 0.1)])


c("Cnr1   Reln Lypd6b  Lypd6") 

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
WT.indices <- colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = t.type == "WT"))
KO.indices <- colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = t.type == "KO"))

#data table of cells and their types
cellTypes_dt <- data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[c("predicted.id", "subclass_label", "neighborhood_label", "t.type", "rep")]], keep.rownames="index")

#combining metadata of all mecp2 Het experiments
cell_metadata_gene_counts_mecp2Het <- data.table(rbind(cell_metadata_gene_counts, cell_metadata_gene_counts2, cell_metadata_gene_counts6, cell_metadata_gene_counts7))

cellTypes_meta <- left_join(x=cellTypes_dt, y=cell_metadata_gene_counts_mecp2Het, by=c("index"="X"))

#Euclidean distance function
euc <- function(a, b) sqrt(sum((a - b)^2))

euc(cellTypes_meta[1, .(center_x, center_y)], cellTypes_meta[2, .(center_x, center_y)])

#finding close cells; try 100 microns away or less
rep1_indices <- cellTypes_meta[rep=="Rep1", index]
rep2_indices <- cellTypes_meta[rep=="Rep2", index]

rep6_indices <- cellTypes_meta[rep=="Rep6", index]


All_Vis_layerdepth <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/layer_depth/All_Vis_layerdepth_fixed.csv")
#add class label
All_Vis_layerdepth <- left_join(x=All_Vis_layerdepth, y=CTX_HIP_annot[, .(cluster_label, subclass_label, neighborhood_label, class_label)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))

All_Vis_layerdepth_Pvalb <- All_Vis_layerdepth[subclass_label=="Pvalb"]

rep2_Vis_Pvalb_indices <- cellTypes_meta[(rep=="Rep2") & (index %in% All_Vis_layerdepth_Pvalb$cell.id), index]
rep2_Vis_Pvalb_pairs <- t(combn(rep2_Vis_Pvalb_indices,2))


euc(cellTypes_meta[index==rep2_Vis_Pvalb_pairs[1,1], .(center_x, center_y)], cellTypes_meta[index==rep2_Vis_Pvalb_pairs[1,2], .(center_x, center_y)])

cellTypes.dt <- data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[c("rep", "t.type", "predicted.id", "subclass_label", "neighborhood_label")]], keep.rownames="index")

rep2_Vis_Pvalb_dist <- data.table(index1=character(), predicted.id1=character(), t.type1=character(), index2=character(), predicted.id2=character(), t.type2=character(), dist=numeric())

for(i in 1:nrow(rep2_Vis_Pvalb_pairs)){
  index1 <- rep2_Vis_Pvalb_pairs[i,1]
  predicted.id1 <- cellTypes.dt[index==index1, predicted.id]
  t.type1 <- cellTypes.dt[index==index1, t.type]
  
  index2 <- rep2_Vis_Pvalb_pairs[i,2]
  predicted.id2 <- cellTypes.dt[index==index2, predicted.id]
  t.type2 <- cellTypes.dt[index==index2, t.type]
  
  dist <- euc(cellTypes_meta[index==index1, .(center_x, center_y)], cellTypes_meta[index==index2, .(center_x, center_y)])
  
  line <- cbind(index1=index1, predicted.id1=predicted.id1, t.type1=t.type1, 
                index2=index2, predicted.id2=predicted.id2, t.type2=t.type2,
                dist=dist)
  rep2_Vis_Pvalb_dist  <- rbind(rep2_Vis_Pvalb_dist, line)
}

head(rep2_Vis_Pvalb_dist)

rep2_Vis_Pvalb_dist$dist <- as.numeric(rep2_Vis_Pvalb_dist$dist)

rep2_Vis_Pvalb_dist <- left_join(x=rep2_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, center_x, center_y)], by=c("index1"="index"))
rep2_Vis_Pvalb_dist <- left_join(x=rep2_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, center_x, center_y)], by=c("index2"="index"))
names(rep2_Vis_Pvalb_dist)[8:11] = c("center_x1", "center_y1", "center_x2", "center_y2")

rep2_Vis_Pvalb_dist <- left_join(x=rep2_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, Mecp2)], by=c("index1"="index"))
rep2_Vis_Pvalb_dist <- left_join(x=rep2_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, Mecp2)], by=c("index2"="index"))
names(rep2_Vis_Pvalb_dist)[12:13] = c("Mecp2.1", "Mecp2.2")

join_col_names = c("index", names(cellTypes_meta)[15:564])
head(cellTypes_meta[, ..join_col_names])

rep2_Vis_Pvalb_dist2 <- left_join(x=rep2_Vis_Pvalb_dist, y=cellTypes_meta[, ..join_col_names], by=c("index1"="index"))
rep2_Vis_Pvalb_dist2 <- left_join(x=rep2_Vis_Pvalb_dist2, y=cellTypes_meta[, ..join_col_names], by=c("index2"="index"))



#pairs of cells with different transcriptotypes
rep2_Vis_Pvalb_dist_filt <- rep2_Vis_Pvalb_dist2[t.type1 != t.type2]
#same subtype
rep2_Vis_Pvalb_dist_filt <- rep2_Vis_Pvalb_dist_filt[predicted.id1 == predicted.id2]

View(rep2_Vis_Pvalb_dist_filt[(Mecp2.1 >5) | (Mecp2.2 >5)])

rep2_Vis_Pvalb_dist_filt[(Mecp2.1 >5) | (Mecp2.2 >5)]
#
#rep2_Vis_Pvalb_dist_filt_close <- 

head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(subclass_label=="Pvalb") & (dir=="Low") & (logFC > 0) & (p.adj.BH.all <0.1)])

head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Any recurrent MR") & (subclass_label=="Pvalb") & (dir=="Low") & (logFC > 0) & (p.adj.BH.all <0.1)])
head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Any recurrent MR") & (subclass_label=="Pvalb") & (dir=="Low") & (logFC > 1.3)])
head(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[(gene_class=="Any recurrent MR") & (subclass_label=="Pvalb") & (dir=="Low") & (logFC > 0) & (FDR < 0.1)])

View(rep2_Vis_Pvalb_dist_filt[(Mecp2.1 >3) & (Sdk1.x < 3) & (Mecp2.2 <3) & (Sdk1.y > 3)])

View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="116_Pvalb") & (predicted.id2=="116_Pvalb") & (t.type1=="KO") & (t.type2=="WT") & (Sdk1.x > 3) & (Sdk1.y < 2)])
View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="116_Pvalb") & (predicted.id2=="116_Pvalb") & (t.type1=="KO") & (t.type2=="WT") & (dist<200) & (Dcc.x==3) & (Dcc.y==0)])
View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="113_Pvalb") & (predicted.id2=="113_Pvalb") & (t.type1=="KO") & (t.type2=="WT") & (dist<200) & (Fat1.x>=5) & (Fat1.y<=1)])
View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="113_Pvalb") & (predicted.id2=="113_Pvalb") & (t.type1=="WT") & (t.type2=="KO") & (dist<200) & (Fat1.y>=5) & (Fat1.x<=1)])

View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="119_Pvalb") & (predicted.id2=="119_Pvalb") & (t.type1=="WT") & (t.type2=="KO") & (dist<200) & (Fat1.y>=5) & (Fat1.x<=1)])
View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="119_Pvalb") & (predicted.id2=="119_Pvalb") & (t.type1=="WT") & (t.type2=="KO") & (dist<200) & (Cdh6.y>=3) & (Cdh6.x<=1)])

View(rep2_Vis_Pvalb_dist_filt[(predicted.id1=="116_Pvalb") & (predicted.id2=="116_Pvalb") & (t.type1=="KO") & (t.type2=="WT") & (dist<200) & (Col25a1.x>=1) & (Col25a1.y<=1)])
Col25a1


hist(cellTypes_meta[rep=="Rep2", Dcc])
hist(cellTypes_meta[rep=="Rep2", Lypd6])
hist(cellTypes_meta[rep=="Rep2", Flrt3])
hist(cellTypes_meta[rep=="Rep2", Fat1])
hist(cellTypes_meta[rep=="Rep2", Cdh6])

###
rep6_Vis_Pvalb_indices <- cellTypes_meta[(rep=="Rep6") & (index %in% All_Vis_layerdepth_Pvalb$cell.id), index]
rep6_Vis_Pvalb_pairs <- t(combn(rep6_Vis_Pvalb_indices,2))

rep6_Vis_Pvalb_dist <- data.table(index1=character(), predicted.id1=character(), t.type1=character(), index2=character(), predicted.id2=character(), t.type2=character(), dist=numeric())

for(i in 1:nrow(rep6_Vis_Pvalb_pairs)){
  index1 <- rep6_Vis_Pvalb_pairs[i,1]
  predicted.id1 <- cellTypes.dt[index==index1, predicted.id]
  t.type1 <- cellTypes.dt[index==index1, t.type]
  
  index2 <- rep6_Vis_Pvalb_pairs[i,2]
  predicted.id2 <- cellTypes.dt[index==index2, predicted.id]
  t.type2 <- cellTypes.dt[index==index2, t.type]
  
  dist <- euc(cellTypes_meta[index==index1, .(center_x, center_y)], cellTypes_meta[index==index2, .(center_x, center_y)])
  
  line <- cbind(index1=index1, predicted.id1=predicted.id1, t.type1=t.type1, 
                index2=index2, predicted.id2=predicted.id2, t.type2=t.type2,
                dist=dist)
  rep6_Vis_Pvalb_dist  <- rbind(rep6_Vis_Pvalb_dist, line)
}

head(rep6_Vis_Pvalb_dist)

rep6_Vis_Pvalb_dist$dist <- as.numeric(rep6_Vis_Pvalb_dist$dist)

rep6_Vis_Pvalb_dist <- left_join(x=rep6_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, center_x, center_y)], by=c("index1"="index"))
rep6_Vis_Pvalb_dist <- left_join(x=rep6_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, center_x, center_y)], by=c("index2"="index"))
names(rep6_Vis_Pvalb_dist)[8:11] = c("center_x1", "center_y1", "center_x2", "center_y2")

rep6_Vis_Pvalb_dist <- left_join(x=rep6_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, Mecp2)], by=c("index1"="index"))
rep6_Vis_Pvalb_dist <- left_join(x=rep6_Vis_Pvalb_dist, y=cellTypes_meta[, .(index, Mecp2)], by=c("index2"="index"))
names(rep6_Vis_Pvalb_dist)[12:13] = c("Mecp2.1", "Mecp2.2")

rep6_Vis_Pvalb_dist2 <- left_join(x=rep6_Vis_Pvalb_dist, y=cellTypes_meta[, ..join_col_names], by=c("index1"="index"))
rep6_Vis_Pvalb_dist2 <- left_join(x=rep6_Vis_Pvalb_dist2, y=cellTypes_meta[, ..join_col_names], by=c("index2"="index"))

####

all.nhoods.func4 <- function(gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table){
  scaled.DEG.table.L2.3.IT  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "L2/3 IT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  L2.3.IT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.L2.3.IT)
  
  scaled.DEG.table.L4.5.6.IT.Car3  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "L4/5/6 IT Car3", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  L4.5.6.IT.Car3.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.L4.5.6.IT.Car3)
  
  scaled.DEG.table.MGE  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "MGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  MGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.MGE)
  
  scaled.DEG.table.CGE  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "CGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  CGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.CGE)
  
  scaled.DEG.table.PT  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "PT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  PT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.PT)
  
  scaled.DEG.table.NP.CT.L6b  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "NP/CT/L6b", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  NP.CT.L6b.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.NP.CT.L6b)
  
  scaled.DEG.table.DG.SUB.CA  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "DG/SUB/CA", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  DG.SUB.CA.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.DG.SUB.CA)
  
  scaled.DEG.table.Other  = group.scaled.DEG.func2(grouping="neighborhood", nhood_label = "Other", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table, avg.exp.table=avg.exp.table)
  Other.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR_table(scaled.DEG.table = scaled.DEG.table.Other)
  
  all.nhoods.dt = data.table(rbind(
    cbind(L2.3.IT.dt, neighborhood="L2/3 IT"),
    cbind(L4.5.6.IT.Car3.dt, neighborhood="L4/5/6 IT Car3"),
    cbind(MGE.dt, neighborhood="MGE"),
    cbind(CGE.dt, neighborhood="CGE"),
    cbind(PT.dt, neighborhood="PT"),
    cbind(NP.CT.L6b.dt, neighborhood="NP/CT/L6b"),
    cbind(DG.SUB.CA.dt, neighborhood="DG/SUB/CA"),
    cbind(Other.dt, neighborhood="Other")))
  
  all.nhoods.dt = all.nhoods.dt %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE","MGE","L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA", "Other")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "PV MR", "SST MR", "L4 MR", "L5 MR", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(all.nhoods.dt)
}

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2 = all.nhoods.func4(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                                                                                                     de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, 
                                                                                                     binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                                                                                                     avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)
write.table(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2, file="HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.withOther.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.txt", quote=F, row.names=F, sep="\t")
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2=fread("HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.withOther.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.txt")

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 <- rbind(
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG")]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2[(gene_class=="Any recurrent MR") & (dir=="Low")]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table2[(gene_class=="Any recurrent MR") & (dir=="High")])
)

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2 %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE","MGE","L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA", "Other")))



all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR2 <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2, neighborhood, dir, gene_class) %>%
  summarise(
    IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
    outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
    outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[2] +1.5*IQR
  ) %>% data.table

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR2 <-left_join(x=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2,y=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR2, by=c("neighborhood", "dir", "gene_class"))

#only using the subtypes that appear in every comparison of a neighborhood
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 <- all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR2[(as.numeric(logFC) >= outlier_min) & (as.numeric(logFC) <= outlier_max)]

subtype_list_all_nhoods = unique(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[, subtype])
subtype_list_all_nhoods_reduced = c()
for(i in subtype_list_all_nhoods){
  if(nrow(unique(
    all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(subtype==i), .(gene_class, dir)]))==3){
    subtype_list_all_nhoods_reduced = c(subtype_list_all_nhoods_reduced,i)
  }
}


all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3  = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[subtype %in% subtype_list_all_nhoods_reduced]
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3 %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[, subtype], cluster_label]))))

ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3, aes(x = neighborhood, y = as.numeric(logFC), color=subtype))+
  ggtitle("")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  #separate by gene class and direction of change
  facet_grid(.~neighborhood + gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  #trim=TRUE means the data is trimmed to fit the range of observations
  geom_violin(trim=TRUE, color="black", fill=NA)+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table[, subtype], cluster_label]))]) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nhoodsWithOther.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.png", width = 16, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nhoodsWithOther.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.eps", width = 16, height = 6, dpi = 300, units = "in", device=cairo_ps, fallback_resolution = 300)

###trying to make a heatmap
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3, neighborhood) %>%
  summarise(
    low = mean(as.numeric(logFC[dir == 'Low']), na.rm=TRUE) -  mean(as.numeric(logFC[dir == 'Non-DEG']), na.rm=TRUE),
    high = mean(as.numeric(logFC[dir == 'High']), na.rm=TRUE) -  mean(as.numeric(logFC[dir == 'Non-DEG']), na.rm=TRUE)) %>% data.frame(row.names="neighborhood")


all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix = t(as.matrix(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-0.16,-0.011,length=100),  # for blue
               seq(-0.01, 0.01,length=100),  # for white
               seq(0.011,0.430,length=100)) # for red


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/mecp2.CTX.HC.100vol.300counts.pred0.2.summary.heatmap.low.vs.nonDEG.high.vs.nonDEG.logFoldDiff.BlueToRed.eps")
heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, symkey=FALSE) 
dev.off()

col_breaks2 = c(seq(-0.153,0.185,length=100),  # for blue
               seq(0.186, 0.2,length=100),  # for white
               seq(0.201,0.430,length=100)) # for red

heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, symkey=FALSE) 
#symmetrical colors

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FC.counts <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3, neighborhood) %>%
  summarise(
    n.nonDEG = sum((dir=="Non-DEG")),
    n.low = sum((dir=="Low")),
    n.high = sum((dir=="High"))) %>% data.frame(row.names="neighborhood") %>% as.matrix

wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="PT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="PT"), as.numeric(logFC)])$p.value


nhoods = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix2 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix[, nhoods]


col_breaks3 = c(seq(-0.43,-0.11,length=100),  # for blue
                seq(-0.1, 0.1,length=100),  # for white
                seq(0.11,0.430,length=100))

heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix2,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=col_breaks3,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, symkey=TRUE)

my_palette2 <- colorRampPalette(c("white", "red"))(n = 299)
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/mecp2.CTX.HC.100vol.300counts.pred0.2.summary.heatmap.low.vs.nonDEG.high.vs.nonDEG.min10Genes.eps")
heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix2,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette2,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE)
dev.off()

all.nhoods.anyMetaMR.pval.matrix=c(
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Non-DEG") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2)
)

all.nhoods.anyMetaMR.pval.matrix = matrix(all.nhoods.anyMetaMR.pval.matrix, nrow=2, ncol=6)
rownames(all.nhoods.anyMetaMR.pval.matrix) = c("Low", "High")
colnames(all.nhoods.anyMetaMR.pval.matrix) = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")


all.nhoods.anyMetaMR.pval.matrix.LowVsHigh =c(
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="Low") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3[(dir=="High") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2)
)
all.nhoods.anyMetaMR.pval.matrix.LowVsHigh = matrix(all.nhoods.anyMetaMR.pval.matrix.LowVsHigh, nrow=1, ncol=6)
rownames(all.nhoods.anyMetaMR.pval.matrix.LowVsHigh) = "Low_vs_High"
colnames(all.nhoods.anyMetaMR.pval.matrix.LowVsHigh) = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")


###saving average expression somewhere
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, file="HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv", quote=F, row.names=T)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt, file="HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv", quote=F, row.names=F)

#average expressions in 100vol, 300 count filter KO
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, subset = t.type == "KO")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg.melt) = c("subtype", "gene", "avgExp")

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg, file="HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv", quote=F, row.names=T)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.avg.melt, file="HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv", quote=F, row.names=F)


####using correct outlier removal 
##color dots in violin plot by subtype, any meta MR genes
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2 <- function(subclass_list, custom_subtype_list=FALSE, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_width=3, plot_height=5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  
  subtype_logfc_data_IQR <- group_by(subtype_logfc_data_table, dir, gene_class) %>%
    summarise(
      IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
      outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
      outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[4] +1.5*IQR
    ) %>% data.table
  
  subtype_logfc_data_table <- left_join(x=subtype_logfc_data_table,y=subtype_logfc_data_IQR, by=c("dir", "gene_class"))
  subtype_logfc_data_table <- subtype_logfc_data_table[(as.numeric(logFC) >= outlier_min) & (as.numeric(logFC) <= outlier_max)]
  
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #print(unique(subtype_logfc_data_table[, subtype]))
  plot_color2 = plot_color[subtype_reduced_list]
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table, aes(x = dir, y = as.numeric(logFC), color=subtype))+
    #ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observation
    #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
    scale_color_manual(name="", values = plot_color2)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    geom_violin(trim=TRUE, color="black", fill=NA)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())+
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  
  pval1 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig1 <- sig_function(pval1)
  pval2 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig2 <- sig_function(pval2)
  pval3 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig3 <- sig_function(pval3)
  sigs <- data.table(rbind(cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR Low"), sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar=paste("Any recurrent MR Low;", "Any recurrent MR High"), sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR High"), sig_symbol=sig3, wilcox.pval=pval3)))
  return(sigs)
  return(vplot)
}
subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=c("Pvalb"),
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/correctOutliers.Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=c("Sst"),
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/correctOutliers.Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)

subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=Rbp4_L5_targets,
                                                                  custom_subtype_list=TRUE,
                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                  ymin=-1.2, ymax=2,
                                                                  plot_color=colors_subtype,
                                                                  plot_width=2.5, plot_height=5,
                                                                  plot_save=TRUE,
                                                                  plot_title="",
                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/correctOutliers.Rbp4.L5.targets.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.viol.dotPlot.subtypeColors.smaller",
                                                                  return_table=FALSE)

####
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR.correct <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2, neighborhood, dir, gene_class) %>%
  summarise(
    IQR = quantile(as.numeric(logFC), na.rm=TRUE)[4] - quantile(as.numeric(logFC), na.rm=TRUE)[2],
    outlier_min = quantile(as.numeric(logFC), na.rm=TRUE)[2] - 1.5*IQR,
    outlier_max = quantile(as.numeric(logFC), na.rm=TRUE)[4] +1.5*IQR
  ) %>% data.table

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct2 <-left_join(x=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table2,y=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.IQR.correct, by=c("neighborhood", "dir", "gene_class"))

#only using the subtypes that appear in every comparison of a neighborhood
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 <- all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct2[(as.numeric(logFC) >= outlier_min) & (as.numeric(logFC) <= outlier_max)]

subtype_list_all_nhoods_correct = unique(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[, subtype])
subtype_list_all_nhoods_correct_reduced = c()
for(i in subtype_list_all_nhoods_correct){
  if(nrow(unique(
    all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(subtype==i), .(gene_class, dir)]))==3){
    subtype_list_all_nhoods_correct_reduced = c(subtype_list_all_nhoods_correct_reduced,i)
  }
}


all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3  = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[subtype %in% subtype_list_all_nhoods_correct_reduced]
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3 %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[, subtype], cluster_label]))))

ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3, aes(x = neighborhood, y = as.numeric(logFC), color=subtype))+
  ggtitle("")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  #separate by gene class and direction of change
  facet_grid(.~neighborhood + gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.4, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
  #trim=TRUE means the data is trimmed to fit the range of observations
  geom_violin(trim=TRUE, color="black", fill=NA)+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[, subtype], cluster_label]))]) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nhoodsWithOther.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.png", width = 16, height = 6, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nhoodsWithOther.nonDEG.lowHigh.anyMetaMR.logFC.outliersRemoved.viol.dotPlot.eps", width = 16, height = 6, dpi = 300, units = "in", device=cairo_ps, fallback_resolution = 300)

###trying to make a heatmap
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3, neighborhood) %>%
  summarise(
    low = mean(as.numeric(logFC[dir == 'Low']), na.rm=TRUE) -  mean(as.numeric(logFC[dir == 'Non-DEG']), na.rm=TRUE),
    high = mean(as.numeric(logFC[dir == 'High']), na.rm=TRUE) -  mean(as.numeric(logFC[dir == 'Non-DEG']), na.rm=TRUE)) %>% data.frame(row.names="neighborhood")


all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix = t(as.matrix(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-0.16,-0.011,length=100),  # for blue
               seq(-0.01, 0.01,length=100),  # for white
               seq(0.011,0.430,length=100)) # for red


setEPS()
#postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/correctOutliers.mecp2.CTX.HC.100vol.300counts.pred0.2.summary.heatmap.low.vs.nonDEG.high.vs.nonDEG.logFoldDiff.BlueToRed.eps")
heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=col_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, symkey=FALSE) 
dev.off()

col_breaks2 = c(seq(-0.153,0.185,length=100),  # for blue
                seq(0.186, 0.2,length=100),  # for white
                seq(0.201,0.430,length=100)) # for red

heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, symkey=FALSE) 
#symmetrical colors

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FC.counts <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3, neighborhood) %>%
  summarise(
    n.nonDEG = sum((dir=="Non-DEG")),
    n.low = sum((dir=="Low")),
    n.high = sum((dir=="High"))) %>% data.frame(row.names="neighborhood") %>% as.matrix

wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="PT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="PT"), as.numeric(logFC)])$p.value


nhoods = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix2 = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix[, nhoods]


#col_breaks3 = c(seq(-0.43,-0.11,length=100),  # for blue
 #               seq(-0.1, 0.1,length=100),  # for white
  #              seq(0.11,0.430,length=100))

#heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR3.FCsumm.matrix2,
#          main = "", # heat map title
#          notecol="black",      # change font color of cell labels to black
 #         breaks=col_breaks3,    # enable color transition at specified limits
  #        density.info="none",  # turns off density plot inside color legend
   #       trace="none",         # turns off trace lines inside the heat map
    #      margins =c(12,9),     # widens margins around plot
     #     col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
      #    dendrogram="none",     # only draw a row dendrogram
       #   Rowv="NA",
        #  Colv="NA",key=TRUE, symkey=TRUE)

my_palette2 <- colorRampPalette(c("white", "red"))(n = 299)
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/correctOutliers.mecp2.CTX.HC.100vol.300counts.pred0.2.summary.heatmap.low.vs.nonDEG.high.vs.nonDEG.min10Genes.eps")
heatmap.2(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3.FCsumm.matrix2,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette2,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE)
dev.off()

correct.all.nhoods.anyMetaMR.pval.matrix=c(
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Non-DEG") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2)
)

correct.all.nhoods.anyMetaMR.pval.matrix = matrix(correct.all.nhoods.anyMetaMR.pval.matrix, nrow=2, ncol=6)
rownames(correct.all.nhoods.anyMetaMR.pval.matrix) = c("Low", "High")
colnames(correct.all.nhoods.anyMetaMR.pval.matrix) = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")


correct.all.nhoods.anyMetaMR.pval.matrix.LowVsHigh =c(
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="CGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="CGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="MGE"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="MGE"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="L2/3 IT"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="L2/3 IT"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="L4/5/6 IT Car3"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="NP/CT/L6b"), as.numeric(logFC)])$p.value,2),
  signif(wilcox.test(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="Low") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.anyMetaMR.table.IQR.correct3[(dir=="High") & (neighborhood=="DG/SUB/CA"), as.numeric(logFC)])$p.value,2)
)
correct.all.nhoods.anyMetaMR.pval.matrix.LowVsHigh = matrix(correct.all.nhoods.anyMetaMR.pval.matrix.LowVsHigh, nrow=1, ncol=6)
rownames(correct.all.nhoods.anyMetaMR.pval.matrix.LowVsHigh) = "Low_vs_High"
colnames(correct.all.nhoods.anyMetaMR.pval.matrix.LowVsHigh) = c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "NP/CT/L6b", "DG/SUB/CA")


##saving tables that went into correct outlier violin plots
Pvalb_vioplot_input_table <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=c("Pvalb"),
                                                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                   ymin=-1.2, ymax=2,
                                                                   plot_color=colors_subtype,
                                                                   plot_width=2.5, plot_height=5,
                                                                   plot_save=FALSE,
                                                                   plot_title="",
                                                                   plot_file="",
                                                                   return_table=TRUE)

Sst_vioplot_input_table <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=c("Sst"),
                                                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                   ymin=-1.2, ymax=2,
                                                                   plot_color=colors_subtype,
                                                                   plot_width=2.5, plot_height=5,
                                                                   plot_save=FALSE,
                                                                   plot_title="",
                                                                   plot_file="",
                                                                   return_table=TRUE)

L5_vioplot_input_table <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func2(subclass_list=Rbp4_L5_targets,
                                                                   custom_subtype_list=TRUE,
                                                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                   ymin=-1.2, ymax=2,
                                                                   plot_color=colors_subtype,
                                                                   plot_width=2.5, plot_height=5,
                                                                   plot_save=FALSE,
                                                                   plot_title="",
                                                                   plot_file="",
                                                                   return_table=TRUE)

write.table(Pvalb_vioplot_input_table, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.correctOutliers.Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(Sst_vioplot_input_table, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.correctOutliers.Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(L5_vioplot_input_table, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.correctOutliers.L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")



subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func_noOutlierRemoval <- function(subclass_list, custom_subtype_list=FALSE, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_width=3, plot_height=5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  
  
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==3){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #print(unique(subtype_logfc_data_table[, subtype]))
  plot_color2 = plot_color[subtype_reduced_list]
  #violin plot of log fold changes 
  vplot <- ggplot(subtype_logfc_data_table, aes(x = dir, y = as.numeric(logFC), color=subtype))+
    #ggtitle(plot_title)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    #separate by gene class and direction of change
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) +
    #trim=TRUE means the data is trimmed to fit the range of observation
    #scale_color_manual(values = colors_subtype[unique(subtype_logfc_data_table[!(is.na(gene_class)), subtype])]) +
    scale_color_manual(name="", values = plot_color2)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.4, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", color="black", size=0.4, width=0.1)+
    geom_violin(trim=TRUE, color="black", fill=NA)+
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())+
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
  
  pval1 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig1 <- sig_function(pval1)
  pval2 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig2 <- sig_function(pval2)
  pval3 <- wilcox.test(x=subtype_logfc_data_table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
                       y=subtype_logfc_data_table[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value
  sig3 <- sig_function(pval3)
  sigs <- data.table(rbind(cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR Low"), sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar=paste("Any recurrent MR Low;", "Any recurrent MR High"), sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar=paste("Non-MR, non-MA Non-DEG;", "Any recurrent MR High"), sig_symbol=sig3, wilcox.pval=pval3)))
  return(sigs)
  return(vplot)
}

##saving tables that went into no outlier violin plots
Pvalb_vioplot_input_table_noOutlierRemoval <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func_noOutlierRemoval(subclass_list=c("Pvalb"),
                                                                                                low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                                ymin=-1.2, ymax=2,
                                                                                                plot_color=colors_subtype,
                                                                                                plot_width=2.5, plot_height=5,
                                                                                                plot_save=FALSE,
                                                                                                plot_title="",
                                                                                                plot_file="",
                                                                                                return_table=TRUE)

Sst_vioplot_input_table_noOutlierRemoval <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func_noOutlierRemoval(subclass_list=c("Sst"),
                                                                                              low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                              ymin=-1.2, ymax=2,
                                                                                              plot_color=colors_subtype,
                                                                                              plot_width=2.5, plot_height=5,
                                                                                              plot_save=FALSE,
                                                                                              plot_title="",
                                                                                              plot_file="",
                                                                                              return_table=TRUE)

L5_vioplot_input_table_noOutlierRemoval <- subclass_low_high_nonDEG_viol_dotPlot_anyMetaMR_subtypeColor_func_noOutlierRemoval(subclass_list=Rbp4_L5_targets,
                                                                                             custom_subtype_list=TRUE,
                                                                                             low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table,
                                                                                             ymin=-1.2, ymax=2,
                                                                                             plot_color=colors_subtype,
                                                                                             plot_width=2.5, plot_height=5,
                                                                                             plot_save=FALSE,
                                                                                             plot_title="",
                                                                                             plot_file="",
                                                                                             return_table=TRUE)

write.table(Pvalb_vioplot_input_table_noOutlierRemoval, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(Sst_vioplot_input_table_noOutlierRemoval, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(L5_vioplot_input_table_noOutlierRemoval, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv", quote=FALSE, row.names=FALSE, sep="\t")
