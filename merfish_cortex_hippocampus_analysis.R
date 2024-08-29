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

#cortex and hippocampus ROIs
exp1_L_Ctx_HC = fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/exported regions/exp1_L_Ctx_HC.csv", col.names="index")
exp1_R_Ctx_HC = fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/exported regions/exp1_R_Ctx_HC.csv", col.names="index")
exp1_Ctx_HC = rbind(exp1_L_Ctx_HC, exp1_R_Ctx_HC)

exp2_L_Ctx_HC = fread("HG_lab/Vizgen/experiment2_2022-09-08/exported ROIs from Visualizer/exp2_L_Ctx_HC.csv", col.names="index")
exp2_R_Ctx_HC = fread("HG_lab/Vizgen/experiment2_2022-09-08/exported ROIs from Visualizer/exp2_R_Ctx_HC.csv", col.names="index")
exp2_Ctx_HC = rbind(exp2_L_Ctx_HC, exp2_R_Ctx_HC)

exp6_L_Ctx_HC = fread("HG_lab/Vizgen/experiment6_2023-01-23/exported brain regions/exp6_L_cortex_HC.csv", col.names="index")
exp6_R_Ctx_HC = fread("HG_lab/Vizgen/experiment6_2023-01-23/exported brain regions/exp6_R_cortex_HC.csv", col.names="index")
exp6_Ctx_HC = rbind(exp6_L_Ctx_HC, exp6_R_Ctx_HC)

exp7_L_Ctx_HC = fread("HG_lab/Vizgen/experiment7_2023-01-27/exported regions from visualizer/exp7_L_ctx_HC.csv", col.names="index")
exp7_R_Ctx_HC = fread("HG_lab/Vizgen/experiment7_2023-01-27/exported regions from visualizer/exp7_R_ctx_HC.csv", col.names="index")
exp7_Ctx_HC = rbind(exp7_L_Ctx_HC, exp7_R_Ctx_HC)

#head(cell_metadata_gene_counts[(X %in% exp1_Ctx_HC[, index]), ])

Seurat_obj_index_filter_func <- function(cell_metadata_gene_counts, vol.min, vol.max, count.min, label_subset, gene_names, rep_label){
  exp = cell_metadata_gene_counts[(volume <= vol.max) & (volume >= vol.min) & (total_counts>=count.min) & (X %in% label_subset), ]
  
  count_cols_exp  = c("X", gene_names)
  exp.df = t(data.frame(exp[, ..count_cols_exp], row.names="X"))
  
  exp.obj = CreateSeuratObject(counts = exp.df, project="MERFISH")
  
  ##preparing for analysis through Seurat
  #normalize and scale
  exp.obj <- SCTransform(exp.obj, clip.range = c(-10, 10), )
  exp.obj <- RunPCA(exp.obj, npcs = 30, features = rownames(exp.obj))
  exp.obj <- RunUMAP(exp.obj, dims = 1:30)
  exp.obj <- FindNeighbors(exp.obj, reduction = "pca", dims = 1:30)
  exp.obj <- FindClusters(exp.obj, resolution = 0.3)
  exp.obj$rep <- rep_label
  return(exp.obj)
}

#cortex and hippocampus,  100 vol, 300 count filter
exp1.100vol.300counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts, label_subset=exp1_Ctx_HC[, index], vol.min=100, vol.max=3*median(cell_metadata_gene_counts[, volume]), count.min=300, gene_names=gene_names, rep_label="Rep1")
saveRDS(exp1.100vol.300counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_minVol100_maxVol3Med_min300counts_CTX_HC_Seurat_SCT.rda")
gc()
exp2.100vol.300counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts2, label_subset=exp2_Ctx_HC[, index], vol.min=100, vol.max=3*median(cell_metadata_gene_counts2[, volume]), count.min=300, gene_names=gene_names2, rep_label="Rep2")
saveRDS(exp2.100vol.300counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp2_minVol100_maxVol3Med_min300counts_CTX_HC_Seurat_SCT.rda")
gc()
exp6.100vol.300counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts6, label_subset=exp6_Ctx_HC[, index], vol.min=100, vol.max=3*median(cell_metadata_gene_counts6[, volume]), count.min=300, gene_names=gene_names6, rep_label="Rep6")
saveRDS(exp6.100vol.300counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp6_minVol100_maxVol3Med_min300counts_CTX_HC_Seurat_SCT.rda")
gc()
exp7.100vol.300counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts7, label_subset=exp7_Ctx_HC[, index], vol.min=100, vol.max=3*median(cell_metadata_gene_counts7[, volume]), count.min=300, gene_names=gene_names7, rep_label="Rep7")
saveRDS(exp7.100vol.300counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp7_minVol100_maxVol3Med_min300counts_CTX_HC_Seurat_SCT.rda")
gc()

#cortex and hippocampus,  300 vol, 600 count filter
exp1.300vol.600counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts, label_subset=exp1_Ctx_HC[, index], vol.min=300, vol.max=3*median(cell_metadata_gene_counts[, volume]), count.min=600, gene_names=gene_names, rep_label="Rep1")
saveRDS(exp1.300vol.600counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_minVol300_maxVol3Med_min600counts_CTX_HC_Seurat_SCT.rda")
gc()
exp2.300vol.600counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts2, label_subset=exp2_Ctx_HC[, index], vol.min=300, vol.max=3*median(cell_metadata_gene_counts2[, volume]), count.min=600, gene_names=gene_names2, rep_label="Rep2")
saveRDS(exp2.300vol.600counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp2_minVol300_maxVol3Med_min600counts_CTX_HC_Seurat_SCT.rda")
gc()
exp6.300vol.600counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts6, label_subset=exp6_Ctx_HC[, index], vol.min=300, vol.max=3*median(cell_metadata_gene_counts6[, volume]), count.min=600, gene_names=gene_names6, rep_label="Rep6")
saveRDS(exp6.300vol.600counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp6_minVol300_maxVol3Med_min600counts_CTX_HC_Seurat_SCT.rda")
gc()
exp7.300vol.600counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts7, label_subset=exp7_Ctx_HC[, index], vol.min=300, vol.max=3*median(cell_metadata_gene_counts7[, volume]), count.min=600, gene_names=gene_names7, rep_label="Rep7")
saveRDS(exp7.300vol.600counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp7_minVol300_maxVol3Med_min600counts_CTX_HC_Seurat_SCT.rda")
gc()

#cortex and hippocampus,  500 vol, 900 or 1200 count filter
exp1.500vol.900counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts, label_subset=exp1_Ctx_HC[, index], vol.min=500, vol.max=3*median(cell_metadata_gene_counts[, volume]), count.min=900, gene_names=gene_names, rep_label="Rep1")
saveRDS(exp1.500vol.900counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_minVol500_maxVol3Med_min900counts_CTX_HC_Seurat_SCT.rda")
gc()
exp2.500vol.900counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts2, label_subset=exp2_Ctx_HC[, index], vol.min=500, vol.max=3*median(cell_metadata_gene_counts2[, volume]), count.min=900, gene_names=gene_names2, rep_label="Rep2")
saveRDS(exp2.500vol.900counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp2_minVol500_maxVol3Med_min900counts_CTX_HC_Seurat_SCT.rda")
gc()
exp6.500vol.1200counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts6, label_subset=exp6_Ctx_HC[, index], vol.min=500, vol.max=3*median(cell_metadata_gene_counts6[, volume]), count.min=1200, gene_names=gene_names6, rep_label="Rep6")
saveRDS(exp6.500vol.1200counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp6_minVol500_maxVol3Med_min1200counts_CTX_HC_Seurat_SCT.rda")
gc()
exp7.500vol.900counts.CTX.HC.obj = Seurat_obj_index_filter_func(cell_metadata_gene_counts=cell_metadata_gene_counts7, label_subset=exp7_Ctx_HC[, index], vol.min=500, vol.max=3*median(cell_metadata_gene_counts7[, volume]), count.min=900, gene_names=gene_names7, rep_label="Rep7")
saveRDS(exp7.500vol.900counts.CTX.HC.obj, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp7_minVol500_maxVol3Med_min900counts_CTX_HC_Seurat_SCT.rda")
gc()

#cortex and hippocampus seurat objects
mecp2Het.CTX.HC.100vol.300counts.obj <- readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.300vol.600counts.obj <- readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.500vol.900or1200counts.obj <- readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol500_maxVol3Med_min900or1200counts_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
#randomizing cell order for UMAP plotting of cells (so as not to have one replicate all stacked below another)
mecp2Het.CTX.HC.100vol.300counts.cellOrder = sample(colnames(mecp2Het.CTX.HC.100vol.300counts.obj))
mecp2Het.CTX.HC.300vol.600counts.cellOrder = sample(colnames(mecp2Het.CTX.HC.300vol.600counts.obj))
mecp2Het.CTX.HC.500vol.900or1200counts.cellOrder = sample(colnames(mecp2Het.CTX.HC.500vol.900or1200counts.obj))
#confirming how well replicates integrate through UMAP
DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj, reduction = "umap", group.by = "rep", order = mecp2Het.CTX.HC.100vol.300counts.cellOrder)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300")+
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "bottom")+
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

DimPlot(mecp2Het.CTX.HC.300vol.600counts.obj, reduction = "umap", group.by = "rep", order = mecp2Het.CTX.HC.300vol.600counts.cellOrder)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600")+
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "bottom")+
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.300vol.600counts.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.300vol.600counts.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

DimPlot(mecp2Het.CTX.HC.500vol.900or1200counts.obj, reduction = "umap", group.by = "rep", order = mecp2Het.CTX.HC.500vol.900or1200counts.cellOrder)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n500 <= cell volume <= 3*median(cell volume)\ncounts >= 900 or 1200")+
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "bottom")+
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.500vol.900or1200counts.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.500vol.900or1200counts.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.yao2021.smartseq.predScore.hist.png")
hist(mecp2Het.CTX.HC.100vol.300counts.obj$prediction.score.max, xlab="Seurat max prediction score", main="Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300")
abline(v=median(mecp2Het.CTX.HC.100vol.300counts.obj$prediction.score.max, na.rm=TRUE), col="red", lwd=2, lty=2)
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.300vol.600counts.yao2021.smartseq.predScore.hist.png")
hist(mecp2Het.CTX.HC.300vol.600counts.obj$prediction.score.max, xlab="Seurat max prediction score", main="Mecp2 KO/+ cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600")
abline(v=median(mecp2Het.CTX.HC.300vol.600counts.obj$prediction.score.max, na.rm=TRUE), col="red", lwd=2, lty=2)
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.500vol.900or1200counts.yao2021.smartseq.predScore.hist.png")
hist(mecp2Het.CTX.HC.500vol.900or1200counts.obj$prediction.score.max, xlab="Seurat max prediction score", main="Mecp2 KO/+ cortex and hippocampus\n500 <= cell volume <= 3*median(cell volume)\ncounts >= 900 or 1200")
abline(v=median(mecp2Het.CTX.HC.500vol.900or1200counts.obj$prediction.score.max, na.rm=TRUE), col="red", lwd=2, lty=2)
dev.off()



mecp2Het.CTX.HC.100vol.300counts.sctCounts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj[["SCT"]], slot = "counts")
mecp2Het.CTX.HC.100vol.300counts.sctCounts.mecp2 <- mecp2Het.CTX.HC.100vol.300counts.sctCounts["Mecp2",]

mecp2Het.CTX.HC.300vol.600counts.sctCounts = GetAssayData(mecp2Het.CTX.HC.300vol.600counts.obj[["SCT"]], slot = "counts")
mecp2Het.CTX.HC.300vol.600counts.sctCounts.mecp2 <- mecp2Het.CTX.HC.300vol.600counts.sctCounts["Mecp2",]

mecp2Het.CTX.HC.500vol.900or1200counts.sctCounts = GetAssayData(mecp2Het.CTX.HC.500vol.900or1200counts.obj[["SCT"]], slot = "counts")
mecp2Het.CTX.HC.500vol.900or1200counts.sctCounts.mecp2 <- mecp2Het.CTX.HC.500vol.900or1200counts.sctCounts["Mecp2",]


ttypes_assign_func <- function(seurat_obj, mecp2_counts, KO_filter, WT_filter){
  KO_cells <- names(mecp2_counts[mecp2_counts < KO_filter])
  WT_cells <- names(mecp2_counts[mecp2_counts > WT_filter])
  ttypes <- data.frame(rbind(
    cbind(index=KO_cells , t.type="KO"),
    cbind(index=WT_cells , t.type="WT")), row.names="index")
  seurat_obj_ttypes = AddMetaData(object=seurat_obj, metadata=ttypes, col.name="t.type")
  return(seurat_obj_ttypes)
}

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.100vol.300counts.obj, mecp2_counts = mecp2Het.CTX.HC.100vol.300counts.sctCounts.mecp2, KO_filter=1, WT_filter=1)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.100vol.300counts.obj, mecp2_counts = mecp2Het.CTX.HC.100vol.300counts.sctCounts.mecp2, KO_filter=1, WT_filter=2)

saveRDS(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_Under1Over1_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
saveRDS(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_Under1Over2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.300vol.600counts.obj, mecp2_counts = mecp2Het.CTX.HC.300vol.600counts.sctCounts.mecp2, KO_filter=1, WT_filter=1)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.300vol.600counts.obj, mecp2_counts = mecp2Het.CTX.HC.300vol.600counts.sctCounts.mecp2, KO_filter=1, WT_filter=2)

mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.500vol.900or1200counts.obj, mecp2_counts = mecp2Het.CTX.HC.500vol.900or1200counts.sctCounts.mecp2, KO_filter=1, WT_filter=1)
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2 = ttypes_assign_func(seurat_obj=mecp2Het.CTX.HC.500vol.900or1200counts.obj, mecp2_counts = mecp2Het.CTX.HC.500vol.900or1200counts.sctCounts.mecp2, KO_filter=1, WT_filter=2)

#function for calculated fold changes between KO and WT cells 
WT.KO.de.noshrink.func <- function(seurat_obj){
  sce <- subset(seurat_obj, subset = ((t.type == "WT") | (t.type == "KO")))
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                              colData = sce@meta.data)
  groups.sce <- colData(sce)[, c("predicted.id", "t.type", "rep")]
  aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)
  aggr_counts.filt <- aggr_counts[,aggr_counts$ncells >= 10]
  col.data = colData(aggr_counts.filt)
  col.data$t.type = factor(col.data$t.type, levels=c("WT", "KO"))
  de.results.noShrink <- pseudoBulkDGE(aggr_counts.filt, 
                                       label=aggr_counts.filt$predicted.id,
                                       col.data=col.data,
                                       design=~ rep + t.type,
                                       coef="t.typeKO",
                                       condition=aggr_counts.filt$t.type,
                                       robust=FALSE)
  return(de.results.noShrink)
}
#processing DEG result file to a data table
subtype_de_results_func <- function(de.results){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = "subtype"))
  names(de.results.dt)[2] = "gene"
  return(de.results.dt )
}

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2)
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1)
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2)


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results)
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results)
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results)

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)

#read
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

###DEG results for meta MR and non-mR, non-mA genes in each filtering regime
##100 vol, 300 count
#Under1Over1
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA[, subtype])

#Under1Over2
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA[, subtype])


##300 vol, 600 count
#Under1Over1
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA[, subtype])

#Under1Over2
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA[, subtype])

##500 vol, 900 or 1200 count
#Under1Over1
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA[, subtype])

#Under1Over2
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt[gene %in% unchanged_genes], gene_class="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt[gene %in% gene_panel_classes[(MR_Meta==TRUE), Gene]], gene_class="Meta MR")
))
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA = mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR")))
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.subtypes = unique(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA[, subtype])

#summary
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary <- group_by(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA, subtype, gene_class) %>%
  summarise(
    med.logfc = median(as.numeric(logFC), na.rm=TRUE),
  ) %>% data.table

medlogFC_Under1Over1_compar <- rbind(
 cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary, filter="100 vol, 300 counts"),
 cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary, filter="300 vol, 600 counts"),
 cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt.metaMR.nonMRMA.summary, filter="500 vol, 900 or 1200 counts")
)

medlogFC_Under1Over2_compar <- rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary, filter="100 vol, 300 counts"),
  cbind(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary, filter="300 vol, 600 counts"),
  cbind(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt.metaMR.nonMRMA.summary, filter="500 vol, 900 or 1200 counts")
)
ggplot(medlogFC_Under1Over1_compar, aes(x = gene_class, y = as.numeric(med.logfc), fill=gene_class))+
  ggtitle("Under 1, Over 1")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Median log2 fold changes (KO/WT)\nof subtypes") + xlab("")+
  scale_fill_manual(name = "", values = c("Non-MR, non-MA"="gray", "Meta MR"="red")) +
  facet_grid(.~filter,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/mecp2Het.CTX.HC.filterCompar.yao2021.smartseqLabels.subtypeMedLogFC.metaMR.genes.Under1Over1.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/mecp2Het.CTX.HC.filterCompar.yao2021.smartseqLabels.subtypeMedLogFC.metaMR.genes.Under1Over1.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')


ggplot(medlogFC_Under1Over2_compar, aes(x = gene_class, y = as.numeric(med.logfc), fill=gene_class))+
  ggtitle("Under 1, Over 2")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Median log2 fold changes (KO/WT)\nof subtypes") + xlab("")+
  scale_fill_manual(name = "", values = c("Non-MR, non-MA"="gray", "Meta MR"="red")) +
  facet_grid(.~filter,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/mecp2Het.CTX.HC.filterCompar.yao2021.smartseqLabels.subtypeMedLogFC.metaMR.genes.Under1Over2.boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/mecp2Het.CTX.HC.filterCompar.yao2021.smartseqLabels.subtypeMedLogFC.metaMR.genes.Under1Over2.boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')

#Yao 2021 metadata
yao_meta_dt = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/metadata.csv", na.strings=c("", "NA"))
yao_meta_dt = yao_meta_dt[!is.na(cluster_label), ]

#subclass assignments of cluster labels from Yao 2021
cluster_subclass = unique(yao_meta_dt[, .(cluster_label, subclass_label)])
yao2021_DE = readRDS("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/comb.de.genes.Zizhen.Yao.rda")


CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg.melt) = c("subtype", "gene", "avgExp")


mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes.dt = data.table(mecp2Het.CTX.HC.100vol.300counts.obj[["predicted.id"]], keep.rownames = "index")
mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes.dt = data.table(left_join(x=mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes.dt, y=CTX_HIP_annot[, .(cluster_label, subclass_label, neighborhood_label)], by=c("predicted.id"="cluster_label")))

mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes = AddMetaData(object=mecp2Het.CTX.HC.100vol.300counts.obj, metadata=data.frame(mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes.dt, row.names="index"))

mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes.dt = data.table(mecp2Het.CTX.HC.300vol.600counts.obj[["predicted.id"]], keep.rownames = "index")
mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes.dt = data.table(left_join(x=mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes.dt, y=CTX_HIP_annot[, .(cluster_label, subclass_label, neighborhood_label)], by=c("predicted.id"="cluster_label")))

mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes = AddMetaData(object=mecp2Het.CTX.HC.300vol.600counts.obj, metadata=data.frame(mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes.dt, row.names="index"))

mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes, subset = prediction.score.max > 0.2)
mecp2Het.CTX.HC.300vol.600counts.pred0.2.obj.cellTypes <- subset(mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes, subset = prediction.score.max > 0.2)

saveRDS(mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

yao_binary_dt = data.table(comparison=names(yao2021_DE))
yao_binary_dt = yao_binary_dt[, c("cl.x", "cl.y") := tstrsplit(comparison, "_",fixed=TRUE, type.convert=TRUE)]
yao_binary_dt = left_join(x=yao_binary_dt, y=CTX_HIP_annot[, .(cl, cluster_label, subclass_label)], by=c("cl.x"="cl"))
yao_binary_dt = left_join(x=yao_binary_dt, y=CTX_HIP_annot[, .(cl, cluster_label, subclass_label)], by=c("cl.y"="cl"))

yao_binary_dt$same_subclass = 0
yao_binary_dt[subclass_label.x == subclass_label.y, same_subclass := 1]
write.table(yao_binary_dt, file="HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/yao2021_binary_comparisons_data_table.txt", quote=F, row.names=F, col.names=T, sep="\t")

yao.binary.mecp2Het.CTX.HC.100vol.300counts =  yao_binary_dt[(same_subclass==1) & 
                                        (cluster_label.x %in% mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes$predicted.id) & 
                                        (cluster_label.y %in% mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes$predicted.id) &
                                        (subclass_label.x %in% mecp2Het.CTX.HC.100vol.300counts.obj.cellTypes$subclass_label),]

yao.binary.mecp2Het.CTX.HC.300vol.600counts =  yao_binary_dt[(same_subclass==1) & 
                                                               (cluster_label.x %in% mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes$predicted.id) & 
                                                               (cluster_label.y %in% mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes$predicted.id) &
                                                               (subclass_label.x %in% mecp2Het.CTX.HC.300vol.600counts.obj.cellTypes$subclass_label),]

yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2 =  yao_binary_dt[(same_subclass==1) & 
                                                               (cluster_label.x %in% mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes$predicted.id) & 
                                                               (cluster_label.y %in% mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes$predicted.id) &
                                                               (subclass_label.x %in% mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes$subclass_label),]
write.csv(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Yao2021_pairwise_DEGs/yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.comparisons.csv", quote=F, row.names=F)
yao.binary.mecp2Het.CTX.HC.300vol.600counts.pred0.2 =  yao_binary_dt[(same_subclass==1) & 
                                                               (cluster_label.x %in% mecp2Het.CTX.HC.300vol.600counts.pred0.2.obj.cellTypes$predicted.id) & 
                                                               (cluster_label.y %in% mecp2Het.CTX.HC.300vol.600counts.pred0.2.obj.cellTypes$predicted.id) &
                                                               (subclass_label.x %in% mecp2Het.CTX.HC.300vol.600counts.pred0.2.obj.cellTypes$subclass_label),]


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


CTX_HIP_annot[, cl := as.character(cl)]
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


##function for making scaled (z-score) table of Yao 2021 DEGs for binary comparisons of subtypes in the same subclass
scaled.DEG.func <- function(subclass_label, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table){
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
    
    scaled.DEG.table= left_join(x=cluster.DEG.dt, y=gene.by.subtype.avg.counts.scale, by=c("subtype", "gene")) 
  }
  else{
    scaled.DEG.table = logical(0)
  }
  return(scaled.DEG.table)
}

scaled.DEG.Pvalb <- scaled.DEG.func(subclass_label="Pvalb", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2)
scaled.DEG.Pvalb <- left_join(scaled.DEG.Pvalb, mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt, by=c("subtype", "gene"))
scaled.DEG.Pvalb$avgExp.med <- median(scaled.DEG.Pvalb[, as.numeric(avgExp)], na.rm=TRUE)


##function of plotting gene fold changes aggregated in neighborhoods or a custom subclass list
group.scaled.DEG.func <- function(grouping="neighborhood", nhood_label=NULL, subclass_list=NULL, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table){
  scaled.DEG.table.group = setNames(data.table(matrix(nrow = 0, ncol = 10)), c("cl", "gene", "dir", "subtype", "zScore", "logFC", "logCPM", "F", "PValue", "FDR"))
  if(grouping=="neighborhood"){
    for(i in unique(CTX_HIP_annot[neighborhood_label==nhood_label, subclass_label])){
      scaled.DEG.table.subclass = scaled.DEG.func(subclass_label=i, gene.by.subtype.avg.counts, de.results.table, binary_comparison_table)
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
      scaled.DEG.table.subclass = scaled.DEG.func(i)
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

acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs <- function(scaled.DEG.table, subclass_abbr, avg.exp.table, plot_title, sub_dir, output_file_prefix, logFC_ymin=-2.5, logFC_ymax=2.5, avgExp_ymin=0, avgExp_ymax=70, save_plot=FALSE){
  subtype_zscore_MR_dt = data.table(rbind(
    cbind(scaled.DEG.table[(dir=="Non-DEG") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[get(paste0(subclass_abbr, "_MR"))==TRUE,Gene])], gene_class=paste(subclass_abbr, "MR")),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[get(paste0(subclass_abbr, "_MR"))==TRUE,Gene])], gene_class=paste(subclass_abbr, "MR")),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[get(paste0(subclass_abbr, "_OTHERCELLTYPE"))==TRUE,Gene])], gene_class=paste(subclass_abbr, "other-cell-type MR")),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[get(paste0(subclass_abbr, "_OTHERCELLTYPE"))==TRUE,Gene])], gene_class=paste(subclass_abbr, "other-cell-type MR")),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% short_lowmC_genes)], gene_class="Short low mC"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% short_lowmC_genes)], gene_class="Short low mC")))
  subtype_zscore_MR_dt = left_join(x=subtype_zscore_MR_dt, y=avg.exp.table, by=c("subtype", "gene"))
  subtype_zscore_MR_dt_nonDEG = subtype_zscore_MR_dt[(dir=="Non-DEG")]
  subtype_zscore_MR_dt_up = subtype_zscore_MR_dt[(dir=="High") & (avgExp > quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt_down = subtype_zscore_MR_dt[(dir=="Low") & (avgExp < quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt = rbind(subtype_zscore_MR_dt_nonDEG, subtype_zscore_MR_dt_up, subtype_zscore_MR_dt_down)
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR", paste(subclass_abbr, "MR"), paste(subclass_abbr, "other-cell-type MR"), "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(subtype_zscore_MR_dt[!(is.na(gene_class))])
  ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(logFC_ymin,logFC_ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if (file.exists(sub_dir)){
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 15, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 15, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 15, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 15, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
  
  ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(avgExp), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(avgExp_ymin,avgExp_ymax))+
    ylab("Average expression") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if (file.exists(sub_dir)){
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 15, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 15, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 15, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 15, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
}
#version of above without INTACT MR gene classes
acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs <- function(scaled.DEG.table, avg.exp.table, plot_title, sub_dir, output_file_prefix, logFC_ymin=-2.5, logFC_ymax=2.5, avgExp_ymin=0, avgExp_ymax=70, save_plot=FALSE){
  subtype_zscore_MR_dt = data.table(rbind(
    cbind(scaled.DEG.table[(dir=="Non-DEG") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% short_lowmC_genes)], gene_class="Short low mC"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% short_lowmC_genes)], gene_class="Short low mC")))
  subtype_zscore_MR_dt = left_join(x=subtype_zscore_MR_dt, y=avg.exp.table, by=c("subtype", "gene"))
  subtype_zscore_MR_dt_nonDEG = subtype_zscore_MR_dt[(dir=="Non-DEG")]
  subtype_zscore_MR_dt_up = subtype_zscore_MR_dt[(dir=="High") & (avgExp > quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt_down = subtype_zscore_MR_dt[(dir=="Low") & (avgExp < quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt = rbind(subtype_zscore_MR_dt_nonDEG, subtype_zscore_MR_dt_up, subtype_zscore_MR_dt_down)
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  
  p1 <- ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(logFC_ymin,logFC_ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if(file.exists(sub_dir)){
      p1
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      p1
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
  
  p2 <- ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(avgExp), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(avgExp_ymin,avgExp_ymax))+
    ylab("Average expression") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if(file.exists(sub_dir)){
      p2
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      p2
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
  return(subtype_zscore_MR_dt[!(is.na(gene_class))])
  ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = subtype, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    #scale_color_manual(values = colors_nhood[levels(factor(switched.all.nhoods.dt[, neighborhood]))]) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3)+
    #coord_cartesian(ylim=c(0,18))+
    ylab("LogFC") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(color=guide_legend(nrow=3, byrow=TRUE))
}

mecp2Het.CTX.HC.100vol.300counts.Under1Over1.scaled.DEG.table.Pvalb = scaled.DEG.func(subclass_label="Pvalb", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts)
mecp2Het.CTX.HC.100vol.300counts.Under1Over1.scaled.DEG.table.Sst = scaled.DEG.func(subclass_label="Sst", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts)


mecp2Het.CTX.HC.100vol.300counts.Under1Over1.Pvalb.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = mecp2Het.CTX.HC.100vol.300counts.Under1Over1.scaled.DEG.table.Pvalb,
                                                                               avg.exp.table = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.WT.avg.melt, 
                                                                               plot_title="Pvalb", 
                                                                               sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/",
                                                                               output_file_prefix = "Pvalb.mecp2Het.CTX.HC.100vol.300counts.Under1Over1", save_plot = TRUE)

ggplot(mecp2Het.CTX.HC.100vol.300counts.Under1Over1.Pvalb.dt[!(is.na(gene_class))], aes(x = subtype, y = as.numeric(logFC), color=subtype))+
  ggtitle("Pvalb, Under 1 Over 1")+
  #scale_color_manual(values = colors_nhood[levels(factor(switched.all.nhoods.dt[, neighborhood]))]) +
  stat_summary(fun = mean, 
               geom = "point", size=0.6) + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3)+
  #coord_cartesian(ylim=c(0,18))+
  ylab("LogFC") + xlab("")+
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(color=guide_legend(nrow=3, byrow=TRUE))
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/allNeighbors_exp6_exp7_subclass_yao2021_DEG_MR_avgExp_mean_se_only.png", width = 12, height = 6, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/allNeighbors_exp6_exp7_subclass_yao2021_DEG_MR_avgExp_mean_se_only.eps", width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)

#summary tables 
sctCount_summary_table_func <- function(seurat_obj_full, seurat_obj_Under1Over1, seurat_obj_Under1Over2, replicate, predScoreFilter, cell_metadata_gene_count_table, annotation_file=CTX_HIP_annot, output_file=NULL, write_file=FALSE, return_file=TRUE){
  Under1Over1.labels = data.table(seurat_obj_Under1Over1[[c("rep", "t.type", "predicted.id", "prediction.score.max")]], keep.rownames="index")
  Under1Over2.labels = data.table(seurat_obj_Under1Over2[[c("rep", "t.type", "predicted.id", "prediction.score.max")]], keep.rownames="index")
  names(Under1Over1.labels)[3] = "t.type.Under1Over1"
  names(Under1Over2.labels)[3] = "t.type.Under1Over2"
  
  ttype.labels = left_join(x=Under1Over1.labels, Under1Over2.labels, by=c("index", "rep", "predicted.id", "prediction.score.max"))
  ttype.labels = ttype.labels[, .(index, rep, t.type.Under1Over1, t.type.Under1Over2, prediction.score.max, predicted.id)]

  obj = subset(seurat_obj_full, subset = ((rep == replicate) & (prediction.score.max > predScoreFilter)))
  obj.sctTable = data.table(t(as.matrix(GetAssayData(obj[["SCT"]], slot = "counts"))), keep.rownames = "index")
  obj.sctTable$total_sct_counts = rowSums(t(as.matrix(GetAssayData(obj[["SCT"]], slot = "counts"))))
  obj.sctTable = left_join(x=obj.sctTable, y=cell_metadata_gene_count_table[, .(X, total_counts, fov, volume, center_x, center_y, min_x, max_x, min_y, max_y)], by=c('index'='X'))
  names(obj.sctTable)[553] = "total_raw_counts"
  obj.sctTable = left_join(x=obj.sctTable, y=ttype.labels, by=c("index"))
  obj.sctTable = left_join(x=obj.sctTable, y=annotation_file[, .(cluster_label, subclass_label, neighborhood_label)], by=c("predicted.id"="cluster_label"))
  if(write_file==TRUE){
    write.csv(obj.sctTable, file=output_file, row.names=FALSE, quote=F)
  }
  if(return_file==TRUE){
    return(obj.sctTable)
  }
}

##summary tables of sctransform-corrected counts and other metadata
#100 vol, 300 count
exp1.CTX.HC.100vol.300counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.100vol.300counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2,
                            replicate="Rep1", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)

exp2.CTX.HC.100vol.300counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.100vol.300counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2,
                            replicate="Rep2", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts2,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)

exp6.CTX.HC.100vol.300counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.100vol.300counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2,
                            replicate="Rep6", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts6,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp6.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)

exp7.CTX.HC.100vol.300counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.100vol.300counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2,
                            replicate="Rep7", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts7,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp7.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)


exp3.100vol.300counts.obj=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/exp3_WT_coronal_minVol100_maxVol3Med_min300counts_Seurat_SCT_labelTransfer.rda")
exp3.100vol.300counts.sctCounts = GetAssayData(exp3.100vol.300counts.obj[["SCT"]], slot = "counts")
exp3.100vol.300counts.sctCounts.mecp2 <- exp3.100vol.300counts.sctCounts["Mecp2",]
exp3.100vol.300counts.obj.Under1Over1 = ttypes_assign_func(seurat_obj=exp3.100vol.300counts.obj, mecp2_counts = exp3.100vol.300counts.sctCounts.mecp2, KO_filter=1, WT_filter=1)
exp3.100vol.300counts.obj.Under1Over2 = ttypes_assign_func(seurat_obj=exp3.100vol.300counts.obj, mecp2_counts = exp3.100vol.300counts.sctCounts.mecp2, KO_filter=1, WT_filter=2)
exp3.100vol.300counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=exp3.100vol.300counts.obj, seurat_obj_Under1Over1 = exp3.100vol.300counts.obj.Under1Over1, seurat_obj_Under1Over2 = exp3.100vol.300counts.obj.Under1Over2,
                                                                             replicate="Exp3 WT coronal", predScoreFilter=0.2, cell_metadata_gene_count_table=WT_exp_cell_metadata_gene_counts,
                                                                             annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp3.100vol.300counts.pred0.2.sctCount.celltype.summary.table.csv",
                                                                             return_file=TRUE,
                                                                             write_file=TRUE)

#300 vol, 600 count
exp1.CTX.HC.300vol.600counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.300vol.600counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2,
                            replicate="Rep1", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp1.300vol.600counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)
exp2.CTX.HC.300vol.600counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.300vol.600counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2,
                            replicate="Rep2", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts2,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp2.300vol.600counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)
exp6.CTX.HC.300vol.600counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.300vol.600counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2,
                            replicate="Rep6", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts6,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp6.300vol.600counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)
exp7.CTX.HC.300vol.600counts.pred0.2.sctTable <- sctCount_summary_table_func(seurat_obj_full=mecp2Het.CTX.HC.300vol.600counts.obj, seurat_obj_Under1Over1 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1, seurat_obj_Under1Over2 = mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2,
                            replicate="Rep7", predScoreFilter=0.2, cell_metadata_gene_count_table=cell_metadata_gene_counts7,
                            annotation_file=CTX_HIP_annot, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp7.300vol.600counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv",
                            return_file=TRUE,
                            write_file=FALSE)

####prediction score filter of 0.2
mecp2Het.CTX.HC.100vol.300counts.pred0.2 <- subset(mecp2Het.CTX.HC.100vol.300counts.obj, subset = prediction.score.max > 0.2)
mecp2Het.CTX.HC.300vol.600counts.pred0.2 <- subset(mecp2Het.CTX.HC.300vol.600counts.obj, subset = prediction.score.max > 0.2)

saveRDS(mecp2Het.CTX.HC.100vol.300counts.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
saveRDS(mecp2Het.CTX.HC.300vol.600counts.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_pred0.2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2 <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1, subset = prediction.score.max > 0.2)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2 <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2, subset = prediction.score.max > 0.2)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2 <- subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1, subset = prediction.score.max > 0.2)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2 <- subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2, subset = prediction.score.max > 0.2)

saveRDS(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
saveRDS(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
saveRDS(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_pred0.2_Under1Over1_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
saveRDS(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol300_maxVol3Med_min600counts_pred0.2_Under1Over2_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results <- WT.KO.de.noshrink.func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2)


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results)
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt <- subtype_de_results_func(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results)

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

#average expressions in 100vol, 300 count filter
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

#average expressions in 300vol, 100 count filter
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT = subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT = subset(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2, subset = t.type == "WT")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT

mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg), keep.rownames="subtype"), id.vars="subtype"))
names(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg.melt) = c("subtype", "gene", "avgExp")


subclass_count_hist_func <- function(current_subclass, gene, cell_metadata_gene_counts_table, counts_matrix, min.vol, min.count, title, count_type="Sctransform", xmin=0, xmax=20){
  cell_metadata_gene_counts_filt_subclass <- cell_metadata_gene_counts_table[X %in% all.exp.types.100vol.300counts[(subclass_label==current_subclass) & (prediction.score.max>0.2),index] & (volume <= 3*median(cell_metadata_gene_counts_table[, volume]) & (volume >= min.vol) & (total_counts >= min.count)),]
  gene.counts = counts_matrix[gene,cell_metadata_gene_counts_filt_subclass[, X]]
  hist(gene.counts, breaks=100, main="", xlab=paste(count_type, gene, "counts"), xlim=c(xmin,xmax))
  abline(v=median(gene.counts), col="red", lwd=2, lty=2)
  title(paste0(title), adj = 0.5, line = 0)
}

subclasses4_count_hist_func <- function(subclass1, subclass2, subclass3, subclass4, gene, cell_metadata_gene_counts_table, counts_matrix=exp1.100vol.300counts.sctCounts, min.vol=100, min.count=300, count_type, plot_title, output_file, xmin=xmin, xmax=xmax){
  png(output_file, width=1600, height=3400, res=300)
  par(mfcol=c(4,1))
  subclass_count_hist_func(current_subclass=subclass1, gene=gene, cell_metadata_gene_counts_table=cell_metadata_gene_counts_table, counts_matrix=counts_matrix, min.vol=min.vol, min.count=min.count, title=paste0(plot_title,"\n",subclass1), count_type=count_type, xmin=xmin, xmax=xmax) 
  subclass_count_hist_func(current_subclass=subclass2, gene=gene, cell_metadata_gene_counts_table=cell_metadata_gene_counts_table, counts_matrix=counts_matrix, min.vol=min.vol, min.count=min.count, title=subclass2, count_type=count_type, xmin=xmin, xmax=xmax) 
  subclass_count_hist_func(current_subclass=subclass3, gene=gene, cell_metadata_gene_counts_table=cell_metadata_gene_counts_table, counts_matrix=counts_matrix, min.vol=min.vol, min.count=min.count, title=subclass3, count_type=count_type, xmin=xmin, xmax=xmax) 
  subclass_count_hist_func(current_subclass=subclass4, gene=gene, cell_metadata_gene_counts_table=cell_metadata_gene_counts_table, counts_matrix=counts_matrix, min.vol=min.vol, min.count=min.count, title=subclass4, count_type=count_type, xmin=xmin, xmax=xmax) 
  dev.off()
}

#sctransform-corrected counts
exp1.CTX.HC.100vol.300counts.pred0.2.sctCounts = GetAssayData(exp1.100vol.300counts.obj[["SCT"]], slot = "counts")

hist(exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label=="Astro", Mecp2], main="Exp1, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nAstro", xlab="Sctransform Mecp2 counts")

hist(exp1.CTX.HC.300vol.600counts.pred0.2.sctTable[, Mecp2], main="Exp1, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2", xlab="Sctransform Mecp2 counts", breaks=100)
hist(exp1.CTX.HC.300vol.600counts.pred0.2.sctTable[subclass_label=="Astro", Mecp2], main="Exp1, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2\nAstro", xlab="Sctransform Mecp2 counts", breaks=100)


subclass_sctCount_hist_func <- function(current_subclass, gene, sct_count_table, title, xmin=0, xmax=20){
  gene.counts = sct_count_table[subclass_label==current_subclass, get(gene)]
  hist(gene.counts, breaks=100, main="", xlab=paste("Sctransform", gene, "counts"), xlim=c(xmin,xmax))
  title(paste0(title), adj = 0.5, line = -1)
}

subclass_sctCount_hist_func(current_subclass="Astro", gene="Mecp2", sct_count_table=exp1.CTX.HC.300vol.600counts.pred0.2.sctTable, title="Exp1, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2\nAstro")


subclasses12_sctCount_hist_func <- function(subclass1, subclass2, subclass3, subclass4, subclass5, subclass6, subclass7, subclass8, subclass9, subclass10, subclass11, subclass12, gene, sct_count_table, plot_title, output_file, xmin=xmin, xmax=xmax){
  png(output_file, width=3400, height=3400, res=300)
  par(mfrow=c(4,3))
  subclass_sctCount_hist_func(current_subclass=subclass1, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass1), xmin=xmin, xmax=xmax) 
  subclass_sctCount_hist_func(current_subclass=subclass2, gene=gene, sct_count_table=sct_count_table, title=paste0(plot_title,"\n",subclass2), xmin=xmin, xmax=xmax)   
  subclass_sctCount_hist_func(current_subclass=subclass3, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass3), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass4, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass4), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass5, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass5), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass6, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass6), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass7, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass7), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass8, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass8), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass9, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass9), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass10, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass10), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass11, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass11), xmin=xmin, xmax=xmax)
  subclass_sctCount_hist_func(current_subclass=subclass12, gene=gene, sct_count_table=sct_count_table, title=paste0(subclass12), xmin=xmin, xmax=xmax)
  dev.off()
}

subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp1, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp1.CTX.HC.100vol.300counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp1, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp1.CTX.HC.100vol.300counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")


subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp2.CTX.HC.100vol.300counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp2.CTX.HC.100vol.300counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")


subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp6.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp6, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp6.CTX.HC.100vol.300counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp6.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp6, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp6.CTX.HC.100vol.300counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")

subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp7.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp7, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp7.CTX.HC.100vol.300counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp7.CTX.HC.100vol.300counts.pred0.2.sctTable, plot_title ="Exp7, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp7.CTX.HC.100vol.300counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")


subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp3.100vol.300counts.pred0.2.sctTable, plot_title ="Exp3, WT coronal\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/exp3.WT.coronal.100vol.300counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")

##300vol, 600count
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp1.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp1, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp1.CTX.HC.300vol.600counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp1.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp1, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp1.CTX.HC.300vol.600counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")


subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp2.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp2, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp2.CTX.HC.300vol.600counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp2.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp2, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp2.CTX.HC.300vol.600counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")


subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp6.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp6, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp6.CTX.HC.300vol.600counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp6.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp6, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp6.CTX.HC.300vol.600counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")

subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Mecp2", sct_count_table=exp7.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp7, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=15,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp7.CTX.HC.300vol.600counts.pred0.2.Mecp2.sctransformCounts.hist.subclasses12.png")
subclasses12_sctCount_hist_func(subclass1="Astro", subclass2="DG", subclass3="CA1-ProS", subclass4="Pvalb", subclass5="Sst", subclass6="Vip", subclass7="L2/3 IT CTX", subclass8="L4/5 IT CTX", subclass9="L5 IT CTX", subclass10="L5 PT CTX", subclass11="L6 IT CTX", subclass12="L6 CT CTX",
                                gene="Blank.37", sct_count_table=exp7.CTX.HC.300vol.600counts.pred0.2.sctTable, plot_title ="Exp7, cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2)", xmin=0, xmax=5,
                                output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/histograms/CTX_HC/exp7.CTX.HC.300vol.600counts.pred0.2.Blank.37.sctransformCounts.hist.subclasses12.png")



###
DimPlot(mecp2Het.CTX.HC.100vol.300counts.pred0.2, reduction = "umap", group.by = "predicted.id", order = mecp2Het.CTX.HC.100vol.300counts.cellOrder, raster=FALSE)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2")+
  scale_color_manual(values = colors_subtype) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "None")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.UMAP.by.subtype.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.UMAP.by.subtype.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.UMAP.by.subtype.cairops.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


DimPlot(mecp2Het.CTX.HC.300vol.600counts.pred0.2, reduction = "umap", group.by = "predicted.id", order = mecp2Het.CTX.HC.300vol.600counts.cellOrder)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2")+
  scale_color_manual(values = colors_subtype) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "None")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.300vol.600counts.predScore0.2.UMAP.by.subtype.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.300vol.600counts.predScore0.2.UMAP.by.subtype.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

DimPlot(mecp2Het.CTX.HC.100vol.300counts.pred0.2, reduction = "umap", group.by = "rep", order = mecp2Het.CTX.HC.100vol.300counts.cellOrder)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2")+
  #scale_color_manual(values = colors_subtype) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position = "bottom")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


###
#neighborhood-level data table
#all neighborhoods function
all.nhoods.func <- function(gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table, output_prefix, plot_save=FALSE){
  scaled.DEG.table.L2.3.IT  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "L2/3 IT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  L2.3.IT.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.L2.3.IT,
                                                                                 avg.exp.table = avg.exp.table, 
                                                                                 plot_title="L2/3 IT", 
                                                                                 sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/L2.3.IT/",
                                                                                 output_file_prefix = paste0("L2.3.IT","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.L4.5.6.IT.Car3  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "L4/5/6 IT Car3", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  L4.5.6.IT.Car3.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.L4.5.6.IT.Car3,
                                                                       avg.exp.table = avg.exp.table, 
                                                                       plot_title="L4/5/6 IT Car3", 
                                                                       sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/L4.5.6.IT.Car3/",
                                                                       output_file_prefix = paste0("L4.5.6.IT.Car3","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.MGE  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "MGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  MGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.MGE,
                                                                              avg.exp.table = avg.exp.table, 
                                                                              plot_title="MGE", 
                                                                              sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/MGE/",
                                                                              output_file_prefix = paste0("MGE","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.CGE  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "CGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  CGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.CGE,
                                                                   avg.exp.table = avg.exp.table, 
                                                                   plot_title="CGE", 
                                                                   sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/CGE/",
                                                                   output_file_prefix = paste0("CGE","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.PT  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "PT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  PT.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.PT,
                                                                   avg.exp.table = avg.exp.table, 
                                                                   plot_title="PT", 
                                                                   sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/PT/",
                                                                   output_file_prefix = paste0("PT","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.NP.CT.L6b  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "NP/CT/L6b", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  NP.CT.L6b.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.NP.CT.L6b,
                                                                  avg.exp.table = avg.exp.table, 
                                                                  plot_title="NP/CT/L6b", 
                                                                  sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/NP.CT.L6b/",
                                                                  output_file_prefix = paste0("NP.CT.L6b","_",output_prefix), save_plot = plot_save)
  
  scaled.DEG.table.DG.SUB.CA  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "DG/SUB/CA", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  DG.SUB.CA.dt = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.DG.SUB.CA,
                                                                         avg.exp.table = avg.exp.table, 
                                                                         plot_title="DG/SUB/CA", 
                                                                         sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/DG.SUB.CA/",
                                                                         output_file_prefix = paste0("DG.SUB.CA","_",output_prefix), save_plot = plot_save)
  
  all.nhoods.dt = data.table(rbind(
    cbind(L2.3.IT.dt, neighborhood="L2/3 IT"),
    cbind(L4.5.6.IT.Car3.dt, neighborhood="L4/5/6 IT Car3"),
    cbind(MGE.dt, neighborhood="MGE"),
    cbind(CGE.dt, neighborhood="CGE"),
    cbind(PT.dt, neighborhood="PT"),
    cbind(NP.CT.L6b.dt, neighborhood="NP/CT/L6b"),
    cbind(DG.SUB.CA.dt, neighborhood="DG/SUB/CA")))
  
  all.nhoods.dt = all.nhoods.dt %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE","MGE","L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  all.nhoods.dt = all.nhoods.dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(all.nhoods.dt)
}

function(scaled.DEG.table, avg.exp.table, plot_title, sub_dir, output_file_prefix, logFC_ymin=-2.5, logFC_ymax=2.5, avgExp_ymin=0, avgExp_ymax=70, save_plot=FALSE)
  
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1 = all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, 
                binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt,
                output_prefix="mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1", plot_save=FALSE)

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2 = all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg,
                                                                                  de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt, 
                                                                                  binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                                                                                  avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.WT.avg.melt,
                                                                                  output_prefix="mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2", plot_save=FALSE)

all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1 = all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg,
                                                                                  de.results.table=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt, 
                                                                                  binary_comparison_table=yao.binary.mecp2Het.CTX.HC.300vol.600counts.pred0.2,
                                                                                  avg.exp.table=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.WT.avg.melt,
                                                                                  output_prefix="mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1", plot_save=FALSE)

all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2 = all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg,
                                                                                  de.results.table=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt, 
                                                                                  binary_comparison_table=yao.binary.mecp2Het.CTX.HC.300vol.600counts.pred0.2,
                                                                                  avg.exp.table=mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.WT.avg.melt,
                                                                                  output_prefix="mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2", plot_save=FALSE)


all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt <- data.table(rbind(
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA") & is.finite(as.numeric(logFC))]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(dir=="Low") & (gene_class=="Any recurrent MR") & is.finite(as.numeric(logFC))]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(dir=="High") & (gene_class=="Any recurrent MR") & is.finite(as.numeric(logFC))])
))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA")))


ggplot(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt, aes(x = neighborhood, y = as.numeric(logFC), color=neighborhood))+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1")+
  scale_color_manual(values = colors_nhood[levels(factor(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt[, neighborhood]))]) +
  stat_summary(fun = mean, 
               geom = "point", size=0.6) + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3)+
  coord_cartesian(ylim=c(-0.5,1.2))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only.png", width = 6, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only.eps", width = 6, height = 6, dpi = 300, units = "in", device='eps')

all.nhoods.plot.func <- function(all.nhoods.table, plot_title, ymin=-0.5, ymax=1.2, output_file, plot_save=FALSE){
  all.nhoods.table.nonDEG.and.recurrent.dt <- data.table(rbind(
    cbind(all.nhoods.table[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA") & is.finite(as.numeric(logFC))]),
    cbind(all.nhoods.table[(dir=="Low") & (gene_class=="Any recurrent MR") & is.finite(as.numeric(logFC))]),
    cbind(all.nhoods.table[(dir=="High") & (gene_class=="Any recurrent MR") & is.finite(as.numeric(logFC))])
  ))
  all.nhoods.table.nonDEG.and.recurrent.dt = all.nhoods.table.nonDEG.and.recurrent.dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Any recurrent MR")))
  all.nhoods.table.nonDEG.and.recurrent.dt = all.nhoods.table.nonDEG.and.recurrent.dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  all.nhoods.table.nonDEG.and.recurrent.dt = all.nhoods.table.nonDEG.and.recurrent.dt %>% mutate(neighborhood = factor(neighborhood, levels=c("CGE", "MGE", "L2/3 IT", "L4/5/6 IT Car3", "PT", "NP/CT/L6b", "DG/SUB/CA")))
  
  ggplot(all.nhoods.table.nonDEG.and.recurrent.dt, aes(x = neighborhood, y = as.numeric(logFC), color=neighborhood))+
    ggtitle(plot_title)+
    scale_color_manual(values = colors_nhood[levels(factor(all.nhoods.table.nonDEG.and.recurrent.dt[, neighborhood]))]) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3)+
    coord_cartesian(ylim=c(ymin,ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(output_file ,".png"), width = 6, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(output_file ,".eps"), width = 6, height = 6, dpi = 300, units = "in", device='eps')
  }
  return(all.nhoods.table.nonDEG.and.recurrent.dt)
}

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt <- all.nhoods.plot.func(all.nhoods.table = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1, 
                     plot_title="Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                     ymin=-0.5, ymax=1.2, plot_save=FALSE,
                     output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only")

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.nonDEG.and.recurrent.dt <- all.nhoods.plot.func(all.nhoods.table = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2, 
                     plot_title="Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over2",
                     ymin=-0.5, ymax=1.2, plot_save=FALSE,
                     output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.nonDEG.and.recurrent.MR.logFC.mean.se.only")

all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt <- all.nhoods.plot.func(all.nhoods.table = all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1, 
                     plot_title="Mecp2 KO/+ cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2\nUnder1Over1",
                     ymin=-0.5, ymax=1.2, plot_save=FALSE,
                     output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only")

all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.nonDEG.and.recurrent.dt <- all.nhoods.plot.func(all.nhoods.table = all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2, 
                     plot_title="Mecp2 KO/+ cortex and hippocampus\n300 <= cell volume <= 3*median(cell volume)\ncounts >= 600, predScore > 0.2\nUnder1Over2",
                     ymin=-0.5, ymax=1.2, plot_save=FALSE,
                     output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/all_neighborhood_plots/allNeighbors.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.nonDEG.and.recurrent.MR.logFC.mean.se.only")


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

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt, neighborhood, dir) %>%
  summarise(
    count = n(),
    mean = mean(as.numeric(logFC), na.rm = TRUE),
    sd = sd(as.numeric(logFC), na.rm = TRUE)
  ) %>% data.table

all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary <- group_by(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.dt, neighborhood, dir) %>%
  summarise(
    count = n(),
    mean = mean(as.numeric(logFC), na.rm = TRUE),
    sd = sd(as.numeric(logFC), na.rm = TRUE)
  ) %>% data.table

#100vol, 300 counts
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Non-DEG", "Low")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Low", "High")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Non-DEG", "High")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)

#300vol, 600 counts
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Non-DEG", "Low")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Low", "High")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)
sig_function(t.test(formula = mean ~ dir,
                    data = all.nhoods.mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.nonDEG.and.recurrent.summary[dir %in% c("Non-DEG", "High")],
                    alternative = "two.sided",
                    mu = 0, 
                    paired = TRUE,   
                    var.equal = TRUE,
                    conf.level = 0.95)$p.value)

Pvalb_subtype_list = CTX_HIP_annot[subclass_label %in% c("Pvalb"), cluster_label]
Sst_subtype_list = CTX_HIP_annot[subclass_label %in% c("Sst"), cluster_label]


Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(gene %in% gene_panel_classes[PV_MR==TRUE, Gene]) & (subtype %in% Pvalb_subtype_list) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC)], gene_class="PV MR", dir="All"),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(subtype %in% Pvalb_subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(subtype %in% Pvalb_subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
))

Pvalb_subtype_reduced_list = c()
for(i in Pvalb_subtype_list){
  if(nrow(unique(Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1[subtype==i, .(gene_class, dir)]))==4){
    Pvalb_subtype_reduced_list = c(Pvalb_subtype_reduced_list,i)
  }
}
Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 = Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1[subtype %in% Pvalb_subtype_reduced_list]
Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 = Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "PV MR", "Any recurrent MR")))
Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 = Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "All", "Low", "High")))
Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 = Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1[, subtype], cluster_label]))))

#same as above but with jitter
ggplot(Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1, aes(x = subtype, y = as.numeric(logFC), color=subtype))+
  ggtitle("Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1")+
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.6, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3, color="black")+
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1[, subtype], cluster_label]))]) +
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only.png", width = 6, height = 6, dpi = 300, units = "in", device='png')


Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt[(gene %in% gene_panel_classes[SST_MR==TRUE, Gene]) & (subtype %in% Sst_subtype_list) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC)], gene_class="SST MR", dir="All"),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(subtype %in% Sst_subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
  cbind(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[(subtype %in% Sst_subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
))
Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 = Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "SST MR", "Any recurrent MR")))
Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 = Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "All", "Low", "High")))
Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 = Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1 %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% Sst_nonDEG_allPvMR_anyMetaMR_genes_dt_Under1Over1[, subtype], cluster_label]))))

#same as above but with jitter
ggplot(Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1, aes(x = subtype, y = as.numeric(logFC), color=subtype))+
  ggtitle("Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1")+
  geom_jitter(size=0.6, alpha=0.9) +
  stat_summary(fun = mean, 
               geom = "point", size=0.6, color="black") + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3, color="black")+
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_Under1Over1[, subtype], cluster_label]))]) +
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only.png", width = 6, height = 6, dpi = 300, units = "in", device='png')

subclass_low_high_nonDEG_plot_func <- function(de_results_table, subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(de_results_table[(gene %in% gene_panel_classes[get(paste0(subclass_abbr,"_MR"))==TRUE, Gene]) & (subtype %in% subtype_list) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC)], gene_class=paste(subclass_abbr, "MR"), dir="All"),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==4){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"), "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "All", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #jitter plot of log fold changes 
  ggplot(subtype_logfc_data_table, aes(x = subtype, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.6, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3, color="black")+
    coord_cartesian(ylim=c(-1.5,1.5))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))]) +
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) + 
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = 6, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = 6, height = 6, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
}

Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_plot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                   subclass_list=c("Pvalb"),
                                   subclass_abbr="PV",
                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1,
                                   plot_save=FALSE,
                                   plot_title="Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                   plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only",
                                   return_table=TRUE)

Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_plot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                   subclass_list=c("Sst"),
                                   subclass_abbr="SST",
                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1,
                                   plot_save=FALSE,
                                   plot_title="Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                   plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only",
                                   return_table=TRUE)

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


L5_nonDEG_allL5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_plot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                   subclass_list=Rbp4_L5_targets,
                                   custom_subtype_list=TRUE,
                                   subclass_abbr="L5",
                                   low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1,
                                   plot_save=FALSE,
                                   plot_title="L5\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                   plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Rbp4.L5.targets/L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.and.recurrent.MR.logFC.mean.se.only",
                                   return_table=TRUE)

subtype_nonDEG_anyMeta_sigs <- function(subtype_nonDEG_anyMeta_table){
  subtype_nonDEG_anyMeta_table_summary <- group_by(subtype_nonDEG_anyMeta_table, subtype, dir) %>%
    summarise(
      count = n(),
      mean = mean(as.numeric(logFC), na.rm = TRUE),
      sd = sd(as.numeric(logFC), na.rm = TRUE)
    ) %>% data.table
  pval1 <- t.test(formula = mean ~ dir,
              data = subtype_nonDEG_anyMeta_table_summary[dir %in% c("Non-DEG", "All")],
              alternative = "two.sided",
              mu = 0, 
              paired = TRUE,   
              var.equal = TRUE,
              conf.level = 0.95)$p.value
  sig1 <- sig_function(pval1)
  pval2 <- t.test(formula = mean ~ dir,
                 data = subtype_nonDEG_anyMeta_table_summary[dir %in% c("Low", "High")],
                 alternative = "two.sided",
                 mu = 0, 
                 paired = TRUE,   
                 var.equal = TRUE,
                 conf.level = 0.95)$p.value
  sig2 <- sig_function(pval2)
  pval3 <- t.test(formula = mean ~ dir,
                  data = subtype_nonDEG_anyMeta_table_summary[dir %in% c("Non-DEG", "Low")],
                  alternative = "two.sided",
                  mu = 0, 
                  paired = TRUE,   
                  var.equal = TRUE,
                  conf.level = 0.95)$p.value
  sig3 <- sig_function(pval3)
  pval4 <- t.test(formula = mean ~ dir,
                  data = subtype_nonDEG_anyMeta_table_summary[dir %in% c("Non-DEG", "High")],
                  alternative = "two.sided",
                  mu = 0, 
                  paired = TRUE,   
                  var.equal = TRUE,
                  conf.level = 0.95)$p.value
  sig4 <- sig_function(pval4)
  sigs <- data.table(rbind(cbind(compar="Non-DEG, All", sig=sig1, pval=pval1),
        cbind(compar="Low, High", sig=sig2, pval=pval2),
        cbind(compar="Non-DEG, Low", sig=sig3, pval=pval3),
        cbind(compar="Non-DEG, High", sig=sig4, pval=pval4)))
  return(sigs)
}

subtype_nonDEG_anyMeta_sigs(subtype_nonDEG_anyMeta_table = Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1)
subtype_nonDEG_anyMeta_sigs(subtype_nonDEG_anyMeta_table = Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1)
subtype_nonDEG_anyMeta_sigs(subtype_nonDEG_anyMeta_table = L5_nonDEG_allL5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1)


##looking at scatterplots of the cells in space
library(data.table)
exp1.CTX.HC.100vol.300counts.pred0.2.sctTable <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
plot(exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="118_Pvalb, exp1 cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore>0.2",
     xlab="Exp1, x", ylab="Exp1, y", pch=16, col="gray", cex=0.2)
points(x=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="WT"), center_x], 
       y=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="WT"), center_y],
       pch=16, col="blue", cex=0.4)
points(x=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="KO"), center_x], 
       y=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="KO"), center_y],
       pch=16, col="red", cex=0.4)
###

head(p.adjust(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt$PValue, method="BH"))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt$PValue, method="BH")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt$PValue, method="BH")

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over1.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.500vol.900or1200counts.obj.Under1Over2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.500vol.900or1200counts.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt$PValue, method="BH")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt$PValue, method="BH")
mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt$p.adj.BH.all = p.adjust(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt$PValue, method="BH")

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over2.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over1.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.300vol.600counts.obj.Under1Over2.pred0.2.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/mecp2Het.CTX.HC.300vol.600counts.pred0.2.Under1Over2.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)


View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="114_Pvalb") & (p.adj.BH.all<0.1), "gene"])
View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="116_Pvalb") & (p.adj.BH.all<0.1), "gene"])
View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="118_Pvalb") & (p.adj.BH.all<0.1), "gene"])
View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="119_Pvalb") & (p.adj.BH.all<0.1), "gene"])


View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="99_Sst") & (p.adj.BH.all<0.1), "gene"])
View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="98_Sst") & (p.adj.BH.all<0.1), "gene"])
View(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(subtype=="119_Pvalb") & (p.adj.BH.all<0.1), "gene"])


Pvalb_subtype_list = unique(CTX_HIP_annot[subclass_label=="Pvalb", cluster_label])
Sst_subtype_list = unique(CTX_HIP_annot[subclass_label=="Sst", cluster_label])
#log fold-change of MR and unchanged genes for subclasses matching our INTACT subclasses; under 1 for KO, over 1 for WT
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(gene %in% gene_panel_classes[PV_MR==TRUE, Gene]) & (subtype %in% Pvalb_subtype_list) & !(is.na(as.numeric(logFC)))], INTACT_label="PV", genetype="PV MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(gene %in% unchanged_genes) & (subtype %in% Pvalb_subtype_list) & !(is.na(as.numeric(logFC)))], INTACT_label="PV", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(gene %in% gene_panel_classes[SST_MR==TRUE, Gene]) & (subtype %in% Sst_subtype_list) & !(is.na(as.numeric(logFC)))], INTACT_label="SST", genetype="SST MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[(gene %in% unchanged_genes) & (subtype %in% Sst_subtype_list) & !(is.na(as.numeric(logFC)))], INTACT_label="SST", genetype="Non-MR, non-MA")
))


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt %>% mutate(INTACT_label = factor(INTACT_label, levels=c("PV", "SST")))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt %>% mutate(genetype = factor(genetype, levels=c("PV MR", "SST MR", "Non-MR, non-MA")))

ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.logfc.agg.PV.SST.dt, aes(x = genetype, y = as.numeric(logFC), fill=genetype))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(aes(middle = mean(as.numeric(logFC))), outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-2,2))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "Gene classes:", values = c("PV MR"=unique(CTX_HIP_annot[subclass_label=="Pvalb", subclass_color]), 
                                                       "SST MR"=unique(CTX_HIP_annot[subclass_label=="Sst", subclass_color]), 
                                                       "Non-MR, non-MA"="gray")) +
  facet_grid(.~INTACT_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.PV.SST.MR.genes.in.PV.SST.cells.smartseqLabels.logFC.sctCounts.Under1Over1.boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.PV.SST.MR.genes.in.PV.SST.cells.smartseqLabels.logFC.sctCounts.Under1Over1.boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2[["predicted.id"]], keep.rownames = "index")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt = data.table(left_join(x=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt, y=CTX_HIP_annot[, .(cluster_label, subclass_label, neighborhood_label)], by=c("predicted.id"="cluster_label")))

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt.celltypes = AddMetaData(object=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, metadata=data.frame(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt, row.names="index"))

subclass_umap_func <- function(Seurat_obj, subclass){
  subclass_obj = subset(Seurat_obj, subset =  subclass_label == subclass)
  subclass_obj <- RunUMAP(subclass_obj, dims = 1:30)
  return(subclass_obj)
}


DimPlot(exp_vizgen.combined.sct.500vol.900counts.Mecp2Reg.cellType.genes.labelTransfer.celltypes, reduction = "umap", group.by = "t.type")+
  ggtitle("500 <= cell volume <= 3*median(cell volume)\n >= 900 counts\nMecp2-regulated genes, cell-type-specific genes")+
  scale_color_manual(values = c("WT"="purple", "KO"="orange", "NA"="grey")) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(size=2), nrow=6))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/exp1_exp2_minVol500_maxVol3Med_900counts_Mecp2Reg_cellType_genes_umap_plot_by_transcriptotype.png",width = 6, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/exp1_exp2_minVol500_maxVol3Med_900counts_Mecp2Reg_cellType_genes_umap_plot_by_transcriptotype.eps",width = 6, height = 5, dpi = 300, units = "in", device='eps')


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb = subclass_umap_func(Seurat_obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt.celltypes, subclass="Pvalb")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst = subclass_umap_func(Seurat_obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.dt.celltypes, subclass="Sst")

DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb, reduction = "umap", group.by = "t.type")+
  ggtitle("Pvalb, labeled by transcriptotype")+
  scale_color_manual(values = c("WT"="purple", "KO"="orange", "NA"="grey")) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(size=2), nrow=6))
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/exp1_exp2_minVol500_maxVol3Med_900counts_Mecp2Reg_cellType_genes_umap_plot_by_transcriptotype_Pvalb.png",width = 6, height = 5, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/exp1_exp2_minVol500_maxVol3Med_900counts_Mecp2Reg_cellType_genes_umap_plot_by_transcriptotype_Pvalb.eps",width = 6, height = 5, dpi = 300, units = "in", device='eps')

DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb, reduction = "umap", group.by = "predicted.id")+
  ggtitle("Pvalb, labeled by transcriptotype")+
  #scale_color_manual(values = c("WT"="purple", "KO"="orange", "NA"="grey")) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(size=2), nrow=6))

DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst, reduction = "umap", group.by = "predicted.id")+
  ggtitle("Sst, labeled by transcriptotype")+
  #scale_color_manual(values = c("WT"="purple", "KO"="orange", "NA"="grey")) +
  theme(plot.title=element_text(hjust=0.5, size=10), legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(size=2), nrow=6))

##overlap of 
fisher.test()
#subtype specific genes defined by Allen Institute in 2018/2019
allen18_pv = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_pv_n_gu_short_1.csv")
allen18_sst = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_sst_n_gu_short_1.csv")
allen18_L5 = fread("HG_lab/Mati/GabelLab/genesets/subtype_specific_genes/allen18_L5_n_gu_short.csv")


mecp2Het.CTX.HC.100vol.300counts.Under1Over1.pred0.2.scaled.DEG.table.Pvalb = scaled.DEG.func(subclass_label="Pvalb", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2)
mecp2Het.CTX.HC.100vol.300counts.Under1Over1.pred0.2.scaled.DEG.table.Sst = scaled.DEG.func(subclass_label="Sst", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2)


Pvalb_test <- acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs(scaled.DEG.table=mecp2Het.CTX.HC.100vol.300counts.Under1Over1.pred0.2.scaled.DEG.table.Pvalb,
                                               subclass_abbr="PV", 
                                               avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt, 
                                               plot_title="", 
                                               sub_dir="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/", 
                                               output_file_prefix="Pvalb_mecp2Het_CTX_HC_100vol_300counts_pred0.2_100vol_300counts_", 
                                               logFC_ymin=-2.5, logFC_ymax=2.5, avgExp_ymin=0, avgExp_ymax=70, save_plot=TRUE)
subclass_abbr="PV"
Pvalb_test2 = data.table(rbind(
  Pvalb_test[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"),],
  Pvalb_test[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")),],
  Pvalb_test[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")),],
  Pvalb_test[(dir=="High") & (gene_class=="Any recurrent MR"),],
  Pvalb_test[(dir=="Low") & (gene_class=="Any recurrent MR"),]))
Pvalb_test2 = Pvalb_test2 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"), "Any recurrent MR")))
Pvalb_test2 = Pvalb_test2 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))

write.table(Pvalb_test2, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/Pvalb.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.PvMR.recurrent.MR.logFC.table.txt", sep="\t", quote=F, row.names=F)

ggplot(Pvalb_test2[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
  ggtitle("Pvalb, area")+
  geom_violin(fill=colors_subclass["Pvalb"], trim=TRUE, scale="area")+
  stat_summary(fun = mean, 
               geom = "point", color="black", size=0.6) + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
  coord_cartesian(ylim=c(-2.5,2.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())

ggplot(Pvalb_test2[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
  ggtitle("Pvalb, count")+
  coord_cartesian(ylim=c(-2.5,2.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  facet_grid(.~gene_class + dir,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  geom_violin(fill=colors_subclass["Pvalb"], trim=TRUE, scale="count")+
  stat_summary(fun = mean, 
               geom = "point", color="black", size=0.6) + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.PvMR.recurrent.MR.logFC.violPlot.mean.se.only.png", width = 6, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.PvMR.recurrent.MR.logFC.violPlot.mean.se.only.eps", width = 6, height = 6, dpi = 300, units = "in", device='eps')


subclass_low_high_nonDEG_plot_func <- function(de_results_table, subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
  if(custom_subtype_list==FALSE){
    subtype_list = CTX_HIP_annot[subclass_label %in% subclass_list, cluster_label]
  } 
  if(custom_subtype_list==TRUE){
    subtype_list = subclass_list
  }
  subtype_logfc_data_table <- data.table(rbind(
    cbind(de_results_table[(gene %in% gene_panel_classes[get(paste0(subclass_abbr,"_MR"))==TRUE, Gene]) & (subtype %in% subtype_list) & is.finite(as.numeric(logFC)), .(subtype, gene, logFC)], gene_class=paste(subclass_abbr, "MR"), dir="All"),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Non-MR, non-MA")) & (dir=="Non-DEG") & is.finite(as.numeric(logFC)), .(subtype, gene, logFC, gene_class, dir)]),
    cbind(low_high_nonDEG_logfc_table[(subtype %in% subtype_list) & (gene_class %in% c("Any recurrent MR") & is.finite(as.numeric(logFC))), .(subtype, gene, logFC, gene_class, dir)])
  ))
  subtype_reduced_list = c()
  for(i in subtype_list){
    if(nrow(unique(subtype_logfc_data_table[subtype==i, .(gene_class, dir)]))==4){
      subtype_reduced_list = c(subtype_reduced_list,i)
    }
  }
  subtype_logfc_data_table = subtype_logfc_data_table[subtype %in% subtype_reduced_list]
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"), "Any recurrent MR")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(dir = factor(dir, levels=c("Non-DEG", "All", "Low", "High")))
  subtype_logfc_data_table = subtype_logfc_data_table %>% mutate(subtype = factor(subtype, levels=c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))))
  
  #jitter plot of log fold changes 
  ggplot(subtype_logfc_data_table, aes(x = subtype, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.6, color="black") + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3, color="black")+
    coord_cartesian(ylim=c(-1.5,1.5))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    scale_color_manual(values = colors_subtype[c(unique(CTX_HIP_annot[cluster_label %in% subtype_logfc_data_table[, subtype], cluster_label]))]) +
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) + 
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(plot_save==TRUE){
    ggsave(paste0(plot_file,".png"), width = 6, height = 6, dpi = 300, units = "in", device='png')
    ggsave(paste0(plot_file,".eps"), width = 6, height = 6, dpi = 300, units = "in", device=cairo_ps)
  }
  if(return_table==TRUE){
    return(subtype_logfc_data_table)
  }
}

head(zscore_exp_acrossSubclass_func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, subclass_names=Rbp4_L5_targets))

#adding
acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR <- function(scaled.DEG.table, avg.exp.table){
  scaled.DEG.table = left_join(scaled.DEG.table, avg.exp.table, by=c("subtype", "gene"))
  gene.high.min =  as.numeric(quantile(scaled.DEG.table[, as.numeric(avgExp)], na.rm=TRUE)[3])
  gene.low.max =  as.numeric(quantile(scaled.DEG.table[, as.numeric(avgExp)], na.rm=TRUE)[3])
  
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
  subtype_zscore_MR_dt_up = subtype_zscore_MR_dt[(dir=="High") & (avgExp > gene.high.min),]
  subtype_zscore_MR_dt_down = subtype_zscore_MR_dt[(dir=="Low") & (avgExp < gene.low.max),]
  subtype_zscore_MR_dt = rbind(subtype_zscore_MR_dt_nonDEG, subtype_zscore_MR_dt_up, subtype_zscore_MR_dt_down)
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "PV MR", "SST MR", "L4 MR", "L5 MR", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC")))
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  return(subtype_zscore_MR_dt[!(is.na(gene_class))])
}


all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt, 
                binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt,
                output_prefix="mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1", plot_save=FALSE)

scaled.DEG.table.L2.3.IT.test  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "L2/3 IT", gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg, de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[,1:7], binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2)

L2.3.IT.dt.test = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.L2.3.IT.test,
                                                                 avg.exp.table = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)

L2.3.IT.dt.test2 = acrossSubclass_yaoDEG_MR_dotplot_func2_reduced_moreArgs(scaled.DEG.table = scaled.DEG.table.L2.3.IT.test,
                                                                           avg.exp.table = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt, 
                                                                           plot_title="L2/3 IT", 
                                                                           sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/L2.3.IT/",
                                                                           output_file_prefix = paste0("L2.3.IT","_",output_prefix), save_plot = FALSE)



all.nhoods.func2 <- function(gene.by.subtype.avg.counts, de.results.table, binary_comparison_table, avg.exp.table){
  scaled.DEG.table.L2.3.IT  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "L2/3 IT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  L2.3.IT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.L2.3.IT,
                                                                       avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.L4.5.6.IT.Car3  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "L4/5/6 IT Car3", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  L4.5.6.IT.Car3.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.L4.5.6.IT.Car3,
                                                                              avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.MGE  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "MGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  MGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.MGE,
                                                                   avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.CGE  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "CGE", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  CGE.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.CGE,
                                                                   avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.PT  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "PT", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  PT.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.PT,
                                                                  avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.NP.CT.L6b  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "NP/CT/L6b", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  NP.CT.L6b.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.NP.CT.L6b,
                                                                         avg.exp.table = avg.exp.table)
  
  scaled.DEG.table.DG.SUB.CA  = group.scaled.DEG.func(grouping="neighborhood", nhood_label = "DG/SUB/CA", gene.by.subtype.avg.counts=gene.by.subtype.avg.counts, de.results.table=de.results.table, binary_comparison_table=binary_comparison_table)
  DG.SUB.CA.dt = acrossSubclass_yaoDEG_MR_dotplot_func_moreArgs_allINTACTMR(scaled.DEG.table = scaled.DEG.table.DG.SUB.CA,
                                                                         avg.exp.table = avg.exp.table)
  
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

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.2 = all.nhoods.func(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                                                                                  de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[, 1:7], 
                                                                                  binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                                                                                  avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt,
                                                                                  output_prefix="mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1", plot_save=FALSE)

all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR = all.nhoods.func2(gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg,
                                                                                  de.results.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.de.results.dt[, 1:7], 
                                                                                  binary_comparison_table=yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2,
                                                                                  avg.exp.table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt)

Pvalb_test2 = data.table(rbind(
  Pvalb_test[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"),],
  Pvalb_test[(dir=="High") & (gene_class==paste(subclass_abbr, "MR")),],
  Pvalb_test[(dir=="Low") & (gene_class==paste(subclass_abbr, "MR")),],
  Pvalb_test[(dir=="High") & (gene_class=="Any recurrent MR"),],
  Pvalb_test[(dir=="Low") & (gene_class=="Any recurrent MR"),]))
Pvalb_test2 = Pvalb_test2 %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", paste(subclass_abbr, "MR"), "Any recurrent MR")))
Pvalb_test2 = Pvalb_test2 %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))


subclass_low_high_nonDEG_violPlot_func <- function(de_results_table, subclass_list, custom_subtype_list=FALSE, subclass_abbr, low_high_nonDEG_logfc_table, ymin=-2.5, ymax=2.5, plot_color="white", plot_title=NULL, plot_save=FALSE, plot_file=NULL, return_table=TRUE){
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
  
  #violin plot of log fold changes 
  ggplot(subtype_logfc_data_table[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC)))+
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
}

Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_violPlot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                                                                                                  subclass_list=c("Pvalb"),
                                                                                                                  subclass_abbr="PV",
                                                                                                                  low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR,
                                                                                                                  plot_color=colors_subclass["Pvalb"],
                                                                                                                  plot_save=TRUE,
                                                                                                                  plot_title="Pvalb\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Pvalb/Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.PvMR.recurrent.MR.logFC.violPlot",
                                                                                                                  return_table=TRUE)

Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_violPlot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                                                                                                           subclass_list=c("Sst"),
                                                                                                                           subclass_abbr="SST",
                                                                                                                           low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR,
                                                                                                                           ymin=-2.5, ymax=2.5,
                                                                                                                           plot_color=colors_subclass["Sst"],
                                                                                                                           plot_save=TRUE,
                                                                                                                           plot_title="Sst\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                                           plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Sst/Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.SstMR.recurrent.MR.logFC.violPlot",
                                                                                                                           return_table=TRUE)

L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 = subclass_low_high_nonDEG_violPlot_func(de_results_table=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.de.results.dt,
                                                                                                               subclass_list=Rbp4_L5_targets,
                                                                                                               custom_subtype_list=TRUE,
                                                                                                               subclass_abbr="L5",
                                                                                                               low_high_nonDEG_logfc_table=all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR,
                                                                                                               ymin=-2.5, ymax=2.5,
                                                                                                               plot_color="#50B2AD",
                                                                                                               plot_save=TRUE,
                                                                                                               plot_title="L5\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2\nUnder1Over1",
                                                                                                               plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/subclass_agg/Rbp4.L5.targets/L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.L5MR.recurrent.MR.logFC.violPlot",
                                                                                                               return_table=TRUE)
L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1_summary <- L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1 %>% group_by(gene_class, dir) %>%
  summarise(
    count = n(),
    mean = mean(as.numeric(logFC), na.rm = TRUE),
    sd = sd(as.numeric(logFC), na.rm = TRUE)
  ) %>% data.table

t.test(x=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="Non-DEG") & (gene_class=="Non-MR, non-MA"), logFC],
       y=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="Low") & (gene_class=="L5 MR"), logFC],
       alternative = "two.sided",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)$p.value

t.test(x=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="Low") & (gene_class=="L5 MR"), logFC],
       y=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="High") & (gene_class=="L5 MR"), logFC],
       alternative = "two.sided",
       conf.level = 0.95)$p.value

t.test(x=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
       y=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="High") & (gene_class=="Any recurrent MR"), logFC],
       alternative = "two.sided",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)$p.value

wilcox.test(x=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="Low") & (gene_class=="Any recurrent MR"), logFC],
            y=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(dir=="High") & (gene_class=="Any recurrent MR"), logFC])$p.value

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

View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1, subclass_abbr="PV"))
View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1, subclass_abbr="SST"))
View(subtype_nonDEG_lowHigh_INTACTMR_anyMeta_sigs(subtype_table=L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1, subclass_abbr="L5"))

nrow(Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="PV MR") & (dir=="Low"),])
nrow(Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="PV MR") & (dir=="High"),])

nrow(Pvalb_nonDEG_lowHigh_PvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="Any recurrent MR") & (dir=="Low"),])
nrow(Pvalb_nonDEG_allPvMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="Any recurrent MR") & (dir=="Low"),])

nrow(Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="Any recurrent MR") & (dir=="Low"),])
nrow(Sst_nonDEG_allSstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="Any recurrent MR") & (dir=="Low"),])



nrow(Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="SST MR") & (dir=="Low"),])
nrow(Sst_nonDEG_lowHigh_SstMR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="SST MR") & (dir=="High"),])

nrow(L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="L5 MR") & (dir=="Low"),])
nrow(L5_nonDEG_lowHigh_L5MR_anyMetaMR_genes_dt_100vol_300counts_pred0.2_Under1Over1[(gene_class=="L5 MR") & (dir=="High"),])

nrow(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[gene_class=="Any recurrent MR"])
nrow(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR[gene_class=="Any recurrent MR"])


setdiff(all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR[gene_class=="Any recurrent MR", paste0(cl,":",gene,":",dir,":",subtype)], all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1[gene_class=="Any recurrent MR", paste0(cl,":",gene,":",dir,":",subtype)])


acrossSubclass_yaoDEG_MR_dotplot_func2_moreArgs <- function(scaled.DEG.table, avg.exp.table, plot_title, sub_dir, output_file_prefix, logFC_ymin=-2.5, logFC_ymax=2.5, avgExp_ymin=0, avgExp_ymax=70, save_plot=FALSE){
  subtype_zscore_MR_dt = data.table(rbind(
    cbind(scaled.DEG.table[(dir=="Non-DEG") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% unchanged_genes)], gene_class="Non-MR, non-MA"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[MR_Meta==TRUE,Gene])], gene_class="Meta MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[Rec_cellspec_MR==TRUE,Gene])], gene_class="Recurrent cell-type-specific MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% any_meta_MR_genes)], gene_class="Any recurrent MR"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% short_lowmC_genes)], gene_class="Short low mC"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% short_lowmC_genes)], gene_class="Short low mC"),
    cbind(scaled.DEG.table[(dir=="High") & (gene %in% gene_panel_classes[PV_MR==TRUE,Gene])], gene_class="PV MR"),
    cbind(scaled.DEG.table[(dir=="Low") & (gene %in% gene_panel_classes[PV_MR==TRUE,Gene])], gene_class="PV MR")))
  subtype_zscore_MR_dt = left_join(x=subtype_zscore_MR_dt, y=avg.exp.table, by=c("subtype", "gene"))
  subtype_zscore_MR_dt_nonDEG = subtype_zscore_MR_dt[(dir=="Non-DEG")]
  subtype_zscore_MR_dt_up = subtype_zscore_MR_dt[(dir=="High") & (avgExp > quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt_down = subtype_zscore_MR_dt[(dir=="Low") & (avgExp < quantile(subtype_zscore_MR_dt[, as.numeric(avgExp)], na.rm=TRUE)[3]),]
  subtype_zscore_MR_dt = rbind(subtype_zscore_MR_dt_nonDEG, subtype_zscore_MR_dt_up, subtype_zscore_MR_dt_down)
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(gene_class = factor(gene_class, levels=c("Non-MR, non-MA", "Meta MR", "Recurrent cell-type-specific MR", "Any recurrent MR", "Short low mC", "PV MR")))
  subtype_zscore_MR_dt = subtype_zscore_MR_dt %>% mutate(dir = factor(dir, levels=c("Non-DEG", "Low", "High")))
  
  p1 <- ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(logFC_ymin,logFC_ymax))+
    ylab("Log2 fold change (KO/WT)") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if(file.exists(sub_dir)){
      p1
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      p1
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_logFC.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
  
  p2 <- ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = dir, y = as.numeric(avgExp), color=subtype))+
    ggtitle(plot_title)+
    geom_jitter(size=0.6, alpha=0.9) +
    stat_summary(fun = mean, 
                 geom = "point", color="black", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, color="black", size=0.3)+
    coord_cartesian(ylim=c(avgExp_ymin,avgExp_ymax))+
    ylab("Average expression") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE))
  if(save_plot==TRUE){
    if(file.exists(sub_dir)){
      p2
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    } else {
      # create a new sub directory inside
      dir.create(sub_dir)
      p2
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.png"), width = 12, height = 6, dpi = 300, units = "in", device='png')
      ggsave(filename=paste0(sub_dir,output_file_prefix,"_subclassAgg_yao2021_DEG_MR_avgExp.eps"), width = 12, height = 6, dpi = 300, units = "in", device=cairo_ps)
    }
  }
  return(subtype_zscore_MR_dt[!(is.na(gene_class))])
  ggplot(subtype_zscore_MR_dt[!(is.na(gene_class))], aes(x = subtype, y = as.numeric(logFC), color=subtype))+
    ggtitle(plot_title)+
    #scale_color_manual(values = colors_nhood[levels(factor(switched.all.nhoods.dt[, neighborhood]))]) +
    stat_summary(fun = mean, 
                 geom = "point", size=0.6) + 
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.4, size=0.3)+
    #coord_cartesian(ylim=c(0,18))+
    ylab("LogFC") + xlab("")+
    facet_grid(.~gene_class + dir,
               switch = "x",
               labeller = label_wrap_gen(width=8),
               scales = "free_x"
    ) + 
    theme_bw()+
    theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
    guides(color=guide_legend(nrow=3, byrow=TRUE))
}

L2.3.IT.dt.test3 = acrossSubclass_yaoDEG_MR_dotplot_func2_moreArgs(scaled.DEG.table = scaled.DEG.table.L2.3.IT.test,
                                                                           avg.exp.table = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.melt, 
                                                                           plot_title="L2/3 IT", 
                                                                           sub_dir ="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/subtype_DE_analysis/yao2021_DEG_plots/neighborhood_agg/L2.3.IT/",
                                                                           output_file_prefix = paste0("L2.3.IT","_",output_prefix), save_plot = FALSE)
nrow(L2.3.IT.dt.test2[gene_class=="Any recurrent MR"])
nrow(L2.3.IT.dt.test3[gene_class=="Any recurrent MR"])

gene.by.subtype.avg.counts=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg
subclass_names="Pvalb"
subclass_subtypes = cluster_subclass[subclass_label%in%subclass_names, cluster_label]
subclass_subtypes_in_data = intersect(colnames(gene.by.subtype.avg.counts), subclass_subtypes)
gene.by.subtype.avg.counts.subclass = gene.by.subtype.avg.counts[, c(subclass_subtypes_in_data)]
gene.by.subtype.avg.counts.subclass.scale = scale(t(gene.by.subtype.avg.counts.subclass))

#subclass-level pseudobulkDGE
WT.KO.subclass.de.noshrink.func <- function(seurat_obj){
  sce <- subset(seurat_obj, subset = ((t.type == "WT") | (t.type == "KO")))
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                              colData = sce@meta.data)
  groups.sce <- colData(sce)[, c("subclass_label", "t.type", "rep")]
  aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)
  aggr_counts.filt <- aggr_counts[,aggr_counts$ncells >= 10]
  col.data = colData(aggr_counts.filt)
  col.data$t.type = factor(col.data$t.type, levels=c("WT", "KO"))
  de.results.noShrink <- pseudoBulkDGE(aggr_counts.filt, 
                                       label=aggr_counts.filt$subclass_label,
                                       col.data=col.data,
                                       design=~ rep + t.type,
                                       coef="t.typeKO",
                                       condition=aggr_counts.filt$t.type,
                                       robust=FALSE)
  return(de.results.noShrink)
}
#processing DEG result file to a data table
subclass_de_results_func <- function(de.results){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = "subclass"))
  names(de.results.dt)[2] = "gene"
  de.results.dt$p.adj.BH.all = p.adjust(de.results.dt$PValue, method="BH")
  return(de.results.dt)
}
#
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2[["predicted.id"]], keep.rownames = "index")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt = data.table(left_join(x=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt, y=CTX_HIP_annot[, .(cluster_label, subclass_label, neighborhood_label)], by=c("predicted.id"="cluster_label")))

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = AddMetaData(object=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2, metadata=data.frame(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt, row.names="index"))
saveRDS(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results <- WT.KO.subclass.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.dt <- subclass_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results)


write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/subclass.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.summ <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.dt %>% group_by(subclass) %>%
  summarise(
    num_DEGs_padj0.1 = sum((as.numeric(p.adj.BH.all) < 0.1), na.rm = TRUE),
  ) %>% data.table


ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.summ, aes(x = subclass, y = as.numeric(num_DEGs_padj0.1)))+
  ggtitle("100 <= cell volume <= 3*median(cell volume)\ncounts >= 300\nUnder 1, Over 1")+
  geom_col()+
  #coord_cartesian(ylim=c(-1,1))+
  ylab("Number of DEGs") + xlab("")+
  #scale_fill_manual(name = "Genes:", values = c("Non-MR, non-MA"="gray", "Meta MR"="red")) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10, angle=90), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/barplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.subclasses.numDEGs.padj0.1.boxplot.png", width = 10, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/barplots/mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.subclasses.numDEGs.padj0.1.boxplot.eps", width = 10, height = 5, dpi = 300, units = "in", device='eps')


#custom label for aggregating for pseudobulkDGE
WT.KO.custom.de.noshrink.func <- function(seurat_obj){
  sce <- subset(seurat_obj, subset = ((t.type == "WT") | (t.type == "KO")))
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                              colData = sce@meta.data)
  groups.sce <- colData(sce)[, c("custom_label", "t.type", "rep")]
  aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)
  aggr_counts.filt <- aggr_counts[,aggr_counts$ncells >= 10]
  col.data = colData(aggr_counts.filt)
  col.data$t.type = factor(col.data$t.type, levels=c("WT", "KO"))
  de.results.noShrink <- pseudoBulkDGE(aggr_counts.filt, 
                                       label=aggr_counts.filt$custom_label,
                                       col.data=col.data,
                                       design=~ rep + t.type,
                                       coef="t.typeKO",
                                       condition=aggr_counts.filt$t.type,
                                       robust=FALSE)
  return(de.results.noShrink)
}
#processing DEG result file to a data table
custom_de_results_func <- function(de.results){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = "custom_label"))
  names(de.results.dt)[2] = "gene"
  de.results.dt$p.adj.BH.all = p.adjust(de.results.dt$PValue, method="BH")
  return(de.results.dt)
}
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = predicted.id %in% Rbp4_L5_targets)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets$custom_label <- "Rbp4_L5_target"

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results <- WT.KO.custom.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt <- custom_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = predicted.id %in% Nr5a1_L4_targets)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets$custom_label <- "Nr5a1_L4_target"

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results <- WT.KO.custom.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt <- custom_de_results_func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.dt[subclass=="Pvalb"]
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.de.results.dt[subclass=="Sst"]
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt)[1] <- "custom_label"
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt)[1] <- "custom_label"

write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Rbp4.L5.targets.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Nr5a1.L4.targets.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Pvalb.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Sst.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)

#read the subclass-level pseudobulkDGE results in
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Rbp4.L5.targets.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Nr5a1.L4.targets.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Pvalb.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/Sst.INTACT.agg.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")



#log fold-change of MR and unchanged genes for subclasses matching our INTACT subclasses
logfc_agg_intact_subclass_dt_CTX_HC <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt[(gene %in% gene_panel_classes[PV_MR==TRUE, Gene])], INTACT_label="PV", genetype="PV MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="PV", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt[(gene %in% gene_panel_classes[SST_MR==TRUE, Gene])], INTACT_label="SST", genetype="SST MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="SST", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt[(gene %in% gene_panel_classes[L4_MR==TRUE, Gene])], INTACT_label="L4", genetype="L4 MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="L4", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt[(gene %in% gene_panel_classes[L5_MR==TRUE, Gene])], INTACT_label="L5", genetype="L5 MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="L5", genetype="Non-MR, non-MA")
))
logfc_agg_intact_subclass_dt_CTX_HC = logfc_agg_intact_subclass_dt_CTX_HC %>% mutate(INTACT_label = factor(INTACT_label, levels=c("PV", "SST","L4", "L5")))
logfc_agg_intact_subclass_dt_CTX_HC = logfc_agg_intact_subclass_dt_CTX_HC %>% mutate(genetype = factor(genetype, levels=c("PV MR", "SST MR", "L4 MR", "L5 MR", "Non-MR, non-MA")))

ggplot(logfc_agg_intact_subclass_dt_CTX_HC[!is.na(logFC)], aes(x = genetype, y = as.numeric(logFC), fill=genetype))+
  ggtitle("INTACT-level pseudobulkDGE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(aes(middle = mean(as.numeric(logFC))), outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "", values = c("PV MR"=unique(CTX_HIP_annot[subclass_label=="Pvalb", subclass_color]), 
                                          "SST MR"=unique(CTX_HIP_annot[subclass_label=="Sst", subclass_color]), 
                                          "L4 MR"="#00E5E5",
                                          "L5 MR"="#50B2AD",
                                          "Non-MR, non-MA"="gray")) +
  facet_grid(.~INTACT_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.in.PV.SST.L4.L5.cells.smartseqLabels.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.png", width = 5, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.in.PV.SST.L4.L5.cells.smartseqLabels.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.eps", width = 5, height = 6, dpi = 300, units = "in", device='eps')

write.table(logfc_agg_intact_subclass_dt_CTX_HC, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/logfc_agg_intact_subclass_CTX_HC_table.txt", quote=F, row.names=F, sep="\t")
ggplot(logfc_agg_intact_subclass_dt_CTX_HC[!is.na(logFC)], aes(x = genetype, y = as.numeric(logFC), fill=genetype))+
  ggtitle("INTACT-level pseudobulkDGE")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(aes(middle = mean(as.numeric(logFC))), outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1,1))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "", values = c("PV MR"=unique(CTX_HIP_annot[subclass_label=="Pvalb", subclass_color]), 
                                          "SST MR"=unique(CTX_HIP_annot[subclass_label=="Sst", subclass_color]), 
                                          "L4 MR"="#00E5E5",
                                          "L5 MR"="#50B2AD",
                                          "Non-MR, non-MA"="gray")) +
  facet_grid(.~INTACT_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.eps", width = 4, height = 6, dpi = 300, units = "in", device='eps')

head(logfc_agg_intact_subclass_dt_CTX_HC)
pval1 <- wilcox.test(x=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="PV") & (genetype=="PV MR"), logFC],
                     y=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="PV") & (genetype=="Non-MR, non-MA"), logFC])$p.value
sig1 <- sig_function(pval1)

INTACTMR_nonMR_sigs <- function(logfc_agg_intact_subclass_dt_CTX_HC){
  pval1 <- wilcox.test(x=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="PV") & (genetype=="PV MR"), logFC],
                       y=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="PV") & (genetype=="Non-MR, non-MA"), logFC])$p.value
  sig1 <- sig_function(pval1)
  
  pval2 <- wilcox.test(x=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="SST") & (genetype=="SST MR"), logFC],
                       y=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="SST") & (genetype=="Non-MR, non-MA"), logFC])$p.value
  sig2 <- sig_function(pval2)
  
  pval3 <- wilcox.test(x=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="L4") & (genetype=="L4 MR"), logFC],
                       y=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="L4") & (genetype=="Non-MR, non-MA"), logFC])$p.value
  sig3 <- sig_function(pval3)
  
  pval4 <- wilcox.test(x=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="L5") & (genetype=="L5 MR"), logFC],
                       y=logfc_agg_intact_subclass_dt_CTX_HC[(INTACT_label=="L5") & (genetype=="Non-MR, non-MA"), logFC])$p.value
  sig4 <- sig_function(pval4)
  
  sigs <- data.table(rbind(cbind(compar="PV MR; Non-MR, non-MA", sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(compar="SST MR; Non-MR, non-MA", sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(compar="L4 MR; Non-MR, non-MA", sig_symbol=sig3, wilcox.pval=pval3),
                           cbind(compar="L5 MR; Non-MR, non-MA", sig_symbol=sig4, wilcox.pval=pval4)))
  return(sigs)
}

INTACTMR_nonMR_sig_table <- INTACTMR_nonMR_sigs(logfc_agg_intact_subclass_dt_CTX_HC=logfc_agg_intact_subclass_dt_CTX_HC)

#meta MR genes
metaMR_logfc_agg_intact_subclass_dt_CTX_HC <- data.table(rbind(
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt[(gene %in% gene_panel_classes[MR_Meta==TRUE, Gene])], INTACT_label="PV", genetype="Meta MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Pvalb.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="PV", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt[(gene %in% gene_panel_classes[MR_Meta==TRUE, Gene])], INTACT_label="SST", genetype="Meta MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Sst.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="SST", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt[(gene %in% gene_panel_classes[MR_Meta==TRUE, Gene])], INTACT_label="L4", genetype="Meta MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Nr5a1.L4.targets.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="L4", genetype="Non-MR, non-MA"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt[(gene %in% gene_panel_classes[MR_Meta==TRUE, Gene])], INTACT_label="L5", genetype="Meta MR"),
  cbind(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.Rbp4.L5.targets.de.results.dt[(gene %in% unchanged_genes)], INTACT_label="L5", genetype="Non-MR, non-MA")
))
metaMR_logfc_agg_intact_subclass_dt_CTX_HC = metaMR_logfc_agg_intact_subclass_dt_CTX_HC %>% mutate(INTACT_label = factor(INTACT_label, levels=c("PV", "SST","L4", "L5")))
metaMR_logfc_agg_intact_subclass_dt_CTX_HC = metaMR_logfc_agg_intact_subclass_dt_CTX_HC %>% mutate(genetype = factor(genetype, levels=c("Non-MR, non-MA", "Meta MR")))

ggplot(metaMR_logfc_agg_intact_subclass_dt_CTX_HC[!is.na(logFC)], aes(x = genetype, y = as.numeric(logFC), fill=genetype))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(aes(middle = mean(as.numeric(logFC))), outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-1.5,1.5))+
  ylab("Log2 fold change (KO/WT)") + xlab("")+
  scale_fill_manual(name = "", values = c("Non-MR, non-MA"="gray", "Meta MR"="red")) +
  facet_grid(.~INTACT_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank()) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.metaMR.genes.in.INTACT.cells.INTACTLevel.pseudobulkDGE.logFC.boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.metaMR.genes.in.INTACT.cells.INTACTLevel.pseudobulkDGE.logFC.boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


#reps 1 and 2 separated from reps 6 and 7, subtype pseudobulkDGE
exp1.exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = (rep == "Rep1") | (rep == "Rep2"))
exp6.exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = (rep == "Rep6") | (rep == "Rep7"))


exp1.exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results <- WT.KO.de.noshrink.func(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = (rep == "Rep1") | (rep == "Rep2")))
exp1.exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results.dt <- subtype_de_results_func(exp1.exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results)

exp6.exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results <- WT.KO.de.noshrink.func(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = (rep == "Rep6") | (rep == "Rep7")))
exp6.exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results.dt <- subtype_de_results_func(exp6.exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results)

write.csv(exp1.exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/subtype.agg.exp1.exp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)
write.csv(exp6.exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subtype.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/subtype.agg.exp6.exp7.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv", quote=F, row.names=F)


