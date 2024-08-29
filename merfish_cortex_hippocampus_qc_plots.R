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
library(patchwork)
#cpm
Pv_raw = data.frame(read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/pv_ko_nondedup_exon_all_rawcounts.tsv', header=TRUE), row.names="Gene")
Pv_raw = Pv_raw[, 7:10]
Pv_cpm = cpm(Pv_raw)
Pv_cpm_means = melt(rowMeans(Pv_cpm))
Pv_cpm_means = data.table(Pv_cpm_means, keep.rownames = "Gene")

#Sst
Sst_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/sst_ko_nondedup_exon_rawcounts_011722.tsv', header=TRUE)
Sst_raw = Sst_raw[,5:8]
Sst_cpm = cpm(Sst_raw)
Sst_cpm_means = melt(rowMeans(Sst_cpm))
Sst_cpm_means = data.table(Sst_cpm_means, keep.rownames = "Gene")

#L4
L4_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/Nr5a1_ko_nondedup_exon_rawcounts.tsv', header=TRUE)
L4_raw = L4_raw[,5:8]
L4_cpm = cpm(L4_raw)
L4_cpm_means = melt(rowMeans(L4_cpm))
L4_cpm_means = data.table(L4_cpm_means, keep.rownames = "Gene")

#L5
L5_raw = read.table('HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/rbp4_all_exon_nondedup_counts.tsv')
#use WT1, WT2, WT4, and WT5 (not WT3)
L5_raw = L5_raw[, c(6,7,9,10)]
L5_cpm = cpm(L5_raw)
L5_cpm_means = melt(rowMeans(L5_cpm))
L5_cpm_means = data.table(L5_cpm_means, keep.rownames = "Gene")

#Yao 2021 metadata
yao_meta_dt = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/metadata.csv", na.strings=c("", "NA"))
yao_meta_dt = yao_meta_dt[!is.na(cluster_label), ]
#subclass assignments of cluster labels from Yao 2021
cluster_subclass = unique(yao_meta_dt[, .(cluster_label, subclass_label)])

CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
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


#mecp2Het.CTX.HC.100vol.300counts.pred0.2.obj.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")


exp1.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep=="Rep1"))
exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep=="Rep2"))
exp6.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep=="Rep6"))
exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = colnames(subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep=="Rep7"))


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.ttypes = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[c("t.type", "rep")]], keep.rownames="index")
WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.ttypes[t.type=="WT", index]
KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.ttypes[t.type=="KO", index]
NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.ttypes[is.na(t.type), index]


#order the cells so that half the exp1 cells will be plotted, then half of exp2, then remaining half of exp1, then remaining half of exp2
rep_cell_order = c(exp1.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:9933], 
               exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:17133],
               exp6.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:16349], 
               exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:14777],
               exp1.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[9934:19867], 
               exp2.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[17134:34267],
               exp6.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[16350:32698], 
               exp7.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[14778:29555])

ttype_cell_order = c(NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:6413],
                     NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[6414:12826],
                     KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:24013],
                     KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[24014:48027],
                     WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:27767],
                     WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[27768:55534])

#ttype_cell_order2 = c(WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:52757],
 #                     KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:45625],
  #                    NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:12184],
   #                   WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[52758:55534],
    #                  KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[45626:48027],
     #                 NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[12185:12826])

ttype_cell_order2 = c(WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:52757],
                      KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:45625],
                      NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[1:12820],
                      WT.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[52758:55534],
                      KO.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[45626:48027],
                      NA.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[12821:12826])
#set.seed(2)
#mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.cellOrder = sample(colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes))

DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, reduction = "umap", group.by = "rep", raster=FALSE)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2")+
  #scale_color_manual(values = c("KO"="orange", "WT"="purple")) +
  theme(plot.title=element_text(hjust=0.5, size=5), legend.position = "bottom")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)



mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.nonNA <- subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = ((t.type=="WT") | (t.type=="KO")))

DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.nonNA, reduction = "umap", group.by = "t.type", order = ttype_cell_order2, raster=FALSE)+
  ggtitle("Mecp2 KO/+ cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore > 0.2")+
  scale_color_manual(values = c("KO"="orange", "WT"="purple")) +
  theme(plot.title=element_text(hjust=0.5, size=5), legend.position = "bottom")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.transcriptotype.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/umaps/CTX_HC/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.transcriptotype.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, reduction = "umap", group.by = "t.type", raster=FALSE, shuffle=TRUE)+
  ggtitle("")+
  scale_color_manual(values = c("KO"="orange", "WT"="purple", "NA"="gray")) +
  theme(plot.title=element_text(hjust=0.5, size=5), legend.position = "bottom")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.transcriptotype.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.transcriptotype.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


DimPlot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, reduction = "umap", group.by = "rep", raster=FALSE, shuffle=TRUE)+
  ggtitle("")+
  theme(plot.title=element_text(hjust=0.5, size=5), legend.position = "bottom")
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.rep.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.Under1Over1.UMAP.by.rep.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


#all reps

DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.allReps.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.allReps.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.allReps.cpmNorm, subset=t.type=="WT")
sce.allReps.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.allReps.cpmNorm.WT[["RNA"]], slot = "data")), 
                                    colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.allReps.cpmNorm.WT@meta.data)



#rep1
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1 = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep == "Rep1")
DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1.cpmNorm, subset=t.type=="WT")

#rep2
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2 = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep == "Rep2")
DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2.cpmNorm, subset=t.type=="WT")


#rep6
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6 = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep == "Rep6")
DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6.cpmNorm, subset=t.type=="WT")


#rep7
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7 = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = rep == "Rep7")
DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7.cpmNorm, subset=t.type=="WT")




sce.exp1.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1.cpmNorm.WT[["RNA"]], slot = "data")), 
                                    colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep1.cpmNorm.WT@meta.data)
groups.exp1.WT <- colData(sce.exp1.WT)[, c("predicted.id")]
agg.exp1.WT.avgCPM <- aggregateAcrossCells(sce.exp1.WT , ids=groups.exp1.WT, statistics="mean")
agg.exp1.WT.avgCPM.dt = data.table(assay(agg.exp1.WT.avgCPM), keep.rownames="Gene")
#exp2 averaging cpm
sce.exp2.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2.cpmNorm.WT[["RNA"]], slot = "data")), 
                                    colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep2.cpmNorm.WT@meta.data)
groups.exp2.WT <- colData(sce.exp2.WT)[, c("predicted.id")]
agg.exp2.WT.avgCPM <- aggregateAcrossCells(sce.exp2.WT , ids=groups.exp2.WT, statistics="mean")
agg.exp2.WT.avgCPM.dt = data.table(assay(agg.exp2.WT.avgCPM), keep.rownames="Gene")
#exp6 averaging cpm
sce.exp6.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6.cpmNorm.WT[["RNA"]], slot = "data")), 
                                    colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep6.cpmNorm.WT@meta.data)
groups.exp6.WT <- colData(sce.exp6.WT)[, c("predicted.id")]
agg.exp6.WT.avgCPM <- aggregateAcrossCells(sce.exp6.WT , ids=groups.exp6.WT, statistics="mean")
agg.exp6.WT.avgCPM.dt = data.table(assay(agg.exp6.WT.avgCPM), keep.rownames="Gene")
#exp7 averaging cpm
sce.exp7.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7.cpmNorm.WT[["RNA"]], slot = "data")), 
                                    colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.rep7.cpmNorm.WT@meta.data)
groups.exp7.WT <- colData(sce.exp7.WT)[, c("predicted.id")]
agg.exp7.WT.avgCPM <- aggregateAcrossCells(sce.exp7.WT , ids=groups.exp7.WT, statistics="mean")
agg.exp7.WT.avgCPM.dt = data.table(assay(agg.exp7.WT.avgCPM), keep.rownames="Gene")



##rep1 vs rep2
#
#averaging CPM over all cells, rep1
groups.exp1.WT.allCells <- colData(sce.exp1.WT)[, c("rep")]
agg.exp1.WT.avgCPM.allCells <- aggregateAcrossCells(sce.exp1.WT, ids=groups.exp1.WT.allCells, statistics="mean")
agg.exp1.WT.avgCPM.allCells.dt = data.table(assay(agg.exp1.WT.avgCPM.allCells), keep.rownames="Gene")
#averaging CPM over all cells, rep2
groups.exp2.WT.allCells <- colData(sce.exp2.WT)[, c("rep")]
agg.exp2.WT.avgCPM.allCells <- aggregateAcrossCells(sce.exp2.WT, ids=groups.exp2.WT.allCells, statistics="mean")
agg.exp2.WT.avgCPM.allCells.dt = data.table(assay(agg.exp2.WT.avgCPM.allCells), keep.rownames="Gene")
#averaging CPM over all cells, rep6
groups.exp6.WT.allCells <- colData(sce.exp6.WT)[, c("rep")]
agg.exp6.WT.avgCPM.allCells <- aggregateAcrossCells(sce.exp6.WT, ids=groups.exp6.WT.allCells, statistics="mean")
agg.exp6.WT.avgCPM.allCells.dt = data.table(assay(agg.exp6.WT.avgCPM.allCells), keep.rownames="Gene")
#averaging CPM over all cells, rep7
groups.exp7.WT.allCells <- colData(sce.exp7.WT)[, c("rep")]
agg.exp7.WT.avgCPM.allCells <- aggregateAcrossCells(sce.exp7.WT, ids=groups.exp7.WT.allCells, statistics="mean")
agg.exp7.WT.avgCPM.allCells.dt = data.table(assay(agg.exp7.WT.avgCPM.allCells), keep.rownames="Gene")



agg.exp1.exp2.WT.avgCPM.allCells = inner_join(x=agg.exp1.WT.avgCPM.allCells.dt, y=agg.exp2.WT.avgCPM.allCells.dt, by="Gene")
agg.exp1.exp6.WT.avgCPM.allCells = inner_join(x=agg.exp1.WT.avgCPM.allCells.dt, y=agg.exp6.WT.avgCPM.allCells.dt, by="Gene")
agg.exp1.exp7.WT.avgCPM.allCells = inner_join(x=agg.exp1.WT.avgCPM.allCells.dt, y=agg.exp7.WT.avgCPM.allCells.dt, by="Gene")

agg.exp2.exp6.WT.avgCPM.allCells = inner_join(x=agg.exp2.WT.avgCPM.allCells.dt, y=agg.exp6.WT.avgCPM.allCells.dt, by="Gene")
agg.exp2.exp7.WT.avgCPM.allCells = inner_join(x=agg.exp2.WT.avgCPM.allCells.dt, y=agg.exp7.WT.avgCPM.allCells.dt, by="Gene")

agg.exp6.exp7.WT.avgCPM.allCells = inner_join(x=agg.exp6.WT.avgCPM.allCells.dt, y=agg.exp7.WT.avgCPM.allCells.dt, by="Gene")


agg.exp6.exp7.WT.avgCPM.allCells = inner_join(x=agg.exp6.WT.avgCPM.allCells.dt, y=agg.exp7.WT.avgCPM.allCells.dt, by="Gene")
agg.exp1.exp6.WT.avgCPM.allCells = inner_join(x=agg.exp1.WT.avgCPM.allCells.dt, y=agg.exp6.WT.avgCPM.allCells.dt, by="Gene")


cor12 = cor(agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep2)], method="spearman", use="complete.obs")
cor16 = cor(agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep6)], method="spearman", use="complete.obs")
cor17 = cor(agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp7.WT.avgCPM.allCells[, log2(Rep7)], method="spearman", use="complete.obs")

cor26 = cor(agg.exp2.exp6.WT.avgCPM.allCells[, log2(Rep2)], agg.exp2.exp6.WT.avgCPM.allCells[, log2(Rep6)], method="spearman", use="complete.obs")
cor27 = cor(agg.exp2.exp7.WT.avgCPM.allCells[, log2(Rep2)], agg.exp2.exp7.WT.avgCPM.allCells[, log2(Rep7)], method="spearman", use="complete.obs")

cor67 = cor(agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep6)], agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep7)], method="spearman", use="complete.obs")

data=c(1, cor12, cor16, cor17, cor12, 1, cor26, cor27, cor16, cor26, 1, cor67, cor17, cor27, cor67, 1)
data_matrix <- matrix(data=data, nrow=4, ncol=4)

library(gplots)
corr_palette <- colorRampPalette(c("white", "red", "darkred"))(n = 20)
#col_breaks = c(seq(0.90,0.949,length=10),  # for white
#               seq(0.95,1,length=10))

library(RColorBrewer)
brewer.pal(5, "Reds")
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/mecp2.CTX.HC.100vol.300counts.pred0.2.experiment.log2cpm.corr.matrix.heatmap.eps")
heatmap.2(data_matrix,
          main = "Correlation matrix", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=seq(0.9,1,0.005),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=corr_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, 
          asp=1) 
dev.off()


#scatterplot showing average CPM in each gene in rep1 and rep2
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.exp1.vs.exp2.avgCPMs.scatterplot.png",  width=1600, height=1600, res=300)
plot(agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep2)], xlab="MERFISH rep1 log2(average CPM)", ylab="MERFISH rep2 log2(average CPM)", pch=16)
legend("topleft", legend = paste0("r=", round(cor(agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep2)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.exp1.vs.exp2.avgCPMs.scatterplot.eps")
plot(agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep2)], xlab="MERFISH rep1 log2(average CPM)", ylab="MERFISH rep2 log2(average CPM)", pch=16)
legend("topleft", legend = paste0("r=", round(cor(agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp2.WT.avgCPM.allCells[, log2(Rep2)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.exp1.vs.exp6.avgCPMs.scatterplot.eps")
plot(agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep6)], xlab="MERFISH rep1 log2(average CPM)", ylab="MERFISH rep6 log2(average CPM)", pch=16)
legend("topleft", legend = paste0("r=", round(cor(agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep1)], agg.exp1.exp6.WT.avgCPM.allCells[, log2(Rep6)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.exp6.vs.exp7.avgCPMs.scatterplot.eps")
plot(agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep6)], agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep7)], xlab="MERFISH rep6 log2(average CPM)", ylab="MERFISH rep7 log2(average CPM)", pch=16)
legend("topleft", legend = paste0("r=", round(cor(agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep6)], agg.exp6.exp7.WT.avgCPM.allCells[, log2(Rep7)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

##subset SingleCellExperiment for chosen L5 subtypes


#Rbp4 L5 targets
sce.allReps.WT.Rbp4.targets = sce.allReps.WT[, sce.allReps.WT$predicted.id %in% Rbp4_L5_targets]
sce.allReps.WT.Rbp4.targets$target = "Rbp4_L5_target"
groups.allReps.WT.Rbp4.targets <- colData(sce.allReps.WT.Rbp4.targets)[, c("target")]
agg.allReps.WT.avgCPM.Rbp4.targets <- aggregateAcrossCells(sce.allReps.WT.Rbp4.targets , ids=groups.allReps.WT.Rbp4.targets, statistics="mean")
agg.allReps.WT.avgCPM.Rbp4.targets.dt = data.table(assay(agg.allReps.WT.avgCPM.Rbp4.targets), keep.rownames="Gene")


L5_cpm_means_join = inner_join(x=L5_cpm_means, y=agg.allReps.WT.avgCPM.Rbp4.targets.dt, by="Gene")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/L5.INTACT.vs.merfish.allReps.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.avgCPMs.scatterplot.eps")
plot(L5_cpm_means_join[, log2(value)], L5_cpm_means_join[, log2(Rbp4_L5_target)],xlab="L5 INTACT log2(mean CPM)", ylab="L5 MERFISH log2(mean CPM)", pch=16, col="#50B2AD")
legend("topleft", legend = paste0("r=", round(cor(L5_cpm_means_join[, log2(value)], L5_cpm_means_join[, log2(Rbp4_L5_target)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#Nr5a1 L4 targets
sce.allReps.WT.Nr5a1.targets = sce.allReps.WT[, sce.allReps.WT$predicted.id %in% Nr5a1_L4_targets]
sce.allReps.WT.Nr5a1.targets$target = "Nr5a1_L4_target"
groups.allReps.WT.Nr5a1.targets <- colData(sce.allReps.WT.Nr5a1.targets)[, c("target")]
agg.allReps.WT.avgCPM.Nr5a1.targets <- aggregateAcrossCells(sce.allReps.WT.Nr5a1.targets , ids=groups.allReps.WT.Nr5a1.targets, statistics="mean")
agg.allReps.WT.avgCPM.Nr5a1.targets.dt = data.table(assay(agg.allReps.WT.avgCPM.Nr5a1.targets), keep.rownames="Gene")


L4_cpm_means_join = inner_join(x=L4_cpm_means, y=agg.allReps.WT.avgCPM.Nr5a1.targets.dt, by="Gene")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/L4.INTACT.vs.merfish.allReps.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.avgCPMs.scatterplot.eps")
plot(L4_cpm_means_join[, log2(value)], L4_cpm_means_join[, log2(Nr5a1_L4_target)],xlab="L4 INTACT log2(mean CPM)", ylab="L4 MERFISH log2(mean CPM)", pch=16, col="#00E5E5")
legend("topleft", legend = paste0("r=", round(cor(L4_cpm_means_join[, log2(value)], L4_cpm_means_join[, log2(Nr5a1_L4_target)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#Pvalb targets
sce.allReps.WT.Pvalb = sce.allReps.WT[,sce.allReps.WT$subclass_label=="Pvalb"]
sce.allReps.WT.Pvalb$target = "Pvalb_target"
groups.allReps.WT.Pvalb <- colData(sce.allReps.WT.Pvalb)[, c("target")]
agg.allReps.WT.avgCPM.Pvalb <- aggregateAcrossCells(sce.allReps.WT.Pvalb , ids=groups.allReps.WT.Pvalb, statistics="mean")
agg.allReps.WT.avgCPM.Pvalb.dt = data.table(assay(agg.allReps.WT.avgCPM.Pvalb), keep.rownames="Gene")

Pv_cpm_means_join = inner_join(x=Pv_cpm_means, y=agg.allReps.WT.avgCPM.Pvalb.dt, by="Gene")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/PV.INTACT.vs.merfish.allReps.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.avgCPMs.scatterplot.eps")
plot(Pv_cpm_means_join[, log2(value)], Pv_cpm_means_join[, log2(Pvalb_target)],xlab="PV INTACT log2(mean CPM)", ylab="PV MERFISH log2(mean CPM)", pch=16, col=unique(CTX_HIP_annot[subclass_label=="Pvalb", subclass_color]))
legend("topleft", legend = paste0("r=", round(cor(Pv_cpm_means_join[, log2(value)], Pv_cpm_means_join[, log2(Pvalb_target)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

#Sst targets
sce.allReps.WT.Sst = sce.allReps.WT[,sce.allReps.WT$subclass_label=="Sst"]
sce.allReps.WT.Sst$target = "Sst_target"
groups.allReps.WT.Sst <- colData(sce.allReps.WT.Sst)[, c("target")]
agg.allReps.WT.avgCPM.Sst <- aggregateAcrossCells(sce.allReps.WT.Sst , ids=groups.allReps.WT.Sst, statistics="mean")
agg.allReps.WT.avgCPM.Sst.dt = data.table(assay(agg.allReps.WT.avgCPM.Sst), keep.rownames="Gene")

Sst_cpm_means_join = inner_join(x=Sst_cpm_means, y=agg.allReps.WT.avgCPM.Sst.dt, by="Gene")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/scatterplots/CTX_HC/SST.INTACT.vs.merfish.allReps.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.avgCPMs.scatterplot.eps")
plot(Sst_cpm_means_join[, log2(value)], Sst_cpm_means_join[, log2(Sst_target)],xlab="SST INTACT log2(mean CPM)", ylab="SST MERFISH log2(mean CPM)", pch=16, col=unique(CTX_HIP_annot[subclass_label=="Sst", subclass_color]))
legend("topleft", legend = paste0("r=", round(cor(Sst_cpm_means_join[, log2(value)], Sst_cpm_means_join[, log2(Sst_target)], method="spearman", use="complete.obs"), 3)), bty = "n", text.col="black")
dev.off()

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[c("predicted.id", "subclass_label", "neighborhood_label")]], keep.rownames = "index")

#summarizing the data table for calculating proportions of each cell type
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt %>% 
                                                 group_by(subclass_label, ) %>%
                                                 summarise(n.val = n()) %>%
                                                 mutate(prop = n.val / sum(n.val)))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[, prop:=as.numeric(prop)]
#plotting the proportions of cell types
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/cluster_proportion_plots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary.barplot.png", width=3000, height=2000, res=300)
par(mfcol=c(3,1))
barplot(t(as.matrix(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[1:14, c("subclass_label", "prop")], rownames="subclass_label")),
        names.arg=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[1:14, subclass_label],
        ylim=c(0,0.12), col=c("gray"), las=2, cex.axis=0.5, cex.names=0.6, cex.main=0.8,
        main="",
        ylab="Proportion",
        legend = FALSE)
barplot(t(as.matrix(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[15:28, c("subclass_label", "prop")], rownames="subclass_label")),
        names.arg=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[15:28, subclass_label],
        ylim=c(0,0.12), col=c("gray"), las=2, cex.axis=0.5, cex.names=0.6, cex.main=0.8,
        main="",
        ylab="Proportion",
        legend = FALSE)
barplot(t(as.matrix(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[29:42, c("subclass_label", "prop")], rownames="subclass_label")),
        names.arg=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.subclass.summary[29:42, subclass_label],
        ylim=c(0,0.12), col=c("gray"), las=2, cex.axis=0.5, cex.names=0.6, cex.main=0.8,
        main="",
        ylab="Proportion",
        legend = FALSE)
dev.off()





##summary tables of sctransform-corrected counts and other metadata
#100 vol, 300 count
exp1.CTX.HC.100vol.300counts.pred0.2.sctTable = fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp2.CTX.HC.100vol.300counts.pred0.2.sctTable = fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp6.CTX.HC.100vol.300counts.pred0.2.sctTable = fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp6.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp7.CTX.HC.100vol.300counts.pred0.2.sctTable = fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp7.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")

#WT
exp3.100vol.300counts.pred0.2.sctTable <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/exp3.100vol.300counts.pred0.2.sctCount.celltype.summary.table.csv")

#combining the replicates
WTandHet.100vol.300counts.pred0.2.sctTable <- rbind(exp1.CTX.HC.100vol.300counts.pred0.2.sctTable,
                                                           exp2.CTX.HC.100vol.300counts.pred0.2.sctTable,
                                                           exp6.CTX.HC.100vol.300counts.pred0.2.sctTable,
                                                           exp7.CTX.HC.100vol.300counts.pred0.2.sctTable,
                                                           exp3.100vol.300counts.pred0.2.sctTable)
#histograms of
#hist(mecp2Het.CTX.HC.100vol.300counts.pred0.2.sctTable)
#hist(mecp2Het.CTX.HC.100vol.300counts.pred0.2.sctTable[, Mecp2], breaks=100, main="")


dWT_Mecp2 = density(WTandHet.100vol.300counts.pred0.2.sctTable[rep=="Exp3 WT coronal", Mecp2], cut=0, bw=0.5)
dHet_Mecp2 = density(WTandHet.100vol.300counts.pred0.2.sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Mecp2], cut=0, bw=0.5)
dHet_Blank.37 = density(WTandHet.100vol.300counts.pred0.2.sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Blank.37], cut=0, bw=0.5)
#x and y values for plot
d_genes_x = c(dWT_Mecp2$x, dHet_Mecp2$x, dHet_Blank.37$x)  
d_genes_y = c(dWT_Mecp2$y, dHet_Mecp2$y, dHet_Blank.37$y)  

#png("HG_lab/Mati/GabelLab/Vizgen/qc/cell_count_Mecp2_someGenes_Gabel_20220818_densityPlot.png")
plot(dWT_Mecp2, lwd = 1, col = "black",
     main = "Gene count density", xlab = "Gene count per cell",
     xlim = c(min(d_genes_x), c(max(d_genes_x))),  # Min and Max X-axis limits
     ylim = c(min(d_genes_y), c(max(d_genes_y))))  # Min and Max Y-axis limits
lines(dHet_Mecp2, col = "darkred", lwd = 1)
lines(dHet_Blank.37, col = "gray", lwd = 1)
legend("topright", legend=c("WT MeCP2","Mecp2 KO/+ MeCP2","Mecp2 KO/+ Blank-37"), col = c("black","darkred","gray"), lty=1, bty="n") 
#dev.off()

WTandHet.Mecp2.blank.density.func <- function(sctTable, return_table=FALSE, plot_title="", plot_width, plot_height, plot_save=FALSE, plot_file=NULL){
  WTandHet.Mecp2.blank.density.dt <- data.table(rbind(cbind(gene.count=sctTable[rep=="Exp3 WT coronal", Mecp2], id.var="WT Mecp2"),
                                                      cbind(gene.count=sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Mecp2], id.var="Mecp2 KO/+ Mecp2"),
                                                      cbind(gene.count=sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Blank.37], id.var="Mecp2 KO/+ Blank-37")))
  
  WTandHet.Mecp2.blank.density.dt = WTandHet.Mecp2.blank.density.dt %>% mutate(id.var = factor(id.var, levels=c("WT Mecp2", "Mecp2 KO/+ Mecp2", "Mecp2 KO/+ Blank-37")))
  ggplot(WTandHet.Mecp2.blank.density.dt)+
    ggtitle(plot_title)+
    geom_density(aes(x=as.numeric(gene.count), color=id.var))+
    scale_color_manual(name="Gene:", values = c("WT Mecp2"="blue", "Mecp2 KO/+ Mecp2"="red", "Mecp2 KO/+ Blank-37"="brown")) +
    geom_vline(xintercept = 1, color="gray", linetype="longdash")+
    #coord_cartesian(ylim=c(0,2))+
    ylab("Density") + xlab("Sctransform counts")+
    theme_bw()+
    theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  if(plot_save){
    ggsave(filename=paste0(plot_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
    ggsave(filename=paste0(plot_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')
  }
  if(return_table){
    return(WTandHet.Mecp2.blank.density.dt)
  }
}

WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="Pvalb"], 
                                  return_table=FALSE,
                                  plot_title="Pvalb",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/Pvalb.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot")

WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="Sst"], 
                                  return_table=FALSE,
                                  plot_title="Sst",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/Sst.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot")

WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="Astro"], 
                                  return_table=FALSE,
                                  plot_title="Astro",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/Astro.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot")

WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="DG"], 
                                  return_table=FALSE,
                                  plot_title="DG",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/DG.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot")

WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="L5 IT CTX"], 
                                  return_table=FALSE,
                                  plot_title="L5 IT CTX",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/L5.IT.CTX.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot")




WTandHet.Mecp2.blank.density.dt <- data.table(rbind(cbind(gene.count=WTandHet.100vol.300counts.pred0.2.sctTable[rep=="Exp3 WT coronal", Mecp2], id.var="WT Mecp2"),
                                                    cbind(gene.count=WTandHet.100vol.300counts.pred0.2.sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Mecp2], id.var="Mecp2 KO/+ Mecp2"),
                                                    cbind(gene.count=WTandHet.100vol.300counts.pred0.2.sctTable[rep %in% c("Rep1","Rep2","Rep6","Rep7"), Blank.37], id.var="Mecp2 KO/+ Blank-37")))

WTandHet.Mecp2.blank.density.dt = WTandHet.Mecp2.blank.density.dt %>% mutate(id.var = factor(id.var, levels=c("WT Mecp2", "Mecp2 KO/+ Mecp2", "Mecp2 KO/+ Blank-37")))


ggplot(WTandHet.100vol.300counts.pred0.2.sctTable) + geom_density(aes(x=quality, fill=sweetnes))

ggplot(WTandHet.Mecp2.blank.density.dt)+
  geom_density(aes(x=as.numeric(gene.count), color=id.var))+
  scale_color_manual(name="Gene:", values = c("WT Mecp2"="blue", "Mecp2 KO/+ Mecp2"="red", "Mecp2 KO/+ Blank-37"="brown")) +
  geom_vline(xintercept = 1, color="gray", linetype="longdash")+
  #coord_cartesian(ylim=c(ymin,ymax))+
  ylab("Density") + xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


##scatterplots of positions of cells
exp2.CTX.HC.100vol.300counts.pred0.2.sctTable <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore>0.2",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000))
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subclass[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore>0.2",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000))
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subclass[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.png", width=1000, height=1000, res=300)
par(pty="s")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore>0.2",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000))
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subclass[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()




setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Exp2, cortex and hippocampus\n100 <= cell volume <= 3*median(cell volume)\ncounts >= 300, predScore>0.2",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000))
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subclass[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()

#ggplot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable)+
#  ggtitle("")+
#  geom_point(aes(x=center_x, y=center_y, color="gray"))+
#  geom_point(data=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label %in% CTX_HIP_annot[class_label=="GABAergic", subclass_label]], 
#             aes(x=center_x, y=center_y, color=subclass_label))+
#  scale_color_manual(name="Gene:", values = colors_subclass) +
#  coord_cartesian(xlim=c(0,10000), ylim=c(0,10000))+
#  xlab("Exp2, x") + ylab("Exp2, y") +
#  theme_bw()+
#  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))


#points(x=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="WT"), center_x], 
#       y=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="WT"), center_y],
#       pch=16, col="blue", cex=0.4)
#points(x=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="KO"), center_x], 
#       y=exp1.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id=="118_Pvalb") & (t.type.Under1Over1=="KO"), center_y],
#       pch=16, col="red", cex=0.4)


WTandHet.Mecp2.blank.density.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable[subclass_label=="Pvalb"], 
                                  return_table=FALSE,
                                  plot_title="Pvalb",
                                  plot_width=5, plot_height=5, 
                                  plot_save=TRUE, 
                                  plot_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/densityplots/Pvalb.allExps.100vol.300count.pred0.2.Under1Over1.sctCount.densityPlot2")

#mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes
#sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                            #colData = sce@meta.data)
#groups.sce <- colData(sce)[, c("predicted.id", "t.type", "rep")]
#aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)

library(patchwork)

WT.Het.blank37.Pvalb.dt <- data.table(rbind(
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep=="Exp3 WT coronal") & (subclass_label=="Pvalb"), Blank.37], id.var="WT"),
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label=="Pvalb"), Blank.37], id.var="Mecp2 KO/+")))

WT.Het.mecp2.Pvalb.dt <- data.table(rbind(
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep=="Exp3 WT coronal") & (subclass_label=="Pvalb"), Mecp2], id.var="WT"),
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label=="Pvalb"), Mecp2], id.var="Mecp2 KO/+")))


p1 <- ggplot(WT.Het.blank37.Pvalb.dt, aes(x = as.numeric(val), fill = id.var)) + 
  ggtitle("Pvalb")+
  geom_histogram(position = "dodge", binwidth = 1) +
  geom_vline(xintercept=1,col="lightgray")+
  coord_cartesian(xlim=c(0,15))+
  xlab("")+ylab("Blank counts")+
  scale_fill_manual(name="",
                    values = c("WT"="blue",
                               "Mecp2 KO/+"="red"))+
  theme_bw()+
  theme(plot.title =element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))

p2 <- ggplot(WT.Het.mecp2.Pvalb.dt, aes(x = as.numeric(val), fill = id.var)) + 
  ggtitle("")+
  geom_histogram(position = "dodge", binwidth = 1) +
  geom_vline(xintercept=1,col="lightgray")+
  coord_cartesian(xlim=c(0,15))+
  xlab("")+ylab("Mecp2 counts")+
  scale_fill_manual(name="",
                    values = c("WT"="blue",
                               "Mecp2 KO/+"="red"))+
  theme_bw()+
  theme(plot.title =element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))

p1 + p2 + plot_layout(ncol=1)


WT.Het.subclass.blank.mecp2.plot.func <-  function(sctTable, subclass, plot_width=3.5, plot_height=5, output_file){
  WT.Het.blank37.dt <- data.table(rbind(
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Blank.37], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Blank.37], id.var="Mecp2 KO/+")))
  
  WT.Het.mecp2.dt <- data.table(rbind(
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Mecp2], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Mecp2], id.var="Mecp2 KO/+")))
  
  
  p1 <- ggplot(WT.Het.blank37.dt, aes(x = as.numeric(val), fill = id.var)) + 
    ggtitle(subclass)+
    geom_histogram(position = "dodge", binwidth = 1) +
    geom_vline(xintercept=1,col="lightgray", linetype = "longdash")+
    coord_cartesian(xlim=c(0,15))+
    xlab("")+ylab("Blank-37 counts")+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"))+
    theme_bw()+
    theme(plot.title =element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  p2 <- ggplot(WT.Het.mecp2.dt, aes(x = as.numeric(val), fill = id.var)) + 
    ggtitle("")+
    geom_histogram(position = "dodge", binwidth = 1) +
    geom_vline(xintercept=1,col="lightgray", linetype = "longdash")+
    coord_cartesian(xlim=c(0,15))+
    xlab("")+ylab("Mecp2 counts")+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"))+
    theme_bw()+
    theme(plot.title =element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  p1 + p2 + plot_layout(ncol=1)
  ggsave(filename=paste0(output_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
  ggsave(filename=paste0(output_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')
}

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="Pvalb",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Pvalb.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="Sst",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Sst.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="L5 IT CTX",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/L5.IT.CTX.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="Astro",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Astro.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="DG",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/DG.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="L5 PT CTX",
                                      plot_width=2.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/L5.PT.CTX.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.hist")

###



WT.Het.blank37.Pvalb.dt.summary = data.table(WT.Het.blank37.Pvalb.dt %>% 
                                             group_by(id.var, val) %>%
                                             summarise(n.val = n()) %>%
                                             mutate(prop = n.val / sum(n.val)))

WT.Het.mecp2.Pvalb.dt.summary = data.table(WT.Het.mecp2.Pvalb.dt %>% 
                                               group_by(id.var, val) %>%
                                               summarise(n.val = n()) %>%
                                               mutate(prop = n.val / sum(n.val)))

p1 <- ggplot(WT.Het.blank37.Pvalb.dt.summary, aes(x = val, y=as.numeric(prop), fill = id.var)) + 
  ggtitle("Pvalb")+
  geom_col(position="dodge") +
  #coord_cartesian(xlim=c(0,15))+
  xlab("")+ylab("Blank counts")+
  scale_fill_manual(name="",
                    values = c("WT"="blue",
                               "Mecp2 KO/+"="red"))+
  theme_bw()+
  theme(plot.title =element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
####
prop.WT.Het.subclass.blank.mecp2.plot.func1 <-  function(sctTable, subclass, x_max=20, plot_width=3.5, plot_height=5, output_file){
  WT.Het.blank37.dt <- data.table(rbind(
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Blank.37], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Blank.37], id.var="Mecp2 KO/+")))
  
  WT.Het.mecp2.dt <- data.table(rbind(                                               
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Mecp2], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Mecp2], id.var="Mecp2 KO/+")))
  
  
  #num_order <- as.character(seq(0, max(WT.Het.mecp2.dt.summary[,as.numeric(val)], na.rm=TRUE), 1))
  num_order <- as.character(seq(0, x_max, 1))
  
  #blank
  WT.Het.blank37.dt <- WT.Het.blank37.dt %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt <- WT.Het.blank37.dt %>% mutate(val = factor(val, levels=num_order))
  
  
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt %>% 
    dplyr::count(id.var, val, .drop = FALSE) %>%  
    group_by(id.var, val) %>%
    summarise(n.val = n) %>%
    mutate(prop = n.val / sum(n.val)) %>% data.table
  
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  #mecp2
  WT.Het.mecp2.dt <- WT.Het.mecp2.dt %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.mecp2.dt <- WT.Het.mecp2.dt %>% mutate(val = factor(val, levels=num_order))
  
  
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt %>% 
    dplyr::count(id.var, val, .drop = FALSE) %>%  
    group_by(id.var, val) %>%
    summarise(n.val = n) %>%
    mutate(prop = n.val / sum(n.val)) %>% data.table
  
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  
  
  p1 <- ggplot(WT.Het.blank37.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    coord_cartesian(ylim=c(0,1))+
    xlab("")+ylab("Blank count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  
  p2 <- ggplot(WT.Het.mecp2.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    coord_cartesian(ylim=c(0,0.45))+
    xlab("")+ylab("Mecp2 count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title =element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  p1 + p2 + plot_layout(ncol=1)
  ggsave(filename=paste0(output_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
  ggsave(filename=paste0(output_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')
}
###

prop.WT.Het.subclass.blank.mecp2.plot.func <-  function(sctTable, subclass, x_max=20, plot_width=3.5, plot_height=5, output_file){
  WT.Het.blank37.dt <- data.table(rbind(
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Blank.37], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Blank.37], id.var="Mecp2 KO/+")))
  
  WT.Het.mecp2.dt <- data.table(rbind(                                               
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Mecp2], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Mecp2], id.var="Mecp2 KO/+")))
  
  
  WT.Het.blank37.dt.summary = data.table(WT.Het.blank37.dt %>% 
                                                 group_by(id.var, val) %>%
                                                 summarise(n.val = n()) %>%
                                                 mutate(prop = n.val / sum(n.val)))
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(val = factor(val, levels=as.character(unique(sort(WT.Het.blank37.dt.summary[,as.numeric(val)])))))
  
  WT.Het.mecp2.dt.summary = data.table(WT.Het.mecp2.dt %>% 
                                               group_by(id.var, val) %>%
                                               summarise(n.val = n()) %>% mutate(prop = n.val / sum(n.val)))
  
  #num_order <- as.character(seq(0, max(WT.Het.mecp2.dt.summary[,as.numeric(val)], na.rm=TRUE), 1))
  num_order <- as.character(seq(0, x_max, 1))
  
  
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  
  
  p1 <- ggplot(WT.Het.blank37.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    #coord_cartesian(ylim=c(0,1), xlim=c(0,20))+
    xlab("")+ylab("Blank count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  
  p2 <- ggplot(WT.Het.mecp2.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    #coord_cartesian(ylim=c(0,1), xlim=c(0,20))+
    xlab("")+ylab("Mecp2 count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title =element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  p1 + p2 + plot_layout(ncol=1)
  ggsave(filename=paste0(output_file,".png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
  ggsave(filename=paste0(output_file,".eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')
}



prop.WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                      subclass="Pvalb",
                                      x_max=20,
                                      plot_width=3.5, plot_height=5,
                                      output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Pvalb.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist")

prop.WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                           subclass="L5 IT CTX",
                                           x_max=20,
                                           plot_width=3.5, plot_height=5,
                                           output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/L5.IT.CTX.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist")

prop.WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                           subclass="Sst",
                                           x_max=20,
                                           plot_width=3.5, plot_height=5,
                                           output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Sst.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist")

prop.WT.Het.subclass.blank.mecp2.plot.func(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                           subclass="Astro",
                                           x_max=20,
                                           plot_width=3.5, plot_height=5,
                                           output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Astro.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist")


ggplot(WT.Het.blank37.Pvalb.dt.summary, aes(x = val, y=as.numeric(prop), fill = id.var)) + 
  ggtitle("Pvalb")+
  geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
  coord_cartesian(ylim=c(0,15))+
  xlab("")+ylab("Blank counts")+
  scale_fill_manual(name="",
                    values = c("WT"="blue",
                               "Mecp2 KO/+"="red"))+
  theme_bw()+
  theme(plot.title =element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))


###
WT.Het.blank37.Pvalb.dt <- data.table(rbind(
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep=="Exp3 WT coronal") & (subclass_label=="Pvalb"), Blank.37], id.var="WT"),
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label=="Pvalb"), Blank.37], id.var="Mecp2 KO/+")))

WT.Het.mecp2.Pvalb.dt <- data.table(rbind(
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep=="Exp3 WT coronal") & (subclass_label=="Pvalb"), Mecp2], id.var="WT"),
  cbind(val=WTandHet.100vol.300counts.pred0.2.sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label=="Pvalb"), Mecp2], id.var="Mecp2 KO/+")))


WT.Het.blank37.Pvalb.dt.summary = data.table(WT.Het.blank37.Pvalb.dt %>% 
                                         group_by(id.var, val) %>%
                                         summarise(n.val = n()) %>%
                                         mutate(prop = n.val / sum(n.val)))
num_order = seq(0, 20, 1)

WT.Het.blank37.Pvalb.dt <- WT.Het.blank37.Pvalb.dt %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
WT.Het.blank37.Pvalb.dt <- WT.Het.blank37.Pvalb.dt %>% mutate(val = factor(val, levels=num_order))


WT.Het.blank37.Pvalb.dt.summary1 <- WT.Het.blank37.Pvalb.dt %>% 
  dplyr::count(id.var, val, .drop = FALSE) %>%  
  group_by(id.var, val) %>%
  summarise(n.val = n) %>%
  mutate(prop = n.val / sum(n.val)) %>% data.table


WT.Het.blank37.Pvalb.dt.summary <- WT.Het.blank37.Pvalb.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
WT.Het.blank37.Pvalb.dt.summary <- WT.Het.blank37.Pvalb.dt.summary %>% mutate(val = factor(val, levels=num_order))


WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(val = factor(val, levels=num_order))



prop.WT.Het.subclass.blank.mecp2.plot.func1(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                           subclass="Astro",
                                           x_max=20,
                                           plot_width=3.5, plot_height=5,
                                           output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Astro.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist2")

prop.WT.Het.subclass.blank.mecp2.plot.func1(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                            subclass="Pvalb",
                                            x_max=20,
                                            plot_width=3.5, plot_height=5,
                                            output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Pvalb.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist2")


#make color key
head(CTX_HIP_annot)

ggplot(data=CTX_HIP_annot, aes(x=cluster_label, y = 0.5, fill=cluster_label))+
  geom_tile() +
  scale_fill_manual(values = colors_subtype) +
  theme_void()+
  theme(legend.position="none") 

ggplot(data=CTX_HIP_annot, aes(x=cluster_label, y = 1, fill=cluster_label))+
  geom_tile() +
  scale_fill_manual(values = colors_subtype) +
  facet_grid(.~subclass_label,
             switch = "x",
             labeller = label_wrap_gen(width=8),
             scales = "free_x"
  ) +
  theme(legend.position="none") 



used_clusters <- unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes$predicted.id)
CTX_HIP_annot_used <- CTX_HIP_annot[cluster_label %in% used_clusters]

#ordering purposes
subclass_order <- unique(CTX_HIP_annot_used$subclass_label)
subclass_order_dt <- data.table(subclass_order)
names(subclass_order_dt) = "subclass_label"
subclass_order_dt <- left_join(x=subclass_order_dt, y=CTX_HIP_annot_used[, .(subclass_label, subclass_color, cluster_label, cluster_color)], by=("subclass_label"))

cluster_order <- unique(subclass_order_dt$cluster_label)

CTX_HIP_annot_ordered <- subclass_order_dt  %>% mutate(cluster_label = factor(cluster_label, levels=cluster_order))
CTX_HIP_annot_ordered <- CTX_HIP_annot_ordered %>% mutate(subclass_label = factor(subclass_label, levels=subclass_order))

ggplot(data=CTX_HIP_annot_ordered, aes(x=1, y = cluster_label, fill = cluster_label))+
  geom_tile() +
  scale_fill_manual(values = colors_subtype[used_clusters]) +
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.png", width = 1, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.eps", width = 1, height = 5, dpi = 300, units = "in", device='eps')

#splitting heatmap

ggplot(data=CTX_HIP_annot_ordered[1:135], aes(x=1, y = cluster_label, fill = cluster_label))+
  geom_tile() +
  scale_fill_manual(values = colors_subtype[used_clusters]) +
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.png", width = 1.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.eps", width = 1.5, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=CTX_HIP_annot_ordered[136:273], aes(x=1, y = cluster_label, fill = cluster_label))+
  geom_tile() +
  scale_fill_manual(values = colors_subtype[used_clusters]) +
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.SecondHalf.png", width = 1.5, height = 5.074075, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.SecondHalf.eps", width = 1.5, height = 5.074075, dpi = 300, units = "in", device='eps')


subclass_order_dt[, final_label := paste0(subclass_label,", ", cluster_label)]
final_order <- unique(subclass_order_dt$final_label)

name_vec_final = unique(subclass_order_dt$final_label)
col_vec_final = unique(subclass_order_dt$cluster_color)
colors_final = setNames(col_vec_final, name_vec_final)


CTX_HIP_annot_final <- subclass_order_dt  %>% mutate(final_label = factor(final_label, levels=final_order))

ggplot(data=CTX_HIP_annot_final[1:135], aes(x=1, y = final_label, fill = final_label))+
  geom_tile() +
  scale_fill_manual(values = colors_final) +
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.png", width = 1.5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.eps", width = 1.5, height = 5, dpi = 300, units = "in", device='eps')


ggplot(data=CTX_HIP_annot_final[136:273], aes(x=1, y = final_label, fill = final_label))+
  geom_tile() +
  scale_fill_manual(values = colors_final) +
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.SecondHalf.png", width = 1.5, height = 5.074075, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.SecondHalf.eps", width = 1.5, height = 5.074075, dpi = 300, units = "in", device='eps')

###remaking scatterplots of locations of cells
#left visual cortex and other areas cells 
exp2Liu_L_VIS <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS.csv")
names(exp2Liu_L_VIS) <- "index"
#right
exp2Liu_R_VIS <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS.csv")
names(exp2Liu_R_VIS) <- "index"

exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS <- exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[index %in% exp2Liu_L_VIS$index]
exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.R.VIS <- exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[index %in% exp2Liu_R_VIS$index]


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_y], main="Glutamatergic",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.15, xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.25, col=colors_subclass[i], xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_y], main="GABAergic",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.15, xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.25, col=colors_subclass[i], xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[, center_y], main="Non-neuronal",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.15, xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.L.VIS[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.25, col=colors_subclass[i], xlim=c(2000,5000), ylim=c(7000,10000), asp=1)
}
dev.off()



plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.R.VIS[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable.R.VIS[, center_y], main="Glutamatergic",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.15, xlim=c(0,10000), ylim=c(0,10000), asp=1)


###
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.WT.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.KO.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.WT.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.KO.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.WT.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.KO.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

###other side
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAaergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subclassColors.scatterplot.png")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,2500), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1000,2500), ylim=c(0,6500), asp=1)
}
dev.off()

###
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2 <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes
DefaultAssay(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2) <- "RNA"
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2.allReps.cpmNorm = NormalizeData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2, normalization.method="RC", scale.factor = 1e6)
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2.allReps.cpmNorm.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2.allReps.cpmNorm, subset=t.type=="WT")
sce.cellTypes2.allReps.WT <- SingleCellExperiment(assays = list(counts = GetAssayData(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2.allReps.cpmNorm.WT[["RNA"]], slot = "data")), 
                                       colData = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes2.allReps.cpmNorm.WT@meta.data)
groups.cellTypes2.allReps.WT <- colData(sce.cellTypes2.allReps.WT)[, c("predicted.id")]
agg.cellTypes2.allReps.WT.avgCPM <- aggregateAcrossCells(sce.cellTypes2.allReps.WT , ids=groups.cellTypes2.allReps.WT , statistics="mean")
agg.exp7.WT.avgCPM.dt = data.table(assay(agg.exp7.WT.avgCPM), keep.rownames="Gene")
####


prop.WT.Het.subclass.blank.mecp2.plot.func2 <-  function(sctTable, subclass, x_max=14, plot_width=5, plot_height=5, output_file){
  WT.Het.blank37.dt <- data.table(rbind(
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Blank.37], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Blank.37], id.var="Mecp2 KO/+")))
  
  WT.Het.mecp2.dt <- data.table(rbind(                                               
    cbind(val=sctTable[(rep=="Exp3 WT coronal") & (subclass_label==subclass), Mecp2], id.var="WT"),
    cbind(val=sctTable[(rep %in% c("Rep1","Rep2","Rep6","Rep7")) & (subclass_label==subclass), Mecp2], id.var="Mecp2 KO/+")))
  
  
  #num_order <- as.character(seq(0, max(WT.Het.mecp2.dt.summary[,as.numeric(val)], na.rm=TRUE), 1))
  num_order <- as.character(seq(0, x_max, 1))
  
  #blank
  WT.Het.blank37.dt <- WT.Het.blank37.dt %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt <- WT.Het.blank37.dt %>% mutate(val = factor(val, levels=num_order))
  
  
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt %>% 
    dplyr::count(id.var, val, .drop = FALSE) %>%  
    group_by(id.var, val) %>%
    summarise(n.val = n) %>%
    mutate(prop = n.val / sum(n.val)) %>% data.table
  
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.blank37.dt.summary <- WT.Het.blank37.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  #mecp2
  WT.Het.mecp2.dt <- WT.Het.mecp2.dt %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.mecp2.dt <- WT.Het.mecp2.dt %>% mutate(val = factor(val, levels=num_order))
  
  
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt %>% 
    dplyr::count(id.var, val, .drop = FALSE) %>%  
    group_by(id.var, val) %>%
    summarise(n.val = n) %>%
    mutate(prop = n.val / sum(n.val)) %>% data.table
  
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(id.var = factor(id.var, levels=c("WT", "Mecp2 KO/+")))
  WT.Het.mecp2.dt.summary <- WT.Het.mecp2.dt.summary %>% mutate(val = factor(val, levels=num_order))
  
  
  
  ggplot(WT.Het.blank37.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,1))+
    xlab("")+ylab("Blank count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme(aspect.ratio = 1)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  ggsave(filename=paste0(output_file,"_blank.png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
  ggsave(filename=paste0(output_file,"_blank.eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')

  ggplot(WT.Het.mecp2.dt.summary, aes(x = factor(val), y=as.numeric(prop), fill = id.var)) + 
    ggtitle(subclass)+
    geom_bar(position=position_dodge(preserve = "single"), stat="identity") +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,0.5))+
    xlab("")+ylab("Mecp2 count proportion")+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(name="",
                      values = c("WT"="blue",
                                 "Mecp2 KO/+"="red"), drop=FALSE)+
    theme(aspect.ratio = 1)+
    theme_bw()+
    #theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8))+
    theme(plot.title =element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x=element_text(size=10))
  
  ggsave(filename=paste0(output_file,"_mecp2.png"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='png')
  ggsave(filename=paste0(output_file,"_mecp2.eps"), width = plot_width, height = plot_height, dpi = 300, units = "in", device='eps')
  
}

prop.WT.Het.subclass.blank.mecp2.plot.func2(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                            subclass="Astro",
                                            x_max=15,
                                            plot_width=4.9, plot_height=5,
                                            output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Astro.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist3")

prop.WT.Het.subclass.blank.mecp2.plot.func2(sctTable=WTandHet.100vol.300counts.pred0.2.sctTable, 
                                            subclass="Pvalb",
                                            x_max=15,
                                            plot_width=4.9, plot_height=5,
                                            output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/histograms/Pvalb.100vol.300count.pred0.2.Under1Over1.sctransform.blank37.mecp2.counts.prop.hist3")

###counting cell number through stacked bar plot

#summarize counts
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary = group_by(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt, predicted.id, subclass_label) %>%
  summarise(counts.total = n()) %>% data.table

#subclass order
subclass_order = unique(CTX_HIP_annot[subclass_label %in% mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary$subclass_label,subclass_label])
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary = mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary %>% mutate(subclass_label = factor(subclass_label, levels=subclass_order))


ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary, aes(fill=predicted.id, y=counts.total, x=subclass_label)) + 
  geom_bar(position='stack', stat="identity")+
  scale_fill_manual(values = colors_subtype[levels(factor(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[, predicted.id]))]) +
  ylab("Number of cells")+xlab("Subclass")+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=6, angle=90))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.png", width = 10, height = 6, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.eps", width = 10, height = 6, dpi = 300, units = "in", device='eps')

ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[subclass_label %in% subclass_order[1:21]], aes(fill=predicted.id, y=counts.total, x=subclass_label)) + 
  geom_bar(position='stack', stat="identity")+
  scale_fill_manual(values = colors_subtype[levels(factor(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[, predicted.id]))]) +
  ylab("Number of cells")+xlab("Subclass")+
  coord_cartesian(ylim=c(0,15000))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=6, angle=90))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.FirstHalf.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.FirstHalf.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')



ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[subclass_label %in% subclass_order[22:42]], aes(fill=predicted.id, y=counts.total, x=subclass_label)) + 
  geom_bar(position='stack', stat="identity")+
  scale_fill_manual(values = colors_subtype[levels(factor(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[, predicted.id]))]) +
  ylab("Number of cells")+xlab("Subclass")+
  coord_cartesian(ylim=c(0,15000))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=6, angle=90))
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.SecondHalf.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.SecondHalf.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

##
ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[subclass_label %in% subclass_order[1:21]], aes(fill=predicted.id, y=counts.total, x=subclass_label)) + 
  geom_bar(position='stack', stat="identity")+
  scale_fill_manual(values = colors_subtype[levels(factor(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[, predicted.id]))]) +
  ylab("Number of cells")+xlab("")+
  coord_cartesian(ylim=c(0,15000))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.FirstHalf.noLabels.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.FirstHalf.noLabels.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')



ggplot(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[subclass_label %in% subclass_order[22:42]], aes(fill=predicted.id, y=counts.total, x=subclass_label)) + 
  geom_bar(position='stack', stat="identity")+
  scale_fill_manual(values = colors_subtype[levels(factor(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.dt.summary[, predicted.id]))]) +
  ylab("Number of cells")+xlab("")+
  coord_cartesian(ylim=c(0,15000))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.x=element_blank())
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.SecondHalf.noLabels.png", width = 5, height = 2.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/barplots/mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.clusterCountsTotal.stackedBarplot.SecondHalf.noLabels.eps", width = 5, height = 2.5, dpi = 300, units = "in", device='eps')

##
nrow(WTandHet.100vol.300counts.pred0.2.sctTable[(rep %in% c("Rep1", "Rep2", "Rep6", "Rep7")) & (Blank.37 < 2)])

##
##scatterplots of positions of cells, subtype colors
exp2.CTX.HC.100vol.300counts.pred0.2.sctTable <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")

plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Subclass",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000))
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[subclass_label==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subclass[i], xlim=c(0,10000), ylim=c(0,10000))
}

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[predicted.id==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subtype[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[predicted.id==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subtype[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.2, xlim=c(0,10000), ylim=c(0,10000), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[predicted.id==i, .(center_x, center_y)], pch=19, cex=0.4, col=colors_subtype[i], xlim=c(0,10000), ylim=c(0,10000))
}
dev.off()

###
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Liu2023.L.VIS.WT.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subclassColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", subclass_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(subclass_label==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subclass[i], xlim=c(1500,5500), ylim=c(6500,9500), asp=1)
}
dev.off()

#right non-neuronal WT
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()
#right non-neuronal KO
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.nonneuronal.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Non-Neuronal, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Non-Neuronal", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

#right glutamatergic WT
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}

#right glutamatergic KO
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.glutamatergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="Glutamatergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="Glutamatergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()

#right gabagergic WT
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.WT.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, WT",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="WT"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()
#right gabaergic KO
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/zoomed_in/right.KO.exp2.CTX.HC.100vol.300counts.pred0.2.gabaergic.cells.subtypeColors.scatterplot.eps")
plot(exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_x], exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[, center_y], main="GABAergic, KO",
     xlab="Exp2, x", ylab="Exp2, y", pch=16, col="gray", cex=0.25, xlim=c(1000,4000), ylim=c(0,6500), asp=1)
for(i in unique(CTX_HIP_annot[class_label=="GABAergic", cluster_label])){
  points(x=exp2.CTX.HC.100vol.300counts.pred0.2.sctTable[(predicted.id==i) & (t.type.Under1Over1=="KO"), .(center_x, center_y)], pch=19, cex=0.35, col=colors_subtype[i], xlim=c(1000,4000), ylim=c(0,6500), asp=1)
}
dev.off()
