library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
library(Seurat)
library(SingleCellExperiment)
library(reshape2)
library(scuttle)
library(scran)
library(tibble)
options(scipy=999)

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
#Seurat object 
mecp2Het.CTX.HC.seu=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_allTaxLabels_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

#file with layer depth values annotated for each cell in V1 L2/3
all_VIS_L23_depths=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_VIS_L23_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")

all_VIS_L23_depths[depth_name=="Q1", index]
#dividing Seurat objects into quintiles of VIS L2/3
#mecp2Het.CTX.HC.seu.VIS.L23.Q1 <- subset(mecp2Het.CTX.HC.seu, cells=)


WT.KO.type.de.noshrink.func <- function(seurat_obj){
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

de_results_func <- function(de.results, idcol){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = idcol))
  names(de.results.dt)[2] = "gene"
  de.results.dt$p.adj.BH.all = p.adjust(de.results.dt$PValue, method="BH")
  return(de.results.dt)
}

###from Rinaldo
all_V1_L23_IT_DEG <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/All_V1.L23IT_layerdepth.5quantiles.noNA.de.results.dt.csv")

#types_of_interest = c("166_L2/3 IT CTX", "171_L2/3 IT CTX", "168_L2/3 IT CTX", "179_L4 IT CTX")
#L23_V1_logFC_TOI <- L23_V1_logFC[subtype %in% types_of_interest]

#L23_V1_logFC_TOI = L23_V1_logFC_TOI %>% mutate(subtype = factor(subtype, levels=types_of_interest))


example_genes_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
#example genes using L2/3 V1 pseudobulkDGE, nCells=10
all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper <- dcast(all_V1_L23_IT_DEG[gene %in% c("Cdh13", "Epha6", "Tox", "Bmper")], layer_quantile ~ gene, value.var = "logFC") %>% data.table
all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper = all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper[, .(layer_quantile, Cdh13, Epha6, Tox, Bmper)]
all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper_matrix = t(as.matrix(all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper , rownames="layer_quantile"))


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells10_L23_V1_L23IT_179L4IT_CTX_Bmper_Cdh13_Epha6_Tox_logFC_merfish_heatmap.eps")
heatmap.2(all_V1_L23_IT_DEG_Cdh13_Epha6_Tox_Bmper_matrix,
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, col=example_genes_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

all_V1_L23_IT_DEG_Cdh13_Tox <- all_V1_L23_IT_DEG[gene %in% c("Cdh13", "Tox")]
all_V1_L23_IT_DEG_Cdh13_Tox$sig = sapply(all_V1_L23_IT_DEG_Cdh13_Tox $p.adj.BH.all, sig_function)

# Define the order of genes and layer_quantile
gene_order <- c("Tox", "Cdh13")
layer_order <- c(1, 2, 3, 4, 5)

# Create factor columns
all_V1_L23_IT_DEG_Cdh13_Tox$gene_factor <- factor(all_V1_L23_IT_DEG_Cdh13_Tox$gene, levels = gene_order)
all_V1_L23_IT_DEG_Cdh13_Tox$layer_factor <- factor(all_V1_L23_IT_DEG_Cdh13_Tox$layer_quantile, levels = layer_order)

# Set the key for sorting
setkey(all_V1_L23_IT_DEG_Cdh13_Tox, gene_factor, layer_factor)

# Reorder the data.table
all_V1_L23_IT_DEG_Cdh13_Tox_reordered <- all_V1_L23_IT_DEG_Cdh13_Tox[]


##Bmper and Epha6
all_V1_L23_IT_DEG_Epha6_Bmper <- all_V1_L23_IT_DEG[gene %in% c("Epha6", "Bmper")]
all_V1_L23_IT_DEG_Epha6_Bmper$sig = sapply(all_V1_L23_IT_DEG_Epha6_Bmper $p.adj.BH.all, sig_function)

# Define the order of genes and layer_quantile
gene_order_Bmper_Epha6 <- c("Bmper", "Epha6")
layer_order <- c(1, 2, 3, 4, 5)

# Create factor columns
all_V1_L23_IT_DEG_Epha6_Bmper$gene_factor <- factor(all_V1_L23_IT_DEG_Epha6_Bmper$gene, levels = gene_order_Bmper_Epha6)
all_V1_L23_IT_DEG_Epha6_Bmper$layer_factor <- factor(all_V1_L23_IT_DEG_Epha6_Bmper$layer_quantile, levels = layer_order)

# Set the key for sorting
setkey(all_V1_L23_IT_DEG_Epha6_Bmper, gene_factor, layer_factor)

# Reorder the data.table
all_V1_L23_IT_DEG_Epha6_Bmper_reordered <- all_V1_L23_IT_DEG_Epha6_Bmper[]

###
WT.KO.layer.de.noshrink.func <- function(seurat_obj){
  sce <- subset(seurat_obj, subset = ((t.type == "WT") | (t.type == "KO")))
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                              colData = sce@meta.data)
  groups.sce <- colData(sce)[, c("layer_quantile", "t.type", "rep")]
  aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)
  aggr_counts.filt <- aggr_counts[,aggr_counts$ncells >= 10]
  col.data = colData(aggr_counts.filt)
  col.data$t.type = factor(col.data$t.type, levels=c("WT", "KO"))
  de.results.noShrink <- pseudoBulkDGE(aggr_counts.filt, 
                                       label=aggr_counts.filt$layer_quantile,
                                       col.data=col.data,
                                       design=~ rep + t.type,
                                       coef="t.typeKO",
                                       condition=aggr_counts.filt$t.type,
                                       robust=FALSE)
  return(de.results.noShrink)
}
#processing DEG result file to a data table
layer_de_results_func <- function(de.results){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = "layer_quantile"))
  names(de.results.dt)[2] = "gene"
  de.results.dt$p.adj.BH.all = p.adjust(de.results.dt$PValue, method="BH")
  return(de.results.dt)
}



###read in layer depth information
All_V1.L23_layerdepth <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")

All_V1.L23_layerdepth_glut <- All_V1.L23_layerdepth[class_label=="Glutamatergic"]


colnames(All_V1.L23_layerdepth_glut)[1] <- "cell.id"

#annotate 5 quantiles of layer depth
All_V1.L23_layerdepth_glut$layer_quantile <- All_V1.L23_layerdepth_glut$layer_quintile

#subset Seurat object to only superficial L2/3 IT cells
mecp2Het.CTX.HC.seu.V1.L23.glut <- subset(mecp2Het.CTX.HC.seu, cells=All_V1.L23_layerdepth_glut$cell.id)

#turn relevant layer metadata into data frames with cell ids as row names
All_V1.L23_layerdepth_glut_meta_df <- data.frame(All_V1.L23_layerdepth_glut[, .(cell.id, norm_depth, layer_quantile)], row.names="cell.id")

#add layer metadata
mecp2Het.CTX.HC.seu.V1.L23.glut = AddMetaData(object=mecp2Het.CTX.HC.seu.V1.L23.glut, metadata=All_V1.L23_layerdepth_glut_meta_df)


#V1 L2/3 glutamatergic pseudobulkDGE
All_V1.L23_layerdepth_glut.de.results <- WT.KO.layer.de.noshrink.func(mecp2Het.CTX.HC.seu.V1.L23.glut)
All_V1.L23_layerdepth_glut.de.results.dt <- layer_de_results_func(All_V1.L23_layerdepth_glut.de.results)

#write out files
write.csv(All_V1.L23_layerdepth_glut.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/pseudobulkDGE.glut.cells.in.V1.L23.layerdepth.5quantiles.noNA.de.results.csv", quote=F, row.names=F)

###
#example genes using L2/3 V1 pseudobulkDGE, nCells=10
all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper <- dcast(All_V1.L23_layerdepth_glut.de.results.dt[gene %in% c("Cdh13", "Epha6", "Tox", "Bmper")], layer_quantile ~ gene, value.var = "logFC") %>% data.table
all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper = all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper[, .(layer_quantile, Cdh13, Epha6, Tox, Bmper)]
all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper_matrix = t(as.matrix(all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper , rownames="layer_quantile"))


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells10_glut_cells_in_V1_L23_Bmper_Cdh13_Epha6_Tox_logFC_merfish_heatmap.eps")
heatmap.2(all_V1_L23_glut_DEG_Cdh13_Epha6_Tox_Bmper_matrix,
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, col=example_genes_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper <- All_V1.L23_layerdepth_glut.de.results.dt[gene %in% c("Cdh13", "Epha6", "Tox", "Bmper")]
all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$sig = sapply(all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$p.adj.BH.all, sig_function)

# Define the order of genes and layer_quantile
gene_order_Cdh13_Tox_Epha6_Bmper <- c("Cdh13", "Tox", "Epha6", "Bmper")
layer_order <- c(1, 2, 3, 4, 5)

# Create factor columns
all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$gene_factor <- factor(all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$gene, levels = gene_order_Cdh13_Tox_Epha6_Bmper)
all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$layer_factor <- factor(all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper$layer_quantile, levels = layer_order)

# Set the key for sorting
setkey(all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper, gene_factor, layer_factor)

# Reorder the data.table
all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper_reordered <- all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper[]

write.csv(all_V1_L23_glut_DEG_dt_Cdh13_Tox_Epha6_Bmper_reordered, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/sig_stats_nCells10_glut_cells_in_V1_L23_Bmper_Cdh13_Epha6_Tox_logFC_merfish_heatmap.csv", row.names=F)
