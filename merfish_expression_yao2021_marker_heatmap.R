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
#Allen19_cluster_exp_medians = fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Allen19_cluster_exp_medians.csv")
#annotation
CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
CTX_HIP_annot[, cl := as.character(cl)]

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

#expression table
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg=fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
#MGE markers
MGE_markers = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/marker_genes/MGE_markers.txt", header=FALSE)$V1
MGE_markers_panel <- intersect(MGE_markers, gene_panel_classes[,Gene])

#CGE markers
CGE_markers = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/marker_genes/CGE_markers.txt", header=FALSE)$V1
CGE_markers_panel <- intersect(CGE_markers, gene_panel_classes[,Gene])


MGE_CGE_subtypes <- intersect(unique(CTX_HIP_annot[(neighborhood_label=="CGE") | (neighborhood_label=="MGE"), cluster_label]), colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg))

CGE_MGE_markers_panel <- c(CGE_markers_panel, MGE_markers_panel)

head(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg[CGE_MGE_markers_panel, MGE_CGE_subtypes])

CGE_MGE_palette <- colorRampPalette(c("white", "goldenrod"))(n = 299)

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/averageExpression.sctransform.CGE.MGE.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.png", width=2000, height=4000, res=300)
heatmap.2(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg[CGE_MGE_markers_panel, MGE_CGE_subtypes]),
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=CGE_MGE_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, keysize=1)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/averageExpression.sctransform.CGE.MGE.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=8, height=20)
heatmap.2(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg[CGE_MGE_markers_panel, MGE_CGE_subtypes]),
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=CGE_MGE_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, keysize=1)
dev.off()

###exc
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT, assays="SCT", slot="data", group.by = c("subclass_label"))$SCT
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT


mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.melt = data.table(melt(data.table(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg), keep.rownames="subclass_label"), id.vars="subclass_label"))
names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.melt) = c("subclass_label", "gene", "avgExp")

exc_markers <- fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/marker_genes/glutamatergic_subclass_markers.txt", header=FALSE)$V1
exc_markers_panel <- intersect(exc_markers, gene_panel_classes[,Gene])

exc_subclass_list <- c("CR", "L2/3 IT PPP", "L3 IT ENT", "L2 IT ENTm", "L2 IT ENTl", "L2/3 IT ENTl", "L2/3 IT CTX",
                       "L2/3 IT RHP", "L4/5 IT CTX", "L5 IT CTX", "L5/6 IT TPE−ENT", "L6 IT CTX",
                       "L6 IT ENTl",
                       "Car3",
                       "L5 PT CTX",
                       "L4 RSP−ACA",
                       "L5 PPP",
                       "L5/6 NP CTX",
                       "NP SUB",
                       "NP PPP",
                       "L6 CT CTX",
                       "CT SUB",
                       "L6b/CT ENT",
                       "L6b CTX",
                       "SUB−ProS",
                       "CA1−ProS",
                       "CA3",
                       "CA2-IG−FC",
                       "DG")
exc_subclass_list <- intersect(exc_subclass_list, colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg))

exc_palette <- colorRampPalette(c("white", "goldenrod"))(n = 299)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.exc <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg[exc_markers_panel, exc_subclass_list]



png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/averageExpression.sctransform.glutamatergic.subclass.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.png", width=2000, height=4000, res=300)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg[exc_markers_panel, exc_subclass_list],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=exc_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, keysize=1)
dev.off()

#z-score
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.exc.scale <- t(scale(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.exc)))
exc_palette2 <- colorRampPalette(c("blue", "white", "red"))(n = 199)
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/zscore.averageExpression.sctransform.glutamatergic.subclass.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps")
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.exc.scale,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=exc_palette2,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, keysize=1, key.title = "Z-score")
dev.off()

inh_subclass_list <- unique(CTX_HIP_annot[(class_label=="GABAergic") & (neighborhood_label!="Other"), subclass_label])

#GABAaergic subclasses, except Meis2
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.inh <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg[CGE_MGE_markers_panel, inh_subclass_list]

#z-score
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.inh.scale <- t(scale(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.inh)))

inh_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
#inh_breaks = c(seq(-2.1,-0.5,length=100),  # for blue,
#               seq(-0.49,0.49, length=100), #for white
#               seq(0.5,2.1,length=100))

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=5, height=10)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subclass.avg.inh.scale,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=inh_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=inh_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, keysize=1, key.title = "Z-score", symkey=T)
dev.off()


#ggplot(data=CTX_HIP_annot_ordered[1:135], aes(x=1, y = cluster_label, fill = cluster_label))+
#  geom_tile() +
#  scale_fill_manual(values = colors_subtype[used_clusters]) +
#  theme_bw()+
#  theme(legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        axis.text.y=element_text(size=3), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title=element_blank())
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.png", width = 1.5, height = 5, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/umaps/mecp2Het.CTX.HC.exp1.exp2.exp6.exp7.100vol.300counts.predScore0.2.colorKey.FirstHalf.eps", width = 1.5, height = 5, dpi = 300, units = "in", device='eps')

##

#ordering purposes

used_clusters_WT <- unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT$predicted.id)
CTX_HIP_annot_used_WT <- CTX_HIP_annot[cluster_label %in% used_clusters_WT]

subclass_order_WT <- unique(CTX_HIP_annot_used_WT$subclass_label)
subclass_order_WT_dt <- data.table(subclass_order_WT)
names(subclass_order_WT_dt) = "subclass_label"
subclass_order_WT_dt <- left_join(x=subclass_order_WT_dt, y=CTX_HIP_annot_used_WT[, .(subclass_label, subclass_color, cluster_label, cluster_color)], by=("subclass_label"))

cluster_order_WT <- unique(subclass_order_WT_dt$cluster_label)

CTX_HIP_annot_ordered_WT <- subclass_order_WT_dt  %>% mutate(cluster_label = factor(cluster_label, levels=cluster_order_WT))
CTX_HIP_annot_ordered_WT <- CTX_HIP_annot_ordered_WT %>% mutate(subclass_label = factor(subclass_label, levels=subclass_order_WT))


setdiff(cluster_order_WT, colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg))
#z-score
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale <- scale(t(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg))

setdiff(gene_panel_classes[,Gene], colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale))
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490]

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.avg.scale[cluster_order]

all_palette <- colorRampPalette(c("blue", "white", "red"))(n = 199)



png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.png", width=10, height=5, units="in", res=300)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a column dendrogram
          Rowv=NA,
          Colv=NA,key=TRUE,
          cexCol=0.1)
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=5)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a column dendrogram
          Rowv=NA,
          Colv=NA,key=TRUE,
          cexCol=0.1,
          cexRow=0.1)
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/subclassColors.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=5)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="none",     # only draw a column dendrogram
          Rowv=NA,
          Colv=NA,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$subclass_color)
dev.off()


#clustered
setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/colClust.subclassColors.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=5)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$subclass_color)
dev.off()

out <- heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
                 main = "", # heat map title
                 notecol="black",      # change font color of cell labels to black
                 #breaks=col_breaks2,    # enable color transition at specified limits
                 density.info="none",  # turns off density plot inside color legend
                 trace="none",         # turns off trace lines inside the heat map
                 #margins =c(10,10),     # widens margins around plot
                 col=all_palette,       # use on color palette defined earlier
                 #col=brewer.pal(9,"YlOrRd"),
                 dendrogram="column",     # only draw a column dendrogram
                 Rowv=NA,
                 Colv=TRUE)
View(cbind(colnames(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490])[out$colInd], "hi"))

#ggplot(test_pval_and_fc,aes(variable,cellType)) + geom_tile(aes(fill=pvalMax),color = "white") +
#  #Creating legend
#  guides(fill=guide_colorbar(title="-log10 p-value", title.position="top")) +
#  #Creating color range
#  scale_fill_gradientn(limits = c(0,20), colors=c("white",color),guide="colorbar") +
#  geom_text(aes(label = fc), color = "black", size = 6)+
#  theme_bw()+
  ##Rotating labels
  #theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
  #      axis.text=element_text(size=18), axis.ticks = element_blank(), axis.title=element_blank())

all_palette2 <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(n = 499)
all_breaks = c(seq(-4,-2.26,length=100),  # for darkblue,
               seq(-2.25,-0.51, length=100),
               seq(-0.5,0.5, length=100), #for white
               seq(0.51, 2.25, length=100), #for red,
               seq(2.26, 4, length=100)) #for darkred

all_breaks2 = c(seq(-2,-1.26,length=100),  # for darkblue,
               seq(-1.25,-0.51, length=100),
               seq(-0.5,0.5, length=100), #for white
               seq(0.51, 1.25, length=100), #for red,
               seq(1.26, 2, length=100)) #for darkred


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/subclassColors.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap2.eps", width=15, height=5)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$subclass_color)
dev.off()


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/subclassColors.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap4.eps", width=15, height=8)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          #col=brewer.pal(9,"YlOrRd"),
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$subclass_color)
dev.off()

CTX_HIP_annot_ordered_WT <- left_join(x=CTX_HIP_annot_ordered_WT, y=CTX_HIP_annot[, .(subclass_label, subclass_color, cluster_label, cluster_color, neighborhood_label, neighborhood_color, cluster_id, subclass_id, neighborhood_id)])


CTX_HIP_annot_ordered_WT[order(subclass_id), cluster_label]
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$cluster_color)

heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[cluster_order_WT, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT$neighborhood_color)

###order
CTX_HIP_annot_ordered_WT_2 <- CTX_HIP_annot_ordered_WT[order(subclass_id)]


setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/subtypeColors.ordered.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=8)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[CTX_HIP_annot_ordered_WT_2$cluster_label, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT_2$cluster_color)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/nhoodColors.ordered.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=8)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[CTX_HIP_annot_ordered_WT_2$cluster_label, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT_2$neighborhood_color)
dev.off()


heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[CTX_HIP_annot_ordered_WT_2$cluster_label, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT_2$cluster_color)

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/subclassColors.ordered.zscore.averageExpression.sctransform.all.subtypes.markers.mecp2.CTX.HC.100vol.300counts.pred0.2.Under1Over1.WT.heatmap.eps", width=15, height=8)
heatmap.2(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.subtype.avg.scale[CTX_HIP_annot_ordered_WT_2$cluster_label, 1:490],
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          breaks=all_breaks2,    # enable color transition at specified limits
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(10,10),     # widens margins around plot
          col=all_palette2,       # use on color palette defined earlier
          dendrogram="column",     # only draw a column dendrogram
          Rowv=NA,
          Colv=TRUE,key=TRUE,
          cexCol=0.1,
          cexRow=0.1,
          RowSideColors = CTX_HIP_annot_ordered_WT_2$subclass_color, key.title = "Z-score")
dev.off()