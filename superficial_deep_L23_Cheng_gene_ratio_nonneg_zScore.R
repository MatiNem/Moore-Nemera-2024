library(data.table)
library(dplyr)
library(ggplot2)

#annotation table from Yao et al, 2021
CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
#makes the cl column a character, useful for some data table joins later
CTX_HIP_annot[, cl := as.character(cl)]

##
#summary cell data for each experiment passing filters
exp1_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp2_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp6_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp6.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp7_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp7.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")

exp1_sct_counts2 <- left_join(x=exp1_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp2_sct_counts2 <- left_join(x=exp2_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp6_sct_counts2 <- left_join(x=exp6_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp7_sct_counts2 <- left_join(x=exp7_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))

#no NA
exp1_sct_counts2_noNA <- exp1_sct_counts2[t.type.Under1Over1 %in% c("WT", "KO")]
exp2_sct_counts2_noNA <- exp2_sct_counts2[t.type.Under1Over1 %in% c("WT", "KO")]
exp6_sct_counts2_noNA <- exp6_sct_counts2[t.type.Under1Over1 %in% c("WT", "KO")]
exp7_sct_counts2_noNA <- exp7_sct_counts2[t.type.Under1Over1 %in% c("WT", "KO")]

#####expression of cells in L2/3 of V1 only
L23_V1_WT_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
L23_V1_KO_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

L23_V1_WT_avgExp <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
L23_V1_KO_avgExp <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

#L23 locations
exp1_L23_L <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_L_VIS_L23.csv")
exp1_L23_R <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_R_VIS_L23.csv")
names(exp1_L23_L)="index"
names(exp1_L23_R)="index"
exp1_L23_cells <- c(exp1_L23_L$index, exp1_L23_R$index)


exp2_L23_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS_L23.csv")
exp2_L23_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS_L23.csv")
names(exp2_L23_L)="index"
names(exp2_L23_R)="index"
exp2_L23_cells <- c(exp2_L23_L$index, exp2_L23_R$index)

exp6_L23_L <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_L_VIS_L23.csv")
exp6_L23_R <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_R_VIS_L23.csv")
names(exp6_L23_L)="index"
names(exp6_L23_R)="index"
exp6_L23_cells <- c(exp6_L23_L$index, exp6_L23_R$index)

exp7_L23_L <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_L_VIS_L23.csv")
exp7_L23_R <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_R_VIS_L23.csv")
names(exp7_L23_L)="index"
names(exp7_L23_R)="index"
exp7_L23_cells <- c(exp7_L23_L$index, exp7_L23_R$index)

#MERFISH gene panel table
gene_panel_classes = fread("HG_lab/Vizgen/Gene_lists_for_analysis/MERFISH_gene_panel_classes.txt")
#coding genes
ensgene_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed", header=FALSE)$V4
#L2/3 genes
L23_genes <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.all.genes.txt", header=FALSE)$V1
L23_genes_mm9 <- intersect(ensgene_mm9, L23_genes)
non_L23_genes <- setdiff(ensgene_mm9, L23_genes_mm9)

#L2/3 A, B, and C genes
L23_genes_A <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.A.genes.txt", header=FALSE)$V1
L23_genes_B <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.B.genes.txt", header=FALSE)$V1
L23_genes_C <- fread("HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.C.genes.txt", header=FALSE)$V1

#MERFISH gene panel A/B/C genes
L23_all_merfish <- intersect(gene_panel_classes$Gene, L23_genes)
L23_A_merfish <- intersect(gene_panel_classes$Gene, L23_genes_A)
L23_B_merfish <- intersect(gene_panel_classes$Gene, L23_genes_B)
L23_C_merfish <- intersect(gene_panel_classes$Gene, L23_genes_C)



#formatting L23 gene names to fit data.frame standards
L23_genes_newNames <- make.names(L23_genes)
L23_A_merfish_newNames <- make.names(L23_A_merfish)
L23_B_merfish_newNames <- make.names(L23_B_merfish)
L23_C_merfish_newNames <- make.names(L23_C_merfish)

"A/C ratio"

L23_V1_WT_avgExp_matrix[Gene %in% L23_A_merfish_newNames, Gene]
#WT
L23_V1_WT_avgExp_A_genes_summ <- group_by(L23_V1_WT_avgExp[gene %in% L23_A_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table

L23_V1_WT_avgExp_B_genes_summ <- group_by(L23_V1_WT_avgExp[gene %in% L23_B_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table

L23_V1_WT_avgExp_C_genes_summ <- group_by(L23_V1_WT_avgExp[gene %in% L23_C_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table
#KO
L23_V1_KO_avgExp_A_genes_summ <- group_by(L23_V1_KO_avgExp[gene %in% L23_A_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table

L23_V1_KO_avgExp_B_genes_summ <- group_by(L23_V1_KO_avgExp[gene %in% L23_B_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table

L23_V1_KO_avgExp_C_genes_summ <- group_by(L23_V1_KO_avgExp[gene %in% L23_C_merfish_newNames], subtype) %>%
  summarise(mean.Exp = mean(avgExp)) %>% data.table

##
L23_V1_WT_KO_avgExp_A_genes_summ <- inner_join(x=L23_V1_WT_avgExp_A_genes_summ, y=L23_V1_KO_avgExp_A_genes_summ, by="subtype", suffix = c(".WT", ".KO"))
L23_V1_WT_KO_avgExp_B_genes_summ <- inner_join(x=L23_V1_WT_avgExp_B_genes_summ, y=L23_V1_KO_avgExp_B_genes_summ, by="subtype", suffix = c(".WT", ".KO"))
L23_V1_WT_KO_avgExp_C_genes_summ <- inner_join(x=L23_V1_WT_avgExp_C_genes_summ, y=L23_V1_KO_avgExp_C_genes_summ, by="subtype", suffix = c(".WT", ".KO"))

L23_V1_WT_KO_avgExp_A_genes_summ$mean.Exp.WT.KO.ratio <- L23_V1_WT_KO_avgExp_A_genes_summ[,]

All_Vis_layerdepth <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/layer_depth/All_Vis_layerdepth_fixed.csv")

V1_layerdepth <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/layer_depth/All_V1_layerdepth_fixed.csv")



#L2/3 cells
allExp_L23_cells= fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/exp1_exp2_exp6_exp7_L23_cells.csv")
#L2/3 cells in VIS only
allExp_L23_V1_cells <- intersect(allExp_L23_cells$cell.id, V1_layerdepth$cell.id)
allExp_L23_Vis_cells <- intersect(allExp_L23_cells$cell.id, All_Vis_layerdepth $cell.id)


V1_layerdepth_df <- data.frame(V1_layerdepth, row.names="cell.id")

V1_layerdepth_df$avg_A_sctCount <-  rowMeans(V1_layerdepth_df[, L23_A_merfish_newNames])
V1_layerdepth_df$avg_C_sctCount <-  rowMeans(V1_layerdepth_df[, L23_C_merfish_newNames])

V1_layerdepth_df$avg_A_C_ratio <- V1_layerdepth_df$avg_A_sctCount / V1_layerdepth_df$avg_C_sctCount

L23_V1_layerdepth_dt <- data.table(V1_layerdepth_df[allExp_L23_V1_cells,], keep.rownames="cell.id")
L23_V1_layerdepth_dt_WTKO <- L23_V1_layerdepth_dt[t.type.Under1Over1 %in% c("WT", "KO")]



Vis_layerdepth_df <- data.frame(All_Vis_layerdepth, row.names="cell.id")

Vis_layerdepth_df$avg_A_sctCount <-  rowMeans(Vis_layerdepth_df[, L23_A_merfish_newNames])
Vis_layerdepth_df$avg_C_sctCount <-  rowMeans(Vis_layerdepth_df[, L23_C_merfish_newNames])

Vis_layerdepth_df$avg_A_C_ratio <- Vis_layerdepth_df$avg_A_sctCount / Vis_layerdepth_df$avg_C_sctCount

L23_Vis_layerdepth_dt <- data.table(Vis_layerdepth_df[allExp_L23_Vis_cells,], keep.rownames="cell.id")
L23_Vis_layerdepth_dt_WTKO <- L23_Vis_layerdepth_dt[t.type.Under1Over1 %in% c("WT", "KO")]


ggplot(L23_V1_layerdepth_dt_WTKO[(predicted.id=="166_L2/3 IT CTX") & (rep=="Rep2")], aes(x = as.numeric(NDR), y = as.numeric(avg_A_C_ratio), color=t.type.Under1Over1))+
  ggtitle("166_L2/3 IT CTX")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("SCTransform-corrected count ratio\n(superficial/deep)") + xlab("Layer depth")+
  scale_fill_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
#ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23.Cheng2022.nonA.vs.A.genes.L23.V1.subtypeAgg.nCells10.pseudobulkDGE.logfc.boxplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
#ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23.Cheng2022.nonA.vs.A.genes.L23.V1.subtypeAgg.nCells10.pseudobulkDGE.logfc.boxplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


exp2_L23_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS_L23.csv")
exp2_L23_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS_L23.csv")
names(exp2_L23_L)="index"
names(exp2_L23_R)="index"
exp2_L23_cells <- c(exp2_L23_L$index, exp2_L23_R$index)

plot(L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index),center_x], L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index),center_y],pch = 20,cex=.3,col = "gray",xlab="x", ylab="y", asp=1)
points(L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & NDR>0.8,center_x],L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & NDR>0.8,center_y],pch = 20,cex=0.5, col="red")
points(L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & NDR>0.7,center_x],L23_V1_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & NDR>0.7,center_y],pch = 20,cex=0.5, col="blue")

plot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="166_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="166_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "purple", xlab="NDR", ylab="SCTransform-corrected count ratio\n(superficial/deep)", main="166_L2/3 IT CTX, L2/3 Vis, Exp2")
points(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="166_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="166_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "orange")

plot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "purple", xlab="NDR", ylab="SCTransform-corrected count ratio\n(superficial/deep)", main="171_L2/3 IT CTX, L2/3 Vis, Exp2")
points(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "orange")

plot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "purple", xlab="NDR", ylab="SCTransform-corrected count ratio\n(superficial/deep)", main="171_L2/3 IT CTX, L2/3 Vis, Exp2")
points(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "orange")

plot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="179_L4 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="179_L4 IT CTX") & (t.type.Under1Over1=="WT"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "purple", xlab="NDR", ylab="SCTransform-corrected count ratio\n(superficial/deep)", main="179_L4 IT CTX, L2/3 Vis, Exp2")
points(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="179_L4 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(NDR)], L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="179_L4 IT CTX") & (t.type.Under1Over1=="KO"), as.numeric(avg_A_C_ratio)], pch = 19, cex=1, col = "orange")


ggplot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="166_L2/3 IT CTX")], aes(x = as.numeric(NDR), y = as.numeric(avg_A_C_ratio), color=t.type.Under1Over1))+
  ggtitle("166_L2/3 IT CTX, Exp2 L, layer 2/3, Vis")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("SCTransform-corrected count ratio\n(superficial/deep)") + xlab("Layer depth")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.166.L23.IT.CTX.L23.Vis.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.166.L23.IT.CTX.L23.Vis.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="171_L2/3 IT CTX")], aes(x = as.numeric(NDR), y = as.numeric(avg_A_C_ratio), color=t.type.Under1Over1))+
  ggtitle("171_L2/3 IT CTX, Exp2 L, layer 2/3, Vis")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("SCTransform-corrected count ratio\n(superficial/deep)") + xlab("Layer depth")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.171.L23.IT.CTX.L23.Vis.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.171.L23.IT.CTX.L23.Vis.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="168_L2/3 IT CTX")], aes(x = as.numeric(NDR), y = as.numeric(avg_A_C_ratio), color=t.type.Under1Over1))+
  ggtitle("168_L2/3 IT CTX, Exp2 L, layer 2/3, Vis")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("SCTransform-corrected count ratio\n(superficial/deep)") + xlab("Layer depth")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.168.L23.IT.CTX.L23.Vis.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.168.L23.IT.CTX.L23.Vis.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(L23_Vis_layerdepth_dt_WTKO[(cell.id %in% exp2_L23_L$index) & (predicted.id=="179_L4 IT CTX")], aes(x = as.numeric(NDR), y = as.numeric(avg_A_C_ratio), color=t.type.Under1Over1))+
  ggtitle("179_L4 IT CTX, Exp2 L, layer 2/3, Vis")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("SCTransform-corrected count ratio\n(superficial/deep)") + xlab("Layer depth")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.179.L4.IT.CTX.L23.Vis.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.ratio.sctCounts.179.L4.IT.CTX.L23.Vis.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

#####
V1_layerdepth2 <- left_join(V1_layerdepth, CTX_HIP_annot[, .(cluster_label, class_label)], by=c("predicted.id"="cluster_label"))
#exp2
exp2_L23_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS_L23.csv")
exp2_L23_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS_L23.csv")
names(exp2_L23_L)="index"
names(exp2_L23_R)="index"
exp2_L23_cells <- c(exp2_L23_L$index, exp2_L23_R$index)

exp2_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp2_sct_counts2 <- left_join(x=exp2_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp2_sct_counts_V1 <- exp2_sct_counts2[index %in% V1_layerdepth2$cell.id,]
exp2_sct_counts_V1_glut <- exp2_sct_counts_V1[class_label=="Glutamatergic",]

exp2_sct_counts_V1_glut_L23 <- exp2_sct_counts_V1_glut[index %in% exp2_L23_cells]
exp2_sct_counts_V1_glut_L23_genesOnly <- data.frame(exp2_sct_counts_V1_glut_L23[, 1:551], row.names="index")
exp2_sct_counts_V1_glut_L23_scaled <- scale(exp2_sct_counts_V1_glut_L23_genesOnly)
join_cols_exp2_sct_counts_V1_glut_L23 <- c("index", names(exp2_sct_counts_V1_glut_L23)[552:573])

exp2_sct_counts_V1_glut_L23_scaled_dt <- data.frame(exp2_sct_counts_V1_glut_L23_scaled) %>% data.table(keep.rownames="index")
exp2_sct_counts_V1_glut_L23_scaled_dt <- inner_join(x=exp2_sct_counts_V1_glut_L23_scaled_dt, y=exp2_sct_counts_V1_glut_L23[, ..join_cols_exp2_sct_counts_V1_glut_L23], by="index")

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio <- exp2_sct_counts_V1_glut_L23_scaled_dt
exp2_sct_counts_V1_glut_L23_scaled_AC_ratio$A_avg <- rowMeans(exp2_sct_counts_V1_glut_L23_scaled_dt[, ..L23_A_merfish_newNames])
exp2_sct_counts_V1_glut_L23_scaled_AC_ratio$C_avg <- rowMeans(exp2_sct_counts_V1_glut_L23_scaled_dt[, ..L23_C_merfish_newNames])
exp2_sct_counts_V1_glut_L23_scaled_AC_ratio$AC_ratio <- exp2_sct_counts_V1_glut_L23_scaled_AC_ratio[, A_avg/C_avg]

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio <- left_join(x=exp2_sct_counts_V1_glut_L23_scaled_AC_ratio, y=V1_layerdepth2[, .(cell.id, NDR)], by=c("index"="cell.id"))

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_L <- exp2_sct_counts_V1_glut_L23_scaled_AC_ratio[index %in% exp2_L23_L$index]
exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_R <- exp2_sct_counts_V1_glut_L23_scaled_AC_ratio[index %in% exp2_L23_R$index]

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_L$layer_decile <- ntile(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_L$NDR, 10)
exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_R$layer_decile <- ntile(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_R$NDR, 10)

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both <- rbind(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_L, exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_R)

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_ACordered <- exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both[order(AC_ratio)]

exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both_summ <- group_by(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both[t.type.Under1Over1 %in% c("WT", "KO")], layer_decile, t.type.Under1Over1, predicted.id) %>%
  summarise(avg_AC_ratio = mean(AC_ratio)) %>% data.table

ggplot(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both_summ[predicted.id=="166_L2/3 IT CTX"], aes(x = as.numeric(layer_decile), y = as.numeric(avg_AC_ratio), color=t.type.Under1Over1))+
  ggtitle("166_L2/3 IT CTX, L2/3 V1 glutamatergic cells, Exp2")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("RNA SCTransform-corrected count z-score ratio\n(superficial/deep)") + xlab("Layer depth decile")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.166.L23.IT.CTX.L23.V1.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.166.L23.IT.CTX.L23.V1.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

ggplot(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both_summ[predicted.id=="171_L2/3 IT CTX"], aes(x = as.numeric(layer_decile), y = as.numeric(avg_AC_ratio), color=t.type.Under1Over1))+
  ggtitle("171_L2/3 IT CTX, L2/3 V1 glutamatergic cells, Exp2")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("RNA SCTransform-corrected count z-score ratio\n(superficial/deep)") + xlab("Layer depth decile")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.171.L23.IT.CTX.L23.V1.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.171.L23.IT.CTX.L23.V1.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


ggplot(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both_summ[predicted.id=="168_L2/3 IT CTX"], aes(x = as.numeric(layer_decile), y = as.numeric(avg_AC_ratio), color=t.type.Under1Over1))+
  ggtitle("168_L2/3 IT CTX, L2/3 V1 glutamatergic cells, Exp2")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("RNA SCTransform-corrected count z-score ratio\n(superficial/deep)") + xlab("Layer depth decile")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.168.L23.IT.CTX.L23.V1.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.168.L23.IT.CTX.L23.V1.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')


ggplot(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both_summ[predicted.id=="179_L4 IT CTX"], aes(x = as.numeric(layer_decile), y = as.numeric(avg_AC_ratio), color=t.type.Under1Over1))+
  ggtitle("179_L4 IT CTX, L2/3 V1 glutamatergic cells, Exp2")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("RNA SCTransform-corrected count z-score ratio\n(superficial/deep)") + xlab("Layer depth decile")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.179.L4.IT.CTX.L23.V1.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.zScore.ratio.179.L4.IT.CTX.L23.V1.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')

library(zoo)


ggplot(exp2_sct_counts_V1_glut_L23_scaled_AC_ratio_both, aes(x=layer_decile, y=AC_ratio)) + 
  geom_point(position=position_jitter(1,3), pch=21, fill="#FF0000AA") +
  geom_line(aes(y=rollmean(AC_ratio, 7, na.pad=TRUE))) +
  theme_bw()
####

#####expression of cells in whole cortex and hippocampus
Vis_WT_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
Vis_KO_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

Vis_WT_avgExp <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
Vis_KO_avgExp <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/Vis.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")


Vis_layerdepth2 <- left_join(All_Vis_layerdepth, CTX_HIP_annot[, .(cluster_label, class_label)], by=c("predicted.id"="cluster_label"))
#exp1
exp1_Vis_L <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_L_VIS.csv")
exp1_Vis_R <- fread("HG_lab/Vizgen/experiment1_2022-08-18/Vizualizer_experiments/Liu 2023 exported regions/exp1_Liu_R_VIS.csv")
names(exp1_Vis_L)="index"
names(exp1_Vis_R)="index"
exp1_Vis_cells <- c(exp1_Vis_L$index, exp1_Vis_R$index)

#exp2
exp2_Vis_L <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_L_VIS.csv")
exp2_Vis_R <- fread("HG_lab/Vizgen/experiment2_2022-09-08/Vizualizer_experiments/Liu 2023 exported regions/exp2Liu_R_VIS.csv")
names(exp2_Vis_L)="index"
names(exp2_Vis_R)="index"
exp2_Vis_cells <- c(exp2_Vis_L$index, exp2_Vis_R$index)

#exp6
exp6_Vis_L <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_L_VIS.csv")
exp6_Vis_R <- fread("HG_lab/Vizgen/experiment6_2023-01-23/Vizualizer_experiments/Liu 2023 exported regions/exp6_Liu_R_VIS.csv")
names(exp6_Vis_L)="index"
names(exp6_Vis_R)="index"
exp6_Vis_cells <- c(exp6_Vis_L$index, exp6_Vis_R$index)

#exp7
exp7_Vis_L <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_L_VIS.csv")
exp7_Vis_R <- fread("HG_lab/Vizgen/experiment7_2023-01-27/Vizualizer_experiments/Liu 2023 exported regions/exp7_Liu2023_R_VIS.csv")
names(exp7_Vis_L)="index"
names(exp7_Vis_R)="index"
exp7_Vis_cells <- c(exp7_Vis_L$index, exp7_Vis_R$index)

#exp1
exp1_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp1.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp1_sct_counts2 <- left_join(x=exp1_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp1_sct_counts_Vis <- exp1_sct_counts2[index %in% Vis_layerdepth2$cell.id,]
exp1_sct_counts_Vis_glut <- exp1_sct_counts_Vis[class_label=="Glutamatergic"]
exp1_sct_counts_Vis_glut_genesOnly <- data.frame(exp1_sct_counts_Vis_glut[, 1:551], row.names="index")
exp1_sct_counts_Vis_glut_scaled <- scale(exp1_sct_counts_Vis_glut_genesOnly)
join_cols_exp1_sct_counts_Vis_glut <- c("index", names(exp1_sct_counts_Vis_glut)[552:573])
exp1_sct_counts_Vis_glut_scaled_dt <- data.frame(exp1_sct_counts_Vis_glut_scaled) %>% data.table(keep.rownames="index")
exp1_sct_counts_Vis_glut_scaled_dt <- inner_join(x=exp1_sct_counts_Vis_glut_scaled_dt, y=exp1_sct_counts_Vis_glut[, ..join_cols_exp1_sct_counts_Vis_glut], by="index")

exp1_sct_counts_Vis_glut_scaled_AC_ratio <- left_join(x=exp1_sct_counts_Vis_glut_scaled_dt, y=Vis_layerdepth2[, .(cell.id, NDR)], by=c("index"="cell.id"))
exp1_sct_counts_Vis_glut_scaled_AC_ratio$A_avg <- rowMeans(exp1_sct_counts_Vis_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp1_sct_counts_Vis_glut_scaled_AC_ratio$C_avg <- rowMeans(exp1_sct_counts_Vis_glut_scaled_dt[, ..L23_C_merfish_newNames])



#exp2
exp2_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp2.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp2_sct_counts2 <- left_join(x=exp2_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp2_sct_counts_Vis <- exp2_sct_counts2[index %in% Vis_layerdepth2$cell.id,]
exp2_sct_counts_Vis_glut <- exp2_sct_counts_Vis[class_label=="Glutamatergic"]
exp2_sct_counts_Vis_glut_genesOnly <- data.frame(exp2_sct_counts_Vis_glut[, 1:551], row.names="index")
exp2_sct_counts_Vis_glut_scaled <- scale(exp2_sct_counts_Vis_glut_genesOnly)
join_cols_exp2_sct_counts_Vis_glut <- c("index", names(exp2_sct_counts_Vis_glut)[552:573])
exp2_sct_counts_Vis_glut_scaled_dt <- data.frame(exp2_sct_counts_Vis_glut_scaled) %>% data.table(keep.rownames="index")
exp2_sct_counts_Vis_glut_scaled_dt <- inner_join(x=exp2_sct_counts_Vis_glut_scaled_dt, y=exp2_sct_counts_Vis_glut[, ..join_cols_exp2_sct_counts_Vis_glut], by="index")

exp2_sct_counts_Vis_glut_scaled_AC_ratio <- left_join(x=exp2_sct_counts_Vis_glut_scaled_dt, y=Vis_layerdepth2[, .(cell.id, NDR)], by=c("index"="cell.id"))
exp2_sct_counts_Vis_glut_scaled_AC_ratio$A_avg <- rowMeans(exp2_sct_counts_Vis_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp2_sct_counts_Vis_glut_scaled_AC_ratio$C_avg <- rowMeans(exp2_sct_counts_Vis_glut_scaled_dt[, ..L23_C_merfish_newNames])


#exp6
exp6_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp6.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp6_sct_counts2 <- left_join(x=exp6_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp6_sct_counts_Vis <- exp6_sct_counts2[index %in% Vis_layerdepth2$cell.id,]
exp6_sct_counts_Vis_glut <- exp6_sct_counts_Vis[class_label=="Glutamatergic"]
exp6_sct_counts_Vis_glut_genesOnly <- data.frame(exp6_sct_counts_Vis_glut[, 1:551], row.names="index")
exp6_sct_counts_Vis_glut_scaled <- scale(exp6_sct_counts_Vis_glut_genesOnly)
join_cols_exp6_sct_counts_Vis_glut <- c("index", names(exp6_sct_counts_Vis_glut)[552:573])
exp6_sct_counts_Vis_glut_scaled_dt <- data.frame(exp6_sct_counts_Vis_glut_scaled) %>% data.table(keep.rownames="index")
exp6_sct_counts_Vis_glut_scaled_dt <- inner_join(x=exp6_sct_counts_Vis_glut_scaled_dt, y=exp6_sct_counts_Vis_glut[, ..join_cols_exp6_sct_counts_Vis_glut], by="index")

exp6_sct_counts_Vis_glut_scaled_AC_ratio <- left_join(x=exp6_sct_counts_Vis_glut_scaled_dt, y=Vis_layerdepth2[, .(cell.id, NDR)], by=c("index"="cell.id"))
exp6_sct_counts_Vis_glut_scaled_AC_ratio$A_avg <- rowMeans(exp6_sct_counts_Vis_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp6_sct_counts_Vis_glut_scaled_AC_ratio$C_avg <- rowMeans(exp6_sct_counts_Vis_glut_scaled_dt[, ..L23_C_merfish_newNames])


#exp7
exp7_sct_counts <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/exp7.100vol.300counts.CTX.HC.pred0.2.sctCount.celltype.summary.table.csv")
exp7_sct_counts2 <- left_join(x=exp7_sct_counts, y=CTX_HIP_annot[,.(cluster_label, subclass_label, neighborhood_label, class_label, cluster_color, subclass_color, neighborhood_color, class_color)], by=c("predicted.id"="cluster_label", "subclass_label", "neighborhood_label"))
exp7_sct_counts_Vis <- exp7_sct_counts2[index %in% Vis_layerdepth2$cell.id,]
exp7_sct_counts_Vis_glut <- exp7_sct_counts_Vis[class_label=="Glutamatergic"]
exp7_sct_counts_Vis_glut_genesOnly <- data.frame(exp7_sct_counts_Vis_glut[, 1:551], row.names="index")
exp7_sct_counts_Vis_glut_scaled <- scale(exp7_sct_counts_Vis_glut_genesOnly)
join_cols_exp7_sct_counts_Vis_glut <- c("index", names(exp7_sct_counts_Vis_glut)[552:573])
exp7_sct_counts_Vis_glut_scaled_dt <- data.frame(exp7_sct_counts_Vis_glut_scaled) %>% data.table(keep.rownames="index")
exp7_sct_counts_Vis_glut_scaled_dt <- inner_join(x=exp7_sct_counts_Vis_glut_scaled_dt, y=exp7_sct_counts_Vis_glut[, ..join_cols_exp7_sct_counts_Vis_glut], by="index")

exp7_sct_counts_Vis_glut_scaled_AC_ratio <- left_join(x=exp7_sct_counts_Vis_glut_scaled_dt, y=Vis_layerdepth2[, .(cell.id, NDR)], by=c("index"="cell.id"))
exp7_sct_counts_Vis_glut_scaled_AC_ratio$A_avg <- rowMeans(exp7_sct_counts_Vis_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp7_sct_counts_Vis_glut_scaled_AC_ratio$C_avg <- rowMeans(exp7_sct_counts_Vis_glut_scaled_dt[, ..L23_C_merfish_newNames])



####combine
allExps_sct_counts_Vis_glut_scaled_AC_ratio <- rbind(exp1_sct_counts_Vis_glut_scaled_AC_ratio, exp2_sct_counts_Vis_glut_scaled_AC_ratio, exp6_sct_counts_Vis_glut_scaled_AC_ratio, exp7_sct_counts_Vis_glut_scaled_AC_ratio)
#separate by WT and KO, order by layer depth (NDR)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered <- allExps_sct_counts_Vis_glut_scaled_AC_ratio[(t.type.Under1Over1=="WT")][order(NDR)]
allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered <- allExps_sct_counts_Vis_glut_scaled_AC_ratio[(t.type.Under1Over1=="KO")][order(NDR)]

#subtract minimum A z-score from each A zscore (and same with C z-score) to make new z-score minimum 0 to avoid negatives
allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$nonneg_A <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$A_avg - min(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$A_avg)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$nonneg_C <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$C_avg - min(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$C_avg)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$nonneg_AC_ratio <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered[, nonneg_A/nonneg_C]

allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$nonneg_A <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$A_avg - min(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$A_avg)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$nonneg_C <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$C_avg - min(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$C_avg)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$nonneg_AC_ratio <- allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered[, nonneg_A/nonneg_C]

#calculating rolling medians
window_size <- 150

# WT, Create a sequence of starting indices for the rolling median calculation
start_indices_WT <- seq(1, nrow(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered) - window_size + 1)

# Calculate rolling median of NDR and non-negative A/C z-score ratio using base R and data.table
nonneg_AC_rolling_medians_WT <- sapply(start_indices_WT, function(i) {
  median(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$nonneg_AC_ratio[i:(i + window_size - 1)])
})

NDR_rolling_medians_WT <- sapply(start_indices_WT, function(i) {
  median(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered$NDR[i:(i + window_size - 1)])
})

# KO, Create a sequence of starting indices for the rolling median calculation
start_indices_KO <- seq(1, nrow(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered) - window_size + 1)

# KO, Calculate rolling median of NDR and non-negative A/C z-score ratio using base R and data.table
nonneg_AC_rolling_medians_KO <- sapply(start_indices_KO, function(i) {
  median(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$nonneg_AC_ratio[i:(i + window_size - 1)])
})

NDR_rolling_medians_KO <- sapply(start_indices_KO, function(i) {
  median(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered$NDR[i:(i + window_size - 1)])
})

allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered_rollMed=cbind(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered[1:length(nonneg_AC_rolling_medians_WT)], nonneg_AC_ratio_rollMed=nonneg_AC_rolling_medians_WT, NDR_rollMed=NDR_rolling_medians_WT)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered_rollMed=cbind(allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered[1:length(nonneg_AC_rolling_medians_KO)], nonneg_AC_ratio_rollMed=nonneg_AC_rolling_medians_KO, NDR_rollMed=NDR_rolling_medians_KO)


allExps_sct_counts_Vis_glut_scaled_AC_ratio_both_rollMed <- rbind(allExps_sct_counts_Vis_glut_scaled_AC_ratio_WT_ordered_rollMed, allExps_sct_counts_Vis_glut_scaled_AC_ratio_KO_ordered_rollMed)
allExps_sct_counts_Vis_glut_scaled_AC_ratio_both_rollMed = allExps_sct_counts_Vis_glut_scaled_AC_ratio_both_rollMed %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))



ggplot(allExps_sct_counts_Vis_glut_scaled_AC_ratio_both_rollMed, aes(x = as.numeric(NDR_rollMed), y = as.numeric(nonneg_AC_ratio_rollMed), color=t.type.Under1Over1))+
  ggtitle("Glutamatergic\nMecp2 KO/+ MERFISH Vis")+
  geom_line()+
  #coord_cartesian(ylim=c(-1.3, 1.3))+
  ylab("Rolling median RNA SCTransform-corrected count\nnon-negative z-score ratio (superficial/deep)") + xlab("Rolling median NDR (deep to superficial)")+
  scale_color_manual(name="Transcriptotype:", values = c("WT"="purple", "KO"="orange")) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.nonneg.zScore.ratio.by.NDR.Liu2023.Vis.mecp2Het.CTX.HC.lineplot.png", width = 5, height = 5 , dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/lineplots/L23.Cheng2022.A.C.sctCount.nonneg.zScore.ratio.by.NDR.Liu2023.Vis.mecp2Het.CTX.HC.lineplot.eps", width = 5, height = 5 , dpi = 300, units = "in", device='eps')



###
#file with layer depth values annotated for each cell in V1 L2/3
all_V1_L23_depths=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_V1_L23_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")
AC_avgZscore_func <- function(exp_sct_count_table, depth_table){
  #subset of cells in depth table
  exp_sct_count_table_dep <- exp_sct_count_table[index %in% depth_table$X,]
  #Glutamatergic cells
  exp_sct_count_table_dep_glut <- exp_sct_count_table_dep[class_label=="Glutamatergic"]
  exp_sct_count_table_dep_glut_genesOnly <- data.frame(exp_sct_count_table_dep_glut[, 1:551], row.names="index")
  exp_sct_count_table_dep_glut_scaled <- scale(exp_sct_count_table_dep_glut_genesOnly)
  
  exp_sct_count_table_dep_glut_scaled_dt <- data.frame(exp_sct_count_table_dep_glut_scaled) %>% data.table(keep.rownames="index")
  exp_sct_count_table_dep_glut_scaled_dt <- left_join(x=exp_sct_count_table_dep_glut_scaled_dt, y=depth_table, by=c("index"="X"))
  
  exp_sct_count_table_dep_glut_scaled_AC_avgZscore <- copy(exp_sct_count_table_dep_glut_scaled_dt)
  exp_sct_count_table_dep_glut_scaled_AC_avgZscore$A_avg <- rowMeans(exp_sct_count_table_dep_glut_scaled_dt[, ..L23_A_merfish_newNames])
  exp_sct_count_table_dep_glut_scaled_AC_avgZscore$C_avg <- rowMeans(exp_sct_count_table_dep_glut_scaled_dt[, ..L23_C_merfish_newNames])
  return(exp_sct_count_table_dep_glut_scaled_AC_avgZscore)
}
#AC ratio for each experiment with new layer quintiles in V1 L2/3
exp1_sct_counts_V1_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp1_sct_counts2, depth_table=all_V1_L23_depths)
exp2_sct_counts_V1_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp2_sct_counts2, depth_table=all_V1_L23_depths)
exp6_sct_counts_V1_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp6_sct_counts2, depth_table=all_V1_L23_depths)
exp7_sct_counts_V1_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp7_sct_counts2, depth_table=all_V1_L23_depths)

#
all_sct_counts_V1_L23_glut_scaled_AC_avgZscore <- rbind(exp1_sct_counts_V1_L23_glut_scaled_AC_avgZscore, exp2_sct_counts_V1_L23_glut_scaled_AC_avgZscore, exp6_sct_counts_V1_L23_glut_scaled_AC_avgZscore, exp7_sct_counts_V1_L23_glut_scaled_AC_avgZscore)
all_sct_counts_V1_L23_glut_scaled_AC_avgZscore = all_sct_counts_V1_L23_glut_scaled_AC_avgZscore %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_V1_L23_glut_scaled_AC_avgZscore = all_sct_counts_V1_L23_glut_scaled_AC_avgZscore %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_V1_L23_glut_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(A_avg), fill=t.type.Under1Over1))+
  ggtitle("Superficial genes in V1 L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\nsuperficial gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_V1_L23_glut_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(C_avg), fill=t.type.Under1Over1))+
  ggtitle("Deep genes in V1 L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\ndeep gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

#Seurat object
#mecp2Het.CTX.HC.seu <-  readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_allTaxLabels_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

#mecp2Het.CTX.HC.seu.avgExp = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT


exp7_sct_counts_V1_L23 <- exp7_sct_counts2[index %in% all_V1_L23_depths$X,]
exp7_sct_counts_V1_L23_glut <- exp7_sct_counts_V1_L23[class_label=="Glutamatergic"]
exp7_sct_counts_V1_L23_glut_dep <- left_join(x=exp7_sct_counts_V1_L23_glut[, 1:551], y=all_V1_L23_depths, by=c("index"="X"))

exp7_sct_counts_V1_L23_glut_dep$A_avg_sctCount <- rowMeans(exp7_sct_counts_V1_L23_glut_dep[, ..L23_A_merfish_newNames])
exp7_sct_counts_V1_L23_glut_dep$C_avg_sctCount <- rowMeans(exp7_sct_counts_V1_L23_glut_dep[, ..L23_C_merfish_newNames])
exp7_sct_counts_V1_L23_glut_dep$AC_sctCount_ratio <- exp7_sct_counts_V1_L23_glut_dep[, A_avg_sctCount/C_avg_sctCount]



#function for calculating ratio of average SCTransform-corrected count of L23 Cheng A genes to L23 Cheng C genes
AC_sctCount_ratio_func <- function(exp_sct_count_table, depth_table){
  exp_sct_count_table_dep <- exp_sct_count_table[index %in% depth_table$X,]
  exp_sct_count_table_dep_glut <- exp_sct_count_table_dep[class_label=="Glutamatergic"]
  exp_sct_count_table_dep_glut <- left_join(x=exp_sct_count_table_dep_glut[, 1:551], y=depth_table, by=c("index"="X"))
  
  exp_sct_count_table_dep_glut$A_avg_sctCount <- rowMeans(exp_sct_count_table_dep_glut[, ..L23_A_merfish_newNames])
  exp_sct_count_table_dep_glut$C_avg_sctCount <- rowMeans(exp_sct_count_table_dep_glut[, ..L23_C_merfish_newNames])
  exp_sct_count_table_dep_glut$AC_sctCount_ratio <- exp_sct_count_table_dep_glut[, A_avg_sctCount/C_avg_sctCount]
  return(exp_sct_count_table_dep_glut)
}

exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp1_sct_counts2, depth_table=all_V1_L23_depths)
exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp2_sct_counts2, depth_table=all_V1_L23_depths)
exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp6_sct_counts2, depth_table=all_V1_L23_depths)
exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp7_sct_counts2, depth_table=all_V1_L23_depths)
all_sct_counts_V1_L23_glut_AC_sctCount_ratio <- rbind(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio, exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio, exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio, exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio)

all_sct_counts_V1_L23_glut_AC_sctCount_ratio = all_sct_counts_V1_L23_glut_AC_sctCount_ratio %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_V1_L23_glut_AC_sctCount_ratio = all_sct_counts_V1_L23_glut_AC_sctCount_ratio %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_V1_L23_glut_AC_sctCount_ratio[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(AC_sctCount_ratio), fill=t.type.Under1Over1))+
  ggtitle("Superficial/deep gene ratio in V1 L2/3 MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,8))+
  ylab("Average SCTransform-corrected RNA count ratio\n(superficial/deep genes)") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/L23_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

####
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

AC_pal <- colorRampPalette(c('blue','red'))

exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered <- exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio[order(AC_sctCount_ratio)]
exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered <- exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio[order(AC_sctCount_ratio)]
exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered <- exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio[order(AC_sctCount_ratio)]
exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered <- exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio[order(AC_sctCount_ratio)]


exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[, get("AC_sctCount_ratio")], breaks=20))]
exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[, get("AC_sctCount_ratio")], breaks=20))]
exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[, get("AC_sctCount_ratio")], breaks=20))]
exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[, get("AC_sctCount_ratio")], breaks=20))]




plot(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L"),center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L"),center_y],pch = 20,cex=.3,col = "white",main="WT", xlab="x", ylab="y", asp=1)
for(i in exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L") &(t.type.Under1Over1=="WT"), index]){
  points(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,plot_color])
}

plot(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L"),center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L"),center_y],pch = 20,cex=.3,col = "white",main="KO", xlab="x", ylab="y", asp=1)
for(i in exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[(hemisphere=="L") &(t.type.Under1Over1=="KO"), index]){
  points(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered[index==i,plot_color])
}

AC_sctCount_ratio_outlier_func <- function(data_table){
  IQR = quantile(data_table$AC_sctCount_ratio)[4] - quantile(data_table$AC_sctCount_ratio)[2] 
  outlier_min <- quantile(data_table$AC_sctCount_ratio)[2] - 1.5*IQR
  outlier_max <- quantile(data_table$AC_sctCount_ratio)[4] + 1.5*IQR
  data_table$AC_sctCount_ratio_minMax <- data_table$AC_sctCount_ratio
  data_table[AC_sctCount_ratio > outlier_max, AC_sctCount_ratio_minMax := outlier_max]
  data_table[AC_sctCount_ratio < outlier_min, AC_sctCount_ratio_minMax := outlier_min]
  return(data_table)
}

exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax = AC_sctCount_ratio_outlier_func(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered)
exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax = AC_sctCount_ratio_outlier_func(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered)
exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax = AC_sctCount_ratio_outlier_func(exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered)
exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax = AC_sctCount_ratio_outlier_func(exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered)

exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")], breaks=20))]
exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")], breaks=20))]
exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")], breaks=20))]
exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax$plot_color <- get("AC_pal")(20)[as.numeric(cut(exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")], breaks=20))]


plot(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L"),center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L"),center_y],pch = 20,cex=.3,col = "white",main="Exp2, V1 L2/3 WT", xlab="x", ylab="y", asp=1)
for(i in exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L") &(t.type.Under1Over1=="WT"), index]){
  points(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,plot_color])
}

plot(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L"),center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L"),center_y],pch = 20,cex=.3,col = "white",main="Exp2, V1 L2/3 KO", xlab="x", ylab="y", asp=1)
for(i in exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[(hemisphere=="L") &(t.type.Under1Over1=="KO"), index]){
  points(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,center_x],exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[index==i,plot_color])
}


color.bar(get("AC_pal")(20), min=signif(min(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), max=signif(max(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), nticks=2, title="Exp2 V1 L2/3 superficial/deep gene\nSCTransform-corrected count ratio")


AC_sctCount_ratio_scatterplot_func <- function(data_table, plot_title, output_file, hemi, WT=TRUE, KO=TRUE){
  if(WT==TRUE){
    setEPS()
    postscript(paste0(output_file,"_WT.eps"))
    plot(data_table[(hemisphere==hemi),center_x],data_table[(hemisphere==hemi),center_y],pch = 20,cex=.3,col = "white",main=paste(plot_title, "WT"), xlab="x", ylab="y", asp=1)
    for(i in data_table[(hemisphere==hemi) &(t.type.Under1Over1=="WT"), index]){
      points(data_table[index==i,center_x],data_table[index==i,center_y],pch = 20,cex=2, col = data_table[index==i,plot_color])
    }
    dev.off()
  }

  if(KO==TRUE){
    setEPS()
    postscript(paste0(output_file,"_KO.eps"))
    plot(data_table[(hemisphere==hemi),center_x],data_table[(hemisphere==hemi),center_y],pch = 20,cex=.3,col = "white",main=paste(plot_title, "KO"), xlab="x", ylab="y", asp=1)
    for(i in data_table[(hemisphere==hemi) &(t.type.Under1Over1=="KO"), index]){
      points(data_table[index==i,center_x],data_table[index==i,center_y],pch = 20,cex=2, col = data_table[index==i,plot_color])
    }
    dev.off()
  }
}

AC_sctCount_ratio_scatterplot_func(data_table=exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax, hemi="L", WT=TRUE, KO=TRUE, plot_title="Exp1, V1 L2/3 glutamatergic\nleft hemisphere",
                                   output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/AC_ratio_avg_sctCount_plots/exp1_L_V1_L23_glut_AC_avg_sctCount_ratio_ordered_cex2_scatterplot")

AC_sctCount_ratio_scatterplot_func(data_table=exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax, hemi="L", WT=TRUE, KO=TRUE, plot_title="Exp2, V1 L2/3 glutamatergic\nleft hemisphere",
                                   output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/AC_ratio_avg_sctCount_plots/exp2_L_V1_L23_glut_AC_avg_sctCount_ratio_ordered_cex2_scatterplot")

AC_sctCount_ratio_scatterplot_func(data_table=exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax, hemi="L", WT=TRUE, KO=TRUE, plot_title="Exp6, V1 L2/3 glutamatergic\nleft hemisphere",
                                   output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/AC_ratio_avg_sctCount_plots/exp6_L_V1_L23_glut_AC_avg_sctCount_ratio_ordered_cex2_scatterplot")

AC_sctCount_ratio_scatterplot_func(data_table=exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax, hemi="L", WT=TRUE, KO=TRUE, plot_title="Exp7, V1 L2/3 glutamatergic\nleft hemisphere",
                                   output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/AC_ratio_avg_sctCount_plots/exp7_L_V1_L23_glut_AC_avg_sctCount_ratio_ordered_cex2_scatterplot")

color.bar(get("AC_pal")(20), min=signif(min(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), max=signif(max(exp1_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), nticks=2, title="Exp1 V1 L2/3 glutamatergic, left\nsuperficial/deep gene\nSCTransform-corrected count ratio")
color.bar(get("AC_pal")(20), min=signif(min(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), max=signif(max(exp2_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), nticks=2, title="Exp2 V1 L2/3 glutamatergic, left\nsuperficial/deep gene\nSCTransform-corrected count ratio")
color.bar(get("AC_pal")(20), min=signif(min(exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), max=signif(max(exp6_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), nticks=2, title="Exp6 V1 L2/3 glutamatergic, left\nsuperficial/deep gene\nSCTransform-corrected count ratio")
color.bar(get("AC_pal")(20), min=signif(min(exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), max=signif(max(exp7_sct_counts_V1_L23_glut_AC_sctCount_ratio_ACordered_minMax[, get("AC_sctCount_ratio_minMax")]),2), nticks=2, title="Exp7 V1 L2/3 glutamatergic, left\nsuperficial/deep gene\nSCTransform-corrected count ratio")


###VIS version of superficial/deep gene ratio plots
#file with layer depth values annotated for each cell in V1 L2/3
all_VIS_L23_depths=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_VIS_L23_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")

#AC ratio for each experiment with new layer quintiles in VIS L2/3
exp1_sct_counts_VIS_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp1_sct_counts2, depth_table=all_VIS_L23_depths)
exp2_sct_counts_VIS_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp2_sct_counts2, depth_table=all_VIS_L23_depths)
exp6_sct_counts_VIS_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp6_sct_counts2, depth_table=all_VIS_L23_depths)
exp7_sct_counts_VIS_L23_glut_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp7_sct_counts2, depth_table=all_VIS_L23_depths)

#
all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore <- rbind(exp1_sct_counts_VIS_L23_glut_scaled_AC_avgZscore, exp2_sct_counts_VIS_L23_glut_scaled_AC_avgZscore, exp6_sct_counts_VIS_L23_glut_scaled_AC_avgZscore, exp7_sct_counts_VIS_L23_glut_scaled_AC_avgZscore)
all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore = all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore = all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(A_avg), fill=t.type.Under1Over1))+
  ggtitle("Superficial genes in VIS L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\nsuperficial gene counts") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_VIS_L23_glut_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_VIS_L23_glut_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_VIS_L23_glut_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(C_avg), fill=t.type.Under1Over1))+
  ggtitle("Deep genes in VIS L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\ndeep gene counts") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_VIS_L23_glut_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_VIS_L23_glut_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

#
exp1_sct_counts_VIS_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp1_sct_counts2, depth_table=all_VIS_L23_depths)
exp2_sct_counts_VIS_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp2_sct_counts2, depth_table=all_VIS_L23_depths)
exp6_sct_counts_VIS_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp6_sct_counts2, depth_table=all_VIS_L23_depths)
exp7_sct_counts_VIS_L23_glut_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp7_sct_counts2, depth_table=all_VIS_L23_depths)
all_sct_counts_VIS_L23_glut_AC_sctCount_ratio <- rbind(exp1_sct_counts_VIS_L23_glut_AC_sctCount_ratio, exp2_sct_counts_VIS_L23_glut_AC_sctCount_ratio, exp6_sct_counts_VIS_L23_glut_AC_sctCount_ratio, exp7_sct_counts_VIS_L23_glut_AC_sctCount_ratio)

all_sct_counts_VIS_L23_glut_AC_sctCount_ratio = all_sct_counts_VIS_L23_glut_AC_sctCount_ratio %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_VIS_L23_glut_AC_sctCount_ratio = all_sct_counts_VIS_L23_glut_AC_sctCount_ratio %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_VIS_L23_glut_AC_sctCount_ratio[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(AC_sctCount_ratio), fill=t.type.Under1Over1))+
  ggtitle("Superficial/deep gene ratio in VIS L2/3 MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,8))+
  ylab("Average SCTransform-corrected RNA count ratio\n(superficial/deep genes)") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_VIS_L23_glut_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/VIS/L23_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_VIS_L23_glut_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

###removing NAs before layer quintile calculation and then performing the rest of the above plots in V1 L2/3
#file with layer depth values annotated for each cell in V1 L2/3
all_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")

#AC ratio for each experiment with new layer quintiles in V1 L2/3
exp1_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp1_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp2_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp2_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp6_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp6_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp7_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- AC_avgZscore_func(exp_sct_count_table=exp7_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)

#
all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- rbind(exp1_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, exp2_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, exp6_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, exp7_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore)
all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore = all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore = all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(A_avg), fill=t.type.Under1Over1))+
  ggtitle("Superficial genes in V1 L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\nsuperficial gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(C_avg), fill=t.type.Under1Over1))+
  ggtitle("Deep genes in V1 L2/3, MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  #coord_cartesian(ylim=c(-0.3,0.3))+
  ylab("Average z-score of SCTransform-corrected\ndeep gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


exp1_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp1_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp2_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp2_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp6_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp6_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp7_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- AC_sctCount_ratio_func(exp_sct_count_table=exp7_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- rbind(exp1_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio, exp2_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio, exp6_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio, exp7_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio)

all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio = all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio = all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio$log2_AC_sctCount_ratio <- log2(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio$AC_sctCount_ratio)

ggplot(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(AC_sctCount_ratio), fill=t.type.Under1Over1))+
  ggtitle("Superficial/deep gene ratio in V1 L2/3 MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,8))+
  ylab("Average SCTransform-corrected RNA count ratio\n(superficial/deep genes)") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

#log2 version of above plot
ggplot(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = log2(as.numeric(AC_sctCount_ratio)), fill=t.type.Under1Over1))+
  ggtitle("Superficial/deep gene ratio in V1 L2/3 MeCP2 KO/+")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-3.5,3.5))+
  ylab("Log2 average SCTransform-corrected RNA count ratio\n(superficial/deep genes)") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_AC_log2_avg_sctCount_ratio_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_glut_noNA_Cheng2022_AC_log2_avg_sctCount_ratio_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


#wilcox.test(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[]

#presentation version
ggplot(all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(A_avg), fill=t.type.Under1Over1))+
  ggtitle("Superficial genes in glutamatergic cells\nof V1 L2/3")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.65,1.8))+
  ylab("Average expression z-score\nsuperficial genes") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/presentation_V1_L23_glut_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/presentation_V1_L23_glut_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(C_avg), fill=t.type.Under1Over1))+
  ggtitle("Deep genes in glutamatergic cells\nof V1 L2/3")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.65,1.8))+
  ylab("Average expression z-score\ndeep genes") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/presentation_V1_L23_glut_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/presentation_V1_L23_glut_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(AC_sctCount_ratio), fill=t.type.Under1Over1))+
  ggtitle("Superficial/deep gene ratio in glutamatergic cells\nof V1 L2/3")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(0,8))+
  ylab("Average gene expression ratio\n(superficial/deep genes)") + xlab("Layer depth quintile (superficial to deep)")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/presentation_V1_L23_glut_noNA_Cheng2022_AC_avg_sctCount_ratio_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')

###GABAergic average z-score
gaba_AC_avgZscore_func <- function(exp_sct_count_table, depth_table){
  #subset of cells in depth table
  exp_sct_count_table_dep <- exp_sct_count_table[index %in% depth_table$X,]
  #Glutamatergic cells
  exp_sct_count_table_dep_gaba <- exp_sct_count_table_dep[class_label=="GABAergic"]
  exp_sct_count_table_dep_gaba_genesOnly <- data.frame(exp_sct_count_table_dep_gaba[, 1:551], row.names="index")
  exp_sct_count_table_dep_gaba_scaled <- scale(exp_sct_count_table_dep_gaba_genesOnly)
  
  exp_sct_count_table_dep_gaba_scaled_dt <- data.frame(exp_sct_count_table_dep_gaba_scaled) %>% data.table(keep.rownames="index")
  exp_sct_count_table_dep_gaba_scaled_dt <- left_join(x=exp_sct_count_table_dep_gaba_scaled_dt, y=depth_table, by=c("index"="X"))
  
  exp_sct_count_table_dep_gaba_scaled_AC_avgZscore <- copy(exp_sct_count_table_dep_gaba_scaled_dt)
  exp_sct_count_table_dep_gaba_scaled_AC_avgZscore$A_avg <- rowMeans(exp_sct_count_table_dep_gaba_scaled_dt[, ..L23_A_merfish_newNames])
  exp_sct_count_table_dep_gaba_scaled_AC_avgZscore$C_avg <- rowMeans(exp_sct_count_table_dep_gaba_scaled_dt[, ..L23_C_merfish_newNames])
  return(exp_sct_count_table_dep_gaba_scaled_AC_avgZscore)
}

#
#AC ratio for each experiment with new layer quintiles in V1 L2/3
exp1_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore <- gaba_AC_avgZscore_func(exp_sct_count_table=exp1_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp2_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore <- gaba_AC_avgZscore_func(exp_sct_count_table=exp2_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp6_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore <- gaba_AC_avgZscore_func(exp_sct_count_table=exp6_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)
exp7_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore <- gaba_AC_avgZscore_func(exp_sct_count_table=exp7_sct_counts2_noNA, depth_table=all_V1_L23_depths_noNA)

#
all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore <- rbind(exp1_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore, exp2_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore, exp6_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore, exp7_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore)
all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore = all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore %>% mutate(t.type.Under1Over1 = factor(t.type.Under1Over1, levels=c("WT", "KO")))
all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore = all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore %>% mutate(depth_name = factor(depth_name, levels=c("Q1", "Q2", "Q3", "Q4", "Q5")))

ggplot(all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(A_avg), fill=t.type.Under1Over1))+
  ggtitle("Superficial genes in GABAergic cells\nof V1 L2/3")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.65,1.8))+
  ylab("Average z-score of SCTransform-corrected\nsuperficial gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_gaba_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_gaba_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

ggplot(all_sct_counts_V1_L23_gaba_noNA_scaled_AC_avgZscore[t.type.Under1Over1 %in% c("WT", "KO")], aes(x = t.type.Under1Over1, y = as.numeric(C_avg), fill=t.type.Under1Over1))+
  ggtitle("Deep genes in GABAergic cells\nof V1 L2/3")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE)+
  coord_cartesian(ylim=c(-0.65,1.8))+
  ylab("Average z-score of SCTransform-corrected\ndeep gene counts") + xlab("Layer depth quintile")+
  scale_fill_manual(name = "", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~depth_name,
             switch = "x",
             labeller = label_wrap_gen(width=15)
  ) + 
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(angle=0, size=10))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.y=element_text(size=10))
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_gaba_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/V1_L23_gaba_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')

#p-values
wilcox.test(all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[(t.type.Under1Over1=="WT") & (layer_quintile==1), C_avg], all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore[(t.type.Under1Over1=="KO") & (layer_quintile==1), C_avg])$p.value

wilcox.test(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="WT") & (layer_quintile==1), AC_sctCount_ratio], all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="KO") & (layer_quintile==1), AC_sctCount_ratio])$p.value

wilcox.test(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="WT") & (layer_quintile==1), log2(AC_sctCount_ratio)], all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="KO") & (layer_quintile==1), log2(AC_sctCount_ratio)])$p.value
wilcox.test(all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="WT") & (layer_quintile==1), log2_AC_sctCount_ratio], all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio[(t.type.Under1Over1=="KO") & (layer_quintile==1), log2_AC_sctCount_ratio])$p.value


exp_zScores_WT_vs_KO_layerDepth_sigs_func <- function(layerDepth_exp_table, value_column, Cheng_gene_type){
  sig1 <- pval1 <- NA
  sig2 <- pval2 <- NA
  sig3 <- pval3 <- NA
  sig4 <- pval4 <- NA
  sig5 <- pval5 <- NA
  
  #number of genes that go into Wilcoxon rank-sum tests
  n.WT.cells.layerQ1 <- length(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==1), as.numeric(get(value_column))])
  n.KO.cells.layerQ1 <- length(layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==1), as.numeric(get(value_column))])
  
  n.WT.cells.layerQ2 <- length(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==2), as.numeric(get(value_column))])
  n.KO.cells.layerQ2 <- length(layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==2), as.numeric(get(value_column))])
  
  n.WT.cells.layerQ3 <- length(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==3), as.numeric(get(value_column))])
  n.KO.cells.layerQ3 <- length(layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==3), as.numeric(get(value_column))])
  
  n.WT.cells.layerQ4 <- length(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==4), as.numeric(get(value_column))])
  n.KO.cells.layerQ4 <- length(layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==4), as.numeric(get(value_column))])
  
  n.WT.cells.layerQ5 <- length(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==5), as.numeric(get(value_column))])
  n.KO.cells.layerQ5 <- length(layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==5), as.numeric(get(value_column))])
  
  tryCatch({
    pval1 <- wilcox.test(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==1), as.numeric(get(value_column))], layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==1), as.numeric(get(value_column))])$p.value
    sig1 <- sig_function(pval1)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval2 <- wilcox.test(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==2), as.numeric(get(value_column))], layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==2), as.numeric(get(value_column))])$p.value
    sig2 <- sig_function(pval2)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval3 <- wilcox.test(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==3), as.numeric(get(value_column))], layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==3), as.numeric(get(value_column))])$p.value
    sig3 <- sig_function(pval3)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval4 <- wilcox.test(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==4), as.numeric(get(value_column))], layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==4), as.numeric(get(value_column))])$p.value
    sig4 <- sig_function(pval4)
  }, error = function(e) {
    # Handle the error
  })
  
  tryCatch({
    pval5 <- wilcox.test(layerDepth_exp_table[(t.type.Under1Over1=="WT") & (layer_quintile==5), as.numeric(get(value_column))], layerDepth_exp_table[(t.type.Under1Over1=="KO") & (layer_quintile==5), as.numeric(get(value_column))])$p.value
    sig5 <- sig_function(pval5)
  }, error = function(e) {
    # Handle the error
  })
  
  sigs <- data.table(rbind(cbind(Cheng_gene_type=Cheng_gene_type, layer_quintile=1, n.WT.cells = n.WT.cells.layerQ1, n.KO.cells = n.KO.cells.layerQ1, sig_symbol=sig1, wilcox.pval=pval1),
                           cbind(Cheng_gene_type=Cheng_gene_type, layer_quintile=2, n.WT.cells = n.WT.cells.layerQ2, n.KO.cells = n.KO.cells.layerQ2, sig_symbol=sig2, wilcox.pval=pval2),
                           cbind(Cheng_gene_type=Cheng_gene_type, layer_quintile=3, n.WT.cells = n.WT.cells.layerQ3, n.KO.cells = n.KO.cells.layerQ3, sig_symbol=sig3, wilcox.pval=pval3),
                           cbind(Cheng_gene_type=Cheng_gene_type, layer_quintile=4, n.WT.cells = n.WT.cells.layerQ4, n.KO.cells = n.KO.cells.layerQ4, sig_symbol=sig4, wilcox.pval=pval4),
                           cbind(Cheng_gene_type=Cheng_gene_type, layer_quintile=5, n.WT.cells = n.WT.cells.layerQ5, n.KO.cells = n.KO.cells.layerQ5, sig_symbol=sig5, wilcox.pval=pval5)))
  return(sigs)
}

A_sigs_all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- exp_zScores_WT_vs_KO_layerDepth_sigs_func(layerDepth_exp_table=all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, value_column="A_avg", Cheng_gene_type="A")
C_sigs_all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore <- exp_zScores_WT_vs_KO_layerDepth_sigs_func(layerDepth_exp_table=all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, value_column="C_avg", Cheng_gene_type="C")
AC_sigs_all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio <- exp_zScores_WT_vs_KO_layerDepth_sigs_func(layerDepth_exp_table=all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio, value_column="log2_AC_sctCount_ratio", Cheng_gene_type="A/C ratio")

write.csv(A_sigs_all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/wilcoxPvals_V1_L23_glut_noNA_Cheng2022_A_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.csv", row.names=F)
write.csv(C_sigs_all_sct_counts_V1_L23_glut_noNA_scaled_AC_avgZscore, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/wilcoxPvals_V1_L23_glut_noNA_Cheng2022_C_genes_avgZscore_sct_counts_by_layerQuintile_boxplot.csv", row.names=F)
write.csv(AC_sigs_all_sct_counts_V1_L23_glut_noNA_AC_sctCount_ratio, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/boxplots/V1/wilcoxPvals_V1_L23_glut_noNA_Cheng2022_AC_log2_avg_sctCount_ratio_by_layerQuintile_boxplot.csv", row.names=F)

#sigs_nCells10_L23_V1_logFC_TOI_A <- L23_vs_nonL23_sigs_func(logfc_table=L23_V1_logFC_TOI_A, geneList1="A", geneList2="non-A")