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
#library(fields)
options(scipy=999)
gene_panel_classes = fread("HG_lab/Vizgen/Gene_lists_for_analysis/MERFISH_gene_panel_classes.txt")
#gene overlap function, Fisher's test
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

gene_overlap_func3 <- function(full_geneList, geneList1, geneList2){
  geneList1 <- intersect(geneList1, full_geneList)
  geneList2 <- intersect(geneList2, full_geneList)
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
  #Fisher odds ratio
  OR = test_contingency_table$estimate
  return(list(p_values_contingency_table, OR, contingency_table))
}

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

nonL23_A_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_A)
nonL23_B_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_B)
nonL23_C_merfish <- setdiff(gene_panel_classes$Gene, L23_genes_C)

#formatting L23 gene names to fit data.frame standards
L23_genes_newNames <- make.names(L23_genes)
L23_A_merfish_newNames <- make.names(L23_A_merfish)
L23_B_merfish_newNames <- make.names(L23_B_merfish)
L23_C_merfish_newNames <- make.names(L23_C_merfish)

write.table(L23_A_merfish_newNames, file="HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.A.genes.merfishPanel.newNames.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(L23_B_merfish_newNames, file="HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.B.genes.merfishPanel.newNames.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(L23_C_merfish_newNames, file="HG_lab/Mati/GabelLab/genesets/Cheng2022/L2.3.C.genes.merfishPanel.newNames.txt", row.names=F, col.names=F, sep="\t", quote=F)

nonL23_A_merfish_newNames <- make.names(nonL23_A_merfish)
nonL23_B_merfish_newNames <- make.names(nonL23_B_merfish)
nonL23_C_merfish_newNames <- make.names(nonL23_C_merfish)
#
#annotation table from Yao et al, 2021
CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
#makes the cl column a character, useful for some data table joins later
CTX_HIP_annot[, cl := as.character(cl)]
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


exp1_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")
exp2_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp2_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")
exp6_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp6_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")
exp7_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp7_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")
all_V1_L23_depths_noNA=fread("HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/V1_layerdepth_files/L23/exp1_exp2_exp6_exp7_V1_L23_noNA_minVol100_maxVol3Med_minCount300_layerDepths_table.csv")

#make separate data tables for each hemisphere
exp1_V1L23.L <- exp1_V1_L23_depths_noNA[hemisphere == "L"]
exp2_V1L23.L <- exp2_V1_L23_depths_noNA[hemisphere == "L"]
exp6_V1L23.L <- exp6_V1_L23_depths_noNA[hemisphere == "L"]
exp7_V1L23.L <- exp7_V1_L23_depths_noNA[hemisphere == "L"]

exp1_V1L23.R <- exp1_V1_L23_depths_noNA[hemisphere == "R"]
exp2_V1L23.R <- exp2_V1_L23_depths_noNA[hemisphere == "R"]
exp6_V1L23.R <- exp6_V1_L23_depths_noNA[hemisphere == "R"]
exp7_V1L23.R <- exp7_V1_L23_depths_noNA[hemisphere == "R"]

#getting z-scores of sctransform-corrected RNA counts for glutamatergic cells in V1 L2/3 of experiment 2
exp2_sct_counts2_noNA_V1_L23 <- exp2_sct_counts2_noNA[index %in% exp2_V1_L23_depths_noNA$X,]
exp2_sct_counts2_noNA_V1_L23_glut <- exp2_sct_counts2_noNA_V1_L23 [class_label=="Glutamatergic",]
exp2_sct_counts2_noNA_V1_L23_glut_genesOnly <- data.frame(exp2_sct_counts2_noNA_V1_L23_glut[, 1:551], row.names="index")
exp2_sct_counts2_noNA_V1_L23_glut_scaled <- scale(exp2_sct_counts2_noNA_V1_L23_glut_genesOnly )
join_cols_exp2_sct_counts2_noNA_V1_L23_glut <- c("index", names(exp2_sct_counts2_noNA_V1_L23_glut)[552:573])

exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt <- data.frame(exp2_sct_counts2_noNA_V1_L23_glut_scaled) %>% data.table(keep.rownames="index")
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt <- inner_join(x=exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt, y=exp2_sct_counts2_noNA_V1_L23_glut[, ..join_cols_exp2_sct_counts2_noNA_V1_L23_glut], by="index")

A_pal3 <- colorRampPalette(c('lightblue','purple3'))
B_pal3 <- colorRampPalette(c('yellow', 'orangered1'))
C_pal3 <- colorRampPalette(c('yellow','darkgreen'))

limx_L23 <- c(2000,4200)
limy_L23 <- c(8000, 9200)


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

example_gene_plot_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_V1_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_V1_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6", dot_size=2)

example_gene_plot_zScoreOrder_eps_func(example_gene="Spon1", exp_data_table=exp2_sct_counts_V1_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="A_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/zScoreOrder_cex2_exp2_V1_L23_L_ZscoreAcrossL23_scatterplot_Cheng_A_gene_Spon1", dot_size=2)


###ordering cell plotting by expression
#example_gene_plot_eps_func(example_gene="Trpc6", exp_data_table=exp2_sct_counts_V1_glut_L23_scaled_dt, cell.id=exp2_L23_L$index, color_pal="B_pal3", color_column="plot_color", xlim=limx_L23, ylim=limy_L23, output_file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/cex2_exp2_L23_L_ZscoreAcrossL23_scatterplot_Cheng_B_gene_Trpc6", dot_size=2)

#exp2_zScoreOrder_Trpc6 = exp2_sct_counts_V1_glut_L23_scaled_dt[order(Trpc6), index]
exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Cdh13  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="WT")][(order(Cdh13)), index]
exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Cdh13  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="KO")][(order(Cdh13)), index]


exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Tox  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="WT")][(order(Tox)), index]
exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Tox  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="KO")][(order(Tox)), index]

exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Bmper  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="WT")][(order(Bmper)), index]
exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Bmper  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="KO")][(order(Bmper)), index]

exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Epha6  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="WT")][(order(Epha6)), index]
exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Epha6  = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp2_V1L23.L$X) & (t.type.Under1Over1=="KO")][(order(Epha6)), index]



exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13 <- copy(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox <- copy(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6 <- copy(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper <- copy(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt)

exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$plot_color <- get("A_pal3")(20)[as.numeric(cut(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")], breaks=20))]
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Tox")], breaks=20))]
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6$plot_color <- get("A_pal3")(20)[as.numeric(cut(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Epha6")], breaks=20))]
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Bmper")], breaks=20))]


#Cdh13 spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Cdh13_WT.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L WT", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Cdh13){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Cdh13_KO.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L KO", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Cdh13){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,plot_color])
}
dev.off()

#Tox spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Tox_WT.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L WT", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Tox){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Tox_KO.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L KO", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Tox){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,plot_color])
}
dev.off()

#
#Bmper spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Bmper_WT.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L WT Bmper", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Bmper){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Bmper_KO.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L KO Bmper", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Bmper){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,plot_color])
}
dev.off()

#Epha6 spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Epha6_WT.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L WT Epha6", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_WT_Epha6){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp2_V1_L23_glut_L_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Epha6_KO.png", width=2000, height=2000, res=300)
plot(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp2_V1L23.L$X),center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp2_V1L23.L$X),center_y],pch = 20,cex=.3,col = "white",xlim=limx_L23, ylim=limy_L23,main="Exp2 V1 L2/3 L KO Epha6", xlab="x", ylab="y", asp=1)
for(i in exp2_zScoreOrder_noNA_V1_L23_glut_L_KO_Epha6){
  points(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_x],exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_y],pch = 20,cex=2, col = exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,plot_color])
}
dev.off()




color.bar(get("A_pal3")(20), min=signif(min(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")]),2), max=signif(max(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")]),2), nticks=2, title="Cdh13")
#color.bar(get("B_pal3")(20), min=signif(min(exp2_sct_counts_V1_glut_L23_scaled_dt_Trpc6[, get("Trpc6")]),2), max=signif(max(exp2_sct_counts_V1_glut_L23_scaled_dt_Trpc6[, get("Trpc6")]),2), nticks=2, title="Trpc6")
color.bar(get("C_pal3")(20), min=signif(min(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Tox")]),2), max=signif(max(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Tox")]),2), nticks=2, title="Tox")

#color bars for Epha6 and Bmper
color.bar(get("A_pal3")(20), min=signif(min(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[, get("Epha6")]),2), max=signif(max(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Epha6")]),2), nticks=2, title="Epha6")
color.bar(get("C_pal3")(20), min=signif(min(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Bmper")]),2), max=signif(max(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Bmper")]),2), nticks=2, title="Bmper")



#####cells in L2/3 only

L23_V1_WT_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.WT.csv")
L23_V1_KO_avgExp_matrix <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/expression_tables/L23.V1.matrix.averageExpression.sctransform.subtypes.mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.KO.csv")

L23_V1_WT_and_KO_avgExp_matrix_new_TOI <- left_join(x=L23_V1_WT_avgExp_matrix[, c("Gene","164_L2/3 IT CTX", "166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX", "180_L4 IT CTX")],
                                                y=L23_V1_KO_avgExp_matrix[, c("Gene","164_L2/3 IT CTX", "166_L2/3 IT CTX","171_L2/3 IT CTX","168_L2/3 IT CTX","179_L4 IT CTX", "180_L4 IT CTX")],
                                                by=c("Gene"), suffix = c(".WT", ".KO"))



L23_V1_WT_and_KO_avgExp_matrix_new_TOI <- data.frame(L23_V1_WT_and_KO_avgExp_matrix_new_TOI[Gene %in% L23_genes_newNames], row.names="Gene")

L23_V1_WT_and_KO_avgExp_matrix_new_TOI_ABC_zscore <- t(scale(t(as.matrix(L23_V1_WT_and_KO_avgExp_matrix_new_TOI))))


L23_V1_palette <- colorRampPalette(c("blue", "white",  "red" ))(n = 299)


new_type_order <- c("X164_L2.3.IT.CTX.WT", "X164_L2.3.IT.CTX.KO",
                "X166_L2.3.IT.CTX.WT", "X166_L2.3.IT.CTX.KO",
                "X171_L2.3.IT.CTX.WT", "X171_L2.3.IT.CTX.KO",
                "X168_L2.3.IT.CTX.WT", "X168_L2.3.IT.CTX.KO",
                "X179_L4.IT.CTX.WT", "X179_L4.IT.CTX.KO",
                "X180_L4.IT.CTX.WT", "X180_L4.IT.CTX.KO")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/L23_V1_L23_ABC_164_166_171_168_179_180_AverageExpression_zscores_merfish_heatmap.eps")
heatmap.2(L23_V1_WT_and_KO_avgExp_matrix_new_TOI_ABC_zscore[c(L23_A_merfish_newNames, L23_B_merfish_newNames, L23_C_merfish_newNames), new_type_order],
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.2, 
          cexCol = 0.6, labRow = NULL, labCol = NULL, col=L23_V1_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()


#pseudobulkDGE of glutamatergic cells in V1 L2/3
V1_L23_glut_layerDepth_logFC =fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/pseudobulkDGE.glut.cells.in.V1.L23.layerdepth.5quantiles.noNA.de.results.csv")



#layer quantile significantly MR gene overlap with L2/3 A/B/C genes
V1_L23_glut_Q1_MR_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q1_MR_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q1_MR_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q2_MR_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q2_MR_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q2_MR_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q3_MR_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q3_MR_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q3_MR_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q4_MR_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q4_MR_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q4_MR_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q5_MR_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q5_MR_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q5_MR_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC > 0), gene], L23_C_merfish_newNames)

#layer quantile significantly MA gene overlap with L2/3 A/B/C genes
V1_L23_glut_Q1_MA_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q1_MA_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q1_MA_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==1) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q2_MA_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q2_MA_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q2_MA_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==2) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q3_MA_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q3_MA_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q3_MA_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==3) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q4_MA_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q4_MA_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q4_MA_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==4) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

V1_L23_glut_Q5_MA_L23_A <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_A_merfish_newNames)
V1_L23_glut_Q5_MA_L23_B <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_B_merfish_newNames)
V1_L23_glut_Q5_MA_L23_C <-  gene_overlap_func3(gene_panel_classes[,Gene], V1_L23_glut_layerDepth_logFC[(layer_quantile==5) & (p.adj.BH.all < 0.05) & (logFC < 0), gene], L23_C_merfish_newNames)

#matrix of p-values of overlaps between statistically Mecp2-regulated genes and L2/3 A/B/C genes
V1_L23_glut_layerQuantiles_sig_ABC_pvals_matrix <- matrix(c(unlist(V1_L23_glut_Q1_MR_L23_A[1]), unlist(V1_L23_glut_Q1_MR_L23_B[1]), unlist(V1_L23_glut_Q1_MR_L23_C[1]),
                                                            unlist(V1_L23_glut_Q2_MR_L23_A[1]), unlist(V1_L23_glut_Q2_MR_L23_B[1]), unlist(V1_L23_glut_Q2_MR_L23_C[1]),
                                                            unlist(V1_L23_glut_Q3_MR_L23_A[1]), unlist(V1_L23_glut_Q3_MR_L23_B[1]), unlist(V1_L23_glut_Q3_MR_L23_C[1]),
                                                            unlist(V1_L23_glut_Q4_MR_L23_A[1]), unlist(V1_L23_glut_Q4_MR_L23_B[1]), unlist(V1_L23_glut_Q4_MR_L23_C[1]),
                                                            unlist(V1_L23_glut_Q5_MR_L23_A[1]), unlist(V1_L23_glut_Q5_MR_L23_B[1]), unlist(V1_L23_glut_Q5_MR_L23_C[1]),
                                                            unlist(V1_L23_glut_Q1_MA_L23_A[1]), unlist(V1_L23_glut_Q1_MA_L23_B[1]), unlist(V1_L23_glut_Q1_MA_L23_C[1]),
                                                            unlist(V1_L23_glut_Q2_MA_L23_A[1]), unlist(V1_L23_glut_Q2_MA_L23_B[1]), unlist(V1_L23_glut_Q2_MA_L23_C[1]),
                                                            unlist(V1_L23_glut_Q3_MA_L23_A[1]), unlist(V1_L23_glut_Q3_MA_L23_B[1]), unlist(V1_L23_glut_Q3_MA_L23_C[1]),
                                                            unlist(V1_L23_glut_Q4_MA_L23_A[1]), unlist(V1_L23_glut_Q4_MA_L23_B[1]), unlist(V1_L23_glut_Q4_MA_L23_C[1]),
                                                            unlist(V1_L23_glut_Q5_MA_L23_A[1]), unlist(V1_L23_glut_Q5_MA_L23_B[1]), unlist(V1_L23_glut_Q5_MA_L23_C[1])), 
                                              nrow=3, ncol=10)
rownames(V1_L23_glut_layerQuantiles_sig_ABC_pvals_matrix) = c("A", "B", "C")
colnames(V1_L23_glut_layerQuantiles_sig_ABC_pvals_matrix) = c("MR_Q1", "MR_Q2", "MR_Q3", "MR_Q4", "MR_Q5",
                                                              "MA_Q1", "MA_Q2", "MA_Q3", "MA_Q4", "MA_Q5")

#write.table(nCells5_L23_V1_sig_ABC_pvals_matrix, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells5_L23_V1_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_merfish_pvals.txt", row.names=TRUE, col.names=TRUE, quote=F, sep="\t")
V1_L23_glut_layerQuantiles_sig_ABC_log10pvals_matrix<- -log10(V1_L23_glut_layerQuantiles_sig_ABC_pvals_matrix + 1 - 1)
write.table(V1_L23_glut_layerQuantiles_sig_ABC_log10pvals_matrix, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells10_glut_cells_in_V1_L23_Cheng_correctOverlaps_stat_sig_MeCP2reg_genes_by_layerDepth_merfish_log10pvals.txt", row.names=TRUE, col.names=TRUE, quote=F, sep="\t")
#log2 odds ratio
V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix <- matrix(c(unlist(V1_L23_glut_Q1_MR_L23_A[2]), unlist(V1_L23_glut_Q1_MR_L23_B[2]), unlist(V1_L23_glut_Q1_MR_L23_C[2]),
                                                            unlist(V1_L23_glut_Q2_MR_L23_A[2]), unlist(V1_L23_glut_Q2_MR_L23_B[2]), unlist(V1_L23_glut_Q2_MR_L23_C[2]),
                                                            unlist(V1_L23_glut_Q3_MR_L23_A[2]), unlist(V1_L23_glut_Q3_MR_L23_B[2]), unlist(V1_L23_glut_Q3_MR_L23_C[2]),
                                                            unlist(V1_L23_glut_Q4_MR_L23_A[2]), unlist(V1_L23_glut_Q4_MR_L23_B[2]), unlist(V1_L23_glut_Q4_MR_L23_C[2]),
                                                            unlist(V1_L23_glut_Q5_MR_L23_A[2]), unlist(V1_L23_glut_Q5_MR_L23_B[2]), unlist(V1_L23_glut_Q5_MR_L23_C[2]),
                                                            unlist(V1_L23_glut_Q1_MA_L23_A[2]), unlist(V1_L23_glut_Q1_MA_L23_B[2]), unlist(V1_L23_glut_Q1_MA_L23_C[2]),
                                                            unlist(V1_L23_glut_Q2_MA_L23_A[2]), unlist(V1_L23_glut_Q2_MA_L23_B[2]), unlist(V1_L23_glut_Q2_MA_L23_C[2]),
                                                            unlist(V1_L23_glut_Q3_MA_L23_A[2]), unlist(V1_L23_glut_Q3_MA_L23_B[2]), unlist(V1_L23_glut_Q3_MA_L23_C[2]),
                                                            unlist(V1_L23_glut_Q4_MA_L23_A[2]), unlist(V1_L23_glut_Q4_MA_L23_B[2]), unlist(V1_L23_glut_Q4_MA_L23_C[2]),
                                                            unlist(V1_L23_glut_Q5_MA_L23_A[2]), unlist(V1_L23_glut_Q5_MA_L23_B[2]), unlist(V1_L23_glut_Q5_MA_L23_C[2])), 
                                                          nrow=3, ncol=10)
V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix <- log2(V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix)
rownames(V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix) = c("A", "B", "C")
colnames(V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix) = c("MR_Q1", "MR_Q2", "MR_Q3", "MR_Q4", "MR_Q5",
                                                              "MA_Q1", "MA_Q2", "MA_Q3", "MA_Q4", "MA_Q5")

round(V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix,2)
write.table(round(V1_L23_glut_layerQuantiles_sig_ABC_log2OR_matrix,2), file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells10_glut_cells_in_V1_L23_Cheng_correctOverlaps_stat_sig_MeCP2reg_genes_by_layerDepth_merfish_log2OR.txt", row.names=TRUE, col.names=TRUE, quote=F, sep="\t")

col_palette2 <- colorRampPalette(c("white", "red", "darkred"))(n = 299)
#col_breaks = c(seq(-3, -1,length=100),  # for blue
#               seq(-0.99,0.99,length=100), #for white
#               seq(1, 3, length=100)) #for red



setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/nCells10_glut_cells_in_V1_L23_ABC_Cheng_overlaps_stat_sig_MeCP2reg_genes_by_layerDepth_merfish_heatmap.eps")
heatmap.2(V1_L23_glut_layerQuantiles_sig_ABC_log10pvals_matrix,
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







###
#getting z-scores of sctransform-corrected RNA counts for glutamatergic cells in V1 L2/3 of experiment 1
exp1_sct_counts2_noNA_V1_L23 <- exp1_sct_counts2_noNA[index %in% exp1_V1_L23_depths_noNA$X,]
exp1_sct_counts2_noNA_V1_L23_glut <- exp1_sct_counts2_noNA_V1_L23 [class_label=="Glutamatergic",]
exp1_sct_counts2_noNA_V1_L23_glut_genesOnly <- data.frame(exp1_sct_counts2_noNA_V1_L23_glut[, 1:551], row.names="index")
exp1_sct_counts2_noNA_V1_L23_glut_scaled <- scale(exp1_sct_counts2_noNA_V1_L23_glut_genesOnly )
join_cols_exp1_sct_counts2_noNA_V1_L23_glut <- c("index", names(exp1_sct_counts2_noNA_V1_L23_glut)[552:573])

exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt <- data.frame(exp1_sct_counts2_noNA_V1_L23_glut_scaled) %>% data.table(keep.rownames="index")
exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt <- inner_join(x=exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt, y=exp1_sct_counts2_noNA_V1_L23_glut[, ..join_cols_exp1_sct_counts2_noNA_V1_L23_glut], by="index")

#getting z-scores of sctransform-corrected RNA counts for glutamatergic cells in V1 L2/3 of experiment 6
exp6_sct_counts2_noNA_V1_L23 <- exp6_sct_counts2_noNA[index %in% exp6_V1_L23_depths_noNA$X,]
exp6_sct_counts2_noNA_V1_L23_glut <- exp6_sct_counts2_noNA_V1_L23 [class_label=="Glutamatergic",]
exp6_sct_counts2_noNA_V1_L23_glut_genesOnly <- data.frame(exp6_sct_counts2_noNA_V1_L23_glut[, 1:551], row.names="index")
exp6_sct_counts2_noNA_V1_L23_glut_scaled <- scale(exp6_sct_counts2_noNA_V1_L23_glut_genesOnly )
join_cols_exp6_sct_counts2_noNA_V1_L23_glut <- c("index", names(exp6_sct_counts2_noNA_V1_L23_glut)[552:573])

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt <- data.frame(exp6_sct_counts2_noNA_V1_L23_glut_scaled) %>% data.table(keep.rownames="index")
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt <- inner_join(x=exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt, y=exp6_sct_counts2_noNA_V1_L23_glut[, ..join_cols_exp6_sct_counts2_noNA_V1_L23_glut], by="index")

#getting z-scores of sctransform-corrected RNA counts for glutamatergic cells in V1 L2/3 of experiment 7
exp7_sct_counts2_noNA_V1_L23 <- exp7_sct_counts2_noNA[index %in% exp7_V1_L23_depths_noNA$X,]
exp7_sct_counts2_noNA_V1_L23_glut <- exp7_sct_counts2_noNA_V1_L23 [class_label=="Glutamatergic",]
exp7_sct_counts2_noNA_V1_L23_glut_genesOnly <- data.frame(exp7_sct_counts2_noNA_V1_L23_glut[, 1:551], row.names="index")
exp7_sct_counts2_noNA_V1_L23_glut_scaled <- scale(exp7_sct_counts2_noNA_V1_L23_glut_genesOnly )
join_cols_exp7_sct_counts2_noNA_V1_L23_glut <- c("index", names(exp7_sct_counts2_noNA_V1_L23_glut)[552:573])

exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt <- data.frame(exp7_sct_counts2_noNA_V1_L23_glut_scaled) %>% data.table(keep.rownames="index")
exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt <- inner_join(x=exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt, y=exp7_sct_counts2_noNA_V1_L23_glut[, ..join_cols_exp7_sct_counts2_noNA_V1_L23_glut], by="index")



exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Cdh13  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="WT")][(order(Cdh13)), index]
exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Cdh13  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="KO")][(order(Cdh13)), index]


exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Tox  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="WT")][(order(Tox)), index]
exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Tox  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="KO")][(order(Tox)), index]

exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Bmper  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="WT")][(order(Bmper)), index]
exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Bmper  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="KO")][(order(Bmper)), index]

exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Epha6  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="WT")][(order(Epha6)), index]
exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Epha6  = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[(index %in% exp6_V1L23.R$X) & (t.type.Under1Over1=="KO")][(order(Epha6)), index]



exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13 <- copy(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox <- copy(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6 <- copy(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper <- copy(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt)

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$plot_color <- get("A_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")], breaks=20))]
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Tox")], breaks=20))]
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6$plot_color <- get("A_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[, get("Epha6")], breaks=20))]
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[, get("Bmper")], breaks=20))]


#Cdh13 spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Cdh13_WT.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white", xlim=c(0, 2200), ylim=c(3200,5400), main="Exp6 V1 L2/3 R WT Cdh13", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Cdh13){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Cdh13_KO.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R KO Cdh13", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Cdh13){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[index==i,plot_color])
}
dev.off()

#Tox spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Tox_WT.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R WT Tox", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Tox){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Tox_KO.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R KO Tox", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Tox){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[index==i,plot_color])
}
dev.off()

#
#Bmper spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Bmper_WT.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R WT Bmper", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Bmper){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_C_gene_Bmper_KO.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R KO Bmper", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Bmper){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[index==i,plot_color])
}
dev.off()

#Epha6 spatial scatter plots
png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Epha6_WT.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R WT Epha6", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_WT_Epha6){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,plot_color])
}
dev.off()

png("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/scatterplots/Cheng_spatialMaps/example_genes/WT_KO_zScoring_only/exp6_V1_L23_glut_R_sctCount_ZscoreAcrossV1L23_scatterplot_Cheng_A_gene_Epha6_KO.png", width=2000, height=2000, res=300)
plot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp6_V1L23.R$X),center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[(index %in% exp6_V1L23.R$X),center_y],pch = 20,cex=.3,col = "white",xlim=c(0, 2200), ylim=c(3200,5400),main="Exp6 V1 L2/3 R KO Epha6", xlab="x", ylab="y", asp=1)
for(i in exp6_zScoreOrder_noNA_V1_L23_glut_R_KO_Epha6){
  points(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_x],exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,center_y],pch = 20,cex=2, col = exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[index==i,plot_color])
}
dev.off()

#Tox and Cdh13
color.bar(get("A_pal3")(20), min=signif(min(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")]),2), max=signif(max(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")]),2), nticks=2, title="Cdh13")
color.bar(get("C_pal3")(20), min=signif(min(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Tox")]),2), max=signif(max(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Tox")]),2), nticks=2, title="Tox")

#color bars for Epha6 and Bmper
color.bar(get("A_pal3")(20), min=signif(min(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Epha6[, get("Epha6")]),2), max=signif(max(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Epha6")]),2), nticks=2, title="Epha6")
color.bar(get("C_pal3")(20), min=signif(min(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Bmper")]),2), max=signif(max(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Tox[, get("Bmper")]),2), nticks=2, title="Bmper")

####average expression by sublayer depth
#Seurat object 
mecp2Het.CTX.HC.seu=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_allTaxLabels_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

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


#subset Seurat object to only superficial L2/3 IT cells
mecp2Het.CTX.HC.seu.V1.L23.glut <- subset(mecp2Het.CTX.HC.seu, cells=All_V1.L23_layerdepth_glut$cell.id)

#turn relevant layer metadata into data frames with cell ids as row names
All_V1.L23_layerdepth_glut_meta_df <- data.frame(All_V1.L23_layerdepth_glut[, .(cell.id, norm_depth, layer_quantile)], row.names="cell.id")

#add layer metadata
mecp2Het.CTX.HC.seu.V1.L23.glut = AddMetaData(object=mecp2Het.CTX.HC.seu.V1.L23.glut, metadata=All_V1.L23_layerdepth_glut_meta_df)

#WT subset of Seurat object
mecp2Het.CTX.HC.seu.V1.L23.glut.WT = subset(mecp2Het.CTX.HC.seu.V1.L23.glut, subset = t.type == "WT")
#KO susbset of Seurat object
mecp2Het.CTX.HC.seu.V1.L23.glut.KO = subset(mecp2Het.CTX.HC.seu.V1.L23.glut, subset = t.type == "KO")

#average expression of WT and KO glutamatergic cells in V1 L2/3 by sublayer quantile
mecp2Het.CTX.HC.seu.V1.L23.glut.WT.avgExp.per.layerDepth <- data.table(AverageExpression(obj=mecp2Het.CTX.HC.seu.V1.L23.glut.WT, assays="SCT", slot="data", group.by = c("layer_quantile"))$SCT, keep.rownames="Gene")
mecp2Het.CTX.HC.seu.V1.L23.glut.KO.avgExp.per.layerDepth <- data.table(AverageExpression(obj=mecp2Het.CTX.HC.seu.V1.L23.glut.KO, assays="SCT", slot="data", group.by = c("layer_quantile"))$SCT, keep.rownames="Gene")


L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth <- left_join(x=mecp2Het.CTX.HC.seu.V1.L23.glut.WT.avgExp.per.layerDepth, y=mecp2Het.CTX.HC.seu.V1.L23.glut.KO.avgExp.per.layerDepth,
                                                               by=c("Gene"), suffix = c(".WT", ".KO"))



L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df <- data.frame(L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth, row.names="Gene")

L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df_zscore <- t(scale(t(as.matrix(L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df))))


L23_V1_palette <- colorRampPalette(c("blue", "white",  "red" ))(n = 299)

#order for columns for heatmap of scaled average gene expression heatmap by sublayer quantile
#quantile_ttype_order <- c("X1.WT", "X2.WT",
#                    "X166_L2.3.IT.CTX.WT", "X166_L2.3.IT.CTX.KO",
 #                   "X171_L2.3.IT.CTX.WT", "X171_L2.3.IT.CTX.KO",
  #                  "X168_L2.3.IT.CTX.WT", "X168_L2.3.IT.CTX.KO",
   #                 "X179_L4.IT.CTX.WT", "X179_L4.IT.CTX.KO",
    #                "X180_L4.IT.CTX.WT", "X180_L4.IT.CTX.KO")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/glut_cells_in_V1_L23_AverageExpression_zscores_by_layerDepth_Cdh13_Bmper_merfish_heatmap.eps")
heatmap.2(L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df_zscore[c("Cdh13", "Bmper"),],
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, labRow = NULL, labCol = NULL, col=L23_V1_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

###using raw RNA counts fo average gene expression
#average expression of WT and KO glutamatergic cells in V1 L2/3 by sublayer quantile
raw.mecp2Het.CTX.HC.seu.V1.L23.glut.WT.avgExp.per.layerDepth <- data.table(AverageExpression(obj=mecp2Het.CTX.HC.seu.V1.L23.glut.WT, assays="RNA", slot="data", group.by = c("layer_quantile"))$RNA, keep.rownames="Gene")
raw.mecp2Het.CTX.HC.seu.V1.L23.glut.KO.avgExp.per.layerDepth <- data.table(AverageExpression(obj=mecp2Het.CTX.HC.seu.V1.L23.glut.KO, assays="RNA", slot="data", group.by = c("layer_quantile"))$RNA, keep.rownames="Gene")


raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth <- left_join(x=raw.mecp2Het.CTX.HC.seu.V1.L23.glut.WT.avgExp.per.layerDepth, y=raw.mecp2Het.CTX.HC.seu.V1.L23.glut.KO.avgExp.per.layerDepth,
                                                               by=c("Gene"), suffix = c(".WT", ".KO"))



raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df <- data.frame(raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth, row.names="Gene")

raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df_zscore <- t(scale(t(as.matrix(raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df))))


L23_V1_palette <- colorRampPalette(c("blue", "white",  "red" ))(n = 299)

#order for columns for heatmap of scaled average gene expression heatmap by sublayer quantile
#quantile_ttype_order <- c("X1.WT", "X2.WT",
#                    "X166_L2.3.IT.CTX.WT", "X166_L2.3.IT.CTX.KO",
#                   "X171_L2.3.IT.CTX.WT", "X171_L2.3.IT.CTX.KO",
#                  "X168_L2.3.IT.CTX.WT", "X168_L2.3.IT.CTX.KO",
#                 "X179_L4.IT.CTX.WT", "X179_L4.IT.CTX.KO",
#                "X180_L4.IT.CTX.WT", "X180_L4.IT.CTX.KO")

setEPS()
postscript("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/glut_cells_in_V1_L23_rawCounts_AverageExpression_zscores_by_layerDepth_Cdh13_Bmper_merfish_heatmap.eps")
heatmap.2(raw_L23_V1_glut_WT_and_KO_avgExp_matrix_by_layerDepth_df_zscore[c("Cdh13", "Bmper"),],
          notecol = "black", offsetRow=0.1, offsetCol=0.1, density.info = "none", 
          trace ="none", margins = c(5,10), ColSideColors, RowSideColors, cexRow = 0.6, 
          cexCol = 0.6, labRow = NULL, labCol = NULL, col=L23_V1_palette, dendrogram="none", 
          Rowv="NA", Colv="NA", key=TRUE, key.title = NA, key.xlab = NA, key.ylab = NA, keysize = 0.8,
          srtCol=25, lmat = rbind(c(0, 3), c(2,1), c(0,4)), lhei = c(0.85, 3.5, 0.5))
dev.off()

head(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper)

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths <- left_join(x=exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper, All_V1.L23_layerdepth_glut[, .(cell.id, layer_quantile)], by=c("index"="cell.id"))
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths <- left_join(x=exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13, All_V1.L23_layerdepth_glut[, .(cell.id, layer_quantile)], by=c("index"="cell.id"))


#exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore <- exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths[, .(avg_Bmper = mean(Bmper)), by = .(t.type.Under1Over1, layer_quantile)]
#exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore <- exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths[, .(avg_Chd13 = mean(Cdh13)), by = .(t.type.Under1Over1, layer_quantile)]

layer_order <- c(1, 2, 3, 4, 5)
# Set the key columns for sorting
setkey(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths, t.type.Under1Over1, layer_quantile)
setkey(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths, t.type.Under1Over1, layer_quantile)

# Reorder the data.table based on the custom order
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore <- exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths[, .(avg_Bmper = mean(Bmper)), by = .(t.type.Under1Over1, layer_quantile)]
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore <- exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths[, .(avg_Cdh13 = mean(Cdh13)), by = .(t.type.Under1Over1, layer_quantile)]

exp6_Bmper_palette <- get("C_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, breaks = 20))]
exp6_Cdh13_palette <- get("A_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$Cdh13, breaks = 20))]

#map average z-score values onto existing color palettes for Bmper and Cdh13 z-scores
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_color <- color_palette[
  findInterval(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, vec = quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, probs = seq(0, 1, length.out = 20 + 1)))
]

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_color <- exp6_Cdh13_palette[
  findInterval(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_Cdh13, vec = quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$Cdh13, probs = seq(0, 1, length.out = 20 + 1)))
]

breaks <- 20
break_values <- quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, probs = seq(0, 1, length.out = breaks + 1))

# Create a color mapping based on unique 'Bmper' values
Bmper_color_mapping <- setNames(exp6_Bmper_palette[findInterval(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, vec = break_values)], exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper)

# Map colors to 'avg_Bmper' in the second data.table using the color mapping
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_color <- Bmper_color_mapping[as.character(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper)]



#exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_matrix <- dcast(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths, layer_quantile ~ t.type.Under1Over1, value.var = "avg_Bmper")
#exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_matrix <- dcast(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths, layer_quantile ~ t.type.Under1Over1, value.var = "avg_Cdh13")

#exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$plot_color <- get("C_pal3")(20)[as.numeric(cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[, get("Bmper")], breaks=20))]

##
# Create a color mapping based on all 'Bmper' values in the first data.table
Bmper_color_mapping <- setNames(exp6_Bmper_palette[findInterval(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, vec = break_values)], exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper)

# Discretize avg_Bmper to match breaks and use it as a key for indexing the color mapping
discretized_avg_Bmper <- cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, breaks = break_values)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_color <- Bmper_color_mapping[as.numeric(discretized_avg_Bmper)]



setorder(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths, match(layer_quantile, layer_order))
setorder(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths, match(layer_quantile, layer_order))


# Use dcast to convert the data.table to a matrix
result_matrix <- dcast(result, t.type.Under1Over1 ~ layer_quantile, value.var = "avg_Bmper")


###
color_palette <- get("C_pal3")(20)
unique_values <- unique(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper)
breaks <- quantile(unique_values, probs = seq(0, 1, length.out = length(color_palette) + 1))

# Function to map numeric values to colors based on the color palette and breaks
map_to_color <- function(values, palette, breaks) {
  color_mapping <- setNames(palette, cut(values, breaks = breaks, include.lowest = TRUE))
  return(color_mapping[as.character(cut(values, breaks = breaks, include.lowest = TRUE))])
}

# Map colors to 'avg_Bmper' in the second data.table using the color mapping function
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_color <- map_to_color(
  exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper,
  color_palette,
  breaks
)


###
intervals_Bmper <- cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[, get("Bmper")], breaks = 20)
intervals_Cdh13 <- cut(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13[, get("Cdh13")], breaks = 20)

hist_values_Bmper <- hist(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper[, get("Bmper")], breaks = 20, plot = FALSE)
breaks_used_Bmper <- hist_values_Bmper$breaks

# Find intervals for 'avg_Bmper' in the second data.table using the same breaks
interval_avg_Bmper <- findInterval(
  exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper,
  breaks_used_Bmper
)

# Map colors based on the intervals
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_color <- get("C_pal3")(20)[interval_avg_Bmper]

##
find_interval_info_vector <- function(values, intervals) {
  result_list <- list()
  
  for (val in values) {
    for (j in seq_along(intervals)) {
      interval_values <- strsplit(gsub("\\[|\\(|\\]|\\)", "", intervals[j]), ",")[[1]]
      lower_bound <- as.numeric(interval_values[1])
      upper_bound <- as.numeric(interval_values[2])
      
      if (val > lower_bound && val <= upper_bound) {
        result_list[[length(result_list) + 1]] <- list(avg_val = val, index = j, interval = intervals[j])
        break  # Move to the next value
      }
    }
  }
  
  return(result_list)
}

results_Bmper <- find_interval_info_vector(values=exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, intervals=levels(intervals_Bmper))
results_Cdh13 <- find_interval_info_vector(values=exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_Cdh13, intervals=levels(intervals_Cdh13))

#results <- find_interval_info_vector(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, intervals_Bmper)

# Convert the results to a data.table
result_Bmper_dt <- data.table(
  avg_val = sapply(results_Bmper, function(x) x$avg_val),
  index = sapply(results_Bmper, function(x) x$index),
  interval = sapply(results_Bmper, function(x) x$interval)
)

result_Cdh13_dt <- data.table(
  avg_val = sapply(results_Cdh13, function(x) x$avg_val),
  index = sapply(results_Cdh13, function(x) x$index),
  interval = sapply(results_Cdh13, function(x) x$interval)
)


exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore <- cbind(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore, result_Bmper_dt)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper_color <- get("C_pal3")(20)[exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$index]

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore <- cbind(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore, result_Cdh13_dt)
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_Cdh13_color <- get("A_pal3")(20)[exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$index]

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_matrix <- dcast(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore, layer_quantile ~ t.type.Under1Over1, value.var = "avg_Bmper")
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_matrix <- dcast(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore, layer_quantile ~ t.type.Under1Over1, value.var = "avg_Cdh13")


ggplot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore, aes(x = layer_quantile, y = t.type.Under1Over1, fill = avg_Bmper_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.text=element_text(size=18), axis.ticks = element_blank())+
  labs(title = "Bmper", x = "Layer Quantile", y = "")
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/Bmper_Cheng_C_gene_V1_L23_glut_noNA_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/Bmper_Cheng_C_gene_V1_L23_glut_noNA_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')
  

ggplot(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore, aes(x = layer_quantile, y = t.type.Under1Over1, fill = avg_Cdh13_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.text=element_text(size=18), axis.ticks = element_blank())+
  labs(title = "Cdh13", x = "Layer Quantile", y = "")
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/Cdh13_Cheng_C_gene_V1_L23_glut_noNA_avgZscore_sct_counts_by_layerQuintile_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/Cdh13_Cheng_C_gene_V1_L23_glut_noNA_avgZscore_sct_counts_by_layerQuintile_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, seq(0, 1, 0.1))
quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, seq(0, 1, 0.1))

quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$Cdh13, seq(0, 1, 0.1))
quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_Cdh13, seq(0, 1, 0.1))

zScore_dist_Bmper_Cdh13 <- data.table(cbind(Bmper_zScore_dist=quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper$Bmper, seq(0, 1, 0.1)),
      Bmper_avg_zScore_layerDepth_dist=quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Bmper_layerDepths_avgZscore$avg_Bmper, seq(0, 1, 0.1)),
      Cdh13_zScore_dist=quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13$Cdh13, seq(0, 1, 0.1)),
      Cdh13_avg_zScore_layerDepth_dist=quantile(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt_Cdh13_layerDepths_avgZscore$avg_Cdh13, seq(0, 1, 0.1))), keep.rownames="percentile")

write.csv(zScore_dist_Bmper_Cdh13, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/heatmaps/V1/zScore_distribution_Bmper_Cdh13_sctCount_and_avgZscore_by_layerDepth.csv", row.names=F)



###

exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt$A_avg <- rowMeans(exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt$B_avg <- rowMeans(exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_B_merfish_newNames])
exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt$C_avg <- rowMeans(exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_C_merfish_newNames])

exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt$A_avg <- rowMeans(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt$B_avg <- rowMeans(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_B_merfish_newNames])
exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt$C_avg <- rowMeans(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_C_merfish_newNames])

exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt$A_avg <- rowMeans(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt$B_avg <- rowMeans(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_B_merfish_newNames])
exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt$C_avg <- rowMeans(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_C_merfish_newNames])


exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt$A_avg <- rowMeans(exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_A_merfish_newNames])
exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt$B_avg <- rowMeans(exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_B_merfish_newNames])
exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt$C_avg <- rowMeans(exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt[, ..L23_C_merfish_newNames])

write.csv(exp1_sct_counts2_noNA_V1_L23_glut_scaled_dt, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/Cheng_gene_zScores/exp1_V1_L23_glut_noNA_minVol100_maxVol3Med_minCount300_ChengABC_gene_zScores.csv", row.names=F)
write.csv(exp2_sct_counts2_noNA_V1_L23_glut_scaled_dt, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/Cheng_gene_zScores/exp2_V1_L23_glut_noNA_minVol100_maxVol3Med_minCount300_ChengABC_gene_zScores.csv", row.names=F)
write.csv(exp6_sct_counts2_noNA_V1_L23_glut_scaled_dt, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/Cheng_gene_zScores/exp6_V1_L23_glut_noNA_minVol100_maxVol3Med_minCount300_ChengABC_gene_zScores.csv", row.names=F)
write.csv(exp7_sct_counts2_noNA_V1_L23_glut_scaled_dt, file="HG_lab/Vizgen/analysis_files_from_all_exps/summary_tables/Cheng_gene_zScores/exp7_V1_L23_glut_noNA_minVol100_maxVol3Med_minCount300_ChengABC_gene_zScores.csv", row.names=F)
