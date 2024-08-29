library(data.table)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(reshape2)
library(scuttle)
library(scran)
library(tibble)
options(scipen=999)

#Yao 2021 metadata
yao_meta_dt = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/metadata.csv", na.strings=c("", "NA"))
yao_meta_dt = yao_meta_dt[!is.na(cluster_label), ]
#subclass assignments of cluster labels from Yao 2021
cluster_subclass = unique(yao_meta_dt[, .(cluster_label, subclass_label)])
yao2021_DE = readRDS("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/comb.de.genes.Zizhen.Yao.rda")


CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
CTX_HIP_annot[, cl := as.character(cl)]

#read in Seurat object containing your data
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")
#you can see what available metadata there is using mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes@meta.data or names(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[]])

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.meta <- mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes@meta.data %>% data.table

#average expressions in 100vol, 300 count filter
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT = subset(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, subset = t.type == "WT")
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.avg = AverageExpression(obj=mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, assays="SCT", slot="data", group.by = c("predicted.id"))$SCT


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

Pvalb_targets <- unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.meta[subclass_label=="Pvalb", predicted.id])
Sst_targets <- unique(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.meta[subclass_label=="Sst", predicted.id])


#type DEGs
type1.type2.de.noshrink.func <- function(seurat_obj, type1, type2){
  sce <- subset(seurat_obj, subset = ((predicted.id == type1) | (predicted.id == type2)))
  sce <- subset(sce, subset = ((t.type == "WT")))
  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(sce[["RNA"]], slot = "counts")), 
                              colData = sce@meta.data)
  #here, predicted.id is the metadata name for subtype;
  groups.sce <- colData(sce)[, c("predicted.id", "t.type", "rep")]
  aggr_counts <- aggregateAcrossCells(sce, ids=groups.sce)
  aggr_counts.filt <- aggr_counts[,aggr_counts$ncells >= 10]
  col.data = colData(aggr_counts.filt)
  col.data$predicted.id = factor(col.data$predicted.id, levels=c(type2, type1))
  de.results.noShrink <- pseudoBulkDGE(x=aggr_counts.filt,
                                       #all the labels of data
                                       label=aggr_counts.filt$t.type,
                                       col.data=col.data,
                                       #design of the model matrix
                                       design=~ rep + predicted.id,
                                       #String or character vector containing the coefficients to drop from the design matrix to form the null hypothesis; when comparing t.type=KO and t.type=WT where you want the output to be KO/WT, it will be t.typeKO
                                       coef=paste0("predicted.id", type1),
                                       #A vector or factor of length equal to ncol(x), specifying the experimental condition for each column of x. Only used for abundance-based filtering of genes.
                                       condition=aggr_counts.filt$predicted.id,
                                       robust=FALSE)
  return(de.results.noShrink)
}

#processing DEG result file to a data table
type_de_results_func <- function(de.results, id, type1_col, type2_col){
  de.results.list = as.list(de.results)
  de.results.list = lapply(de.results.list, as.data.frame)
  de.results.list = lapply(de.results.list, rownames_to_column)
  de.results.dt = data.table(rbindlist(de.results.list, idcol = id))
  names(de.results.dt)[2] = "gene"
  de.results.dt$p.adj.BH.all = p.adjust(de.results.dt$PValue, method="BH")
  #adding columns saying what is type 1 and what is type 2
  de.results.dt$type1  = type1_col
  de.results.dt$type2  = type2_col
  return(de.results.dt)
}

Pvalb.116.vs.114.de.results <- type1.type2.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, type1="116_Pvalb", type2="114_Pvalb")

Pvalb.116.vs.114.de.results.dt <- type_de_results_func(de.results=Pvalb.116.vs.114.de.results, id="t.type", type1_col="116_Pvalb", type2_col="114_Pvalb")


execution_time <- system.time({
  Pvalb.116.vs.114.de.results <- type1.type2.de.noshrink.func(
    mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, 
    type1 = "116_Pvalb", 
    type2 = "114_Pvalb"
  )
})

print(execution_time)

#checking to see if directionality is correct; Serpinf1 should be higher in 116_Pvalb than 114_Pvalb in expression; and it is.
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes.WT.avg["Serpinf1", c("116_Pvalb", "114_Pvalb")]

#binary comparisons from Yao 2021 of types that are in our merfish data
yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2 <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Yao2021_pairwise_DEGs/yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.comparisons.csv")

yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(subclass_label.x=="Pvalb") & (subclass_label.y=="Pvalb")] 


yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb <- yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(cluster_label.x %in% Pvalb_targets) & (cluster_label.y %in% Pvalb_targets)]
yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst <- yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(cluster_label.x %in% Sst_targets) & (cluster_label.y %in% Sst_targets)]
yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4 <- yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(cluster_label.x %in% Nr5a1_L4_targets) & (cluster_label.y %in% Nr5a1_L4_targets)]
yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5 <- yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2[(cluster_label.x %in% Rbp4_L5_targets) & (cluster_label.y %in% Rbp4_L5_targets)]

# Define the column names
DEG_column_names <- c("t.type", "gene", "logFC", "logCPM", "F", "PValue", "FDR", "p.adj.BH.all", "type1", "type2")

# Create an empty data.table with the specified column names
Pvalb.de.results.dt <- data.table(matrix(ncol = length(DEG_column_names), nrow = 0))
setnames(Pvalb.de.results.dt, DEG_column_names)

# Create an empty data.table with the specified column names
Sst.de.results.dt <- data.table(matrix(ncol = length(DEG_column_names), nrow = 0))
setnames(Sst.de.results.dt, DEG_column_names)


# Create an empty data.table with the specified column names
L4.de.results.dt <- data.table(matrix(ncol = length(DEG_column_names), nrow = 0))
setnames(L4.de.results.dt, DEG_column_names)


# Create an empty data.table with the specified column names
L5.de.results.dt <- data.table(matrix(ncol = length(DEG_column_names), nrow = 0))
setnames(L5.de.results.dt, DEG_column_names)

for(i in 1:nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb)){
  #define types
  Pvalb_type_1 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb[i, cluster_label.x]
  Pvalb_type_2 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb[i, cluster_label.y]
  #pseudobulkDGE
  Pvalb.de.results <- type1.type2.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, type1=Pvalb_type_1, type2=Pvalb_type_2)
  if(length(Pvalb.de.results) > 0){
    #make a readable table
    Pvalb.de.results.dataTable <- type_de_results_func(de.results=Pvalb.de.results, id="t.type", type1_col=Pvalb_type_1, type2_col=Pvalb_type_2)
  } else{
    Pvalb.de.results.dataTable <- cbind(t.type="WT", gene=NA, logFC=NA, logCPM=NA, F=NA, PValue=NA, FDR=NA, p.adj.BH.all=NA, type1=Pvalb_type_1, type2=Pvalb_type_2) %>% data.table
  }
  #add to existing table of same subclass
  Pvalb.de.results.dt <- rbind(Pvalb.de.results.dt, Pvalb.de.results.dataTable)
  print(paste0(i, "/", nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb)))
}

write.csv(Pvalb.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/type_vs_type_DEG_tables/Pvalb_type_vs_type_WT_pseudobulkDGE_exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_table.csv")


#Sst
for(i in 1:nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst)){
  #define types
  Sst_type_1 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst[i, cluster_label.x]
  Sst_type_2 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst[i, cluster_label.y]
  #pseudobulkDGE
  Sst.de.results <- type1.type2.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, type1=Sst_type_1, type2=Sst_type_2)
  if(length(Sst.de.results) > 0){
    #make a readable table
    Sst.de.results.dataTable <- type_de_results_func(de.results=Sst.de.results, id="t.type", type1_col=Sst_type_1, type2_col=Sst_type_2)
  } else{
    Sst.de.results.dataTable <- cbind(t.type="WT", gene=NA, logFC=NA, logCPM=NA, F=NA, PValue=NA, FDR=NA, p.adj.BH.all=NA, type1=Sst_type_1, type2=Sst_type_2) %>% data.table
  }
  #add to existing table of same subclass
  Sst.de.results.dt <- rbind(Sst.de.results.dt, Sst.de.results.dataTable)
  print(paste0(i, "/", nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst)))
}

write.csv(Sst.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/type_vs_type_DEG_tables/Sst_type_vs_type_WT_pseudobulkDGE_exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_table.csv")

#L4
for(i in 1:nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4)){
  #define types
  L4_type_1 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4[i, cluster_label.x]
  L4_type_2 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4[i, cluster_label.y]
  #pseudobulkDGE
  L4.de.results <- type1.type2.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, type1=L4_type_1, type2=L4_type_2)
  if(length(L4.de.results) > 0){
    #make a readable table
    L4.de.results.dataTable <- type_de_results_func(de.results=L4.de.results, id="t.type", type1_col=L4_type_1, type2_col=L4_type_2)
  } else{
    L4.de.results.dataTable <- cbind(t.type="WT", gene=NA, logFC=NA, logCPM=NA, F=NA, PValue=NA, FDR=NA, p.adj.BH.all=NA, type1=L4_type_1, type2=L4_type_2) %>% data.table
  }
  #add to existing table of same subclass
  L4.de.results.dt <- rbind(L4.de.results.dt, L4.de.results.dataTable)
  print(paste0(i, "/", nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4)))
}

write.csv(L4.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/type_vs_type_DEG_tables/Nr5a1_L4_type_vs_type_WT_pseudobulkDGE_exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_table.csv")

#L5
for(i in 1:nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5)){
  #define types
  L5_type_1 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5[i, cluster_label.x]
  L5_type_2 = yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5[i, cluster_label.y]
  #pseudobulkDGE
  L5.de.results <- type1.type2.de.noshrink.func(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes, type1=L5_type_1, type2=L5_type_2)
  if(length(L5.de.results) > 0){
    #make a readable table
    L5.de.results.dataTable <- type_de_results_func(de.results=L5.de.results, id="t.type", type1_col=L5_type_1, type2_col=L5_type_2)
  } else{
    L5.de.results.dataTable <- cbind(t.type="WT", gene=NA, logFC=NA, logCPM=NA, F=NA, PValue=NA, FDR=NA, p.adj.BH.all=NA, type1=L5_type_1, type2=L5_type_2) %>% data.table
  }
  #add to existing table of same subclass
  L5.de.results.dt <- rbind(L5.de.results.dt, L5.de.results.dataTable)
  print(paste0(i, "/", nrow(yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5)))
}

write.csv(L5.de.results.dt, file="HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/type_vs_type_DEG_tables/Rbp4_L5_type_vs_type_WT_pseudobulkDGE_exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_CTX_HC_table.csv")

#
yao2021_DE[[]]

head(Pvalb.de.results.dt)


Pvalb.de.results.dt2 <- left_join(x=Pvalb.de.results.dt, yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Pvalb, by=c("type1"="cluster_label.x", "type2"="cluster_label.y"))
Sst.de.results.dt2 <- left_join(x=Sst.de.results.dt, yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Sst, by=c("type1"="cluster_label.x", "type2"="cluster_label.y"))
L4.de.results.dt2 <- left_join(x=L4.de.results.dt, yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L4, by=c("type1"="cluster_label.x", "type2"="cluster_label.y"))
L5.de.results.dt2 <- left_join(x=L5.de.results.dt, yao.binary.mecp2Het.CTX.HC.100vol.300counts.pred0.2.L5, by=c("type1"="cluster_label.x", "type2"="cluster_label.y"))


names(yao2021_DE[[]][1])


names(yao2021_DE[[c(Pvalb.de.results.dt2[1, comparison], "up.genes")]])

#for(i in unique(binary_comparison_table[(subclass_label.x%in%subclass_label), cl.x])){
#  for(j in binary_comparison_table[(subclass_label.x%in%subclass_label) & (cl.x==i), comparison]){
#   up.genes = names(yao2021_DE[[c(j, "up.genes")]])
#    down.genes= names(yao2021_DE[[c(j, "down.genes")]])
#    non_DEGs = setdiff(gene.by.subtype.avg.counts.scale[, gene], c(up.genes, down.genes))


#type comparisons of interest
all_compar_of_interest <- unique(c(Pvalb.de.results.dt2$comparison, Sst.de.results.dt2$comparison, L4.de.results.dt2$comparison, L5.de.results.dt2$comparison))
# Initialize an empty list to store data.tables
up_down_yao2021_DE <- list()

# Loop through each comparison in yao2021_DE
for (comparison in all_compar_of_interest) {
  # Extract up and down genes
  up_genes <- names(yao2021_DE[[c(comparison, "up.genes")]])
  down_genes <- names(yao2021_DE[[c(comparison, "down.genes")]])
  
  # Create data.tables for up and down genes
  dt_up <- data.table(comparison = comparison, dir = "High", gene = up_genes)
  dt_down <- data.table(comparison = comparison, dir = "Low", gene = down_genes)
  
  # Append to the result list
  up_down_yao2021_DE <- c(up_down_yao2021_DE, list(dt_up, dt_down))
}

# Combine all data.tables in the list into a single data.table
up_down_yao2021_DE_dt <- rbindlist(up_down_yao2021_DE)


#adding high and low gene identities
Pvalb.de.results.dt.geneDir <- left_join(x=Pvalb.de.results.dt2, y=up_down_yao2021_DE_dt, by=c("comparison", "gene"))
Sst.de.results.dt.geneDir <- left_join(x=Sst.de.results.dt2, y=up_down_yao2021_DE_dt, by=c("comparison", "gene"))
L4.de.results.dt.geneDir <- left_join(x=L4.de.results.dt2, y=up_down_yao2021_DE_dt, by=c("comparison", "gene"))
L5.de.results.dt.geneDir <- left_join(x=L5.de.results.dt2, y=up_down_yao2021_DE_dt, by=c("comparison", "gene"))


#pseudobulkDGE of KO/WT
typeAgg.KO.WT.de.results.dt <- fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/DE_genes/CTX_HC/typeAgg.mecp2Het.CTX.only.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.csv")

#renaming some columns
typeAgg.KO.WT.de.results.dt2 <- copy(typeAgg.KO.WT.de.results.dt)
typeAgg.KO.WT.de.results.dt2 <- typeAgg.KO.WT.de.results.dt2[type %in% c(Pvalb_targets, Sst_targets, Rbp4_L5_targets)]
typeAgg.KO.WT.de.results.dt2$dir ="KO/WT"
typeAgg.KO.WT.de.results.dt2$INTACT_subclass <- "NA"
typeAgg.KO.WT.de.results.dt2[type %in% Pvalb_targets, INTACT_subclass := "PV"]
typeAgg.KO.WT.de.results.dt2[type %in% Sst_targets, INTACT_subclass := "SST"]
typeAgg.KO.WT.de.results.dt2[type %in% Rbp4_L5_targets, INTACT_subclass := "L5"]

#combining data tables containing fold changes of DEGs between types and also KO/WT fold changes for each type
highLow_KOWT_logFC_dt <- rbind(
  cbind(Pvalb.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="PV"),
  cbind(Sst.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="SST"),
  cbind(L5.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="L5"),
  cbind(typeAgg.KO.WT.de.results.dt2[, .(type, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir, INTACT_subclass)])
)

highLow_KOWT_logFC_dt <- highLow_KOWT_logFC_dt[dir %in% c("High", "Low", "KO/WT")]

highLow_KOWT_logFC_dt = highLow_KOWT_logFC_dt %>% mutate(dir = factor(dir, levels=c("High", "Low", "KO/WT")))
highLow_KOWT_logFC_dt = highLow_KOWT_logFC_dt %>% mutate(INTACT_subclass = factor(INTACT_subclass, levels=c("L5", "PV", "SST")))

ggplot(highLow_KOWT_logFC_dt, aes(x = INTACT_subclass, y = as.numeric(logFC), fill=INTACT_subclass))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("L5"="orchid", "PV"="forestgreen", "SST"="orange"))+
  coord_cartesian(ylim=c(-5.5,5.5))+
  ylab("Log2 fold difference") + xlab("")+
  facet_grid(.~dir,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWTall_logFC_INTACT_subclasses_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWTall_logFC_INTACT_subclasses_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')



##
highLow_dt <- rbind(
  cbind(Pvalb.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="PV"),
  cbind(Sst.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="SST"),
  cbind(L5.de.results.dt.geneDir[, .(type=type1, gene, logFC, logCPM, F, PValue, FDR, p.adj.BH.all, dir)], INTACT_subclass="L5")
)

highLow_dt <- highLow_dt[dir %in% c("High", "Low")]

highLow_dt_matched_KOWT <- left_join(x=highLow_dt, y=typeAgg.KO.WT.de.results.dt, by=c("type", "gene"), suffix=c(".highLow", ".KOWT"))

highLow_dt_matched_KOWT_melt <- melt(highLow_dt_matched_KOWT, id.vars=c("type", "gene", "dir", "INTACT_subclass"), measure.vars=c("logFC.highLow", "logFC.KOWT")) %>% data.table

ggplot(highLow_dt_matched_KOWT_melt, aes(x = INTACT_subclass, y = as.numeric(value), fill=INTACT_subclass))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("L5"="orchid", "PV"="forestgreen", "SST"="orange"))+
  coord_cartesian(ylim=c(-5.5,5.5))+
  ylab("Log2 fold difference") + xlab("")+
  facet_grid(.~variable+dir,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_sameGenes_logFC_INTACT_subclasses_boxplot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_sameGenes_INTACT_subclasses_boxplot.eps", width = 4, height = 5, dpi = 300, units = "in", device='eps')



###
all.nhoods.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.allINTACTMR.table <- fread("HG_lab/Vizgen/analysis_files_from_all_exps/DE_genes/CTX_HC/all.neighborhoods.lowHigh.geneClasses.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.sctransformMecp2Counts.pseudoBulkDGE.noShrink.txt")

Pvalb_logFC_table=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.Pvalb.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv")
Sst_logFC_table=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.Sst.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv")
L5_logFC_table=fread("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/manuscript_plots/violinplots/table.noOutlierRemoval.L5.mecp2Het.CTX.HC.100vol.300counts.pred0.2.Under1Over1.nonDEG.lowHigh.recurrent.MR.logFC.csv")

#doing opposite value of log2  fold change 
highLow_dt2 <- copy(highLow_dt)
highLow_dt2[, logFC := as.numeric(logFC)]
highLow_dt2[dir == "Low", logFC := -logFC]

highLow_dt2



highDEGs_MR_KOWT_logFC <- rbind(
  cbind(Pvalb.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="PV"),
  cbind(Pvalb.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="PV"),
  cbind(Sst.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="SST"),
  cbind(Sst.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="SST"),
  cbind(L5.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="L5"),
  cbind(L5.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="L5"),
  cbind(Pvalb_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="PV"),
  cbind(Pvalb_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="PV"),
  cbind(Sst_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="SST"),
  cbind(Sst_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="SST"),
  cbind(L5_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="L5"),
  cbind(L5_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="L5")
)

highDEGs_MR_KOWT_logFC = highDEGs_MR_KOWT_logFC %>% mutate(INTACT_subclass = factor(INTACT_subclass, levels=c("L5", "PV", "SST")))
highDEGs_MR_KOWT_logFC = highDEGs_MR_KOWT_logFC %>% mutate(group.ID = factor(group.ID, levels=c("High type/low type, all DEGs", "KO/WT, low MR DEGs", "KO/WT, high MR DEGs")))


ggplot(highDEGs_MR_KOWT_logFC, aes(x = INTACT_subclass, y = as.numeric(logFC), fill=INTACT_subclass))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  scale_fill_manual(name = "", values = c("L5"="orchid", "PV"="forestgreen", "SST"="orange"))+
  coord_cartesian(ylim=c(-2,5))+
  ylab("Log2 fold difference") + xlab("")+
  facet_grid(.~group.ID,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text = element_text(size = 8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_anyMetaMR_DEGs_logFC_INTACT_subclasses_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_anyMetaMR_DEGs_logFC_INTACT_subclasses_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')


##include non DEGs
highDEGs_MR_KOWT_all_and_nonDEGs_logFC <- rbind(
  cbind(Pvalb.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="PV"),
  cbind(Pvalb.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="PV"),
  cbind(Sst.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="SST"),
  cbind(Sst.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="SST"),
  cbind(L5.de.results.dt.geneDir[dir=="High", .(type=type1, gene, logFC=as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="L5"),
  cbind(L5.de.results.dt.geneDir[dir=="Low", .(type=type2, gene, logFC=-as.numeric(logFC))], group.ID="High type/low type, all DEGs", INTACT_subclass="L5"),
  cbind(Pvalb_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="PV"),
  cbind(Pvalb_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="PV"),
  cbind(Pvalb_logFC_table[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG"), .(type=subtype, gene, logFC)], group.ID="KO/WT, non-DEGs", INTACT_subclass="PV"),
  cbind(Sst_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="SST"),
  cbind(Sst_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="SST"),
  cbind(Sst_logFC_table[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG"), .(type=subtype, gene, logFC)], group.ID="KO/WT, non-DEGs", INTACT_subclass="SST"),
  cbind(L5_logFC_table[(gene_class=="Any recurrent MR") & (dir=="High"), .(type=subtype, gene, logFC)], group.ID="KO/WT, high MR DEGs", INTACT_subclass="L5"),
  cbind(L5_logFC_table[(gene_class=="Any recurrent MR") & (dir=="Low"), .(type=subtype, gene, logFC)], group.ID="KO/WT, low MR DEGs", INTACT_subclass="L5"),
  cbind(L5_logFC_table[(gene_class=="Non-MR, non-MA") & (dir=="Non-DEG"), .(type=subtype, gene, logFC)], group.ID="KO/WT, non-DEGs", INTACT_subclass="L5")
)

highDEGs_MR_KOWT_all_and_nonDEGs_logFC = highDEGs_MR_KOWT_all_and_nonDEGs_logFC %>% mutate(INTACT_subclass = factor(INTACT_subclass, levels=c("L5", "PV", "SST")))
highDEGs_MR_KOWT_all_and_nonDEGs_logFC = highDEGs_MR_KOWT_all_and_nonDEGs_logFC %>% mutate(group.ID = factor(group.ID, levels=c("High type/low type, all DEGs", "KO/WT, low MR DEGs", "KO/WT, high MR DEGs", "KO/WT, non-DEGs")))

ggplot(highDEGs_MR_KOWT_all_and_nonDEGs_logFC, aes(x = INTACT_subclass, y = as.numeric(logFC), fill=INTACT_subclass))+
  ggtitle("")+
  stat_boxplot(geom='errorbar', width=0.25)+
  geom_boxplot(outlier.shape = NA, notch=TRUE) + 
  #scale_fill_manual(name = "", values = c("L5"="orchid", "PV"="forestgreen", "SST"="orange"))+
  scale_fill_manual(name = "", values = c("PV"=unique(CTX_HIP_annot[subclass_label=="Pvalb", subclass_color]), 
                                          "SST"=unique(CTX_HIP_annot[subclass_label=="Sst", subclass_color]), 
                                          "L5"="#50B2AD"))+ 
  coord_cartesian(ylim=c(-2,5))+
  ylab("Log2 fold difference") + xlab("")+
  facet_grid(.~group.ID,
             switch = "x", # Moves the labels from the top to the bottom
             labeller = label_value # Adds the labels
  ) +
  theme_bw()+
  theme(strip.background =element_rect(fill="white"), strip.text = element_text(size = 8))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.title.x=element_blank(), axis.text.y=element_text(size=14), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_anyMetaMR_and_non_DEGs_logFC_INTACT_subclasses_boxplot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/Yao2021_DEGs_highType_lowType_KOWT_anyMetaMR_and_non_DEGs_logFC_INTACT_subclasses_boxplot.eps", width = 5, height = 5, dpi = 300, units = "in", device='eps')
