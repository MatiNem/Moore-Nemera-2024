library(data.table)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(reshape2)
library(scuttle)
library(scran)
library(tibble)
options(scipen=999)
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
#annotation
CTX_HIP_annot = fread("HG_lab/Mati/GabelLab/RNAseq/Yao2021_mouse/CTX_HIP_Annotation_20190820_annotation_20200913.csv")
CTX_HIP_annot[, cl := as.character(cl)]

#Seurat object
mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2 = readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

mecp2Het.some.meta = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2[[c("rep", "t.type", "predicted.id", "subclass_label", "neighborhood_label")]], keep.rownames = "index")
mecp2Het.some.meta <- left_join(x=mecp2Het.some.meta, y=CTX_HIP_annot[, .(cluster_label, supertype_label)], by=c("predicted.id"="cluster_label"))



mecp2Het.some.meta[!(t.type %in% c("WT", "KO")), t.type:="NA"]

mecp2Het.some.meta.noNA <- mecp2Het.some.meta[t.type %in% c("WT", "KO")]

#summarizing the data table for calculating proportions of each cell type
mecp2Het.some.meta.subclass.summary = group_by(mecp2Het.some.meta , subclass_label, t.type, rep) %>%
   summarise(n.val = n()) %>%
   mutate(prop = n.val / sum(n.val)) %>% data.table
mecp2Het.some.meta.subclass.summary[, prop:=as.numeric(prop)]

#no NA
mecp2Het.some.meta.noNA.subclass.summary = group_by(mecp2Het.some.meta.noNA , subclass_label, t.type, rep) %>%
  summarise(n.val = n()) %>%
  mutate(prop = n.val / sum(n.val)) %>% data.table
mecp2Het.some.meta.noNA.subclass.summary[, prop:=as.numeric(prop)]
##

#colors types, subclasses, neighborhoods for ggplots
name_vec_type = unique(CTX_HIP_annot$cluster_label)
col_vec_type = unique(CTX_HIP_annot$cluster_color)
colors_type = setNames(col_vec_type, name_vec_type)

name_vec_subclass = unique(CTX_HIP_annot$subclass_label)
col_vec_subclass = unique(CTX_HIP_annot$subclass_color)
colors_subclass = setNames(col_vec_subclass, name_vec_subclass)

used_clusters <- unique(mecp2Het.some.meta$predicted.id)
CTX_HIP_annot_used <- CTX_HIP_annot[cluster_label %in% used_clusters]

#ordering purposes
subclass_order <- unique(CTX_HIP_annot_used$subclass_label)
subclass_order_dt <- data.table(subclass_order)
names(subclass_order_dt) = "subclass_label"
subclass_order_dt <- left_join(x=subclass_order_dt, y=CTX_HIP_annot_used[, .(subclass_label, subclass_color, cluster_label, cluster_color)], by=("subclass_label"))

cluster_order <- unique(subclass_order_dt$cluster_label)


mecp2Het.some.meta.subclass.summary <- mecp2Het.some.meta.subclass.summary  %>% mutate(subclass_label = factor(subclass_label, levels=subclass_order))
mecp2Het.some.meta.subclass.summary <- mecp2Het.some.meta.subclass.summary  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO")))

mecp2Het.some.meta.noNA.subclass.summary <- mecp2Het.some.meta.noNA.subclass.summary  %>% mutate(subclass_label = factor(subclass_label, levels=subclass_order))
mecp2Het.some.meta.noNA.subclass.summary <- mecp2Het.some.meta.noNA.subclass.summary  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO")))



#cell type numbers
ggplot(mecp2Het.some.meta.subclass.summary[(rep=="Rep1") & (t.type=="WT") ], aes(x=subclass_label, y=prop, fill=subclass_label)) + 
  ggtitle("Exp1, CTX and HC")+
  geom_bar(position='stack', stat="identity")+
  #coord_cartesian(ylim=c(0,20000))+
  scale_fill_manual(values = colors_subclass[levels(factor(mecp2Het.some.meta.subclass.summary[, subclass_label]))]) +
  ylab("Proportion")+xlab("Subclass")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5), legend.position = "None", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y=element_text(size=10), axis.text.x=element_text(size=6, angle=90))
#ggsave("HG_lab/Mati/GabelLab/Vizgen/barplots/D3A.mutant.exp8.minVol100.minCount300.clusterCountsTotal.stackedBarplot.png", width = 10, height = 6, dpi = 300, units = "in", device='png')
#ggsave("HG_lab/Mati/GabelLab/Vizgen/barplots/D3A.mutant.exp8.minVol100.minCount300.clusterCountsTotal.stackedBarplot.eps", width = 10, height = 6, dpi = 300, units = "in", device='eps')


ggplot(mecp2Het.some.meta.noNA.subclass.summary[rep=="Rep1"], 
       aes(x=subclass_label, y=prop, fill=subclass_label)) + 
  ggtitle("Exp1, CTX and HC") +
  geom_bar(position=position_dodge(width=0.9), stat="identity") +
  scale_fill_manual(values = colors_subclass[levels(factor(mecp2Het.some.meta.noNA.subclass.summary[rep=="Rep1",subclass_label]))]) +
  ylab("Proportion") +
  xlab("Subclass") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5), 
        legend.position = "none", 
        legend.margin=margin(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), 
        axis.ticks.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.text.x=element_text(size=6, angle=90))



ggplot(mecp2Het.some.meta.noNA.subclass.summary[rep=="Rep1"], 
       aes(x=subclass_label, y=prop, fill=t.type)) + 
  ggtitle("Exp1, CTX and HC") +
  geom_bar(aes(fill=t.type), position=position_dodge(width=0.8), stat="identity") +
  #scale_fill_manual(values = colors_subclass[levels(factor(mecp2Het.some.meta.noNA.subclass.summary[rep=="Rep1",subclass_label]))]) +
  scale_fill_manual(name="Transcriptotype:", values = c("WT"="black", "KO"="red")) + # Use subclass colors if needed
  ylab("Proportion") +
  xlab("Subclass") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5), 
        legend.position = "bottom", 
        legend.margin=margin(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line.y=element_line(color="black"), 
        axis.ticks.x = element_blank(), 
        axis.text.y=element_text(size=10), 
        axis.text.x=element_text(size=6, angle=90))

mecp2Het.some.meta.noNA.subclass.summary[rep=="Rep1"]


#Rbp4-Cre and Nr5a-Cre targets
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

Pvalb_types <- CTX_HIP_annot[subclass_label=="Pvalb", cluster_label]
Sst_types <- CTX_HIP_annot[subclass_label=="Sst", cluster_label]

#
mecp2Het.some.meta.INTACTLabs <- rbind(
  cbind(mecp2Het.some.meta[predicted.id %in% Nr5a1_L4_targets], INTACT_label="L4"),
  cbind(mecp2Het.some.meta[predicted.id %in% Rbp4_L5_targets], INTACT_label="L5"),
  cbind(mecp2Het.some.meta[predicted.id %in% Pvalb_types], INTACT_label="PV"),
  cbind(mecp2Het.some.meta[predicted.id %in% Sst_types], INTACT_label="SST")
)

mecp2Het.some.meta.INTACTLabs$animal_label <- "NA"
mecp2Het.some.meta.INTACTLabs[rep=="Rep1", animal_label := "B1"]
mecp2Het.some.meta.INTACTLabs[rep=="Rep2", animal_label := "B2"]
mecp2Het.some.meta.INTACTLabs[rep=="Rep6", animal_label := "B3"]
mecp2Het.some.meta.INTACTLabs[rep=="Rep7", animal_label := "B3"]

#summarizing the data table for calculating proportions of each INTACT subclass
mecp2Het.some.meta.INTACTLabs.summary = group_by(mecp2Het.some.meta.INTACTLabs, INTACT_label, t.type, animal_label) %>%
  summarise(n.val = n()) %>%
  mutate(prop = n.val / sum(n.val)) %>% data.table
mecp2Het.some.meta.INTACTLabs.summary[, prop:=as.numeric(prop)]

###
mecp2Het.some.meta.INTACTLabs.summary2 = group_by(mecp2Het.some.meta.INTACTLabs, INTACT_label, t.type, animal_label) %>%
  summarise(n.val = n(),) %>%
  mutate(sum.val = sum(n.val)) %>%
  mutate(prop = n.val / sum(n.val)) %>% data.table
mecp2Het.some.meta.INTACTLabs.summary2[, prop:=as.numeric(prop)]
##


diff_L4_WT_KO <- mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L4"), prop] - mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L4"), prop]
diff_L5_WT_KO <- mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L5"), prop] - mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L5"), prop]
diff_PV_WT_KO <- mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="PV"), prop] - mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="PV"), prop]
diff_SST_WT_KO <- mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="SST"), prop] - mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="SST"), prop]


shapiro.test(diff_L4_WT_KO)$p.value
shapiro.test(diff_L5_WT_KO)$p.value
shapiro.test(diff_PV_WT_KO)$p.value
shapiro.test(diff_SST_WT_KO)$p.value
#

t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1

t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop], paired=TRUE)$p.value #1


t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop], paired=TRUE)$p.value #1


t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop], paired=TRUE)$p.value #1


t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop], paired=TRUE)$p.value #1


#numbers
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, n.val], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, n.val], paired=TRUE)$p.value #0.2077018
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, n.val], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, n.val], paired=TRUE)$p.value #0.09698513
t.test(mecp2Het.some.meta.INTACTLabs.summary[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, n.val], mecp2Het.some.meta.INTACTLabs.summary[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, n.val], paired=TRUE)$p.value #0.1072143

mecp2Het.some.meta.INTACTLabs.summary <- mecp2Het.some.meta.INTACTLabs.summary  %>% mutate(animal_label = factor(animal_label, levels=c("B1", "B2", "B3")))
mecp2Het.some.meta.INTACTLabs.summary <- mecp2Het.some.meta.INTACTLabs.summary  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO", "NA")))
mecp2Het.some.meta.INTACTLabs.summary <- mecp2Het.some.meta.INTACTLabs.summary  %>% mutate(INTACT_label = factor(INTACT_label, levels=c("L4", "L5", "PV", "SST")))

ggplot(mecp2Het.some.meta.INTACTLabs.summary, aes(x = t.type, y = as.numeric(prop), color=t.type))+
  ggtitle("")+
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1))+
  ylab("Proportion of cells") + xlab("")+
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="gray")) +
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
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/dark_rear/exp21.region.0.1.2.mecp2Het.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
#ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/boxplots/dark_rear/exp21.region.0.1.2.mecp2Het.minVol100.maxVol3Med.min300counts.pred0.2.INTACT.MR.genes.INTACTLevel.pseudobulkDGE.logFC.sctCounts.Under1Over1.boxplot.eps", width = 4, height = 6, dpi = 300, units = "in", device='eps')



###
mecp2Het.some.meta2 <- copy(mecp2Het.some.meta)
mecp2Het.some.meta2$animal_label <- "NA"
mecp2Het.some.meta2[rep=="Rep1", animal_label := "B1"]
mecp2Het.some.meta2[rep=="Rep2", animal_label := "B2"]
mecp2Het.some.meta2[rep=="Rep6", animal_label := "B3"]
mecp2Het.some.meta2[rep=="Rep7", animal_label := "B3"]


#calculate number of cells with combination of predicted.id, animal_label, and t.type, as well as proportion of cells of that t.type and animal_label with the same predicted.id
mecp2Het.some.meta2.type.summary <- mecp2Het.some.meta2 %>%
  group_by(predicted.id, t.type, animal_label) %>%
  summarise(n.val = n(), .groups = 'drop') %>%
  group_by(t.type, animal_label) %>%
  mutate(prop.val = n.val / sum(n.val)) %>%
  ungroup() %>%
  data.table()


mecp2Het.some.meta2.type.summary.INTACTLabs <- rbind(
  cbind(mecp2Het.some.meta2.type.summary[predicted.id %in% Nr5a1_L4_targets], INTACT_label="L4"),
  cbind(mecp2Het.some.meta2.type.summary[predicted.id %in% Rbp4_L5_targets], INTACT_label="L5"),
  cbind(mecp2Het.some.meta2.type.summary[predicted.id %in% Pvalb_types], INTACT_label="PV"),
  cbind(mecp2Het.some.meta2.type.summary[predicted.id %in% Sst_types], INTACT_label="SST")
)

mecp2Het.some.meta2.type.summary.INTACTLabs <- group_by(mecp2Het.some.meta2.type.summary.INTACTLabs, INTACT_label, t.type, animal_label) %>%
  summarise(INTACT_total = sum(n.val)) %>% data.table

t.type.animal.count.summary <- mecp2Het.some.meta2 %>%
  group_by(t.type, animal_label) %>%
  summarise(total_cells = n()) %>% data.table

mecp2Het.some.meta.INTACTLabs.join <- left_join(x=mecp2Het.some.meta.INTACTLabs.summary, t.type.animal.count.summary, by=c("t.type", "animal_label"))

mecp2Het.some.meta.INTACTLabs.join$prop_of_total = mecp2Het.some.meta.INTACTLabs.join$n.val / mecp2Het.some.meta.INTACTLabs.join$total_cells

mecp2Het.some.meta.INTACTLabs.join <- mecp2Het.some.meta.INTACTLabs.join  %>% mutate(animal_label = factor(animal_label, levels=c("B1", "B2", "B3")))
mecp2Het.some.meta.INTACTLabs.join <- mecp2Het.some.meta.INTACTLabs.join  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO", "NA")))
mecp2Het.some.meta.INTACTLabs.join <- mecp2Het.some.meta.INTACTLabs.join  %>% mutate(INTACT_label = factor(INTACT_label, levels=c("L4", "L5", "PV", "SST")))


ggplot(mecp2Het.some.meta.INTACTLabs.join, aes(x = t.type, y = as.numeric(prop_of_total), color=t.type))+
  ggtitle("")+
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.15))+
  ylab("Proportion of cells") + xlab("")+
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
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
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT.total.transcriptotype.proportions.by.exp.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/INTACT.total.transcriptotype.proportions.by.exp.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)




diff_L4_WT_KO <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4"), prop_of_total]
diff_L5_WT_KO <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5"), prop_of_total]
diff_PV_WT_KO <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV"), prop_of_total]
diff_SST_WT_KO <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST"), prop_of_total]

diff_L4_WT_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4"), prop_of_total]
diff_L5_WT_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5"), prop_of_total]
diff_PV_WT_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV"), prop_of_total]
diff_SST_WT_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST"), prop_of_total]

diff_L4_KO_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4"), prop_of_total]
diff_L5_KO_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5"), prop_of_total]
diff_PV_KO_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV"), prop_of_total]
diff_SST_KO_NA <- mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST"), prop_of_total] - mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST"), prop_of_total]



shapiro.test(diff_L4_WT_KO)$p.value #p=0.7274099
shapiro.test(diff_L5_WT_KO)$p.value #p=0.02651254
shapiro.test(diff_PV_WT_KO)$p.value #p=0.3514851
shapiro.test(diff_SST_WT_KO)$p.value #p=0.1897221

shapiro.test(diff_L4_WT_NA)$p.value #p=0.4190694
shapiro.test(diff_L5_WT_NA)$p.value #p=0.3677651
shapiro.test(diff_PV_WT_NA)$p.value #p=0.1901459
shapiro.test(diff_SST_WT_NA)$p.value #p=0.0006013113

shapiro.test(diff_L4_KO_NA)$p.value #p=0.6081674
shapiro.test(diff_L5_KO_NA)$p.value #p=0.5033553
shapiro.test(diff_PV_KO_NA)$p.value #p=0.8940575
shapiro.test(diff_SST_KO_NA)$p.value #p=0.7909877

#paired t-tests
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.6326021
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.1311154
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.08820824

t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.8120726
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.4117509
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.3133418


t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.2810344
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.469442
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.2283557


t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.01832333
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.20975
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.2757628


t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.02259281
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.02357964
t.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.0447392

#paired wilcoxon
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.75
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25

wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #1
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.5
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25


wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.5
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.75
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25


wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.25
wilcox.test(mecp2Het.some.meta.INTACTLabs.join[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], mecp2Het.some.meta.INTACTLabs.join[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_total], paired=TRUE)$p.value #0.5



##proportion of each INTACT subclass that is each t.type
mecp2Het.some.meta2.type.summary.INTACTLabs.sum <- group_by(mecp2Het.some.meta2.type.summary.INTACTLabs, INTACT_label, animal_label) %>%
  summarise(
    animal_INTACT_total = sum(INTACT_total),
    .groups = 'drop' # This ensures the result is not grouped by default
  ) %>% data.table

#join INTACT subclass-animal sum table with INTACT subclass-animal-t.type sum table
mecp2Het.some.meta2.type.summary.INTACTLabs.prop <- left_join(x=mecp2Het.some.meta2.type.summary.INTACTLabs, y=mecp2Het.some.meta2.type.summary.INTACTLabs.sum, by=c("INTACT_label", "animal_label"))
#calculate ratio of the number cells with a given INTACT_label, animal_label, and t.type to the number of all cells of the same INTACT_label and animal_label
mecp2Het.some.meta2.type.summary.INTACTLabs.prop$prop_of_subclass_with_ttype <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop[, (INTACT_total)/(animal_INTACT_total)]
#c

mecp2Het.some.meta2.type.summary.INTACTLabs.prop <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop  %>% mutate(animal_label = factor(animal_label, levels=c("B1", "B2", "B3")))
mecp2Het.some.meta2.type.summary.INTACTLabs.prop <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO", "NA")))
mecp2Het.some.meta2.type.summary.INTACTLabs.prop <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop  %>% mutate(INTACT_label = factor(INTACT_label, levels=c("L4", "L5", "PV", "SST")))


ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type))+
  ggtitle("")+
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  #coord_cartesian(ylim=c(0,0.15))+
  ylab("Proportion of cells") + xlab("")+
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
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
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type, shape=animal_label)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  #coord_cartesian(ylim=c(0,0.15)) +
  ylab("Proportion of cells") + xlab("") +
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
  scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~INTACT_label, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE), shape = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


###all cells summary
# Group by t.type and animal_label, then summarize
mecp2Het.some.meta2.all.t.type.summary <- mecp2Het.some.meta2 %>%
  group_by(t.type, animal_label) %>%
  summarise(
    INTACT_total = n(),  # Count the number of rows for each t.type and animal_label
  ) %>%
  ungroup()

# Calculate total_cells for each animal_label
mecp2Het.some.meta2.all.t.type.summary <- mecp2Het.some.meta2.all.t.type.summary %>%
  left_join(
    mecp2Het.some.meta2 %>%
      group_by(animal_label) %>%
      summarise(animal_INTACT_total = n()),  # Calculate total cells for each animal_label
    by = "animal_label"
  ) %>% data.table

mecp2Het.some.meta2.all.t.type.summary$prop_of_subclass_with_ttype <- mecp2Het.some.meta2.all.t.type.summary[, (INTACT_total)/(animal_INTACT_total)]

mecp2Het.some.meta2.all.t.type.summary$INTACT_label <- "All"
mecp2Het.some.meta2.all.t.type.summary <- mecp2Het.some.meta2.all.t.type.summary[, .(INTACT_label, t.type, animal_label, INTACT_total, animal_INTACT_total, prop_of_subclass_with_ttype)]

#join all with INTACT subclass proportion table
mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll <-  rbind(mecp2Het.some.meta2.all.t.type.summary, mecp2Het.some.meta2.type.summary.INTACTLabs.prop)

mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll  %>% mutate(animal_label = factor(animal_label, levels=c("B1", "B2", "B3")))
mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll  %>% mutate(t.type = factor(t.type, levels=c("WT", "KO", "NA")))
mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll  %>% mutate(INTACT_label = factor(INTACT_label, levels=c("All", "L4", "L5", "PV", "SST")))



ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type))+
  ggtitle("")+
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1))+
  ylab("Proportion of cells") + xlab("")+
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
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
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type, shape=animal_label)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("Proportion of cells") + xlab("") +
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
  scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~INTACT_label, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE), shape = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/proportion.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)


ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type, shape=animal_label)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("Proportion of cells") + xlab("") +
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
  scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~INTACT_label, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE), shape = guide_legend(nrow=3, byrow=TRUE))


mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[t.type %in% c("WT", "KO")]

new_diff_All_WT_KO <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="All"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="All"), prop_of_subclass_with_ttype]
new_diff_All_KO_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="All"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="All"), prop_of_subclass_with_ttype]
new_diff_All_WT_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="All"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="All"), prop_of_subclass_with_ttype]


new_diff_L4_WT_KO <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L4"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L4"), prop_of_subclass_with_ttype]
new_diff_L4_KO_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L4"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L4"), prop_of_subclass_with_ttype]
new_diff_L4_WT_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L4"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L4"), prop_of_subclass_with_ttype]


new_diff_L5_WT_KO <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L5"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L5"), prop_of_subclass_with_ttype]
new_diff_L5_KO_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L5"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L5"), prop_of_subclass_with_ttype]
new_diff_L5_WT_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L5"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L5"), prop_of_subclass_with_ttype]


new_diff_PV_WT_KO <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="PV"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="PV"), prop_of_subclass_with_ttype]
new_diff_PV_KO_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="PV"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="PV"), prop_of_subclass_with_ttype]
new_diff_PV_WT_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="PV"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="PV"), prop_of_subclass_with_ttype]

new_diff_SST_WT_KO <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="SST"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="SST"), prop_of_subclass_with_ttype]
new_diff_SST_KO_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="SST"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="SST"), prop_of_subclass_with_ttype]
new_diff_SST_WT_NA <- mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="SST"), prop_of_subclass_with_ttype] - mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="SST"), prop_of_subclass_with_ttype]


shapiro.test(new_diff_All_WT_KO)$p.value #p=0.7096474
shapiro.test(new_diff_All_KO_NA)$p.value #p=0.898721
shapiro.test(new_diff_All_WT_NA)$p.value #p=0.5535426

shapiro.test(new_diff_L4_WT_KO)$p.value #p=0.5239123
shapiro.test(new_diff_L4_KO_NA)$p.value #p=0.9248865
shapiro.test(new_diff_L4_WT_NA)$p.value #p=0.3602667

shapiro.test(new_diff_L5_WT_KO)$p.value #p=0.7321901
shapiro.test(new_diff_L5_KO_NA)$p.value #p=0.95621
shapiro.test(new_diff_L5_WT_NA)$p.value #p=0.6171516

shapiro.test(new_diff_PV_WT_KO)$p.value #p=0.3038825
shapiro.test(new_diff_PV_KO_NA)$p.value #p=0.1144727
shapiro.test(new_diff_PV_WT_NA)$p.value #p=0.485119

shapiro.test(new_diff_SST_WT_KO)$p.value #p=0.8804252
shapiro.test(new_diff_SST_KO_NA)$p.value #p=0.693635
shapiro.test(new_diff_SST_WT_NA)$p.value #p=0.9835502

t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.6266634, ns
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.01147746, *
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="All")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.01301677, *

t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.4454429, ns
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.004880728, **
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L4")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.01888321, *

t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.7549387, ns
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.00263715, ** 
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="L5")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.008888044, **

t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.4240908, ns
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.04915351, *
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="PV")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.0258688, *

t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.2629006, ns
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="KO") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.0113638, *
t.test(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="WT") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll[(t.type=="NA") & (INTACT_label=="SST")][order(animal_label)][, prop_of_subclass_with_ttype], paired=TRUE)$p.value #0.01104687, *




#bigger
ggplot(mecp2Het.some.meta2.type.summary.INTACTLabs.prop.withAll, aes(x = t.type, y = as.numeric(prop_of_subclass_with_ttype), color=t.type, shape=animal_label)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 4, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("Proportion of cells") + xlab("") +
  scale_color_manual(name = "Transcriptotype", values = c("WT"="purple", "KO"="orange", "NA"="darkgray")) +
  scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~INTACT_label, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE), shape = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/big.prop.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/big.prop.of.allCells.and.INTACT.subclass.per.animal.with.transcriptotype.mecp2Het.CTX.HC.minVol100.maxVol3Med.min300counts.pred0.2.Under1Over1.dotplot.animalShapes.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)
