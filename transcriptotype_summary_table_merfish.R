library(data.table)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(reshape2)
library(scuttle)
library(scran)
library(tibble)
options(scipen=999)

mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes=readRDS("HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/Seurat_objects/CTX_HC/exp1_exp2_exp6_exp7_minVol100_maxVol3Med_min300counts_pred0.2_Under1Over1_cellTypes_CTX_HC_Seurat_SCT_integrated_yao2021_smartseq_labelTransfer.rda")

#metadata for Seurat object, specifically transcriptotypes
mecp2Het.t.type = data.table(mecp2Het.CTX.HC.100vol.300counts.obj.Under1Over1.pred0.2.cellTypes[[c("rep", "t.type", "predicted.id")]], keep.rownames="index")

mecp2Het.t.type$exp_label = "NA"
mecp2Het.t.type[rep == "Rep1", exp_label := "Exp1"]
mecp2Het.t.type[rep == "Rep2", exp_label := "Exp2"]
mecp2Het.t.type[rep == "Rep6", exp_label := "Exp3"]
mecp2Het.t.type[rep == "Rep7", exp_label := "Exp4"]

mecp2Het.t.type$animal_label <- "B1"
mecp2Het.t.type[rep == "Rep1", animal_label := "B1"]
mecp2Het.t.type[rep == "Rep2", animal_label := "B2"]
mecp2Het.t.type[rep == "Rep6", animal_label := "B3"]
mecp2Het.t.type[rep == "Rep7", animal_label := "B3"]

mecp2Het.t.type[!(t.type %in% c("WT", "KO")), t.type:="NA"]


#summarize transcriptotype numbers by experiment
mecp2Het.t.type.exp.summ <- group_by(mecp2Het.t.type, exp_label) %>%
  summarise(
    number.WT = sum(t.type == "WT", na.rm = TRUE),
    number.KO = sum(t.type == "KO", na.rm = TRUE),
    number.NA = sum(is.na(t.type) | t.type == "NA")
  ) %>%
  mutate(
    total = number.WT + number.KO + number.NA,
    proportion.WT = number.WT / total,
    proportion.KO = number.KO / total,
    proportion.NA = number.NA / total
  ) %>% data.table

mecp2Het.t.type.exp.summ$animal_label <- "NA"
mecp2Het.t.type.exp.summ[exp_label == "Exp1", animal_label :="B1"]
mecp2Het.t.type.exp.summ[exp_label == "Exp2", animal_label :="B2"]
mecp2Het.t.type.exp.summ[exp_label == "Exp3", animal_label :="B3"]
mecp2Het.t.type.exp.summ[exp_label == "Exp4", animal_label :="B3"]

write.csv(mecp2Het.t.type.exp.summ, "HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/t.type.number.by.experiment.mecp2Het.CTX.HC.minVol100.maxVol3Med.minCount300.predScore0.2.Under1Over1.csv", row.names=F)

#summarize transcriptotype numbers by animal
mecp2Het.t.type.animal.summ <- group_by(mecp2Het.t.type, animal_label) %>%
  summarise(
    number.WT = sum(t.type == "WT", na.rm = TRUE),
    number.KO = sum(t.type == "KO", na.rm = TRUE),
    number.NA = sum(is.na(t.type) | t.type == "NA")
  ) %>%
  mutate(
    total = number.WT + number.KO + number.NA,
    proportion.WT = number.WT / total,
    proportion.KO = number.KO / total,
    proportion.NA = number.NA / total
  ) %>% data.table

mecp2Het.t.type.animal.summ$exp_label <- "NA"

mecp2Het.t.type.animal.summ[animal_label == "B1", exp_label := "Exp1"]
mecp2Het.t.type.animal.summ[animal_label == "B2", exp_label := "Exp2"]
mecp2Het.t.type.animal.summ[animal_label == "B3", exp_label := c("Exp3_Exp4")]

write.csv(mecp2Het.t.type.animal.summ, "HG_lab/Mati/GabelLab/Vizgen/mecp2_MERFISH/summary_tables/t.type.number.by.animal.mecp2Het.CTX.HC.minVol100.maxVol3Med.minCount300.predScore0.2.Under1Over1.csv", row.names=F)
