library(data.table)
library(dplyr)
library(ggplot2)

#global average methylation values
PV_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041642_summary.tsv")
PV_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041644_summary.tsv")

PV_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041641_summary.tsv")
PV_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041643_summary.tsv")

L5_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041645_summary.tsv")
L5_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041648_summary.tsv")

L5_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041646_summary.tsv")
L5_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041647_summary.tsv")

L4_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041634_summary.tsv")
L4_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042981_summary.tsv")

L4_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB041633_summary.tsv")
L4_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042982_summary.tsv")

SST_WT1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042979_summary.tsv")
SST_WT2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042980_summary.tsv")
SST_WT3 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB045574_22CGWFLT4_S263_summary.tsv")

SST_KO1 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB042978_summary.tsv")
SST_KO2 <- fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/LIB045575_22CGWFLT4_S264_summary.tsv")


#data table of global CA methylation values, background non-conversion rate subtracted
global_meth_dt <- rbind(
  cbind(global_meth_corrected = PV_WT1[(V1=="CA"), V2] - PV_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="PV", gtac_id="LIB041642"),
  cbind(global_meth_corrected = PV_WT2[(V1=="CA"), V2] - PV_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="PV", gtac_id="LIB041644"),
  cbind(global_meth_corrected = PV_KO1[(V1=="CA"), V2] - PV_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="PV", gtac_id="LIB041641"),
  cbind(global_meth_corrected = PV_KO2[(V1=="CA"), V2] - PV_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="PV", gtac_id="LIB041643"),
  cbind(global_meth_corrected = L5_WT1[(V1=="CA"), V2] - L5_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="L5", gtac_id="LIB041645"),
  cbind(global_meth_corrected = L5_WT2[(V1=="CA"), V2] - L5_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="L5", gtac_id="LIB041648"),
  cbind(global_meth_corrected = L5_KO1[(V1=="CA"), V2] - L5_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="L5", gtac_id="LIB041646"),
  cbind(global_meth_corrected = L5_KO2[(V1=="CA"), V2] - L5_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="L5", gtac_id="LIB041647"),
  cbind(global_meth_corrected = L4_WT1[(V1=="CA"), V2] - L4_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="L4", gtac_id="LIB041634"),
  cbind(global_meth_corrected = L4_WT2[(V1=="CA"), V2] - L4_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="L4", gtac_id="LIB042981"),
  cbind(global_meth_corrected = L4_KO1[(V1=="CA"), V2] - L4_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="L4", gtac_id="LIB041633"),
  cbind(global_meth_corrected = L4_KO2[(V1=="CA"), V2] - L4_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="L4", gtac_id="LIB042982"),
  cbind(global_meth_corrected = SST_WT1[(V1=="CA"), V2] - SST_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB042979"),
  cbind(global_meth_corrected = SST_WT2[(V1=="CA"), V2] - SST_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB042980"),
  cbind(global_meth_corrected = SST_WT3[(V1=="CA"), V2] - SST_WT3[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB045574"),
  cbind(global_meth_corrected = SST_KO1[(V1=="CA"), V2] - SST_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="SST", gtac_id="LIB042978"),
  cbind(global_meth_corrected = SST_KO2[(V1=="CA"), V2] - SST_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="SST", gtac_id="LIB045575")
) %>% data.table 

global_meth_dt <- global_meth_dt %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
global_meth_dt <- global_meth_dt %>% mutate(subclass = factor(subclass, levels=c("L4", "L5", "PV", "SST")))

ggplot(global_meth_dt, aes(x = genotype, y = as.numeric(global_meth_corrected), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.06)) +
  ylab("Global mean mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype", values = c("WT"="purple", "KO"="orange")) +
  #scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~subclass, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCA_backgroundSub_stripplot.png", width = 4, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCA_backgroundSub_stripplot.eps", width = 4, height = 6, dpi = 300, units = "in", device=cairo_ps)

wilcox.test(global_meth_dt[(genotype=="WT") & (subclass=="L4"), as.numeric(global_meth_corrected)], global_meth_dt[(genotype=="KO") & (subclass=="L4"), as.numeric(global_meth_corrected)])$p.value #p.value=0.3333333

wilcox.test(global_meth_dt[(genotype=="WT") & (subclass=="L5"), as.numeric(global_meth_corrected)], global_meth_dt[(genotype=="KO") & (subclass=="L5"), as.numeric(global_meth_corrected)])$p.value #p.value=0.6666667

wilcox.test(global_meth_dt[(genotype=="WT") & (subclass=="PV"), as.numeric(global_meth_corrected)], global_meth_dt[(genotype=="KO") & (subclass=="PV"), as.numeric(global_meth_corrected)])$p.value #p.value=0.3333333

wilcox.test(global_meth_dt[(genotype=="WT") & (subclass=="SST"), as.numeric(global_meth_corrected)], global_meth_dt[(genotype=="KO") & (subclass=="SST"), as.numeric(global_meth_corrected)])$p.value #p.value=0.8

#bigger
ggplot(global_meth_dt, aes(x = genotype, y = as.numeric(global_meth_corrected), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,0.06)) +
  ylab("Global mean mCA/CA") + xlab("") +
  scale_color_manual(name = "Genotype", values = c("WT"="purple", "KO"="orange")) +
  #scale_shape_manual(name = "Animal Label", values = c("B1"=16, "B2"=17, "B3"=15)) + # Customize shapes
  facet_grid(.~subclass, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCA_backgroundSub_stripplot_bigger.png", width = 3, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCA_backgroundSub_stripplot_bigger.eps", width = 3, height = 6, dpi = 300, units = "in", device=cairo_ps)

  

####mCG
global_mCG_dt <- rbind(
  cbind(global_mCG_corrected = PV_WT1[(V1=="CG"), V2] - PV_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="PV", gtac_id="LIB041642"),
  cbind(global_mCG_corrected = PV_WT2[(V1=="CG"), V2] - PV_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="PV", gtac_id="LIB041644"),
  cbind(global_mCG_corrected = PV_KO1[(V1=="CG"), V2] - PV_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="PV", gtac_id="LIB041641"),
  cbind(global_mCG_corrected = PV_KO2[(V1=="CG"), V2] - PV_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="PV", gtac_id="LIB041643"),
  cbind(global_mCG_corrected = L5_WT1[(V1=="CG"), V2] - L5_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="L5", gtac_id="LIB041645"),
  cbind(global_mCG_corrected = L5_WT2[(V1=="CG"), V2] - L5_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="L5", gtac_id="LIB041648"),
  cbind(global_mCG_corrected = L5_KO1[(V1=="CG"), V2] - L5_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="L5", gtac_id="LIB041646"),
  cbind(global_mCG_corrected = L5_KO2[(V1=="CG"), V2] - L5_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="L5", gtac_id="LIB041647"),
  cbind(global_mCG_corrected = L4_WT1[(V1=="CG"), V2] - L4_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="L4", gtac_id="LIB041634"),
  cbind(global_mCG_corrected = L4_WT2[(V1=="CG"), V2] - L4_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="L4", gtac_id="LIB042981"),
  cbind(global_mCG_corrected = L4_KO1[(V1=="CG"), V2] - L4_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="L4", gtac_id="LIB041633"),
  cbind(global_mCG_corrected = L4_KO2[(V1=="CG"), V2] - L4_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="L4", gtac_id="LIB042982"),
  cbind(global_mCG_corrected = SST_WT1[(V1=="CG"), V2] - SST_WT1[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB042979"),
  cbind(global_mCG_corrected = SST_WT2[(V1=="CG"), V2] - SST_WT2[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB042980"),
  cbind(global_mCG_corrected = SST_WT3[(V1=="CG"), V2] - SST_WT3[(V1=="Lambda"), V2], genotype="WT", subclass="SST", gtac_id="LIB045574"),
  cbind(global_mCG_corrected = SST_KO1[(V1=="CG"), V2] - SST_KO1[(V1=="Lambda"), V2], genotype="KO", subclass="SST", gtac_id="LIB042978"),
  cbind(global_mCG_corrected = SST_KO2[(V1=="CG"), V2] - SST_KO2[(V1=="Lambda"), V2], genotype="KO", subclass="SST", gtac_id="LIB045575")
) %>% data.table 

global_mCG_dt <- global_mCG_dt %>% mutate(genotype = factor(genotype, levels=c("WT", "KO")))
global_mCG_dt <- global_mCG_dt %>% mutate(subclass = factor(subclass, levels=c("L4", "L5", "PV", "SST")))


wilcox.test(global_mCG_dt[(genotype=="WT") & (subclass=="L4"), as.numeric(global_mCG_corrected)], global_mCG_dt[(genotype=="KO") & (subclass=="L4"), as.numeric(global_mCG_corrected)])$p.value #p.value=1

wilcox.test(global_mCG_dt[(genotype=="WT") & (subclass=="L5"), as.numeric(global_mCG_corrected)], global_mCG_dt[(genotype=="KO") & (subclass=="L5"), as.numeric(global_mCG_corrected)])$p.value #p.value=0.6666667

wilcox.test(global_mCG_dt[(genotype=="WT") & (subclass=="PV"), as.numeric(global_mCG_corrected)], global_mCG_dt[(genotype=="KO") & (subclass=="PV"), as.numeric(global_mCG_corrected)])$p.value #p.value=0.6666667

wilcox.test(global_mCG_dt[(genotype=="WT") & (subclass=="SST"), as.numeric(global_mCG_corrected)], global_mCG_dt[(genotype=="KO") & (subclass=="SST"), as.numeric(global_mCG_corrected)])$p.value #p.value=1

#bigger
ggplot(global_mCG_dt, aes(x = genotype, y = as.numeric(global_mCG_corrected), color=genotype)) +
  ggtitle("") +
  geom_jitter(width = 0.2, size = 5, alpha = 0.6) +
  coord_cartesian(ylim=c(0,1)) +
  ylab("Global mean mCG/CG") + xlab("") +
  scale_color_manual(name = "Genotype", values = c("WT"="purple", "KO"="orange")) +
  facet_grid(.~subclass, switch = "x", labeller = label_wrap_gen(width=8), scales = "free_x") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), strip.text.x = element_text(angle=0, size=8)) +
  theme(plot.title = element_text(hjust=0.5), legend.position = "bottom", legend.margin = margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(color="black"), axis.ticks.x = element_blank(), axis.text.y = element_text(size=10), axis.text.x = element_blank()) +
  guides(color = guide_legend(nrow=3, byrow=TRUE))
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCG_backgroundSub_stripplot_bigger.png", width = 3, height = 6, dpi = 300, units = "in", device='png')
ggsave(filename="HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/dotplots/L5_L4_PV_SST_INTACT_global_mCG_backgroundSub_stripplot_bigger.eps", width = 3, height = 6, dpi = 300, units = "in", device=cairo_ps)
