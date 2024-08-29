library(data.table)
library(dplyr)
library(ggplot2)
library(vioplot)
library(ggridges)
library(pheatmap)
library(viridis)
#Protein-coding genes
gene_coord = fread('HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed')
names(gene_coord)[1:6] = c("Chrom","Start","End","Gene","Num_transcripts","Strand")

#List of relevant chromosomes
chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

#gene body and flank methylation tables, using genic methylation from TSS+3kb to TES
PV_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L4_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
SST_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")
L5_INTACT_gene_body_TSSplus3kb_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt")

PV_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
L4_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
SST_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")
L5_INTACT_gene_body_TSSplus3kb_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt")

#replacing negative methylation values with zeros
PV_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
PV_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCA[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCA[flank_methylation_corrected < 0, flank_methylation_corrected := 0]

PV_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
PV_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
SST_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L4_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCG[gene_methylation_corrected < 0, gene_methylation_corrected := 0]
L5_INTACT_gene_body_TSSplus3kb_flank_mCG[flank_methylation_corrected < 0, flank_methylation_corrected := 0]



#cCREs
mousebrain_union_cCREs_PV_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_PV_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_SST_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_SST_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L4_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_L4_mCA[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L5_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L5_WT_KO_deep_INTACT_mCA_mm9.bed")
mousebrain_union_cCREs_L5_mCA[, cCRE_methylation := V5/V6]

mousebrain_union_cCREs_PV_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_PV_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_SST_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_SST_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L4_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_L4_mCG[, cCRE_methylation := V5/V6]
mousebrain_union_cCREs_L5_mCG = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/INTACT_meth_cCREs/mousebrain_union_cCRE_L5_WT_KO_deep_INTACT_mCG_mm9.bed")
mousebrain_union_cCREs_L5_mCG[, cCRE_methylation := V5/V6]

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

mousebrain_union_cCREs_PV_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCA$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_SST_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_SST_mCA$cCRE_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L4_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_L4_mCA$cCRE_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L5_mCA$cCRE_methylation_corrected <- mousebrain_union_cCREs_L5_mCA$cCRE_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

mousebrain_union_cCREs_PV_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_PV_mCG$cCRE_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_SST_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_SST_mCG$cCRE_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L4_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_L4_mCG$cCRE_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
mousebrain_union_cCREs_L5_mCG$cCRE_methylation_corrected <- mousebrain_union_cCREs_L5_mCG$cCRE_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

#turning negative methylation values into zeros
mousebrain_union_cCREs_PV_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_SST_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L4_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L5_mCA[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]

mousebrain_union_cCREs_PV_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_SST_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L4_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]
mousebrain_union_cCREs_L5_mCG[cCRE_methylation_corrected < 0, cCRE_methylation_corrected := 0]



mPv_mL4_snmcseq_gene_and_flank = data.table(full_join(x=mPv_snmcseq_gene_and_flank, y=mL4_snmcseq_gene_and_flank, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mPv_mL4_snmcseq_gene_and_flank[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mPv_mL4_snmcseq_gene_and_flank[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]

mPv_mL4_snmcseq_gene_and_flank_mCG = data.table(full_join(x=mPv_snmcseq_gene_and_flank_mCG, y=mL4_snmcseq_gene_and_flank_mCG, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mPv_mL4_snmcseq_gene_and_flank_mCG[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mPv_mL4_snmcseq_gene_and_flank_mCG[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]


#Calculate z-score of mL4-mPv methylation difference of gene bodies and flanks
mPv_mL4_snmcseq_gene_and_flank[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mPv_mL4_snmcseq_gene_and_flank[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mPv_mL4_snmcseq_gene_and_flank = data.table(mPv_mL4_snmcseq_gene_and_flank%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mPv_mL4_snmcseq_gene_and_flank = data.table(mPv_mL4_snmcseq_gene_and_flank%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

mPv_mL4_snmcseq_gene_and_flank_mCG[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mPv_mL4_snmcseq_gene_and_flank_mCG[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mPv_mL4_snmcseq_gene_and_flank_mCG = data.table(mPv_mL4_snmcseq_gene_and_flank_mCG%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mPv_mL4_snmcseq_gene_and_flank_mCG = data.table(mPv_mL4_snmcseq_gene_and_flank_mCG%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

#Combining mPv and mL4 cCRE mCA/CA or mCG/CG tables
mPv_mL4_cCRE_mCA = full_join(x=mousebrain_union_cCREs_mPv_mCA, y=mousebrain_union_cCREs_mL4_mCA, by=c("V1", "V2", "V3", "V4"))
mPv_mL4_cCRE_mCG = full_join(x=mousebrain_union_cCREs_mPv_mCG, y=mousebrain_union_cCREs_mL4_mCG, by=c("V1", "V2", "V3", "V4"))

#cCRE mCA/CA fold differences
mPv_mL4_cCRE_mCA$cCRE_meth_fd =  mPv_mL4_cCRE_mCA[, .(cCRE_methylation.x)]/mPv_mL4_cCRE_mCA[, .(cCRE_methylation.y)]
mPv_mL4_cCRE_mCG$cCRE_meth_fd =  mPv_mL4_cCRE_mCG[, .(cCRE_methylation.x)]/mPv_mL4_cCRE_mCG[, .(cCRE_methylation.y)]

#cCRE file with genic location of cCRE
mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_Cicero_linked_genes_mm9_genicBooleans.txt")
mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9 = unique(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene==1), .(cCRE_chrom, cCRE_start, cCRE_end, cCRE_label)])
#write.table(mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9, file="HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/gene_cCRE_positively_correlated_connections/mousebrain_union_nonPromoter_cCREs_intragenicLinked_genes_cCREcoords_mm9.bed", quote=F, row.names=F, col.names=F, sep="\t")

mPv_mL4_snmcseq_gene_and_flank_bothFinite =  data.table(mPv_mL4_snmcseq_gene_and_flank[is.finite(gene_meth_fd) & is.finite(flank_meth_fd), ] %>% tidyr::gather("element", "value", c("gene_meth_fd","flank_meth_fd")))
ggplot(data=mPv_mL4_snmcseq_gene_and_flank_bothFinite, aes(x=element, y=value))+
  geom_violin()+
  ylab("(mPv mCA/CA)/(mL4 mCA/CA)") + 
  coord_cartesian(ylim=c(0,20))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

mPv_mL4_snmcseq_gene_flank_cCRE_mCA = data.table(rbind(
  cbind(val=mPv_mL4_snmcseq_gene_and_flank[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mPv_mL4_snmcseq_gene_and_flank[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mPv_mL4_cCRE_mCA[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mPv_mL4_snmcseq_gene_flank_cCRE_mCA[, val := as.numeric(val)]

mPv_mL4_snmcseq_gene_flank_cCRE_mCG = data.table(rbind(
  cbind(val=mPv_mL4_snmcseq_gene_and_flank_mCG[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mPv_mL4_snmcseq_gene_and_flank_mCG[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mPv_mL4_cCRE_mCG[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mPv_mL4_snmcseq_gene_flank_cCRE_mCG[, val := as.numeric(val)]

ggplot(mPv_mL4_snmcseq_gene_flank_cCRE_mCA, aes(x=val)) + 
  geom_line(stat="function", fun=dnorm, args = list(mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="Regional"), val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="Regional"), val])), mapping = aes(colour="Regional"), show.legend = TRUE) + 
  geom_line(stat="function", fun=dnorm, args = list(mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="Gene body"), val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="Gene body"), val])), mapping = aes(colour="Gene body"), show.legend = TRUE) + 
  geom_line(stat="function", fun=dnorm, args = list(mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="cCRE"), val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[(element=="cCRE"), val])), mapping = aes(colour="cCRE"), show.legend = TRUE) +
  scale_colour_manual(limits = c("Regional", "Gene body", "cCRE"), labels=c("Regional", "Gene body", "cCRE"), values = c("black", "darkred", "mediumpurple2"))+
  xlab("(mPv mCA/CA)/(mL4 mCA/CA)") + ylab("Probability density") + 
  coord_cartesian(xlim=c(0,15), ylim=c(0,1.4))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_cCRE_mCAperCA_foldDiff_pdf.eps")
plot(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val], dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="darkred", xlim=c(0,12), ylim=c(0, 1.4), pch=16)
par(new=TRUE)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="black", xlim=c(0,12), ylim=c(0,1.4), pch=16)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="mediumpurple3", xlim=c(0,12), ylim=c(0,1.4), pch=16)
legend("topright", legend=c('Regional', 'Gene body', 'cCRE'), cex=0.8, pt.cex=1,
       col=c('darkred', 'black', 'mediumpurple3'), pch=c(16,16, 16), bty="n")
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_cCRE_mCAperCA_foldDiff_pdf.png")
plot(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val], dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Regional", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="darkred", xlim=c(0,12), ylim=c(0, 1.4), pch=16)
par(new=TRUE)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="Gene body", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="black", xlim=c(0,12), ylim=c(0,1.4), pch=16)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCA[element=="cCRE", val])), xlab="(mPv mCA/CA)/(mL4 mCA/CA)", ylab="Probability density", col="mediumpurple3", xlim=c(0,12), ylim=c(0,1.4), pch=16)
legend("topright", legend=c('Regional', 'Gene body', 'cCRE'), cex=0.8, pt.cex=1,
       col=c('darkred', 'black', 'mediumpurple3'), pch=c(16,16, 16), bty="n")
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_cCRE_mCGperCG_foldDiff_pdf.eps")
plot(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val], dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="darkred", xlim=c(0,5), ylim=c(0, 6), pch=16)
par(new=TRUE)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="black", xlim=c(0,5), ylim=c(0,6), pch=16)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="mediumpurple3", xlim=c(0,5), ylim=c(0,6), pch=16)
legend("topright", legend=c('Regional', 'Gene body', 'cCRE'), cex=0.8, pt.cex=1,
       col=c('darkred', 'black', 'mediumpurple3'), pch=c(16,16, 16), bty="n")
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_cCRE_mCGperCG_foldDiff_pdf.png")
plot(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val], dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Regional", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="darkred", xlim=c(0,5), ylim=c(0, 6), pch=16)
par(new=TRUE)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="Gene body", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="black", xlim=c(0,5), ylim=c(0,6), pch=16)
points(x=mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val], y=dnorm(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val], mean=mean(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val]), sd=sd(mPv_mL4_snmcseq_gene_flank_cCRE_mCG[element=="cCRE", val])), xlab="(mPv mCG/CG)/(mL4 mCG/CG)", ylab="Probability density", col="mediumpurple3", xlim=c(0,5), ylim=c(0,6), pch=16)
legend("topright", legend=c('Regional', 'Gene body', 'cCRE'), cex=0.8, pt.cex=1,
       col=c('darkred', 'black', 'mediumpurple3'), pch=c(16,16, 16), bty="n")
dev.off()



setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_allLinkedcCRE_mCAperCA_foldDiff_violinPlot.eps")
with(mPv_mL4_snmcseq_gene_flank_cCRE_mCA , vioplot( 
  val[element=="Regional"] , val[element=="Gene body"], val[element=="cCRE"],  
  names=c("Regional","Gene body","cCRE"),
  ylab="(mPv mCA/CA)/(mL4 mCA/CA)"
))
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_allLinkedcCRE_mCAperCA_foldDiff_violinPlot.png")
with(mPv_mL4_snmcseq_gene_flank_cCRE_mCA , vioplot( 
  val[element=="Regional"] , val[element=="Gene body"], val[element=="cCRE"],  
  names=c("Regional","Gene body","cCRE"),
  ylab="(mPv mCA/CA)/(mL4 mCA/CA)"
))
dev.off()

png("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_regional_geneBody_allLinkedcCRE_mCGperCG_foldDiff_violinPlot.png")
with(mPv_mL4_snmcseq_gene_flank_cCRE_mCG , vioplot( 
  val[element=="Regional"] , val[element=="Gene body"], val[element=="cCRE"],  
  names=c("Regional","Gene body","cCRE"),
  ylab="(mPv mCG/CG)/(mL4 mCG/CG)"
))
dev.off()



mPv_mL4_cCRE_mCA[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mPv_mL4_cCRE_mCA = data.table(mPv_mL4_cCRE_mCA%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))
mPv_mL4_cCRE_mCG[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mPv_mL4_cCRE_mCG = data.table(mPv_mL4_cCRE_mCG%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))

pv_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

sst_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
sst_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
sst_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

L4_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L4_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L4_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")

L5_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_unchanged_genes_p0.5_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed")
L5_otherCellType_mr_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")
L5_otherCellType_ma_genes_q0.1_nondedup_mm9 = fread("HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed")


pv_mr_L4_unchanged_genes = intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], L4_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])
L4_mr_pv_unchanged_genes = intersect(L4_mr_genes_q0.1_nondedup_mm9[,V4], pv_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])

View(mPv_mL4_cCRE_mCA[(zScore_cCRE_meth_diff>=2) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_mr_L4_unchanged_genes) & (Intragenic_to_linked_gene==1),cCRE_label]),])

View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_mr_L4_unchanged_genes) &
                                                                   (Gene %in% mPv_mL4_snmcseq_gene_and_flank[zScore_gene_meth_diff>=1, gene]) & 
                                                                   (Intragenic_to_linked_gene==1) &
                                                                   (cCRE_label %in% mPv_mL4_cCRE_mCA[(zScore_cCRE_meth_diff>=1), V4]),])

mPv_10kb_mCA = fread("HG_lab/Mati/GabelLab/Luo2017_data/mPv_methylation/mPv_Luo2017_snmcseq_CA_merged_mm9_10kb_windows.bed")
mSst_10kb_mCA = fread("HG_lab/Mati/GabelLab/Luo2017_data/mSst_all_methylation/mSst_all_Luo2017_snmcseq_CA_merged_mm9_10kb_windows.bed")
mL4_10kb_mCA = fread("HG_lab/Mati/GabelLab/Luo2017_data/mL4_methylation/mL4_Luo2017_snmcseq_CA_merged_mm9_10kb_windows.bed")
mL5_10kb_mCA = fread("HG_lab/Mati/GabelLab/Luo2017_data/mL5_all_methylation/mL5_all_Luo2017_snmcseq_CA_merged_mm9_10kb_windows.bed")

names(mPv_10kb_mCA) = c("chrom", "start", "end", "Pv_meth", "Pv_cov")
mPv_10kb_mCA[, Pv_unmeth := Pv_cov - Pv_meth]
mPv_10kb_mCA[, Pv_methylation := Pv_meth/Pv_cov]
mPv_10kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
mPv_10kb_mCA = mPv_10kb_mCA[chrom %in% chrom_list]

names(mSst_10kb_mCA) = c("chrom", "start", "end", "Sst_meth", "Sst_cov")
mSst_10kb_mCA[, Sst_unmeth := Sst_cov - Sst_meth]
mSst_10kb_mCA[, Sst_methylation := Sst_meth/Sst_cov]
mSst_10kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
mSst_10kb_mCA = mSst_10kb_mCA[chrom %in% chrom_list]

names(mL4_10kb_mCA) = c("chrom", "start", "end", "L4_meth", "L4_cov")
mL4_10kb_mCA[, L4_unmeth := L4_cov - L4_meth]
mL4_10kb_mCA[, L4_methylation := L4_meth/L4_cov]
mL4_10kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
mL4_10kb_mCA = mL4_10kb_mCA[chrom %in% chrom_list]

names(mL5_10kb_mCA) = c("chrom", "start", "end", "L5_meth", "L5_cov")
mL5_10kb_mCA[, L5_unmeth := L5_cov - L5_meth]
mL5_10kb_mCA[, L5_methylation := L5_meth/L5_cov]
mL5_10kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
mL5_10kb_mCA = mL5_10kb_mCA[chrom %in% chrom_list]

mm9_10kb_windows_genicRegions = fread("HG_lab/Mati/GabelLab/databases/mm9_10kb_windows_genicRegions.bed")

mPv_mSst_10kb_mCA = unique(data.table(inner_join(x=mPv_10kb_mCA, y=mSst_10kb_mCA, by=c("chrom","start","end","coord"))))
mPv_mSst_mL4_10kb_mCA = unique(data.table(inner_join(x=mPv_mSst_10kb_mCA, y=mL4_10kb_mCA, by=c("chrom","start","end","coord"))))
mPv_mSst_mL4_mL5_10kb_mCA = unique(data.table(inner_join(x=mPv_mSst_mL4_10kb_mCA, y=mL5_10kb_mCA, by=c("chrom","start","end","coord"))))


mPv_mSst_mL4_mL5_10kb_mCA[(coord %in% mm9_10kb_windows_genicRegions[, V4]), genic_location := "Intragenic"]
mPv_mSst_mL4_mL5_10kb_mCA[!(coord %in% mm9_10kb_windows_genicRegions[, V4]), genic_location := "Extragenic"]

mPv_mL4_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Pv_methylation/L4_methylation], element="Extragenic regions"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Pv_methylation/L4_methylation], element="Intragenic regions")))
mPv_mL4_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

with(mPv_mL4_10kb_windows_mCA_methFoldDiff[is.finite(val)] , vioplot( 
  val[element=="Extragenic regions"] , val[element=="Intragenic regions"],  
  names=c("Extragenic regions", "Intragenic regions"),
  ylab="(mPv mCA/CA)/(mL4 mCA/CA)"
))
#10kb windows centered on nonpromoter cCREs
mPv_10kb_cCRE_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_10kbWindows/mousebrain_union_nonPromoter_cCREs_10kbWindows_mPv_Luo2017_snmcseq_CA_merged_mm9.bed")
mSst_10kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_10kbWindows/mousebrain_union_nonPromoter_cCREs_10kbWindows_mSst_all_Luo2017_snmcseq_CA_merged_mm9.bed")
mL4_10kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_10kbWindows/mousebrain_union_nonPromoter_cCREs_10kbWindows_mL4_Luo2017_snmcseq_CA_merged_mm9.bed")
mL5_10kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_10kbWindows/mousebrain_union_nonPromoter_cCREs_10kbWindows_mL5_all_Luo2017_snmcseq_CA_merged_mm9.bed")

names(mPv_10kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "Pv_meth", "Pv_cov")
names(mSst_10kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "Sst_meth", "Sst_cov")
names(mL4_10kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "L4_meth", "L4_cov")
names(mL5_10kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "L5_meth", "L5_cov")

mPv_10kb_cCRE_mCA[, Pv_methylation := Pv_meth/Pv_cov]
mSst_10kb_cCRE_mCA[, Sst_methylation := Sst_meth/Sst_cov]
mL4_10kb_cCRE_mCA[, L4_methylation := L4_meth/L4_cov]
mL5_10kb_cCRE_mCA[, L5_methylation := L5_meth/L5_cov]

mPv_mSst_10kb_cCRE_mCA = unique(data.table(inner_join(x=mPv_10kb_cCRE_mCA, y=mSst_10kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))
mPv_mSst_mL4_10kb_cCRE_mCA = unique(data.table(inner_join(x=mPv_mSst_10kb_cCRE_mCA, y=mL4_10kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA = unique(data.table(inner_join(x=mPv_mSst_mL4_10kb_cCRE_mCA, y=mL5_10kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))

#intragenic non-promoter cCREs
intragenic_nonPromoter_cCREs_mm9 =fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/mousebrain_union_nonPromoter_cCREs_ensgene_mm9.txt")

mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = mPv_mSst_mL4_mL5_10kb_cCRE_mCA[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label],]

mPv_mL4_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Pv_methylation/L4_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Pv_methylation/L4_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Pv_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mPv_mL4_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/mPv_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_violinPlot.eps")
with(mPv_mL4_10kb_windows_mCA_methFoldDiff[is.finite(val)] , vioplot( 
  val[element=="Extragenic"] , val[element=="Intragenic"], val[element=="Intragenic,\n cCRE-centered"], 
  names=c("Extragenic", "Intragenic", "Intragenic,\n cCRE-centered"),
  xlab="10kb regions", ylab="(mPv mCA/CA)/(mL4 mCA/CA)"
))


##
mPv_mL5_snmcseq_gene_and_flank = data.table(full_join(x=mPv_snmcseq_gene_and_flank, y=mL5_snmcseq_gene_and_flank, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mPv_mL5_snmcseq_gene_and_flank[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mPv_mL5_snmcseq_gene_and_flank[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]

mPv_mL5_snmcseq_gene_and_flank_mCG = data.table(full_join(x=mPv_snmcseq_gene_and_flank_mCG, y=mL5_snmcseq_gene_and_flank_mCG, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mPv_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mPv_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]


#Calculate z-score of mPv-mL5 methylation difference of gene bodies and flanks
mPv_mL5_snmcseq_gene_and_flank[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mPv_mL5_snmcseq_gene_and_flank[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mPv_mL5_snmcseq_gene_and_flank = data.table(mPv_mL5_snmcseq_gene_and_flank%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mPv_mL5_snmcseq_gene_and_flank = data.table(mPv_mL5_snmcseq_gene_and_flank%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

mPv_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mPv_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mPv_mL5_snmcseq_gene_and_flank_mCG = data.table(mPv_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mPv_mL5_snmcseq_gene_and_flank_mCG = data.table(mPv_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

#Combining mPv and mL4 cCRE mCA/CA or mCG/CG tables
mPv_mL5_cCRE_mCA = full_join(x=mousebrain_union_cCREs_mPv_mCA, y=mousebrain_union_cCREs_mL5_mCA, by=c("V1", "V2", "V3", "V4"))
mPv_mL5_cCRE_mCG = full_join(x=mousebrain_union_cCREs_mPv_mCG, y=mousebrain_union_cCREs_mL5_mCG, by=c("V1", "V2", "V3", "V4"))

#cCRE mCA/CA fold differences
mPv_mL5_cCRE_mCA$cCRE_meth_fd =  mPv_mL5_cCRE_mCA[, .(cCRE_methylation.x)]/mPv_mL5_cCRE_mCA[, .(cCRE_methylation.y)]
mPv_mL5_cCRE_mCG$cCRE_meth_fd =  mPv_mL5_cCRE_mCG[, .(cCRE_methylation.x)]/mPv_mL5_cCRE_mCG[, .(cCRE_methylation.y)]

mPv_mL5_cCRE_mCA[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mPv_mL5_cCRE_mCA = data.table(mPv_mL5_cCRE_mCA%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))
mPv_mL5_cCRE_mCG[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mPv_mL5_cCRE_mCG = data.table(mPv_mL5_cCRE_mCG%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))


mPv_mL5_snmcseq_gene_flank_cCRE_mCA = data.table(rbind(
  cbind(val=mPv_mL5_snmcseq_gene_and_flank[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mPv_mL5_snmcseq_gene_and_flank[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mPv_mL5_cCRE_mCA[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mPv_mL5_snmcseq_gene_flank_cCRE_mCA[, val := as.numeric(val)]

mPv_mL5_snmcseq_gene_flank_cCRE_mCG = data.table(rbind(
  cbind(val=mPv_mL5_snmcseq_gene_and_flank_mCG[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mPv_mL5_snmcseq_gene_and_flank_mCG[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mPv_mL5_cCRE_mCG[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mPv_mL5_snmcseq_gene_flank_cCRE_mCG[, val := as.numeric(val)]

##Sst vs L5
mSst_mL5_snmcseq_gene_and_flank = data.table(full_join(x=mSst_snmcseq_gene_and_flank, y=mL5_snmcseq_gene_and_flank, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mSst_mL5_snmcseq_gene_and_flank[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mSst_mL5_snmcseq_gene_and_flank[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]

mSst_mL5_snmcseq_gene_and_flank_mCG = data.table(full_join(x=mSst_snmcseq_gene_and_flank_mCG, y=mL5_snmcseq_gene_and_flank_mCG, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mSst_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mSst_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]


#Calculate z-score of mSst-mL5 methylation difference of gene bodies and flanks
mSst_mL5_snmcseq_gene_and_flank[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mSst_mL5_snmcseq_gene_and_flank[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mSst_mL5_snmcseq_gene_and_flank = data.table(mSst_mL5_snmcseq_gene_and_flank%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mSst_mL5_snmcseq_gene_and_flank = data.table(mSst_mL5_snmcseq_gene_and_flank%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

mSst_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mSst_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mSst_mL5_snmcseq_gene_and_flank_mCG = data.table(mSst_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mSst_mL5_snmcseq_gene_and_flank_mCG = data.table(mSst_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

#Combining mSst and mL5 cCRE mCA/CA or mCG/CG tables
mSst_mL5_cCRE_mCA = full_join(x=mousebrain_union_cCREs_mSst_mCA, y=mousebrain_union_cCREs_mL5_mCA, by=c("V1", "V2", "V3", "V4"))
mSst_mL5_cCRE_mCG = full_join(x=mousebrain_union_cCREs_mSst_mCG, y=mousebrain_union_cCREs_mL5_mCG, by=c("V1", "V2", "V3", "V4"))

#cCRE mCA/CA fold differences
mSst_mL5_cCRE_mCA$cCRE_meth_fd =  mSst_mL5_cCRE_mCA[, .(cCRE_methylation.x)]/mSst_mL5_cCRE_mCA[, .(cCRE_methylation.y)]
mSst_mL5_cCRE_mCG$cCRE_meth_fd =  mSst_mL5_cCRE_mCG[, .(cCRE_methylation.x)]/mSst_mL5_cCRE_mCG[, .(cCRE_methylation.y)]

mSst_mL5_cCRE_mCA[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mSst_mL5_cCRE_mCA = data.table(mSst_mL5_cCRE_mCA%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))
mSst_mL5_cCRE_mCG[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mSst_mL5_cCRE_mCG = data.table(mSst_mL5_cCRE_mCG%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))


mSst_mL5_snmcseq_gene_flank_cCRE_mCA = data.table(rbind(
  cbind(val=mSst_mL5_snmcseq_gene_and_flank[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mSst_mL5_snmcseq_gene_and_flank[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mSst_mL5_cCRE_mCA[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mSst_mL5_snmcseq_gene_flank_cCRE_mCA[, val := as.numeric(val)]

mSst_mL5_snmcseq_gene_flank_cCRE_mCG = data.table(rbind(
  cbind(val=mSst_mL5_snmcseq_gene_and_flank_mCG[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mSst_mL5_snmcseq_gene_and_flank_mCG[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mSst_mL5_cCRE_mCG[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mSst_mL5_snmcseq_gene_flank_cCRE_mCG[, val := as.numeric(val)]

##L4 vs L5
mL4_mL5_snmcseq_gene_and_flank = data.table(full_join(x=mL4_snmcseq_gene_and_flank, y=mL5_snmcseq_gene_and_flank, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mL4_mL5_snmcseq_gene_and_flank[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mL4_mL5_snmcseq_gene_and_flank[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]

mL4_mL5_snmcseq_gene_and_flank_mCG = data.table(full_join(x=mL4_snmcseq_gene_and_flank_mCG, y=mL5_snmcseq_gene_and_flank_mCG, by=c('gene', 'chrom', 'start', 'end', 'transcripts', 'strand')))
mL4_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_fd := as.numeric(gene_methylation.x)/as.numeric(gene_methylation.y)]
mL4_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_fd := as.numeric(flank_methylation.x)/as.numeric(flank_methylation.y)]


#Calculate z-score of mL4-mL5 methylation difference of gene bodies and flanks
mL4_mL5_snmcseq_gene_and_flank[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mL4_mL5_snmcseq_gene_and_flank[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mL4_mL5_snmcseq_gene_and_flank = data.table(mL4_mL5_snmcseq_gene_and_flank%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mL4_mL5_snmcseq_gene_and_flank = data.table(mL4_mL5_snmcseq_gene_and_flank%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

mL4_mL5_snmcseq_gene_and_flank_mCG[, gene_meth_diff := gene_methylation.x - gene_methylation.y]
mL4_mL5_snmcseq_gene_and_flank_mCG[, flank_meth_diff := flank_methylation.x - flank_methylation.y]
mL4_mL5_snmcseq_gene_and_flank_mCG = data.table(mL4_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_gene_meth_diff = (gene_meth_diff - mean(gene_meth_diff, na.rm=TRUE))/sd(gene_meth_diff, na.rm=TRUE)))
mL4_mL5_snmcseq_gene_and_flank_mCG = data.table(mL4_mL5_snmcseq_gene_and_flank_mCG%>% mutate(zScore_flank_meth_diff = (flank_meth_diff - mean(flank_meth_diff, na.rm=TRUE))/sd(flank_meth_diff, na.rm=TRUE)))

#Combining mL4 and mL5 cCRE mCA/CA or mCG/CG tables
mL4_mL5_cCRE_mCA = full_join(x=mousebrain_union_cCREs_mL4_mCA, y=mousebrain_union_cCREs_mL5_mCA, by=c("V1", "V2", "V3", "V4"))
mL4_mL5_cCRE_mCG = full_join(x=mousebrain_union_cCREs_mL4_mCG, y=mousebrain_union_cCREs_mL5_mCG, by=c("V1", "V2", "V3", "V4"))

#cCRE mCA/CA fold differences
mL4_mL5_cCRE_mCA$cCRE_meth_fd =  mL4_mL5_cCRE_mCA[, .(cCRE_methylation.x)]/mL4_mL5_cCRE_mCA[, .(cCRE_methylation.y)]
mL4_mL5_cCRE_mCG$cCRE_meth_fd =  mL4_mL5_cCRE_mCG[, .(cCRE_methylation.x)]/mL4_mL5_cCRE_mCG[, .(cCRE_methylation.y)]

mL4_mL5_cCRE_mCA[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mL4_mL5_cCRE_mCA = data.table(mL4_mL5_cCRE_mCA%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))
mL4_mL5_cCRE_mCG[, cCRE_meth_diff := cCRE_methylation.x - cCRE_methylation.y]
mL4_mL5_cCRE_mCG = data.table(mL4_mL5_cCRE_mCG%>% mutate(zScore_cCRE_meth_diff = (cCRE_meth_diff - mean(cCRE_meth_diff, na.rm=TRUE))/sd(cCRE_meth_diff, na.rm=TRUE)))


mL4_mL5_snmcseq_gene_flank_cCRE_mCA = data.table(rbind(
  cbind(val=mL4_mL5_snmcseq_gene_and_flank[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mL4_mL5_snmcseq_gene_and_flank[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mL4_mL5_cCRE_mCA[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mL4_mL5_snmcseq_gene_flank_cCRE_mCA[, val := as.numeric(val)]

mL4_mL5_snmcseq_gene_flank_cCRE_mCG = data.table(rbind(
  cbind(val=mL4_mL5_snmcseq_gene_and_flank_mCG[is.finite(flank_meth_fd), (flank_meth_fd)], element="Regional"),
  cbind(val=mL4_mL5_snmcseq_gene_and_flank_mCG[is.finite(gene_meth_fd), (gene_meth_fd)], element="Gene body"),
  cbind(val=mL4_mL5_cCRE_mCG[is.finite(cCRE_meth_fd) & (V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[, cCRE_label]), cCRE_meth_fd], element="cCRE")))
mL4_mL5_snmcseq_gene_flank_cCRE_mCG[, val := as.numeric(val)]




#genes that are unchanged in L5 but MR in other cell types
pv_mr_L5_unchanged_genes = intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], L5_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])
sst_mr_L5_unchanged_genes = intersect(sst_mr_genes_q0.1_nondedup_mm9[,V4], L5_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])
L4_mr_L5_unchanged_genes = intersect(L4_mr_genes_q0.1_nondedup_mm9[,V4], L5_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])

pv_sst_L4_mr_L5_unchanged_genes = intersect(intersect(pv_mr_L5_unchanged_genes, sst_mr_L5_unchanged_genes), L4_mr_L5_unchanged_genes)

pv_sst_L4_mr_L5_unchanged_subthreshold_genes = setdiff(intersect(intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], sst_mr_genes_q0.1_nondedup_mm9[,V4]), L4_mr_genes_q0.1_nondedup_mm9[,V4]), L5_mr_genes_q0.1_nondedup_mm9[,V4])

View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_mr_L5_unchanged_genes) &
                                                                     (Gene %in% mPv_mL5_snmcseq_gene_and_flank[zScore_gene_meth_diff>=1, gene]) & 
                                                                     (Intragenic_to_linked_gene==1) &
                                                                     (cCRE_label %in% mPv_mL5_cCRE_mCA[(zScore_cCRE_meth_diff>=1), V4]),])

View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_sst_L4_mr_L5_unchanged_subthreshold_genes) &
                                                                     (Intragenic_to_linked_gene==1) &
                                                                     (cCRE_label %in% mPv_mL5_cCRE_mCA[(zScore_cCRE_meth_diff>=1), V4]),])

mPv_mL5_snmcseq_gene_and_flank[gene %in% pv_sst_L4_mr_L5_unchanged_subthreshold_genes, .(gene, zScore_gene_meth_diff, zScore_flank_meth_diff)]
View(mPv_mL5_cCRE_mCA[(V4 %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[Gene %in% pv_sst_L4_mr_L5_unchanged_subthreshold_genes, cCRE_label]) & (zScore_cCRE_meth_diff>=1),])



pv_L4_L5_mr_sst_unchanged_genes = intersect(intersect(intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], L4_mr_genes_q0.1_nondedup_mm9[,V4]), L5_mr_genes_q0.1_nondedup_mm9[,V4]), sst_otherCellType_mr_genes_q0.1_nondedup_mm9[,V4])

mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = cbind(mousebrain_union_cCREs_mPv_mCA[,.(V1,V2,V3,V4,cCRE_methylation)], mousebrain_union_cCREs_mSst_mCA[,.(cCRE_methylation)], mousebrain_union_cCREs_mL4_mCA[,.(cCRE_methylation)], mousebrain_union_cCREs_mL5_mCA[,.(cCRE_methylation)])
names(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA) = c("cCRE_chrom", "cCRE_start", "cCRE_end", "cCRE_label", "Pv_cCRE_methylation", "Sst_cCRE_methylation", "L4_cCRE_methylation", "L5_cCRE_methylation")
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, Sst_Pv_cCRE_meth_diff := Sst_cCRE_methylation - Pv_cCRE_methylation]
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, Sst_L5_cCRE_meth_diff := Sst_cCRE_methylation - L5_cCRE_methylation]
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, Sst_L4_cCRE_meth_diff := Sst_cCRE_methylation - L4_cCRE_methylation]
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, Pv_L4_cCRE_meth_diff := Pv_cCRE_methylation - L4_cCRE_methylation]
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, Pv_L5_cCRE_meth_diff := Pv_cCRE_methylation - L5_cCRE_methylation]
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[, L5_L4_cCRE_meth_diff := L5_cCRE_methylation - L4_cCRE_methylation]

mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_Sst_Pv_cCRE_meth_diff = (Sst_Pv_cCRE_meth_diff - mean(Sst_Pv_cCRE_meth_diff, na.rm=TRUE))/sd(Sst_Pv_cCRE_meth_diff, na.rm=TRUE)))
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_Sst_L5_cCRE_meth_diff = (Sst_L5_cCRE_meth_diff - mean(Sst_L5_cCRE_meth_diff, na.rm=TRUE))/sd(Sst_L5_cCRE_meth_diff, na.rm=TRUE)))
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_Sst_L4_cCRE_meth_diff = (Sst_L4_cCRE_meth_diff - mean(Sst_L4_cCRE_meth_diff, na.rm=TRUE))/sd(Sst_L4_cCRE_meth_diff, na.rm=TRUE)))
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_Pv_L4_cCRE_meth_diff = (Pv_L4_cCRE_meth_diff - mean(Pv_L4_cCRE_meth_diff, na.rm=TRUE))/sd(Pv_L4_cCRE_meth_diff, na.rm=TRUE)))
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_Pv_L5_cCRE_meth_diff = (Pv_L5_cCRE_meth_diff - mean(Pv_L5_cCRE_meth_diff, na.rm=TRUE))/sd(Pv_L5_cCRE_meth_diff, na.rm=TRUE)))
mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA = data.table(mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA %>% mutate(zScore_L5_L4_cCRE_meth_diff = (L5_L4_cCRE_meth_diff - mean(L5_L4_cCRE_meth_diff, na.rm=TRUE))/sd(L5_L4_cCRE_meth_diff, na.rm=TRUE)))


View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_L4_L5_mr_sst_unchanged_genes) &
                                                                     (Intragenic_to_linked_gene==1) &
                                                                     (cCRE_label %in% mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[(zScore_Sst_Pv_cCRE_meth_diff>=1), cCRE_label]),])

View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_L4_L5_mr_sst_unchanged_genes) &
                                                                     (Intragenic_to_linked_gene==1) &
                                                                     (cCRE_label %in% mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[(zScore_Sst_Pv_cCRE_meth_diff<=-1) & (zScore_Sst_L4_cCRE_meth_diff<=-1) & (zScore_Sst_L5_cCRE_meth_diff<=-1), cCRE_label]),])


mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = cbind(mPv_snmcseq_gene_and_flank[, .(chrom, start, end, gene, gene_methylation, flank_methylation)], mSst_snmcseq_gene_and_flank[, .(gene_methylation, flank_methylation)], mL4_snmcseq_gene_and_flank[, .(gene_methylation, flank_methylation)], mL5_snmcseq_gene_and_flank[, .(gene_methylation, flank_methylation)])
names(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank)[5:12] = c("Pv_gene_methylation", "Pv_flank_methylation", "Sst_gene_methylation", "Sst_flank_methylation", "L4_gene_methylation", "L4_flank_methylation", "L5_gene_methylation", "L5_flank_methylation")

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_Pv_gene_meth_diff := Sst_gene_methylation - Pv_gene_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_L4_gene_meth_diff := Sst_gene_methylation - L4_gene_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_L5_gene_meth_diff := Sst_gene_methylation - L5_gene_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Pv_L4_gene_meth_diff := Pv_gene_methylation - L4_gene_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Pv_L5_gene_meth_diff := Pv_gene_methylation - L5_gene_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, L5_L4_gene_meth_diff := L5_gene_methylation - L4_gene_methylation]

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_Pv_flank_meth_diff := Sst_flank_methylation - Pv_flank_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_L4_flank_meth_diff := Sst_flank_methylation - L4_flank_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Sst_L5_flank_meth_diff := Sst_flank_methylation - L5_flank_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Pv_L4_flank_meth_diff := Pv_flank_methylation - L4_flank_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, Pv_L5_flank_meth_diff := Pv_flank_methylation - L5_flank_methylation]
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[, L5_L4_flank_meth_diff := L5_flank_methylation - L4_flank_methylation]

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_Pv_gene_meth_diff = (Sst_Pv_gene_meth_diff - mean(Sst_Pv_gene_meth_diff, na.rm=TRUE))/sd(Sst_Pv_gene_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_L4_gene_meth_diff = (Sst_L4_gene_meth_diff - mean(Sst_L4_gene_meth_diff, na.rm=TRUE))/sd(Sst_L4_gene_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_L5_gene_meth_diff = (Sst_L5_gene_meth_diff - mean(Sst_L5_gene_meth_diff, na.rm=TRUE))/sd(Sst_L5_gene_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Pv_L4_gene_meth_diff = (Pv_L4_gene_meth_diff - mean(Pv_L4_gene_meth_diff, na.rm=TRUE))/sd(Pv_L4_gene_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Pv_L5_gene_meth_diff = (Pv_L5_gene_meth_diff - mean(Pv_L5_gene_meth_diff, na.rm=TRUE))/sd(Pv_L5_gene_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_L5_L4_gene_meth_diff = (L5_L4_gene_meth_diff - mean(L5_L4_gene_meth_diff, na.rm=TRUE))/sd(L5_L4_gene_meth_diff, na.rm=TRUE)))

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_Pv_flank_meth_diff = (Sst_Pv_flank_meth_diff - mean(Sst_Pv_flank_meth_diff, na.rm=TRUE))/sd(Sst_Pv_flank_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_L4_flank_meth_diff = (Sst_L4_flank_meth_diff - mean(Sst_L4_flank_meth_diff, na.rm=TRUE))/sd(Sst_L4_flank_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Sst_L5_flank_meth_diff = (Sst_L5_flank_meth_diff - mean(Sst_L5_flank_meth_diff, na.rm=TRUE))/sd(Sst_L5_flank_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Pv_L4_flank_meth_diff = (Pv_L4_flank_meth_diff - mean(Pv_L4_flank_meth_diff, na.rm=TRUE))/sd(Pv_L4_flank_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_Pv_L5_flank_meth_diff = (Pv_L5_flank_meth_diff - mean(Pv_L5_flank_meth_diff, na.rm=TRUE))/sd(Pv_L5_flank_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_snmcseq_gene_and_flank = data.table(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank %>% mutate(zScore_L5_L4_flank_meth_diff = (L5_L4_flank_meth_diff - mean(L5_L4_flank_meth_diff, na.rm=TRUE))/sd(L5_L4_flank_meth_diff, na.rm=TRUE)))

View(mPv_mSst_mL4_mL5_snmcseq_gene_and_flank)

View(mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Gene %in% pv_mr_L5_unchanged_genes) &
                                                                     (Intragenic_to_linked_gene==1) &
                                                                     (Gene %in% mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[zScore_Pv_L5_gene_meth_diff>=1, gene]) &
                                                                     (Gene %in% mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[abs(zScore_Pv_L5_flank_meth_diff)<1, gene]) &
                                                                     (cCRE_label %in% mousebrain_union_cCREs_mPv_mSst_mL4_mL5_mCA[(zScore_Pv_L5_cCRE_meth_diff>=1), cCRE_label]),])


mPv_mL5_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Pv_methylation/L5_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Pv_methylation/L5_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Pv_methylation/L5_methylation], element="Intragenic,\n cCRE-centered")))
mPv_mL5_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mSst_mPv_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Sst_methylation/Pv_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Sst_methylation/Pv_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_methylation/Pv_methylation], element="Intragenic,\n cCRE-centered")))
mSst_mPv_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mSst_mL5_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Sst_methylation/L5_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Sst_methylation/L5_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_methylation/L5_methylation], element="Intragenic,\n cCRE-centered")))
mSst_mL5_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mSst_mL4_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", Sst_methylation/L4_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", Sst_methylation/L4_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mSst_mL4_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mL5_mL4_10kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Extragenic", L5_methylation/L4_methylation], element="Extragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_mCA[genic_location=="Intragenic", L5_methylation/L4_methylation], element="Intragenic"),
  cbind(val=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, L5_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mL5_mL4_10kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]


ggplot(data=mPv_mL5_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mPv mCA/CA)/(mL5 mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,4))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPv_mL5_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPv_mL5_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSst_mPv_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSst mCA/CA)/(mPv mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,3))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mPv_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png",  width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mPv_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps",  width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)


ggplot(data=mSst_mL5_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSst mCA/CA)/(mL5 mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,4))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mL5_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png",  width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mL5_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps",  width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSst_mL4_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSst mCA/CA)/(mL4 mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,6))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png",  width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSst_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps",  width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mPv_mL4_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mPv mCA/CA)/(mL4 mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,5))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPv_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png",  width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPv_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps",  width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mL5_mL4_10kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mL5 mCA/CA)/(mL4 mCA/CA)") + ylab("10kb regions")+
  coord_cartesian(xlim=c(0,3))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png",  width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_10kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps",  width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)


pv_sst_mr_L4_L5_unchanged_genes = intersect(intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], sst_mr_genes_q0.1_nondedup_mm9[,V4]), intersect(L4_unchanged_genes_p0.5_nondedup_mm9[,V4], L5_unchanged_genes_p0.5_nondedup_mm9[,V4]))
sst_L5_mr_pv_L4_unchanged_genes = intersect(intersect(sst_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]), intersect(pv_unchanged_genes_p0.5_nondedup_mm9[,V4], L4_unchanged_genes_p0.5_nondedup_mm9[,V4]))
pv_L4_mr_sst_L5_unchanged_genes = intersect(intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], L4_mr_genes_q0.1_nondedup_mm9[,V4]), intersect(sst_unchanged_genes_p0.5_nondedup_mm9[,V4], L5_unchanged_genes_p0.5_nondedup_mm9[,V4]))
L4_L5_mr_pv_sst_unchanged_genes = intersect(intersect(L4_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]), intersect(pv_unchanged_genes_p0.5_nondedup_mm9[,V4], sst_unchanged_genes_p0.5_nondedup_mm9[,V4]))


pv_sst_mr_L4_L5_subThresholdUnchanged_genes = setdiff(intersect(pv_mr_genes_q0.1_nondedup_mm9[,V4], sst_mr_genes_q0.1_nondedup_mm9[,V4]), c(L4_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]))

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[(gene %in% pv_sst_mr_L4_L5_unchanged_genes) &
                                        (abs(zScore_Sst_Pv_flank_meth_diff)<1) &
                                        (abs(zScore_Sst_L4_flank_meth_diff)<1) &
                                        (abs(zScore_Sst_L5_flank_meth_diff)<1) &
                                        (abs(zScore_Pv_L4_flank_meth_diff)<1) &
                                        (abs(zScore_Pv_L5_flank_meth_diff)<1) &
                                        (abs(zScore_L5_L4_flank_meth_diff)<1),
                                        gene]

mPv_mSst_mL4_mL5_snmcseq_gene_and_flank[(gene %in% sst_L5_mr_pv_L4_unchanged_genes) &
                                          (abs(zScore_Sst_Pv_flank_meth_diff)<1) &
                                          (abs(zScore_Sst_L4_flank_meth_diff)<1) &
                                          (abs(zScore_Sst_L5_flank_meth_diff)<1) &
                                          (abs(zScore_Pv_L4_flank_meth_diff)<1) &
                                          (abs(zScore_Pv_L5_flank_meth_diff)<1) &
                                          (abs(zScore_L5_L4_flank_meth_diff)<1),
                                        gene]

#calculating methylation differences in the 10kb window regions
mPv_mSst_mL4_mL5_10kb_mCA[, Sst_Pv_meth_diff :=  Sst_methylation - Pv_methylation]
mPv_mSst_mL4_mL5_10kb_mCA[, Sst_L4_meth_diff :=  Sst_methylation - L4_methylation]
mPv_mSst_mL4_mL5_10kb_mCA[, Sst_L5_meth_diff :=  Sst_methylation - L5_methylation]
mPv_mSst_mL4_mL5_10kb_mCA[, Pv_L4_meth_diff :=  Pv_methylation - L4_methylation]
mPv_mSst_mL4_mL5_10kb_mCA[, Pv_L5_meth_diff :=  Pv_methylation - L5_methylation]
mPv_mSst_mL4_mL5_10kb_mCA[, L5_L4_meth_diff :=  L5_methylation - L4_methylation]

mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_Pv_meth_diff :=  Sst_methylation - Pv_methylation]
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_L4_meth_diff :=  Sst_methylation - L4_methylation]
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Sst_L5_meth_diff :=  Sst_methylation - L5_methylation]
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Pv_L4_meth_diff :=  Pv_methylation - L4_methylation]
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, Pv_L5_meth_diff :=  Pv_methylation - L5_methylation]
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked[, L5_L4_meth_diff :=  L5_methylation - L4_methylation]


#calculating z-scores of methylation differences
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_Sst_Pv_meth_diff = (Sst_Pv_meth_diff - mean(Sst_Pv_meth_diff, na.rm=TRUE))/sd(Sst_Pv_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_Sst_L4_meth_diff = (Sst_L4_meth_diff - mean(Sst_L4_meth_diff, na.rm=TRUE))/sd(Sst_L4_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_Sst_L5_meth_diff = (Sst_L5_meth_diff - mean(Sst_L5_meth_diff, na.rm=TRUE))/sd(Sst_L5_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_Pv_L4_meth_diff = (Pv_L4_meth_diff - mean(Pv_L4_meth_diff, na.rm=TRUE))/sd(Pv_L4_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_Pv_L5_meth_diff = (Pv_L5_meth_diff - mean(Pv_L5_meth_diff, na.rm=TRUE))/sd(Pv_L5_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_mCA = data.table(mPv_mSst_mL4_mL5_10kb_mCA %>% mutate(zScore_L5_L4_meth_diff = (L5_L4_meth_diff - mean(L5_L4_meth_diff, na.rm=TRUE))/sd(L5_L4_meth_diff, na.rm=TRUE)))

mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_Sst_Pv_meth_diff = (Sst_Pv_meth_diff - mean(Sst_Pv_meth_diff, na.rm=TRUE))/sd(Sst_Pv_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_Sst_L4_meth_diff = (Sst_L4_meth_diff - mean(Sst_L4_meth_diff, na.rm=TRUE))/sd(Sst_L4_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_Sst_L5_meth_diff = (Sst_L5_meth_diff - mean(Sst_L5_meth_diff, na.rm=TRUE))/sd(Sst_L5_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_Pv_L4_meth_diff = (Pv_L4_meth_diff - mean(Pv_L4_meth_diff, na.rm=TRUE))/sd(Pv_L4_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_Pv_L5_meth_diff = (Pv_L5_meth_diff - mean(Pv_L5_meth_diff, na.rm=TRUE))/sd(Pv_L5_meth_diff, na.rm=TRUE)))
mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked = data.table(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked %>% mutate(zScore_L5_L4_meth_diff = (L5_L4_meth_diff - mean(L5_L4_meth_diff, na.rm=TRUE))/sd(L5_L4_meth_diff, na.rm=TRUE)))

#plot of distribution of all extragenic and intragenic regions
ggplot(mPv_mSst_mL4_mL5_10kb_mCA, aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("All 10kb regions")+
  xlab("Genic location of 10kb regions") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/all_10kbWindows_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/all_10kbWindows_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_Sst_Pv_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mSst vs mPv differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mSst-mPv mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mPv_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mPv_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')


ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_Sst_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mSst vs mL4 differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mSst-mL4 mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')


ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_Sst_L5_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mSst vs mL5 differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mSst-mL5 mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')


ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_Pv_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mPv vs mL4 differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mPv-mL4 mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')


ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_Pv_L5_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mPv vs mL5 differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mPv-mL5 mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA[(abs(zScore_L5_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  ggtitle("mL5 vs mL4 differential mCA/CA regions")+
  xlab("Genic location of 10kb regions,\n |z-score of mL5-mL4 mCA/CA| >= 2") + ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mL5_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.png", width=3.7, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mL5_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_barPlot.eps", width=3.7, height=5, dpi = 300, units = "in", device='eps')

setkeyv(mPv_mSst_mL4_mL5_10kb_mCA, cols=c("chrom", "start", "end"))
setkeyv(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked, cols=c("chrom", "start", "end"))
#determining the 10kb genomic regions that overlap the 10kb regions previously defined as centered at intragenic cognate-linked cCREs
mPv_mSst_mL4_mL5_10kb_mCA_coords_within5kb_intragenicLinkedcCRE_centers = unique(foverlaps(y=mPv_mSst_mL4_mL5_10kb_mCA, x=mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked, by.y=key(mPv_mSst_mL4_mL5_10kb_mCA), by.x=key(mPv_mSst_mL4_mL5_10kb_cCRE_mCA_intragenicLinked), type="any", nomatch=NULL)[, coord])

#z-score column names
zScore_colnames = c("zScore_Sst_Pv_meth_diff", "zScore_Sst_L4_meth_diff", "zScore_Sst_L5_meth_diff", "zScore_Pv_L4_meth_diff", "zScore_Pv_L5_meth_diff", "zScore_L5_L4_meth_diff")
mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity = data.table(rbind(cbind(mPv_mSst_mL4_mL5_10kb_mCA[!(coord %in% mPv_mSst_mL4_mL5_10kb_mCA_coords_within5kb_intragenicLinkedcCRE_centers) & (genic_location=="Extragenic"), ..zScore_colnames], genic_location = "Extragenic, >5kb from cognate-linked cCREs"),
                                                           cbind(mPv_mSst_mL4_mL5_10kb_mCA[!(coord %in% mPv_mSst_mL4_mL5_10kb_mCA_coords_within5kb_intragenicLinkedcCRE_centers) & (genic_location=="Intragenic"), ..zScore_colnames], genic_location = "Intragenic, >5kb from cognate-linked cCREs"),
                                                           cbind(mPv_mSst_mL4_mL5_10kb_mCA[(coord %in% mPv_mSst_mL4_mL5_10kb_mCA_coords_within5kb_intragenicLinkedcCRE_centers) & (genic_location=="Extragenic"), ..zScore_colnames], genic_location = "Extragenic, <=5kb from cognate-linked cCREs"),
                                                           cbind(mPv_mSst_mL4_mL5_10kb_mCA[(coord %in% mPv_mSst_mL4_mL5_10kb_mCA_coords_within5kb_intragenicLinkedcCRE_centers) & (genic_location=="Intragenic"), ..zScore_colnames], genic_location = "Intragenic, <=5kb from cognate-linked cCREs")))

mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity = mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity %>% mutate(genic_location = factor(genic_location, levels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs")))

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_Sst_Pv_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mSst vs mPv differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mPv_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mPv_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_Sst_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mSst vs mL4 differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_Sst_L5_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mSst vs mL5 differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mSst_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_Pv_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mPv vs mL4 differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_Pv_L5_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mPv vs mL5 differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mPv_mL5_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')

ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_L5_L4_meth_diff) >= 2),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mL5 vs mL4 differential mCA/CA regions,\n |z-score| >=2")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mL5_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.png", width=4.5, height=5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/genic_distribution_plots/mL5_mL4_10kbWindows_diffmCAperCA_zScore2_extragenic_intragenic_dist_cognateLinkedcCRE_proximity_barPlot.eps", width=4.5, height=5, dpi = 300, units = "in", device='eps')


ggplot(mPv_mSst_mL4_mL5_10kb_mCA_cCREproximity[(abs(zScore_L5_L4_meth_diff) >= 3),], aes(x = genic_location, y = ..prop.., fill=factor(..x..), group = 1)) +
  geom_bar()+
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.15) +
  scale_fill_manual(guide=guide_legend(nrow=4, byrow=TRUE), name="Genic location:", values = c("darkred", "darksalmon", "dodgerblue4", "dodgerblue1"),
                    labels=c("Extragenic, >5kb from cognate-linked cCREs", "Extragenic, <=5kb from cognate-linked cCREs", "Intragenic, >5kb from cognate-linked cCREs", "Intragenic, <=5kb from cognate-linked cCREs"))+
  ggtitle("mL5 vs mL4 differential mCA/CA regions,\n |z-score| >= 3")+
  ylab("Proportion") +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=12), legend.position = "bottom", legend.margin=margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_text(size=10), axis.text.x = element_blank())

###using 1kb windows for examining methylation variation by genic location of regions
PV_1kb_mCA = fread("HG_lab/Mati/GabelLab/INTACT_methylation/PV/PV_WT_KO_deep_INTACT_mCA_mm9_1kb_windows.bed")
SST_1kb_mCA = fread("HG_lab/Mati/GabelLab/INTACT_methylation/SST/SST_WT_KO_deep_INTACT_mCA_mm9_1kb_windows.bed")
L4_1kb_mCA = fread("HG_lab/Mati/GabelLab/INTACT_methylation/L4/L4_WT_KO_deep_INTACT_mCA_mm9_1kb_windows.bed")
L5_1kb_mCA = fread("HG_lab/Mati/GabelLab/INTACT_methylation/L5/L5_WT_KO_deep_INTACT_mCA_mm9_1kb_windows.bed")


names(PV_1kb_mCA) = c("chrom", "start", "end", "PV_meth", "PV_cov")
PV_1kb_mCA[, PV_unmeth := PV_cov - PV_meth]
PV_1kb_mCA[, PV_methylation := PV_meth/PV_cov]
PV_1kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
PV_1kb_mCA = PV_1kb_mCA[chrom %in% chrom_list]

names(SST_1kb_mCA) = c("chrom", "start", "end", "SST_meth", "SST_cov")
SST_1kb_mCA[, SST_unmeth := SST_cov - SST_meth]
SST_1kb_mCA[, SST_methylation := SST_meth/SST_cov]
SST_1kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
SST_1kb_mCA = SST_1kb_mCA[chrom %in% chrom_list]

names(L4_1kb_mCA) = c("chrom", "start", "end", "L4_meth", "L4_cov")
L4_1kb_mCA[, L4_unmeth := L4_cov - L4_meth]
L4_1kb_mCA[, L4_methylation := L4_meth/L4_cov]
L4_1kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
L4_1kb_mCA = L4_1kb_mCA[chrom %in% chrom_list]

names(L5_1kb_mCA) = c("chrom", "start", "end", "L5_meth", "L5_cov")
L5_1kb_mCA[, L5_unmeth := L5_cov - L5_meth]
L5_1kb_mCA[, L5_methylation := L5_meth/L5_cov]
L5_1kb_mCA[, coord := paste0(chrom,":",start,"-",end)]
L5_1kb_mCA = L5_1kb_mCA[chrom %in% chrom_list]

#average nonconversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")

#subtracting nonconversion rate
PV_1kb_mCA$PV_methylation_corrected <- PV_1kb_mCA$PV_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
PV_1kb_mCA[PV_methylation_corrected < 0, PV_methylation_corrected := 0]

SST_1kb_mCA$SST_methylation_corrected <- SST_1kb_mCA$SST_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
SST_1kb_mCA[SST_methylation_corrected < 0, SST_methylation_corrected := 0]

L4_1kb_mCA$L4_methylation_corrected <- L4_1kb_mCA$L4_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L4_1kb_mCA[L4_methylation_corrected < 0, L4_methylation_corrected := 0]

L5_1kb_mCA$L5_methylation_corrected <- L5_1kb_mCA$L5_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
L5_1kb_mCA[L5_methylation_corrected < 0, L5_methylation_corrected := 0]



mm9_1kb_windows_genicRegions = fread("HG_lab/Mati/GabelLab/databases/mm9_1kb_windows_genicRegions.bed")

PV_SST_1kb_mCA = unique(data.table(inner_join(x=PV_1kb_mCA, y=SST_1kb_mCA, by=c("chrom","start","end","coord"))))
PV_SST_L4_1kb_mCA = unique(data.table(inner_join(x=PV_SST_1kb_mCA, y=L4_1kb_mCA, by=c("chrom","start","end","coord"))))
PV_SST_L4_L5_1kb_mCA = unique(data.table(inner_join(x=PV_SST_L4_1kb_mCA, y=L5_1kb_mCA, by=c("chrom","start","end","coord"))))


PV_SST_L4_L5_1kb_mCA[(coord %in% mm9_1kb_windows_genicRegions[, V4]), genic_location := "Intragenic"]
PV_SST_L4_L5_1kb_mCA[!(coord %in% mm9_1kb_windows_genicRegions[, V4]), genic_location := "Extragenic"]

#1kb windows centered on nonpromoter cCREs
PV_1kb_cCRE_mCA = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_1kbWindows/mousebrain_union_nonPromoter_cCREs_1kbWindows_PV_WT_KO_deep_INTACT_CA_mm9.bed")
SST_1kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_1kbWindows/mousebrain_union_nonPromoter_cCREs_1kbWindows_SST_WT_KO_deep_INTACT_CA_mm9.bed")
L4_1kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_1kbWindows/mousebrain_union_nonPromoter_cCREs_1kbWindows_L4_WT_KO_deep_INTACT_CA_mm9.bed")
L5_1kb_cCRE_mCA  = fread("HG_lab/Mati/GabelLab/ATACseq/Li_Ren_2021/Union/cCRE_centered_1kbWindows/mousebrain_union_nonPromoter_cCREs_1kbWindows_L5_WT_KO_deep_INTACT_CA_mm9.bed")

names(PV_1kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "PV_meth", "PV_cov")
names(SST_1kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "SST_meth", "SST_cov")
names(L4_1kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "L4_meth", "L4_cov")
names(L5_1kb_cCRE_mCA) = c("chrom", "start", "end", "cCRE_label", "L5_meth", "L5_cov")

PV_1kb_cCRE_mCA[, PV_methylation := PV_meth/PV_cov]
SST_1kb_cCRE_mCA[, SST_methylation := SST_meth/SST_cov]
L4_1kb_cCRE_mCA[, L4_methylation := L4_meth/L4_cov]
L5_1kb_cCRE_mCA[, L5_methylation := L5_meth/L5_cov]

#subtracting nonconversion rate
PV_1kb_cCRE_mCA$PV_methylation_corrected <- PV_1kb_cCRE_mCA$PV_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
PV_1kb_cCRE_mCA[PV_methylation_corrected < 0, PV_methylation_corrected := 0]

SST_1kb_cCRE_mCA$SST_methylation_corrected <- SST_1kb_cCRE_mCA$SST_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
SST_1kb_cCRE_mCA[SST_methylation_corrected < 0, SST_methylation_corrected := 0]

L4_1kb_cCRE_mCA$L4_methylation_corrected <- L4_1kb_cCRE_mCA$L4_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L4_1kb_cCRE_mCA[L4_methylation_corrected < 0, L4_methylation_corrected := 0]

L5_1kb_cCRE_mCA$L5_methylation_corrected <- L5_1kb_cCRE_mCA$L5_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
L5_1kb_cCRE_mCA[L5_methylation_corrected < 0, L5_methylation_corrected := 0]



PV_SST_1kb_cCRE_mCA = unique(data.table(inner_join(x=PV_1kb_cCRE_mCA, y=SST_1kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))
PV_SST_L4_1kb_cCRE_mCA = unique(data.table(inner_join(x=PV_SST_1kb_cCRE_mCA, y=L4_1kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))
PV_SST_L4_L5_1kb_cCRE_mCA = unique(data.table(inner_join(x=PV_SST_L4_1kb_cCRE_mCA, y=L5_1kb_cCRE_mCA, by=c("chrom","start","end","cCRE_label"))))

PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked = PV_SST_L4_L5_1kb_cCRE_mCA[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Intragenic==1) & (Intragenic_to_linked_gene==1), cCRE_label],]
PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked = PV_SST_L4_L5_1kb_cCRE_mCA[cCRE_label %in% mousebrain_union_nonPromoter_cCREs_linked_genes_genicBooleans[(Intragenic==0) & (Intragenic_to_linked_gene==0), cCRE_label],]


mSST_mPV_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/PV_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/PV_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/PV_methylation], element="Intragenic,\n cCRE-centered")))
mSST_mPV_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mSST_mL4_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mSST_mL4_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mSST_mL5_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/L5_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/L5_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/L5_methylation], element="Intragenic,\n cCRE-centered")))
mSST_mL5_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mPV_mL4_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", PV_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", PV_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, PV_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mPV_mL4_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

mPV_mL5_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", PV_methylation/L5_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", PV_methylation/L5_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, PV_methylation/L5_methylation], element="Intragenic,\n cCRE-centered")))
mPV_mL5_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]


mL5_mL4_1kb_windows_mCA_methFoldDiff = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", L5_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", L5_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, L5_methylation/L4_methylation], element="Intragenic,\n cCRE-centered")))
mL5_mL4_1kb_windows_mCA_methFoldDiff[, val := as.numeric(val)]

ggplot(data=mSST_mPV_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mSST_mPV_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mPV_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mPV_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mSST mCA/CA)/(mPV mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,6))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mPV_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mPV_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSST_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mSST_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mSST mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSST_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mSST_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                               round(sd(mSST_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mSST mCA/CA)/(mL5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,8))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mPV_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mPV_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                               round(sd(mPV_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                               round(sd(mPV_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mPV mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mPV_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mPV_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                               round(sd(mPV_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                               round(sd(mPV_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mPV mCA/CA)/(mL5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,7))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mPV_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mPV mCA/CA)/(mL5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,7))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


ggplot(data=mL5_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow"),
                    name="Standard\n deviation:",
                    labels = c(round(sd(mL5_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),2),
                             round(sd(mL5_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),2),
                             round(sd(mL5_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),2)),
                    guide=guide_legend(nrow=3, byrow=TRUE)) +
  xlab("(mL5 mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,5))+
  theme_bw()+
  theme(legend.position = "right", legend.margin=margin(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 4, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 4, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mL5_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mL5 mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,5))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mL5_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mPV_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mPV mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mPV_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSST_mPV_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSST mCA/CA)/(mPV mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,6))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mPV_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mPV_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSST_mL4_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSST mCA/CA)/(mL4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL4_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=mSST_mL5_1kb_windows_mCA_methFoldDiff[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow")) +
  xlab("(mSST mCA/CA)/(mL5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,8))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/mSST_mL5_1kbRegions_extragenic_intragenic_intragenicLinkedcCRE_mCAperCA_foldDiff_density_ridgelinePlot2.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

#function for standard deviation of methylation difference
genicLoc_sd_methFoldDiff = function(data_table){
c(sd(data_table[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),
  sd(data_table[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),
  sd(data_table[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE))
}

cellType_1kbWindows_mCA_methFoldDiff_matrix = rbind(genicLoc_sd_methFoldDiff(mSST_mPV_1kb_windows_mCA_methFoldDiff),
                                genicLoc_sd_methFoldDiff(mSST_mL4_1kb_windows_mCA_methFoldDiff),
                                genicLoc_sd_methFoldDiff(mSST_mL5_1kb_windows_mCA_methFoldDiff),
                                genicLoc_sd_methFoldDiff(mPV_mL4_1kb_windows_mCA_methFoldDiff),
                                genicLoc_sd_methFoldDiff(mPV_mL5_1kb_windows_mCA_methFoldDiff),
                                genicLoc_sd_methFoldDiff(mL5_mL4_1kb_windows_mCA_methFoldDiff))
row.names(cellType_1kbWindows_mCA_methFoldDiff_matrix) = c("SST/PV", "SST/L4", "SST/L5", "PV/L4", "PV/L5", "L5/L4")
colnames(cellType_1kbWindows_mCA_methFoldDiff_matrix) = c("Extragenic", "Intragenic", "Intragenic,\n cCRE-centered")

pheatmap(cellType_1kbWindows_mCA_methFoldDiff_matrix,
         color = colorRampPalette(c("white","goldenrod2"))(99),
         display_numbers = TRUE,
         #labels_col=c("Region", "Gene body", "cCRE"),
         fontsize=20,
         #fontsize_row = 15,
         #fontsize_col = 15,
         fontsize_number = 20,
         cluster_cols=FALSE,
         cluster_rows=FALSE
)
heatmap_func <- function(color_matrix, number_matrix, column_names, color, color_max, title, x_size=18, y_size=18, x_angle=0, y_angle=0, tile_text_size=6){
  test_color <- as.data.frame(color_matrix)
  test_number <- as.data.frame(number_matrix)
  colnames(test_color) = column_names
  colnames(test_number) = column_names
  test_color$cellType = row.names(test_color)
  test_number$cellType = row.names(test_number)
  #Melting data so we can plot it with GGplot
  test_color.m <- suppressWarnings(melt(test_color,id.vars = c("cellType")))
  test_number.m <- suppressWarnings(melt(test_number,id.vars = c("cellType")))
  #Resetting factors
  test_color.m = test_color.m %>% mutate(cellType = factor(cellType, levels=c("L5/L4", "PV/L5", "PV/L4", "SST/L5", "SST/L4", "SST/PV")))
  test_number.m = test_number.m %>% mutate(cellType = factor(cellType, levels=c("L5/L4", "PV/L5", "PV/L4", "SST/L5", "SST/L4", "SST/PV")))
  test_color_and_number = data.table(cbind(test_color.m, number = test_number.m$value))
  test_color_and_number$colorMax = test_color_and_number$value
  test_color_and_number[value > color_max, colorMax := color_max]
  
  #Creating the plot itself
  ggplot(test_color_and_number,aes(variable,cellType)) + geom_tile(aes(fill=colorMax),color = "white") +
    ggtitle(title)+
    #Creating legend
    guides(fill=guide_colorbar(title="Standard deviation", title.position="top")) +
    #Creating color range
    scale_fill_gradientn(limits = c(0,color_max), colors=c("white",color),guide="colorbar") +
    geom_text(aes(label = colorMax), color = "black", size = tile_text_size)+
    theme_bw()+
    #Rotating labels
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
          axis.text.x=element_text(size=x_size, angle=x_angle), axis.text.y=element_text(size=y_size, angle=y_angle), axis.ticks = element_blank(), axis.title=element_blank())
}
heatmap_func(round(cellType_1kbWindows_mCA_methFoldDiff_matrix,2), round(cellType_1kbWindows_mCA_methFoldDiff_matrix,2), c("Extragenic", "Intragenic", "Intragenic,\n cCRE-centered"), "goldenrod2", 3, "mCA/CA fold difference\n standard deviations",  x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/cellType_1kbWindows_mCA_methFoldDiff_stdevs_heatmap.png", width=3.5, height=5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/cellType_1kbWindows_mCA_methFoldDiff_stdevs_heatmap.eps", width=3.5, height=5.5, dpi = 300, units = "in", device='eps')


heatmap_func_noText <- function(color_matrix, number_matrix, column_names, color, color_max, title, x_size=18, y_size=18, x_angle=0, y_angle=0, tile_text_size=6){
  test_color <- as.data.frame(color_matrix)
  test_number <- as.data.frame(number_matrix)
  colnames(test_color) = column_names
  colnames(test_number) = column_names
  test_color$cellType = row.names(test_color)
  test_number$cellType = row.names(test_number)
  #Melting data so we can plot it with GGplot
  test_color.m <- suppressWarnings(melt(test_color,id.vars = c("cellType")))
  test_number.m <- suppressWarnings(melt(test_number,id.vars = c("cellType")))
  #Resetting factors
  test_color.m = test_color.m %>% mutate(cellType = factor(cellType, levels=c("L5/L4", "PV/L5", "PV/L4", "SST/L5", "SST/L4", "SST/PV")))
  test_number.m = test_number.m %>% mutate(cellType = factor(cellType, levels=c("L5/L4", "PV/L5", "PV/L4", "SST/L5", "SST/L4", "SST/PV")))
  test_color_and_number = data.table(cbind(test_color.m, number = test_number.m$value))
  test_color_and_number$colorMax = test_color_and_number$value
  test_color_and_number[value > color_max, colorMax := color_max]
  
  #Creating the plot itself
  ggplot(test_color_and_number,aes(variable,cellType)) + geom_tile(aes(fill=colorMax),color = "white") +
    ggtitle(title)+
    #Creating legend
    guides(fill=guide_colorbar(title="Standard deviation", title.position="top")) +
    #Creating color range
    scale_fill_viridis()+
    #scale_fill_gradientn(limits = c(0,color_max), colors=c("white",color),guide="colorbar") +
    #geom_text(aes(label = colorMax), color = "black", size = tile_text_size)+
    theme_bw()+
    #Rotating labels
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
          axis.text.x=element_text(size=x_size, angle=x_angle), axis.text.y=element_text(size=y_size, angle=y_angle), axis.ticks = element_blank(), axis.title=element_blank())
}


##ridgeline plots and heatmaps including extragenic cCRE-centered regions
SST_PV_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/PV_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/PV_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/PV_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, SST_methylation/PV_methylation], element="Extragenic,\n cCRE-centered")))
SST_PV_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
SST_PV_1kb_windows_mCA_methFoldDiff_extra = SST_PV_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))

SST_L4_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/L4_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, SST_methylation/L4_methylation], element="Extragenic,\n cCRE-centered")))
SST_L4_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
SST_L4_1kb_windows_mCA_methFoldDiff_extra = SST_L4_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))


SST_L5_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", SST_methylation/L5_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", SST_methylation/L5_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, SST_methylation/L5_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, SST_methylation/L5_methylation], element="Extragenic,\n cCRE-centered")))
SST_L5_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
SST_L5_1kb_windows_mCA_methFoldDiff_extra = SST_L5_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))


PV_L4_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", PV_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", PV_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, PV_methylation/L4_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, PV_methylation/L4_methylation], element="Extragenic,\n cCRE-centered")))
PV_L4_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
PV_L4_1kb_windows_mCA_methFoldDiff_extra = PV_L4_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))


PV_L5_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", PV_methylation/L5_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", PV_methylation/L5_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, PV_methylation/L5_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, PV_methylation/L5_methylation], element="Extragenic,\n cCRE-centered")))
PV_L5_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
PV_L5_1kb_windows_mCA_methFoldDiff_extra = PV_L5_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))


L5_L4_1kb_windows_mCA_methFoldDiff_extra = data.table(rbind(
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Extragenic", L5_methylation/L4_methylation], element="Extragenic"),
  cbind(val=PV_SST_L4_L5_1kb_mCA[genic_location=="Intragenic", L5_methylation/L4_methylation], element="Intragenic"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_intragenicLinked[, L5_methylation/L4_methylation], element="Intragenic,\n cCRE-centered"),
  cbind(val=PV_SST_L4_L5_1kb_cCRE_mCA_extragenicLinked[, L5_methylation/L4_methylation], element="Extragenic,\n cCRE-centered")))
L5_L4_1kb_windows_mCA_methFoldDiff_extra[, val := as.numeric(val)]
L5_L4_1kb_windows_mCA_methFoldDiff_extra = L5_L4_1kb_windows_mCA_methFoldDiff_extra %>% mutate(element = factor(element, levels=c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")))


ggplot(data=PV_L5_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(PV mCA/CA)/(L5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,7))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_PV_L5_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_PV_L5_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)


ggplot(data=L5_L4_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(L5 mCA/CA)/(L4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,5))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_L5_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_L5_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=PV_L4_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(PV mCA/CA)/(L4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_PV_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_PV_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=SST_PV_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(SST mCA/CA)/(PV mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,6))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_PV_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_PV_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=SST_L4_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(SST mCA/CA)/(L4 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,10))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_L4_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

ggplot(data=SST_L5_1kb_windows_mCA_methFoldDiff_extra[is.finite(val)], aes(x = val, y = element, fill=element)) + 
  geom_density_ridges(rel_min_height=0.01, scale = 7, alpha=0.4)+
  scale_fill_manual(values = c("Extragenic"="darkred", "Intragenic"="lightblue", "Intragenic,\n cCRE-centered"="yellow", "Extragenic,\n cCRE-centered"="darkgreen")) +
  xlab("(SST mCA/CA)/(L5 mCA/CA)") + ylab("1kb regions")+
  coord_cartesian(xlim=c(0,8))+
  theme_bw()+
  theme(legend.position = "None", legend.margin=margin(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(color="black"), axis.text.y=element_text(size=10), axis.text.x = element_text(size=10))
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_L5_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.png", width = 5, height = 5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_SST_L5_1kbRegions_extragenic_intragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCAperCA_foldDiff_density_ridgelinePlot.eps", width = 5, height = 5, dpi = 300, units = "in", device=cairo_ps)

#function for standard deviation of methylation difference
genicLoc_sd_methFoldDiff_extra = function(data_table){
  c(sd(data_table[is.finite(val) & (element=="Intragenic"), val], na.rm=TRUE),
    sd(data_table[is.finite(val) & (element=="Extragenic"), val], na.rm=TRUE),
    sd(data_table[is.finite(val) & (element=="Intragenic,\n cCRE-centered"), val], na.rm=TRUE),
    sd(data_table[is.finite(val) & (element=="Extragenic,\n cCRE-centered"), val], na.rm=TRUE))
}

cellType_1kbWindows_mCA_methFoldDiff_extra_matrix = rbind(genicLoc_sd_methFoldDiff_extra(SST_PV_1kb_windows_mCA_methFoldDiff_extra),
                                                    genicLoc_sd_methFoldDiff_extra(SST_L4_1kb_windows_mCA_methFoldDiff_extra),
                                                    genicLoc_sd_methFoldDiff_extra(SST_L5_1kb_windows_mCA_methFoldDiff_extra),
                                                    genicLoc_sd_methFoldDiff_extra(PV_L4_1kb_windows_mCA_methFoldDiff_extra),
                                                    genicLoc_sd_methFoldDiff_extra(PV_L5_1kb_windows_mCA_methFoldDiff_extra),
                                                    genicLoc_sd_methFoldDiff_extra(L5_L4_1kb_windows_mCA_methFoldDiff_extra))
row.names(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix) = c("SST/PV", "SST/L4", "SST/L5", "PV/L4", "PV/L5", "L5/L4")
colnames(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix) = c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered")

heatmap_func_noText(round(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix,2), round(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix,2), c("Intragenic", "Extragenic", "Intragenic,\n cCRE-centered", "Extragenic,\n cCRE-centered"), "green", 3, "mCA/CA fold difference\n standard deviations",  x_size=16, y_size=16, x_angle=90, y_angle=0, tile_text_size = 8)
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCA_methFoldDiff_stdevs_heatmap_noText.png", width=3.5, height=5.5, dpi = 300, units = "in", device='png')
ggsave("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/backgroundSub_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCA_methFoldDiff_stdevs_heatmap_noText.eps", width=3.5, height=5.5, dpi = 300, units = "in", device='eps')


cellType_1kbWindows_mCA_methFoldDiff_extra_matrix_dt <- data.table(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix, keep.rownames="Comparison")
names(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix_dt) <- c("Comparison", "Intragenic", "Extragenic", "Intragenic_cCRE-centered", "Extragenic_cCRE-centered")
write.csv(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix_dt, "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/table_of_backgroundSub_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCA_methFoldDiff_stdevs.csv", quote=F, row.names=F)
write.table(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix_dt, "HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/table_of_backgroundSub_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_WT_KO_INTACT_mCA_methFoldDiff_stdevs.txt", quote=F, row.names=F, sep="\t")


#
# Get the number of columns
# Get the number of columns
num_columns <- ncol(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix)
columns <- colnames(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix)

# Initialize a matrix to store p-values
p_value_matrix <- matrix(NA, nrow = num_columns, ncol = num_columns)
rownames(p_value_matrix) <- columns
colnames(p_value_matrix) <- columns

# Perform paired Wilcoxon rank sum tests between each pair of columns
for (i in 1:(num_columns - 1)) {
  for (j in (i + 1):num_columns) {
    col1 <- columns[i]
    col2 <- columns[j]
    test_result <- wilcox.test(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col1], 
                               cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col2], 
                               paired = TRUE)
    p_value_matrix[i, j] <- test_result$p.value
    p_value_matrix[j, i] <- test_result$p.value  # Since the test is symmetric
  }
}

# Print the p-value matrix
print(p_value_matrix)

rownames(p_value_matrix) <- c("Intragenic", "Extragenic", "Intragenic_cCRE-centered", "Extragenic_cCRE_centered")
colnames(p_value_matrix) <- c("Intragenic", "Extragenic", "Intragenic_cCRE-centered", "Extragenic_cCRE_centered")


log10_p_value_matrix <- -log10(p_value_matrix)

wilcox_col_palette <- colorRampPalette(c("pink", "red", "darkred"))(n = 299)

png("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/wilcox_negLog10_pval_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_INTACT_mCA_methFoldDiff_stdevs_heatmap.png", width=1500, height=1500, res=300)
heatmap.2(log10_p_value_matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=wilcox_col_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, key.title="-Log10 p-value",
          cexRow=0.6,            # adjust row label text size (smaller value for smaller text)
          cexCol=0.6             # adjust column label text size (smaller value for smaller text))
)
dev.off()

###t.test
# Initialize a matrix to store p-values
t_test_p_value_matrix <- matrix(NA, nrow = num_columns, ncol = num_columns)
rownames(t_test_p_value_matrix) <- columns
colnames(t_test_p_value_matrix) <- columns

# Perform paired Wilcoxon rank sum tests between each pair of columns
for (i in 1:(num_columns - 1)) {
  for (j in (i + 1):num_columns) {
    col1 <- columns[i]
    col2 <- columns[j]
    t_test_result <- t.test(cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col1], 
                               cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col2], 
                               paired = TRUE)
    t_test_p_value_matrix[i, j] <- t_test_result$p.value
    t_test_p_value_matrix[j, i] <- t_test_result$p.value  # Since the test is symmetric
  }
}

##shapiro test
# Initialize a matrix to store p-values
shapiro_p_value_matrix <- matrix(NA, nrow = num_columns, ncol = num_columns)
rownames(shapiro_p_value_matrix) <- columns
colnames(shapiro_p_value_matrix) <- columns

# Perform paired Wilcoxon rank sum tests between each pair of columns
for (i in 1:(num_columns - 1)) {
  for (j in (i + 1):num_columns) {
    col1 <- columns[i]
    col2 <- columns[j]
    #differences in mCA fold difference standard deviations between groups-
    differences <- cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col1] - cellType_1kbWindows_mCA_methFoldDiff_extra_matrix[, col2]
    shapiro_test_result <- shapiro.test(differences)
    shapiro_p_value_matrix[i, j] <- shapiro_test_result$p.value
    shapiro_p_value_matrix[j, i] <- shapiro_test_result$p.value  # Since the test is symmetric
  }
}

shapiro_p_value_matrix[1:4, 1:4]
min(shapiro_p_value_matrix[1:4, 1:4], na.rm=TRUE) #0.4939623

#Conclusion, these differences are normally distributed; paired t-test should be appropriate
p_value_matrix_sub <- p_value_matrix[1:4, 1:4]
p_value_matrix_sub_dt <- data.table(p_value_matrix_sub, keep.row.names="Group")
t_test_p_value_matrix_sub<- t_test_p_value_matrix[1:4, 1:4]
t_test_p_value_matrix_sub_dt <- data.table(t_test_p_value_matrix_sub, keep.row.names="Group")

write.csv(p_value_matrix_sub_dt, "HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/paired.wilcox.rank.sum.test_pval_table_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_mCA_methFoldDiff_stdevs.csv")
write.csv(t_test_p_value_matrix_sub_dt, "HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/paired.t.test_pval_table_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_mCA_methFoldDiff_stdevs.csv")


log10_t_test_p_value_matrix_sub <- -log10(t_test_p_value_matrix_sub)

ttest_col_palette <- colorRampPalette(c("pink", "red", "darkred"))(n = 299)

png("HG_lab/Mati/GabelLab/cell_confusion_summaryPlots/ridgeline_plots/paired.t.test_negLog10_pval_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_mCA_methFoldDiff_stdevs_heatmap.png", width=1500, height=1500, res=300)
heatmap.2(log10_t_test_p_value_matrix_sub,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=wilcox_col_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, key.title="-Log10 t-test p-value",
          key.par=list(cex=0.5),
          cexRow=0.6,            # adjust row label text size (smaller value for smaller text)
          cexCol=0.6             # adjust column label text size (smaller value for smaller text))
)
dev.off()

setEPS()
postscript("HG_lab/Mati/GabelLab/cell_confusion_revision_plots_and_data/ridgeline_plots/paired.t.test_negLog10_pval_cellType_1kbWindows_intragenic_extragenic_IntraAndExtraLinkedcCRE_INTACT_mCA_methFoldDiff_stdevs_heatmap.eps", width=1500, height=1500)
heatmap.2(log10_t_test_p_value_matrix,
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          #breaks=col_breaks,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          #col=my_palette,       # use on color palette defined earlier
          #col=brewer.pal(9,"Reds"),
          col=wilcox_col_palette,
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",key=TRUE, key.title="-Log10 t-test p-value",
          key.par=list(cex=0.5),
          cexRow=0.6,            # adjust row label text size (smaller value for smaller text)
          cexCol=0.6             # adjust column label text size (smaller value for smaller text))
)
dev.off()