library(data.table)
library(dplyr)

##mCA equivalent of combined flank and body methylation tables
PV_INTACT_gene_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
L4_INTACT_gene_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
SST_INTACT_gene_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
L5_INTACT_gene_flank_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_L5_WT_KO_deep_INTACT_mCA_mm9.bed")


names(PV_INTACT_gene_flank_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(L4_INTACT_gene_flank_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(SST_INTACT_gene_flank_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(L5_INTACT_gene_flank_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")


#add up reads of each pair of flanks per gene; remove flanks with missing values before adding
PV_INTACT_gene_flank_mCA = PV_INTACT_gene_flank_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
L4_INTACT_gene_flank_mCA = L4_INTACT_gene_flank_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
SST_INTACT_gene_flank_mCA = SST_INTACT_gene_flank_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
L5_INTACT_gene_flank_mCA = L5_INTACT_gene_flank_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]

PV_INTACT_gene_flank_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
L4_INTACT_gene_flank_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
SST_INTACT_gene_flank_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
L5_INTACT_gene_flank_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]


#gene body (TSS + 3kb to TES) INTACT mCA
PV_gene_body_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_WT_KO_deep_INTACT_mCA_mm9.bed")
L4_gene_body_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_WT_KO_deep_INTACT_mCA_mm9.bed")
SST_gene_body_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_WT_KO_deep_INTACT_mCA_mm9.bed")
L5_gene_body_TSSplus3kb_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_WT_KO_deep_INTACT_mCA_mm9.bed")

names(PV_gene_body_TSSplus3kb_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L4_gene_body_TSSplus3kb_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(SST_gene_body_TSSplus3kb_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L5_gene_body_TSSplus3kb_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")

PV_gene_body_TSSplus3kb_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L4_gene_body_TSSplus3kb_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
SST_gene_body_TSSplus3kb_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L5_gene_body_TSSplus3kb_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

PV_gene_and_flank_mCA = data.table(left_join(x=PV_gene_body_TSSplus3kb_mCA, y=PV_INTACT_gene_flank_mCA, by=c("gene")))
L4_gene_and_flank_mCA = data.table(left_join(x=L4_gene_body_TSSplus3kb_mCA, y=L4_INTACT_gene_flank_mCA, by=c("gene")))
SST_gene_and_flank_mCA = data.table(left_join(x=SST_gene_body_TSSplus3kb_mCA, y=SST_INTACT_gene_flank_mCA, by=c("gene")))
L5_gene_and_flank_mCA = data.table(left_join(x=L5_gene_body_TSSplus3kb_mCA, y=L5_INTACT_gene_flank_mCA, by=c("gene")))

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")
#subtracting nonconversion rate from gene body and flanking methylation values
PV_gene_and_flank_mCA$gene_methylation_corrected <- PV_gene_and_flank_mCA$gene_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
PV_gene_and_flank_mCA$flank_methylation_corrected <- PV_gene_and_flank_mCA$flank_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]

L4_gene_and_flank_mCA$gene_methylation_corrected <- L4_gene_and_flank_mCA$gene_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L4_gene_and_flank_mCA$flank_methylation_corrected <- L4_gene_and_flank_mCA$flank_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]

L5_gene_and_flank_mCA$gene_methylation_corrected <- L5_gene_and_flank_mCA$gene_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
L5_gene_and_flank_mCA$flank_methylation_corrected <- L5_gene_and_flank_mCA$flank_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

SST_gene_and_flank_mCA$gene_methylation_corrected <- SST_gene_and_flank_mCA$gene_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
SST_gene_and_flank_mCA$flank_methylation_corrected <- SST_gene_and_flank_mCA$flank_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]


write.table(PV_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L4_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(SST_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L5_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")



##mCG equivalent of combined flank and body methylation tables
PV_INTACT_gene_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
L4_INTACT_gene_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
SST_INTACT_gene_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
L5_INTACT_gene_flank_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/all_genes_flanks200kb_geneMasked_L5_WT_KO_deep_INTACT_mCG_mm9.bed")


names(PV_INTACT_gene_flank_mCG) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(L4_INTACT_gene_flank_mCG) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(SST_INTACT_gene_flank_mCG) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(L5_INTACT_gene_flank_mCG) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")


#add up reads of each pair of flanks per gene; remove flanks with missing values before adding
PV_INTACT_gene_flank_mCG = PV_INTACT_gene_flank_mCG[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
L4_INTACT_gene_flank_mCG = L4_INTACT_gene_flank_mCG[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
SST_INTACT_gene_flank_mCG = SST_INTACT_gene_flank_mCG[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
L5_INTACT_gene_flank_mCG = L5_INTACT_gene_flank_mCG[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]

PV_INTACT_gene_flank_mCG[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
L4_INTACT_gene_flank_mCG[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
SST_INTACT_gene_flank_mCG[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
L5_INTACT_gene_flank_mCG[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]


#gene body (TSS + 3kb to TES) INTACT mCG
PV_gene_body_TSSplus3kb_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_PV_WT_KO_deep_INTACT_mCG_mm9.bed")
L4_gene_body_TSSplus3kb_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L4_WT_KO_deep_INTACT_mCG_mm9.bed")
SST_gene_body_TSSplus3kb_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_SST_WT_KO_deep_INTACT_mCG_mm9.bed")
L5_gene_body_TSSplus3kb_mCG = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/ensgene_mm9_TSS_plus_3kb_L5_WT_KO_deep_INTACT_mCG_mm9.bed")

names(PV_gene_body_TSSplus3kb_mCG) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L4_gene_body_TSSplus3kb_mCG) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(SST_gene_body_TSSplus3kb_mCG) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(L5_gene_body_TSSplus3kb_mCG) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")

PV_gene_body_TSSplus3kb_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L4_gene_body_TSSplus3kb_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
SST_gene_body_TSSplus3kb_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
L5_gene_body_TSSplus3kb_mCG[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

PV_gene_and_flank_mCG = data.table(left_join(x=PV_gene_body_TSSplus3kb_mCG, y=PV_INTACT_gene_flank_mCG, by=c("gene")))
L4_gene_and_flank_mCG = data.table(left_join(x=L4_gene_body_TSSplus3kb_mCG, y=L4_INTACT_gene_flank_mCG, by=c("gene")))
SST_gene_and_flank_mCG = data.table(left_join(x=SST_gene_body_TSSplus3kb_mCG, y=SST_INTACT_gene_flank_mCG, by=c("gene")))
L5_gene_and_flank_mCG = data.table(left_join(x=L5_gene_body_TSSplus3kb_mCG, y=L5_INTACT_gene_flank_mCG, by=c("gene")))

#bisulfite non-conversion rates
avg_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/lambda_average_nonconversion_table.tsv")
#subtracting nonconversion rate from gene body and flanking methylation values
PV_gene_and_flank_mCG$gene_methylation_corrected <- PV_gene_and_flank_mCG$gene_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]
PV_gene_and_flank_mCG$flank_methylation_corrected <- PV_gene_and_flank_mCG$flank_methylation - avg_nonconv[label=="PV_WT_KO_lambda", nonconversion_rate]

L4_gene_and_flank_mCG$gene_methylation_corrected <- L4_gene_and_flank_mCG$gene_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]
L4_gene_and_flank_mCG$flank_methylation_corrected <- L4_gene_and_flank_mCG$flank_methylation - avg_nonconv[label=="L4_WT_KO_lambda", nonconversion_rate]

L5_gene_and_flank_mCG$gene_methylation_corrected <- L5_gene_and_flank_mCG$gene_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]
L5_gene_and_flank_mCG$flank_methylation_corrected <- L5_gene_and_flank_mCG$flank_methylation - avg_nonconv[label=="L5_WT_KO_lambda", nonconversion_rate]

SST_gene_and_flank_mCG$gene_methylation_corrected <- SST_gene_and_flank_mCG$gene_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]
SST_gene_and_flank_mCG$flank_methylation_corrected <- SST_gene_and_flank_mCG$flank_methylation - avg_nonconv[label=="SST_WT_KO_lambda", nonconversion_rate]


write.table(PV_gene_and_flank_mCG, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L4_gene_and_flank_mCG, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L4_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(SST_gene_and_flank_mCG, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/SST_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L5_gene_and_flank_mCG, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_WT_KO_deep_INTACT_gene_body_TSSplus3kb_flank_mCG_mm9.txt", quote=F, row.names=F, sep="\t")



###WT and KO replicates
##mCA flank and body methylation tables
#PV
all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA_mm9.bed")

#L5
all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA_mm9.bed")
all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA_mm9.bed")


names(all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")

names(all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")
names(all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA) = c("chrom", "flank_start", "flank_end", "gene", "strand", "flank_meth_reads", "flank_cyto_reads")


#add up reads of each pair of flanks per gene; remove flanks with missing values before adding
#PV
all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
#L5
all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]
all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA = all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA[is.finite(flank_meth_reads) & is.finite(flank_cyto_reads), lapply(.SD, sum), by = gene, .SDcols = c("flank_meth_reads", "flank_cyto_reads")]


#flank methylation level calculation
all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]

all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]
all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA[, flank_methylation := as.integer(flank_meth_reads)/as.integer(flank_cyto_reads)]


#gene body (TSS + 3kb to TES) INTACT mCA
#PV
ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041642_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041642_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041644_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041644_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041641_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041641_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041643_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041643_deep_INTACT_mCA_mm9.bed")

#L5
ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041645_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041645_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041648_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041648_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041646_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041646_deep_INTACT_mCA_mm9.bed")
ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041647_deep_INTACT_mCA = fread("HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041647_deep_INTACT_mCA_mm9.bed")


names(ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041642_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041644_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041641_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041643_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")

names(ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041645_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041648_deep_INTACT_mCA ) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041646_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")
names(ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041647_deep_INTACT_mCA) = c("chrom", "start", "end", "gene", "transcripts", "strand", "meth_reads", "cyto_reads")


ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041642_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041644_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041641_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041643_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]

ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041645_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041648_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041646_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]
ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041647_deep_INTACT_mCA[, gene_methylation := as.integer(meth_reads)/as.integer(cyto_reads)]



PV_WT_LIB041642_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041642_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_PV_WT_LIB041642_deep_INTACT_mCA, by=c("gene")))
PV_WT_LIB041644_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_PV_WT_LIB041644_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_PV_WT_LIB041644_deep_INTACT_mCA, by=c("gene")))
PV_KO_LIB041641_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041641_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_PV_KO_LIB041641_deep_INTACT_mCA, by=c("gene")))
PV_KO_LIB041643_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_PV_KO_LIB041643_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_PV_KO_LIB041643_deep_INTACT_mCA, by=c("gene")))


L5_WT_LIB041645_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041645_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_L5_WT_LIB041645_deep_INTACT_mCA, by=c("gene")))
L5_WT_LIB041648_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_L5_WT_LIB041648_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_L5_WT_LIB041648_deep_INTACT_mCA, by=c("gene")))
L5_KO_LIB041646_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041646_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_L5_KO_LIB041646_deep_INTACT_mCA, by=c("gene")))
L5_KO_LIB041647_gene_and_flank_mCA = data.table(left_join(x=ensgene_mm9_TSS_plus_3kb_L5_KO_LIB041647_deep_INTACT_mCA, y=all_genes_flanks200kb_geneMasked_L5_KO_LIB041647_deep_INTACT_mCA, by=c("gene")))

#bisulfite non-conversion rates, all
lambda_nonconv=fread("HG_lab/Mati/GabelLab/bisulfite/summary_tsv/intact_lambda_nonconversion_table.csv")

#subtracting nonconversion rate from gene body and flanking methylation values
#PV
PV_WT_LIB041642_gene_and_flank_mCA$gene_methylation_corrected <- PV_WT_LIB041642_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041642", nonconversion_rate]
PV_WT_LIB041642_gene_and_flank_mCA$flank_methylation_corrected <- PV_WT_LIB041642_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041642", nonconversion_rate]

PV_WT_LIB041644_gene_and_flank_mCA$gene_methylation_corrected <- PV_WT_LIB041644_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041644", nonconversion_rate]
PV_WT_LIB041644_gene_and_flank_mCA$flank_methylation_corrected <- PV_WT_LIB041644_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041644", nonconversion_rate]

PV_KO_LIB041641_gene_and_flank_mCA$gene_methylation_corrected <- PV_KO_LIB041641_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041641", nonconversion_rate]
PV_KO_LIB041641_gene_and_flank_mCA$flank_methylation_corrected <- PV_KO_LIB041641_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041641", nonconversion_rate]

PV_KO_LIB041643_gene_and_flank_mCA$gene_methylation_corrected <- PV_KO_LIB041643_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041643", nonconversion_rate]
PV_KO_LIB041643_gene_and_flank_mCA$flank_methylation_corrected <- PV_KO_LIB041643_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041643", nonconversion_rate]

#L5
L5_WT_LIB041645_gene_and_flank_mCA$gene_methylation_corrected <- L5_WT_LIB041645_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041645", nonconversion_rate]
L5_WT_LIB041645_gene_and_flank_mCA$flank_methylation_corrected <- L5_WT_LIB041645_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041645", nonconversion_rate]

L5_WT_LIB041648_gene_and_flank_mCA$gene_methylation_corrected <- L5_WT_LIB041648_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041648", nonconversion_rate]
L5_WT_LIB041648_gene_and_flank_mCA$flank_methylation_corrected <- L5_WT_LIB041648_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041648", nonconversion_rate]

L5_KO_LIB041646_gene_and_flank_mCA$gene_methylation_corrected <- L5_KO_LIB041646_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041646", nonconversion_rate]
L5_KO_LIB041646_gene_and_flank_mCA$flank_methylation_corrected <- L5_KO_LIB041646_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041646", nonconversion_rate]

L5_KO_LIB041647_gene_and_flank_mCA$gene_methylation_corrected <- L5_KO_LIB041647_gene_and_flank_mCA$gene_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041647", nonconversion_rate]
L5_KO_LIB041647_gene_and_flank_mCA$flank_methylation_corrected <- L5_KO_LIB041647_gene_and_flank_mCA$flank_methylation - lambda_nonconv[GTAC_ESP_ID=="LIB041647", nonconversion_rate]



write.table(PV_WT_LIB041642_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_WT_LIB041642_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(PV_WT_LIB041644_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_WT_LIB041644_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(PV_KO_LIB041641_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_KO_LIB041641_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(PV_KO_LIB041643_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/PV_reps/PV_KO_LIB041643_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")

write.table(L5_WT_LIB041645_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_WT_LIB041645_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L5_WT_LIB041648_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_WT_LIB041648_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L5_KO_LIB041646_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_KO_LIB041646_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")
write.table(L5_KO_LIB041647_gene_and_flank_mCA, file="HG_lab/Mati/GabelLab/genesets/INTACT_meth_genes/L5_reps/L5_KO_LIB041647_deep_INTACT_gene_body_TSSplus3kb_flank_mCA_mm9.txt", quote=F, row.names=F, sep="\t")

