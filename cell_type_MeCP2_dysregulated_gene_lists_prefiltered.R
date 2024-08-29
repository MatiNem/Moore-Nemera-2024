library(data.table)
library(dplyr)
all_genes_mm9 = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_flat.txt")
chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
protein_coding_genes = fread("HG_lab/Mati/GabelLab/genesets/ensgene_mm9_coding_2.bed")
protein_coding_genes_chromList = protein_coding_genes[V1 %in% chrom_list,]

pv_deseq = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/pv_ko_exon_full_dseq_rep_prefilt5_res_df_all_090821.tsv"), keep.rownames = "Gene")
sst_deseq = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/sst_ko_exon_full dseq_rep_prefilt5_res_df_all_090821.tsv"), keep.rownames = "Gene")
L4_deseq = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nr5a1_ko_exon_full_dseq_rep_prefilt5_res_df_all_090821.tsv"), keep.rownames = "Gene")



pv_mr_apeglm = pv_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange > 0), Gene]
pv_mr_ashr = pv_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]
pv_ma_apeglm = pv_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange < 0), Gene]
pv_ma_ashr = pv_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]
pv_unchanged_apeglm = pv_deseq[(apeglm_gWK_pvalue > 0.5), Gene]
pv_unchanged_ashr = pv_deseq[(ashr_pvalue > 0.5), Gene]


sst_mr_apeglm = sst_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange > 0), Gene]
sst_mr_ashr = sst_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]
sst_ma_apeglm = sst_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange < 0), Gene]
sst_ma_ashr = sst_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]
sst_unchanged_apeglm = sst_deseq[(apeglm_gWK_pvalue > 0.5), Gene]
sst_unchanged_ashr = sst_deseq[(ashr_pvalue > 0.5), Gene]

L4_mr_apeglm = L4_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange > 0), Gene]
L4_mr_ashr = L4_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]
L4_ma_apeglm = L4_deseq[(apeglm_gWK_padj <= 0.1) & (apeglm_gWK_log2FoldChange < 0), Gene]
L4_ma_ashr = L4_deseq[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]
L4_unchanged_apeglm = L4_deseq[(apeglm_gWK_pvalue > 0.5), Gene]
L4_unchanged_ashr = L4_deseq[(ashr_pvalue > 0.5), Gene]


Pv_mr_genes_q0.1_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% pv_mr_ashr), ]
Pv_ma_genes_q0.1_proteinCoding_prefilt5_mm9 =  protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% pv_ma_ashr), ]
Pv_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% pv_unchanged_ashr), ]

Sst_mr_genes_q0.1_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% sst_mr_ashr), ]
Sst_ma_genes_q0.1_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% sst_ma_ashr), ]
Sst_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% sst_unchanged_ashr), ]

L4_mr_genes_q0.1_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% L4_mr_ashr), ]
L4_ma_genes_q0.1_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% L4_ma_ashr), ]
L4_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9 = protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% L4_unchanged_ashr), ]

write.table(Pv_mr_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Pv/Pv_mr_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(Pv_ma_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Pv/Pv_ma_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(Pv_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Pv/Pv_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(Sst_mr_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Sst/Sst_mr_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(Sst_ma_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Sst/Sst_ma_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(Sst_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Sst/Sst_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(L4_mr_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/L4/L4_mr_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_ma_genes_q0.1_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/L4/L4_ma_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/L4/L4_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")


Pv_unchanged_SstorL4_mr_genes = intersect(Pv_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9[, V4], c(Sst_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4], L4_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4]))
Sst_unchanged_PvorL4_mr_genes = intersect(Sst_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9[, V4], c(Pv_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4], L4_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4]))
L4_unchanged_PvorSst_mr_genes = intersect(L4_unchanged_genes_p0.5_proteinCoding_prefilt5_mm9[, V4], c(Pv_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4], Sst_mr_genes_q0.1_proteinCoding_prefilt5_mm9[, V4]))


write.table(protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% Pv_unchanged_SstorL4_mr_genes), ], file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Pv/Pv_unchanged_SstorL4_MR_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% Sst_unchanged_PvorL4_mr_genes), ], file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/Sst/Sst_unchanged_PvorL4_MR_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(protein_coding_genes[(V1 %in% chrom_list) & (V4 %in% L4_unchanged_PvorSst_mr_genes), ], file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/L4/L4_unchanged_PvorSst_MR_genes_q0.1_proteinCoding_prefilt5_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")



###Using non-deduplicated files
pv_deseq_nondedup = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/pv_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv"), keep.rownames = "Gene")
sst_deseq_nondedup = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/sst_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv"), keep.rownames = "Gene")
L4_deseq_nondedup = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/nr5a1_ko_exon_nondedup_coding_dseq_rep_prefilt5_res_df_all_110921.tsv"), keep.rownames = "Gene")
L5_deseq_nondedup = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/rbp4_ko_exon_nondedup_coding_d3_s12_dseq_rep_prefilt5_res_df_all_110921.tsv"), keep.rownames = "Gene")

pv_mr_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
pv_ma_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
pv_unchanged_genes_p0.5_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup[(ashr_pvalue > 0.5), Gene]]

sst_mr_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
sst_ma_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
sst_unchanged_genes_p0.5_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup[(ashr_pvalue > 0.5), Gene]]

L4_mr_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
L4_ma_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
L4_unchanged_genes_p0.5_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup[(ashr_pvalue > 0.5), Gene]]

L5_mr_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
L5_ma_genes_q0.1_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
L5_unchanged_genes_p0.5_nondedup_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup[(ashr_pvalue > 0.5), Gene]]
#Other-cell-type MR genes
pv_otherCellType_mr_genes_q0.1_nondedup_mm9 = pv_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_mm9[,V4], L4_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]), ]
sst_otherCellType_mr_genes_q0.1_nondedup_mm9 = sst_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(pv_mr_genes_q0.1_nondedup_mm9[,V4], L4_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]), ]
L4_otherCellType_mr_genes_q0.1_nondedup_mm9 = L4_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_mm9[,V4], pv_mr_genes_q0.1_nondedup_mm9[,V4], L5_mr_genes_q0.1_nondedup_mm9[,V4]), ]
L5_otherCellType_mr_genes_q0.1_nondedup_mm9 = L5_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_mm9[,V4], pv_mr_genes_q0.1_nondedup_mm9[,V4], L4_mr_genes_q0.1_nondedup_mm9[,V4]), ]

#Other-cell-type MA genes
pv_otherCellType_ma_genes_q0.1_nondedup_mm9 = pv_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_mm9[,V4], L4_ma_genes_q0.1_nondedup_mm9[,V4], L5_ma_genes_q0.1_nondedup_mm9[,V4]), ]
sst_otherCellType_ma_genes_q0.1_nondedup_mm9 = sst_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(pv_ma_genes_q0.1_nondedup_mm9[,V4], L4_ma_genes_q0.1_nondedup_mm9[,V4], L5_ma_genes_q0.1_nondedup_mm9[,V4]), ]
L4_otherCellType_ma_genes_q0.1_nondedup_mm9 = L4_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_mm9[,V4], pv_ma_genes_q0.1_nondedup_mm9[,V4], L5_ma_genes_q0.1_nondedup_mm9[,V4]), ]
L5_otherCellType_ma_genes_q0.1_nondedup_mm9 = L5_unchanged_genes_p0.5_nondedup_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_mm9[,V4], pv_ma_genes_q0.1_nondedup_mm9[,V4], L4_ma_genes_q0.1_nondedup_mm9[,V4]), ]


write.table(pv_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_unchanged_genes_p0.5_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_otherCellType_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_otherCellType_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(sst_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_unchanged_genes_p0.5_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_otherCellType_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_otherCellType_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/Sst/Sst_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(L4_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_unchanged_genes_p0.5_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_otherCellType_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_otherCellType_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L4/L4_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(L5_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_unchanged_genes_p0.5_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_otherCellType_mr_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_otherCellType_ma_genes_q0.1_nondedup_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/L5/L5_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")


###using new RUVG normalization
pv_deseq_nondedup_new_ruvg = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/new_ruvg/pv_ko_nondedup_deseq_cod_pfilt5_ruvg_emp10o_011722.tsv"), keep.rownames = "Gene")
sst_deseq_nondedup_new_ruvg = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/new_ruvg/sst_ko_nondedup_deseq_cod_pfilt5_ruvg_emp10o_011722.tsv"), keep.rownames = "Gene")
L4_deseq_nondedup_new_ruvg = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/new_ruvg/nr5a1_ko_nondedup_deseq_cod_pfilt5_ruvg_emp10o_011722.tsv"), keep.rownames = "Gene")
L5_deseq_nondedup_new_ruvg = data.table(read.table("HG_lab/Mati/GabelLab/genesets/deseq_files/nondedup/new_ruvg/rbp4_ko_nondedup_deseq_cod_pfilt5_ruvg_emp10o_011722.tsv"), keep.rownames = "Gene")

pv_mr_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
pv_ma_genes_q0.1_nondedup_new_ruvg_mm9 =  protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
pv_unchanged_genes_p0.5_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% pv_deseq_nondedup_new_ruvg[(ashr_pvalue > 0.5), Gene]]

sst_mr_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
sst_ma_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
sst_unchanged_genes_p0.5_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% sst_deseq_nondedup_new_ruvg[(ashr_pvalue > 0.5), Gene]]

L4_mr_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
L4_ma_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
L4_unchanged_genes_p0.5_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L4_deseq_nondedup_new_ruvg[(ashr_pvalue > 0.5), Gene]]

L5_mr_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange > 0), Gene]]
L5_ma_genes_q0.1_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup_new_ruvg[(ashr_padj <= 0.1) & (ashr_log2FoldChange < 0), Gene]]
L5_unchanged_genes_p0.5_nondedup_new_ruvg_mm9 = protein_coding_genes_chromList[V4 %in% L5_deseq_nondedup_new_ruvg[(ashr_pvalue > 0.5), Gene]]
#Other-cell-type MR genes
pv_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9 = pv_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
sst_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9 = sst_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(pv_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
L4_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9 = L4_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], pv_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
L5_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9 = L5_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], pv_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_mr_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]

#Other-cell-type MA genes
pv_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9 = pv_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
sst_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9 = sst_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(pv_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
L4_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9 = L4_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], pv_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L5_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]
L5_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9 = L5_unchanged_genes_p0.5_nondedup_new_ruvg_mm9[V4 %in% c(sst_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], pv_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4], L4_ma_genes_q0.1_nondedup_new_ruvg_mm9[,V4]), ]


write.table(pv_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Pv/Pv_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Pv/Pv_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_unchanged_genes_p0.5_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Pv/Pv_unchanged_genes_p0.5_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Pv/Pv_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(pv_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Pv/Pv_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(sst_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Sst/Sst_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Sst/Sst_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_unchanged_genes_p0.5_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Sst/Sst_unchanged_genes_p0.5_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Sst/Sst_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sst_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/Sst/Sst_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(L4_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L4/L4_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L4/L4_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_unchanged_genes_p0.5_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L4/L4_unchanged_genes_p0.5_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L4/L4_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L4_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L4/L4_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

write.table(L5_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L5/L5_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L5/L5_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_unchanged_genes_p0.5_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L5/L5_unchanged_genes_p0.5_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_otherCellType_mr_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L5/L5_otherCellType_mr_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")
write.table(L5_otherCellType_ma_genes_q0.1_nondedup_new_ruvg_mm9, file="HG_lab/Mati/GabelLab/genesets/prefilt5_genes/nondedup/new_ruvg/L5/L5_otherCellType_ma_genes_q0.1_coding_prefilt5_nondedup_new_ruvg_mm9.bed", quote=F, col.names=F, row.names=F, sep="\t")

