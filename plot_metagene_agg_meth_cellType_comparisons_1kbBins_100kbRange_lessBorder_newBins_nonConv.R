library(data.table)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
group1 = fread(args[1], header=FALSE)
group1_name = args[2]
group1_col = args[3]

group2 = fread(args[4], header=FALSE)
group2_name = args[5]
group2_col = args[6]

group3 = fread(args[7], header=FALSE)
group3_name = args[8]
group3_col = args[9]

methylation_type = args[10]
min_y = as.numeric(args[11]) 
max_y = as.numeric(args[12])

output = args[13]

#non-conversion rate table
nonconv_table <- fread(args[14])
nonconv_label <- args[15]
nonconv_rate <- nonconv_table[label==nonconv_label, as.numeric(nonconversion_rate)]


names(group1) = c("bin_number", "bin_label", "mean_methylation")
group1[, gene_list := group1_name]
group1$corrected_mean_methylation <- group1$mean_methylation - nonconv_rate 

names(group2) = c("bin_number", "bin_label", "mean_methylation")
group2[, gene_list := group2_name]
group2$corrected_mean_methylation <- group2$mean_methylation - nonconv_rate

names(group3) = c("bin_number", "bin_label", "mean_methylation")
group3[, gene_list := group3_name]
group3$corrected_mean_methylation <- group3$mean_methylation - nonconv_rate

full = rbind(group1, group2, group3)

name_vec = c(group1_name, group2_name, group3_name)
col_vec = c(group1_col, group2_col, group3_col)
colors = setNames(col_vec, name_vec)

filename1 = paste0(as.character(args[1]), "_nonConvSub")
filename2 = paste0(as.character(args[4]), "_nonConvSub")
filename3 = paste0(as.character(args[7]), "_nonConvSub")

write.table(group1, file=filename1, quote=FALSE, row.names=FALSE)
write.table(group2, file=filename2, quote=FALSE, row.names=FALSE)
write.table(group3, file=filename3, quote=FALSE, row.names=FALSE)

h <- ggplot(full ,aes(x = bin_number, y = corrected_mean_methylation, group = gene_list))
h +
  coord_cartesian(ylim=c(min_y,max_y))+
  geom_rect(aes(xmin = -22.5, xmax = 27.5, ymin = -Inf, ymax = Inf), fill="gray", alpha = 0.5) +
  geom_text(x=1.5, y=0.01, label="Metagene") +
  geom_line(aes(color = gene_list)) +
  scale_color_manual(values = colors) +
  xlab("") + ylab(paste0("Mean m",methylation_type,"/",methylation_type)) +
  #geom_vline(xintercept=-22.5,col="lightgray") + geom_vline(xintercept=27.5,col="gray") +
  scale_x_continuous(breaks=c(-127.5,-27.5,-22.5,27.5,127.5),labels=c("-100kb","TSS","+5kb","TES","+100kb"), limits = c(-127.5,127.5)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.margin=margin(), legend.title=element_blank(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(color="black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

ggsave(paste0(output,".png"), device="png")
ggsave(paste0(output,".eps"), device=cairo_ps)


#write.table(group2, file=paste0(as.character(args[4]),"_sub"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#write.table(group3, file=paste0(as.character(args[7]),"_sub"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


