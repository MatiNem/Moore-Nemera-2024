library(data.table)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
group1 = fread(args[1], header=FALSE)
group1_name = args[2]
group1_col = args[3]

group2 = fread(args[4], header=FALSE)
group2_name = args[5]
group2_col = args[6]

#group3 = fread(args[7], header=FALSE)
#group3_name = args[8]
#group3_col = args[9]

methylation_type = args[7]
min_y = as.numeric(args[8]) 
max_y = as.numeric(args[9])

output = args[10]


names(group1) = c("bin_number", "bin_label", "mean_methylation")
group1[, gene_list := group1_name]

names(group2) = c("bin_number", "bin_label", "mean_methylation")
group2[, gene_list := group2_name]

#names(group3) = c("bin_number", "bin_label", "mean_methylation")
#group3[, gene_list := group3_name]

full = rbind(group1, group2)

name_vec = c(group1_name, group2_name)
col_vec = c(group1_col, group2_col)
colors = setNames(col_vec, name_vec)


h <- ggplot(full ,aes(x = bin_number, y = mean_methylation, group = gene_list))
h +
  coord_cartesian(ylim=c(min_y,max_y))+
  geom_rect(aes(xmin = -23, xmax = 27, ymin = -Inf, ymax = Inf), fill="gray", alpha = 0.5) +
  geom_text(x=1.5, y=0.01, label="Metagene") +
  geom_line(aes(color = gene_list)) +
  scale_color_manual(values = colors) +
  xlab("") + ylab(paste0("Mean m",methylation_type,"/",methylation_type)) +
  geom_vline(xintercept=-23,col="lightgray") + geom_vline(xintercept=26,col="gray") +
  scale_x_continuous(breaks=c(-128,-28,-23,27,127),labels=c("-100kb","TSS","+5kb","TES","+100kb"), limits = c(-128,128)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.margin=margin(), legend.title=element_blank(), panel.border=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(color="black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

ggsave(paste0(output,".png"), device="png")
ggsave(paste0(output,".eps"), device=cairo_ps)
