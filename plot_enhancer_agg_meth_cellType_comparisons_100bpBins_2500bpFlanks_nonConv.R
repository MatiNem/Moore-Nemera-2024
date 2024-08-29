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

y_axis_title = args[10]
min_y = as.numeric(args[11])
max_y = as.numeric(args[12])
#upper and lower bounds of rectangle, usually (rectangle_max - rectangle_min) is the median length of the enhancer
rectangle_min = as.numeric(args[13])
rectangle_max = as.numeric(args[14])
#output file
output = args[15]

#non-conversion rate table
nonconv_table <- fread(args[16])
nonconv_label <- args[17]
nonconv_rate <- nonconv_table[label==nonconv_label, as.numeric(nonconversion_rate)]


names(group1) = c("bin_number", "bin_label", "mean_methylation")
group1[, gene_list := group1_name]
group1$mean_methylation <- as.numeric(group1$mean_methylation)
group1$corrected_mean_methylation <- group1$mean_methylation - nonconv_rate
print(median(group1$corrected_mean_methylation, na.rm=TRUE))

names(group2) = c("bin_number", "bin_label", "mean_methylation")
group2[, gene_list := group2_name]
group2$mean_methylation <- as.numeric(group2$mean_methylation)
group2$corrected_mean_methylation <- group2$mean_methylation - nonconv_rate
print(median(group2$corrected_mean_methylation, na.rm=TRUE))

names(group3) = c("bin_number", "bin_label", "mean_methylation")
group3[, gene_list := group3_name]
group3$mean_methylation <- as.numeric(group3$mean_methylation)
group3$corrected_mean_methylation <- group3$mean_methylation - nonconv_rate
print(median(group3$corrected_mean_methylation, na.rm=TRUE))

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
  geom_rect(aes(xmin = rectangle_min, xmax = rectangle_max, ymin = -Inf, ymax = Inf), fill="gray", alpha = 0.5) +
  geom_line(aes(color = gene_list)) +
  #scale_color_manual(name="Enhancers linked to:", values = colors) +
  scale_color_manual(name="cCREs linked to:", values = colors) +
  xlab("") + ylab(paste0(y_axis_title)) +
  scale_x_continuous(breaks=c(-25,-12.5,0,12.5,25),labels=c("-2.5kb","","0","","+2.5kb"), limits = c(-25,25)) +
  theme_bw() +
  theme(legend.text= element_text(size = 10), legend.position = "bottom", legend.margin=margin(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(color="black"), axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_text(size=15), axis.text.y=element_text(size=15))+
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

ggsave(paste0(output,".png"), device="png")
ggsave(paste0(output,".eps"), device=cairo_ps)
