library(data.table)
library(dplyr)
library(stringr)

tissue_meths1 <- fread("HG_lab/NGS/SIC_Submission_194_YR/reports/ATCGmap/194_YR_script_copy.tsv")

tissue_meths2 <- fread("HG_lab/NGS/SIC_Submission_196_YR/Gabel_196A/reports/ATCGmap/Gabel_196A_script.tsv")


tissue_meths_all <- left_join(tissue_meths1, tissue_meths2, by="V1")

names(tissue_meths_all)[1]= "nucleotides"

#
update_colnames <- function(dt) {
  new_names <- sapply(names(dt), function(col_name) {
    if (str_detect(col_name, "(506-2|598-2)_(CB|HYPO|STR)")) {
      return(str_replace(col_name, "^(.*)_(506-2|598-2)_(CB|HYPO|STR)(.*)$", "\\3_KO_\\1_\\2_\\3\\4"))
    } else if (str_detect(col_name, "(506-3|598-1)_(CB|HYPO|STR)")) {
      return(str_replace(col_name, "^(.*)_(506-3|598-1)_(CB|HYPO|STR)(.*)$", "\\3_WT_\\1_\\2_\\3\\4"))
    } else {
      return(col_name)
    }
  })
  setnames(dt, old = names(dt), new = new_names)
}

# Example data.table for testing
dt <- data.table(
  "Gabel_3_506-2_HYPO_D703_D503CGCTCATTAT_AGGATAGG_S56_summary" = 1:5,
  "Gabel_17_598-1_HYPO_D703_D503_CGCTCATTAT_AGGATAGG_S99_summary" = 6:10,
  "Other_Column" = 11:15
)

# Apply the function
update_colnames(dt)

# Print updated data.table
print(dt)



update_colnames(tissue_meths_all)


write.table(tissue_meths_all, file="HG_lab/Mati/GabelLab/bisulfite/summary_tsv/brain_tissue_summary_global_meth.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
