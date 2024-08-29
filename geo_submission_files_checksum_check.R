library(data.table)
library(dplyr)

mismatched_checks <- fread("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/mismatched_files", header=FALSE)


##
new_checks <- fread("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/new_checksums.md5", header=FALSE)
old_checks <- fread("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/checksums.md5", header=FALSE)


remove_leading_backslashes <- function(x) {
  sub("^\\\\+", "", x)
}

# Test the function with the example string
test_string <- "\\1824ad06dc893b838b1d8fe29d8baa2b"
cleaned_string <- remove_leading_backslashes(test_string)

# Create a new column with the double backslashes removed
new_checks[, V2_clean := remove_leading_backslashes(V2)]

#
# Function to replace backslashes with forward slashes
replace_backslashes <- function(path) {
  gsub("\\\\", "/", path)
}

# Test the function with the example string
test_string <- ".\\CB\\CB_KO_Gabel_20_598-2_CB_D706_D506_GAATTCGTAT_TAAGATTA_S102_R1_001.fastq.gz"
modified_string <- replace_backslashes(test_string)


new_checks[, V1_clean := replace_backslashes(V1)]


#
old_new_checks_join <- inner_join(x=old_checks, y=new_checks[, .(V1_clean, V2_clean)], by=c("V1"="V1_clean", "V2"="V2_clean")) %>% data.table

old_checks$V2 %in% new_checks$V2_clean


remove_leading_path <- function(path) {
  sub("^\\.\\/[A-Za-z0-9]+\\/", "", path)
}

remove_leading_path("./STR/STR_WT_Gabel_2_506-3_STR_D702_D502TCCGGAGAAT_GCCTCTAT_S55_R2_001.fastq.gz")
remove_leading_path("./L4/L4_KO_LIB041633_227WGJLT4_S206_L004_R2_001.fastq.gz")

old_new_checks_join[, V1 := remove_leading_path(V1)]
