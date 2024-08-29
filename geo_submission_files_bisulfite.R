library(dplyr)

# Set the directory path
#PV
folder <- "HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/PV"

# List all files in the directory
all_files <- list.files(folder, full.names = FALSE)

# Filter the files to get those ending with "R1_001.fastq"
r1_files <- all_files[grep("R1_001.fastq", all_files)]

# For each R1 file, find the corresponding R2 file
paired_files <- lapply(r1_files, function(r1_file) {
  r2_file <- sub("R1_001.fastq", "R2_001.fastq", r1_file)
  if (r2_file %in% all_files) {
    return(data.frame(R1 = r1_file, R2 = r2_file))
  } else {
    return(NULL)
  }
})

# Combine the pairs into a single data frame
paired_files_df <- do.call(rbind, paired_files)

# Print the two-column table
print(paired_files_df)

#
r1_r2_function <- function(folder_path){
  folder <- folder_path
  
  # List all files in the directory
  all_files <- list.files(folder, full.names = FALSE)
  
  # Filter the files to get those ending with "R1_001.fastq"
  r1_files <- all_files[grep("R1_001.fastq", all_files)]
  
  # For each R1 file, find the corresponding R2 file
  paired_files <- lapply(r1_files, function(r1_file) {
    r2_file <- sub("R1_001.fastq", "R2_001.fastq", r1_file)
    if (r2_file %in% all_files) {
      return(data.frame(R1 = r1_file, R2 = r2_file))
    } else {
      return(NULL)
    }
  })
  
  paired_files_df <- do.call(rbind, paired_files)
  
  return(paired_files_df)
}

View(r1_r2_function("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/PV"))
View(r1_r2_function("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/SST"))
View(r1_r2_function("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/L4"))
View(r1_r2_function("HG_lab/Manuscripts/Cellular_confusion/GEO_submission/bisulfite/L5"))