#######################
##### Extract significant SNPs
#######################

# Specify the main directory path
main_directory <- "220_phenotype_download/unzip"

# Get a list of sub-directories within the main directory
subdirectories <- list.dirs(main_directory, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store the significant SNPs
BBJ_screen_POS <- list()
count=0

# Loop through each sub-directory
for (subdir in subdirectories) {
  if (length(grep("BBJ",subdir)) > 0){
    # Get a list of all files with .gz extension within the subdirectory
    file_list <- list.files(subdir, pattern = "\\.gz$", full.names = TRUE)
    
    # Loop through each .gz file in the sub-directory
    for (file in file_list) {
      # Check if the file name contains "chrX"
      if (length(grep("auto",file)) > 0) {
        # Unzip the file and read its contents
        unzipped_data <- data.table::fread(file, sep="\t", header=TRUE, select=c("CHR","POS", "p.value.NA", "SE"))
        unzipped_data <- tidyr::unite(unzipped_data, CHR_POS, c(CHR, POS))
        
        # Combine the chromosome and the position information
        CHR_POS <- unzipped_data$CHR_POS
        
        # Extract the significant SNPs
        ##################
        ##### p-value < 1e-3 && standard error < 0.2
        ##################
        
        POS_screen <- unzipped_data[unzipped_data$p.value.NA < 1e-3 & unzipped_data$SE < 0.2,]$CHR_POS
        BBJ_screen_POS <- append(BBJ_screen_POS, list(POS_screen))
        
        count = count + 1
        print(count)
      }
    }
  }
}

file_path <- "220_phenotype_download/processed/BBJ_auto_screen_POS_e3.RDS"
saveRDS(BBJ_screen_POS, file = file_path)

####### Note: here we only provide the example of Japanese population dataset, European population dataset is handled in the same way


#######################
##### Find SNPs that are significant in both European and Japanese datasets
#######################

BBJ_auto_screen <- readRDS(file="BBJ_auto_screen_POS_e3.RDS")

EUR_auto_screen <- readRDS(file="EUR_auto_screen_POS_e3.RDS")

UKB_auto_screen <- readRDS(file="UKB_auto_screen_POS_e3.RDS")

##### Japanese dataset
# 1624000 SNPs
BBJ_unlist <- unlist(BBJ_auto_screen)
# 1217828 SNPs
BBJ_unlist_unique <- unique(BBJ_unlist)

##### European dataset
# 2674267 SNPs
EUR_unlist <- unlist(EUR_auto_screen)
# 1875832 SNPs
EUR_unlist_unique <- unique(EUR_unlist)

# 261515 SNPs
UKB_unlist <- unlist(UKB_auto_screen)
# 233156 SNPs
UKB_unlist_unique <- unique(UKB_unlist)

# 2043115 SNPs
EUR_all_unlist_unique <- union(EUR_unlist_unique, UKB_unlist_unique)

SNP_auto_list <- intersect(BBJ_unlist_unique, EUR_all_unlist_unique)

# Extract the chromosome and position
snp_info = strsplit(SNP_auto_list, "_")
chromosome = sapply(snp_info, "[[", 1)
position = as.numeric(sapply(snp_info, "[[", 2))

# Filter the SNPs in the MHC region
mhc_snps = SNP_auto_list[chromosome == "6" & position >= 25000000 & position <= 34000000]

# Remove the MHC SNPs from the original list
SNP_auto_EHC_rm_list <- SNP_auto_list[!SNP_auto_list %in% mhc_snps]

saveRDS(SNP_auto_EHC_rm_list, file="SNP_auto_EHC_rm_POS_e3.RDS")


#######################
##### Ancestry Information Extraction
#######################

ancestry_info <- read.delim("220_phenotype/reference/Ancestry_info/1000-Geno-ancestry-info.tsv", header=TRUE, sep="\t")
Eastern_Asian_info <- ancestry_info[ancestry_info$Superpopulation.code=="EAS" | grepl("East Asia", ancestry_info$Superpopulation.name), ]
Eastern_Asian_number <- Eastern_Asian_info[,1]

EAS_ids <- data.frame(FID=Eastern_Asian_number, IID=Eastern_Asian_number)
write.table(Eastern_Asian_number, "220_phenotype/reference/Ancestry_info/Eastern_Asian_ancestry.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(EAS_ids, "220_phenotype/reference/Ancestry_info/EAS_id.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")


European_info <- ancestry_info[ancestry_info$Superpopulation.code=="EUR" | grepl("European", ancestry_info$Superpopulation.name), ]
European_number <- European_info[,1]
# 670 in total

EUR_ids <- data.frame(FID = European_number, IID = European_number)
write.table(EUR_ids, "/n/holyscratch01/duan_lab/ch_zhu/220_phenotype/reference/Ancestry_info/EUR_id.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


#####################
###### Subset data for Ancestry Groups
#####################

# Construct files that only contain Eastern Asian subjects
# Define the path to the plink executable
plink_path <- "/n/home01/czhu/bin/plink"
# Define the directory path for input and output files
base_path <- "/220_phenotype/reference/bfile"
keep_file <- "/220_phenotype/reference/Ancestry_info/EAS_id.txt"

# Loop
for (chr_num in 1:22){
  input_file <- paste0(base_path, "/chr", chr_num)
  output_file <- paste0(base_path, "/EAS_chr", chr_num)
  
  cmd <- paste0(plink_path, " --bfile ", input_file, " --keep ", keep_file, " --make-bed", " --out ", output_file)
  
  system(cmd)
}

# Construct files that only contain European subjects
# Define the path to the plink executable
plink_path <- "/n/home01/czhu/bin/plink"
# Define the directory path for input and output files
base_path <- "/220_phenotype/reference/bfile"
keep_file_EUR <- "/220_phenotype/reference/Ancestry_info/EUR_id.txt"

# Loop
for (chr_num in 1:22){
  input_file <- paste0(base_path, "/chr", chr_num)
  output_file <- paste0(base_path, "/EUR_chr", chr_num)
  
  cmd <- paste0(plink_path, " --bfile ", input_file, " --keep ", keep_file_EUR, " --make-bed", " --out ", output_file)
  
  system(cmd)
}

####################
##### Minor Allele Frequency (MAF) Computation
####################

# Computation of MAF
for (chr_num in 1:22){
  input_file <- paste0(base_path, "/EAS_chr", chr_num)
  output_file <- paste0(base_path, "/EAS_chr", chr_num)
  
  cmd <- paste0(plink_path, " --bfile ", input_file, " --freq", " --out ", output_file)
  system(cmd)
}


for (chr_num in 1:22){
  input_file <- paste0(base_path, "/EUR_chr", chr_num)
  output_file <- paste0(base_path, "/EUR_chr", chr_num)
  
  cmd <- paste0(plink_path, " --bfile ", input_file, " --freq", " --out ", output_file)
  system(cmd)
}

count_low_maf_snps <- vector("list", 22)
count_snps <- vector("list", 22)
prop_snps <- vector("list", 22)

high_maf_snps <- vector("list", 22)

count <- 0

# Filter SNPs based on MAF thresholds (>0.005)

for (chr_num in 1:22){
  file_path <- paste0(base_path, "/EAS_chr", chr_num, ".frq")
  frq_data <- data.table::fread(file_path)
  high_maf_snps[chr_num] <- frq_data[frq_data$MAF > 0.005, 2]
  count_low_maf_snps[chr_num] <- sum(frq_data$MAF < 0.005, na.rm=TRUE)
  count_snps[chr_num] <- nrow(frq_data)
  prop_snps[chr_num] <- sum(frq_data$MAF < 0.005, na.rm=TRUE)/nrow(frq_data)
  
  # Clean up memory for each iteration
  rm(frq_data)
  gc()
  
  count = count+1
  print(count)
}

write.table(prop_snps, "/220_phenotype/reference/bfile/EAS_low_maf_rate.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
saveRDS(high_maf_snps, "/220_phenotype/reference/bfile/EAS_high_maf_list.rds")


EUR_count_low_maf_snps <- vector("list", 22)
EUR_count_snps <- vector("list", 22)
EUR_prop_snps <- vector("list", 22)

EUR_high_maf_snps <- vector("list", 22)

count <- 0

for (chr_num in 1:22){
  file_path <- paste0(base_path, "/EUR_chr", chr_num, ".frq")
  frq_data <- data.table::fread(file_path)
  EUR_high_maf_snps[chr_num] <- frq_data[frq_data$MAF > 0.005, 2]
  EUR_count_low_maf_snps[chr_num] <- sum(frq_data$MAF < 0.005, na.rm=TRUE)
  EUR_count_snps[chr_num] <- nrow(frq_data)
  EUR_prop_snps[chr_num] <- sum(frq_data$MAF < 0.005, na.rm=TRUE)/nrow(frq_data)
  
  # Clean up memory for each iteration
  rm(frq_data)
  gc()
  
  count = count+1
  print(count)
}

write.table(EUR_prop_snps, "/220_phenotype/reference/bfile/EUR_low_maf_rate.txt", quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
saveRDS(EUR_high_maf_snps, "/220_phenotype/reference/bfile/EUR_high_maf_list.rds")


##################
##### identify SNPs shared between EAS and EUR population
#################

# Extracted snps of EUR dataset
EUR_beta <- readRDS("220_phenotype_download/processed/SNP_auto_beta_e3.RDS")
EUR_POS <- readRDS("/220_phenotype/EUR_auto_screen_POS_e3.RDS")
EUR_POS_all <- unlist(EUR_POS)
EUR_POS_unique <- unique(EUR_POS_all)


# Extracted snps of UKB dataset
UKB_beta <- readRDS("220_phenotype_download/processed/UKB_auto_beta_e3.RDS")
UKB_POS <- readRDS("/220_phenotype/UKB_auto_screen_POS_e3.RDS")
UKB_POS_all <- unlist(UKB_POS)
UKB_POS_unique <- unique(UKB_POS_all)

EUR_meta_POS_unique <- c(UKB_POS_unique, EUR_POS_unique)
EUR_meta_POS_unique <- unique(EUR_meta_POS_unique)


# Extracted snps of BBJ dataset
BBJ_beta <- readRDS("220_phenotype_download/processed/BBJ_auto_beta_e3.RDS")
BBJ_POS <- readRDS("/220_phenotype/BBJ_auto_screen_POS_e3.RDS")
BBJ_POS_all <- unlist(BBJ_POS)
BBJ_POS_unique <- unique(BBJ_POS_all)

sig_POS_all <- intersect(EUR_meta_POS_unique, BBJ_POS_unique)


EAS_rate <- read.csv("/220_phenotype/reference/bfile/EAS_low_maf_rate.txt", sep="\t", header=F)
EUR_rate <- read.csv("/220_phenotype/reference/bfile/EUR_low_maf_rate.txt", sep="\t", header=F)
EUR_high_snps <- readRDS("/220_phenotype/reference/bfile/EUR_high_maf_list.rds")
EAS_high_snps <- readRDS("/220_phenotype/reference/bfile/EAS_high_maf_list.rds")

# 11053424 SNPs in EUR population
EUR_high_snps_unlist <- unlist(EUR_high_snps)
# 9774425 SNPs in EAS population
EAS_high_snps_unlist <- unlist(EAS_high_snps)
# 7348804 snps that pass both the thresholds
high_snps_list <- intersect(EUR_high_snps_unlist, EAS_high_snps_unlist)
saveRDS(high_snps_list,"/220_phenotype/reference/bfile/high_maf_all_list.rds")

sig_POS_avail_all <- intersect(high_snps_list, sig_POS_all)
saveRDS(sig_POS_avail_all," /220_phenotype/reference/bfile/avail_sig_POS.rds")



####################
##### Linkage Disequilibrium (LD) Pruning
####################

# conduct LD pruning

library(data.table)
library(genio)

sig_POS_avail_all <- readRDS(file="/220_phenotype/reference/bfile/avail_sig_POS.rds")

vec.POS <- c()
vec.intersect.POS <- c()

# Loop through all 22 chromosomes
for (i in 1:22) {
  
  # Filter SNPs based on prefix
  SNP_auto_char_list <- sig_POS_avail_all[grepl(paste0("^", i, ":"), sig_POS_avail_all)]
  
  # Read the BIM file into a data.table
  bim_data <- fread(paste0("/220_phenotype/reference/bfile/chr", i, ".bim"), select = c(4))
  bim_data <- as.vector(bim_data$V4)
  
  # Extract positions from the filtered SNPs
  SNP_char_POS <- sub(paste0("^", i, ":"), "", SNP_auto_char_list)
  vec.POS[i] <- length(SNP_char_POS)
  
  # Intersect positions
  SNP_char_POS_intersect <- SNP_char_POS[SNP_char_POS %in% bim_data]
  vec.intersect.POS[i] <- length(SNP_char_POS_intersect)
  
  # Add "i:" in front of every position
  SNP_char_POS_intersect_complete <- paste0(i, ":", SNP_char_POS_intersect)
  
  # Write to file
  write.table(SNP_char_POS_intersect_complete, paste0("/220_phenotype/reference/filtered/22_chromosomes/SNP_new_char",i, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Make bfile that only contains the selected SNPs
  system(paste0("/n/home01/czhu/bin/plink --bfile /220_phenotype/reference/bfile/chr", i, " --extract /220_phenotype/reference/filtered/22_chromosomes/SNP_new_char", i, ".txt --make-bed --out /220_phenotype/reference/bfile/chr", i, "_1000_geno_filtered_new"))
  
  # Run PLINK command for pruning
  system(paste0("/n/home01/czhu/bin/plink --bfile /220_phenotype/reference/bfile/chr", i, "_1000_geno_filtered_new --indep-pairwise 50 5 0.1 --out /220_phenotype/reference/filtered/22_chromosomes/chr", i, "_1000_geno_pruned_new"))
}


pruned_snps_all <- list()

for (i in 1:22)
{
  # Read the pruned snps
  pruned_snps <- read.table(paste0("/220_phenotype/reference/filtered/22_chromosomes/chr", i, "_1000_geno_pruned_new.prune.in"))
  pruned_snps <- as.vector(pruned_snps$V1)
  pruned_snps_all[[i]] <- pruned_snps
}

vector_lengths <- lengths(pruned_snps_all)

# 27696 SNPs in total after pruning
total_elements <- sum(vector_lengths)

pruned_snps_unlist <- unlist(pruned_snps_all)
write.csv(pruned_snps_unlist, quote=FALSE, row.names=FALSE, file="/220_phenotype/reference/filtered/22_chromosomes/new_pruned_snps_all.txt")

pruned_snps_unlist <-  read.csv("/220_phenotype/reference/filtered/22_chromosomes/new_pruned_snps_all.txt")
pruned_snps_unlist <- as.vector(pruned_snps_unlist$x)


# Get the MAF of the selected SNPs, for later modification use

maf_EAS <- vector("list", 22)
count = 0

# Define the path to the plink executable
plink_path <- "/n/home01/czhu/bin/plink"

# Define the directory path for input and output files
base_path <- "/220_phenotype/reference/bfile"

for (chr_num in 1:22){
  file_path <- paste0(base_path, "/EAS_chr", chr_num, ".frq")
  frq_data <- data.table::fread(file_path)
  maf_EAS[[chr_num]] <- frq_data[frq_data$SNP %in% pruned_snps_unlist, c(2,5)]
  
  # Clean up memory for each iteration
  rm(frq_data)
  gc()
  
  count = count+1
  print(count)
}

maf_EAS_unlist <- do.call(rbind, maf_EAS)
saveRDS(maf_EAS_unlist, "/220_phenotype/reference/Ancestry_info/EAS_filtered_maf.rds")

maf_EUR <- vector("list", 22)
count = 0

# Define the path to the plink executable
plink_path <- "/n/home01/czhu/bin/plink"

# Define the directory path for input and output files
base_path <- "/220_phenotype/reference/bfile"

for (chr_num in 1:22){
  file_path <- paste0(base_path, "/EUR_chr", chr_num, ".frq")
  frq_data <- data.table::fread(file_path)
  maf_EUR[[chr_num]] <- frq_data[frq_data$SNP %in% pruned_snps_unlist, c(2,5)]
  
  # Clean up memory for each iteration
  rm(frq_data)
  gc()
  
  count = count+1
  print(count)
}

maf_EUR_unlist <- do.call(rbind, maf_EUR)
saveRDS(maf_EUR_unlist, "/220_phenotype/reference/Ancestry_info/EUR_filtered_maf.rds")


# 146 phenotypes
BBJ_beta <- readRDS("220_phenotype_download/processed/BBJ_auto_beta_e3.RDS")
# 118 phenotypes
EUR_beta <- readRDS("220_phenotype_download/processed/EUR_auto_beta_e3.RDS")
# 28 phenotypes
UKB_beta <- readRDS("220_phenotype_download/processed/UKB_auto_beta_e3.RDS")


# Specify the main directory path
main_directory <- "220_phenotype_download/unzip"

# Get a list of sub-directories within the main directory
subdirectories <- list.dirs(main_directory, full.names = TRUE, recursive = FALSE)

# Eliminate the hampered file
subdirectories <- subdirectories[-106]

# Initialize an empty list to store the BETA columns
file_name <- vector("character", 147)

count=1

# Loop through each sub-directory
for (subdir in subdirectories) {
  if (length(grep("BBJ",subdir)) > 0){
    # Get a list of all files with .gz extension within the subdirectory
    file_list <- list.files(subdir, pattern = "\\.gz$", full.names = TRUE)
    
    # Loop through each .gz file in the sub-directory
    for (file in file_list) {
      # Check if the file name contains "chrX"
      if (length(grep("auto",file)) > 0) {
        
        file_name[[count]] <- file
        count = count + 1
        print(count)
      }
    }
  }
}

BBJ_results <- gsub(".*\\.BBJ\\.(.*?)\\.v1.*", "\\1", file_name)
BBJ_results <- BBJ_results[-147]
BBJ_beta <- BBJ_beta[-147]
names(BBJ_beta) <- BBJ_results

saveRDS(BBJ_beta, "220_phenotype_download/processed/BBJ_auto_beta_withname_e3.RDS")


test <- readRDS("220_phenotype_download/processed/BBJ_auto_beta_withname_e3.RDS")


file.list <- c()
count = 1 

# Loop through each sub-directory
for (subdir in subdirectories) {
  if (length(grep("EUR",subdir)) > 0){
    # Get a list of all files with .gz extension within the subdirectory
    file_list <- list.files(subdir, pattern = "\\.gz$", full.names = TRUE)
    for (file in file_list) {
      # Check if the file name contains "chrX"
      if (length(grep("auto",file)) > 0) {
        file.list <- append(file.list, file)
        count = count+1
        print(count)
      }
    }
  }
}

# Use grep with invert = TRUE to find entries without "UKB" in the vector
non_matching_entries <- grep("UKB", file.list, value = TRUE, invert = TRUE)
EUR_results <- gsub(".*\\.EUR\\.(.*?)\\.v1.*", "\\1", non_matching_entries)
names(EUR_beta) <- EUR_results

saveRDS(EUR_beta, "220_phenotype_download/processed/EUR_auto_beta_withname_e3.RDS")

matching_entries <- grep("UKB", file.list, value = TRUE, invert = FALSE)
UKB_results <- gsub(".*\\.EUR\\.(.*?)\\.v1.*", "\\1", matching_entries)
names(UKB_beta) <- UKB_results

saveRDS(UKB_beta, "220_phenotype_download/processed/UKB_auto_beta_withname_e3.RDS")

library(dplyr)
filtered_BBJ_list <- lapply(BBJ_beta, function(df){
  df %>% filter(CHR_POS %in% pruned_snps_unlist)
})


# Each contains all 20k+ snps
# 146 phenotypes
BBJ_beta_filtered <- readRDS("220_phenotype_download/processed/BBJ_new_filtered_beta.RDS")
# 118 phenotypes
EUR_beta_filtered <- readRDS("220_phenotype_download/processed/EUR_new_filtered_beta.RDS")
# 28 phenotypes
UKB_beta_filtered <- readRDS("220_phenotype_download/processed/UKB_new_filtered_beta.RDS")



names(EUR_beta_filtered) <- EUR_results
names(UKB_beta_filtered) <- UKB_results

saveRDS(EUR_beta_filtered, "220_phenotype_download/processed/EUR_new_filtered_beta.RDS") 
saveRDS(UKB_beta_filtered, "220_phenotype_download/processed/UKB_new_filtered_beta.RDS")

EUR_beta_filtered <- readRDS("220_phenotype_download/processed/EUR_new_filtered_beta.RDS")
UKB_beta_filtered <- readRDS("220_phenotype_download/processed/UKB_new_filtered_beta.RDS")
BBJ_beta_filtered <- readRDS("220_phenotype_download/processed/BBJ_auto_beta_withname_e3.RDS")

BBJ_pheno <- names(BBJ_beta_filtered)
EUR_pheno <- names(EUR_beta_filtered)
UKB_pheno <- names(UKB_beta_filtered)
EUR_combine_pheno <- c(EUR_pheno, UKB_pheno)

intersect_pheno <- intersect(BBJ_pheno, EUR_combine_pheno)
setdiff(BBJ_pheno, EUR_combine_pheno)
setdiff(EUR_combine_pheno, BBJ_pheno)


# remove the duplicate SNPs in each list
library(dplyr)
remove_duplicates <- function(df){
  df <- as.data.frame(df)
  df_filtered <- df[!duplicated(df$CHR_POS),]
  return(df_filtered)
}

# 27696 SNPs
cleaned_EUR_beta_filtered <- lapply(EUR_beta_filtered, remove_duplicates)
saveRDS(cleaned_EUR_beta_filtered, "220_phenotype_download/processed/EUR_new_filtered_beta_cleaned.RDS") 

# 27124 SNPs
cleaned_UKB_beta_filtered <- lapply(UKB_beta_filtered, remove_duplicates)
saveRDS(cleaned_UKB_beta_filtered, "220_phenotype_download/processed/UKB_new_filtered_beta_cleaned.RDS") 


cleaned_EUR_beta_filtered <- readRDS("220_phenotype_download/processed/EUR_new_filtered_beta_cleaned.RDS")
cleaned_UKB_beta_filtered <- readRDS("220_phenotype_download/processed/UKB_new_filtered_beta_cleaned.RDS")


#####################
##### Beta Coefficient modification
#####################

# Combine the data frames from both lists into one list
combined_list <- c(cleaned_EUR_beta_filtered, cleaned_UKB_beta_filtered)

# Create a list of unique IDs
unique_ids <- unique(unlist(lapply(combined_list, function(df) df$CHR_POS)))

# Create an empty data frame with columns as unique IDs
merged_df <- data.frame(matrix(0, nrow = length(combined_list), ncol = length(unique_ids)))
colnames(merged_df) <- unique_ids

# Fill in the merged data frame with coefficients
for (i in seq_along(combined_list)) {
  df <- combined_list[[i]]
  id_indices <- match(df$CHR_POS, unique_ids)
  merged_df[i, id_indices] <- df$BETA
}

# Set row names to be the names of the original data frames
rownames(merged_df) <- names(combined_list)
saveRDS(merged_df, "220_phenotype_download/processed/EUR_population_original_matrix.RDS")

merged_df <- readRDS( "220_phenotype_download/processed/EUR_population_original_matrix.RDS")

# read the minor allele frequency data
EUR_filtered <- readRDS("220_phenotype/reference/Ancestry_info/EUR_filtered_maf.rds")

EUR_frq <- EUR_filtered$MAF
EUR_sd <- sqrt(EUR_frq * (1-EUR_frq))

merged_df_EUR_maf_modified <- as.matrix(merged_df) %*% diag(EUR_sd)
merged_df_EUR_maf_modified <- as.data.frame(merged_df_EUR_maf_modified)
colnames(merged_df_EUR_maf_modified) <- colnames(merged_df)
saveRDS(merged_df_EUR_maf_modified, "220_phenotype_download/processed/EUR_population_maf_modified_matrix.RDS")

merged_df_EUR_maf_modified <- readRDS("220_phenotype_download/processed/EUR_population_maf_modified_matrix.RDS")


# 146 phenotypes
BBJ_beta_filtered <- readRDS("220_phenotype_download/processed/BBJ_new_filtered_beta.RDS")
BBJ_beta_filtered <- BBJ_beta_filtered[-c(147:150)]

names(BBJ_beta_filtered) <- BBJ_results
saveRDS(BBJ_beta_filtered,"220_phenotype_download/processed/BBJ_new_filtered_beta.RDS")

cleaned_BBJ_beta_filtered <- lapply(BBJ_beta_filtered, remove_duplicates)
saveRDS(cleaned_BBJ_beta_filtered, "220_phenotype_download/processed/BBJ_new_filtered_beta_cleaned.RDS") 



# Create a list of unique IDs
BBJ_unique_ids <- unique(unlist(lapply(cleaned_BBJ_beta_filtered, function(df) df$CHR_POS)))

# Create an empty data frame with columns as unique IDs
BBJ_merged_df <- data.frame(matrix(0, nrow = length(cleaned_BBJ_beta_filtered), ncol = length(BBJ_unique_ids)))
colnames(BBJ_merged_df) <- BBJ_unique_ids

# Fill in the merged data frame with coefficients
for (i in seq_along(cleaned_BBJ_beta_filtered)) {
  df <- cleaned_BBJ_beta_filtered[[i]]
  id_indices <- match(df$CHR_POS, BBJ_unique_ids)
  BBJ_merged_df[i, id_indices] <- df$beta_column
}

# Set row names to be the names of the original data frames
rownames(BBJ_merged_df) <- names(cleaned_BBJ_beta_filtered)
saveRDS(BBJ_merged_df, "220_phenotype_download/processed/BBJ_population_original_matrix.RDS")

# read the minor allele frequency data
EAS_filtered <- readRDS("220_phenotype/reference/Ancestry_info/EAS_filtered_maf.rds")

EAS_frq <- EAS_filtered$MAF
EAS_sd <- sqrt(EAS_frq * (1-EAS_frq))

merged_df_BBJ_maf_modified <- matrix(as.numeric(unlist(BBJ_merged_df)), nrow=nrow(BBJ_merged_df)) %*% diag(EAS_sd)
merged_df_BBJ_maf_modified <- as.data.frame(merged_df_BBJ_maf_modified)
colnames(merged_df_BBJ_maf_modified) <- colnames(BBJ_merged_df)
rownames(merged_df_BBJ_maf_modified) <- rownames(BBJ_merged_df)

saveRDS(merged_df_BBJ_maf_modified, "220_phenotype_download/processed/BBJ_population_maf_modified_matrix.RDS")
