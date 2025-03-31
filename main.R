# Project: Microbiome Contamination Analysis - Reorganized
# Description: Analysis of contaminants in microbiome sequencing data
# Author: QIANG DONG (Original code)
# Last Updated: 2025-03-31

#==============================================================================
# Package Management and Environment Setup
#==============================================================================
## Required packages
required_packages <- list(
  cran = c("reshape2", "ggplot2", "vegan", "stringr", "gridExtra", "ape",
           "RColorBrewer", "dplyr", "knitr", "cowplot", "openxlsx", 
           "circlize", "plotly", "here", "openxlsx", "UpSetR", "reshape2", 
           "ggrepel", "httr", "R.utils"),
  bioc = c("phyloseq", "decontam", "ComplexHeatmap", "dada2", "pheatmap")
)

## Install and load packages
install_and_load <- function(packages, type = "cran") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (type == "bioc") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
    library(pkg, character.only = TRUE)
    message(sprintf("%s version: %s", pkg, packageVersion(pkg)))
  }
}

## Initialize environment
install_and_load(required_packages$cran)
install_and_load(required_packages$bioc, "bioc")

## Project directory structure - Reorganized
project_dirs <- list(
  root = here(),
  data = here("data"),
  data_raw = here("data/raw"),
  data_processed = here("data/processed"),
  results = here("results"),
  figures = here("results/Figures"),
  tables = here("results/Tables"),
  ref_db = here("ref_databases")
)

## Create directories
sapply(project_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

## Color schemes
color_schemes <- list(
  main = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
           "#0072B2", "#D55E00", "#CC79A7", "#000000"),
  methods = c(Frequency = "darkorange2", 
              Prevalence = "turquoise3", 
              Combined = "darkorchid4"),
  prevalence = c("2" = "#edf8e9", "3-5" = "#bae4b3", 
                 "6-10" = "#74c476", "11+" = "#238b45")
)

#==============================================================================
# 0. Reference Database Management
#==============================================================================
## Define reference database files with direct Dropbox links
ref_db_links <- list(
  "silva_nr99_v138.2_toSpecies_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/erhd1gcwuncxkkd5nxeku/silva_nr99_v138.2_toSpecies_trainset.fa.gz?rlkey=ilo06ftnjl9we6pfcgaypy5a1&st=bcywl890&dl=1",
  "silva_nr99_v138.2_toGenus_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/lnwwdbqm5eyqgj8xktw3f/silva_nr99_v138.2_toGenus_trainset.fa.gz?rlkey=tjwi2bb4lkejemvie93dsxeu9&st=slcdc9h3&dl=1",
  "rdp_19_toSpecies_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/7cfeowmxfw19ftouomzth/rdp_19_toSpecies_trainset.fa.gz?rlkey=sjex7itpw766df6gn82t23f3q&st=9n0sp2ye&dl=1",
  "rdp_19_toGenus_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/uvr2u3nwzxunt72cfnjyl/rdp_19_toGenus_trainset.fa.gz?rlkey=kdlmuuxcu3hdyhmcoe0l2rbfz&st=u8m032tg&dl=1",
  "gg2_2024_09_toSpecies_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/zcd24frb4idegpe3r0fuv/gg2_2024_09_toSpecies_trainset.fa.gz?rlkey=7fxzvqmoego7v24051zgouhbv&st=2svmowvv&dl=1",
  "gg2_2024_09_toGenus_trainset.fa.gz" = "https://www.dropbox.com/scl/fi/qc82i4q8ocfjapjoytw6f/gg2_2024_09_toGenus_trainset.fa.gz?rlkey=0j0lhzf5wnarv8xq2lb1wkmyf&st=lghz38r0&dl=1"
)

## Function to download reference databases from individual Dropbox links
download_reference_dbs <- function() {
  message("Starting download of reference databases from Dropbox...")
  
  # Track downloaded files
  downloaded_files <- list()
  
  for (file_name in names(ref_db_links)) {
    dest_path <- file.path(project_dirs$ref_db, file_name)
    
    # Check if file already exists
    if (file.exists(dest_path)) {
      message(sprintf("File already exists: %s", file_name))
      downloaded_files[[file_name]] <- dest_path
      next
    }
    
    # Get direct download URL for this file (adding dl=1 ensures direct download)
    file_url <- ref_db_links[[file_name]]
    
    message(sprintf("Downloading %s...", file_name))
    
    # Use httr to download the file
    tryCatch({
      response <- httr::GET(file_url, httr::write_disk(dest_path, overwrite = TRUE),
                            httr::progress())
      
      # Check if download was successful
      if (httr::status_code(response) == 200) {
        message(sprintf("Successfully downloaded: %s", file_name))
        downloaded_files[[file_name]] <- dest_path
      } else {
        warning(sprintf("Failed to download %s. Status code: %d", 
                        file_name, httr::status_code(response)))
      }
    }, error = function(e) {
      warning(sprintf("Error downloading %s: %s", file_name, e$message))
    })
  }
  
  # Check which files were successfully downloaded
  if (length(downloaded_files) > 0) {
    message(sprintf("Successfully downloaded %d reference database files.", length(downloaded_files)))
  } else {
    warning("No reference database files were downloaded successfully.")
  }
  
  return(downloaded_files)
}

## Configure reference databases
configure_ref_dbs <- function(downloaded_files) {
  ref_dbs <- list()
  
  # Add Silva databases
  if ("silva_nr99_v138.2_toGenus_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$silva_genus <- list(
      name = "Silva Genus",
      version = "v138.2",
      path = downloaded_files[["silva_nr99_v138.2_toGenus_trainset.fa.gz"]],
      description = "Silva rRNA database for bacterial and archaeal taxonomy (genus level)"
    )
  }
  
  if ("silva_nr99_v138.2_toSpecies_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$silva_species <- list(
      name = "Silva Species",
      version = "v138.2",
      path = downloaded_files[["silva_nr99_v138.2_toSpecies_trainset.fa.gz"]],
      description = "Silva rRNA database for bacterial and archaeal taxonomy (species level)"
    )
  }
  
  # Add GreenGenes databases
  if ("gg2_2024_09_toGenus_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$gg_genus <- list(
      name = "GreenGenes Genus",
      version = "2024_09",
      path = downloaded_files[["gg2_2024_09_toGenus_trainset.fa.gz"]],
      description = "GreenGenes database for 16S rRNA taxonomy (genus level)"
    )
  }
  
  if ("gg2_2024_09_toSpecies_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$gg_species <- list(
      name = "GreenGenes Species",
      version = "2024_09",
      path = downloaded_files[["gg2_2024_09_toSpecies_trainset.fa.gz"]],
      description = "GreenGenes database for 16S rRNA taxonomy (species level)"
    )
  }
  
  # Add RDP databases
  if ("rdp_19_toGenus_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$rdp_genus <- list(
      name = "RDP Genus",
      version = "19",
      path = downloaded_files[["rdp_19_toGenus_trainset.fa.gz"]],
      description = "RDP (Ribosomal Database Project) for bacterial taxonomy (genus level)"
    )
  }
  
  if ("rdp_19_toSpecies_trainset.fa.gz" %in% names(downloaded_files)) {
    ref_dbs$rdp_species <- list(
      name = "RDP Species",
      version = "19",
      path = downloaded_files[["rdp_19_toSpecies_trainset.fa.gz"]],
      description = "RDP (Ribosomal Database Project) for bacterial taxonomy (species level)"
    )
  }
  
  # Save database configurations for later use
  saveRDS(ref_dbs, file.path(project_dirs$data_processed, "reference_db_configs.rds"))
  
  return(ref_dbs)
}

# Download and configure reference databases
downloaded_files <- download_reference_dbs()
reference_dbs <- configure_ref_dbs(downloaded_files)

#===============================================================================
# 1. Load Raw Data and Add Metadata
#===============================================================================
## Load the initial DADA2 phyloseq object
raw_phyloseq <- readRDS(file.path(project_dirs$data_raw, "dada2_phyloseq.rds"))

## Load and process metadata
metadata <- read.csv(file.path(project_dirs$data_raw, "metadata.csv"), row.names = 1)
sample_data(raw_phyloseq) <- sample_data(metadata)

## Save phyloseq object with new metadata
total <- raw_phyloseq  # Rename for consistency with later code
saveRDS(total, file.path(project_dirs$data_processed, "DADA2_with_new_metadata.rds"))


#===============================================================================
# 2. Assign Taxonomy
#===============================================================================

## Get the tax table and create initial TAXannotate
tax_data <- tax_table(total)
TAXannotate <- data.frame(
  ASV_ID = rownames(tax_data),
  sequence = tax_data[, "sequence"],
  stringsAsFactors = FALSE
)
rownames(TAXannotate) <- TAXannotate$ASV_ID

## Set random seed
set.seed(100)

## Function to assign taxonomy
assign_taxonomy <- function(sequences, db_path, output_name) {
  # Check if sequences is a vector or single sequence
  if (!is.vector(sequences)) {
    sequences <- as.character(sequences)
  }
  
  # Remove any NA or empty sequences
  sequences <- sequences[!is.na(sequences) & nchar(sequences) > 0]
  
  # Proceed with taxonomy assignment
  tax <- tryCatch({
    assignTaxonomy(sequences, db_path, 
                   multithread = TRUE, 
                   taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  }, error = function(e) {
    print(paste("Error in taxonomy assignment:", e$message))
    return(NULL)
  })
  
  if (!is.null(tax)) {
    saveRDS(tax, file.path(here("results/Tables"), paste0(output_name, ".rds")))
    write.csv(tax, file.path(here("results/Tables"), paste0(output_name, ".csv")), row.names = TRUE)
  }
  
  return(tax)
}

## Assign taxonomy using different databases
t138 <- assign_taxonomy(TAXannotate$sequence,
                        file.path(here("ref_databases"), "silva_nr99_v138.2_toSpecies_trainset.fa.gz"),
                        "Table1A_silva_taxonomy")

rdp19 <- assign_taxonomy(TAXannotate$sequence,
                         file.path(here("ref_databases"), "rdp_19_toSpecies_trainset.fa.gz"),
                         "Table1B_rdp_taxonomy")

gg2 <- assign_taxonomy(TAXannotate$sequence,
                       file.path(here("ref_databases"), "gg2_2024_09_toSpecies_trainset.fa.gz"),
                       "Table1C_greengenes_taxonomy")


## Add all taxonomic levels assignments
TAXannotate <- cbind(TAXannotate, 
                     silva_kingdom = t138[, "Kingdom"],
                     silva_phylum = t138[, "Phylum"],
                     silva_class = t138[, "Class"],
                     silva_order = t138[, "Order"],
                     silva_family = t138[, "Family"],
                     silva_genus = t138[, "Genus"],
                     silva_species = t138[, "Species"],
                     rdp_kingdom = rdp19[, "Kingdom"],
                     rdp_phylum = rdp19[, "Phylum"],
                     rdp_class = rdp19[, "Class"],
                     rdp_order = rdp19[, "Order"],
                     rdp_family = rdp19[, "Family"],
                     rdp_genus = rdp19[, "Genus"],
                     rdp_species = rdp19[, "Species"],
                     # Modified gsub to handle NAs properly
                     gg_kingdom = ifelse(!is.na(gg2[, "Kingdom"]), gsub("d__", "", gg2[, "Kingdom"]), NA),
                     gg_phylum = ifelse(!is.na(gg2[, "Phylum"]), gsub("p__", "", gg2[, "Phylum"]), NA),
                     gg_class = ifelse(!is.na(gg2[, "Class"]), gsub("c__", "", gg2[, "Class"]), NA),
                     gg_order = ifelse(!is.na(gg2[, "Order"]), gsub("o__", "", gg2[, "Order"]), NA),
                     gg_family = ifelse(!is.na(gg2[, "Family"]), gsub("f__", "", gg2[, "Family"]), NA),
                     gg_genus = ifelse(!is.na(gg2[, "Genus"]), gsub("g__", "", gg2[, "Genus"]), NA),
                     gg_species = ifelse(!is.na(gg2[, "Species"]), gsub("s__", "", gg2[, "Species"]), NA))

# Clean up function to remove duplicate columns
clean_TAXannotate <- function(df) {
  unique_cols <- unique(colnames(df))
  df <- df[, match(unique_cols, colnames(df))]
  return(df)
}

# Clean up any duplicate columns
TAXannotate <- clean_TAXannotate(TAXannotate)

# Define all taxonomic levels
tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Function to count non-NA annotations in lower levels for a database
count_lower_annotations <- function(df, row_idx, current_level, database) {
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  current_idx <- which(tax_levels == current_level)
  if (current_idx == length(tax_levels)) return(0)
  
  lower_levels <- tax_levels[(current_idx + 1):length(tax_levels)]
  count <- 0
  
  for (level in lower_levels) {
    level_lower <- tolower(level)
    if (!is.na(df[row_idx, paste0(database, "_", level_lower)])) {
      count <- count + 1
    }
  }
  
  return(count)
}

# Function to determine best database based on lower level completeness
get_best_database <- function(df, row_idx, current_level, agreeing_dbs) {
  # Count lower level annotations for each database
  counts <- sapply(agreeing_dbs, function(db) {
    count_lower_annotations(df, row_idx, current_level, db)
  })
  
  # If counts are different, use the database with more annotations
  max_count <- max(counts)
  if (sum(counts == max_count) == 1) {
    return(agreeing_dbs[which.max(counts)])
  }
  
  # If counts are tied, use priority order
  if ("silva" %in% agreeing_dbs) return("silva")
  if ("rdp" %in% agreeing_dbs) return("rdp")
  return("gg")
}

# Create recommendation function
create_recommendation <- function(df, level) {
  level_lower <- tolower(level)
  cols <- c(
    paste0("silva_", level_lower),
    paste0("rdp_", level_lower),
    paste0("gg_", level_lower)
  )
  
  df[[paste0(level, "_recommendation")]] <- case_when(
    (!is.na(df[[cols[1]]]) & !is.na(df[[cols[2]]]) & !is.na(df[[cols[3]]]) &
       df[[cols[1]]] == df[[cols[2]]] & df[[cols[2]]] == df[[cols[3]]]) ~ "3DB",
    (!is.na(df[[cols[1]]]) & !is.na(df[[cols[2]]]) & 
       df[[cols[1]]] == df[[cols[2]]]) ~ "2DB_SR",
    (!is.na(df[[cols[1]]]) & !is.na(df[[cols[3]]]) & 
       df[[cols[1]]] == df[[cols[3]]]) ~ "2DB_SG",
    (!is.na(df[[cols[2]]]) & !is.na(df[[cols[3]]]) & 
       df[[cols[2]]] == df[[cols[3]]]) ~ "2DB_RG",
    !is.na(df[[cols[1]]]) ~ "1DB_S",
    !is.na(df[[cols[2]]]) ~ "1DB_R",
    !is.na(df[[cols[3]]]) ~ "1DB_G",
    TRUE ~ "0DB"
  )
  
  return(df)
}

# Modified assign_final_recommendation function
assign_final_recommendation <- function(df, level, database_assignments) {
  level_lower <- tolower(level)
  rec_col <- paste0(level, "REC")
  df[[rec_col]] <- NA_character_
  
  recommendation_col <- paste0(level, "_recommendation")
  cols <- c(
    paste0("silva_", level_lower),
    paste0("rdp_", level_lower),
    paste0("gg_", level_lower)
  )
  
  # For each row
  for (i in 1:nrow(df)) {
    current_rec <- df[i, recommendation_col]
    prev_db <- database_assignments[i]
    
    # Determine which databases to consider
    if (!is.na(prev_db)) {
      # Use the same database as before if it has a value
      db_idx <- match(prev_db, c("silva", "rdp", "gg"))
      if (!is.na(df[i, cols[db_idx]])) {
        df[i, rec_col] <- df[i, cols[db_idx]]
        next
      }
    }
    
    # Determine which databases agree at this level
    agreeing_dbs <- switch(current_rec,
                           "3DB" = c("silva", "rdp", "gg"),
                           "2DB_SR" = c("silva", "rdp"),
                           "2DB_SG" = c("silva", "gg"),
                           "2DB_RG" = c("rdp", "gg"),
                           "1DB_S" = "silva",
                           "1DB_R" = "rdp",
                           "1DB_G" = "gg",
                           character(0)
    )
    
    if (length(agreeing_dbs) > 0) {
      if (length(agreeing_dbs) == 2) {
        # For 2DB conditions, look at lower level completeness
        chosen_db <- get_best_database(df, i, level, agreeing_dbs)
        database_assignments[i] <- chosen_db
        db_idx <- match(chosen_db, c("silva", "rdp", "gg"))
        df[i, rec_col] <- df[i, cols[db_idx]]
      } else if (length(agreeing_dbs) == 1) {
        # For 1DB conditions, use that database
        database_assignments[i] <- agreeing_dbs
        db_idx <- match(agreeing_dbs, c("silva", "rdp", "gg"))
        df[i, rec_col] <- df[i, cols[db_idx]]
      } else {
        # For 3DB conditions, prioritize silva > rdp > gg
        for (db in c("silva", "rdp", "gg")) {
          db_idx <- match(db, c("silva", "rdp", "gg"))
          if (!is.na(df[i, cols[db_idx]])) {
            df[i, rec_col] <- df[i, cols[db_idx]]
            database_assignments[i] <- db
            break
          }
        }
      }
    }
  }
  
  return(list(df = df, database_assignments = database_assignments))
}

# Initialize database assignments
database_assignments <- rep(NA_character_, nrow(TAXannotate))

# Process all taxonomic levels
for (level in tax_levels) {
  # Create recommendations
  TAXannotate <- create_recommendation(TAXannotate, level)
  
  # Assign taxonomy with the new strategy
  result <- assign_final_recommendation(TAXannotate, level, database_assignments)
  TAXannotate <- result$df
  database_assignments <- result$database_assignments
}

# Update the phyloseq object
tax_table(total) <- tax_table(as.matrix(TAXannotate))

## Save results
write.csv(TAXannotate, 
          file.path(here("results/Tables"), "Table1D_TAXannotate_full.csv"),
          row.names = FALSE)

saveRDS(total, here("data/processed/DADA2_with_new_metadata.rds"))

## Create subset based on Decontam analysis
ps <- subset_samples(total, Analysis_2 == "Decontam")
saveRDS(ps, here("data/processed/DADA2_with_new_metadata_subset.rds"))


#===============================================================================
# 4. Decontamination analysis
#===============================================================================

# Create necessary directories
dir.create(here("results/Tables/contamination_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results/Figures/contamination_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results/Tables/extracted_asvs"), recursive = TRUE, showWarnings = FALSE)

# Calculate per-sample contamination metrics
calculate_per_sample_metrics <- function(ps_obj, tax2_data = NULL) {
  # Extract OTU table
  otu_mat <- as(otu_table(ps_obj), 'matrix')
  if (!taxa_are_rows(ps_obj)) {
    otu_mat <- t(otu_mat)
  }
  
  # Get metadata
  meta <- data.frame(sample_data(ps_obj))
  neg_samples <- meta$Analysis_1 == "GF_control"
  
  # Get negative control data
  neg_data <- otu_mat[, neg_samples]
  
  # Calculate negative control frequencies
  neg_freqs <- sweep(neg_data, 2, colSums(neg_data), '/')
  neg_mean_freqs <- rowMeans(neg_freqs)
  neg_spans <- rowSums(neg_data > 0) / ncol(neg_data)
  
  # Get decontam values if TAX2 is provided
  if (!is.null(tax2_data)) {
    decontam_cols <- c("prob.f", "prob.p", "prob.c", "if_GF")
    tax2_data$ASV_ID <- rownames(tax2_data)
    decontam_values <- tax2_data[, c("ASV_ID", intersect(decontam_cols, colnames(tax2_data)))]
  }
  
  # Process each non-negative sample
  sample_results <- list()
  for(sample in colnames(otu_mat)[!neg_samples]) {
    # Get sample data and convert to frequencies
    sample_counts <- otu_mat[, sample]
    sample_freqs <- sample_counts / sum(sample_counts)
    
    # Filter for non-zero ASVs in this sample
    non_zero_idx <- which(sample_counts > 0)
    
    # Skip if no non-zero ASVs
    if(length(non_zero_idx) == 0) {
      warning(sprintf("Sample %s has no non-zero ASVs, skipping.", sample))
      next
    }
    
    # Filter sample data for non-zero ASVs
    sample_counts_filtered <- sample_counts[non_zero_idx]
    
    # Calculate frequencies for filtered data
    sample_freqs <- sample_counts_filtered / sum(sample_counts_filtered)
    
    # Get corresponding negative control data for these ASVs
    filtered_neg_mean_freqs <- neg_mean_freqs[non_zero_idx]
    filtered_neg_spans <- neg_spans[non_zero_idx]
    
    # Calculate metrics for this sample
    results <- data.frame(
      ASV_ID = names(sample_counts_filtered),
      Sample = sample,
      sample_freq = sample_freqs,
      neg_mean_freq = filtered_neg_mean_freqs,
      neg_span = filtered_neg_spans,
      sample_reads = sample_counts_filtered,
      neg_mean_reads = rowMeans(neg_data[non_zero_idx, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
    
    # Calculate frequency ratio
    results$freq_ratio <- results$neg_mean_freq / results$sample_freq
    results$freq_ratio[is.infinite(results$freq_ratio)] <- max(results$freq_ratio[is.finite(results$freq_ratio)])
    results$freq_ratio[is.nan(results$freq_ratio)] <- 0
    
    # Add decontam values if available
    if (!is.null(tax2_data)) {
      results <- merge(results, decontam_values, by = "ASV_ID", all.x = TRUE)
    }
    
    # Add taxonomy if available
    if (!is.null(tax_table(ps_obj))) {
      tax_df <- data.frame(tax_table(ps_obj))
      tax_cols <- setdiff(colnames(tax_df), c("prob.f", "prob.p", "prob.c", "if_standard"))
      results <- cbind(
        results,
        tax_df[match(results$ASV_ID, rownames(tax_df)), tax_cols, drop = FALSE]
      )
    }
    
    # Sort by frequency ratio and neg_span
    results <- results[order(results$freq_ratio, results$neg_span, decreasing = TRUE), ]
    
    # Store results for this sample
    sample_results[[sample]] <- results
  }
  
  return(sample_results)
}

# Export results to Excel
export_per_sample_metrics <- function(ps_obj, tax2_data = NULL) {
  # Calculate metrics
  sample_results <- calculate_per_sample_metrics(ps_obj, tax2_data)
  
  # Create Excel workbook
  wb <- createWorkbook()
  
  # Add overview sheet
  addWorksheet(wb, "Overview")
  overview <- data.frame(
    Sample = names(sample_results),
    Total_ASVs = sapply(sample_results, nrow),
    ASVs_in_50pct_neg = sapply(sample_results, function(x) sum(x$neg_span > 0.5)),
    Mean_freq_ratio = sapply(sample_results, function(x) mean(x$freq_ratio, na.rm = TRUE)),
    Possible_contaminants = sapply(sample_results, function(x) sum(x$freq_ratio > 1 & x$neg_span > 0.1)),
    stringsAsFactors = FALSE
  )
  writeData(wb, "Overview", overview)
  
  # Add individual sample sheets
  for(sample_name in names(sample_results)) {
    # Clean sample name for Excel sheet
    sheet_name <- make.names(sample_name)
    sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name length limit
    
    # Add worksheet
    addWorksheet(wb, sheet_name)
    
    # Write data
    writeData(wb, sheet_name, sample_results[[sample_name]])
    
    # Add autofilter and freeze pane
    addFilter(wb, sheet_name, row = 1, cols = 1:ncol(sample_results[[sample_name]]))
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
  
  # Save workbook
  saveWorkbook(wb, here("results/Tables/Table2A_per_sample_contamination_analysis.xlsx"), overwrite = TRUE)
  
  return(sample_results)
}

# Execute per-sample metrics export
export_per_sample_metrics(ps)

# Extract ASVs with specific cutoffs
# Define cutoffs at the beginning
freq_ratio_cutoff <- 0.297
sample_freq_cutoff <- 0.002

extract_asvs <- function() {
  library(openxlsx)
  
  # Read Excel file
  xlsx_path <- here("results/Tables/Table2A_per_sample_contamination_analysis.xlsx")
  wb <- loadWorkbook(xlsx_path)
  sheet_names <- sheets(wb)[sheets(wb) != "Overview"]
  
  # Initialize list to store ASV IDs for each sample
  asv_list <- list()
  
  # Initialize summary dataframe for Excel
  summary_df <- data.frame(
    Sample = character(),
    Total_ASVs = numeric(),
    GF_ASVs = numeric(),
    NonGF_ASVs = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Create a new workbook for output
  wb_out <- createWorkbook()
  
  # Add overview sheet
  addWorksheet(wb_out, "Overview")
  
  # Process each sample sheet
  for (sheet in sheet_names) {
    data <- read.xlsx(wb, sheet = sheet)
    
    # Apply cutoffs to extract ASVs
    extracted_asvs <- data[data$freq_ratio < freq_ratio_cutoff & 
                             data$sample_freq > sample_freq_cutoff, ]
    
    # Store ASV IDs for this sample
    asv_list[[sheet]] <- list(
      all_asvs = extracted_asvs$ASV_ID,
      gf_asvs = extracted_asvs$ASV_ID[extracted_asvs$if_GF == "yes"],
      non_gf_asvs = extracted_asvs$ASV_ID[extracted_asvs$if_GF == "no"]
    )
    
    # Add sample data to Excel workbook
    addWorksheet(wb_out, sheet)
    writeData(wb_out, sheet, extracted_asvs)
    
    # Add to summary dataframe
    summary_df <- rbind(summary_df, data.frame(
      Sample = sheet,
      Total_ASVs = length(asv_list[[sheet]]$all_asvs),
      GF_ASVs = length(asv_list[[sheet]]$gf_asvs),
      NonGF_ASVs = length(asv_list[[sheet]]$non_gf_asvs),
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate total unique ASVs across all samples
  all_unique_asvs <- unique(unlist(lapply(asv_list, function(x) x$all_asvs)))
  all_unique_gf <- unique(unlist(lapply(asv_list, function(x) x$gf_asvs)))
  all_unique_nongf <- unique(unlist(lapply(asv_list, function(x) x$non_gf_asvs)))
  
  # Add summary row to summary_df
  summary_df <- rbind(summary_df, data.frame(
    Sample = "Total_Unique",
    Total_ASVs = length(all_unique_asvs),
    GF_ASVs = length(all_unique_gf),
    NonGF_ASVs = length(all_unique_nongf),
    stringsAsFactors = FALSE
  ))
  
  # Add mean row
  summary_df <- rbind(summary_df, data.frame(
    Sample = "Mean_per_sample",
    Total_ASVs = mean(summary_df$Total_ASVs[summary_df$Sample != "Total_Unique"]),
    GF_ASVs = mean(summary_df$GF_ASVs[summary_df$Sample != "Total_Unique"]),
    NonGF_ASVs = mean(summary_df$NonGF_ASVs[summary_df$Sample != "Total_Unique"]),
    stringsAsFactors = FALSE
  ))
  
  # Write overview to Excel
  writeData(wb_out, "Overview", summary_df)
  
  # Save the Excel file
  excel_filename <- sprintf("Table2B_extracted_asvs_freqRatio%.3f_sampleFreq%.3f.xlsx", 
                            freq_ratio_cutoff, sample_freq_cutoff)
  saveWorkbook(wb_out, 
               here("results/Tables/extracted_asvs", excel_filename), 
               overwrite = TRUE)
  
  # Save individual CSV files
  dir.create(here("results/Tables/extracted_asvs"), showWarnings = FALSE)
  
  for (sample in names(asv_list)) {
    filename <- sprintf("Table2B_%s_freqRatio%.3f_sampleFreq%.3f.csv", 
                        sample, freq_ratio_cutoff, sample_freq_cutoff)
    
    # Create data frame for this sample
    sample_df <- data.frame(
      ASV_ID = asv_list[[sample]]$all_asvs,
      Type = ifelse(asv_list[[sample]]$all_asvs %in% asv_list[[sample]]$gf_asvs, 
                    "GF", "Non-GF")
    )
    
    # Save to CSV
    write.csv(sample_df,
              here("results/Tables/extracted_asvs", filename),
              row.names = FALSE)
  }
  
  return(list(
    asv_data = asv_list,
    summary = summary_df,
    totals = list(
      all_unique = all_unique_asvs,
      gf_unique = all_unique_gf,
      nongf_unique = all_unique_nongf
    )
  ))
}

# Execute the function to extract ASVs
extracted_asvs <- extract_asvs()

#===============================================================================
# 6. Reconstitute phyloseq project
#===============================================================================

process_and_export_data <- function(ps_obj, prefix, save_path = NULL) {
  ## Get full taxonomy information
  full_tax <- read.csv(here("results/Tables/Table1D_TAXannotate_full.csv"))
  
  ## Get reads count matrix and ensure it's samples as columns
  reads_count_matrix <- as.matrix(otu_table(ps_obj))
  if (!taxa_are_rows(ps_obj)) {
    reads_count_matrix <- t(reads_count_matrix)
  }
}

reconstitute_phyloseq <- function(original_phyloseq, extracted_asv_path, freq_ratio = 0.297, sample_freq = 0.002) {
  # Read Excel file and get sheets
  excel_filename <- sprintf("Table2B_extracted_asvs_freqRatio%.3f_sampleFreq%.3f.xlsx", freq_ratio, sample_freq)
  wb <- loadWorkbook(file.path(extracted_asv_path, excel_filename))
  sheet_names <- sheets(wb)[sheets(wb) != "Overview"]
  
  # Filter samples
  keep_samples <- sample_names(original_phyloseq)[sample_data(original_phyloseq)$Analysis_1 == "Sample"]
  filtered_ps <- prune_samples(keep_samples, original_phyloseq)
  
  # Get all ASVs from sheets
  asv_lists <- lapply(sheet_names, function(sheet) {
    read.xlsx(wb, sheet = sheet)$ASV_ID
  })
  all_asvs <- unique(unlist(asv_lists))
  
  # Create presence matrix
  presence_matrix <- matrix(0, nrow = length(all_asvs), ncol = length(sheet_names),
                            dimnames = list(all_asvs, sheet_names))
  for(sheet in sheet_names) {
    presence_matrix[, sheet] <- all_asvs %in% asv_lists[[which(sheet_names == sheet)]]
  }
  
  # Create OTU matrix
  orig_otu <- as(otu_table(filtered_ps), "matrix")
  new_otu <- matrix(0, nrow = length(all_asvs), ncol = ncol(orig_otu),
                    dimnames = list(all_asvs, colnames(orig_otu)))
  common_asvs <- intersect(rownames(orig_otu), all_asvs)
  new_otu[common_asvs, ] <- orig_otu[common_asvs, ]
  
  # Set zeros for non-extracted ASVs
  for(sample in sheet_names) {
    if(sample %in% colnames(new_otu)) {
      non_extracted <- all_asvs[!presence_matrix[, sample]]
      new_otu[non_extracted, sample] <- 0
    }
  }
  
  # Create new phyloseq object
  new_phyloseq <- phyloseq(
    otu_table(new_otu, taxa_are_rows = TRUE),
    tax_table(filtered_ps)[common_asvs, ],
    sample_data(filtered_ps)
  )
  
  # Return results
  list(
    absolute = new_phyloseq,
    relative = transform_sample_counts(new_phyloseq, function(x) x/sum(x)),
    extraction_summary = data.frame(
      ASV = all_asvs,
      Times_Extracted = rowSums(presence_matrix),
      Sample_List = apply(presence_matrix, 1, function(x) paste(sheet_names[x > 0], collapse = ", "))
    ),
    presence_matrix = presence_matrix
  )
}

# Reconstitute the phyloseq object with extracted ASVs
reconstituted_ps <- reconstitute_phyloseq(ps, "results/Tables/extracted_asvs")
recons_ps <- reconstituted_ps$absolute

# Export extraction summary
write.csv(reconstituted_ps$extraction_summary, 
          here("results/Tables/Table3A_extraction_summary.csv"), 
          row.names = FALSE)

# Export processed data
process_and_export_data(recons_ps, "reconstituted_ps", save_path = "results/Tables/reconstituted_analysis")



#===============================================================================
# Figure 2B: Beta Diversity Analysis
#===============================================================================

# Create distance matrix
vegdist <- phyloseq::distance(recons_ps, method = "bray")

# PCoA calculations
CmdScale <- cmdscale(vegdist, k = 10)
vars <- apply(CmdScale, 2, var)
percentVar <- round(100 * (vars/sum(vars)))

# Merge with metadata 
newResults <- merge(x = data.frame(CmdScale), 
                    y = data.frame(sample_data(recons_ps)), 
                    by = "row.names", all.x = TRUE)

colnames(newResults)[c(2,3)] <- c("PC1", "PC2")
colnames(newResults)[1] <- "name"

# Calculate centroids
centroids <- aggregate(cbind(PC1,PC2) ~ Group, data = newResults, mean)
newResults <- merge(newResults, centroids, by = "Group", suffixes = c("",".centroid"))

# Load required package first
library(ggrepel)

# Define the PDF output file
pdf(here("results/Figures/Figure2B_beta_diversity.pdf"), height = 5, width = 5)

# Create the PCoA plot
p_beta <- ggplot(newResults, aes(PC1, PC2, color = Group)) +
  geom_point(size = 5, alpha = 0.8) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("#1a80bb", "#ea801c")) +  
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = Group), size = 0) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold', size = 20),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_line(linetype = "dashed", linewidth = 0.5, colour = "grey80"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.title = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(colour = "grey80", size = rel(0.75)),
    axis.text.y = element_text(colour = "grey80", size = rel(0.75)),
    axis.ticks = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "line"),
    legend.position = "none"
  ) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2, color = Group), alpha = 0.2) +
  geom_label_repel(data = centroids, aes(x = PC1, y = PC2, label = Group), 
                   size = 6, fontface = "bold", point.padding = 0.5, box.padding = 0.5)

# Print the plot to PDF
print(p_beta)

# Close the PDF device
dev.off()

# Save the plot as a JPG file too
ggsave(here("results/Figures/Figure2B_beta_diversity.jpg"), p_beta, height = 5, width = 5)

# Run statistical test for beta diversity
adonis_result <- adonis2(vegdist ~ sample_data(recons_ps)$Group)

# Save statistical results
write.csv(as.data.frame(adonis_result), 
          here("results/Tables/Table3B_beta_diversity_stats.csv"), 
          row.names = TRUE)


#===============================================================================
# Figure 2C: Stacked Bar Plot and Presence/Absence Plot
#===============================================================================

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(tidyr)

# Transform counts to relative abundance
ps_rel = transform_sample_counts(recons_ps, function(x) x/sum(x))

# Get genus-level abundances
genus_abundance <- tax_glom(ps_rel, taxrank="GenusREC")
genus_abundance_df <- psmelt(genus_abundance)

# Calculate mean abundance for each genus
genus_means <- genus_abundance_df %>%
  group_by(GenusREC) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

# Get top 45 genera
top_45_genera <- genus_means$GenusREC[1:45]

# Create a temporary column for grouping and ensure proper aggregation
genus_abundance_df <- genus_abundance_df %>%
  mutate(GenusGroup = if_else(GenusREC %in% top_45_genera, GenusREC, "Other")) %>%
  group_by(Sample, GenusGroup) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  rename(GenusREC = GenusGroup)

# Get the order of genera based on mean abundance (excluding "Other")
genus_order <- genus_abundance_df %>%
  filter(GenusREC != "Other") %>%
  group_by(GenusREC) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(GenusREC)

# Add "Other" to the end of the order
genus_order <- c(genus_order, "Other")

# Convert GenusREC to factor with levels in order of abundance
genus_abundance_df$GenusREC <- factor(genus_abundance_df$GenusREC, 
                                      levels = genus_order)

# Renormalize to ensure sums are exactly 1
genus_abundance_df <- genus_abundance_df %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Generate colors for the known genera and add grey for "Other"
genera_colors <- colorRampPalette(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", 
                                    "#ff7f00", "#ffff33"))(45)
colors <- c(genera_colors, "grey70")
names(colors) <- genus_order  # Assign names to colors to maintain order

# Upper panel: Genus-level barplot (without clustering)
p1 <- ggplot(genus_abundance_df, aes(x = Sample, y = Abundance, fill = GenusREC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm")) +
  labs(fill = "Genus", y = "Relative Abundance") +
  guides(fill = guide_legend(ncol = 1))

# Middle panel: Group information
sample_data_df <- data.frame(sample_data(recons_ps))
sample_data_df$Sample <- rownames(sample_data_df)
group_colors <- c("#ea801c", "#1a80bb")  # Colors for Tumor and Healthy
names(group_colors) <- c("Tumor", "Healthy")

p2 <- ggplot(sample_data_df, aes(x = Sample, y = 1, fill = Group)) +
  geom_tile(height = 0.8) +  # Adjusted height of tiles
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10),
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")) +
  labs(fill = "Group")

# Specified bacteria in specific order
specified_bacteria <- c(
  "Delftia",
  "Reyranella",
  "Afipia",
  "Sediminibacterium",
  "Lactococcus",
  "Lactobacillus",
  "Curvibacter",
  "SYSU-D60015",
  "Streptococcus",
  "Janthinobacterium",
  "Nocardioides",
  "Streptomyces",
  "Lysobacter",
  "Bordetella",
  "Phenylobacterium"
)

# Prepare presence/absence matrix with ordered factors for specified bacteria
presence_matrix <- genus_abundance_df %>%
  filter(GenusREC %in% specified_bacteria) %>%
  mutate(
    present = Abundance > 0,
    # Create a factor with levels in the specified order
    GenusREC = factor(GenusREC, levels = specified_bacteria)
  ) %>%
  select(Sample, GenusREC, present)

# Create p3 with vertical lines for presence (with specified bacteria only)
p3 <- ggplot(presence_matrix, aes(x = Sample, y = GenusREC)) +
  geom_segment(data = subset(presence_matrix, present),
               aes(x = Sample, xend = Sample,
                   y = as.numeric(factor(GenusREC)) - 0.4,
                   yend = as.numeric(factor(GenusREC)) + 0.4),
               color = "#A9A9A9", size = 1.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10, face = "italic", hjust = 1),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 100, unit = "pt")
  ) +
  scale_y_discrete(limits = specified_bacteria)

# Combine all plots
combined_plot <- grid.arrange(p1, p2, p3, 
                              heights = c(3, 0.8, 4),  # Adjusted height ratio
                              ncol = 1,
                              padding = unit(1, "line"))

# Export figure
tryCatch({
  ggsave(here("results/Figures/Figure2C_microbiome_composition.pdf"), 
         combined_plot, width = 12, height = 12)
}, error = function(e) {
  message("Error saving PDF: ", e$message)
})

# Export raw data for figures
tryCatch({
  write.csv(genus_abundance_df %>%
              group_by(Sample, GenusREC) %>%
              summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
              pivot_wider(names_from = GenusREC, values_from = Abundance),
            here("results/Tables/Table3C_genus_abundance.csv"),
            row.names = FALSE)
}, error = function(e) {
  message("Error saving abundance data: ", e$message)
})

# Export group data
tryCatch({
  write.csv(sample_data_df %>%
              select(Sample, Group),
            here("results/Tables/Table3D_sample_groups.csv"), 
            row.names = FALSE)
}, error = function(e) {
  message("Error saving group data: ", e$message)
})

# Export presence/absence data - modified to use only specified bacteria
tryCatch({
  presence_wide <- presence_matrix %>%
    group_by(GenusREC, Sample) %>%
    summarise(present = max(present), .groups = "drop") %>%
    pivot_wider(
      names_from = Sample,
      values_from = present,
      values_fill = FALSE
    )
  
  write.csv(presence_wide, 
            here("results/Tables/Table3E_genus_presence.csv"), 
            row.names = FALSE)
}, error = function(e) {
  message("Error saving presence data: ", e$message)
})


#===============================================================================
# Figure 2D: Phylogenetic plot
#===============================================================================
# Load required libraries
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(phyloseq)
library(Biostrings)
library(phangorn)
library(DECIPHER)

# Generate phylogenetic tree from ASV sequences
message("Generating phylogenetic tree from ASV sequences...")

# Extract sequences from the tax_table of recons_ps
seqs <- tax_table(recons_ps)[,'sequence']
names(seqs) <- taxa_names(recons_ps)

# Convert to DNAStringSet and align
seqs_dna <- DNAStringSet(seqs)
names(seqs_dna) <- names(seqs)  # Explicitly set names for DNAStringSet
alignment <- AlignSeqs(seqs_dna, anchor=NA)
names(alignment) <- names(seqs)  # Set names for alignment

# Convert to matrix while preserving names
align_matrix <- as(alignment, "matrix")
rownames(align_matrix) <- names(seqs)  # Explicitly set row names for matrix

# Convert to phyDat format with preserved names
phang.align <- phyDat(align_matrix, type="DNA")

# Calculate distance matrix and create NJ tree
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

# Set the correct tip labels
treeNJ$tip.label <- names(seqs)

# Initial model fit
fit = pml(treeNJ, data=phang.align)

# Start with a simpler model first (with progress tracking)
message("Optimizing simple JC model...")
fit_simple <- optim.pml(fit, model="JC", rearrangement = "NNI",
                        control = pml.control(trace = 1))

# Proceed to more complex GTR model
message("Updating to GTR model...")
fitGTR <- update(fit_simple, model="GTR", k=4, inv=0.2)
message("Optimizing GTR model (this may take some time)...")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 1))

# Add the tree to your existing phyloseq object
message("Creating final phyloseq object with tree...")
recons_ps_tree <- merge_phyloseq(recons_ps, phy_tree(fitGTR$tree))

# Print summary of the final object
message("Summary of the final phyloseq object:")
print(recons_ps_tree)

# Save the resulting tree to the reorganized directory structure
message("Saving tree to file...")
dir.create(here("results/Tables/phylotree"), recursive = TRUE, showWarnings = FALSE)
saveRDS(fitGTR$tree, here("results/Tables/phylotree/phylogenetic_tree.rds"))
saveRDS(recons_ps_tree, here("results/Tables/phylotree/phyloseq_with_tree.rds"))

# Get tree and data
tree <- phy_tree(recons_ps_tree)
sample_df <- data.frame(sample_data(recons_ps_tree))
otu_mat <- as.matrix(otu_table(recons_ps_tree))
tax_df <- data.frame(tax_table(recons_ps_tree))

# Ensure OTU matrix has taxa as rows
if(!taxa_are_rows(recons_ps_tree)) {
  otu_mat <- t(otu_mat)
}

# Get sample groups
healthy_samples <- rownames(sample_df)[sample_df$Group == "Healthy"]
tumor_samples <- rownames(sample_df)[sample_df$Group == "Tumor"]

# Calculate relative abundance
sample_sums <- colSums(otu_mat)
rel_abund_mat <- sweep(otu_mat, 2, sample_sums, '/')

# Separate by group
healthy_abund <- rel_abund_mat[, healthy_samples, drop=FALSE]
tumor_abund <- rel_abund_mat[, tumor_samples, drop=FALSE]

# Calculate mean abundance
healthy_mean <- rowMeans(healthy_abund) * 100
tumor_mean <- rowMeans(tumor_abund) * 100

# Create long format data for abundance plotting
abundance_data <- data.frame(
  ASV_ID = rep(rownames(healthy_abund), 2),
  value = c(healthy_mean, tumor_mean),
  group = rep(c("Healthy", "Tumor"), each = length(healthy_mean))
)

# Calculate prevalence
prevalence_healthy <- rowSums(otu_mat[, healthy_samples, drop=FALSE] > 0) / length(healthy_samples) * 100
prevalence_tumor <- rowSums(otu_mat[, tumor_samples, drop=FALSE] > 0) / length(tumor_samples) * 100

# Combine tree data
tree_data <- data.frame(
  ASV_ID = rownames(otu_mat),
  prevalence_healthy = prevalence_healthy,
  prevalence_tumor = prevalence_tumor,
  PhylumREC = tax_df$PhylumREC,
  row.names = rownames(otu_mat)
)

# Create comprehensive data export
export_data <- data.frame(
  ASV_ID = rownames(otu_mat),
  Taxonomy = tax_df[, c("PhylumREC", "ClassREC", "OrderREC", "FamilyREC", "GenusREC", "SpeciesREC")],
  Prevalence_Healthy = prevalence_healthy,
  Prevalence_Tumor = prevalence_tumor,
  Mean_Abundance_Healthy = healthy_mean,
  Mean_Abundance_Tumor = tumor_mean,
  Family_Recommendation = tax_df$Family_recommendation
)

# Export raw data
write.csv(export_data, 
          file = here("results/Tables/Table3F_phylogenetic_data.csv"), 
          row.names = FALSE)

# Identify nodes with Family_recommendation = "3DB"
three_db_taxa <- rownames(tax_df)[tax_df$Family_recommendation == "3DB"]
three_db_nodes <- which(tree$tip.label %in% three_db_taxa)

# Create sector dataframe for highlighting
sector_df <- data.frame(
  node = three_db_nodes,
  sector = "3DB",
  stringsAsFactors = FALSE
)

# Create groups based on OrderREC
order_groups <- split(rownames(tax_df), tax_df$OrderREC)

# Create the base tree and bind data first
p <- ggtree(tree, layout='circular', 
            size=0.15,
            branch.length="none") %<+% tax_df

# Add OrderREC grouping and coloring
p <- groupOTU(p, order_groups, 'OrderREC') +
  aes(color=OrderREC) +
  scale_color_discrete(
    name = "OrderREC",
    breaks = sort(unique(tax_df$OrderREC))
  )

# Add 3DB highlighting
p <- p + geom_hilight(data=sector_df,
                      aes(node=node),
                      fill="grey",
                      extend=-3,
                      alpha=0.3)

# Add new scale for PhylumREC points
p <- p + new_scale_color()

# Add tip points with PhylumREC coloring
p <- p + geom_tippoint(aes(color=PhylumREC), 
                       size=1.0) +
  scale_color_manual(
    name = "PhylumREC",
    values = c(
      "Acidobacteriota" = "#800080",    
      "Actinomycetota" = "#B0171F",     
      "Armatimonadota" = "#7B68EE",     
      "Bacillota" = "#9ACD32",          
      "Bacteroidota" = "#87CEFA",       
      "Deinococcota" = "#006400",       
      "Gemmatimonadota" = "#D15FEE",    
      "Myxococcota" = "#FFC125",        
      "Pseudomonadota" = "#EE6A50"      
    ),
    guide=guide_legend(
      keywidth = 0.5, 
      keyheight = 0.5, 
      order=1,
      override.aes=list(starshape=15)
    ),
    na.translate=FALSE
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9)
  )

# Healthy Prevalence with 50% limit
p <- p + new_scale_fill() +
  geom_fruit(
    data = tree_data,
    geom = geom_tile,
    mapping = aes(y = ASV_ID, x = 1, fill = prevalence_healthy),
    offset=0.05,   
    width= 3 
  ) +
  scale_fill_gradient(
    name = "Healthy\nPrevalence (%)",
    low = "white",  
    high = "#0072b5",  
    limits = c(0, 50)
  )

# Tumor Prevalence with 50% limit
p <- p + new_scale_fill() +
  geom_fruit(
    data = tree_data,
    geom = geom_tile,
    mapping = aes(y = ASV_ID, x = 1, fill = prevalence_tumor),
    offset=0.001,   
    width= 3 
  ) +
  scale_fill_gradient(
    name = "Tumor\nPrevalence (%)",
    low = "white",  
    high = "#fc7a00",  
    limits = c(0, 50)
  )

# Add the abundance layer
p <- p + new_scale_fill() +
  # Add central axis line first
  geom_fruit(
    geom = geom_vline,
    mapping = aes(xintercept = 0),
    color = "grey50",
    size = 0.5,
    offset = 0.1
  ) +
  # Add abundance bars
  geom_fruit(
    data = abundance_data,
    geom = geom_col,
    mapping = aes(
      y = ASV_ID,
      x = ifelse(group == "Tumor", -value, value),
      fill = group
    ),
    width = 2,
    offset = 0.02,
    axis.params = list(
      axis = "x",
      text.size = 1.8,
      hjust = 1,
      vjust = 0.5,
      nbreak = 6,
      text.angle = 0
    )
  ) +
  scale_fill_manual(
    name = "Mean Abundance (%)",
    values = c(
      "Healthy" = "#66cdaa",
      "Tumor" = "#eb8c8c"
    )
  )

# Theme customization
p <- p + theme(
  legend.position = "right",
  legend.box = "vertical",
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9, face = "bold"),
  legend.key.size = unit(0.4, "cm")
)

# Add title
p <- p + ggtitle("Microbiome Composition in Healthy and Tumor Samples") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

# Save the plot with increased resolution
ggsave(here("results/Figures/Figure2D_phylogenetic_tree.pdf"), p, width = 16, height = 16, bg = "white")
ggsave(here("results/Figures/Figure2D_phylogenetic_tree.png"), p, width = 16, height = 16, dpi = 600, bg = "white")

# Clean up
detach("package:phangorn", unload=TRUE)
message("Phylogenetic tree generation and plotting completed!")


#===============================================================================
# Figure S2B: Phylum stacked bar plot
#===============================================================================

# Transform counts to relative abundance
ps_rel = transform_sample_counts(recons_ps, function(x) x/sum(x))

# Get phylum-level abundances
phylum_abundance <- tax_glom(ps_rel, taxrank="PhylumREC")
phylum_abundance_df <- psmelt(phylum_abundance)

# Calculate mean abundance for each phylum
phylum_means <- phylum_abundance_df %>%
  group_by(PhylumREC) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))

# Get top 8 phyla
top_8_phyla <- phylum_means$PhylumREC[1:8]

# Step 1: Create the PhylumGroup column
phylum_abundance_df$PhylumGroup <- ifelse(phylum_abundance_df$PhylumREC %in% top_8_phyla, 
                                          as.character(phylum_abundance_df$PhylumREC), 
                                          "Other")

# Step 2: Group and summarize
phylum_summarized <- aggregate(Abundance ~ Sample + PhylumGroup + Group, 
                               data = phylum_abundance_df,
                               FUN = sum)

# Step 3: Rename the column
names(phylum_summarized)[names(phylum_summarized) == "PhylumGroup"] <- "PhylumREC"

# Assign back to original variable name
phylum_abundance_df <- phylum_summarized

# Get the order of phyla based on mean abundance (excluding "Other")
phylum_order <- phylum_abundance_df %>%
  filter(PhylumREC != "Other") %>%
  group_by(PhylumREC) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(PhylumREC)

# Add "Other" to the end of the order
phylum_order <- c(phylum_order, "Other")

# Convert PhylumREC to factor with levels in order of abundance
phylum_abundance_df$PhylumREC <- factor(phylum_abundance_df$PhylumREC, 
                                        levels = phylum_order)

# Make sure sample names are available for ordering
phylum_abundance_df$Sample <- as.character(phylum_abundance_df$Sample)

# Create a sample order based on Group and then by Pseudomonadota abundance
sample_order <- phylum_abundance_df %>%
  filter(PhylumREC == phylum_order[1]) %>%  # Use most abundant phylum for ordering
  group_by(Sample, Group) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  arrange(Group, desc(Abundance)) %>%
  pull(Sample)

# Convert Sample to factor with the determined order
phylum_abundance_df$Sample <- factor(phylum_abundance_df$Sample, 
                                     levels = sample_order)

# Generate colors for the phyla
phyla_colors <- c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", 
  "#ff7f00", "#ffff33", "#a65628", "#f781bf", 
  "grey70"  # For "Other"
)
names(phyla_colors) <- phylum_order

# Create the stacked barplot
p_stacked <- ggplot(phylum_abundance_df, aes(x = Sample, y = Abundance, fill = PhylumREC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyla_colors) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels (too many to show)
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(0.1, "lines"),  # Reduce space between facets
    strip.text = element_text(size = 12, face = "bold"),  # Larger facet labels
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 10)
  ) +
  labs(
    title = "Phylum-Level Composition by Sample",
    x = NULL,  # No x-axis label needed
    y = "Relative Abundance",
    fill = "Phylum"
  )

# Create a group color bar data frame
sample_order_with_groups <- phylum_abundance_df %>%
  select(Sample, Group) %>%
  distinct() %>%
  arrange(factor(Sample, levels = levels(phylum_abundance_df$Sample)))

# Define group colors
group_colors <- c("Tumor" = "#ea801c", "Healthy" = "#1a80bb")

# Create a simpler version of the group annotation bar
p_group <- ggplot(sample_order_with_groups, aes(x = Sample, y = 1, fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = group_colors) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +  # Match faceting of main plot
  theme_void() +  # Simplify the theme completely
  theme(
    strip.text = element_blank(),  # Remove facet labels 
    legend.position = "none",      # Remove legend since it's in the main plot
    plot.margin = margin(0, 0, 0, 0)
  )

# Combine the two plots with the group color bar below the stacked bar plot
library(gridExtra)
combined_plot <- grid.arrange(
  p_stacked + theme(plot.margin = margin(b = 0)),  # Remove bottom margin
  p_group + theme(plot.margin = margin(t = 0)),    # Remove top margin
  heights = c(0.95, 0.05),
  ncol = 1
)

# Export the combined stacked barplot with group annotation
ggsave(here("results/Figures/FigureS2B_phylum_stacked_barplot.pdf"), 
       combined_plot, width = 12, height = 9)

# Create a mean abundance stacked bar plot (just 2 bars - one for each group)
group_means <- phylum_abundance_df %>%
  group_by(Group, PhylumREC) %>%
  summarise(mean_abundance = mean(Abundance), .groups = 'drop')

p_mean_stacked <- ggplot(group_means, aes(x = Group, y = mean_abundance, fill = PhylumREC)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = phyla_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "Mean Phylum Abundance by Group",
    x = NULL,
    y = "Mean Relative Abundance",
    fill = "Phylum"
  )

# Export the mean abundance bar plot
ggsave(here("results/Figures/FigureS2B_phylum_mean_abundance.pdf"), 
       p_mean_stacked, width = 8, height = 6)

# Export the data used for plotting
write.csv(phylum_abundance_df, 
          here("results/Tables/TableS2B_phylum_abundance_data.csv"), 
          row.names = FALSE)

# Export the mean abundance data
write.csv(group_means, 
          here("results/Tables/TableS2B_phylum_group_means.csv"), 
          row.names = FALSE)

#===============================================================================
# Figure S2C: DESeq2 genus level differential abundance
#===============================================================================
# Load required libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Function for calculating geometric mean
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Function to normalize samples
normalizeSample <- function(x) {
  x / sum(x)
}

# Set theme for plots
theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line"),
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
)

# Set alpha threshold for significance
alpha = 0.05

# Get the raw count data and taxonomy
otu_mat <- as.matrix(otu_table(recons_ps))
tax_mat <- as.matrix(tax_table(recons_ps))

# Create a data frame with counts and taxonomy
counts_df <- data.frame(otu_mat)
tax_df <- data.frame(tax_mat)

# Add GenusREC information to counts and standardize genus names
tax_df$GenusREC[is.na(tax_df$GenusREC)] <- "Unknown"
tax_df$GenusREC <- gsub("^\\s+|\\s+$", "", tax_df$GenusREC)  # Remove extra spaces
tax_df$GenusREC <- gsub("\\s+", " ", tax_df$GenusREC)        # Standardize internal spaces

# Aggregate counts by GenusREC
genus_counts <- matrix(0, 
                       nrow = length(unique(tax_df$GenusREC)), 
                       ncol = ncol(counts_df))
colnames(genus_counts) <- colnames(counts_df)
rownames(genus_counts) <- unique(tax_df$GenusREC)

# Sum counts for each genus
for (genus in rownames(genus_counts)) {
  genus_asvs <- rownames(tax_df)[tax_df$GenusREC == genus]
  genus_counts[genus,] <- colSums(counts_df[genus_asvs,, drop = FALSE])
}

# Create new phyloseq object with aggregated counts
sample_data_new <- sample_data(recons_ps)
otu_table_new <- otu_table(genus_counts, taxa_are_rows = TRUE)
tax_table_new <- tax_table(matrix(unique(tax_df$GenusREC), 
                                  ncol = 1, 
                                  dimnames = list(unique(tax_df$GenusREC), "GenusREC")))

genus_ps <- phyloseq(otu_table_new, sample_data_new, tax_table_new)

# Calculate relative abundances for visualization
genus_rel_ps <- transform_sample_counts(genus_ps, normalizeSample)

# Calculate mean relative abundance for each group
tumor_rel <- subset_samples(genus_rel_ps, Group == "Tumor")
healthy_rel <- subset_samples(genus_rel_ps, Group == "Healthy")
tumor_means <- rowMeans(otu_table(tumor_rel))
healthy_means <- rowMeans(otu_table(healthy_rel))

# Perform DESeq2 analysis
diagdds <- phyloseq_to_deseq2(genus_ps, ~ Group)

# Calculate geometric means for size factors
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)

# Estimate dispersions
diagdds <- estimateDispersions(diagdds)

# Set reference level
diagdds$Group <- relevel(diagdds$Group, ref = "Healthy")

# Perform differential abundance analysis
diagdds <- DESeq(diagdds, fitType = "local")
res <- results(diagdds)

# Sort results by adjusted p-value
res <- res[order(res$padj, na.last=NA), ]

# Export complete DESeq2 results with genus names and relative abundances
res_df <- as.data.frame(res)
res_df$GenusREC <- rownames(res_df)

# Add relative abundance information for both groups
res_df$rel_abundance_tumor <- tumor_means[res_df$GenusREC]
res_df$rel_abundance_healthy <- healthy_means[res_df$GenusREC]

# Calculate fold change in relative abundance
res_df$rel_abundance_fold_change <- res_df$rel_abundance_tumor / res_df$rel_abundance_healthy

write.csv(res_df, 
          file = here("results/Tables/TableS2C_deseq2_genus_results.csv"), 
          row.names = FALSE)

# Calculate significance for volcano plot
res_df$sig <- -log10(res_df$padj)
res_df$sig[is.infinite(res_df$sig)] <- max(res_df$sig[!is.infinite(res_df$sig)]) + 1

# Set colors for volcano plot
cols <- rep("black", nrow(res_df))
cols[res_df$padj < alpha & res_df$log2FoldChange > 0] <- "#ea801c"  # Increased in Tumor
cols[res_df$padj < alpha & res_df$log2FoldChange < 0] <- "#1a80bb"  # Increased in Healthy

# Create the volcano plot with x-axis range from -29 to +29
p2 <- ggplot(res_df, aes(x = log2FoldChange, y = sig)) +
  geom_point(aes(size = rel_abundance_tumor), color = cols, alpha = 1, shape = 19) +
  geom_text_repel(
    data = subset(res_df, padj < alpha),
    aes(label = GenusREC),
    size = 3,
    force = 15,
    max.overlaps = 15,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey",
    segment.alpha = 0.5,
    min.segment.length = 0,
    seed = 42
  ) +
  geom_hline(yintercept = -log10(alpha), color = "red", linetype = "dashed") +
  scale_size_continuous(
    name = "Relative\nAbundance\n(Tumor)",
    range = c(0.1, 6),
    limits = c(0.0005, 0.2),
    breaks = c(0.0005, 0.01, 0.2),
    labels = function(x) sprintf("%.3f", x)
  ) +
  xlim(-25, 25) +
  xlab("Fold-change (log2)") +
  ylab("Adj. p-value (-log10)") +
  ggtitle("Differential Abundance: Tumor vs. Healthy\n(GenusREC level)") +
  theme +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(color = "black"),
    legend.key = element_blank()
  ) +
  guides(size = guide_legend(
    override.aes = list(
      color = "#ea801c"
    )
  ))

# Save the plot
ggsave(here("results/Figures/FigureS2C_genus_volcano.pdf"), p2, width = 5.5, height = 4.5, units = "in")

# Export significant results with additional information
sig_res <- subset(res_df, padj < alpha)
sig_res <- sig_res[order(sig_res$padj), ]
write.csv(sig_res[, c("GenusREC", "baseMean", "log2FoldChange", "lfcSE", 
                      "stat", "pvalue", "padj", 
                      "rel_abundance_tumor", "rel_abundance_healthy", "rel_abundance_fold_change")], 
          file = here("results/Tables/TableS2C_significant_genera.csv"), 
          row.names = FALSE)

#===============================================================================
# Figure S2D: DESeq2 ASV level differential abundance
#===============================================================================
# Load required libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set alpha threshold for significance
alpha = 0.05

# Calculate relative abundances for visualization
asv_rel_ps <- transform_sample_counts(recons_ps, normalizeSample)

# Calculate mean relative abundance for each group
tumor_rel <- subset_samples(asv_rel_ps, Group == "Tumor")
healthy_rel <- subset_samples(asv_rel_ps, Group == "Healthy")
tumor_means <- rowMeans(otu_table(tumor_rel))
healthy_means <- rowMeans(otu_table(healthy_rel))

# Calculate standard errors for relative abundances
tumor_se <- apply(otu_table(tumor_rel), 1, function(x) sd(x)/sqrt(length(x)))
healthy_se <- apply(otu_table(healthy_rel), 1, function(x) sd(x)/sqrt(length(x)))

# Perform DESeq2 analysis
diagdds <- phyloseq_to_deseq2(recons_ps, ~ Group)

# Calculate geometric means for size factors
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)

# Estimate dispersions
diagdds <- estimateDispersions(diagdds)

# Set reference level
diagdds$Group <- relevel(diagdds$Group, ref = "Healthy")

# Perform differential abundance analysis
diagdds <- DESeq(diagdds, fitType = "local")
res <- results(diagdds)

# Sort results by adjusted p-value
res <- res[order(res$padj, na.last=NA), ]

# Create results dataframe
res_df <- as.data.frame(res)
res_df$ASV_ID <- rownames(res_df)

# Add taxonomy information
tax_info <- data.frame(ASV_ID = rownames(tax_table(recons_ps)),
                       KingdomREC = tax_table(recons_ps)[,"KingdomREC"],
                       PhylumREC = tax_table(recons_ps)[,"PhylumREC"],
                       ClassREC = tax_table(recons_ps)[,"ClassREC"],
                       OrderREC = tax_table(recons_ps)[,"OrderREC"],
                       FamilyREC = tax_table(recons_ps)[,"FamilyREC"],
                       GenusREC = tax_table(recons_ps)[,"GenusREC"],
                       SpeciesREC = tax_table(recons_ps)[,"SpeciesREC"])

res_df <- merge(res_df, tax_info, by = "ASV_ID", all.x = TRUE)

# Add relative abundance information for both groups
res_df$rel_abundance_tumor <- tumor_means[res_df$ASV_ID]
res_df$rel_abundance_healthy <- healthy_means[res_df$ASV_ID]
res_df$rel_abundance_tumor_se <- tumor_se[res_df$ASV_ID]
res_df$rel_abundance_healthy_se <- healthy_se[res_df$ASV_ID]
res_df$rel_abundance_fold_change <- res_df$rel_abundance_tumor / res_df$rel_abundance_healthy

# Save complete results with all relative abundance information
write.csv(res_df, 
          file = here("results/Tables/TableS2D_deseq2_asv_results.csv"), 
          row.names = FALSE)

# Calculate significance for volcano plot
res_df$sig <- -log10(res_df$padj)
res_df$sig[is.infinite(res_df$sig)] <- max(res_df$sig[!is.infinite(res_df$sig)]) + 1

# Set colors for volcano plot
cols <- rep("black", nrow(res_df))
cols[res_df$padj < alpha & res_df$log2FoldChange > 0] <- "#ea801c"  # Increased in Tumor
cols[res_df$padj < alpha & res_df$log2FoldChange < 0] <- "#1a80bb"  # Increased in Healthy

# Create ASV labels combining taxonomy information (with NA handling)
res_df$plot_label <- sapply(1:nrow(res_df), function(i) {
  genus <- res_df$GenusREC[i]
  species <- res_df$SpeciesREC[i]
  asv_suffix <- substr(res_df$ASV_ID[i], nchar(res_df$ASV_ID[i])-3, nchar(res_df$ASV_ID[i]))
  
  if (is.na(genus) || is.na(species) || 
      genus %in% c(NA, "Not_Assigned") || species %in% c(NA, "Not_Assigned")) {
    paste0("Not_Assigned_", asv_suffix)
  } else {
    paste0(genus, " ", species, "_", asv_suffix)
  }
})

# Create the volcano plot
p_asv_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = sig)) +
  # Add reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
  # Add points
  geom_point(aes(size = rel_abundance_tumor), color = cols, alpha = 1, shape = 19) +
  # Add labels
  geom_text_repel(
    data = subset(res_df, padj < alpha),
    aes(label = plot_label),
    size = 2,
    force = 15,
    max.overlaps = 15,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey",
    segment.alpha = 0.5,
    min.segment.length = 0,
    seed = 42
  ) +
  # Add significance threshold line
  geom_hline(yintercept = -log10(alpha), color = "red", linetype = "dashed") +
  # Configure size scale
  scale_size_continuous(
    name = "Relative\nAbundance\n(Tumor)",
    range = c(0.1, 6),
    limits = c(0.0001, 0.16),
    breaks = c(0.0001, 0.004, 0.16),
    labels = function(x) sprintf("%.3f", x)
  ) +
  # Add labels
  xlab("Fold-change (log2)") +
  ylab("Adj. p-value (-log10)") +
  ggtitle("Differential Abundance: Tumor vs. Healthy\n(ASV level)") +
  # Add theme
  theme +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(color = "black"),
    legend.key = element_blank()
  ) +
  guides(size = guide_legend(
    override.aes = list(
      color = "#ea801c"
    )
  ))

# Save the plot
ggsave(here("results/Figures/FigureS2D_asv_volcano.pdf"), p_asv_volcano, width = 5.5, height = 4.5, units = "in")

# Calculate and save summary statistics with relative abundance info
summary_stats <- data.frame(
  total_ASVs = nrow(res_df),
  significant_ASVs = sum(res_df$padj < alpha, na.rm = TRUE),
  increased_in_tumor = sum(res_df$padj < alpha & res_df$log2FoldChange > 0, na.rm = TRUE),
  increased_in_healthy = sum(res_df$padj < alpha & res_df$log2FoldChange < 0, na.rm = TRUE),
  mean_rel_abundance_tumor = mean(res_df$rel_abundance_tumor, na.rm = TRUE),
  mean_rel_abundance_healthy = mean(res_df$rel_abundance_healthy, na.rm = TRUE)
)

write.csv(summary_stats, 
          file = here("results/Tables/TableS2D_asv_analysis_summary.csv"), 
          row.names = FALSE)

# Export significant results with additional information
sig_res <- subset(res_df, padj < alpha)
sig_res <- sig_res[order(sig_res$padj), ]

# Select columns for significant results export
selected_columns <- c("ASV_ID", "KingdomREC", "PhylumREC", "ClassREC", 
                      "OrderREC", "FamilyREC", "GenusREC", "SpeciesREC",
                      "baseMean", "log2FoldChange", "lfcSE", "stat", 
                      "pvalue", "padj", "rel_abundance_tumor", 
                      "rel_abundance_healthy", "rel_abundance_fold_change",
                      "rel_abundance_tumor_se", "rel_abundance_healthy_se")

write.csv(sig_res[, selected_columns], 
          file = here("results/Tables/TableS2D_significant_asvs.csv"), 
          row.names = FALSE)

#==============================================================================
# FINAL STEP: Print out reorganization summary
#==============================================================================

message("\n==================== FILE REORGANIZATION SUMMARY ====================\n")

message("The code has been reorganized with the following structure:")
message("- Figures are saved to: results/Figures/")
message("- Data tables are saved to: results/Tables/")
message("- File naming follows the section structure in the code")

message("\nMain figures:")
message("- Figure 2B: Beta diversity PCoA plot")
message("- Figure 2C: Stacked bar plot and presence/absence plot")
message("- Figure 2D: Phylogenetic tree visualization")

message("\nSupplementary figures:")
message("- Figure S2B: Phylum stacked bar plot")
message("- Figure S2C: Genus-level differential abundance")
message("- Figure S2D: ASV-level differential abundance")

message("\nCode structure has been cleaned and organized according to figures and analyses.")
message("\n==================================================================\n")