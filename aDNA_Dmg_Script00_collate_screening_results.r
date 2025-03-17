#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#                    aDNA Damage through time -  Script 00                     #
#                         Collating screening results                          #
#                                                                              #
#   Collate results from aDNA screening generated with Latorre et.al, 2020     #
#     "Plant aDNA pipeline" (https://doi.org/10.1002/cppb.20121) protocol:     #
#               https://gitlab.com/smlatorreo/plant-adna-pipeline              #
#                                                                              #
# Note: sample_metadata.txt can be provided in a tab-separated file, ensure    #
#       names under the column "Sample" match the name of the samples          #
#       used during screening.                                                 #
#------------------------------------------------------------------------------#

# Set warning handling and messages
options(warn = 1)  # Print warnings as they occur
print("Starting script execution...")

# Install required packages
install.packages(c("dplyr", "readr", "MASS"))

# Load necessary libraries
library(readr)
library(dplyr)
library(MASS)
print("Libraries loaded successfully")

# Define the base directories
base_dir <- "."  # Change this to your base directory if needed
mapdamage_dir <- file.path(base_dir, "5_aDNA_characteristics")
mapping_dir <- file.path(base_dir, "4_mapping")

#------------------------------------------------------------------------------#
# Helper Functions                                                             #
#------------------------------------------------------------------------------#

# Function to process flagstat files
process_flagstat <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(c(NA, NA))
  }

  lines <- readLines(file_path)
  mapped_line <- lines[grep("mapped \\(", lines)[1]]

  if (length(mapped_line) == 0) {
    warning("No mapped reads line found in flagstat file")
    return(c(NA, NA))
  }

  # Extract number of mapped reads
  no_merged <- as.numeric(strsplit(mapped_line, " ")[[1]][1])

  # Extract endogenous percentage
  percentage_match <- regexpr("\\([0-9.]+%", mapped_line)
  if (percentage_match > 0) {
    endogenous <- as.numeric(gsub("[^0-9.]", "",
                                  substr(mapped_line, percentage_match + 1,
                                         percentage_match +
                                           attr(percentage_match,
                                                "match.length") - 1)))
  } else {
    warning("No percentage found in mapped reads line")
    endogenous <- NA
  }

  return(c(endogenous, no_merged))
}

# Function to extract duplication rate
get_duplication_rate <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NA)
  }

  lines <- readLines(file_path)
  dup_line <- tail(lines, 1)
  dup_rate <- as.numeric(gsub("Duplication Rate: ", "", dup_line))
  return(dup_rate)
}

# Function to calculate lambda (damage fraction per site)
calculate_lambda <- function(lengths, frequencies) {
  # Create data frame for regression
  regression_data <- data.frame(
    length = lengths,
    frequency = frequencies
  ) %>%
    filter(frequency > 0) %>%  # Remove zeros before log transformation
    mutate(log_freq = log(frequency))

  # Perform linear regression
  model <- try(lm(log_freq ~ length, data = regression_data), silent = TRUE)

  if (inherits(model, "try-error") || is.null(model)) {
    warning("Failed to fit linear model for lambda calculation")
    return(list(lambda = NA, r_squared = NA))
  }

  # Lambda is the negative of the slope
  lambda <- -coef(model)[2]
  r_squared <- summary(model)$r.squared

  return(list(
    lambda = lambda,
    r_squared = r_squared
  ))
}

# Function to process length distribution and calculate metrics
process_lgdistribution <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(list(
      log_stats = c(NA, NA),
      lambda_stats = c(NA, NA),
      median_size = NA,
      raw_distribution = NULL
    ))
  }

  # Read the length distribution data
  lgdist_data <- read_delim(file_path, delim = "\t", 
                           comment = "#", col_names = TRUE) %>%
    setNames(trimws(names(.)))

  # Create expanded lengths vector for calculations
  expanded_lengths <- rep(lgdist_data$Length, lgdist_data$Occurences)

  # Fit lognormal distribution
  fit <- try(fitdistr(expanded_lengths, "lognormal"), silent = TRUE)

  if (inherits(fit, "try-error")) {
    warning(paste("Failed to fit lognormal distribution for:", file_path))
    log_mean <- NA
    log_median <- NA
  } else {
    log_mean <- fit$estimate['meanlog']
    log_median <- exp(fit$estimate['meanlog'])
  }

  # Calculate median size
  median_size <- median(expanded_lengths)

  # Calculate lambda
  lambda_results <- calculate_lambda(lgdist_data$Length, lgdist_data$Occurences)

  return(list(
    log_stats = c(log_mean, log_median),
    lambda_stats = c(lambda_results$lambda, lambda_results$r_squared),
    median_size = median_size,
    raw_distribution = lgdist_data
  ))
}

# Function to extract damage frequencies
extract_damage_metrics <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(c(NA, NA, NA))
  }

  data <- try({
    read.table(file_path, header = TRUE, check.names = FALSE)
  }, silent = TRUE)

  if(inherits(data,"try-error") || nrow(data) < 3 || !"5pC>T" %in% names(data)) {
    warning(paste("Invalid data format in:", file_path))
    return(c(NA, NA, NA))
  }

  # Return first three positions
  return(data[1:3, "5pC>T"])
}

#------------------------------------------------------------------------------
# Main Processing Function
#------------------------------------------------------------------------------

process_samples <- function(output_file, metadata_file = NULL) {
  # Read metadata if provided
  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
  } else {
    metadata <- NULL
  }

  # Get all flagstat files
  flagstat_files <- list.files(path = mapping_dir,
                             pattern = "\\.flagstat\\.log$",
                             full.names = TRUE)

  # Initialize results list
  results_list <- list()

  # Process each sample
  for (flagstat_file in flagstat_files) {
    sample_name <- sub("\\.flagstat\\.log$", "", basename(flagstat_file))
    message("Processing sample: ", sample_name)

    # Process flagstat file
    flagstat_results <- process_flagstat(flagstat_file)

    # Get duplication rate
    dup_file <- file.path(mapping_dir, paste0(sample_name, ".mapped.sorted.log"))
    dup_rate <- get_duplication_rate(dup_file)

    # Process MapDamage files
    sample_mapdamage_dir <- file.path(mapdamage_dir, sample_name)

    # Get size metrics and lambda
    lgdist_file <- file.path(sample_mapdamage_dir, "lgdistribution.txt")
    size_results <- process_lgdistribution(lgdist_file)

    # Get damage metrics
    ctot_file <- file.path(sample_mapdamage_dir, "5pCtoT_freq.txt")
    dmg_metrics <- extract_damage_metrics(ctot_file)

    # Create results row
    results_list[[sample_name]] <- data.frame(
      Sample = sample_name,
      Endogenous_fraction = flagstat_results[1],
      NO_MERGED_READS = flagstat_results[2],
      DUP_rate = dup_rate,
      MEDIAN_SIZE = size_results$median_size,
      Log_Mean = size_results$log_stats[1],
      Log_Median = size_results$log_stats[2],
      Lambda = size_results$lambda_stats[1],
      Lambda_R_squared = size_results$lambda_stats[2],
      `5P_DMG_POS1` = dmg_metrics[1],
      `5P_DMG_POS2` = dmg_metrics[2],
      `5P_DMG_POS3` = dmg_metrics[3]
    )

    # Save raw distribution data separately if needed
    if (!is.null(size_results$raw_distribution)) {
      write.table(size_results$raw_distribution,
                 file = file.path(dirname(output_file), 
                                paste0(sample_name, "_length_dist.txt")),
                 sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }

  # Combine all results
  results <- bind_rows(results_list)

  # Merge with metadata if available
  if (!is.null(metadata)) {
    results <- left_join(results, metadata, by = "Sample")
  }

  # Write results
  write.table(results, 
              file = output_file, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)

  message("Analysis complete. Results saved in ", output_file)

  # Cleanup length distribution files
  cleanup_length_files(dirname(output_file))

  return(invisible(results))
}

#------------------------------------------------------------------------------
# Cleanup Function to Remove Length Distribution Files
#------------------------------------------------------------------------------

cleanup_length_files <- function(output_dir) {
  length_files <- list.files(path = output_dir,
                             pattern = "_length_dist.txt$", full.names = TRUE)
  if (length(length_files) > 0) {
    file.remove(length_files)
    message("Removed temporary length distribution files.")
  } else {
    message("No length distribution files found for removal.")
  }
}

#------------------------------------------------------------------------------
# Script Execution
#------------------------------------------------------------------------------

# Define output files
output_file <- "aDNA_damage_screening.txt"
metadata_file <- "sample_metadata.txt"  # Optional, set to NULL if not available

# Run analysis
results <- process_samples(output_file, metadata_file)

print("Data preparation complete!")