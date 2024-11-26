# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# import_data.R

library(dplyr)
library(readr)
library(gt)
library(ggplot2)
library(here)
library(stringr)

import_data <- function(constit, tumoral) {
  
  constit <- read_delim(file = constit, delim = "\t", na = "", trim_ws = TRUE)
  tumoral <- read_delim(file = tumoral, delim = "\t", na = "", trim_ws = TRUE)
  
  # Replace blank spaces in column names with dots
  colnames(constit) <- gsub(" ", ".", colnames(constit))
  colnames(tumoral) <- gsub(" ", ".", colnames(tumoral))
  
  # Set the list of genes and columns to filter out
  excluded_genes <- c("CYP2D6", "CYP1A2", "CYP2C19", "CYP3A4", "CYP3A5", "PMS2", "SDHA")
  retained_columns <- c("Gene", "Transcript", "Pos.", "Type", "Nuc.Change", "Coverage", "AA.Change", "c..HGVS", "p..HGVS")
  
  # Filter out the unwanted columns and genes
  constit_filtered <- constit %>%
    filter(!Gene %in% excluded_genes) %>%
    select(all_of(retained_columns)) %>%
    filter(!grepl('CNV', Coverage))
  
  tumoral_filtered <- tumoral %>%
    filter(!Gene %in% excluded_genes) %>%
    select(all_of(retained_columns))   %>%
    filter(!grepl('CNV', Coverage))
  
  # Table for the new tumoral variants
  unique_tumoral <- tumoral_filtered %>%
    anti_join(constit_filtered, by = c("Pos.")) %>% 
    mutate(Pos. = str_extract(Pos., "chr[XY\\d]+:g\\.\\d+")) %>%
    arrange(desc(Pos.)) %>%
    mutate(
      # Extract the numeric part of the chromosome, handling cases like chrX, chrY, etc.
      chrom_num = as.numeric(str_extract(Pos., "(?<=chr)[0-9]+")),
      # Extract the position number from Pos.
      pos_num = as.numeric(str_extract(Pos., "(?<=\\.g\\.)[0-9]+"))
    )
  
  unique_tumoral$Coverage <- as.character(unique_tumoral$Coverage)
  # Keep only the allele percentage and make a frequency instead
  unique_tumoral <- unique_tumoral %>%
    mutate(
      Coverage = sapply(strsplit(Coverage, "%"), function(x) as.numeric(x[1])),
      Coverage = ifelse(is.na(Coverage), 0, Coverage) / 100 
    ) %>%
    rename(
      VAF.tum = Coverage
    ) %>%
    select(Pos., Gene, c..HGVS, VAF.tum, chrom_num, pos_num) %>%
    # Sort by chrom_num first, then by pos_num
    arrange(chrom_num, pos_num) %>%
    # Remove temporary sorting columns
    select(-chrom_num, -pos_num)
  
  # Merge the tables for the constit/tumoral comparison
  cons_tum = left_join(constit_filtered, tumoral_filtered, by = "Pos.", suffix = c(".cons", ".tum"))
  
  # Keep the genomic position and order by it
  cons_tum <- cons_tum %>% 
    mutate(Pos. = str_extract(Pos., "chr[XY\\d]+:g\\.\\d+")) %>%
    arrange(desc(Pos.)) %>%
    mutate(
      # Extract the numeric part of the chromosome, handling cases like chrX, chrY, etc.
      chrom_num = as.numeric(str_extract(Pos., "(?<=chr)[0-9]+")),
      # Extract the position number from Pos.
      pos_num = as.numeric(str_extract(Pos., "(?<=\\.g\\.)[0-9]+"))
    )
  
  # Convert columns to character if they are not already
  cons_tum$Coverage.cons <- as.character(cons_tum$Coverage.cons)
  cons_tum$Coverage.tum <- as.character(cons_tum$Coverage.tum)
  
  # Keep only the allele percentage and make a frequency instead
  cons_tum <- cons_tum %>%
    mutate(
      Coverage.cons = sapply(strsplit(Coverage.cons, "%"), function(x) as.numeric(x[1])),
      Coverage.tum = sapply(strsplit(Coverage.tum, "%"), function(x) as.numeric(x[1])),
      Coverage.cons = Coverage.cons / 100,  # Divide by 100
      Coverage.tum = Coverage.tum / 100      # Divide by 100
    ) %>%
    rename(
      VAF.cons = Coverage.cons,  # Rename Coverage.cons to VAF.cons
      VAF.tum = Coverage.tum      # Rename Coverage.tum to VAF.tum
    ) %>%
    filter(VAF.cons >= 0.45 & VAF.cons <= 0.55)
  
  # Identify variants from cons_tum with VAF.tum < 0.1
  low_coverage_variants <- cons_tum %>%
    filter(VAF.cons < 0.1 & VAF.tum > 0.1) %>%
    select(Pos., Gene.tum, c..HGVS.tum, VAF.tum) %>%
    rename(
      Gene = Gene.tum,
      c..HGVS = c..HGVS.tum,
      VAF.tum = VAF.tum
    )
  
  # Add low-coverage variants to unique_tumoral
  unique_tumoral <- unique_tumoral %>%
    bind_rows(low_coverage_variants) %>%
    distinct() %>%  # Remove duplicates
    arrange(Pos.)
  
  saveRDS(cons_tum, file = "cons_tum_cleaned.rds")
  saveRDS(unique_tumoral, file = "unique_tumoral.rds")
}