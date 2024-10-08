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
  
  View(constit)
  
  # Set the list of genes and columns to filter out
  excluded_genes <- c("CYP2D6", "CYP1A2", "CYP2C19", "CYP3A4", "CYP3A5", "PMS2", "SDHA")
  retained_columns <- c("Gene", "Transcript", "Pos.", "Type", "Nuc.Change", "Coverage", "AA.Change", "c..HGVS", "p..HGVS")
  
  # Filter out the unwanted columns and genes
  constit_filtered <- constit %>%
    filter(!Gene %in% excluded_genes) %>%
    select(all_of(retained_columns))
  
  tumoral_filtered <- tumoral %>%
    filter(!Gene %in% excluded_genes) %>%
    select(all_of(retained_columns))
  
  # Merge the tables
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
    )
  
  saveRDS(cons_tum, file = "cons_tum_cleaned.rds")
  
}


