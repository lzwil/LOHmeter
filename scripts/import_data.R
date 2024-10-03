# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# import_data.R

library(dplyr)
library(gt)
library(ggplot2)
library(here)
library(stringr)

import_data <- function(constit, tumoral) {

  constit <- read.table(constit, header = TRUE, sep = "\t", fill = TRUE)
  tumoral <- read.table(tumoral, header = TRUE, sep = "\t", fill = TRUE)
  
  
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
    arrange(desc(Pos.))
  
  # Convert columns to character if they are not already
  cons_tum$Coverage.cons <- as.character(cons_tum$Coverage.cons)
  cons_tum$Coverage.tum <- as.character(cons_tum$Coverage.tum)
  
  # Keep only the allele percentage
  cons_tum <- cons_tum %>%
    mutate(
      Coverage.cons = sapply(strsplit(Coverage.cons, "%"), function(x) as.numeric(x[1])),
      Coverage.tum = sapply(strsplit(Coverage.tum, "%"), function(x) as.numeric(x[1])),
    )
  
  # Normalise for the heterozygous in constit
  cons_tum <- cons_tum %>%
    mutate(
      Coverage.tum = if_else(
        Coverage.cons > 40 & Coverage.cons < 60,
        (Coverage.tum * 50) / Coverage.cons,
        Coverage.tum
      )
    )
  
    
  # Classify the Loss of Heterozygousity in TRANS or CIS 
  cons_tum <- cons_tum %>%
    mutate(
      LOH = case_when(
        Coverage.tum >= 10 & Coverage.tum <= 40 ~ "CIS",
        Coverage.tum >= 60 & Coverage.tum <= 90 ~ "TRANS",
        TRUE ~ "" # Default value if none of the above conditions are met
      )
    )
  
  #Estimate the %tumoral
  cons_tum <- cons_tum %>%
    mutate(
      `%tumoral` = case_when(
        LOH == "CIS" ~ abs((2 * (Coverage.tum / 100) - 1) / (1 - (Coverage.tum / 100))),
        LOH == "TRANS" ~ 2 - (1 / Coverage.tum) * 100
        )
    )

  saveRDS(cons_tum, file = "cons_tum_cleaned.rds")
  
}
