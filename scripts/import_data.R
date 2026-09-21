# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# import_data.R

library(dplyr)
library(readr)
library(stringr)

parse_vaf <- function(x) {
  x <- as.character(x)
  value <- suppressWarnings(as.numeric(str_extract(x, "[0-9]+(?:[.,][0-9]+)?")))
  value / 100
}

extract_chr_info <- function(pos_col) {
  chrom <- str_extract(pos_col, "(?<=chr)[A-Za-z0-9]+")
  pos_num <- suppressWarnings(as.numeric(str_extract(pos_col, "(?<=\\.g\\.)[0-9]+")))
  chrom_num <- dplyr::case_when(
    chrom == "X" ~ 23,
    chrom == "Y" ~ 24,
    chrom %in% c("M", "MT") ~ 25,
    TRUE ~ suppressWarnings(as.numeric(chrom))
  )
  
  tibble(chrom = chrom, chrom_num = chrom_num, pos_num = pos_num)
}

import_data <- function(constit, tumoral,
                        output_cons_tum = "cons_tum_cleaned.rds",
                        output_unique_tumoral = "unique_tumoral.rds") {
  
  constit <- read_delim(file = constit, delim = "\t", na = "", trim_ws = TRUE, show_col_types = FALSE)
  tumoral <- read_delim(file = tumoral, delim = "\t", na = "", trim_ws = TRUE, show_col_types = FALSE)
  
  clean_colnames <- function(x) {
    x <- gsub("[\u00A0\uFEFF]", " ", x)  # espaces insécables / BOM -> espace normal
    x <- trimws(x)
    gsub(" +", ".", x)
  }
  
  colnames(constit) <- clean_colnames(colnames(constit))
  colnames(tumoral) <- clean_colnames(colnames(tumoral))
  
  excluded_genes <- c("CYP2D6", "CYP1A2", "CYP2C19", "CYP3A4", "CYP3A5")
  retained_columns <- c("Gene", "Transcript", "Pos.", "Type", "Nuc.Change", "Coverage", "AA.Change", "c..HGVS", "p..HGVS")
  join_keys <- c("Pos.", "Gene", "c..HGVS")
  
  for (rc in retained_columns) {
    if (!(rc %in% names(constit))) {
      stop(
        "Colonne '", rc, "' introuvable dans le fichier constitutionnel. ",
        "Colonnes disponibles : ", paste(names(constit), collapse = ", ")
      )
    }
    if (!(rc %in% names(tumoral))) {
      stop(
        "Colonne '", rc, "' introuvable dans le fichier tumoral. ",
        "Colonnes disponibles : ", paste(names(tumoral), collapse = ", ")
      )
    }
  }
  
  constit_filtered <- constit %>%
    filter(!Gene %in% excluded_genes) %>%
    select(any_of(retained_columns)) %>%
    filter(!grepl("CNV", Coverage)) %>%
    mutate(
      Pos. = str_extract(Pos., "chr[XYMTR0-9]+:g\\.[0-9]+"),
      VAF = parse_vaf(Coverage)
    )
  
  tumoral_filtered <- tumoral %>%
    filter(!Gene %in% excluded_genes) %>%
    select(any_of(retained_columns)) %>%
    filter(!grepl("CNV", Coverage)) %>%
    mutate(
      Pos. = str_extract(Pos., "chr[XYMTR0-9]+:g\\.[0-9]+"),
      VAF = parse_vaf(Coverage)
    )
  
  unique_tumoral <- tumoral_filtered %>%
    anti_join(constit_filtered, by = join_keys) %>%
    bind_cols(extract_chr_info(.$Pos.)) %>%
    transmute(
      Pos.,
      Gene,
      c..HGVS,
      VAF.tum = VAF,
      chrom_num,
      pos_num
    ) %>%
    arrange(chrom_num, pos_num) %>%
    select(-chrom_num, -pos_num)
  
  cons_tum_all <- constit_filtered %>%
    select(all_of(join_keys), Transcript, Type, Nuc.Change, AA.Change, p..HGVS, VAF) %>%
    rename(
      Transcript.cons = Transcript,
      Type.cons = Type,
      Nuc.Change.cons = Nuc.Change,
      AA.Change.cons = AA.Change,
      p..HGVS.cons = p..HGVS,
      VAF.cons = VAF
    ) %>%
    left_join(
      tumoral_filtered %>%
        select(all_of(join_keys), Transcript, Type, Nuc.Change, AA.Change, p..HGVS, VAF) %>%
        rename(
          Transcript.tum = Transcript,
          Type.tum = Type,
          Nuc.Change.tum = Nuc.Change,
          AA.Change.tum = AA.Change,
          p..HGVS.tum = p..HGVS,
          VAF.tum = VAF
        ),
      by = join_keys
    ) %>%
    bind_cols(extract_chr_info(.$Pos.))
  
  low_coverage_variants <- cons_tum_all %>%
    filter(!is.na(VAF.cons), !is.na(VAF.tum), VAF.cons < 0.10, VAF.tum > 0.10) %>%
    transmute(
      Pos.,
      Gene = Gene,
      c..HGVS = c..HGVS,
      VAF.tum = VAF.tum
    )
  
  unique_tumoral <- unique_tumoral %>%
    bind_rows(low_coverage_variants) %>%
    distinct(Pos., Gene, c..HGVS, .keep_all = TRUE) %>%
    bind_cols(extract_chr_info(.$Pos.)) %>%
    arrange(chrom_num, pos_num) %>%
    select(-chrom, -chrom_num, -pos_num)
  
  cons_tum <- cons_tum_all %>%
    filter(!is.na(VAF.cons), VAF.cons >= 0.40, VAF.cons <= 0.60) %>%
    arrange(chrom_num, pos_num)
  
  saveRDS(cons_tum, file = output_cons_tum)
  saveRDS(unique_tumoral, file = output_unique_tumoral)
  
  invisible(list(cons_tum = cons_tum, unique_tumoral = unique_tumoral))
}
