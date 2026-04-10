# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM

library(dplyr)

analyse_data <- function(import_rds) {
  
  cons_tum <- readRDS(file = import_rds)
  
  cons_tum <- cons_tum %>%
    mutate(
      VAF.tum = round(
        if_else(
          !is.na(VAF.cons) & VAF.cons > 0.4 & VAF.cons < 0.6,
          (VAF.tum * 0.5) / VAF.cons,
          VAF.tum
        ),
        2
      ),
      ratio = VAF.tum / VAF.cons,
      LOH = case_when(
        !is.na(ratio) & ratio >= 1.2 ~ "TRANS",
        !is.na(ratio) & ratio <= 0.8 ~ "CIS",
        TRUE ~ NA_character_
      ),
      `%tumoral` = case_when(
        LOH == "CIS" & !is.na(VAF.tum) & VAF.tum != 1 ~ (200 * VAF.tum - 100) / (VAF.tum - 1),
        LOH == "TRANS" & !is.na(VAF.tum) & VAF.tum != 0 ~ (200 * VAF.tum - 100) / VAF.tum,
        TRUE ~ NA_real_
      ),
      `%tumoral` = if_else(`%tumoral` >= 0 & `%tumoral` <= 100, `%tumoral`, NA_real_),
      VAFtheoTRANS = if_else(
        !is.na(`%tumoral`),
        `%tumoral` / (`%tumoral` + 2 * (100 - `%tumoral`)),
        NA_real_
      ),
      VAFtheoPASdeLOH = if_else(
        !is.na(`%tumoral`),
        `%tumoral` / (`%tumoral` * 2 + (100 - `%tumoral`) * 2),
        NA_real_
      )
    ) %>%
    select(
      Pos., Gene, c..HGVS, VAF.cons, VAF.tum, LOH, `%tumoral`,
      VAFtheoTRANS, VAFtheoPASdeLOH, chrom_num, pos_num
    ) %>%
    arrange(chrom_num, pos_num) %>%
    mutate(`%tumoral` = round(`%tumoral`, 2))
  
  cons_tum
}
