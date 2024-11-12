# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM

library(dplyr)
library(gt)
library(here)
library(ggplot2)

analyse_data <- function(import_rds) {
  
  # Load the cleaned data
  cons_tum <- readRDS(file = import_rds)

  # Normalise for the heterozygous in tum relative to constit
  cons_tum <- cons_tum %>%
    mutate(
      VAF.tum = round(
        if_else(
          VAF.cons > 0.4 & VAF.cons < 0.6,
          (VAF.tum * 0.5) / VAF.cons,
          VAF.tum
        ),
        2
      )
    )
  
  
  # Classify the Loss of Heterozygousity in TRANS or CIS 
  cons_tum <- cons_tum %>%
    mutate(
      LOH = case_when(
        VAF.tum / VAF.cons >= 1.2 ~ "TRANS",
        VAF.tum / VAF.cons <= 0.8 ~ "CIS",
        TRUE ~ "" # Default value if none of the above conditions are met
      )
    )
  
  #Estimate the %tumoral
  cons_tum <- cons_tum %>%
    mutate(
      
      `%tumoral` = case_when(
        LOH == "CIS" ~ (200 * VAF.tum - 100) / (VAF.tum - 1),
        
        LOH == "TRANS" ~ (200 * VAF.tum - 100) / VAF.tum
      )
    )
  
  
  
  # Select relevant columns
  cons_tum <- cons_tum %>%
    select(Pos., Gene.cons, c..HGVS.cons, VAF.cons, VAF.tum, LOH, `%tumoral`, chrom_num, pos_num) %>%
    # Sort by chrom_num first, then by pos_num
    arrange(chrom_num, pos_num) %>%
    # Remove temporary sorting columns
    select(-chrom_num, -pos_num) %>%
    mutate(`%tumoral` = round(`%tumoral`, 2))
  
    
  
  # Extract lists of %tumoral values based on LOH classification
  data_for_plot <- cons_tum %>%
    filter(!is.na(`%tumoral`)) %>%
    select(LOH, `%tumoral`)
  
  # Calculate mean and standard deviation for each LOH classification
  summary_stats <- data_for_plot %>%
    group_by(LOH) %>%
    summarise(
      Mean = mean(`%tumoral`, na.rm = TRUE),
      SD = sd(`%tumoral`, na.rm = TRUE)
    )
  
  # Calculate mean of the two classes
  mean_CIS_TRANS <- cons_tum %>%
    filter(LOH %in% c("CIS", "TRANS")) %>%
    summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
    pull(Mean)
  
  #calculate the VAFtheoric
  if (cons_tum$LOH == "CIS") {
    cons_tum %>%
      mutate(VAFtheoric = (100 - `%tumoral`) / (200 - `%tumoral`))
  } else if (cons_tum$LOH == "TRANS") {
    cons_tum %>% 
      mutate(VAFtheoric = 100 / (200 - `%tumoral`))
  }

  
  # Create the boxplot
  boxplot <- ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
    geom_boxplot(varwidth = TRUE) +  # No need for fill=LOH here, it's already set in aes()
    geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
    geom_text(data = summary_stats, aes(x = LOH, y = Mean, label = paste("Mean:", round(Mean, 2), "±", round(SD, 2))), 
              vjust = -0.5, hjust = 1, color = "#4D4D4D", size = 3.5) +  
    labs(title = "Pourcentage estimé de cellules tumorales par classification LOH",
         x = NULL,
         y = "% Tumoral") +
    scale_fill_manual(values = c("CIS" = "#d4f1bc", "TRANS" = "#ffcccb")) +
    theme_minimal()
  
  # Mean of the two classes
  mean_CIS_TRANS <- cons_tum %>%
    filter(LOH %in% c("CIS", "TRANS")) %>%
    summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
    pull(Mean)
  

  
  return(list(cons_tum = cons_tum, boxplot = boxplot))
}


