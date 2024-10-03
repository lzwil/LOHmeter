# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# analyse_data.R

library(dplyr)
library(gt)
library(here)
library(ggplot2)

analyse_data <- function(import_rds) {

  # Load the cleaned data
  cons_tum <- readRDS(file = import_rds)
  
  # Create a gt table with color formatting
  # Apply conditional color formatting row by row
  cons_tum <- cons_tum %>%
    select(Pos., Gene.cons, Coverage.cons, Coverage.tum, LOH, `%tumoral`)

  
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
  
  
  # Create the boxplot
  boxplot <- ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
               geom_boxplot(varwidth = TRUE, fill="plum") +
               geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
               geom_text(data = summary_stats, aes(x = LOH, y = Mean, label = paste("Mean:", round(Mean, 2), "±", round(SD, 2))), 
                         vjust = -0.5, hjust = 1, color = "#4D4D4D", size = 3.5) +  
               labs(title = "Box Plot",
                   subtitle = "% Tumoral by LOH Classification",
                   x = NULL,
                   y = "% Tumoral") +
               theme_minimal()
  
  # Mean of the two classes
  mean_CIS_TRANS <- cons_tum %>%
    filter(LOH %in% c("CIS", "TRANS")) %>%
    summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
    pull(Mean)
  
  return(list(cons_tum = cons_tum, boxplot = boxplot))
}

