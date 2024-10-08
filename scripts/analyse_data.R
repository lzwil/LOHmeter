# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# analyse_data.Rdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAADbElEQVR4XuWXW0gUURjHfez22MyIiTPjQ82M0UtQQRHRjS5iJUmkkC+9daGrZVIqVkYWZWlKRhezVoVMM21bdVdNIYIuGERYYARKEUhCuVkP0/eN7OHM5647u81bP/gQj3O+//+c852LCQmEuanafElOWy2lGGvEVH2tqOrrRFVbLyhpGwRZ3zRX1jdLqVq6pBoZoqJvFWUtU1DTtguqkSUp+g5J0XYmqsYSRVFm0NyOQHHzH2hsajG7n/U1gZllcZnAkdOksYAGgsGg6e/paxZTFixPTk6eSTWmBaedJo0FNICgiY6u7lZB1VfGZALXm+SMiZABBE34urrbEhVtVVLS4llUKyzywhXbWnrfmyequ8zsogdmRp7H+om/Y3tw4g8nNxXeAIImvB3+JzizUU1sOebJTj9SN4qikSK35KHZ8/qTTYSHGkDQRJvX58PZjWgCkhdSsenC0/GW6liEM4CgidZ2b6cgGxsladFsmziOnE/+8fM3WwwMjpiNvjfmruLJJQlF4NUQ1YloAEETj9q9ATxHmHhWYeMcSPaFT4yM//rNDIyO/bTahr9+txnIKWoyf4xP8BqWgdIL5VGDGcjIq9/LJw0ZQGG+DWcBKah8amt/3D/I6ztCULUznAGPz4kB/4sPVnvprYCt/VRNgEvtDDjCz/EGhsMZwGlv631nxfOByarHNvptbkkzn9sRcH+U8QYmaNJw4BLQ6cfIzG+gn0ZFVI2LvIERmhQJLUFo9BX1/VPEMeKZAUHRL/MGpq2BfWWt1o7A6afbEONkjZ+kj44ka1c4A9F3AdYBgoVIv41nF8Abo4IZiHQO8AZw5KGzgK+DcOeAE2AXXGMGEJiFHN4AbjVacLgUtN3/cojmdoSoaNU2A5Mm3LkLnCDK+nWqb4EzkX6odoyK8YG3Ybg7IBbgkXKDajPm6UsP4L2fX9UZ13vACfBwvUl1GfDy3UM7uA0U4W2qyxAUYz/t4DZQA7VUlwHP6YO0g9uAgbtUlwEXxWHawW1ExbhHdRlwSh2lHdwGauA+1WUIinacdnAbuIzqqS5DkvUC2sFtoM4aqC4D/wGFo7IQ1qkYDozTkmKcxRcMFM55vMfB/SXYx+WSbFyF5aqEv1Xh0YqnGx4wuMdxm0HcwWqHp3gdTjn08+DIUVxS9d1U9//lLxSyAwwsz5wdAAAAAElFTkSuQmCC

library(dplyr)
library(gt)
library(here)
library(ggplot2)

analyse_data <- function(import_rds) {
  
  # Load the cleaned data
  cons_tum <- readRDS(file = import_rds)

  # Normalise for the heterozygous in tum relaive to constit
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
  
  
  
  # Create a gt table with color formatting
  # Apply conditional color formatting row by row
  cons_tum <- cons_tum %>%
    select(Pos., Gene.cons, c..HGVS.cons, VAF.cons, VAF.tum, LOH, `%tumoral`, chrom_num, pos_num) %>%
    # Sort by chrom_num first, then by pos_num
    arrange(chrom_num, pos_num) %>%
    # Remove temporary sorting columns
    select(-chrom_num, -pos_num)
  
    
  
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


