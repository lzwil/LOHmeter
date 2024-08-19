# Analyse de la fréquence allélique 
# et déduction du pourcentage tumoral de l'échantillon
# Léo Zwilling
# APHM
# analyse_data.R

library(dplyr)
library(ggplot2)

cons_tum <- readRDS(file = "data/cons_tum_cleaned.rds")

# Color the LOH column depending on the class
# Create a gt table with color formatting
# Apply conditional color formatting row by row
cons_tum %>%
  gt() %>%
  tab_style(
    style = cell_fill(color = "red"),
    locations = cells_body(
      columns = vars(LOH, Coverage.tum, `%tumoral`),
      rows = LOH == "CIS"
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "#228B22"),
    locations = cells_body(
      columns = vars(LOH, Coverage.tum, `%tumoral`),
      rows = LOH == "TRANS"
    )
  )

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


# Create the boxplot
ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
  geom_boxplot(varwidth = T, fill="plum") +
  geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
  geom_text(data = summary_stats, aes(x = LOH, y = Mean, label =  paste(round(Mean, 2), "±", round(SD, 2))), 
            vjust = -0.5, hjust = 1, color = "#4D4D4D", size = 3.5) +  
  labs(title= "Box Plot",
       subtitle = "% Tumoral by LOH Classification",
       x = NULL,
       y = "% Tumoral") +
  theme_minimal()

# Mean of the two classes
mean_CIS_TRANS <- cons_tum %>%
  filter(LOH %in% c("CIS", "TRANS")) %>%
  summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
  pull(Mean)
print(mean_CIS_TRANS)



