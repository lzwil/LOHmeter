# Calculate `VAFtheoric` based on LOH classification
filtered_data <- result$cons_tum%>%
  mutate(
    VAFtheoric = case_when(
      LOH == "CIS" ~ (100 - `%tumoral`) / (200 - `%tumoral`),
      LOH == "TRANS" ~ 100 / (200 - `%tumoral`),
      TRUE ~ NA_real_
    )
  )



#Calculer avec le pourcentage tumoral théorique global
# filtered_data <- result$cons_tum %>%
#   mutate(
#     VAFtheoric = case_when(
#       LOH == "CIS" & !is.na(mean_tumor_percentage()) ~ (100 - mean_tumor_percentage()) / (200 - mean_tumor_percentage()),
#       LOH == "TRANS" & !is.na(mean_tumor_percentage()) ~ 100 / (200 - mean_tumor_percentage()),
#       TRUE ~ NA_real_
#     )
#   )