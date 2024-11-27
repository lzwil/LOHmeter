# Calculate `VAFtheoric` based on LOH classification
filtered_data <- result$cons_tum%>%
  mutate(
    VAFtheoric = case_when(
      LOH == "CIS" ~ (100 - `%tumoral`) / (200 - `%tumoral`),
      LOH == "TRANS" ~ 100 / (200 - `%tumoral`),
      TRUE ~ NA_real_
    )
  )



# # Calculer avec le pourcentage tumoral théorique global
#  filtered_data <- result$cons_tum %>%
#    mutate(
#      VAFtheoric = case_when(
#        LOH == "CIS" & !is.na(mean_tumor_percentage()) ~ (100 - mean_tumor_percentage()) / (200 - mean_tumor_percentage()),
#        LOH == "TRANS" & !is.na(mean_tumor_percentage()) ~ 100 / (200 - mean_tumor_percentage()),
#        TRUE ~ NA_real_
#      )
#    )


# # Reactive expression for the standard deviation of %tumoral for CIS rows
# sd_tumor_percentageCIS <- reactive({
#   req(processed_data())
#   cis_data <- processed_data()$cons_tum %>% filter(LOH == "CIS")
#   
#   # Ensure there are rows with LOH == "CIS"
#   req(nrow(cis_data) > 0)  # Proceed only if there are matching rows
#   
#   cis_data %>%
#     summarise(SD = sd(`%tumoral`, na.rm = TRUE)) %>%
#     pull(SD)
# })
# 
# # Reactive expression for the standard deviation of %tumoral for CIS rows
# sd_tumor_percentageTRANS <- reactive({
#   req(processed_data())
#   trans_data <- processed_data()$cons_tum %>% filter(LOH == "TRANS")
#   
#   # Ensure there are rows with LOH == "CIS"
#   req(nrow(trans_data) > 0)  # Proceed only if there are matching rows
#   
#   trans_data %>%
#     summarise(SD = sd(`%tumoral`, na.rm = TRUE)) %>%
#     pull(SD)
# })



# Identify variants from cons_tum with VAF.tum < 0.1
low_coverage_variants <- cons_tum %>%
  filter(VAF.cons < 0.1) %>%
  select(Pos., Gene.tum, c..HGVS.tum, VAF.tum) %>%
  rename(
    Gene = Gene.tum,
    c..HGVS = c..HGVS.tum,
    VAF = VAF.tum
  )

# Add low-coverage variants to unique_tumoral
unique_tumoral <- unique_tumoral %>%
  bind_rows(low_coverage_variants) %>%
  distinct() %>%  # Remove duplicates
  arrange(Pos.)



VAFtheo_CIS <- reactive({
  req(processed_data())
  
  VAFtheoCIS <- round((100 - mean_tumor_percentage()) / (200 - mean_tumor_percentage()), 2)
  return(VAFtheoCIS)
})

VAFtheo_TRANS <- reactive({
  req(processed_data())
  
  VAFtheoTRANS <- mean_tumor_percentage() / (mean_tumor_percentage() + 2 * (100 - mean_tumor_percentage()))
  return(VAFtheoTRANS)
})