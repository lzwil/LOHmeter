library(DT)
library(dplyr)
library(ggplot2)
library(shiny)
library(bsicons)
library(bslib)
library(tidyr)

source("import_data.R")
source("analyse_data.R")

ui <- navbarPage(
  collapsible = TRUE,
  title = tags$span(style = "font-size: 38px;", "LOHmeter"),
  tabPanel(
    "Evaluation du pourcentage tumoral",
    page_sidebar(
      sidebar = sidebar(
        fileInput(inputId = "constit", label = "Constitutionel"),
        fileInput(inputId = "tum", label = "Tumoral"),
        checkboxInput(inputId = "filter_rows", label = "Afficher uniquement les lignes CIS et TRANS", value = TRUE),
        uiOutput("delete_button_ui"),
        uiOutput("gene_selector")
      ),
      fluidRow(
        fluidRow(
          column(
            width = 12,
            card(
              width = 12,
              style = "height: 500px; overflow-y: auto;",
              DTOutput("table_ui")
            )
          ),
          column(
            width = 4,
            value_box(
              title = "Pourcentage tumoral estimé",
              style = "height: 370px;",
              value = tags$div(
                style = "font-size: 50px;",
                textOutput(outputId = "mean_ui")
              ),
              showcase = tags$img(src = "test-tube.png", height = "130px")
            )
          ),
          column(
            width = 8,
            card(
              width = 12,
              style = "height: 370px",
              full_screen = TRUE,
              plotOutput(outputId = "plot")
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "LOH ou non ?",
    fluidRow(
      column(
        width = 5,
        h4("VAF estimée pour un nouveau variant avec LOH TRANS"),
        card(
          width = 12,
          style = "height: 450px",
          full_screen = TRUE,
          plotOutput(outputId = "conclu_plot")
        )
      ),
      column(
        width = 7,
        style = "display: flex; flex-direction: column; height: 100%;",
        card(
          width = 12,
          style = "flex: 1; overflow-y: auto; padding: 0;",
          DTOutput("table_uiTum")
        ),
        checkboxInput(inputId = "new_variants", label = "Nouveaux Variants Somatiques", value = TRUE)
      )
    )
  )
)

generate_boxplot <- function(data) {
  data_for_plot <- data %>%
    filter(!is.na(`%tumoral`), LOH %in% c("CIS", "TRANS")) %>%
    select(LOH, `%tumoral`)
  
  validate(need(nrow(data_for_plot) > 0, "Aucune donnée exploitable pour le graphique."))
  
  summary_stats <- data_for_plot %>%
    group_by(LOH) %>%
    summarise(
      Mean = mean(`%tumoral`, na.rm = TRUE),
      SD = sd(`%tumoral`, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
    geom_boxplot(varwidth = TRUE, outlier.shape = NA, linetype = 1) +
    geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
    geom_text(
      data = summary_stats,
      aes(x = LOH, y = Mean, label = paste("Mean:", round(Mean, 2), "±", round(SD, 2))),
      vjust = -2, hjust = 1.1, color = "#4D4D4D", size = 5, fontface = "bold"
    ) +
    labs(
      title = "Pourcentage estimé de cellules tumorales par classification LOH",
      x = NULL,
      y = "% Tumoral"
    ) +
    scale_fill_manual(values = c("CIS" = "#d4f1bc", "TRANS" = "#ffcccb")) +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 13),
      plot.title = element_text(size = 18)
    ) +
    coord_cartesian(ylim = c(0, 100))
}

generate_boxplotConclu <- function(data, selected_VAF) {
  data_for_plot <- data %>%
    filter(!is.na(VAFtheoTRANS) & !is.na(VAFtheoPASdeLOH)) %>%
    pivot_longer(
      cols = c(VAFtheoTRANS, VAFtheoPASdeLOH),
      names_to = "Category",
      values_to = "VAF"
    ) %>%
    mutate(Category = recode(Category,
                             VAFtheoTRANS = "LOH TRANS",
                             VAFtheoPASdeLOH = "PAS DE LOH"
    ))
  
  validate(need(nrow(data_for_plot) > 0, "Aucune donnée exploitable pour le graphique."))
  
  plot <- ggplot(data_for_plot, aes(x = Category, y = VAF)) +
    geom_boxplot(aes(fill = Category), color = "#4D4D4D", width = 0.4) +
    scale_fill_manual(values = c("LOH TRANS" = "#FFCCCB", "PAS DE LOH" = "#ADD8E6")) +
    labs(y = "VAF estimée", x = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0, 1))
  
  if (!is.null(selected_VAF) && is.numeric(selected_VAF) && length(selected_VAF) == 1) {
    plot <- plot +
      annotate("point", x = "LOH TRANS", y = selected_VAF, color = "red", size = 4, shape = 17) +
      annotate("text", x = "LOH TRANS", y = selected_VAF, label = round(selected_VAF, 2), vjust = -1, color = "red") +
      annotate("point", x = "PAS DE LOH", y = selected_VAF, color = "red", size = 4, shape = 17) +
      annotate("text", x = "PAS DE LOH", y = selected_VAF, label = round(selected_VAF, 2), vjust = -1, color = "red")
  }
  
  plot
}

server <- function(input, output, session) {
  processed_data <- reactiveVal(NULL)
  result_tumoral <- reactiveVal(NULL)
  selected_VAF <- reactiveVal(NULL)
  
  observeEvent(list(input$constit, input$tum), {
    req(input$constit, input$tum)
    
    import_data(
      constit = input$constit$datapath,
      tumoral = input$tum$datapath,
      output_cons_tum = "cons_tum_cleaned.rds",
      output_unique_tumoral = "unique_tumoral.rds"
    )
    
    req(file.exists("cons_tum_cleaned.rds"), file.exists("unique_tumoral.rds"))
    
    result <- analyse_data("cons_tum_cleaned.rds") %>%
      mutate(.row_id = row_number())
    
    processed_data(result)
    result_tumoral(readRDS(file = "unique_tumoral.rds"))
    selected_VAF(NULL)
  })
  
  filtered_processed_data <- reactive({
    req(processed_data())
    data <- processed_data()
    
    if (isTRUE(input$filter_rows)) {
      data <- data %>% filter(LOH %in% c("CIS", "TRANS"))
    }
    
    if (!is.null(input$selected_gene) && !("Tous les locus" %in% input$selected_gene)) {
      data <- data %>% filter(Gene %in% input$selected_gene)
    }
    
    data
  })
  
  tumor_table_data <- reactive({
    if (isTRUE(input$new_variants)) {
      req(result_tumoral())
      result_tumoral()
    } else {
      filtered_processed_data()
    }
  })
  
  output$delete_button_ui <- renderUI({
    req(processed_data())
    actionButton(inputId = "delete_rows", label = "Supprimer les lignes sélectionnées", icon = icon("trash-alt"))
  })
  
  output$gene_selector <- renderUI({
    req(processed_data())
    genes <- processed_data() %>%
      {if (isTRUE(input$filter_rows)) filter(., LOH %in% c("CIS", "TRANS")) else .} %>%
      pull(Gene) %>%
      unique() %>%
      sort()
    
    selectInput(
      "selected_gene",
      "Sélectionner un ou plusieurs locus:",
      choices = c("Tous les locus", genes),
      selected = "Tous les locus",
      multiple = TRUE
    )
  })
  
  mean_tumor_percentage <- reactive({
    data <- filtered_processed_data()
    validate(need(nrow(data) > 0, NA))
    data %>% summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>% pull(Mean)
  })
  
  output$mean_ui <- renderText({
    value <- mean_tumor_percentage()
    if (is.na(value)) "NA" else paste0(round(value, 2), "%")
  })
  
  selected_columns <- c("Pos.", "Gene", "c..HGVS", "VAF.cons", "VAF.tum", "LOH", "%tumoral")
  
  output$table_ui <- renderDT({
    data <- filtered_processed_data() %>% select(any_of(c(".row_id", selected_columns)))
    display_data <- data %>% select(-.row_id)
    
    datatable(
      display_data,
      options = list(
        pageLength = 50,
        lengthMenu = c(10, 25, 50, 100),
        autowidth = TRUE,
        scrollY = "400px",
        fixedHeader = TRUE,
        order = list(0, "asc")
      ),
      class = "display nowrap compact stripe hover row-border order-column",
      escape = FALSE,
      selection = "multiple"
    ) %>%
      formatStyle(
        "LOH",
        target = "row",
        backgroundColor = styleEqual(
          c("CIS", "TRANS"),
          c("#d4f1bc", "#ffcccb")
        )
      )
  })
  
  output$plot <- renderPlot({
    generate_boxplot(filtered_processed_data())
  })
  
  observeEvent(input$delete_rows, {
    req(processed_data())
    selected_rows <- input$table_ui_rows_selected
    if (length(selected_rows) == 0) return()
    
    current_filtered <- filtered_processed_data() %>% select(.row_id)
    ids_to_remove <- current_filtered$.row_id[selected_rows]
    
    updated <- processed_data() %>% filter(!(.row_id %in% ids_to_remove))
    processed_data(updated)
    selected_VAF(NULL)
  })
  
  output$table_uiTum <- renderDT({
    data <- tumor_table_data()
    
    dt <- datatable(
      data,
      options = list(
        pageLength = 50,
        lengthMenu = c(10, 25, 50, 100),
        autowidth = TRUE,
        scrollY = "400px",
        fixedHeader = TRUE,
        order = list(0, "asc")
      ),
      class = "display nowrap compact stripe hover row-border order-column",
      escape = FALSE,
      selection = "single"
    )
    
    if ("LOH" %in% names(data)) {
      dt <- dt %>%
        formatStyle(
          "LOH",
          target = "row",
          backgroundColor = styleEqual(
            c("CIS", "TRANS"),
            c("#d4f1bc", "#ffcccb")
          )
        )
    }
    
    dt
  })
  
  observeEvent(input$table_uiTum_rows_selected, {
    idx <- input$table_uiTum_rows_selected
    if (length(idx) != 1) {
      selected_VAF(NULL)
      return()
    }
    
    data <- tumor_table_data()
    vaf_value <- data[idx, "VAF.tum", drop = TRUE]
    if (!is.null(vaf_value) && is.numeric(vaf_value) && length(vaf_value) == 1) {
      selected_VAF(vaf_value)
    } else {
      selected_VAF(NULL)
    }
  })
  
  output$conclu_plot <- renderPlot({
    generate_boxplotConclu(filtered_processed_data(), selected_VAF())
  })
}

shinyApp(ui = ui, server = server)
