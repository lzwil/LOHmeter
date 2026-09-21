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
        fileInput(
          inputId = "constit",
          label = "Constitutionnel",
          buttonLabel = "Parcourir...",
          placeholder = "Aucun fichier sélectionné",
          width = "100%"
        ),
        uiOutput("constit_filename"),
        fileInput(
          inputId = "tum",
          label = "Tumoral",
          buttonLabel = "Parcourir...",
          placeholder = "Aucun fichier sélectionné",
          width = "100%"
        ),
        uiOutput("tum_filename"),
        checkboxInput(inputId = "filter_rows", label = "Afficher uniquement les lignes CIS et TRANS", value = FALSE),
        uiOutput("delete_button_ui"),
        uiOutput("gene_selector")
      ),
      fluidRow(
        uiOutput("main_content_ui")
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
          plotOutput(outputId = "conclu_plot", height = "500px")
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
  ),
  tabPanel(
    "Comment utiliser l'outil ?",
    fluidPage(
      br(),
      fluidRow(
        column(
          width = 12,
          tags$div(
            style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            tags$h3("Guide d'utilisation", style = "margin: 0;"),
            tags$a(
              href = "guide_utilisation.pdf",
              target = "_blank",
              "Ouvrir le PDF dans un nouvel onglet"
            )
          )
        )
      ),
      tags$div(
        style = "height: calc(100vh - 190px); width: 100%;",
        tags$iframe(
          src = "guide_utilisation.pdf",
          style = "width: 100%; height: 100%; border: 1px solid #ddd; border-radius: 8px;"
        )
      )
    )
  )
)


theme_lohmeter <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0, color = "#1f2937"),
      plot.subtitle = element_text(size = 11.5, hjust = 0, color = "#6b7280"),
      axis.title = element_text(face = "bold", color = "#374151"),
      axis.text = element_text(color = "#374151"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(12, 18, 12, 12)
    )
}

generate_boxplot <- function(data) {
  data_for_plot <- data %>%
    filter(!is.na(`%tumoral`), LOH %in% c("CIS", "TRANS")) %>%
    mutate(LOH = factor(LOH, levels = c("CIS", "TRANS")))
  
  validate(need(nrow(data_for_plot) > 0, "Aucune donnée exploitable pour le graphique."))
  
  summary_stats <- data_for_plot %>%
    group_by(LOH) %>%
    summarise(
      n = n(),
      mean_value = mean(`%tumoral`, na.rm = TRUE),
      median_value = median(`%tumoral`, na.rm = TRUE),
      sd_value = sd(`%tumoral`, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
    geom_boxplot(
      width = 0.52,
      alpha = 0.85,
      outlier.shape = NA,
      color = "#4b5563"
    ) +
    geom_jitter(
      width = 0.10,
      alpha = 0.55,
      size = 2.2,
      color = "#374151"
    ) +
    geom_text(
      data = summary_stats,
      aes(x = LOH, y = 98, label = paste0("n = ", n)),
      inherit.aes = FALSE,
      size = 4.2,
      fontface = "bold",
      color = "#374151"
    ) +
    scale_fill_manual(values = c("CIS" = "#d9f2c7", "TRANS" = "#ffd6d6")) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(
      title = "Estimation du pourcentage tumoral",
      subtitle = "Distribution des variants classés CIS et TRANS",
      x = NULL,
      y = "% tumoral estimé"
    ) +
    theme_lohmeter() +
    theme(
      axis.text.x = element_text(face = "bold", size = 12, color = "#1f2937")
    )
}


generate_boxplotConclu <- function(data, selected_VAF) {
  data_for_plot <- data %>%
    filter(!is.na(VAFtheoTRANS), !is.na(VAFtheoPASdeLOH)) %>%
    pivot_longer(
      cols = c(VAFtheoTRANS, VAFtheoPASdeLOH),
      names_to = "Category",
      values_to = "VAF"
    ) %>%
    mutate(
      Category = recode(
        Category,
        VAFtheoTRANS = "LOH TRANS",
        VAFtheoPASdeLOH = "PAS DE LOH"
      ),
      Category = factor(Category, levels = c("LOH TRANS", "PAS DE LOH"))
    )
  
  validate(need(nrow(data_for_plot) > 0, "Aucune donnée exploitable pour le graphique."))
  
  summary_stats <- data_for_plot %>%
    group_by(Category) %>%
    summarise(
      n = n(),
      median_value = median(VAF, na.rm = TRUE),
      q1 = quantile(VAF, 0.25, na.rm = TRUE),
      q3 = quantile(VAF, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  subtitle_text <- if (!is.null(selected_VAF) && is.numeric(selected_VAF) && length(selected_VAF) == 1) {
    paste0("Variant sélectionné : VAF observée = ", round(selected_VAF, 3))
  } else {
    "Sélectionner un variant dans le tableau pour l'ajouter au graphique"
  }
  
  plot <- ggplot(data_for_plot, aes(x = Category, y = VAF, fill = Category)) +
    geom_boxplot(
      width = 0.50,
      alpha = 0.90,
      outlier.shape = NA,
      color = "#4b5563"
    ) +
    geom_jitter(
      width = 0.08,
      alpha = 0.28,
      size = 1.8,
      color = "#4b5563"
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 8,
      color = "#111827"
    ) +
    geom_text(
      data = summary_stats,
      aes(x = Category, y = 0.98, label = paste0("n = ", n)),
      inherit.aes = FALSE,
      size = 4.0,
      fontface = "bold",
      color = "#374151"
    ) +
    scale_fill_manual(values = c("LOH TRANS" = "#ffd6d6", "PAS DE LOH" = "#d9ecff")) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(
      title = "Comparaison à la VAF théorique",
      subtitle = subtitle_text,
      x = NULL,
      y = "VAF estimée"
    ) +
    theme_lohmeter()
  
  if (!is.null(selected_VAF) && is.numeric(selected_VAF) && length(selected_VAF) == 1) {
    plot <- plot +
      geom_hline(
        yintercept = selected_VAF,
        linetype = "dashed",
        linewidth = 0.7,
        color = "#dc2626"
      ) +
      annotate(
        "point",
        x = 1,
        y = selected_VAF,
        color = "#dc2626",
        size = 4,
        shape = 18
      ) +
      annotate(
        "point",
        x = 2,
        y = selected_VAF,
        color = "#dc2626",
        size = 4,
        shape = 18
      ) +
      annotate(
        "text",
        x = 1.5,
        y = min(0.99, selected_VAF + 0.05),
        label = paste0("VAF observée = ", round(selected_VAF, 3)),
        color = "#dc2626",
        fontface = "bold",
        size = 4.2
      )
  }
  
  plot
}

server <- function(input, output, session) {
  processed_data <- reactiveVal(NULL)
  result_tumoral <- reactiveVal(NULL)
  selected_VAF <- reactiveVal(NULL)
  import_error <- reactiveVal(NULL)
  
  required_columns_display <- c(
    "Gene", "Pos.", "Coverage", "c. HGVS",
    "Transcript", "Type", "Nuc Change", "AA Change", "p. HGVS"
  )
  
  observeEvent(list(input$constit, input$tum), {
    req(input$constit, input$tum)
    
    import_attempt <- tryCatch({
      import_data(
        constit = input$constit$datapath,
        tumoral = input$tum$datapath,
        output_cons_tum = "cons_tum_cleaned.rds",
        output_unique_tumoral = "unique_tumoral.rds"
      )
      NULL
    }, error = function(e) conditionMessage(e))
    
    if (!is.null(import_attempt)) {
      import_error(import_attempt)
      processed_data(NULL)
      result_tumoral(NULL)
      return()
    }
    
    req(file.exists("cons_tum_cleaned.rds"), file.exists("unique_tumoral.rds"))
    
    analyse_attempt <- tryCatch({
      analyse_data("cons_tum_cleaned.rds") %>% mutate(.row_id = row_number())
    }, error = function(e) conditionMessage(e))
    
    if (is.character(analyse_attempt)) {
      import_error(analyse_attempt)
      processed_data(NULL)
      result_tumoral(NULL)
      return()
    }
    
    import_error(NULL)
    processed_data(analyse_attempt)
    result_tumoral(readRDS(file = "unique_tumoral.rds"))
    selected_VAF(NULL)
  })
  
  output$main_content_ui <- renderUI({
    if (!is.null(import_error())) {
      div(
        style = "height: 500px; display: flex; align-items: center; justify-content: center; text-align: center;",
        div(
          style = "max-width: 650px; padding: 30px;",
          tags$h4("Vérifier le format des données d'entrées.", style = "color: #b91c1c; margin-bottom: 12px;"),
          tags$p(
            style = "font-size: 15px;",
            paste("Colonnes minimales :", paste(required_columns_display, collapse = ", "))
          ),
          tags$p(style = "color: #888; font-size: 13px; margin-top: 15px;", import_error())
        )
      )
    } else {
      fluidRow(
        column(
          width = 12,
          card(
            width = 12,
            style = "height: 500px; overflow-y: auto;",
            uiOutput("table_ui_wrapper")
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
            plotOutput(outputId = "plot", height = "400px")
          )
        )
      )
    }
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
  
  output$constit_filename <- renderUI({
    req(input$constit)
    tags$div(
      style = paste(
        "margin-top: -12px; margin-bottom: 12px;",
        "padding: 8px 10px;",
        "font-size: 13px;",
        "line-height: 1.35;",
        "color: #374151;",
        "background: #f9fafb;",
        "border: 1px solid #e5e7eb;",
        "border-radius: 6px;",
        "white-space: normal;",
        "overflow-wrap: anywhere;",
        "word-break: break-word;"
      ),
      tags$strong("Fichier : "),
      input$constit$name
    )
  })

  output$tum_filename <- renderUI({
    req(input$tum)
    tags$div(
      style = paste(
        "margin-top: -12px; margin-bottom: 12px;",
        "padding: 8px 10px;",
        "font-size: 13px;",
        "line-height: 1.35;",
        "color: #374151;",
        "background: #f9fafb;",
        "border: 1px solid #e5e7eb;",
        "border-radius: 6px;",
        "white-space: normal;",
        "overflow-wrap: anywhere;",
        "word-break: break-word;"
      ),
      tags$strong("Fichier : "),
      input$tum$name
    )
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
    if (nrow(data) == 0) {
      return(NA_real_)
    }
    data %>% summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>% pull(Mean)
  })
  
  output$mean_ui <- renderText({
    value <- mean_tumor_percentage()
    if (is.na(value)) "NA" else paste0(round(value, 2), "%")
  })
  
  selected_columns <- c("Pos.", "Gene", "c..HGVS", "VAF.cons", "VAF.tum", "LOH", "%tumoral")
  
  output$table_ui_wrapper <- renderUI({
    if (nrow(filtered_processed_data()) == 0) {
      div(
        style = "height: 460px; display: flex; align-items: center; justify-content: center; text-align: center; color: #888; font-size: 17px; padding: 20px;",
        "Aucune variation de VAF détectée attestant de perte d'hétérozygotie (LOH)"
      )
    } else {
      DTOutput("table_ui")
    }
  })
  
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
