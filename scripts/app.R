library(DT)  
library(fontawesome)  # Optional: for additional icons
library(here)
library(dplyr)
library(ggplot2)  # Make sure ggplot2 is loaded for plotting
library(shiny)
library(bsicons)
library(bslib)
library(tidyr)
source("import_data.R")
source("analyse_data.R")

# Define UI for application
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
          # Column for the table 
          column(
            width = 12,
            card(
              width = 12,
              style = "height: 500px; overflow-y: auto;",
              DTOutput("table_ui") 
            )
          ),
          # Column for the %tumoral
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
          # Column for the plot content
          column(
            width = 8,  # 2/3 of the row
            card(
              width = 12,
              style = "height: 370px",
              full_screen = TRUE,
              uiOutput(outputId = "plot_ui")
            )
          )  
        )
      )
    )
  ),
  tabPanel("LOH ou non ?", 
           fluidRow(
             column(
               width = 5,
               h4("VAF estimée pour un nouveau variant avec LOH TRANS"),
               card(
                 width = 12,
                 style = "height: 450px",
                 full_screen = TRUE,
                 uiOutput(outputId = "concluPlot_ui")
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


# Function to generate the box plot
generate_boxplot <- function(data) {
  data_for_plot <- data %>%
    filter(!is.na(`%tumoral`)) %>%
    select(LOH, `%tumoral`)
  
  summary_stats <- data_for_plot %>%
    group_by(LOH) %>%
    summarise(
      Mean = mean(`%tumoral`, na.rm = TRUE),
      SD = sd(`%tumoral`, na.rm = TRUE)
    )
  
  ggplot(data_for_plot, aes(x = LOH, y = `%tumoral`, fill = LOH)) +
    geom_boxplot(varwidth = TRUE, outlier.shape = NA, linetype = 1) +  
    geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
    geom_text(data = summary_stats, aes(x = LOH, y = Mean, label = paste("Mean:", round(Mean, 2), "±", round(SD, 2))), 
              vjust = -2, hjust = 1.1, color = "#4D4D4D", size = 5, fontface = "bold") +  
    labs(title = "Pourcentage estimé de cellules tumorales par classification LOH",
         x = NULL,
         y = "% Tumoral") +
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
    ylim(0, 100) 
}

generate_boxplotConclu <- function(data, selected_VAF) {
  # Reshape the data for a combined box plot
  data_for_plot <- data %>%
    filter(!is.na(VAFtheoTRANS) & !is.na(VAFtheoPASdeLOH)) %>%
    pivot_longer(cols = c(VAFtheoTRANS, VAFtheoPASdeLOH),
                 names_to = "Category", values_to = "VAF") %>%
    mutate(Category = recode(Category, 
                             VAFtheoTRANS = "LOH TRANS",
                             VAFtheoPASdeLOH = "PAS de LOH"))  # Rename for clarity
  
  # Initialize the plot
  plot <- ggplot(data_for_plot, aes(x = Category, y = VAF)) +
    geom_boxplot(aes(fill = Category), color = "#4D4D4D", width = 0.4) +
    scale_fill_manual(values = c("LOH TRANS" = "#FFCCCB", "PAS de LOH" = "#ADD8E6")) + # Custom colors
    labs(title = "VAF estimée pour LOH TRANS et PAS de LOH",
         y = "VAF estimée", x = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16),
      legend.position = "none"  # Remove legend as x-axis labels are descriptive
    ) +
    ylim(0, 1)
  
  # Conditionally add the point and text layers for VAF if selected_VAF is not NULL
  if (!is.null(selected_VAF)) {
    plot <- plot +
      annotate("point", x = "LOH TRANS", y = selected_VAF, color = "red", size = 4, shape = 17) +
      annotate("text", x = "LOH TRANS", y = selected_VAF, label = round(selected_VAF, 2), vjust = -1, color = "red")
  }
  
  # Return the final plot
  return(plot)
}





server <- function(input, output) {
  
  
  # Initialize a reactive value to hold the processed data
  processed_data <- reactiveVal(NULL)  
  result_tumoral <- reactiveVal(NULL)
  
  # Load and process data when files are uploaded
  observeEvent(c(input$constit, input$tum), {
    req(input$constit, input$tum)  # Ensure both files are uploaded
    
    # Import and process data
    import_data(constit = input$constit$datapath, tumoral = input$tum$datapath)
    req(file.exists("cons_tum_cleaned.rds"))
    
    # Analyze the data 
    result <- analyse_data("cons_tum_cleaned.rds")
    #print(result)
    processed_data(result)  # Store the result in the reactive value
    
    # Store the unique_tumoral data in reactive value
    req(file.exists("unique_tumoral.rds"))
    result_tumoral(readRDS(file = "unique_tumoral.rds"))
  })
  
  # Delete button
  output$delete_button_ui <- renderUI({
    req(input$constit, input$tum)  # Ensure both files are uploaded
    
    actionButton(inputId = "delete_rows", label = "Supprimer les lignes sélectionnées", icon = icon("trash-alt"))
  })
  
  # Render the gene selection UI based on the processed data
  output$gene_selector <- renderUI({
    req(processed_data())
    
    # Get the processed data
    result <- processed_data()
    filtered_data <- result
    
    # Filter the data if the "Afficher uniquement les lignes CIS et TRANS" checkbox is checked
    if (input$filter_rows) {
      filtered_data <- filtered_data %>% filter(LOH %in% c("CIS", "TRANS"))
    }
    
    genes <- unique(filtered_data$Gene.cons)
    genes <- c("Tous les locus", genes)  # Add "All Genes" as the first option
    selectInput("selected_gene", 
                "Sélectionner un ou plusieurs locus:", 
                choices = genes, 
                selected = "Tous les locus", 
                multiple = TRUE
    )
  })
  
  
  # Reactive expression for the tumor percentage
  mean_tumor_percentage <- reactive({
    req(processed_data())
    # Get the current data
    result <- processed_data()
    data <- result
    
    # Filter the data by the selected gene(s)
    if (!("Tous les locus" %in% input$selected_gene)) {
      data <- data %>% filter(Gene.cons %in% input$selected_gene)
    }
    
    # Calculate the mean tumor percentage
    data %>%
      summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
      pull(Mean)
  })
  
  # Reactive expression for the tumor percentage
  mean_tumor_percentageTRANS <- reactive({
    req(processed_data())
    mean_TRANS_data <- processed_data()$cons_tum %>% filter(LOH == "TRANS")
    
    mean_TRANS_data %>%
      summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
      pull(Mean)
  })
  
  # Reactive expression for the tumor percentage
  mean_tumor_percentageCIS <- reactive({
    req(processed_data())
    mean_cis_data <- processed_data()$cons_tum %>% filter(LOH == "CIS")
    
    mean_cis_data %>%
      summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
      pull(Mean)
  })
  
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
  
  output$concluPlot_ui <- renderUI({
    req(processed_data())
    result <- processed_data()
    
    if (is.null(result)) {
      return(tags$div(
        style = "display: flex; justify-content: center; align-items: center; height: 100%;",
        tags$img(
          src = "box-plot.png",  
          height = "150px",
          width = "150px",
          alt = "Loading plot placeholder"
        )
      ))
    } else {
      # Show the plot when data is ready
      return(plotOutput(outputId = "conclu_plot"))
    }
  })
  
  
  selected_tumoral_data <- reactive({
    req(result_tumoral())  # Ensure tumor data is available
    selected_row <- input$table_uiTum_rows_selected  # Get selected row index
    if (length(selected_row) > 0) {
      result_tumoral()[selected_row, ]  # Extract the data for the selected row
    } else {
      NULL  # Return NULL if no row is selected
    }
  })

  # Render the plot
  output$conclu_plot <- renderPlot({
      req(processed_data(), selected_tumoral_data(), input$selected_gene)  # Ensure both data and selected row are available
      result <- processed_data()
      selected_data <- selected_tumoral_data()
      
      
      # Filter the data by the selected gene(s)
      if (!("Tous les locus" %in% input$selected_gene)) {
        result <- result %>%
          filter(Gene.cons %in% input$selected_gene)
      }
      
      # Pass the selected VAF value from the selected row to the plot
      selected_VAF <- if (!is.null(selected_data)) {
        selected_data$VAF.tum
      } else {
        NULL  # If no row is selected, no point will be shown
      }
      
      # Generate the plot with the reusable function
      
      generate_boxplotConclu(result, selected_VAF)  
    })
  

  output$mean_ui <- renderText({
    paste0(round(mean_tumor_percentage(), 2), "%")  # Display as percentage
  })
  
  # Render the DataTable
  output$table_ui <- renderDT({
    req(processed_data())  # Ensure data is available
    result <- processed_data()
    df <- result
    filtered_data <- result
    
    # Filter the data if the box is checked
    if (input$filter_rows) {
      filtered_data <- filtered_data %>%
        filter(LOH %in% c("CIS", "TRANS"))
    }
    
    return(datatable(
      filtered_data,
      options = list(
        pageLength = 50,
        lengthMenu = c(10, 25, 50, 100),
        autowidth = TRUE,
        scrollY = "400px",
        fixedHeader = TRUE,
        order = list(0, 'asc'),  
        rowCallback = JS(
          "function(row, data, index) {",
          "  if (data[6] == 'CIS') {",  
          "    $('td', row).css('background-color', '#d4f1bc');",  
          "  } else if (data[6] == 'TRANS') {",  
          "    $('td', row).css('background-color', '#ffcccb');",  
          "  }",
          "}"
        )
      ),
      class = 'display nowrap compact stripe hover row-border order-column',
      escape = FALSE
    ))
  }
  )

  # Render the DataTable
  observe({
    if (input$new_variants) {
      output$table_uiTum <- renderDT({
        req(result_tumoral())  # Ensure data is available
        return(datatable(
          result_tumoral(),
          options = list(
            pageLength = 50,
            lengthMenu = c(10, 25, 50, 100),
            autowidth = TRUE,
            scrollY = "400px",
            fixedHeader = TRUE,
            order = list(0, 'asc'),
            searching = FALSE,
            lengthChange = FALSE,
            selection = 'single'
          ),
          class = 'display nowrap compact stripe hover row-border order-column',
          escape = FALSE
        ))
      })
    } else {
      output$table_uiTum <- renderDT({
        req(processed_data())  # Ensure data is available
        result <- processed_data()
        df <- result
        filtered_data <- result
        if (input$filter_rows) {
          filtered_data <- filtered_data %>%
            filter(LOH %in% c("CIS", "TRANS"))
        }
        return(datatable(
          filtered_data,
          options = list(
            pageLength = 50,
            lengthMenu = c(10, 25, 50, 100),
            autowidth = TRUE,
            scrollY = "400px",
            fixedHeader = TRUE,
            order = list(0, 'asc'),
            rowCallback = JS(
              "function(row, data, index) {",
              "  if (data[6] == 'CIS') {",
              "    $('td', row).css('background-color', '#d4f1bc');",
              "  } else if (data[6] == 'TRANS') {",
              "    $('td', row).css('background-color', '#ffcccb');",
              "  }",
              "}"
            )
          ),
          class = 'display nowrap compact stripe hover row-border order-column',
          escape = FALSE
        ))
      })
    }
  })
  
  
  # Render the plot or placeholder
  output$plot_ui <- renderUI({
    req(processed_data())
    result <- processed_data()
    
    if (is.null(result)) {
      return(tags$div(
        style = "display: flex; justify-content: center; align-items: center; height: 100%;",  # Centering styles
        tags$img(
          src = "box-plot.png",  
          height = "150px",
          width = "150px",
          alt = "Loading plot placeholder"
        )
      ))
      
    } else {
      # Show the plot when data is ready
      return(plotOutput(outputId = "plot"))  
      
      
    }
  })
  
  # Render the plot
  output$plot <- renderPlot({
    result <- processed_data()
    req(result)  
    
    # Filter the data by the selected gene
    
    if (any("Tous les locus" %in% input$selected_gene)){
      filtered_data <- result
    } else {
      
      filtered_data <- result %>%
        filter(Gene.cons %in% input$selected_gene)
    }
    # Generate the plot using the reusable function
    generate_boxplot(filtered_data)
  })
  
  
  # Update table and plot when rows are deleted
  observeEvent(input$delete_rows, {
    result <- processed_data()  # Get the current data
    req(result)  
    
    filtered_data <- result
    
    if (input$filter_rows) {
      filtered_data <- filtered_data %>% filter(LOH %in% c("CIS", "TRANS"))
    }
    
    selected_rows <- input$table_ui_rows_selected  # Get selected rows from the filtered table
    
    if (length(selected_rows) > 0) {
      selected_data <- filtered_data[selected_rows, ]  
      
      original_data <- result
      new_data <- original_data %>%
        filter(!Pos. %in% selected_data$Pos.)  
      
      processed_data(new_data)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)