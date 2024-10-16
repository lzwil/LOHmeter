library(DT)  
library(fontawesome)  # Optional: for additional icons
library(here)
library(dplyr)
library(ggplot2)  # Make sure ggplot2 is loaded for plotting
library(shiny)
library(bsicons)
library(bslib)
source("import_data.R")
source("analyse_data.R")

# Define UI for application
ui <- navbarPage(
  collapsible = TRUE,
  title = tags$span(style = "font-size: 38px;", "LOHmeter"),
  tabPanel(
    "Main Page",
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
                style = "font-size: 60px;",
                textOutput(outputId = "mean_ui")
              ),
              showcase = tags$img(src = "test-tube.PNG", height = "150px")
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
  tabPanel("Page 2", 
           fluidRow(
             # Set the fluidRow to use flexbox to fill available vertical space
             style = "display: flex",  # Make the row take the full viewport height
             column(
               width = 4,  # First card will take 4/12 of the width
               style = "display: flex; height: 100%;",  # Make the column take full height
               card(
                 width = 12,  # Make the card take full width of the column
                 style = "flex: 1; overflow-y: auto; padding: 0;",  # Allow the card to grow and remove padding
                 h3("Welcome to Page 2"),
                 p("This is the content for the first card.")
               )
             ),
             column(
               width = 8,  # Second card will take 8/12 of the width
               style = "display: flex; flex-direction: column; height: 100%;",  # Make the column take full height
               card(
                 width = 12,  # Make the card take full width of the column
                 style = "flex: 1; overflow-y: auto; padding: 0;", 
                 DTOutput("table_uiTum") 
               )
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
              vjust = -0.5, hjust = 1, color = "#4D4D4D", size = 5, fontface = "bold") +  
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
    filtered_data <- result$cons_tum
    
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
  
  output$mean_ui <- renderText({
    req(processed_data())
    result <- processed_data()
    
    # Compute the mean of the '%tumoral' column
    mean_value <- result$cons_tum %>%
      summarise(Mean = mean(`%tumoral`, na.rm = TRUE)) %>%
      pull(Mean) # Extract the numeric value of the mean
    
    # Format the mean value as a string with appropriate text
    paste0(round(mean_value, 2), "%")  # Display as percentage
  })
  
  
  # Render the DataTable
  output$table_ui <- renderDT({
    req(processed_data())  # Ensure data is available
    result <- processed_data()
    df <- result$cons_tum
    
    
    # Filter the data if the box is checked
    filtered_data <- result$cons_tum
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
        order = list(0, 'asc')
      ),
      class = 'display nowrap compact stripe hover row-border order-column',
      escape = FALSE
    ))
  }
  )
  
  # Render the plot or placeholder
  output$plot_ui <- renderUI({
    result <- processed_data()
    
    if (is.null(result$boxplot)) {
      tags$div(
        style = "display: flex; justify-content: center; align-items: center; height: 100%;",  # Centering styles
        tags$img(
          src = "box-plot.png",  
          height = "150px",
          width = "150px",
          alt = "Loading plot placeholder"
        )
      )
      
    } else {
      # Show the plot when data is ready
      plotOutput(outputId = "plot")  
    }
  })
  
  # Render the plot
  output$plot <- renderPlot({
    result <- processed_data()
    req(result)  
    req(result$cons_tum)  
    
    # Filter the data by the selected gene
    
    if (any("Tous les locus" %in% input$selected_gene)){
      filtered_data <- result$cons_tum
    } else {
      filtered_data <- result$cons_tum %>%
        filter(Gene.cons %in% input$selected_gene)
    }
    # Generate the plot using the reusable function
    generate_boxplot(filtered_data)
  })
  
  
  # Update table and plot when rows are deleted
  observeEvent(input$delete_rows, {
    result <- processed_data()  # Get the current data
    req(result)  
    
    filtered_data <- result$cons_tum  
    
    if (input$filter_rows) {
      filtered_data <- filtered_data %>% filter(LOH %in% c("CIS", "TRANS"))
    }
    
    selected_rows <- input$table_ui_rows_selected  # Get selected rows from the filtered table
    
    if (length(selected_rows) > 0) {
      selected_data <- filtered_data[selected_rows, ]  
      
      original_data <- result$cons_tum  
      new_data <- original_data %>%
        filter(!Pos. %in% selected_data$Pos.)  
      
      processed_data(list(cons_tum = new_data, boxplot = generate_boxplot(new_data)))
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)