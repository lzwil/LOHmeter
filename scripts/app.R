library(shiny)
library(bslib)
library(DT)
library(fontawesome)  # Optional: for additional icons
library(here)
library(dplyr)
library(ggplot2)  # Make sure ggplot2 is loaded for plotting
source("import_data.R")
source("analyse_data.R")

# Define UI for application
ui <- page_sidebar(
  title = "LOHmeter",
  sidebar = sidebar(
    fileInput(inputId = "constit", label = "Constitutionel"),
    br(),
    fileInput(inputId = "tum", label = "Tumoral"),
    checkboxInput(inputId = "filter_rows", label = "Afficher uniquement les lignes CIS et TRANS", value = TRUE),  
    br(),
    actionButton(inputId = "delete_rows", label = "Supprimer les lignes sélectionnées", icon = icon("trash-alt")),  
    br()
  ),
  # Main content layout with side-by-side cards
  fluidRow(
    # Column for the table UI
    column(
      width = 12,
      card(
        width = 12,
        style = "height: 500px; overflow-y: auto;",
        DTOutput("table_ui")  # Changed from uiOutput to DTOutput
      )
    ),
    # Column for the locus selection
    column(
      width = 4,  # 1/3 of the row
      card(
        width = 12,
        style = "height: 370px;", 
        full_screen = FALSE
      )
    ),
    # Column for the plot content
    column(
      width = 8,  # 2/3 of the row
      card(
        width = 12,
        style = "height: 370px;", 
        full_screen = TRUE,
        uiOutput(outputId = "plot_ui")
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
    geom_boxplot(varwidth = TRUE) +  
    geom_point(data = summary_stats, aes(x = LOH, y = Mean), color = "#4D4D4D", size = 3, shape = 20, show.legend = FALSE) +
    geom_text(data = summary_stats, aes(x = LOH, y = Mean, label = paste("Mean:", round(Mean, 2), "±", round(SD, 2))), 
              vjust = -0.5, hjust = 1, color = "#4D4D4D", size = 3.5) +  
    labs(title = "Pourcentage estimé de cellules tumorales par classification LOH",
         x = NULL,
         y = "% Tumoral") +
    scale_fill_manual(values = c("CIS" = "#d4f1bc", "TRANS" = "#ffcccb")) +
    theme_minimal()
}

server <- function(input, output) {
  # Initialize a reactive value to hold the processed data
  processed_data <- reactiveVal(NULL)  
  
  # Load and process data when files are uploaded
  observeEvent(c(input$constit, input$tum), {
    req(input$constit, input$tum)  # Ensure both files are uploaded
    
    # Import and process data
    import_data(constit = input$constit$datapath, tumoral = input$tum$datapath)
    req(file.exists("cons_tum_cleaned.rds"))
    
    # Analyze the data and set it in reactive value
    result <- analyse_data("cons_tum_cleaned.rds")
    processed_data(result)  # Store the result in the reactive value
  })
  
  # Render the DataTable
  output$table_ui <- renderDT({
    req(processed_data())  # Ensure data is available
    result <- processed_data()
    
    # Filter the data if the box is checked
    filtered_data <- result$cons_tum
    if (input$filter_rows) {
      filtered_data <- filtered_data %>%
        filter(LOH %in% c("CIS", "TRANS"))
    }
    
    datatable(
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
    )
  })
  
  # Render the plot or placeholder
  output$plot_ui <- renderUI({
    result <- processed_data()
    
    if (is.null(result$boxplot)) {
      # Show a placeholder image when the plot is not ready
      tags$img(
        src = "test.txt",  
        style = "max-width: 100%; max-height: 100%; display: block; margin: auto;",
        alt = "Loading plot placeholder"
      )
    } else {
      # Show the plot when data is ready
      plotOutput(outputId = "plot", height = "100%")  
    }
  })
  
  # Render the plot
  output$plot <- renderPlot({
    result <- processed_data()
    req(result)  
    req(result$cons_tum)  
    
    # Generate the plot using the reusable function
    generate_boxplot(result$cons_tum)
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
