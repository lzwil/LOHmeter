
library(shiny)
library(bslib)
library(DT)
source("import_data.R")
source("analyse_data.R")

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "LOHmeter",
  sidebar = sidebar(
    fileInput(inputId = "constit",
              label = "Constit"),
    br(),
    fileInput(inputId = "tum",
              label = "Tumoral")
  ),
  card(
    DTOutput(outputId = "tableau")
    ),
  card(
    plotOutput(outputId = "plot"))
)


server <- function(input, output) {
  
  # Reactive expression to process and analyze data
  processed_data <- reactive({
    
    # Check if both files are uploaded
    req(input$constit, input$tum)
    
    # Import data
    import_data(constit = input$constit$datapath, tumoral = input$tum$datapath)
    
    # Ensure the RDS file exists
    req(file.exists("cons_tum_cleaned.rds"))
    
    # Analyze the data and return the result
    result <- analyse_data("cons_tum_cleaned.rds")
    
    # Return the result
    return(result)
  })
  
  # Render the table
  output$tableau <- renderDT({
    # Ensure data is available before rendering
    result <- processed_data()
    req(result$boxplot)  # Ensure the plot is available
    
    # Render the plot
    result$cons_tum
  })
  
  # Render the plot
  output$plot <- renderPlot({
    # Ensure data is available before rendering
    result <- processed_data()
    req(result$boxplot)  # Ensure the plot is available
    
    # Render the plot
    result$boxplot
  })
  
 
}

# Run the application 
shinyApp(ui = ui, server = server)
