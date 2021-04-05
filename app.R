library(shiny)

ui <- fluidPage(
  
  # App title ----
  titlePanel("OncoSim"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("select", h3("Simulation Style"), 
                  choices = list("Original" = 1, "Evofreeze" = 2), 
                  selected = 1),
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "generations",
                  label = "Number of generations:",
                  min = 100,
                  max = 2000,
                  value = 1000),
      
      sliderInput(inputId = "metastasis",
                  label = "Generation where metastasis occurs:",
                  min = 200,
                  max = 1900,
                  value = 500),
      
      sliderInput(inputId = "metagenerations",
                  label = "Number of metastasis generations:",
                  min = 10,
                  max = 2000,
                  value = 30),
      
      actionButton("start", label = "Go")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "tree")
      
    )
  )
)


server <- function(input, output) {
  num_generations <- (eventReactive(input$start, {input$generations}))
  time_metastasis <- (eventReactive(input$start, {input$metastasis}))
  num_metastasis <- (eventReactive(input$start, {input$metagenerations}))
  
  output$tree <- renderPlot({

    barplot(c(num_generations, time_metastasis, num_metastasis))
    
  })
  
}


shinyApp(ui = ui, server = server)











