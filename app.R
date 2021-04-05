library(shiny)
source("helpers.R")


ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
      )
    )
  ),
  # App title ----
  titlePanel("OncoSim"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("select", h4("Simulation Style"), 
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
                  min = 100,
                  max = 2000,
                  value = 500),
      
      sliderInput(inputId = "metagenerations",
                  label = "Number of metastasis generations:",
                  min = 100,
                  max = 1000,
                  value = 500),
      
      sliderInput(inputId = "purity",
                  label = "Sample purity:",
                  min = 0.1,
                  max = 1,
                  value = 0.9),
      
      actionButton("start", label = "Go")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      p("Oncosim is a clonal evolution tree simulator that allows the CCF of 
        clones and VAFs of mutations on human autosomal chromosomes to be
        modelled across primary and metastatic samples. Set your parameters in
        the sidebar, click go, and wait for the simulation to be complete."),
      br(),
      br(),
      plotOutput(outputId = "chroms"),
      br(),
      p("This is a visualisation of where on the genome the single nucleotide 
      variants have occurred during the simulation. Chromosomes 1-22 are shown
        from left to right with the paternal and maternal copies coloured blue
        and pink respectively. The locations of SNVs are marked by red Xs"),
      br(),
      br(),

      plotOutput(outputId = "tree"),
      br(),
      p("This is the clonal evolution tree of the mutations produced during the
        simulation. The cancer spawns from a grey 'normal' cell, with nodes coloured
        according to their cancer cell fraction, and labelled with their cluster 
        genotype. The metastatic tree is connected to its seed in the primary tree
        with a dotted line"),
      br(),
      br(),

      plotOutput(outputId = "ccf"),
      br(),
      p("This is a scatterplot of the cancer cell fractions of each cluster in
        the primary and metastatic samples. Shared clusters will have two positive
        coordinates, whereas clusters unique to either the primary or metastatic
        tumours will line the x and y axis."),
      br(),
      br(),
      
      plotOutput(outputId = "vaf"),
      br(),
      p("This is a scatterplot of the cancer cell fractions of each cluster in
        the primary and metastatic samples. Shared clusters will have two positive
        coordinates, whereas clusters unique to either the primary or metastatic
        tumours will line the x and y axis."),
    )
  )
)


server <- function(input, output) {
  observeEvent(input$start, {
    if(input$select == 1) {
    withProgress(message = "Running Simulation", value = 0, {
      chromplot <- bigcancersimulator(input$generations, input$metastasis, 
                                      input$metagenerations, input$purity)
    })
    }
    else {
      withProgress(message = "Running Simulation", value = 0, {
        chromplot <- evofreezesimulator(input$generations, input$metastasis, 
                                        input$metagenerations, input$purity)
    })
    }  
    output$chroms <- renderPlot({
      chromplot[[1]]
    })
    output$tree <- renderPlot({
      chromplot[[2]]
    })
    output$ccf <- renderPlot({
      chromplot[[3]]
    })
    output$vaf <- renderPlot({
      chromplot[[4]]
    })
  })


  output$tree <- renderPlot({
    
  })
  
}


shinyApp(ui = ui, server = server)











