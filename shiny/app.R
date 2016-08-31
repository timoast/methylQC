library(shiny)
library(ggplot2)
library(methylQC)
library(plotly)


ui <- fluidPage(theme = "theme.css",
                titlePanel("methylQC"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('data', 'Choose a file'),
                    checkboxInput(inputId = "demo", label = "Use demo data", value = FALSE),
                    tags$hr(),
                    uiOutput("chromSelect"),
                    numericInput(inputId = "position",
                                 label = "Position", min = 0, value = 5000, step = 2000),
                    sliderInput(inputId = "zoom", label = "Zoom", min = 1, max = 100, value = 15),
                    tags$hr(),
                    radioButtons('format', 'Document format', c('PDF', 'HTML'),
                                 inline = TRUE),
                    textInput("filename", label = ("File name"), placeholder = "Enter a file name",
                              value = "methylQC-report"),
                    downloadButton("downloadReport", "Generate report")
                  ),
                  mainPanel(
                    uiOutput("tabs")
                )
          )
)

server <- function(input, output) {
  # change the max file size to 1 GB
  options(shiny.maxRequestSize=1000*1024^2)
  
  # load the data. If user selects demo data, load from the demo file path
  dat <- reactive({
    file1 <- input$data
    if(input$demo == TRUE) {
      loadData("../data/example_data.CGmap.gz")
    } else if(is.null(file1)){
      return()
    } else {
      if(tools::file_ext(file1$name) == "gz"){
        zipped = TRUE
      } else {
        zipped = FALSE
      }
      loadData(file1$datapath, forceZipped = zipped)
    }
  })
  
  rv <- reactiveValues(data = observe(dat()), start = 2000, end = 5000)
  
  # Set chromosome options if input is specified
  output$chromSelect <- renderUI({
    if (is.null(dat())){return()}
    selectInput("chrom", "Chr", choices = unique(dat()$chr), selected = unique(dat()$chr)[[2]])
  })
  
  # Change dataset to user-specified chromosome
  observeEvent(input$chrom, {rv$data <- subset(dat(), dat()$chr == input$chrom)})
  
  # Change start and stop positions for browser based on user-unput position
  observeEvent(input$position, {rv$start <- input$position - (300 * input$zoom)
                                rv$end <- input$position + (300 * input$zoom)})
  
  # Change the dataset to user-specified positions based on zoom level
  observeEvent(input$zoom, {rv$start <- input$position - (300 * input$zoom)
                            rv$end <- input$position + (300 * input$zoom)})
  
  # If data is loaded show the about, plots, statistics, and data tabs
  # If data is not loaded, only show the about tab
  output$tabs <- renderUI({
    if (is.null(dat())){
      tabsetPanel(type = "tabs",
                  tabPanel("About", includeMarkdown("./text/about.md")))
    } else {
      tabsetPanel(type = "tabs",
                  tabPanel("About", includeMarkdown("./text/about.md")),
                  tabPanel("Plots",
                           fluidRow(plotOutput("browser")),
                           fluidRow(plotlyOutput("survival")),
                           fluidRow(
                             column(6,
                                    plotOutput("hist")
                             ),
                             column(6, plotOutput("lambdaCov"),
                                    tableOutput("conversion")
                             )
                           )
                  ),
                  tabPanel("Statistics", h3("Sequencing coverage"),
                           includeMarkdown("./text/statistics_coverage.md"),
                           tableOutput("stats"),
                           h3("Non-conversion rates"),
                           tableOutput("nonconversion")),
                  tabPanel("Data", tableOutput("headData"))
                  )
    }
})
  
  # Create the browser view if data has been loaded
  output$browser <- renderPlot({
    if (is.null(dat())){return()}
    plotBrowser(rv$data, start = rv$start, stop = rv$end)
  })
  
  # If selected chromosome changes, re-sample data to plot the coverage distribution
  # Take 5000 lines randomly to generate the distribution
  d <- eventReactive(rv$data, {
    if (is.null(dat())){return()}
    rv$data[sample(nrow(rv$data), 5000, replace = F),]
  })
  
  # Create density plot showing distibution of sequencing depth
  output$hist <- renderPlot({
    if (is.null(dat())){return()}
    med <- median(d()$depth)
    ggplot(d(), aes(depth)) +  geom_density(fill="blue", alpha = 0.8) + theme_bw() +
      geom_vline(xintercept = med) + ggtitle(paste("Sequencing depth for ", as.character(input$chrom)))
  })
  
  # Create summary statistics based on input data
  output$stats <- renderTable({
    if (is.null(dat())){return()}
    methylomeStats(dat())
  })
  
  # Plot coverage of the lambda chromsome
  output$lambdaCov <- renderPlot({
    if (is.null(dat())){return()}
    ggplot(subset(dat(), dat()$chr == "L"), aes(depth)) + geom_histogram(fill = "darkgreen", color="black") +
      theme_bw() + ggtitle("Lambda chromosome sequencing depth")
  })
  
  # Generate the non-conversion statistics
  output$nonconversion <- renderTable({
    if (is.null(dat())){return()}
    nonConversion(dat())
  })
  
  # Create survival plot for the chosen chromosome
  output$survival <- renderPlotly({
    if (is.null(dat())){return()}
    ggplotly(plotSurvival(rv$data, species = "arabidopsis", chromosome = input$chrom))
  })
  
  # Show the first 50 rows of the input data (shown in data tab)
  output$headData <- renderTable({
    if (is.null(dat())){return()}
    head(rv$data, n = 50)
  })
  
  # Report generation. Uses template at ./report.Rmd
  output$downloadReport <- downloadHandler(
    filename = function() {
      fname <- unlist(strsplit(input$filename, ".", fixed = TRUE))[[1]] # remove extension if present
      paste(fname, sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html'
      ))
    },
    
    content = function(file) {
      src <- normalizePath('report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd')
      
      out <- rmarkdown::render('report.Rmd',
                               params = list(data = dat()),
        switch(
        input$format,
        PDF = rmarkdown::pdf_document(), HTML = rmarkdown::html_document()
      ))
      file.rename(out, file)
    }
  )
}

# Run the app
shinyApp(server = server, ui = ui)