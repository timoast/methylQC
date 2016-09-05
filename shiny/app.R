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
  # read cytosine counts for Arabidopsis. Will need to change later to user input
  cytosines <- read.table("../data/arabidopsis.tsv", header = T)
  
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
                           fluidRow(h3("Browser"),
                                    plotOutput("browser")),
                           fluidRow(h3("Cytosine coverage"),
                                    plotlyOutput("survival")),
                           fluidRow(h3("Sequencing depth"),
                                    column(10,
                                    plotOutput("hist")),
                                    column(2,
                                    numericInput(inputId = "covbreaks",
                                                 label = "Breaks", value = 50, min = 1, max = 1000),
                                    numericInput(inputId = "covlim", label = "xlim", value = NA, min = 1))),
                           fluidRow(h3("Strand bias"),
                                    column(10,
                                    plotOutput("biasplot")),
                                    column(2,
                                    numericInput(inputId = "biasbreaks",
                                                 label = "Breaks", value = 50, min = 1, max = 1000),
                                    numericInput(inputId = "biaslim", label = "xlim", value = NA, min = 1)))
                  ),
                  tabPanel("Statistics",
                           h3("Sequencing depth"),
                           includeMarkdown("./text/statistics_coverage.md"),
                           tableOutput("stats"),
                           h3("Sequencing coverage"),
                           tableOutput("coveragestats"),
                           h3("Non-conversion rates"),
                           tableOutput("nonconversion"),
                           h3("Strand bias"),
                           tableOutput("biastable")),
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
  
  biasChunk <- eventReactive(rv$data, {
    if (is.null(dat())){return()}
    if(nrow(rv$data) < 100000) {
      rv$data[100:nrow(rv$data),]
    } else {
      rv$data[1000:100000,]
    }
  })
  
  bias <- eventReactive(biasChunk(), {
    if (is.null(dat())){return()}
    strandBias(biasChunk())
  })
  
  survival <- eventReactive(rv$data, {
    if (is.null(dat())){return()}
    coverageSurvival(rv$data, cytosines, input$chrom)
  })
  
  # Create summary statistics based on input data
  output$stats <- renderTable({
    if (is.null(dat())){return()}
    methylomeStats(dat())
  })
  
  # Generate the non-conversion statistics
  output$nonconversion <- renderTable({
    if (is.null(dat())){return()}
    nonConversion(dat())
  })
  
  # Create survival plot for the chosen chromosome
  output$survival <- renderPlotly({
    if (is.null(dat())){return()}
    ggplotly(plotCoverage(survival()))
  })
  
  output$coveragestats <- renderTable({
    if (is.null(dat())){return()}
    cc <- dplyr::filter(survival(), depth %in% c(1, 5, 10))
    tidyr::spread(cc, depth, Cytosines)
  })
  
  # Plot histogram of sequencing depth
  output$hist <- renderPlot({
    if (is.null(dat())){return()}
    if(is.na(input$covlim)){
      h <- hist(d()$depth, plot = F)
      xl <- range(h$breaks)
    } else {
      xl <- c(0, input$covlim)
    }
    hist(d()$depth, col = "orangered3",
         main = paste("Sequencing depth for ", as.character(input$chrom)),
         breaks = input$covbreaks, xlim = xl, xlab = "Depth")
  })
  
  output$biasplot <- renderPlot({
    if (is.null(dat())){return()}
    if(is.na(input$biaslim)) {
      h <- hist(bias()$strandBias, plot = FALSE)
      xl <- range(h$breaks)
    } else {
      xl <- c(-input$biaslim, input$biaslim)
    }
    hist(bias()$strandBias, col = "lightblue",
         breaks = input$biasbreaks, xlab = "Strand bias",
         xlim = xl, main = paste("Strand bias for ", as.character(input$chrom)))
  })
  
  output$biastable <- renderTable({
    if (is.null(dat())){return()}
    b <- broom::tidy(summary(bias()$strandBias))
    cbind(b, data.frame(standardDeviation = sd(bias()$strandBias, na.rm = T)))
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
                               params = list(data = dat(), cytosines = cytosines),
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