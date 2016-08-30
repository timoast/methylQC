library(shiny)
library(ggplot2)
library(dplyr)
library(methylQC)
library(shinythemes)
library(plotly)
library(moments)

dat <- loadData(path="../data/example_data.CGmap.gz", datasource = "BSseeker2")
chroms <- unique(dat$chr)
m <- max(dat$position)

ui <- fluidPage(theme = "theme.css",
                titlePanel("methylQC"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('data', 'Choose a file'),
                    radioButtons('type', 'Source',
                                 c(BSseeker2='bsseeker2',
                                   Bismark='bismark',
                                   Bisulfiter='bisulfiter')),
                    checkboxInput(inputId = "demo", label = "Use demo data", value = FALSE),
                    tags$hr(),
                    selectInput(inputId = "chrom", label = "Chr", choices = chroms, selected = "chr1"),
                    numericInput(inputId = "position",
                                 label = "Position", min = 0, max = m, value = 5000, step = 2000),
                    sliderInput(inputId = "zoom", label = "Zoom", min = 1, max = 100, value = 15),
                    tags$hr(),
                    radioButtons('format', 'Document format', c('PDF', 'HTML'),
                                 inline = TRUE),
                    downloadButton("downloadReport", "Generate report")
                  ),
                  mainPanel(
                    tabsetPanel(type = "tabs",
                        tabPanel("About",
                                 includeMarkdown("./text/about.md")),
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
                        tabPanel("Data", tableOutput("dat"))
                    )
                )
          )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2) # max 1 GB files
  
  rv <- reactiveValues(data = dat, start = 2000, end = 5000)  # need to change to input$data, but have set to NULL initially
  observeEvent(input$chrom, {rv$data <- subset(dat, dat$chr == input$chrom)})
  
  observeEvent(input$zoom, {rv$start <- input$position - (300 * input$zoom)
  rv$end <- input$position + (300 * input$zoom)})
  
  observeEvent(input$position, {rv$start <- input$position - (300 * input$zoom)
  rv$end <- input$position + (300 * input$zoom)})
  
  output$browser <- renderPlot({
    plotBrowser(rv$data, start = rv$start, stop = rv$end)
  })
  
  d <- eventReactive(rv$data, {
    rv$data[sample(nrow(rv$data), 5000, replace = F),]
  })
  
  output$hist <- renderPlot({
    med <- median(d()$depth)
    ggplot(d(), aes(depth)) +  geom_density(fill="blue", alpha = 0.8) + theme_bw() +
      geom_vline(xintercept = med) + ggtitle(paste("Sequencing depth for ", as.character(input$chrom)))
  })
  
  output$stats <- renderTable({
    methylomeStats(dat)
  })
  
  l <- subset(dat, dat$chr == "L")
  
  output$lambdaCov <- renderPlot({
    ggplot(l, aes(depth)) + geom_histogram(fill = "darkgreen", color="black") +
      theme_bw() + ggtitle("Lambda chromosome sequencing depth")
  })
  
  output$nonconversion <- renderTable({
    l %>%
      group_by(context) %>%
      mutate(`Non-conversion rate` = sum(mC) / sum(depth) * 100) %>%
      select(context, `Non-conversion rate`) %>%
      unique() %>%
      filter(context %in% c("CG", "CHG", "CHH")) %>%
      arrange(context)
  })
  
  output$survival <- renderPlotly({
    ggplotly(plotSurvival(rv$data, species = "arabidopsis", chromosome = input$chrom))
  })
  
  output$dat <- renderTable({
    head(rv$data, n = 50)
  })
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('methylQC-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
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
                               params = list(l = l,
                                             data = dat),
        switch(
        input$format,
        PDF = rmarkdown::pdf_document(), HTML = rmarkdown::html_document()
      ))
      file.rename(out, file)
    }
  )
}

shinyApp(server = server, ui = ui)