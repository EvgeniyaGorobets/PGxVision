#library(shiny)
#library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "PGxVision"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Biomarkers", tabName = "biomarkers", icon = icon("dna")),
      menuItem("Treatment Response", tabName = "drugResponse",
              icon = icon("capsules"))
    )
  ),

  dashboardBody(
    tabItems(
      # Biomarkers tab
      tabItem(
        tabName = "biomarkers",
        h2("Biomarker Analysis"),
        fluidRow(box(width=12,
          column(width=6, fileInput(
            "biomarkerFile", "Upload Biomarker CSV:",
            accept = c("text/csv", ".csv"), buttonLabel="Browse files")),
          column(width=6, fileInput(
            "genomeFile", "Upload Genome JSON:",
            accept = c("text/csv", ".csv"), buttonLabel="Browse files")),
          p("NOTE: If no biomarker file is uploaded, then a sample biomarker
              dataset will be used by default (see LINK for more).")
        )),
        fluidRow(box(width=12,
          column(width=4, uiOutput("tissueSelect")),
          column(width=4, uiOutput("compoundSelect")),
          column(width=4, uiOutput("mDataTypeSelect"))
        )),
        fluidRow(
          box(plotOutput("manhattanPlot", click = "mouseClick", height = 350)),
          box(plotOutput("volcanoPlot", click = "mouseClick", height = 350))
        ),
        fluidRow(box(width=12,
          column(
            width=9,
            div(tags$b("Gene: "), textOutput("geneName", inline = TRUE)),
            div(tags$b("Genome Position: "),
                textOutput("geneStart", inline = TRUE)),
            div(tags$b("Chromosome: "), textOutput("geneChr", inline = TRUE)),
            div(tags$b("Estimate: "),
                textOutput("geneEstimate", inline = TRUE)),
            div(tags$b("P-Value: "), textOutput("genePVal", inline = TRUE)),
            div(tags$b("FDR: "), textOutput("geneFdr", inline = TRUE))),
          column(width=3, actionButton(
            "runGsAnalysis", "Run Gene Set Analysis!"))
                  #TODO: change styling of button
        ))
      ),

      # Drug Response tab
      tabItem(
        tabName = "drugResponse",
        h2("Drug Response"),
        fluidRow(
          box(plotOutput("waterfallPlot", height = 350)),
          box(plotOutput("forestPlot", height = 350))
        )
      )
    )

  )
)

server <- function(input, output) {
  # Define reactive values biomarkerDf and chromosomeDf
  rv <- reactiveValues(biomarkerDf = Biomarkers,
                       chromosomeDf = GRCh38.p13.Assembly,
                       plottedBiomrkrs = NULL,
                       selectedGene = NULL)

  # Update biomarkerDf and chromosomeDf based on file uploads
  observeEvent(input$biomarkerFile, {
    req(input$biomarkerFile)
    rv$biomarkerDf <- read.csv(input$biomarkerFile$datapath)
  })

  # Get all tissues, compounds, and mDataTypes from biomarkerDf
  # and update dropdown elements accordingly
  output$tissueSelect <- renderUI({
    tissueChoices <- unique(rv$biomarkerDf$tissue) #FIXME: unsafe
    selectInput("tissue", "Tissue",
                c("Select a tissue..." = "", tissueChoices),
                selected = "")
  })

  output$compoundSelect <- renderUI({
    compoundChoices <- unique(rv$biomarkerDf$compound) #FIXME: unsafe
    selectInput("compound", "Compound/Drug",
                c("Select a compound..." = "", compoundChoices),
                selected = "")
  })

  output$mDataTypeSelect <- renderUI({
    mDataTypeChoices <- unique(rv$biomarkerDf$mDataType) #FIXME: unsafe
    selectInput("mDataType", "Molecular Data Type",
                c("Select a molecular data type..." = "", mDataTypeChoices),
                selected = "")
  })

  # Update experiment based on dropdown selections
  experiment <- reactive({
    setNames(c(input$tissue, input$compound, input$mDataType),
             c("tissue", "compound", "mDataType"))
  })

  # Update plots based on experiment
  output$manhattanPlot <- renderPlot({
    result <- buildManhattanPlot(rv$biomarkerDf, rv$chromosomeDf, experiment())
    rv$plottedBiomrkrs <- result$dt
    result$plot
  })

  output$volcanoPlot <- renderPlot({
    buildVolcanoPlot(rv$biomarkerDf, experiment())$plot
  })

  # Update selected gene when users click on plots
  observeEvent(input$mouseClick, {
    req(input$mouseClick)
    req(rv$plottedBiomrkrs)
    rv$selectedGene <- nearPoints(rv$plottedBiomrkrs, input$mouseClick,
                                  threshold = 5, maxpoints = 1)
  })

  output$geneName <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, gene]
  })
  output$geneStart <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, abs_gene_seq_start]
  })
  output$geneChr <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, chr]
  })
  output$geneEstimate <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, estimate]
  })
  output$genePVal <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, pvalue]
  })
  output$geneFdr <- renderText({
    req(rv$selectedGene)
    rv$selectedGene[1, fdr]
  })

  # data <- eventReactive(input$runGsAnalysis, {code to run it})
}

shinyApp(ui, server)
