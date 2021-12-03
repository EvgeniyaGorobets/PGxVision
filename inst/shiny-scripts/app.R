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
          box(plotOutput("manhattanPlot", hover = "mouseHover", height = 350)),
          box(plotOutput("volcanoPlot", hover = "mouseHover", height = 350))
        ),
        fluidRow(
          box(tableOutput("data"))
        )
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
  # Update biomarkerDf based on file uploads
  biomarkerDf <- reactive({
    biomarkerFile <- input$biomarkerFile
    if (!is.null(biomarkerFile)) {
      read.csv(input$biomarkerFile$datapath)
    } else {
      Biomarkers
    }
  })

  # Get all tissues, compounds, and mDataTypes from biomarkerDf
  # and update dropdown elements accordingly
  output$tissueSelect <- renderUI({
    tissueChoices <- unique(biomarkerDf()$tissue) #FIXME: unsafe
    selectInput("tissue", "Tissue",
                c("Select a tissue..." = "", tissueChoices),
                selected = "")
  })

  output$compoundSelect <- renderUI({
    compoundChoices <- unique(biomarkerDf()$compound) #FIXME: unsafe
    selectInput("compound", "Compound/Drug",
                c("Select a compound..." = "", compoundChoices),
                selected = "")
  })

  output$mDataTypeSelect <- renderUI({
    mDataTypeChoices <- unique(biomarkerDf()$mDataType) #FIXME: unsafe
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
    buildManhattanPlot(biomarkerDf(), GRCh38.p13.Assembly, experiment())
  })

  output$volcanoPlot <- renderPlot({
    buildVolcanoPlot(biomarkerDf(), experiment())
  })

  output$data <- renderTable({
    req(input$mouseHover)
    nearPoints(biomarkerDf(), input$mouseHover)
  })
}

shinyApp(ui, server)
