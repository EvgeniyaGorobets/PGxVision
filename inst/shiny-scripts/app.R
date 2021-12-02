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
          column(width=4, uiOutput("tissueSelect")),
          column(width=4, uiOutput("compoundSelect")),
          column(width=4, uiOutput("mDataTypeSelect"))
        )),
        fluidRow(
          box(plotOutput("manhattanPlot", height = 350)),
          box(plotOutput("volcanoPlot", height = 350))
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
  biomarkerDf <- Biomarkers

  # Get all tissues, compounds, and mDataTypes from biomarkerDf
  output$tissueSelect <- renderUI({
    tissueChoices <- unique(biomarkerDf$tissue) #FIXME: unsafe
    selectInput("tissue", "Tissue",
                c("Select a tissue..." = "", tissueChoices),
                selected = "")
  })

  output$compoundSelect <- renderUI({
    compoundChoices <- unique(biomarkerDf$compound) #FIXME: unsafe
    selectInput("compound", "Compound/Drug",
                c("Select a compound..." = "", compoundChoices),
                selected = "")
  })

  output$mDataTypeSelect <- renderUI({
    mDataTypeChoices <- unique(biomarkerDf$mDataType) #FIXME: unsafe
    selectInput("mDataType", "Molecular Data Type",
                c("Select a molecular data type..." = "", mDataTypeChoices),
                selected = "")
  })


  experiment <- c()

  # renderPlot, renderImage, renderDataTable, renderTable, renderText, renderUI, etc.
  getExperiment <- reactive({
    experiment$tissue <- input$tissue
    experiment$compound <- input$compound
    experiment$mDataType <- input$mDataType
    return(experiment)
  })

  output$manhattanPlot <- renderPlot({
    buildManhattanPlot(biomarkerDf, GRCh38.p13.Assembly, getExperiment())
  })

  output$volcanoPlot <- renderPlot({
    buildVolcanoPlot(biomarkerDf, getExperiment())
  })
}

shinyApp(ui, server)
