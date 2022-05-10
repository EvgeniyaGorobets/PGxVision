### INITIALIZE VARIABLES ###

sampleDataDir <- "extdata/sample_data/"
pdxFile <- system.file(paste0(sampleDataDir, "sampleTreatmentResponse.csv"),
                       package="PGxVision")
brcaPdxePaxlitaxelResponse <- read.csv(pdxFile)

### INPUT ###

# Helper functions
sensitivityFileUploadBox <- box(
  width = 12,
  fileInput("drugSensFile", "Upload Drug Sensitivity CSV:",
            accept = c("text/csv", ".csv"), buttonLabel="Browse files")
)

waterfallPlotBox <- box(
  width = 12,
  column(width = 4, uiOutput("wfXSelect")),
  column(width = 4, uiOutput("wfYSelect")),
  column(width = 4, uiOutput("wfColSelect")),
  column(width = 12, plotOutput("waterfallPlot"))
)

# Return all tab input rows 
treatmentTabInputUI = tabItem(
  tabName = "drugResponse",
  h2("Drug Response"),
  fluidRow(sensitivityFileUploadBox),
  fluidRow(waterfallPlotBox),
  fluidRow(uiOutput("wfLabels"))
)

### REACTIVE VALUES AND OBSERVERS ###

treatmentTabCreateRV <- function() {
  return( reactiveValues(sensitivityDf = brcaPdxePaxlitaxelResponse))
}

# Return all reactive variable observers
treatmentTabObservers <- function(input, rv) {
  # Update waterfallDf based on file upload
  observeEvent(input$drugSensFile, {
    req(input$drugSensFile)
    rv$sensitivityDf <- read.csv(input$drugSensFile$datapath)
  })
}

### OUTPUT ###

treatmentTabOutputUI <- function(input, rv, output) {
  # Create waterfall plot dropdowns based on file upload
  output$wfXSelect <- renderUI({
    columns <- colnames(rv$sensitivityDf)
    columnClasses <- sapply(rv$sensitivityDf, class)
    discreteColumns <- columns[columnClasses == "character"]
    selectInput("wfX", "x-Axis", discreteColumns, selectize = F)
  })
  
  output$wfYSelect <- renderUI({
    columns <- colnames(rv$sensitivityDf)
    columnClasses <- sapply(rv$sensitivityDf, class)
    numericColumns <- columns[columnClasses == "numeric"]
    selectInput("wfY", "y-Axis", numericColumns, selectize = F)
  })
  
  output$wfColSelect <- renderUI({
    columns <- colnames(rv$sensitivityDf)
    columnClasses <- sapply(rv$sensitivityDf, class)
    numericColumns <- columns[columnClasses == "numeric"]
    selectInput("wfCol", "Color", numericColumns, selectize = F)
  })
  
  # Create text input elements for custom labeling
  output$wfLabels <- renderUI({
    box(width = 12, title = "Plot Labels",
        column(
          width = 4,
          textInput("wfTitle", "Title", value = "")
        ),
        column(
          width = 4,
          textInput("wfXLabel", "x-Axis Label", value = input$wfX)
        ),
        column(
          width = 4,
          textInput("wfYLabel", "y-Axis Label", value = input$wfY)
        )
    )
  })
  
  # Update plots based on selected columns
  output$waterfallPlot <- renderPlot({
    if (typeof(input$wfX) == "character" &&
        typeof(input$wfY) == "character" &&
        typeof(input$wfCol) == "character") {
      PGxVision::buildWaterfallPlot(
        rv$sensitivityDf, xAxisCol=input$wfX, drugSensitivityCol=input$wfY,
        colorCol=input$wfCol, xLabel=input$wfXLabel, yLabel=input$wfYLabel,
        title=input$wfTitle)
    }
  })
}


