### INPUT ###

# Helper functions
clinicalBiomarkerFileUploadBox <- box(
  width=12,
  column(
    width=6,
    fileInput("patientFile", "Upload Patient Molecular Profile CSV:",
              accept=c("text/csv", ".csv"), buttonLabel="Browse files"
    )
  ),
  column(
    width=6,
    fileInput("referencePopulationFile",
              "Upload Reference Population CSV:",
              accept=c("text/csv", ".csv"), buttonLabel="Browse files"
    )
  )
)

filterClinicalBiomarkersBox <- box(
  width=12,
  column(width=4,
         uiOutput("featureSelect")
  ),
)

clinicalBiomarkerDensityPlot <- box(
  width=12,
  column(width=2),
  column(width=8, align="center", plotOutput("clinicalDensity")),
  column(width=2)
)

pharmacodbBiomarkersTable <- box(
  width=12,
  column(width=8, align="center", dataTableOutput("pdbBiomarkerDfFiltered"))
)

# Return all tab input rows 
clinicalTabInputUI <- tabItem(
  tabName = "clinical",
  h2("Single Sample Predictor"),
  fluidRow(clinicalBiomarkerFileUploadBox),
  fluidRow(filterClinicalBiomarkersBox),
  fluidRow(clinicalBiomarkerDensityPlot),
  fluidRow(pharmacodbBiomarkersTable)
)

### REACTIVE VALUES AND OBSERVERS ###

# object containing all relevant reactive values for the clinical tab
clinicalTabCreateRV <- function() {
  
  return( reactiveValues(patientDf = NULL,
                         referenceDf = NULL,
                         pdbBiomarkersDf=tryCatch(
                           PGxVision::fetchPharmacoDbBiomarkers(),
                           error=function(e) fread(file.path(
                             system.file("extdata", package="PGxVision"),
                             "sample_data",
                             "pharmacodb_biomarker_df.csv"
                           )))))
}

# Return all reactive variable observers
clinicalTabObservers <- function(input, rv) {
  observeEvent(input$patientFile, {
    req(input$patientFile)
    df_ <- data.table::fread(input$patientFile$datapath)
    rv$patientDf <- data.frame(df_[, -1], row.names=df_[[1]])
  })
  
  observeEvent(input$referencePopulationFile, {
    req(input$referencePopulationFile)
    df_ <- data.table::fread(input$referencePopulationFile$datapath)
    rv$referenceDf <- data.frame(df_[, -1], row.names=df_[[1]])
  })
}

### OUTPUT ###

#Return all output objects
clinicalTabOutputUI <- function (input, rv, output) {
  output$featureSelect <- renderUI({
    featureChoices <- unique(row.names(rv$patientDf))
    selectInput("feature", "Feature", c("Select a feature..." = "",
                                        featureChoices), selected="ERBB2")
  })
  
  output$clinicalDensity <- renderPlot(
    PGxVision::densityPlotBiomarkerPercentile(
      input$feature,
      rv$patientDf[input$feature, ],
      rv$referenceDf)
  )
  
  output$pdbBiomarkerDfFiltered <- renderDataTable({
    df_ <- if (input$feature != "") {
      rv$pdbBiomarkersDf[gene_symbol == input$feature]
    } else {
      rv$pdbBiomarkersDf
    }
    df_[
      order(pvalue, -abs(estimate)),
      .(compound_name, correlation=estimate, pvalue, inchikey, pubchem,
        chembl_id, tissue)
    ]
  })
}