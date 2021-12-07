library(shiny)
library(shinydashboard)
library(plotly)
library(visNetwork)

gsTypes <- unique(msigdbr::msigdbr_collections()$gs_subcat)
blankGene <- data.table::data.table(gene="", abs_gene_seq_start="", chr="",
                                    pvalue="", estimate="", fdr="")


biomarkerFileUploadBox <- box(
  width = 12,
  column(
    width = 6,
    fileInput("biomarkerFile", "Upload Biomarker CSV:",
              accept = c("text/csv", ".csv"), buttonLabel="Browse files")
  ),
  column(
    width=6,
    fileInput("genomeFile", "Upload Genome CSV:",
              accept = c("text/csv", ".csv"), buttonLabel="Browse files")
  ),
  p("NOTE: If no biomarker file is uploaded, then a sample biomarker dataset
  will be used by default (see LINK for more). If no genome file is uploaded,
  then the GRCh38.p13.Assembly genome will be used automatically (see LINK for
  more).")
)


sensitivityFileUploadBox <- box(
  width = 12,
  fileInput("drugSensFile", "Upload Drug Sensitivity CSV:",
            accept = c("text/csv", ".csv"), buttonLabel="Browse files")
)


filterBiomarkersBox <- box(
  width = 12, title = "Filter Biomarkers",
  column(width = 4, uiOutput("tissueSelect")),
  column(width = 4, uiOutput("compoundSelect")),
  column(width = 4, uiOutput("mDataTypeSelect"))
)


plotPropertiesBox <- box(
  width = 12, collapsible = T, collapsed = T, title = "Plot Properties",
  column(
    width = 4,
    sliderInput("pValCutoff", "P-Value Cutoff", min=0, max=1, value=0.05)),
  # TODO: something funny happens when pval >= 0.39 ?? consider
  # reducing range of slider
  column(width = 4, p("Color manipulation: under construction")),
  column(width = 4, p("Title/axis label manipulation: under construction"))
)


geneSetAnalysisBox <- box(
  width=12,
  h3("Gene Set Analysis"),
  column(width = 3, uiOutput("geneSelect")),
  column(
    width = 3,
    selectInput("gsType", "Gene Set Type",
                c("Select a gene set type..." = "", gsTypes), selected = "")
  ),
  column(
    width = 3, selectInput("simAlgo", "Similarity Algorithm", c("overlap"))
  ),
  column(
    width = 3, br(),
    actionButton("runGsAnalysis", "Run Gene Set Analysis!")
  ),
  column(width = 9, visNetworkOutput("networkPlot")),
  column(
    width = 3, br(),
    sliderInput("simCutoff", "Similarity Cutoff",
                min = 0, max = 1, value = 0.5),
    br(), h4("Gene Set Info"), uiOutput("gsInfo")
  ),
)


waterfallPlotBox <- box(
  width = 12,
  column(width = 4, uiOutput("wfXSelect")),
  column(width = 4, uiOutput("wfYSelect")),
  column(width = 4, uiOutput("wfColSelect")),
  plotOutput("waterfallPlot")
)


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
        fluidRow(biomarkerFileUploadBox),
        fluidRow(filterBiomarkersBox),
        fluidRow(plotPropertiesBox),
        fluidRow(
          box(plotlyOutput("manhattanPlot")),
          box(plotlyOutput("volcanoPlot"))
        ),
        fluidRow(uiOutput("geneInfoBox")),
        fluidRow(geneSetAnalysisBox),
      ),

      # Drug Response tab
      tabItem(
        tabName = "drugResponse",
        h2("Drug Response"),
        fluidRow(sensitivityFileUploadBox),
        fluidRow(waterfallPlotBox)
      )
    )

  )
)

server <- function(input, output) {
  # Define reactive values
  rv <- reactiveValues(biomarkerDf = Biomarkers,
                       chromosomeDf = GRCh38.p13.Assembly,
                       plottedBiomrkrs = NULL,
                       selectedGene = blankGene,
                       gsSimilarityDf = NULL,
                       sensitivityDf = BRCA.PDXE.paxlitaxel.response)

  # ------------------ BIOMARKER TAB ------------------ #

  # Update biomarkerDf and chromosomeDf based on file uploads
  observeEvent(input$biomarkerFile, {
    req(input$biomarkerFile)
    rv$biomarkerDf <- read.csv(input$biomarkerFile$datapath)
  })

  observeEvent(input$genomeFile, {
    req(input$genomeFile)
    rv$chromosomeDf <- read.csv(input$genomeFile$datapath)
  })

  # Get all gene, tissues, compounds, and mDataTypes from biomarkerDf
  # and update dropdown elements accordingly
  output$geneSelect <- renderUI({
    geneChoices <- unique(rv$biomarkerDf$gene) #FIXME: unsafe
    selectInput("gene", "Gene", c("Select a gene..." = "", geneChoices),
                selected = "")
  })

  output$tissueSelect <- renderUI({
    tissueChoices <- unique(rv$biomarkerDf$tissue) #FIXME: unsafe
    selectInput("tissue", "Tissue", c("Select a tissue..." = "", tissueChoices),
                selected = "")
  })

  output$compoundSelect <- renderUI({
    compoundChoices <- unique(rv$biomarkerDf$compound) #FIXME: unsafe
    selectInput("compound", "Compound/Drug",
                c("Select a compound..." = "", compoundChoices), selected = "")
  })

  output$mDataTypeSelect <- renderUI({
    mDataTypeChoices <- unique(rv$biomarkerDf$mDataType) #FIXME: unsafe
    selectInput("mDataType", "Molecular Data Type",
                c("Select a molecular data type..." = "", mDataTypeChoices),
                selected = "")
  })

  # Update plots based on tissue, compound, mDataType selections
  output$manhattanPlot <- renderPlotly({
    req(rv$biomarkerDf)

    # Wait until everything renders and inputs are not NULL
    if (typeof(input$tissue) == "character" &&
        typeof(input$compound) == "character" &&
        typeof(input$mDataType) == "character") {
      result <- suppressWarnings(
        buildManhattanPlot(rv$biomarkerDf, rv$chromosomeDf, input$tissue,
                           input$compound, input$mDataType, input$pValCutoff))
      rv$plottedBiomrkrs <- result$dt
      ggplotly(result$plot, source = "manhattan")
    }
  })

  # TODO: improve styling by using %>% layout()
  # https://plotly.com/ggplot2/getting-started/
  # Styling options: https://plotly.com/r/font/
  output$volcanoPlot <- renderPlotly({
    req(rv$biomarkerDf)

    # Wait until everything renders and inputs are not NULL
    if (typeof(input$tissue) == "character" &&
        typeof(input$compound) == "character" &&
        typeof(input$mDataType) == "character") {
      p <- suppressWarnings(
        buildVolcanoPlot(rv$biomarkerDf, input$tissue, input$compound,
                         input$mDataType, pValCutoff = input$pValCutoff)$plot)
      ggplotly(p, source = "volcano")
    }
  })

  # Update selected gene when users click on plots
  observe({
    req(rv$plottedBiomrkrs)

    if (!is.null(event_data("plotly_click", source = "manhattan"))) {
      d <- event_data("plotly_click", source = "manhattan")
      rv$selectedGene <- rv$plottedBiomrkrs[abs_gene_seq_start == d$x,]
    }
  })

  observe({
    req(rv$plottedBiomrkrs)

    if (!is.null(event_data("plotly_click", source = "volcano"))) {
      d <- event_data("plotly_click", source = "volcano")
      rv$selectedGene <- rv$plottedBiomrkrs[estimate == d$x,]
    }
  })

  #TODO: I want the equivalent data points on the other plot to highlight

  # Update biomarker info box in response to new selected gene
  output$geneInfoBox <- renderUI({
    req(rv$selectedGene)

    box(
      width = 12,
      column(
        width = 9,
        h3("Biomarker Info"),
        div(tags$b("Gene: "), rv$selectedGene[1, gene]),
        div(tags$b("Genome Position: "),rv$selectedGene[1, abs_gene_seq_start]),
        div(tags$b("Chromosome: "), rv$selectedGene[1, chr]),
        div(tags$b("Estimate: "), rv$selectedGene[1, estimate]),
        div(tags$b("P-Value: "), rv$selectedGene[1, pvalue]),
        div(tags$b("FDR: "), rv$selectedGene[1, fdr])
      ),
      column(
        width = 3, br(), br(),
        p("Click on any point to see more information about the gene."), br()
      )
    )
  })

  # Update & rerender network plot only when the runGsAnalysis button is pressed
  observeEvent(input$runGsAnalysis, {
    rv$gsSimilarityDf <- geneSetAnalysis(
      input$gene, input$gsType, input$simAlgo)
  })

  output$networkPlot <- renderVisNetwork({
    req(rv$gsSimilarityDf)
    p <- buildNetworkPlot(rv$gsSimilarityDf, input$simCutoff) %>%
      # Add JS hook to react to node selection
      # Js code taken from xclotet:
      # xclotet. (2016). Get selected Node data from visNetwork graph without
      # actionButton. StackOverflow.
      # https://stackoverflow.com/questions/41018899/get-selected-node-data-from-visnetwork-graph-without-actionbutton/41020222
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('currentNodeId', nodes.nodes);}")
    p
  })

  # React to node selection
  output$gsInfo <- renderUI({
    print(input$currentNodeId)
    id <- ifelse(length(input$currentNodeId)>0, input$currentNodeId, "")
    div(tags$b("ID: "), id)
  })

  #observeEvent(input$currentNodeId, {
  #  print(input$currentNodeId)
  #})

  # ------------------ END BIOMARKER TAB ------------------ #

  # ------------------ DRUG RESPONSE TAB ------------------ #

  # Update waterfallDf based on file upload
  observeEvent(input$drugSensFile, {
    req(input$drugSensFile)
    rv$sensitivityDf <- read.csv(input$drugSensFile$datapath)
  })

  # Create waterfall plot dropdowns based on file upload
  output$wfXSelect <- renderUI({
    selectInput("wfX", "x-Axis", colnames(rv$sensitivityDf), selectize = F)
  })

  output$wfYSelect <- renderUI({
    selectInput("wfY", "y-Axis", colnames(rv$sensitivityDf), selectize = F)
  })

  output$wfColSelect <- renderUI({
    selectInput("wfCol", "Color", colnames(rv$sensitivityDf), selectize = F)
  })


  # Update plots based on file upload
  output$waterfallPlot <- renderPlot({
    buildWaterfallPlot(
      rv$sensitivityDf, xAxisCol=input$wfX, drugSensitivityCol=input$wfY,
      colorCol=input$wfCol, xLabel="Tumour",
      yLabel="Angle Between Treatment and Control",
      title="Paclitaxel Response in BRCA Tumours")
  })

  # ------------------ END DRUG RESPONSE TAB ------------------ #
}

shinyApp(ui, server)

# [END]
