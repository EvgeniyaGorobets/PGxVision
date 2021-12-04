library(shiny)
library(shinydashboard)

gsTypes <- unique(msigdbr::msigdbr_collections()$gs_subcat)


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
  column(
    width = 6,
    fileInput("drugSensFile", "Upload Drug Sensitivity CSV:",
              accept = c("text/csv", ".csv"), buttonLabel="Browse files")
  ),
  column(
    width=6,
    fileInput("forestFile", "Upload Forest Plot CSV???:",
              accept = c("text/csv", ".csv"), buttonLabel="Browse files")
  ),
  p("NOTE: If no drug sensitivity file is uploaded, then a sample dataset
  will be used by default (see LINK for more). (forest plot???)")
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


geneInfoBox <- box(width = 12,
  column(
    width = 9,
    h3("Biomarker Info"),
    div(tags$b("Gene: "), textOutput("geneName", inline = TRUE)),
    div(tags$b("Genome Position: "), textOutput("geneStart", inline = TRUE)),
    div(tags$b("Chromosome: "), textOutput("geneChr", inline = TRUE)),
    div(tags$b("Estimate: "), textOutput("geneEstimate", inline = TRUE)),
    div(tags$b("P-Value: "), textOutput("genePVal", inline = TRUE)),
    div(tags$b("FDR: "), textOutput("geneFdr", inline = TRUE))
  ),
  column(
    width = 3, br(), br(),
    p("Click on any point to see more information about the gene."), br()
  )
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
  column(width = 9, plotOutput("networkPlot", height = 400)),
  # FIXME: plot doesn't always fully show, overflows
  column(
    width = 3, br(),
    sliderInput("simCutoff", "Similarity Cutoff",
                min = 0, max = 1, value = 0.5),
    br(), h4("Gene Set Info"), p("Under construction")
  ),
)


waterfallPlotBox <- box(
  width = 6,
  column(width = 4, uiOutput("wfXSelect")),
  column(width = 4, uiOutput("wfYSelect")),
  column(width = 4, uiOutput("wfColSelect")),
  plotOutput("waterfallPlot", height = 350)
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
          box(plotOutput("manhattanPlot", click = "mouseClick", height = 350)),
          box(plotOutput("volcanoPlot", click = "mouseClick", height = 350))
        ),
        fluidRow(geneInfoBox),
        fluidRow(geneSetAnalysisBox),
      ),

      # Drug Response tab
      tabItem(
        tabName = "drugResponse",
        h2("Drug Response"),
        fluidRow(sensitivityFileUploadBox),
        fluidRow(
          waterfallPlotBox,
          box(plotOutput("forestPlot"))
        )
      )
    )

  )
)

server <- function(input, output) {
  # Define reactive values
  rv <- reactiveValues(biomarkerDf = Biomarkers,
                       chromosomeDf = GRCh38.p13.Assembly,
                       plottedBiomrkrs = NULL,
                       selectedGene = NULL,
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

  # Update experiment based on dropdown selections
  experiment <- reactive({
    setNames(c(input$tissue, input$compound, input$mDataType),
             c("tissue", "compound", "mDataType"))
  })

  # Update plots based on experiment
  output$manhattanPlot <- renderPlot({
    result <- buildManhattanPlot(rv$biomarkerDf, rv$chromosomeDf, experiment(),
                                 pValCutoff = input$pValCutoff)
    rv$plottedBiomrkrs <- result$dt
    result$plot
  })

  output$volcanoPlot <- renderPlot({
    buildVolcanoPlot(
      rv$biomarkerDf, experiment(), pValCutoff = input$pValCutoff)$plot
  })

  # Update selected gene when users click on plots
  observeEvent(input$mouseClick, {
    req(input$mouseClick)
    req(rv$plottedBiomrkrs)
    rv$selectedGene <- nearPoints(rv$plottedBiomrkrs, input$mouseClick,
                                  threshold = 5, maxpoints = 1)
  })

  # Update biomarker info box in response to new selected gene
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

  # Update & rerender network plot only when the runGsAnalysis button is pressed
  observeEvent(input$runGsAnalysis, {
    rv$gsSimilarityDf <- geneSetAnalysis(
      input$gene, input$gsType, input$simAlgo)
  })

  output$networkPlot <- renderPlot({
    req(rv$gsSimilarityDf)
    buildNetworkPlot(rv$gsSimilarityDf, input$simCutoff)
  })

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
