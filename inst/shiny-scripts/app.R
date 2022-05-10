library(shiny)
library(shinybusy)
library(shinydashboard)
library(plotly)
library(magrittr)
library(jsonlite)
library(data.table)

# Import files for each tab 
source ('clinicalTab.R')
source ('logicalTab.R')
source ('biomarkerTab.R')
source ('treatmentTab.R')

# Change options for this session
opts <- options()
options(shiny.maxRequestSize=250*1024^2) # max upload of 250 MB
on.exit(options(opts))

ui <- dashboardPage(
  dashboardHeader(title = "PGxVision"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Clinical Biomarkers", tabName="clinical",
        icon=icon("stethoscope")),
      menuItem("Logical Model", tabName="logicalModel",
        icon=icon("stethoscope")),
      menuItem("Biomarker Analysis", tabName = "biomarkers", icon = icon("dna")),
      menuItem("Treatment Response", tabName = "drugResponse",
              icon = icon("capsules"))
    )
  ),
  dashboardBody(
    # Custom CSS to put shinybusy spinner in center of network plot
    tags$head(
      tags$style(HTML("
      .shinybusy {
        position: relative !important;
      }
      .shinybusy-busy {
        margin: auto;
      }"))
    ),
    tabItems(
      clinicalTabInputUI,
      logicalTabInputUI,
      biomarkerTabInputUI,
      treatmentTabInputUI
    )
  )
)

server <- function(input, output) {

  # Define reactive values
  clinicalTabRV <- clinicalTabCreateRV()
  logicalTabRV <- logicalTabCreateRV()
  biomarkerTabRV <- biomarkerTabCreateRV()
  treatmentTabRV <- treatmentTabCreateRV()
  
  # Observers
  clinicalTabObservers(input, clinicalTabRV)
  logicalTabObservers(input, logicalTabRV)
  biomarkerTabObservers(input, biomarkerTabRV)
  treatmentTabObservers(input, treatmentTabRV)
  
  # Output UI
  clinicalTabOutputUI(input, clinicalTabRV, output)
  logicalTabOutputUI(input, logicalTabRV, output)
  biomarkerTabOutputUI(input, biomarkerTabRV, output)
  treatmentTabOutputUI(input, treatmentTabRV, output)
}

shinyApp(ui, server)

# [END]
