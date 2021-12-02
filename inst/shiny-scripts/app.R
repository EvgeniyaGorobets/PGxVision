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
          column(width=4, selectInput("tissue", "Tissue",
                      c("Select a tissue..." = "", "Lung", "Breast"),
                      selected = "")),
          column(width=4, selectInput("compound", "Compound/Drug",
                      c("Select a compound..." = "", "Trametinib", "Daporinad", "Dasatinib"),
                      selected = "")),
          column(width=4, selectInput("mDataType", "Molecular Data Type",
                      c("Select a molecular data type..." = "", "RNA"),
                      selected = ""))
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

server <- function(input, output) { }

shinyApp(ui, server)
