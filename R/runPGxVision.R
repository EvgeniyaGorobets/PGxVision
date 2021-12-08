#' Launch the PGxVision Shiny app
#'
#' Launches a shinydashboard that lets users interact with PGxVision through
#' UI elements. The app allows users to upload biomarker, drug sensitivity, and
#' genome files, and will plot the data in the files. Users can dynamically
#' rerender plots by choosing different experiments (tissues, compounds,
#' mDataTypes) and can select points on plots to see more information about the
#' gene. Users can also run gene set analysis on any gene in the dataset, and
#' the app will provide an interactive network plot. Users can click on nodes
#' to see more information about the gene set.
#' The code has been placed in \code{./inst/shiny-scripts/app.R}.
#'
#' @return No return value; open up a shiny page
#'
#' @examples
#' \dontrun{
#' runPGxVision()
#' }
#'
#' @importFrom shiny runApp
#' @export
runPGxVision <- function() {
  appDir <- system.file("shiny-scripts", package = "PGxVision")
  shiny::runApp(appDir, display.mode = "normal")
  return(invisible(NULL))
}

# [END]
