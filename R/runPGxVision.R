#' Launch the shinydashboard for PGxVision
#'
#' TODO: description
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
