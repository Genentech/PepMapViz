# R/runShinyApp.R

#' Launch PepMapViz Shiny Application
#'
#' This function launches a Shiny application that provides an interactive
#' interface for the PepMapViz package functionality.
#'
#' @param ... Additional arguments to pass to shiny::runApp()
#' @return The Shiny application object
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{
#' run_pepmap_app()
#' }
run_pepmap_app <- function(...) {
  app_dir <- system.file("shiny_apps", "PepMapVizApp", package = "PepMapViz")  # Correct path

  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try reinstalling 'PepMapViz'.", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
