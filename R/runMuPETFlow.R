#' Launch the MuPETFlow app
#'
#' This function launches the Shiny app included in this package.
#' @export
#' @import shiny
#' @import shinythemes
#' @import zoo
#' @import ggplot2
#' @import DT
#' @import ggrepel
#' @import gridExtra
#' @import markdown
#' @import BiocManager
#' @import tidyr
#' @import dplyr
runMuPETFlow <- function(){
  appDir <- system.file("shiny", "app", package = "MuPETFlow")
  if(appDir == ""){
    stop("Could not find the directory containing the Shiny app. Try re-installing `MuPETFlow`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
