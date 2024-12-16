#' Run the MuPET-Flow app
#'
#' This function launches the Shiny app included in this package.
#' The first time the tool launches, it will try to install additional **Bioconductor** dependencies that are not managed via `install.packages()`. This process might take a few minutes.
#' Once the application is launched, you can either:
#' 1. Load your experimental data.
#' 2. Run an in-app example by clicking the 'Example' button.
#' <br>
#' <br>For the first case, selecting the channel where the data was acquired is mandatory. For the second case, the example channel `FL4-A` is automatically detected for demonstration purposes.
#'
#' After launching the app, you can follow the app flow, which is divided into three tabs: _**Peaks**_, _**Regression**_ and _**Summary**_.  Below is a general description of the options available in each tab:

#' ## Peaks
#' * **Select a sample (optional):** Allows visual exploration of individual samples if desired.
#' * **Adjust smoothing (optional):** Adjusts the histogram curve for noisy samples.
#' * **Adjust window width (optional):** Defines the interval where the app will look for peaks.
#' * **Select minimum cell count to call a peak (optional):** Useful for samples with a low number of events.
#' * **Select maximum number of peaks to plot (optional):** Useful for samples with heterogeneous populations where more peaks are present.

#' ## Regression
#' * **Select type of analysis:** Choose between "Ploidy" or "Genome size" analysis.
#' * **Select number of standards:** A minimum of two standards is required.
#' * **Select standard samples and values:** This is the ploidy or genome size of your standards.

#' ## Summary
#' * **Results preview:** Creates a compiled figure with histograms for all samples.
#' * **Save plot:** Saves the histograms in either PNG or TIFF format with customizable size and quality. Optionally, you can control the grid layout.
#' * **Save table:** Exports the parameters used and the estimated ploidy or genome size as a CSV file.
#'
#' <br>
#' If any errors are detected, you can go back to the Peaks tab to review the parameters. Note that the regression must also be re-done after parameter adjustments.
#'
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
#' @importFrom dplyr mutate
runMuPETFlow <- function(){
  appDir <- system.file("shiny", "app", package = "MuPETFlow")
  if(appDir == ""){
    stop("Could not find the directory containing the Shiny app. Try re-installing `MuPETFlow`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
