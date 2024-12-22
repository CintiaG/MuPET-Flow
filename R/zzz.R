.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    packageStartupMessage(
      "The 'flowCore' package is required for MuPETFlow to work. ",
      "Please install it manually using:\n",
      "BiocManager::install('flowCore')"
    )
  }
}
