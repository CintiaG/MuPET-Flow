.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    packageStartupMessage(
      "Warning: The 'flowCore' package is required but not installed.\n",
      "Please install it manually using:\n",
      "BiocManager::install('flowCore')\n"
    )
  } else {
    packageStartupMessage("Thank you for using MuPETFlow!")
  }
}
