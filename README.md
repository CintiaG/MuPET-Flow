# MuPET-Flow: Multiple Ploidy Estimation Tool from Flow Cytometry data

Created: January 30, 2023

Updated: January 5, 2024

Authors: Gómez-Muñoz, C.*

*cintia.gomez_munoz@sorbonne-universite.fr

---

## Introduction

MuPET-Flow is a Shiny app intended to analyze flow cytometry data to produce fluorescence histograms, estimate peaks' intensity and calculate ploidy using known standards.

We will update the app description shortly. Feel free to report any bugs you encounter.

## Installation

MuPET-Flow can be installed by downloading the files on this repository manually or via **GitHub CLI**:

```bash
gh repo clone CintiaG/MuPET-Flow
```

`app.R` file needs to be run in **RStudio** by clicking the “Run App” button.


### Pre-requisites

These are the libraries required to run MuPET-Flow. Here we provide links to their source or installation instructions.

* [shiny](https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/index.html)
* [shinythemes](https://rstudio.github.io/shinythemes/)
* [flowCore](https://bioconductor.org/packages/release/bioc/html/flowCore.html)
* [zoo](https://cran.r-project.org/web/packages/zoo/index.html)
* [ggplot2](https://ggplot2.tidyverse.org/)
* [DT](https://rstudio.github.io/DT/)
* [tidyverse](https://www.tidyverse.org/packages/)
* [ggrepel](https://cran.r-project.org/web/packages/ggrepel/readme/README.html)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
* [markdown](https://cran.r-project.org/web/packages/markdown/index.html)

<!--
Pending
A minimum of two different standards is required, but more are recommended.
 -->
