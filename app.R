####################################
# MuPET-Flow                       #
# An app to calculate peaks from   #
# flow cytometer histograms        #
####################################

# Load R packages
library(shiny)
library(shinythemes)
library(flowCore)
library(zoo)
library(ggplot2)
library(DT)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(markdown)

# Define UI
ui <- fluidPage(theme = shinytheme("united"),
                # App name
                navbarPage("MuPET-Flow",
                           # Panel 1
                           tabPanel("Peaks",
                                    # Panel 1 inputs
                                    sidebarPanel(
                                      h3("Inputs"),
                                      # Input files
                                      uiOutput("UiFileUp"),
                                      # Select channel
                                      uiOutput("UiChannelSel"),
                                      # Select sample
                                      uiOutput("UiSampleSel"),
                                      # Adjust smoothing
                                      uiOutput("UiSmoothNum"),
                                      # Adjust window
                                      uiOutput("UiWindowNum"),
                                      # Select minimum number of events
                                      uiOutput("UiEveMin"),
                                      # Select maximum number of peaks
                                      uiOutput("UiPeaksNum"),
                                      # Select calculated peaks to plot
                                      uiOutput("UiPeaksBox"),
                                    ),
                                    # Panel 1 outputs
                                    mainPanel(
                                      # Histogram plot
                                      h3("Histogram"),
                                      verbatimTextOutput("Warn1"),
                                      plotOutput("HisPlot"),
                                      #  Table results 1
                                      h3("Estimated peaks and used parameters"),
                                      DT::dataTableOutput("ResDf1"),
                                    )
                           ),
                           # Panel 2
                           tabPanel("Regression",
                                    # Panel 2 inputs
                                    sidebarPanel(
                                      h3("Inputs"),
                                      # Select type of analysis
                                      # HTML instructions for proper display
                                      div(
                                        # Change variable name
                                        uiOutput("UiType")
                                      ),
                                      # Select standards
                                      uiOutput("UiCtrlsNum"),
                                      # Display boxes
                                      div(style = "display: inline-block;vertical-align:top; width: 150px;",
                                          # Change variable name
                                          uiOutput("UiCtrlsSampleSel")
                                      ),
                                      # Separator
                                      div(style = "display: inline-block;vertical-align:top; width: 50px;", HTML("<br>")),
                                      # Ploidy level box
                                      div(style = "display: inline-block;vertical-align:top; width: 150px;",
                                          uiOutput("UiCtrlsPloNum")
                                      ),
                                      # Perform regression and prediction
                                      HTML("</br>"),
                                      actionButton(inputId = "InReg", label = "Regression"),
                                      # Regression summary
                                      verbatimTextOutput("RegText"),
                                    ),
                                    # Panel 2 outputs
                                    mainPanel(
                                      # Regression plot
                                      h3("Regression"),
                                      verbatimTextOutput("Warn2"),
                                      verbatimTextOutput("Warn3"),
                                      plotOutput("RegPlot"),
                                      # Table results 2
                                      h3("Estimated ploidy"),
                                      DT::dataTableOutput("ResDf2"),
                                    ),
                           ),
                           # Panel 3
                           tabPanel("Summary",
                                    # Panel 3 inputs
                                    fluidPage(
                                      fluidRow(
                                        # Obtain summary
                                        actionButton(inputId = "InSum", label = "Preview"),
                                        # Histogram of all samples
                                        verbatimTextOutput("Warn4"),
                                        verbatimTextOutput("Warn5"),
                                        plotOutput("HisPlotAll"),
                                      ),
                                      HTML("</br>"),
                                      fluidRow(
                                        # Display plot parameters
                                        uiOutput("UiPlotParam"),
                                        # Download plot
                                        downloadButton("downloadPlot", "Save Plot"),
                                      ),
                                      HTML("</br>"),
                                      fluidRow(
                                        # Table results 3
                                        DT::dataTableOutput("ResDf3"),
                                        # Download data
                                        downloadButton("DownloadData", "Save table"),
                                      ),
                                    )
                                    
                           ),
                           # Panel 4
                           tabPanel("Help",
                                    uiOutput('textWithHTML')
                           )
                )
)

# ------------------------------------------------------------------------------------------------------------------- #

# Define server function  
server <- function(input, output, session) {
  # Session
  session$onSessionEnded(function() {
    stopApp()
  })
  # Inputs calculated in server
  # Panel 1 inputs
  # Input files
  output$UiFileUp <- renderUI({
    fileInput(inputId = "InFiles",
              label = "Upload multiple FCS files",
              multiple = TRUE,
              accept = c(".FCS", ".fcs"))
  })
  # Select channel
  output$UiChannelSel <- renderUI({
    selectInput(inputId = "InChan",
                label = "Select a channel",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])),
                selected = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])[1]))
  })
  # Select sample
  output$UiSampleSel <- renderUI({
    selectInput(inputId = "InSample",
                label = "Select a sample",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files)))
  })
  # Adjust smoothing
  output$UiSmoothNum <- renderUI({
    numericInput(inputId = "InSmooth",
                 label = "Adjust smoothing",
                 value = 0.1,
                 min = 0.02,
                 max = 1,
                 step = 0.02)
  })
  # Adjust window
  output$UiWindowNum <- renderUI({
    numericInput(inputId = "InWindow",
                 label = "Adjust window width",
                 value = 50,
                 min = 1,
                 max = 1000,
                 step = 10)
  })
  # Select maximum number of peaks
  output$UiPeaksNum <- renderUI({
    numericInput(inputId = "InMaxPeaks",
                 label = "Select maximun number of peaks to plot",
                 value = 3,
                 min = 1,
                 step = 1)
  })
  # Select minimun number of events
  output$UiEveMin <- renderUI({
    numericInput(inputId = "InMinEve",
                 label = "Select minimun count number for peak",
                 value = 5,
                 min = 0.5,
                 step = 0.5)
  })
  
  # Update UI when when files are uploaded
  observeEvent(input$InFiles, {
    req(InitDf())
    # Display sample names
    output$UiSampleSel <- renderUI({
      selectInput(inputId = "InSample",
                  label = "Select a sample",
                  choices = names(InitDf()$Files))
    })
    # Display channel names
    output$UiChannelSel <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Select a channel",
                  choices = names(InitDf()$Files[[1]]),
                  selected = names(InitDf()$Files[[1]])[1])
    })
  })
  

  # Update UI peaks when new channel is selected
  observeEvent(input$InChan, {
    output$UiPeaksBox <- renderUI({
      checkboxGroupInput(inputId = "InPeaksPlot",
                         label = "Select G1 and G2 peaks",
                         choices = PlotPoints()$MaxIndex,
                         selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
    })
  })
  
  # Display the last used parameters when changing sample 
  # Last used channel
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$UiChannelSel <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Select a channel",
                  choices = names(InitDf()$Files[[1]]),
                  selected = Df$DataPeaks$Channel[SampNum])
    })
  })
  # Last used smoothing
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$UiSmoothNum <- renderUI({
      numericInput(inputId = "InSmooth",
                   label = "Adjust smoothing",
                   value = Df$DataPeaks$Smoothing[SampNum],
                   min = 0.02,
                   max = 1,
                   step = 0.02)
    })
  })
  # Last used window
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$UiWindowNum <- renderUI({
      numericInput(inputId = "InWindow",
                   label = "Adjust window",
                   value = Df$DataPeaks$Window[SampNum],
                   min = 1,
                   max = 1000,
                   step = 10)
    })
  })
  # Last selected peaks
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    LabNum <- match(Df$DataPeaks[SampNum, c(5,6)], PlotPoints()$MaxIndex)
    output$UiPeaksBox <- renderUI({
      checkboxGroupInput(inputId = "InPeaksPlot",
                         label = "Select G1 and G2 peaks",
                         choices = PlotPoints()$MaxIndex,
                         selected = PlotPoints()$MaxIndex[LabNum])
    })
  })
  
  # Panel 2 inputs
  # Select number of standards
  output$UiCtrlsNum <- renderUI({
    numericInput(inputId = "InNumCtrl",
                 label = "Select number of standards",
                 value = ifelse(length(input$InFiles[, 1]) < 4, length(input$InFiles[, 1]), 4),
                 min = 1,
                 step = 1)
  })
  # Create ploidy or genome size type of analysis
  observeEvent(input$InFiles, {
    output$UiType <- renderUI({
      req(input$InNumCtrl)
      radioButtons(inputId = "InType",
                   label = "Select type of analysis",
                   choices = c("Ploidy", "Genome size"),
                   selected = "Ploidy",
                   inline = TRUE)
    })
  })
  # Create select standards
  observeEvent(input$InFiles, {
    output$UiCtrlsSampleSel <- renderUI({
      req(input$InNumCtrl)
      StandardList <- paste("Standard", 1:input$InNumCtrl)
      Ls <- list()
      for (i in 1:length(StandardList)){
        Ls[[i]] <- selectInput(inputId = paste0("InCtrlSample", i),
                               label = StandardList[i],
                               choices = c("", names(InitDf()$Files)),
                               selected = "")
      }
      Ls
    })
  })
  # Create standards ploidy
  observeEvent(input$InFiles, {
    output$UiCtrlsPloNum <- renderUI({
      req(input$InNumCtrl)
      PloidyList <- paste(input$InType, 1:input$InNumCtrl)
      Ls <- list()
      for (i in 1:length(PloidyList)){
        Ls[[i]] <- numericInput(inputId = paste0("InCtrlPlo", i),
                                label = PloidyList[i],
                                value = i,
                                min = 1,
                                step = 1)
      }
      Ls
    })
  })
  
  # Panel 3 inputs
  # Render compilation of all histograms
  output$HisPlotAll <- renderPlot({
    req(Df$Res)
    # Compile and arrange plots
    grid.arrange(grobs = AllPl(), ncol = input$InGrid)
  })
  # Create selection inputs for plot saving
  output$UiPlotParam <- renderUI({
    PlotParamType <- c("Sel", "Sel", "Num", "Num", "Num", "Num")
    PlotParamList <- c("Device", "Units", "Grid", "Width", "Height", "Dpi")
    PlotParamChoi <- list(c("png", "tiff"), c("in", "cm", "mm", "px"), NA, NA, NA, NA)
    PlotParamVal <- c(NA, NA, 4, 300, 200, 300)
    PlotParamStep <- c(NA, NA, 1, 50, 50, 50)
    PlotParamSel <- c("png", "mm")
    Ls <- list()
    for (i in 1:length(PlotParamList)){
      if (PlotParamType[i] == "Sel"){
        Ls[[i]] <- div(
          selectInput(inputId = paste0("In", PlotParamList[i]),
                      label = PlotParamList[i],
                      choices = PlotParamChoi[[i]],
                      selected = PlotParamSel[i]),
          style = "display: inline-block;vertical-align:top; width: 150px;",
        )
      } else {
        Ls[[i]] <- div(
          numericInput(inputId = paste0("In", PlotParamList[i]),
                       label = PlotParamList[i],
                       value = PlotParamVal[i],
                       step = PlotParamStep[i]),
          style = "display: inline-block;vertical-align:top; width: 150px;",
        )
      }
      
    }
    Ls
  })
  # Download figure
  output$downloadPlot <- downloadHandler(
    filename <- function() {
      paste("all_histograms", input$InDevice, sep = ".")
    },
    content <- function(file) {
      ggsave(file, plot = grid.arrange(grobs = AllPl(), ncol = input$InGrid), device = input$InDevice, width = input$InWidth, units = input$InUnits, height = input$InHeight, dpi = input$InDpi)
    }
  )
  
  # Download data
  output$DownloadData <- downloadHandler(
    filename <- function() {
      "peaks.csv"
    },
    content <- function(file) {
      write.csv(Df$Sum, file, row.names = FALSE)
    }
  )
  
  # Define specific functions for peaks calculation
  # Lines
  GetLine <- function(File, ChanNum, Span){
    # Extract and filter expressions
    Exp <- exprs(File[, ChanNum])
    # Estimate the maximum number of bins from the data, and exclude negative values
    MaxBreaks <- length(seq(0, range(File)[,ChanNum][2]))
    # Assumes that lowest possible binning is 256, which allows to exclude low range channels
    MaxBreaks <- ifelse(MaxBreaks < 256, 256, MaxBreaks)
    # Filter low and high values to fit bins in histogram
    Exp <- Exp[Exp >= 1]
    Exp <- Exp[Exp <= MaxBreaks]
    # Calculate histogram
    Hist <- hist(Exp, breaks = 1:MaxBreaks, plot = FALSE)
    # Counts is the y histogram variable
    Counts <- Hist$counts
    # Index is the x histogram variable
    Index <- 1:length(Counts)
    # Fit a line to the data
    Fit <- loess(Counts ~ Index, span = Span)
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create plotting data frame
    data.frame(Fluorescence = Index, Count = FitLine)
  }
  
  # Points
  GetPoints <- function(PlotLine, Width){
    # The window to calculate the maximum
    Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    FlatLine <- rollapply(data = zoo(PlotLine$Count), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    Delta <- FlatLine - PlotLine$Count[-c(1:Width, length(PlotLine$Fluorescence) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    MaxIndex <- which(Delta <= 0) + Width
    # Obtain intensity of the points (y axis)
    Intensity <- PlotLine$Count[MaxIndex]
    # Remove peaks with number of events lower than specified by user
    MaxIndex <- MaxIndex[Intensity > input$InMinEve]
    Intensity <- Intensity[Intensity > input$InMinEve]
    # Create points data frame
    data.frame(MaxIndex = MaxIndex, Intensity = Intensity)
  }
  
  # Peak calculations
  # Get initial files and information
  InitDf <- eventReactive(input$InFiles, {
    # Code used to evaluate run time (START)
    timing_result <- system.time({
      
      FilesLs <- list()
      Channels <- NULL
      Smoothings <- NULL
      Windows <- NULL
      LineLs <- list()
      PointLs <- list()
      G1s <- NULL
      G2s <- NULL
      Names <- NULL
      
      # Read FCS files and calculate initial information
      for(i in 1:length(input$InFiles[, 1])){
        FilesLs[[i]] <- read.FCS(input$InFiles[[i, 'datapath']], emptyValue = FALSE, alter.names = TRUE)
        Name <- sub(" .*", "", sub(".FCS", "", sub(".fcs", "", input$InFiles$name[i])))
        Names <- c(Names, Name)
        Channels[i] <- names(FilesLs[[i]])[1]
        Smoothings[i] <- 0.1
        Windows[i] <- 50
        LineLs[[i]] <- GetLine(File = FilesLs[[i]], ChanNum = 1, Span = 0.1)
        PointLs[[i]] <- GetPoints(PlotLine = LineLs[[i]], Width = 50)
        G1s[i] <- PointLs[[i]]$MaxIndex[1]
        G2s[i] <- PointLs[[i]]$MaxIndex[2]
      }
      
      names(FilesLs) <- Names
    })
    # Code used to evaluate run time (END)
    print("Upload time")
    print(timing_result)
    
    # Return a list with the information
    list(Files = FilesLs, Channels = Channels, Smoothings = Smoothings, Windows = Windows, G1s = G1s, G2s = G2s)
  })
  
  # Create reactive values to modify initially calculated data (It had to be created as reactive value to be modified via proxy)
  Df <- reactiveValues(DataPeaks = NULL)
  # Copy information form initial data frame
  observeEvent(input$InFiles, {
    req(InitDf())
    Df$DataPeaks <- data.frame(Sample = names(InitDf()$Files),
                               Channel = InitDf()$Channels,
                               Smoothing = InitDf()$Smoothings,
                               Window = InitDf()$Windows,
                               G1 = InitDf()$G1s,
                               G2 = InitDf()$G2s)
    
    output$ResDf1 <- DT::renderDataTable(isolate(Df$DataPeaks),
                                         editable = FALSE)
  })
  
  # Create plot line
  PlotLine <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan), {
    req(InitDf())
    # Get sample name
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get file information
    File <- InitDf()$Files[[SampNum]]
    # Get channels information
    ChanNum <- grep(input$InChan, names(File))
    # Return Line
    GetLine(File = File, ChanNum = ChanNum, Span = Df$DataPeaks[SampNum, 3])
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan, input$InMaxPeaks, input$InMinEve), {
    req(InitDf())
    # Get sample name
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get width
    Width <- Df$DataPeaks[SampNum, 4]
    # Return points
    GetPoints(Width = Width, PlotLine = PlotLine())
  })
  # Warning to incorrect execution of regression
  WarnHis1 <- reactive({
    validate(
      need(input$InFiles != "", "Please upload files to obtain histograms"),
    )
  })
  
  output$Warn1 <- renderPrint({
    WarnHis1()
  })
  # Plot data
  output$HisPlot <- renderPlot({
    req(input$InFiles)
    req(InitDf())
    req(input$InSample)
    req(PlotLine())
    req(input$InMaxPeaks)
    req(input$InMinEve)
    
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get name for plotting and table
    Name <- names(InitDf()$Files)[SampNum]
    # Number of labels
    LabNum <- match(input$InPeaksPlot, PlotPoints()$MaxIndex)
    # Labels
    Labs <- paste0("G", 1:length(PlotPoints()$MaxIndex))
    # Plot histogram with points
    ggplot(PlotLine(), aes(x = Fluorescence, y = Count)) +
      geom_line() +
      annotate("point", x = PlotPoints()$MaxIndex[1:input$InMaxPeaks],
               y = PlotPoints()$Intensity[1:input$InMaxPeaks],
               col = "red") +
      annotate("text", x = PlotPoints()$MaxIndex[LabNum],
               y = PlotPoints()$Intensity[LabNum] + 10,
               col = "black",
               label = Labs[1:length(LabNum)]) +
      labs(title = Name) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Proxys for data replacement
  # Channel
  observeEvent(input$InChan, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update channel for all samples
    Df$DataPeaks[, 2] <- input$InChan
    # Recalculate peaks for all samples
    for (i in 1:length(InitDf()$Files)){
      # Get file information
      File <- InitDf()$Files[[i]]
      # Get channels information
      ChanNum <- grep(input$InChan, names(File))
      # Return Line working data frame
      LineWkDf <- GetLine(File = File, ChanNum = ChanNum, Span = Df$DataPeaks[i, 3])
      # Get width
      Width <- Df$DataPeaks[i, 4]
      # Return points
      PointsWkDf <-  GetPoints(Width = Width, PlotLine = LineWkDf)
      # Update peaks
      Df$DataPeaks[i, 5] <- PointsWkDf$MaxIndex[1]
      Df$DataPeaks[i, 6] <- PointsWkDf$MaxIndex[2]
    }
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$DataPeaks)
  })
  # Smoothing
  observeEvent(input$InSmooth, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$DataPeaks[SampNum, 3] <- input$InSmooth
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$DataPeaks)
  })
  # Window
  observeEvent(input$InWindow, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update window
    Df$DataPeaks[SampNum, 4] <- input$InWindow
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$DataPeaks)
  })
  # Peaks
  observeEvent(input$InPeaksPlot, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update G1 and G2 peaks
    Df$DataPeaks[SampNum, 5] <- input$InPeaksPlot[1]
    Df$DataPeaks[SampNum, 6] <- input$InPeaksPlot[2]
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$DataPeaks)
  })
  
  # Regression calculations
  # Warning to incorrect execution of regression
  WarnReg2 <- reactive({
    validate(
      need(input$InFiles != "", "Please upload files to perform regression"),
      need(input$InCtrlSample2 != "", "Please provide at least two standards to perform regression"),
    )
  })
  
  output$Warn2 <- renderPrint({
    WarnReg2()
  })
  
  WarnReg3 <- eventReactive(input$InReg, {
    validate(
      need(input$InCtrlSample2 != "", "Please provide at least two standards to perform regression"),
      need(length(input$InPeaksPlot) != 0, "No peaks detected"),
    )
  })
  
  output$Warn3 <- renderPrint({
    WarnReg3()
  })
  
  # Create linear regression model
  observeEvent(input$InReg, {
    # Code used to evaluate run time (START)
    timing_result <- system.time({
      
      req(input$InCtrlSample2)
      req(input$InPeaksPlot)
      
      # Select required information from DataPeaks
      Df$DataReg <- Df$DataPeaks[,c(1,5,6)]
      # Combine G1 and G2 fluorescent intensity
      Df$DataReg <- gather(data = Df$DataReg,
                           key = "Phase",
                           value = "Intensity",
                           -Sample)
      # Convert intensity to numeric
      Df$DataReg$Intensity <- as.numeric(Df$DataReg$Intensity)
      # Order alphabetically
      Df$DataReg <- Df$DataReg[order(Df$DataReg$Sample),]
      # Initialize variable
      Pattern <- NULL
      # Create control sample strings to evaluate
      CtrlSamples <- paste0("input$InCtrlSample", 1:input$InNumCtrl)
      # Loop over strings
      for (Sample in CtrlSamples){
        # Evaluate each string to obtain the input stored value
        Pattern <- c(Pattern, eval(parse(text = Sample)))
      }
      # Create pattern to search samples
      PatternSearch <- paste(Pattern, collapse = "|")
      # Separate standards and test data frames
      Df$Ctrls <- Df$DataReg[grep(PatternSearch, Df$DataReg$Sample),]
      Df$Tests <- Df$DataReg[-grep(PatternSearch, Df$DataReg$Sample),]
      # Add test and control type information
      Df$Ctrls$Type <- "Standard"
      # To control in case no standards are selected
      if (nrow(Df$Tests) > 0){
        Df$Tests$Type <- "Test"
      }
      # Create control ploidy strings to evaluate
      CtrlPloidies <- paste0("input$InCtrlPlo", 1:input$InNumCtrl)
      # Initialize variable
      AllPloidies <- NULL
      # Loop over strings
      for (Ploidy in CtrlPloidies){
        # Evaluate each string to obtain the input stored value
        AllPloidies <- c(AllPloidies, eval(parse(text = Ploidy)))
      }
      # Change order of all ploidies to match the one in the data frame
      AllPloidies <- AllPloidies[match(unique(Df$Ctrls$Sample), Pattern)]
      # Multiply ploidy for G2 peaks
      AllPloidies <- rep(AllPloidies, each = 2) * 1:2
      # Add ploidy
      Df$Ctrls$Ploidy <- AllPloidies
      # Linear model
      Mod <- lm(Ploidy ~ Intensity, data = Df$Ctrls)
      # Extract adjusted R squared
      Rsquared <- substr(as.character(summary(Mod)$adj.r.squared), 0, 5)
      # Extract p value
      Pval <- pf(summary(Mod)$fstatistic[1], summary(Mod)$fstatistic[2], summary(Mod)$fstatistic[3], lower.tail = FALSE)
      Pval <- formatC(Pval, format = "e", digits = 2)
      # Print summary
      output$RegText <- renderPrint({
        print(summary(Mod))
      })
      # Perform prediction
      Df$Tests$Ploidy <- predict(Mod, Df$Tests) %>% round(2)
      # Combine results
      Df$Res <- rbind(Df$Ctrls, Df$Tests)
      
      # Before rendering data frame, change row names and column order
      rownames(Df$Res) <- 1:nrow(Df$Res)
      Df$Res <- Df$Res[,c(1, 4, 2, 3, 5)]
      Df$ResPrelim <- Df$Res
      # Add size to ploidy column name
      colnames(Df$ResPrelim)[5] <- input$InType
      # Final output data frame
      output$ResDf2 <- DT::renderDataTable(isolate(Df$ResPrelim),
                                           editable = FALSE)
      # Get Rsquared coordinates to plot
      RsquaredXcoord <- Df$Res$Intensity
      RsquaredXcoord <- RsquaredXcoord[!is.na(RsquaredXcoord)]
      RsquaredXcoord <- min(RsquaredXcoord)
      
      RsquaredYcoord <- Df$Res$Ploidy
      RsquaredYcoord <- RsquaredYcoord[!is.na(RsquaredYcoord)]
      RsquaredYcoord <- max(RsquaredYcoord)
      # Plot regression
      output$RegPlot <- renderPlot({
        ggplot(data = Df$Res, mapping = aes(Intensity, Ploidy)) +
          geom_point(color = "red") +
          geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
          geom_text_repel(label = paste(Df$Res$Sample, Df$Res$Phase)) +
          annotate("text", x = RsquaredXcoord + 50, y = RsquaredYcoord - 1, label = paste0("R\u00B2 = ", Rsquared, "\np = ", Pval), size = 5, parse = FALSE) +
          labs(y = input$InType)
      })
    })
    # Code used to evaluate run time (END)
    print("Regression time")
    print(timing_result)
  })
  # Warning to incorrect execution of summary
  WarnSum4 <- reactive({
    validate(
      need(input$InFiles != "", "Please peform regression firts"),
    )
  })
  
  output$Warn4 <- renderPrint({
    WarnSum4()
  })
  
  WarnSum5 <- eventReactive(input$InSum, {
    validate(
      need(Df$Res != "", "No regression found")
    )
  })
  
  output$Warn5 <- renderPrint({
    WarnSum5()
  })
  # Summary
  observeEvent(input$InSum, {
    # Code used to evaluate run time (START)
    timing_result <- system.time({
      
      req(InitDf())
      req(Df$Res)
      # Copy regression results data frame
      Df$Sum <-Df$Res
      # Remove intensity to compare ploidies of different phases
      Df$Sum <- Df$Sum[,-4]
      # Spread data
      Df$Sum <- spread(Df$Sum,
                       key = "Phase",
                       value = "Ploidy")
      # Divide G2 by 2
      Df$Sum$G2 <- round(Df$Sum$G2 / 2, 2)
      # Obtain mean of G1 and G2
      Df$Sum <- mutate(Df$Sum, Mean = round(rowMeans(select(Df$Sum, G1, G2), na.rm = TRUE), 2))
      # Obtain rounded ploidy
      Df$Sum$Rounded <- round(Df$Sum$Mean, 0)
      # Rename phase G1 and G2 in the second data frame
      colnames(Df$Sum)[3] <- paste(input$InType, "G1")
      colnames(Df$Sum)[4] <- paste(input$InType, "G2")
      # Combine with histograms panel data frame
      Df$Sum <- left_join(Df$DataPeaks, Df$Sum)
      # Convert intensity to numeric
      Df$Sum$G1 <- as.numeric(Df$Sum$G1)
      Df$Sum$G2 <- as.numeric(Df$Sum$G2)
      # Rename phase G1 and G2 in the summary data frame
      colnames(Df$Sum)[5] <- "Intensity G1"
      colnames(Df$Sum)[6] <- "Intensity G2"
      # Rearrange data frame order and remove ploidy inferred from G1 and G2 peaks (columns 8 and 9)
      Df$Sum <- Df$Sum[,c(1,7,2:6,10,11)]
      # Final output data frame
      output$ResDf3 <- DT::renderDataTable(isolate(Df$Sum),
                                           editable = FALSE)
    })
    # Code used to evaluate run time (END)
    print("Peaks calculation time")
    print(timing_result)
  })
  
  # Plot all histograms
  AllPl <- eventReactive(input$InSum, {
    # Code used to evaluate run time (START)
    timing_result <- system.time({
      # Create empty list to store plots
      PlotLs <- list()
      # Loop over the samples and create lines
      for (i in 1:length(InitDf()$Files)){
        # Create plot line
        # Get file information
        File <- InitDf()$Files[[i]]
        # Get channels information
        ChanNum <- grep(input$InChan, names(File))
        # Return Line working data frame
        LineWkDf <- GetLine(File = File, ChanNum = ChanNum, Span = Df$DataPeaks[i, 3])
        # Create plot points
        # Get width
        Width <- Df$DataPeaks[i, 4]
        # Return points
        PointsWkDf <-  GetPoints(Width = Width, PlotLine = LineWkDf)
        # Get name for plotting and table
        Name <- names(InitDf()$Files)[i]
        # Number of labels but only the peaks selected in Panel 1
        LabNum <- match(Df$DataPeaks[i, c(5,6)], PointsWkDf$MaxIndex)
        # Labels
        Labs <- paste0("G", 1:length(PointsWkDf$MaxIndex))
        # Plot histogram with points
        Pl <- ggplot(LineWkDf, aes(x = Fluorescence, y = Count)) +
          geom_line() +
          annotate("point", x = PointsWkDf$MaxIndex[LabNum],
                   y = PointsWkDf$Intensity[LabNum],
                   col = "red") +
          annotate("text", x = PointsWkDf$MaxIndex[LabNum],
                   y = PointsWkDf$Intensity[LabNum] + 10,
                   col = "black",
                   label = Labs[1:length(LabNum)]) +
          labs(title = Name) +
          theme(plot.title = element_text(hjust = 0.5))
        # Save plot in list
        PlotLs[[i]] <- Pl
      }
    })
    # Code used to evaluate run time (END)
    print("Plots generation time")
    print(timing_result)
    
    # Return list of plots
    PlotLs
  })
  
  # Panel 4 help
  output$textWithHTML <- renderUI({
    HTML(markdown::markdownToHTML('help.md', template = FALSE))
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)