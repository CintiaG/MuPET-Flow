####################################
# Flow Peaks                       #
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

# Define UI
ui <- fluidPage(theme = shinytheme("united"),
                # App name
                navbarPage("Flow Peaks",
                           # Panel 1
                           tabPanel("Peaks",
                                    # Panel 1 inputs
                                    sidebarPanel(
                                      h3("Inputs"),
                                      # Input files
                                      uiOutput("UiFileUp"),
                                      # Select sample
                                      uiOutput("UiSampleSel"),
                                      # Select channel
                                      uiOutput("UiChannelSel"),
                                      # Adjust smoothing
                                      uiOutput("UiSmoothNum"),
                                      # Adjust window
                                      uiOutput("UiWindowNum"),
                                      # Select maximum number of peaks
                                      uiOutput("UiPeaksNum"),
                                      # Select calculated peaks to plot
                                      uiOutput("UiPeaksBox"),
                                    ),
                                    # Panel 1 outputs
                                    mainPanel(
                                      # Histogram plot
                                      h3("Histogram"),
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
                                      # Select controls
                                      uiOutput("UiCtrlsNum"),
                                      # Display boxes
                                      # HTML instructions for proper display
                                      div(style = "display: inline-block;vertical-align:top; width: 150px;",
                                          # Change variable name
                                          uiOutput("UiCtrlsSampleSel")
                                      ),
                                      # Separator
                                      div(style = "display: inline-block;vertical-align:top; width: 50px;", HTML("<br>")),
                                      # Ploidy level box (Change variable name)
                                      div(style = "display: inline-block;vertical-align:top; width: 150px;",
                                          uiOutput("UiCtrlsPloNum")
                                      ),
                                      # Perform regression and prediction
                                      actionButton(inputId = "InReg", label = "Regression"),
                                    ),
                                    # Panel 2 outputs
                                    mainPanel(
                                      # Regression plot
                                      h3("Regression"),
                                      plotOutput("RegPlot"),
                                      # Table results 2
                                      h3("Estimated ploidy"),
                                      DT::dataTableOutput("ResDf2"),
                                    ),
                           ),
                           # Panel 4
                           tabPanel("Summary",
                                    sidebarPanel(
                                      h3("Inputs?"),
                                      # Obtain summary
                                      actionButton(inputId = "InSum", label = "Summary"),
                                      # Download data
                                      downloadButton("DownloadData", "Download"),
                                    ),
                                    mainPanel(
                                      # Test test (remove)
                                      textOutput("selected_var"),
                                      # Histogram of all samples
                                      plotOutput("HisPlotAll"),
                                      # Table results 3
                                      DT::dataTableOutput("ResDf3"),
                                    ),
                           ),
                           # Panel 3
                           tabPanel("Help",
                                    "In construction..."
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
              label = "Upload FCS files",
              multiple = TRUE,
              accept = c(".FCS", ".fcs"))
  })
  # Select sample
  output$UiSampleSel <- renderUI({
    selectInput(inputId = "InSample",
                label = "Select a sample",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files)))
  })
  # Select channel
  output$UiChannelSel <- renderUI({
    selectInput(inputId = "InChan",
                label = "Select a channel",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])),
                selected = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])[20]))
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
                 label = "Adjust window",
                 value = 10,
                 min = 1,
                 max = 100,
                 step = 1)
  })
  # Select maximum number of peaks
  output$UiPeaksNum <- renderUI({
    numericInput(inputId = "InMaxPeaks",
                 label = "Select maximun number of peaks",
                 value = 3,
                 min = 1,
                 step = 1)
  })
  # Select calculated peaks to plot
  output$UiPeaksBox <- renderUI({
    checkboxGroupInput(inputId = "InPeaksPlot",
                       label = "Select G1 and G2 peaks",
                       choices = PlotPoints()$MaxIndex,
                       selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
  })
  # Download data
  output$DownloadData <- downloadHandler(
    filename <- function() {
      "peaks.csv"
    },
    content <- function(file) {
      write.csv(Df$Sum, file, row.names = FALSE)
    }
  )
  
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
                  selected = names(InitDf()$Files[[1]])[20])
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
                   max = 100,
                   step = 1)
    })
  })
  
  # Panel 2 inputs
  # Select number of controls
  output$UiCtrlsNum <- renderUI({
    numericInput(inputId = "InNumCtrl",
                 label = "Number of controls",
                 value = ifelse(length(input$InFiles[, 1]) < 5, length(input$InFiles[, 1]), 5),
                 min = 1,
                 step = 1)
  })
  # Create select controls
  observeEvent(input$InFiles, {
    req(InitDf())
    output$UiCtrlsSampleSel <- renderUI({
      ControlsList <- paste("Control", 1:input$InNumCtrl)
      Ls <- list()
      for (i in 1:length(ControlsList)){
        Ls[[i]] <- selectInput(inputId = paste0("InCtrlSample", i),
                              label = ControlsList[i],
                              choices = c("", names(InitDf()$Files)),
                              selected = "")
      }
      Ls
    })
  })
  
  # Create controls ploidy
  observeEvent(input$InFiles, {
    req(InitDf())
    output$UiCtrlsPloNum <- renderUI({
      PloidyList <- paste("Ploidy", 1:input$InNumCtrl)
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

  # Define specific functions for peaks calculation
  # Lines
  GetLine <- function(File, ChanNum, Span){
    # Extract and filter expressions
    Exp <- exprs(File[, ChanNum])
    Exp <- Exp[Exp >= 1]
    # Calculate histogram
    Hist <- hist(Exp, breaks = 1:1000, plot = FALSE)
    # Counts is the y histogram variable
    Counts <- Hist$counts
    # Index is the x histogram variable
    Index <- 1:length(Counts)
    # Fit a line to the data
    Fit <- loess(Counts ~ Index, span = Span)
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create plotting data frame
    data.frame(Fluorescence = Index, Freq = FitLine)
  }
  
  # Points
  GetPoints <- function(PlotLine, Width){
    # The window to calculate the maximum
    Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    FlatLine <- rollapply(data = zoo(PlotLine$Freq), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    Delta <- FlatLine - PlotLine$Freq[-c(1:Width, length(PlotLine$Fluorescence) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    MaxIndex <- which(Delta <= 0) + Width
    # Obtain intensity of the points (y axis)
    Intensity <- PlotLine$Freq[MaxIndex]
    # Remove peaks with intensity lower than 10
    MaxIndex <- MaxIndex[Intensity > 10]
    Intensity <- Intensity[Intensity > 10]
    # Create points data frame
    data.frame(MaxIndex = MaxIndex, Intensity = Intensity)
  }
  
  # Peak calculations
  # Get initial files and information
  InitDf <- eventReactive(input$InFiles, {
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
      Name <- sub(" .*", "", input$InFiles$name[i])
      Names <- c(Names, Name)
      Channels[i] <- names(FilesLs[[i]])[20]
      Smoothings[i] <- 0.1
      Windows[i] <- 10
      LineLs[[i]] <- GetLine(File = FilesLs[[i]], ChanNum = 20, Span = 0.1)
      PointLs[[i]] <- GetPoints(PlotLine = LineLs[[i]], Width = 10)
      G1s[i] <- PointLs[[i]]$MaxIndex[1]
      G2s[i] <- PointLs[[i]]$MaxIndex[2]
    }
    
    names(FilesLs) <- Names
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
  PlotPoints <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan, input$InMaxPeaks), {
    req(InitDf())
    # Get sample name
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get width
    Width <- Df$DataPeaks[SampNum, 4]
    # Return points
    GetPoints(Width = Width, PlotLine = PlotLine())
  })
  
  # Plot data
  output$HisPlot <- renderPlot({
    req(input$InFiles)
    req(InitDf())
    req(input$InSample)
    req(PlotLine())
    req(input$InMaxPeaks)
    req(input$InPeaksPlot)
    
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get name for plotting and table
    Name <- names(InitDf()$Files)[SampNum]
    # Number of labels
    LabNum <- match(input$InPeaksPlot, PlotPoints()$MaxIndex)
    # Labels
    Labs <- paste0("G", 1:length(PlotPoints()$MaxIndex))
    # Plot histogram with points
    ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
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
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update channel
    Df$DataPeaks[SampNum, 2] <- input$InChan
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
  observeEvent(input$InReg, {
    req(input$InCtrlSample1)
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
    # Create pattern to serach samples
    PatternSearch <- paste(Pattern, collapse = "|")
    # Seprate controls and test data frames
    Df$Ctrls <- Df$DataReg[grep(PatternSearch, Df$DataReg$Sample),]
    Df$Tests <- Df$DataReg[-grep(PatternSearch, Df$DataReg$Sample),]
    # Add test and control type information
    Df$Ctrls$Type <- "Control"
    # To control in case no controls are selected
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
    
    # Summary
    # render text? or plot in graph???
    #summary(Mod)
    
    # Perform prediction
    Df$Tests$Ploidy <- predict(Mod, Df$Tests) %>% round(2)
    # Combine results
    Df$Res <- rbind(Df$Ctrls, Df$Tests)
    
    # Before rendering data frame, change row names and column order
    rownames(Df$Res) <- 1:nrow(Df$Res)
    Df$Res <- Df$Res[,c(1, 4, 2, 3, 5)]
    # Final output data frame
    output$ResDf2 <- DT::renderDataTable(isolate(Df$Res),
                                         editable = FALSE)
    # Plot regression
    output$RegPlot <- renderPlot({
      ggplot(data = Df$Res, mapping = aes(Intensity, Ploidy)) +
        geom_point(color = "red") +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
        geom_text_repel(label = paste(Df$Res$Sample, Df$Res$Phase))
    })
  })
  
  # Summary
  observeEvent(input$InSum, {
    # Plot
    output$HisPlotAll <- renderPlot({
      # Test if these are necessary
      req(input$InFiles)
      req(InitDf())
      req(input$InSample)
      req(PlotLine())
      req(input$InMaxPeaks)
      req(input$InPeaksPlot)
      
      # Create a mock loop
      
      Things <- c("hose", "door", "table")
      
      Things_to_render <- NULL
      
      for (Element in Things){
        Things_to_render <- c(Things_to_render, Element)
      }
      
      output$selected_var <- renderText({ 
        Things_to_render
      })
      
      # Create plot line
      #PlotLine <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan), {
      #  req(InitDf())
      #  # Get sample name
      #  SampNum <- grep(input$InSample, names(InitDf()$Files))
      #  # Get file information
      #  File <- InitDf()$Files[[SampNum]]
      #  # Get channels information
      #  ChanNum <- grep(input$InChan, names(File))
      #  # Return Line
      #  GetLine(File = File, ChanNum = ChanNum, Span = Df$DataPeaks[SampNum, 3])
      #})
      
      # Create plot points
      #PlotPoints <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan, input$InMaxPeaks), {
      #  req(InitDf())
      #  # Get sample name
      #  SampNum <- grep(input$InSample, names(InitDf()$Files))
      #  # Get width
      #  Width <- Df$DataPeaks[SampNum, 4]
      #  # Return points
      #  GetPoints(Width = Width, PlotLine = PlotLine())
      #})
      
      
      
      
      SampNum <- grep(input$InSample, names(InitDf()$Files))
      # Get name for plotting and table
      Name <- names(InitDf()$Files)[SampNum]
      # Number of labels
      LabNum <- match(input$InPeaksPlot, PlotPoints()$MaxIndex)
      # Labels
      Labs <- paste0("G", 1:length(PlotPoints()$MaxIndex))
      # Plot histogram with points
      ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
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
    # Rename phase G1 and G2 in the second data 
    colnames(Df$Sum)[3] <- "Ploidy G1"
    colnames(Df$Sum)[4] <- "Ploidy G2"
    # Combine with histograms panel data frame
    Df$Sum <- left_join(Df$DataPeaks, Df$Sum)
    # Final output data frame
    output$ResDf3 <- DT::renderDataTable(isolate(Df$Sum),
                                         editable = FALSE)
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)


# Pending
# las modification of name sape
# orde peraks before saving?
# The number of peaks behaviour is weird, it is not updated
# add summary or R2 in plot? with p?
# Obtain average of both peaks and give a final result?
# Help section
# Dowload final data perhaps only and not the peaks?
# It forgets the last seletec peaks