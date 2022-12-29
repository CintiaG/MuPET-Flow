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

# Pending
# las modification of name sape
# orde peraks before saving?
# The number of peaks behaviour is weird, it is not updated

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
                                      # Download data
                                      downloadButton("DownloadData", "Download"),
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
                                      #  Table results 2
                                      h3("Estimated ploidy"),
                                      DT::dataTableOutput("ResDf2"),
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
                 max = 5,
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
      write.csv(Df$data, file, row.names = FALSE)
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
                  selected = Df$data$Channel[SampNum])
    })
  })
  # Last used smoothing
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$UiSmoothNum <- renderUI({
      numericInput(inputId = "InSmooth",
                   label = "Adjust smoothing",
                   value = Df$data$Smoothing[SampNum],
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
                   value = Df$data$Window[SampNum],
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
                              choices = names(InitDf()$Files))
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
                               value = 3,
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
  
  # Create reactive values
  # Get initial files and information
  # Read FCS files and calculate initial information
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

  # Create parameter and results data frame to modify initially calculated data
  # In this case, it had to be created and reactive value in order to be modified via proxy
  Df <- reactiveValues(data = NULL)
  
  observeEvent(input$InFiles, {
    req(InitDf())
    Df$data <- data.frame(Sample = names(InitDf()$Files),
                          Channel = InitDf()$Channels,
                          Smoothing = InitDf()$Smoothings,
                          Window = InitDf()$Windows,
                          G1 = InitDf()$G1s,
                          G2 = InitDf()$G2s)
    
    output$ResDf1 <- DT::renderDataTable(isolate(Df$data),
                                        editable = FALSE)
  })
  
  #Df3 <- reactiveValues(data = NULL)
  # observe more things? correct things?
  # init Df2 that updates with df?
  #observeEvent(c(input$InFiles, input$InPeaksPlot, input$InNumCtrl), {
  #  req(InitDf())
    #Df$data2 <- Df$data
    
  #  Df$data2 <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data
    # Combine G1 and G2 fluorescent intensity
    #Df2$data <- gather(data = Df2$data,
  #  Df$data2 <- gather(data = Df$data2,
  #                key = "Phase",
  #                value = "Intensity",
  #                -Sample)
    # Add ploidy
  #  Df$data2$Ploidy <- NA#input$InNumCtrl
    
    # Order alphabetically
  #  Df$data2 <- Df$data2[order(Df$data2$Sample),]
  #  output$ResDf2 <- DT::renderDataTable(isolate(Df$data2),
  #                                       editable = FALSE)
  #})
  
  # Try to select dataframe per sample
  #observeEvent(input$CtrlPlo, {
  # What to observe InNumCtrl
  #MyThingsToObs <- "casa"#paste0("input$InCtrlSample", 1:input$InNumCtrl)
  
  #observeEvent(c(input$InFiles, input$InCtrlSample1), {
  
  
  observeEvent(input$InReg, {
  #observeEvent(eval(parse(text = paste("c(", paste0("input$InCtrlSample", 1:input$InNumCtrl, collapse = ", "), ")"))), {
    req(input$InCtrlSample1)
    
    # Select required information from Df
    Df$data2 <- Df$data[,c(1,5,6)]
    
    # Combine G1 and G2 fluorescent intensity
    Df$data2 <- gather(data = Df$data2,
                       key = "Phase",
                       value = "Intensity",
                       -Sample)
    
    # Order alphabetically
    Df$data2 <- Df$data2[order(Df$data2$Sample),]
    # Rename rownames 
    #rownames(Df$data2) <- 1:(input$InNumCtrl * 2)
    
    
    
    Pattern <- NULL
    #Sliders <- grep("slider", names(input), value = TRUE)
    #for (Slider in Sliders){
    #  Pattern <- input$Slider
    #}
    #paste0("input$slider", 1:input$InNumCtrl)
    #input$slider1
    #Pattern <- paste(names(input)[grep("slider", names(input))], sep = "|")
    #Pattern
    # Change to sliders name?
    CtrlSamples <- paste0("input$InCtrlSample", 1:input$InNumCtrl)
    #eval(parse(text = Sliders))
    for (Sample in CtrlSamples){
      Pattern <- c(Pattern, eval(parse(text = Sample)))
    }
    
    Pattern2 <- Pattern
    
    Pattern <- paste(Pattern, collapse = "|")
    
    # Here I would have to seprate the dataframes first
    # in orifinal stran patter ins control, kkepp?
    #Df$Ctrls <- Df$data2[grep(Pattern, Df$data2$Strain),]
    Df$Ctrls <- Df$data2[grep(Pattern, Df$data2$Sample),]
    Df$Tests <- Df$data2[-grep(Pattern, Df$data2$Sample),]
    
    # Add test and control type information
    Df$Ctrls$Type <- "Control"
    # To control in case you select no controls, but same for test?
    if (nrow(Df$Tests) > 0){
      Df$Tests$Type <- "Test"
    }
    #Df$Tests$Type <- "Test"
    
   
    
    Ploidies <- paste0("input$InCtrlPlo", 1:input$InNumCtrl)
    
    AllPloidies <- NULL
    
    # Same but to obtain ploidues
    for (Ploidy in Ploidies){
      AllPloidies <- c(AllPloidies, eval(parse(text = Ploidy)))
    }
    
    # Change order of all ploidies to match the one in the data frame
    
    Idx <- match(unique(Df$Ctrls$Sample), Pattern2)
    
    AllPloidies <- AllPloidies[Idx]
    
    AllPloidies <- rep(AllPloidies, each = 2) * 1:2#NA#input$InNumCtrl
    
    
    
    # Add ploidy
    Df$Ctrls$Ploidy <- AllPloidies
    Df$Ctrls$Ploidy2 <- AllPloidies
    # maybe above
    Df$Ctrls$Intensity <- as.numeric(Df$Ctrls$Intensity)
    Df$Tests$Intensity <- as.numeric(Df$Tests$Intensity)
    
    
    
    # Linear model
    Mod <- lm(Ploidy ~ Intensity, data = Df$Ctrls)
    
    # Summary
    # render text?
    #summary(Mod)
    
    # Prediction
    Df$Tests$Ploidy <- predict(Mod, Df$Tests) %>% round(2) # round can be changed
    Df$Tests$Ploidy2 <- predict(Mod, Df$Tests) %>% round(0) # round can be changed
    
    # Combine results
    Df$Res <- rbind(Df$Ctrls, Df$Tests)
    
    # Before rendering data frame, change rownames
    rownames(Df$Res) <- 1:nrow(Df$Res)
    
    # Final output data frame
    output$ResDf2 <- DT::renderDataTable(isolate(Df$Res),
                                         editable = FALSE)
    
    output$RegPlot <- renderPlot({
      ggplot(data = Df$Res, mapping = aes(Intensity, Ploidy)) +
        geom_point(color = "red") +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
      #geom_point(data = DfTests, color = "blue") +
      #geom_text(label = DfRes$Strain)
      #geom_text_repel(label = DfRes$Strain)
        geom_text_repel(label = paste(Df$Res$Sample, Df$Res$Phase))
    })
    
    #Df$Ctrls <- Df$data[grep(Pattern, Df$data$Sample),]
    
    
    #output$ResDf3 <- DT::renderDataTable(isolate(Df$Ctrls),
    #                                     editable = FALSE)
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
    GetLine(File = File, ChanNum = ChanNum, Span = Df$data[SampNum, 3])
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan, input$InMaxPeaks), {
    req(InitDf())
    # Get sample name
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get width
    Width <- Df$data[SampNum, 4]
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
    #Name <- names(InitDf()$Files)
    Name <- names(InitDf()$Files)[SampNum]
    # Number of labels
    LabNum <- match(input$InPeaksPlot, PlotPoints()$MaxIndex)
    Labs <- c("G1", "G2", "G3", "G4", "G5")
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
  observeEvent(input$InChan, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$data[SampNum, 2] <- input$InChan
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$data)
  })
  # Smoothing
  observeEvent(input$InSmooth, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$data[SampNum, 3] <- input$InSmooth
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$data)
  })
  # Window
  observeEvent(input$InWindow, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update window
    Df$data[SampNum, 4] <- input$InWindow
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$data)
  })
  # Peaks
  observeEvent(input$InPeaksPlot, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update G1 and G2 peaks
    Df$data[SampNum, 5] <- input$InPeaksPlot[1]
    Df$data[SampNum, 6] <- input$InPeaksPlot[2]
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf1')
    DT::replaceData(Proxy, Df$data)
  })
  # Controls ploidy
  # is this working?
  observeEvent(input$InCtrlPlo1, {
#    req(Df2)
#    Df2$Ploidy <- input$InNumCtrl
    #SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$data2$Ploidy <- input$InCtrlPlo1
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf2')
    DT::replaceData(Proxy, Df$data2)
  })
  
  
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)