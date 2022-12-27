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

# Define UI
ui <- fluidPage(theme = shinytheme("united"),
                navbarPage("Flow Peaks",
                           tabPanel("Peaks",
                                    sidebarPanel(
                                      # Inputs
                                      tags$h3("Input:"),
                                      # Input files
                                      # Perhaps for style, calculte this in server aswell
                                      fileInput(inputId = "InFiles",
                                                label = "Upload FCS files",
                                                multiple = TRUE,
                                                accept = c(".FCS", ".fcs")),
                                      # Select input file
                                      uiOutput("FileDropdown"),
                                      # Select channel
                                      uiOutput("ChanDropdown"),
                                      # Select smoothing
                                      uiOutput("SmoothSel"),
                                      # Select window
                                      uiOutput("WindowSel"),
                                      # Select maximum number of peaks
                                      numericInput(inputId = "MaxPeaks",
                                                   label = "Select maximun number of peaks",
                                                   value = 3,
                                                   min = 1,
                                                   max = 5,
                                                   step = 1),
                                      # Generate select G1 an G2 peaks selection box
                                      uiOutput("InBox"),
                                      # Save data
                                      downloadButton("DownloadData", "Download"),
                                    ),
                                    mainPanel(
                                      # Outputs
                                      h3("Histogram"),
                                      # Plot histogram
                                      plotOutput("HisPlot"),
                                      h3("Estimated peaks and used parameters"),
                                      # Display table result
                                      DT::dataTableOutput("ResDf"),
                                    )
                           ),
                           tabPanel("Regression",
                                    sidebarPanel(
                                      # Select controls
                                      uiOutput("CtrlSel"),
                                      # Display boxes 
                                      div(style="display: inline-block;vertical-align:top; width: 150px;",
                                          # Change nvariable name
                                          uiOutput("myboxes")
                                      ),
                                      # Separator
                                      div(style="display: inline-block;vertical-align:top; width: 50px;", HTML("<br>")),
                                      # Ploidy level box (Change variable name)
                                      div(style="display: inline-block;vertical-align:top; width: 150px;",
                                          uiOutput("myboxes2")
                                      ),
                                    ),
                                    mainPanel(
                                      #"Some tables and plots"
                                      # Remove selected var, to prove variables
                                      textOutput("selected_var"),
                                      h3("Controls"),
                                      DT::dataTableOutput("ResDf2"),
                                      DT::dataTableOutput("ResDf3"),
                                      #           DT::dataTableOutput("ResDf")
                                    ),
                           ),
                           #"Estimate genome size based on known controls (in construction)."),
                           tabPanel("Help", "In construction...")
                )
)

# ------------------------------------------------------------------------------------------------------------------- #

# Define server function  
server <- function(input, output, session) {
  # Session
  session$onSessionEnded(function() {
    stopApp()
  })
  # Inputs operations calculated in server
  # Panel 1
  # Sample
  output$FileDropdown <- renderUI({
    selectInput(inputId = "InSample",
                label = "Select a sample",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files)))
  })
  # Channels
  output$ChanDropdown <- renderUI({
    selectInput(inputId = "InChan",
                label = "Select a channel",
                choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])),
                selected = ifelse(is.na(input$InFiles), "", names(InitDf()$Files[[1]])[20]))
  })
  # Smoothing
  output$SmoothSel <- renderUI({
    numericInput(inputId = "InSmooth",
               label = "Adjust smoothing",
               value = 0.1,
               min = 0.02,
               max = 1,
               step = 0.02)
  })
  # Window
  output$WindowSel <- renderUI({
    numericInput(inputId = "InWindow",
                label = "Adjust window",
                value = 10,
                min = 1,
                max = 100,
                step = 1)
  })
  # Box peaks
  output$InBox <- renderUI({
    checkboxGroupInput(inputId = "InBox",
                       label = "Select G1 and G2 peaks",
                       choices = PlotPoints()$MaxIndex,
                       selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
  })
  
  # Update UI when when files are uploaded
  observeEvent(input$InFiles, {
    req(InitDf())
    # Displays file names
    output$FileDropdown <- renderUI({
      selectInput(inputId = "InSample",
                  label = "Select a sample",
                  choices = names(InitDf()$Files))
    })
    # Display channel names
    output$ChanDropdown <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Select a channel",
                  choices = names(InitDf()$Files[[1]]),
                  selected = names(InitDf()$Files[[1]])[20])
    })
  })
  
  # Display the last used parameters
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$ChanDropdown <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Select a channel",
                  choices = names(InitDf()$Files[[1]]),
                  selected = Df$data$Channel[SampNum])
    })
  })
  
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$SmoothSel <- renderUI({
      numericInput(inputId = "InSmooth",
                   label = "Adjust smoothing",
                   value = Df$data$Smoothing[SampNum],
                   min = 0.02,
                   max = 1,
                   step = 0.02)
    })
  })
  
  observeEvent(input$InSample, {
    req(InitDf())
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    output$WindowSel <- renderUI({
      numericInput(inputId = "InWindow",
                   label = "Adjust window",
                   value = Df$data$Window[SampNum],
                   min = 1,
                   max = 100,
                   step = 1)
    })
  })
  
  # Panel 2
  # Number of controls
  output$CtrlSel <- renderUI({
    numericInput(inputId = "InCtrl",
                 label = "Number of controls",
                 #value = 5,
                 #value = ifelse(length(input$InFiles) < 5, length(input$InFiles), 5),
                 value = ifelse(length(input$InFiles[, 1]) < 5, length(input$InFiles[, 1]), 5),
                 min = 1,
                 step = 1)
  })
  
  # Add multiple boxes
  #observeEvent(input$InSample, {
  observeEvent(input$InFiles, {
    req(InitDf())
#    BoxesList <- paste("Control", 1:input$InCtrl)
#    v <- list()
#    for (i in 1:length(BoxesList)){
#      v[[i]] <- #box(#width = 3,
#        #title = h4(BoxesList[i]),
#        selectInput(inputId = paste0("slider", i),
#                    label = BoxesList[i],
#                    #choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files)))
#                    choices = names(InitDf()$Files))
#                    ##choices = names(InitDf()$Files))
#                    #choices = list("Not good", "average" , "good"))
#      #)
#    }
    
    #output$myboxes <- renderUI(v)
    output$myboxes <- renderUI({
      BoxesList <- paste("Control", 1:input$InCtrl)
      v <- list()
      for (i in 1:length(BoxesList)){
        v[[i]] <- #box(#width = 3,
          #title = h4(BoxesList[i]),
          selectInput(inputId = paste0("slider", i),
                      label = BoxesList[i],
                      #choices = ifelse(is.na(input$InFiles), "", names(InitDf()$Files)))
                      choices = names(InitDf()$Files))
        ##choices = names(InitDf()$Files))
        #choices = list("Not good", "average" , "good"))
        #)
      }
      v
    })
  })
  
  # For controls ploidy
  observeEvent(input$InFiles, {
    req(InitDf())
    
    #output$myboxes <- renderUI(v)
    output$myboxes2 <- renderUI({
      BoxesList <- paste("Ploidy", 1:input$InCtrl)
      v <- list()
      for (i in 1:length(BoxesList)){
        v[[i]] <- #box(#width = 3,
          #title = h4(BoxesList[i]),
          numericInput(inputId = paste0("CtrlPlo", i),
                       label = BoxesList[i],
                       value = 3,
                       min = 1,
                       #max = 5,
                       step = 1)
      }
      v
    })
  })
  
  
  # Update
#  observeEvent(input$InFiles, {
#    req(InitDf())
#    # Displays file names
#    output$CtrlSel <- renderUI({
#      numericInput(inputId = "InCtrl",
#                   label = "Number of controls",
#                   value = 5,
#                   #value = ifelse(length(input$InFiles) < 5, length(input$InFiles), 5),
#                   #value = ifelse(length(input$InFiles[, 1]) < 5, length(input$InFiles[, 1]), 5),
#                   min = 1,
#                   step = 1)
#    })
#  })

  # Define specific functions
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
      #Name <- sub(" .*", " ", input$InFiles$name[i])
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
  Df <- reactiveValues(data = NULL)
  # In this case, it had to be created and reactive value in order to be modified via proxy
  observeEvent(input$InFiles, {
    req(InitDf())
    Df$data <- data.frame(Sample = names(InitDf()$Files),
                          Channel = InitDf()$Channels,
                          Smoothing = InitDf()$Smoothings,
                          Window = InitDf()$Windows,
                          G1 = InitDf()$G1s,
                          G2 = InitDf()$G2s)
    
    output$ResDf <- DT::renderDataTable(isolate(Df$data),
                                        editable = FALSE)
  #})
  
  #Df2 <- reactiveValues(data = NULL)
  #Df2 <- reactiveValues()
  
  #observeEvent(input$InFiles, {
  #  req(Df)
    #req(input$InCtrl)
    # Generate new data frame for regression calculation
    # Does this have to be together?
   # Df2 <- Df
    #Df2 <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data
    # Combine G1 and G2 fluorescent intensity
    #Df2$data <- gather(data = Df2$data,
    #Df2 <- gather(data = Df2,
    #              key = "Phase",
    #              value = "Intensity",
    #              -Sample)
    # Add ploidy
    #Df2$Ploidy <- NA#input$InCtrl
    #Df2$Ploidy <- 1:input$InCtrl
    
    # Order alphabetically
    #Df2 <- Df2[order(Df2$Sample),]
    
    
    
    #output$ResDf2 <- DT::renderDataTable(isolate(Df2),
    #                                     editable = FALSE)
  })
  
  #Df3 <- reactiveValues(data = NULL)
  # observe more things? correct things?
  # init Df2 that updates with df?
  observeEvent(c(input$InFiles, input$InBox, input$InCtrl), {
    req(InitDf())
    #Df$data2 <- Df$data
    
    Df$data2 <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data[,c(1,5,6)]
    #Df2$data <- Df$data
    # Combine G1 and G2 fluorescent intensity
    #Df2$data <- gather(data = Df2$data,
    Df$data2 <- gather(data = Df$data2,
                  key = "Phase",
                  value = "Intensity",
                  -Sample)
    # Add ploidy
    Df$data2$Ploidy <- NA#input$InCtrl
    
    # Order alphabetically
    Df$data2 <- Df$data2[order(Df$data2$Sample),]
    output$ResDf3 <- DT::renderDataTable(isolate(Df$data2),
                                         editable = FALSE)
  })
  
  # Try to select dataframe per sample
  #observeEvent(input$CtrlPlo, {
  observeEvent(c(input$InFiles, input$slider1), {
    req(input$slider1)
    #Df$Ctrls <- Df$data[grep(input$slider1, Df$data$Sample),]
    output$selected_var <- renderText({ 
      #input$slider1
      grep(input$slider1, Df$data$Sample)
    })
    
    
    
    Df$Ctrls <- Df$data[grep(input$slider1, Df$data$Sample),]
    output$ResDf2 <- DT::renderDataTable(isolate(Df$Ctrls),
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
    GetLine(File = File, ChanNum = ChanNum, Span = Df$data[SampNum, 3])
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$InSample, input$InSmooth, input$InWindow, input$InChan, input$MaxPeaks), {
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
    req(input$MaxPeaks)
    req(input$InBox)
    
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Get name for plotting and table
    #Name <- names(InitDf()$Files)
    Name <- names(InitDf()$Files)[SampNum]
    # Number of labels
    LabNum <- match(input$InBox, PlotPoints()$MaxIndex)
    Labs <- c("G1", "G2", "G3", "G4", "G5")
    # Plot histogram with points
    ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
      geom_line() +
      annotate("point", x = PlotPoints()$MaxIndex[1:input$MaxPeaks],
               y = PlotPoints()$Intensity[1:input$MaxPeaks],
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
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Smoothing
  observeEvent(input$InSmooth, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$data[SampNum, 3] <- input$InSmooth
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Window
  observeEvent(input$InWindow, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update window
    Df$data[SampNum, 4] <- input$InWindow
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Peaks
  observeEvent(input$InBox, {
    SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update G1 and G2 peaks
    Df$data[SampNum, 5] <- input$InBox[1]
    Df$data[SampNum, 6] <- input$InBox[2]
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Controls ploidy
  observeEvent(input$CtrlPlo1, {
#    req(Df2)
#    Df2$Ploidy <- input$InCtrl
    #SampNum <- grep(input$InSample, names(InitDf()$Files))
    # Update smooth
    Df$data2$Ploidy <- input$CtrlPlo1
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf3')
    DT::replaceData(Proxy, Df$data2)
  })
  
  
  # Download CSV
  output$DownloadData <- downloadHandler(
    filename <- function() {
      "peaks.csv"
    },
    content <- function(file) {
      write.csv(Df$data, file, row.names = FALSE)
    }
  )
}

# Create Shiny object
shinyApp(ui = ui, server = server)