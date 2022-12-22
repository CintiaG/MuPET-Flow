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

# List of bugs
# That ghace the display smooth when sample is selected
# Zero thicks
# that g1 and g2 labels are the ones ticked
# Add channel to output
# maybe do not change line and cruces but save plots?

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  "Flow Peaks",
                  tabPanel("Peaks",
                           sidebarPanel(
                             # Inputs
                             tags$h3("Input:"),
                             # Input files
                             fileInput(inputId = "InFiles",
                                       label = "Upload FCS files",
                                       multiple = TRUE,
                                       accept = c(".FCS", ".fcs")),
                             # Select input file
                             uiOutput("FileDropdown"),
                             # Select channel
                             uiOutput("ChanDropdown"),
                             # Adjust smoothing
                             uiOutput("SmoothSel"),
                             #numericInput(inputId = "InSmooth",
                            #              label = "Adjust smoothing",
                            #              value = 0.1,
                            #              min = 0.02,
                            #              max = 1,
                            #              step = 0.02),
                             # Adjust window
                             numericInput(inputId = "InWindow",
                                          label = "Adjust window",
                                          value = 10,
                                          min = 1,
                                          max = 100,
                                          step = 1),
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
                             actionButton(inputId = "InSave",
                                          label = "Save (pending)"),
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
                  tabPanel("Regression", "Estimate genome size based on known controls (in construction)."),
                  tabPanel("Help", "In construction...")
                )
)

# ------------------------------------------------------------------------------------------------------------------- #

# Define server function  
server <- function(input, output) {
  # Inputs operations calculated in server
  # Sample
  output$FileDropdown <- renderUI({
    selectInput(inputId = "Sample",
                label = "Select a sample",
                choices = ifelse(is.na(input$InFiles), "", names(Flow()$Files)))
  })
  # Channels
  output$ChanDropdown <- renderUI({
    selectInput(inputId = "InChan",
                label = "Select a channel",
                choices = ifelse(is.na(input$InFiles), "", names(Flow()$Files[[1]])),
                selected = ifelse(is.na(input$InFiles), "", names(Flow()$Files[[1]])[20]))
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
  # Box peaks
  output$InBox <- renderUI({
    checkboxGroupInput(inputId = "InBox",
                       label = "Select G1 and G2 peaks",
                       choices = PlotPoints()$MaxIndex,
                       selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
  })
  
  # Update UI when when files are uploaded
  observeEvent(input$InFiles, {
    req(Flow())
    # Displays file names
    output$FileDropdown <- renderUI({
      selectInput(inputId = "Sample",
                  label = "Choose a sample",
                  choices = names(Flow()$Files))
    })
    # Display channel names
    output$ChanDropdown <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Choose a channel",
                  choices = names(Flow()$Files[[1]]),
                  selected = names(Flow()$Files[[1]])[20])
    })
  })
  
#  observeEvent(input$Sample, {
#    req(Flow())
#    numericInput(inputId = "InSmooth",
#                 label = "Adjust smoothing",
#                 value = 0.1,
#                 min = 0.02,
#                 max = 1,
#                 step = 0.02),
#  })
  
  # Define specific functions
  # Lines
  GetLine <- function(File, ChanNum, Span){
    # Get sample name
    #SampNum <- grep(input$Sample, names(Flow()$Files))
    # Get channels information
    #ChanNum <- grep(input$InChan, names(Flow()$Files[[SampNum]]))
    # Calculate histogram
    #Hist <- hist(exprs(Flow()$Files[[SampNum]][, ChanNum]), breaks = 1:1000, plot = FALSE)
    Hist <- hist(exprs(File[, ChanNum]), breaks = 1:1000, plot = FALSE)
    # Counts is the y histogram variable
    Counts <- Hist$counts
    # Index is the x histogram variable
    Index <- 1:length(Counts)
    # Fit a line to the data
    #Fit <- loess(Counts ~ Index, span = Df$data[SampNum, 2])
    Fit <- loess(Counts ~ Index, span = Span)
    #Fit <- loess(Counts ~ Index, span = ifelse(is.na(Df$data), 0.1, Df$data[SampNum, 2]))
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create plotting data frame
    data.frame(Fluorescence = Index, Freq = FitLine)
  }
  
  # Points
  GetPoints <- function(PlotLine, Width){
    # Get sample name
    #SampNum <- grep(input$Sample, names(Flow()$Files))
    # Create a second curve with flattened peaks by calculating the maximum in a given window
    # PENDING: The starting width to calculate the window, which can be changed according to needs
    #Width <- Df$data[SampNum, 3]
    #Width <- Df$data[SampNum, 3]
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
  # Read FCS files and calculate inital lineas and peasks??
  Flow <- eventReactive(input$InFiles, {
    FilesLs <- list()
    Spans <- NULL
    Windows <- NULL
    LineLs <- list()
    PointLs <- list()
    G1s <- NULL
    G2s <- NULL
    
    Names <- NULL
    
    for(i in 1:length(input$InFiles[, 1])){
      FilesLs[[i]] <- read.FCS(input$InFiles[[i, 'datapath']], emptyValue = FALSE, alter.names = TRUE)
      Name <- sub(" .*", " ", identifier(FilesLs[[i]]))
      Names <- c(Names, Name)
      Spans[i] <- 0.1
      Windows[i] <- 10
      #LineLs[[i]] <- GetLine(File = FilesLs[[i]], SampNum = 1, ChanNum = 20, Span = 0.1)
      LineLs[[i]] <- GetLine(File = FilesLs[[i]], ChanNum = 20, Span = 0.1)
      PointLs[[i]] <- GetPoints(PlotLine = LineLs[[i]], Width = 10)
      G1s[i] <- PointLs[[i]]$MaxIndex[1]
      G2s[i] <- PointLs[[i]]$MaxIndex[2]
    }
    
    names(FilesLs) <- Names
    #names(SpanLs) <- Names
    # I culd retrieve the line aswell
    list(Files = FilesLs, Spans = Spans, Windows = Windows, G1s = G1s, G2s = G2s)
  })

  # Create parameter and results data frame to modify initially calculated data
  Df <- reactiveValues(data = NULL)
  # In this case, it had to be created and reactive value in order to be modified via proxy
  observeEvent(input$InFiles, {
    req(Flow())
    Df$data <- data.frame(Sample = names(Flow()$Files),
                          Smoothing = Flow()$Spans,#rep(0.1, length(Flow()$Files)), # Jalarlo de FLOW?
                          Window = Flow()$Windows,#rep(10, length(Flow()$Files)),
                          G1 = Flow()$G1s,#PlotPoints()$MaxIndex[1],
                          G2 = Flow()$G2s)#PlotPoints()$MaxIndex[2]) # For some reason doe not recongnize the value of the thir sample...
    
    output$ResDf <- DT::renderDataTable(isolate(Df$data),
                                        editable = FALSE)
  })

  # Create plot line
  PlotLine <- eventReactive(c(input$Sample, input$InSmooth, input$InWindow, input$InChan), {
    req(Flow())
    # Get sample name
    SampNum <- grep(input$Sample, names(Flow()$Files))
    # Get file information
    File <- Flow()$Files[[SampNum]]
    # Get channels information
    ChanNum <- grep(input$InChan, names(File))
    #GetLine(File = File, SampNum = SampNum, ChanNum = ChanNum, Span = Df$data[SampNum, 2])
    GetLine(File = File, ChanNum = ChanNum, Span = Df$data[SampNum, 2])
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$Sample, input$InSmooth, input$InWindow, input$InChan, input$MaxPeaks), {
    req(Flow())
    # Get sample name
    SampNum <- grep(input$Sample, names(Flow()$Files))
    Width <- Df$data[SampNum, 3]
    GetPoints(Width = Width, PlotLine = PlotLine())
  })
  
  # Plot data
  
  output$HisPlot <- renderPlot({
    req(input$InFiles)
    req(Flow())
    req(input$Sample)
    req(PlotLine())
    req(input$MaxPeaks)
    req(input$InBox)
    
    SampNum <- grep(input$Sample, names(Flow()$Files))
    # Get name for plotting and table
    Name <- sub(" .*", "", identifier(Flow()$Files[[SampNum]]))
    # Plot hitogram with points
    ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
      geom_line() +
      annotate("point", x = PlotPoints()$MaxIndex[1:input$MaxPeaks], # bug maybe here an and bug of zero ticks
               y = PlotPoints()$Intensity[1:input$MaxPeaks],
               col = "red") +
      # but here I wat to use the box values and restrict it to two
      annotate("text", x = PlotPoints()$MaxIndex[1:2],
               y = PlotPoints()$Intensity[1:2] + 10,
               col = "black",
               label = c("G1", "G2")) +
      labs(title = Name) +
      theme(plot.title = element_text(hjust = 0.5))
  })

  # Proxys for data replacement
  # Smoothing
  observeEvent(input$InSmooth, {
    SampNum <- grep(input$Sample, names(Flow()$Files))
    # Update smooth
    Df$data[SampNum, 2] <- input$InSmooth
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Window
  observeEvent(input$InWindow, {
    #input$InBox <- ifelse(length(input$InBox) == 0, NA, input$InBox) # Debug minimun number of choices
    SampNum <- grep(input$Sample, names(Flow()$Files))
    # Update windown
    Df$data[SampNum, 3] <- input$InWindow
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Peaks
  observeEvent(input$InBox, {
    #input$InBox <- ifelse(length(input$InBox) == 0, NA, input$InBox) # Debug minimun number of choices
    SampNum <- grep(input$Sample, names(Flow()$Files))
    # Update G1 and G2 peaks
    Df$data[SampNum, 4] <- input$InBox[1]
    Df$data[SampNum, 5] <- input$InBox[2]
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)