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

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  "Flow Peaks",
                  tabPanel("Peaks",
                           sidebarPanel(
                             # Inputs
                             tags$h3("Input:"),
                             # add anotations
                             # Input files
                             fileInput(inputId = "Files",
                                       label = "Upload FCS files",
                                       multiple = TRUE,
                                       accept = c(".FCS", ".fcs")),
                             # Select input file
                             uiOutput('FileDropdown'),
                             # Select channel
                             uiOutput('ChanDropdown'),
                             # Adjust smoothing
                             numericInput(inputId = "InSmooth",
                                          label = "Adjust smoothing",
                                          value = 0.1,
                                          min = 0.02,
                                          max = 1,
                                          step = 0.02),
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
                  tabPanel("Correlation", "Estimate genome size based on known controls (in construction)."),
                  tabPanel("Help", "In construction...")
                  
                )
)

# -------------------------------------- #

# Define server function  
server <- function(input, output) {
  # Inputs operations calculated in server
  # Sample
  output$FileDropdown <- renderUI({
    selectInput(inputId = "Sample",
                label = "Choose a sample",
                choices = ifelse(is.na(input$Files), "", names(Flow())))
  })
  # Channels
  output$ChanDropdown <- renderUI({
    selectInput(inputId = "InChan",
                label = "Choose a channel",
                choices = ifelse(is.na(input$Files), "", names(Flow()[[1]])),
                selected = ifelse(is.na(input$Files), "", names(Flow()[[1]])[20]))
  })
  # Box peaks
  output$InBox <- renderUI({
    checkboxGroupInput(inputId = "InBox", label = "Select G1 and G2 peaks", choices = PlotPoints()$MaxIndex, selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
  })
  
  # Update UI when when files are uploaded
  observeEvent(input$Files, {
    req(Flow())
    # Displays file names
    output$FileDropdown <- renderUI({
      selectInput(inputId = "Sample",
                  label = "Choose a sample",
                  choices = names(Flow()))
    })
    # Display channel names
    output$ChanDropdown <- renderUI({
      selectInput(inputId = "InChan",
                  label = "Choose a channel",
                  choices = names(Flow()[[1]]),
                  selected = names(Flow()[[1]])[20])
    })
  })
  
  # Define specific functions
  # Lines
  GetLine <- function(){
    # Get sample name
    SampNum <- grep(input$Sample, names(Flow()))
    # Get channels information
    ChanNum <- grep(input$InChan, names(Flow()[[SampNum]]))
    # Calculate histogram
    Hist <- hist(exprs(Flow()[[SampNum]][, ChanNum]), breaks = 1:1000, plot = FALSE)
    # Counts is the y histogram variable
    Counts <- Hist$counts
    # Index is the x histogram variable
    Index <- 1:length(Counts)
    # Fit a line to the data
    Fit <- loess(Counts ~ Index, span = Df$data[SampNum, 2])
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create plotting data frame
    data.frame(Fluorescence = Index, Freq = FitLine)
  }
  
  # Points
  GetPoints <- function(){
    # Get sample name
    SampNum <- grep(input$Sample, names(Flow()))
    # Create a second curve with flattened peaks by calculating the maximum in a given window
    # PENDING: The starting width to calculate the window, which can be changed according to needs
    Width <- Df$data[SampNum, 3]
    # The window to calculate the maximum
    Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    FlatLine <- rollapply(data = zoo(PlotLine()$Freq), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    Delta <- FlatLine - PlotLine()$Freq[-c(1:Width, length(PlotLine()$Fluorescence) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    MaxIndex <- which(Delta <= 0) + Width
    # Obtain intensity of the points (y axis)
    Intensity <- PlotLine()$Freq[MaxIndex]
    # Remove peaks with intensity lower than 10
    MaxIndex <- MaxIndex[Intensity > 10]
    Intensity <- Intensity[Intensity > 10]
    # Create points data frame
    data.frame(MaxIndex = MaxIndex, Intensity = Intensity)
  }
  
  
  # Create reactive values
  # Get initial files and information
  # Read FCS files and calculate inital lineas and peasks??
  Flow <- eventReactive(input$Files, {
    FilesLs = list()
    
    Names <- NULL
    
    for(i in 1:length(input$Files[, 1])){
      FilesLs[[i]] <- read.FCS(input$Files[[i, 'datapath']], emptyValue = FALSE, alter.names = TRUE)
      Name <- sub(" .*", " ", identifier(FilesLs[[i]]))
      Names <- c(Names, Name)
    }
    
    names(FilesLs) <- Names
    
    FilesLs
  })
  
#  G1G2Peaks <- eventReactive(input$Files, {
#    PeaksLs = list()
    
#    Names <- NULL
    
#    for(i in 1:length(input$Files[, 1])){
#      FilesLs[[i]] <- read.FCS(input$Files[[i, 'datapath']], emptyValue = FALSE, alter.names = TRUE)
#      Name <- sub(" .*", " ", identifier(FilesLs[[i]]))
#      Names <- c(Names, Name)
#    }
    
#    names(FilesLs) <- Names
    
#    FilesLs
#  })
  
  
  # Create parameter and results data frame to modify initially calculated data
  Df <- reactiveValues(data = NULL)
  # In this case, it had to be created and reactive value in order to be modified via proxy
  observeEvent(input$Files, {
    Df$data <- data.frame(Sample = names(Flow()),
                          Smoothing = rep(0.1, length(Flow())),
                          Window = rep(10, length(Flow())),
                          G1 = NA,#PlotPoints()$MaxIndex[1],
                          G2 = NA)#PlotPoints()$MaxIndex[2])
    
    output$ResDf <- DT::renderDataTable(isolate(Df$data),
                                        editable = FALSE)
  })
  
  
  
  
  
  #NotVisDf <- eventReactive(input$Files, {
  #  CurvesLs <- list()
  #  for (i in i:length(Flow())){
      
  #  }
  #})
  
  # Create plot line
  PlotLine <- eventReactive(c(input$Sample, input$InSmooth, input$InWindow, input$InChan), {
    req(Flow())
    GetLine()
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$Sample, input$InSmooth, input$InWindow, input$InChan, input$MaxPeaks), {
    req(Flow())
    GetPoints()
  })
  
  # Plot data
  
  output$HisPlot <- renderPlot({
    req(input$Files)
    req(Flow())
    req(input$Sample)
    req(PlotLine())
    req(input$MaxPeaks)
    req(input$InBox)
    
    SampNum <- grep(input$Sample, names(Flow()))
    # Get name for plotting and table
    Name <- sub(" .*", "", identifier(Flow()[[SampNum]]))
    # Plot hitogram with points
    ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
      geom_line() +
      annotate("point", x = PlotPoints()$MaxIndex[1:input$MaxPeaks],
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
    SampNum <- grep(input$Sample, names(Flow()))
    # Update smooth
    Df$data[SampNum, 2] <- input$InSmooth
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Window
  observeEvent(input$InWindow, {
    #input$InBox <- ifelse(length(input$InBox) == 0, NA, input$InBox) # Debug minimun number of choices
    SampNum <- grep(input$Sample, names(Flow()))
    # Update windown
    Df$data[SampNum, 3] <- input$InWindow
    # Replace data
    Proxy <- DT::dataTableProxy('ResDf')
    DT::replaceData(Proxy, Df$data)
  })
  # Peaks
  observeEvent(input$InBox, {
    #input$InBox <- ifelse(length(input$InBox) == 0, NA, input$InBox) # Debug minimun number of choices
    SampNum <- grep(input$Sample, names(Flow()))
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