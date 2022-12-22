####################################
# Flow Peaks                       #
# An app to calculate peaks from   #
# flow cytometer histograms         #
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
                  tabPanel("Histograms",
                           sidebarPanel(
                             tags$h3("Input:"),
                             # Inputs
                             fileInput(inputId = "Files", label = "Choose FCS File", multiple = TRUE, accept = c(".FCS", ".fcs")),
                             uiOutput('Dropdown1'),
                             uiOutput('Dropdown2'),
                             numericInput(inputId = "Smooth", label = "Smoothing", value = 0.1, min = 0.02, max = 1, step = 0.02),
                             numericInput(inputId = "Window", label = "Window", value = 10, min = 1, max = 100, step = 1),
                             numericInput(inputId = "Peaks", label = "Peaks", value = 5, min = 1, max = 5, step = 1),
                             actionButton(inputId = "Update", label = "Update"),
                             
                           ),
                           mainPanel(
                             h1("Preview"),
                             
                             h4("Output 1"),
                             plotOutput("plot1"),
                             tableOutput("Table"),
                             DT::dataTableOutput("tablepar"),
                             uiOutput("Box"),
                             
                           )
                           
                  ),
                  tabPanel("Navbar 2", "This panel is intentionally left blank"),
                  tabPanel("Navbar 3", "This panel is intentionally left blank")
                  
                )
)

# Define server function  
server <- function(input, output) {
  Flow <- reactive({
    req(input$Files)
    FilesLs = list()
    
    Names <- NULL
    
    for(i in 1:length(input$Files[, 1])){
      FilesLs[[i]] <- read.FCS(input$Files[[i, 'datapath']], emptyValue = FALSE, alter.names = TRUE)
      Name <- sub(" .*", " ", identifier(FilesLs[[i]]))
      Names <- c(Names, Name)
    }
    names(FilesLs) <- Names
    
    return(FilesLs)
  })
  
  Df <- reactiveValues(data = NULL)
  
  observeEvent(input$Files, {
  
  Df$data <- data.frame(Sample = names(Flow()), Smoothing = rep(0.1, length(Flow())), Window = rep(10, length(Flow())))

  output$tablepar <- DT::renderDataTable(isolate(Df$data), editable = FALSE)
  })
  observeEvent(input$Smooth, {
    SampNum <- grep(input$Sample, names(Flow()))
    Df$data[SampNum, 2] <- input$Smooth
    
    proxy1 <- DT::dataTableProxy('tablepar')
    DT::replaceData(proxy1, Df$data)
  })

  output$Dropdown1 <- renderUI({
    req(Flow())
    selectInput("Channels", label = "Channels", choices = names(Flow()[[1]]), selected = names(Flow()[[1]])[20])
  })
  
  output$Dropdown2 <- renderUI({
    req(Flow())
    selectInput("Sample", label = "Sample", choices = names(Flow()))
  })

  output$plot1 <- renderPlot({
    req(Flow())
    SampNum <- grep(input$Sample, names(Flow()))
    # Get name for ploting and table
    Name <- sub(" .*", "", identifier(Flow()[[SampNum]]))
    PeaksDf <- NULL
    # if there is a for should be here?
    # Get channels information
    ChanNum <- grep(input$Channels, names(Flow()[[SampNum]]))
    # Calculate histogram
    Hist <- hist(exprs(Flow()[[SampNum]][, ChanNum]), breaks = 1:1000, plot = FALSE)
    # Counts is the y histogram variable
    Counts <- Hist$counts
    # Index is the x histogram variable
    Index <- 1:length(Counts)
    # Fit a line to the data
    # Fit is the entire result of the calculation. 25% smoothing span is used by default, but it can be tuned
    Fit <- loess(Counts ~ Index, span = input$Smooth)
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create a second curve with flattened peaks by calculating the maximum in a given window
    # The starting width to calculate the window, which can be changed according to needs
    Width <- input$Window
    # The window to calculate the maximum
    Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    FlatLine <- rollapply(data = zoo(FitLine), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    # In this case, we remove the indexes in fit line that do not exist in the flat line after the sliding window
    Delta <- FlatLine - FitLine[-c(1:Width, length(Counts) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    MaxIndex <- which(Delta <= 0) + Width
    # Remove peaks with intensity lower than 10
    Intensity <- FitLine[MaxIndex]
    MaxIndex <- MaxIndex[Intensity > 10]
    Intensity <- Intensity[Intensity > 10]
    # Keep the number of peaks with the highest intensity
    # Maybe keep all and let the user choose the number of displayed peaks and select peaks
    MaxIndex <- MaxIndex[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    Intensity <- Intensity[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    # Create plotting data frame
    PlotDf <- data.frame(Fluorescence = Index, Freq = FitLine)
    # Partial peaks information
    WkDf <- data.frame("Strain" = Name, "G1" = MaxIndex[1], "G2" = MaxIndex[2])
    # Compile peaks information
    PeaksDf <- rbind(PeaksDf, WkDf)
    
    # Order peaks
    PeaksDf <- PeaksDf[order(PeaksDf$G1),]
    
    # The objective would be to create a table that contains the used parameters per sample
    PeaksDf$Smoothing <- input$Smooth
    PeaksDf$Window <- input$Window
    
    output$Table <- renderTable(expr = PeaksDf)#,
                                         #extensions = "Buttons",
                                         #options = list(dom = "Bfrtip",
                                        #                buttons = c("copy", "csv", "excel", "pdf", "print"),
                                        #               searching = FALSE),
                                        # editable = TRUE)
    output$Box <- renderUI({
      checkboxGroupInput(inputId = "Box", label = "Peaks", choices = MaxIndex)
    })
    
    
    ggplot(PlotDf, aes(x = Fluorescence, y = Freq)) +
      geom_line() +
      annotate("point", x = MaxIndex, y = Intensity, col = "red") +
      labs(title = Name) +
      theme(plot.title = element_text(hjust = 0.5))
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)