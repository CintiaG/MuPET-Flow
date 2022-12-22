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
                             #selectInput("Channels", label = "Channels", choices = NA),
                             uiOutput('Dropdown1'),
                             uiOutput('Dropdown2'),
                             numericInput(inputId = "Smooth", label = "Smoothing", value = 0.1, min = 0.02, max = 1, step = 0.02),
                             numericInput(inputId = "Window", label = "Window", value = 10, min = 1, max = 100, step = 1),
                             numericInput(inputId = "Peaks", label = "Peaks", value = 5, min = 1, max = 5, step = 1),
                             uiOutput("Box"),
                             actionButton(inputId = "Save", label = "Update (pending)"),
                             
                           ),
                           mainPanel(
                             h1("Preview"),
                             
                             h4("Output 1"),
                             plotOutput("plot1"),
                             #tableOutput("Table"),
                             DT::dataTableOutput("tablepar"),
                           )
                           
                  ),
                  tabPanel("Navbar 2", "This panel is intentionally left blank"),
                  tabPanel("Navbar 3", "This panel is intentionally left blank")
                  
                )
)

# Define server function  
server <- function(input, output) {
  # Create parameter and results data frame
  Df <- reactiveValues(data = NULL)
  # In this case, it had to be created and reactive value in order to be modified via proxy
  observeEvent(input$Files, {
    Df$data <- data.frame(Sample = names(Flow()),
                          Smoothing = rep(0.1, length(Flow())),
                          Window = rep(10, length(Flow())),
                          G1 = rep(NA, length(Flow())),
                          G2 = rep(NA, length(Flow())))
    
    output$tablepar <- DT::renderDataTable(isolate(Df$data), editable = FALSE)
  })
  
  # Create reactive values
  # Read FCS files
  Flow <- eventReactive(input$Files, {
    #req(input$Files)
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
  
  # Create plot line
  PlotLine <- eventReactive(c(input$Sample, input$Smooth, input$Window), {
    req(Flow())
    #req(input$Sample)
    SampNum <- grep(input$Sample, names(Flow()))
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
    Fit <- loess(Counts ~ Index, span = Df$data[SampNum, 2])
    # FitLine is the smoothed line
    FitLine <- predict(Fit)
    # Create plotting data frame
    #PlotDf <- data.frame(Fluorescence = Index, Freq = FitLine)
    data.frame(Fluorescence = Index, Freq = FitLine)
  })
  
  # Create plot points
  PlotPoints <- eventReactive(c(input$Sample, input$Smooth, input$Window, input$Peaks), {
    req(Flow())
    SampNum <- grep(input$Sample, names(Flow()))
    # Get name for plotting and table
    #Name <- sub(" .*", "", identifier(Flow()[[SampNum]]))
    PeaksDf <- NULL
    # Create a second curve with flattened peaks by calculating the maximum in a given window
    # The starting width to calculate the window, which can be changed according to needs
    Width <- Df$data[SampNum, 3]
    # The window to calculate the maximum
    Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    FlatLine <- rollapply(data = zoo(PlotLine()$Freq), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    # In this case, we remove the indexes in fit line that do not exist in the flat line after the sliding window
    Delta <- FlatLine - PlotLine()$Freq[-c(1:Width, length(PlotLine()$Fluorescence) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    MaxIndex <- which(Delta <= 0) + Width
    # Remove peaks with intensity lower than 10
    Intensity <- PlotLine()$Freq[MaxIndex]
    # May not be necesary if I am portraying only a few
    MaxIndex <- MaxIndex[Intensity > 10]
    Intensity <- Intensity[Intensity > 10]
    # Keep the number of peaks with the highest intensity
    # Maybe keep all and let the user choose the number of displayed peaks and select peaks
    #MaxIndex <- MaxIndex[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    #Intensity <- Intensity[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    MaxIndex <- MaxIndex[sort(order(Intensity, decreasing = TRUE))]
    Intensity <- Intensity[sort(order(Intensity, decreasing = TRUE))]
    #MaxIndex <- MaxIndex[sort(order(Intensity, decreasing = TRUE))]
    #Intensity <- Intensity[sort(order(Intensity, decreasing = TRUE))]
    
    data.frame(MaxIndex = MaxIndex, Intensity = Intensity)
  })
  
  #Peaks <- eventReactive()
  
  
  # Proxys for replacemnet
  
  # Update smooth
  observeEvent(input$Smooth, {
    SampNum <- grep(input$Sample, names(Flow()))
    Df$data[SampNum, 2] <- input$Smooth
    #Df$data[SampNum, 3] <- input$Window
    
    proxy1 <- DT::dataTableProxy('tablepar')
    DT::replaceData(proxy1, Df$data)
  })
 
  # Update window
   observeEvent(input$Window, {
    SampNum <- grep(input$Sample, names(Flow()))
    Df$data[SampNum, 3] <- input$Window
    
    proxy1 <- DT::dataTableProxy('tablepar')
    DT::replaceData(proxy1, Df$data)
  })
   
   observeEvent(input$Box, {
     SampNum <- grep(input$Sample, names(Flow()))
     Df$data[SampNum, 4] <- input$Box[1]
     Df$data[SampNum, 5] <- input$Box[2]
     
     proxy1 <- DT::dataTableProxy('tablepar')
     DT::replaceData(proxy1, Df$data)
   })

  output$Dropdown1 <- renderUI({
    #req(Flow())
    selectInput("Channels", label = "Channels", choices = ifelse(is.na(input$Files), "", names(Flow()[[1]])), selected = ifelse(is.na(input$Files), "", names(Flow()[[1]])[20]))
    #selectInput("Channels", label = "Channels", choices = names(Flow()[[1]]), selected = names(Flow()[[1]])[20])
  })
  
  # Update select input when files are uploaded
  observeEvent(input$Files, {
    req(Flow())
    output$Dropdown1 <- renderUI({
      selectInput("Channels", label = "Channels", choices = names(Flow()[[1]]), selected = names(Flow()[[1]])[20])
    })
  })
  
  output$Dropdown2 <- renderUI({
    #req(Flow())
    selectInput("Sample", label = "Sample", choices = ifelse(is.na(input$Files), "", names(Flow())))
  })
  
  observeEvent(input$Files, {
    req(Flow())
    output$Dropdown2 <- renderUI({
      selectInput("Sample", label = "Sample", choices = names(Flow()))
    })
  })
  
  output$Box <- renderUI({
    checkboxGroupInput(inputId = "Box", label = "Select peaks", choices = PlotPoints()$MaxIndex, selected = c(PlotPoints()$MaxIndex[1], PlotPoints()$MaxIndex[2]))
  })

  output$plot1 <- renderPlot({
    req(Flow())
    req(input$Sample)
    req(PlotLine())
    req(input$Peaks)
    req(input$Box)
    
    SampNum <- grep(input$Sample, names(Flow()))
    # Get name for plotting and table
    Name <- sub(" .*", "", identifier(Flow()[[SampNum]]))
    #PeaksDf <- NULL
    # Create a second curve with flattened peaks by calculating the maximum in a given window
    # The starting width to calculate the window, which can be changed according to needs
    ###Width <- Df$data[SampNum, 3]
    # The window to calculate the maximum
    ###Window <- 2 * Width + 1
    # Do a sliding window of 2 * Width + 1, and step of 1, and calculate the maximum value of a given window
    ###FlatLine <- rollapply(data = zoo(PlotLine()$Freq), width = Window, FUN = max, align = "center")
    # Calculate the difference between flattened and fitted lines
    # In this case, we remove the indexes in fit line that do not exist in the flat line after the sliding window
    ###Delta <- FlatLine - PlotLine()$Freq[-c(1:Width, length(PlotLine()$Fluorescence) + 1 - 1:Width)]
    # Obtain the indices in which the difference equals to zero, and adjust to the width because of the sliding window step
    ###MaxIndex <- which(Delta <= 0) + Width
    # Remove peaks with intensity lower than 10
    ###Intensity <- PlotLine()$Freq[MaxIndex]
    # May not be necesary if I am portraying only a few
    ###MaxIndex <- MaxIndex[Intensity > 10]
    ###Intensity <- Intensity[Intensity > 10]
    # Keep the number of peaks with the highest intensity
    # Maybe keep all and let the user choose the number of displayed peaks and select peaks
    ###MaxIndex <- MaxIndex[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    ###Intensity <- Intensity[sort(order(Intensity, decreasing = TRUE)[1:input$Peaks])]
    
    ###Res <- data.frame(MaxIndex = MaxIndex, Intensity = Intensity)
    # Partial peaks information
    #WkDf <- data.frame("Strain" = Name, "G1" = PlotPoints()$MaxIndex[1], "G2" = PlotPoints()$MaxIndex[2])
    # Compile peaks information
    #PeaksDf <- data.frame("Strain" = Name, "G1" = input$Box[1], "G2" = input$Box[2])
    
    # Order peaks
    #PeaksDf <- PeaksDf[order(PeaksDf$G1),]
    
    # The objective would be to create a table that contains the used parameters per sample
    # This possible has to be substituted
    #PeaksDf$Smoothing <- input$Smooth
    #PeaksDf$Window <- input$Window
    
    
    
    
    #output$Table <- renderTable(expr = PeaksDf)#,
                                         #extensions = "Buttons",
                                         #options = list(dom = "Bfrtip",
                                        #                buttons = c("copy", "csv", "excel", "pdf", "print"),
                                        #               searching = FALSE),
                                        # editable = TRUE)
    # Para sacar esto, necesito ambie sacar Res
    
    
    
    ggplot(PlotLine(), aes(x = Fluorescence, y = Freq)) +
      geom_line() +
      annotate("point", x = PlotPoints()$MaxIndex[1:input$Peaks], y = PlotPoints()$Intensity[1:input$Peaks], col = "red") +
      labs(title = Name) +
      theme(plot.title = element_text(hjust = 0.5))
  })
}

# Create Shiny object
shinyApp(ui = ui, server = server)