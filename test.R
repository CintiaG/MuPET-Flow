library(shinydashboard)
library(shiny)

QRSList <- c("Box1","Box2","Box3","Box4","Box5")

ui <- dashboardPage(
  dashboardHeader(title = "render Boxes"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Test", tabName = "Test")
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "Test",
              fluidRow(
                tabPanel("Boxes",uiOutput("myboxes"))
              )  
      )
    )   
  )
)


server <- function(input, output) {
  
  v <- list()
  for (i in 1:length(QRSList)){
    v[[i]] <- box(width = 3, background = "blue",
                  title = h3(QRSList[i], style = "display:inline; font-weight:bold"),
                  selectInput(paste0("slider",i), label = NULL,choices = list("Not good" = "danger", "average" = "warning", "good" = "success"))
    )
  }
  output$myboxes <- renderUI(v)
}

shinyApp(ui = ui, server = server)