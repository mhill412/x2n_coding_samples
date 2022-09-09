library(shiny)
library(leaflet)
library(leaflet.minicharts)
library(ggplot2)
library(dplyr)
library(plotly)
library(rsconnect)

a4data <- read.csv(
  "https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv"
)
a4data <- a4data[a4data$country %in% c("China","United States","India","Russia",
                                       "Japan", "Germany", "Iran", "South Korea",
                                       "Saudi Arabia", "Indonesia"),]
a4data[is.na(a4data)] <- 0

countries <- c('China', 'Germany', 'India', 'Indonesia', 'Iran', 'Japan',
               'Russia', 'Saudi Arabia', 'South Korea', 'United States')
color <- c("red", "coral2","orange","yellow","brown","green","cyan",
           "blue4", "blue", "purple")


# Define UI for app ----
ui <- fluidPage(
  
  
  
  pageWithSidebar(
    # App title ----
    headerPanel("The World's Biggest CO2 Contributors"),
    
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the year ----
      sliderInput(inputId = "years",
                  label = "Year (1940-2019)",
                  value = "2010",
                  min = 1940,
                  max = 2019,
                  sep = ""
      ),
      
      selectInput(inputId = "emtype",
                  label = "Type of CO2 Emission",
                  choices = c("All CO2 Emissions" = "co2", "From Cement" = "cement_co2", 
                              "From Coal" = "coal_co2", "From Flaring" = "flaring_co2", 
                              "From Gas" = "gas_co2",
                              "From Oil" = "oil_co2", "Other" = "other_industry_co2"),
                  selected = "All CO2 Emissions")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Introduction", verbatimTextOutput("intro")),           
                  tabPanel("Reactive Data", plotlyOutput("barplot"), width = 8)
      )
    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  reactive_data <- reactive({
    selected_year <- a4data[a4data$year %in% c(input$years),]
    bar_data <- selected_year %>% pull(input$emtype)
    return(bar_data)
  })
  
  react_max <- reactive({
    selected_year <- a4data[a4data$year %in% c(input$years),]
    bar_data <- selected_year %>% pull(input$emtype)
    bar_max <- max(bar_data)
    bar_max <- as.vector(bar_max)
    return(bar_max)
  })
  
  
  
  
  output$barplot <- renderPlotly({
    barplot <- plot_ly(x = ~countries, y = ~reactive_data(), type = "bar",
                       text = reactive_data(), 
                       marker = list(color = color)
    )
    barplot <- barplot %>% layout(title = "CO2 Emissions by Country", 
                                  xaxis = list(title = "Country"),
                                  yaxis = list(title = 
                                                 "CO2 Emitted in a Year (in Millions of Tonnes)")
    )
  })
  
  output$intro <- renderText({
    totalco2 <- a4data$co2
    damtype <- select(a4data, matches("_co2"))
    paste("How much total CO2 has the world emitted into the Earth's atmosphere?",
          sum(totalco2), "million tonnes of CO2", "", 
          "Which form of CO2 emission has done the most damage to the atmosphere?",
          sep="\n"
    )
    
  })
  
  
}
shinyApp(ui = ui, server = server)

