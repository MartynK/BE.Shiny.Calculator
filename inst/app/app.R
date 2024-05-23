
library(shiny)
library(readxl)
library(dplyr)
library(PowerTOST)

power.table <- function(CV, theta0 = .95, design = "2x2", dropout = 10) {

  res <- pa.ABE(CV=CV,theta0=theta0,design=design,targetpower=.9)$paN

  return(
    res %>%
      select(pwr,N) %>%
      `colnames<-`(c("Power","N completer")) %>%
      mutate(`N with dropout` =
               ceiling((`N completer` / (1 - dropout/100))/2)*2)
  )

}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Calculator for BE stuff"),

    # Sidebar
    sidebarLayout(
        sidebarPanel(
            hr(),
            numericInput(
              "CI_lower",
              "Lower limit of reported Conf.Interval:",
              value = 0.81,
              step = .01
            ),
            numericInput(
              "CI_upper",
              "Upper limit of reported Conf.Interval:",
              value = 1.11,
              step = .01
            ),
            hr(),
            sliderInput("N_sub",
                        "Number of subjects:",
                        min = 3,
                        max = 100,
                        value = 24),
            selectInput("design","Reported design:",
                        c("2x2"="2x2",
                          "2x2x4"="2x2x4",
                          "3x3"="3x3")),


            # # Download button
            # downloadButton("downloadData", "Download results"),

            # Add your discreet message at the bottom
            tags$hr(),  # Horizontal line for separation

            # Footer with discreet message and custom link
            tags$footer(

              # Github link
              tags$p(
                tags$a(href = "https://github.com/MartynK/BE.Shiny.Calculator/", "Github repo"),
                style = "font-size: 80%; color: grey; text-align: center; padding-top: 10px;"
              ),

              # Author's link
              tags$p(
                tags$a(href = "https://www.linkedin.com/in/marton-kiss-md-29308a112/", "by Márton Kiss"),
                style = "font-size: 80%; color: grey; text-align: center; padding-top: 10px;"
              )
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(

          tags$p(
            tags$p(textOutput("orig_val")),
            style = "font-size: 180%; color: grey; text-align: center; padding-top: 10px;"
          ),



          hr(),
          tags$p(
            tags$p("The below table is for a standard 2x2 design"),
            style = "font-size: 180%; color: grey; text-align: center; padding-top: 10px;"
          ),
          tableOutput("end_val")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


  iscv <- reactive({
    PowerTOST::CI2CV(lower = input$CI_lower,
                     upper = input$CI_upper,
                     n = input$N_sub,
                     design = input$design)
  })

  power_tab <- reactive({
      power.table(CV = iscv(),
                  theta0 = .95,
                  design = "2x2") #%>%
      # as.data.frame() %>%
      # `colnames<-`(c("n","dropouts","drp","power")) %>%
      #  mutate( n_completers = n - dropouts,
      #          `n with 10% dropout rate` = ceiling((n_completers / 0.9) / 2)*2,
      #          `Power(%)` = power * 100) %>%
      # dplyr::select(`Power(%)`,n_completers, `n with 10% dropout rate`)


  })


  output$orig_val <- renderText({
    paste0( "ISCV based on the provided data: ", round( iscv() * 100, digits = 2), "%")
    })

  output$end_val  <-  renderTable({ power_tab()},
                                  align = "c",
                                  hover = TRUE,
                                  digits = c(2,0,0))

  # # Server response for download
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     "calculate_daily_data.xlsx"
  #   },
  #   content = function(file) {
  #     file.copy("calculate_daily_data.xlsx", file)
  #   }
  # )
}

# Run the application
shinyApp(ui = ui, server = server)
