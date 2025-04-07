library(knitr)
library(tidyverse)
library(shiny)
library(htmltools)
library(deSolve)


col_c  <- c("#FF0066","#7030A0","#1EA7C4" )

lwd_size <- 2
size_plot <- 18
size_strip <- 18
axis_text_size <- 18
axis_title_size <- 18
plot_title_size <- 18
legend_text_size <- 25
  
ui <- fluidPage(
  fluidRow(
    column(6, sliderInput(inputId = "repro",label = "Reproductive number", value = 1.5, min = 1, max = 6)),
    column(6,sliderInput(inputId = "period",label = "Infection period", value = 5, min = 2, max = 7))),
  fluidRow(
    column(6,sliderInput(inputId = "suscep",label = "Initial susceptible", value = 100000, min = 1000, max = 20000000)),
    column(6, sliderInput(inputId = "initial",label = "Initial infected", value = 1, min = 1, max = 100))),
   fluidRow(
     column(6,sliderInput(inputId = "recov",label = "Initial recovered", value = 0, min = 0, max = 10)),
     column(6,sliderInput(inputId = "timemax",label = "Max days", value = 100, min = 10, max = 150))),
  dataTableOutput('values'),
  plotOutput("SIR")
)

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)

  # Initial conditions
  S0 <- reactive({ input$suscep })
  I0 <- reactive({ input$initial })
  R0 <- reactive({ input$recov }) 
  N0 <- reactive({ S0() + I0() + R0()}) 
  

  # Parameters
  params <- reactive({
    list(
      R_0 = input$repro,  # reproductive number
      D = input$period      # Average illness duration in days
    )
  })

  # Time sequence
  times <- reactive({
    seq(0, input$timemax, by = 1)
  })

  # Normalize by population size
  initial_state <- reactive({
    c(S = S0() / N0(), I = I0() / N0(), R = R0() / N0())
  })

  # SIR model function
  sir_model <- function(time, state, params) {
    S <- state[1]
    I <- state[2]
    R <- state[3]
    
    N <- S+I+R

    gamma <- 1 / params$D  # Recovery rate
  beta <- params$R_0 * gamma / N
    

    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I

    list(c(dS, dI, dR))
  }

  # Solve ODE and convert to data frame
  dt1 <- reactive({
    as.data.frame(ode(y = initial_state(), times = times(), func = sir_model, parms = params()))
  })

  # Output for data table
  
 # R_rep <- reactive({
 #    R_rep <- params()$beta*(1/params()$D)
 #     return(R_rep)   
 #  })
 #  
  # output$values <- renderDataTable({
  #   data.frame(dt1())
  #   })

  # Output for plot
  output$SIR <- renderPlot({
    ggplot(dt1(), aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptible"), lwd=lwd_size) +
      geom_line(aes(y = I, color = "Infected") , lwd=lwd_size) +
      geom_line(aes(y = R, color = "Recovered") , lwd=lwd_size) +
      labs(title = "SIR Model for Influenza", x = "Days", y = "Proportion of Population") +
      scale_color_manual("",values = col_c) +
       theme_bw()+
  theme(
    strip.text = element_text(size=size_plot),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_title_size),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())
  })
}

shinyApp(ui = ui, server = server)
