library(knitr)
library(tidyverse)
library(shiny)
library(htmltools)
library(deSolve)

col_c  <- c("#FF0066","#7030A0","#1EA7C4" )

lwd_size <- 2
size_plot <- 15
size_strip <- 15
axis_text_size <- 15
axis_title_size <- 15
plot_title_size <- 15
legend_text_size <- 15
  
ui <- fluidPage(
  fluidRow(
    column(4, sliderInput(inputId = "repro",label = "Reproductive number", value = 1.5, min = 1, max = 6)),
    column(4,sliderInput(inputId = "period",label = "Infection period", value = 5, min = 1, max = 14))),
  fluidRow(
    column(4,sliderInput(inputId = "suscep",label = "Initial susceptible", value = 10000, min = 10000, max = 20000000)),
    column(4, sliderInput(inputId = "initial",label = "Initial infected", value = 1, min = 1, max = 100))),
   fluidRow(
     column(4,sliderInput(inputId = "recov",label = "Initial recovered", value = 0, min = 0, max = 100)),
     column(4,sliderInput(inputId = "timemax",label = "Max days", value = 100, min = 10, max = 150))),
  fluidRow(1,plotOutput("SIR"))
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
      D = input$period    # Average infection period in days
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
    beta <- params$R_0 * gamma / N # transmission rate
    

    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I

    list(c(dS, dI, dR))
  }

  # Solve ODE and convert to data frame
  dt1 <- reactive({
    as.data.frame(ode(y = initial_state(), times = times(), func = sir_model, parms = params()))
  })

  # Output for plot
  output$SIR <- renderPlot({
    ggplot(dt1(), aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptible"), lwd=lwd_size) +
      geom_line(aes(y = I, color = "Infected") , lwd=lwd_size) +
      geom_line(aes(y = R, color = "Recovered") , lwd=lwd_size) +
      labs(title = "SIR Model for Influenza", x = "Days", y = "Proportion of Population") +
      scale_color_manual("",values = col_c) +
       theme_bw()+
  theme(aspect.ratio = 0.4,
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
