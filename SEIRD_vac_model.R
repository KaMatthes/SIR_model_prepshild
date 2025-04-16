# Load required libraries
library(deSolve)
library(tidyverse)
col_c  <- c( "#1A1F2BFF","#7030A0","#1EA7C4","#FF0066","#EB7746FF","#B9CA5DFF")
col_sun <- c("#F6A1CAFF", "#1A1F2BFF","#EB7746FF", "#B9CA5DFF","#FBDB82FF", "#9F272EFF", "#D36B8DFF", "#EC9E58FF", "#24719BFF", "#6A1639FF")

lwd_size <- 2
size_plot <- 18
size_strip <- 18
axis_text_size <- 18
axis_title_size <- 18
plot_title_size <- 18
legend_text_size <- 25

# Time-dependent vaccination rate function
vaccination_rate <- function(t) {
  if (t < 30) {
    return(0)       # No vaccination in the first 30 days
  } else if (t %in% 30:60) {
    return(0.005)   # Moderate vaccination rate between days 30 and 59
  } else {
    return(0.01)    # Higher vaccination rate from day 60 onwards
  }
}
  
# Define SEIRD Model
seird_model <- function(time, state, parameters) {
  
  # 
  v <- vaccination_rate(time)
  
  S <- state[1]
  E <- state[2]
  I <- state[3]
  R <- state[4]
  D <- state[5]
  
  beta <- parameters["beta"]
  D_inf <- parameters["D"]
  sigma <- parameters["sigma"]
  CFR <- parameters["CFR"]
  
  gamma <- 1 / D_inf # Recovery rate
  
  dS <- -beta * S * I - v * S 
  dE <- beta * S * I - sigma * E
  dI <- sigma * E - gamma * I
  dR <- (1 - CFR) * gamma * I + v * S
  dD <- CFR * gamma * I
  dV <- v * S  # Vaccinated individuals
  
  list(c(dS, dE, dI, dR, dD, dV))
}

# Initial conditions
N <- 100000          # Total population
I0 <- 1              # Initial infected
E0 <- 10             # Initial exposed
S0 <- N - I0 - E0    # Initial susceptible
R0 <- 0              # Initial recovered
D0 <- 0              # Initial dead
V0 <- 0              # Initial vaccinated individuals

initial_state <- c(
  S = S0 / N,
  E = E0 / N,
  I = I0 / N,
  R = R0 / N,
  D = D0 / N,
  V = V0 / N
)

# Parameters
params <- c(
  beta = 0.5,      # Transmission rate
  D = 5,           # Duration of infection
  sigma = 1 / 2,   # Incubation rate (1 / incubation period)
  CFR = 0.05   # Case Fatality Rate (1%)

)

# Time vector
times <- seq(0, 300, by = 1)

# Run the model
out <- ode(y = initial_state, times = times, func = seird_model, parms = params)
out_df <- as.data.frame(output) %>%
  gather(., comp, prop,S:V)

# Plot the results
ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "SEIRD Model with Time-Dependent Vaccination",
    x = "Days",
    y = "Proportion of Population")+
  scale_color_manual("",
                     breaks=c("D","I","S","E","R","V"),
                     labels=c("Deceased","Infected","Susceptible", "Exposed","Recovered","Vaccinated"),
                     values = col_c) +
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

