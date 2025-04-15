# Load required libraries
library(deSolve)
library(tidyverse)



# Age groups
age_groups <- c("0_19", "20_64", "65_plus")
n_groups <- length(age_groups)

# Population per age group
population <- c(30000, 50000, 20000)  # Example values

# Time-dependent vaccination rate function
vaccination_rate <- function(t, group) {
  if (t < 30) return(0)
  else if (t < 60) {
    if (group == 3) return(0.01)  # Prioritize 65+
    else if (group == 2) return(0.005)  # Adults
    else return(0.002)  # Young
  } else {
    if (group == 3) return(0.02)
    else if (group == 2) return(0.01)
    else return(0.005)
  }
}


# Define SEIRD Model
seird_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Initialize derivatives
    dS <- numeric(n_groups)
    dE <- numeric(n_groups)
    dI <- numeric(n_groups)
    dR <- numeric(n_groups)
    dD <- numeric(n_groups)
    dV <- numeric(n_groups)
    
    for (i in 1:n_groups) {
      # Retrieve current vaccination rate
      v <- vaccination_rate(time, i)
      
      # Differential equations for each age group
      dS[i] <- -beta[i] * S[i] * I[i] - v * S[i]
      dE[i] <- beta[i] * S[i] * I[i] - sigma[i] * E[i]
      dI[i] <- sigma[i] * E[i] - gamma[i] * I[i]
      dR[i] <- (1 - CFR[i]) * gamma[i] * I[i]
      dD[i] <- CFR[i] * gamma[i] * I[i]
      dV[i] <- v * S[i]
    }
    
    # Combine derivatives
    list(c(dS, dE, dI, dR, dD, dV))
  })
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
  S = population / sum(population) - 0.0001,
  E = rep(0.00005, n_groups),
  I = rep(0.00005, n_groups),
  R = rep(0, n_groups),
  D = rep(0, n_groups),
  V = rep(0, n_groups)
)
# Parameters
params <- list(
  beta = c(0.3, 0.4, 0.2),       # Transmission rates
  sigma = rep(1 / 5.2, n_groups),  # Incubation rates
  gamma = rep(1 / 3.0, n_groups),  # Recovery rates
  CFR = c(0.001, 0.01, 0.08)     # Case Fatality Rates
)

# Time vector
times <- seq(0, 100, by = 1)

# Run the model
output <- ode(y = initial_state, times = times, func = seird_model, parms = params)
output_df <- as.data.frame(output)

# Melt the data for plotting
output_melted <- melt(output_df, id.vars = "time")

# Add age group and compartment information
output_melted$Group <- rep(rep(age_groups, each = length(times)), times = 6)
output_melted$Compartment <- rep(c("S", "E", "I", "R", "D", "V"), each = length(times) * n_groups)



# Map numeric group to age label
group_labels <- setNames(age_groups, c("1", "2", "3"))
long_df$Group <- factor(group_labels[long_df$Group], levels = age_groups)

# Plot the results
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = E, color = "Exposed")) +
  geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = R, color = "Recovered")) +
  geom_line(aes(y = D, color = "Dead")) +
  geom_line(aes(y = V, color = "Vaccinated")) +
  labs(
    title = "SEIRD Model with Time-Dependent Vaccination",
    x = "Days",
    y = "Proportion of Population",
    color = "Compartment"
  ) +
  scale_color_manual(values = c(
    "Susceptible" = "blue",
    "Exposed"     = "orange",
    "Infected"    = "red",
    "Recovered"   = "green",
    "Dead"        = "black",
    "Vaccinated"  = "purple"
  )) +
  theme_minimal()