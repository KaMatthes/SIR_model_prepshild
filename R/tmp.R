library(deSolve)
library(ggplot2)
library(reshape2)

# SEIHRD model function
seihrd_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Transitions
    dS <- -beta * S * I / N
    dE <-  beta * S * I / N - sigma * E
    
    # From E to I, then to R, H, or D (non-hospital)
    dI <-  sigma * E -
      (1 - p_hosp - p_death_I) * gamma * I -     # Recover
      p_hosp * delta * I -                       # Hospital
      p_death_I * alpha * I                      # Death (non-hospital)
    
    dH <-  p_hosp * delta * I - rho * H - mu * H # Recovery & death from hospital
    dR <-  (1 - p_hosp - p_death_I) * gamma * I + rho * H
    dD <-  p_death_I * alpha * I + mu * H
    
    list(c(dS, dE, dI, dH, dR, dD))
  })
}

# Parameters
parameters <- c(
  beta = 0.4,         # Transmission rate
  sigma = 1/5.2,      # Incubation rate
  gamma = 1/7,        # Recovery (no hospital)
  delta = 1/4,        # Hospitalisation rate
  rho = 1/10,         # Recovery from hospital
  mu = 1/15,          # Death after hospital
  alpha = 1/10,       # Death rate from I
  p_hosp = 0.2,       # Proportion of I to hospital
  p_death_I = 0.05    # Proportion of I dying before hospital
)

# Initial states
init_state <- c(
  S = 9990,
  E = 10,
  I = 0,
  H = 0,
  R = 0,
  D = 0
)

# Total population (constant)
N <- sum(init_state)

# Time sequence
times <- seq(0, 160, by = 1)

# Solve ODE
out <- ode(y = init_state, times = times, func = seihrd_model, parms = parameters)
out <- as.data.frame(out)

# Plotting
out_long <- melt(out, id = "time")

ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1.2) +
  labs(title = "SEIHRD Model: Hospitalisations and Deaths (Before & After Hospital)",
       x = "Time (days)",
       y = "Population",
       color = "Compartment") +
  theme_minimal()
