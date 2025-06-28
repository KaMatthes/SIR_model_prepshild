
library(deSolve)
library(ggplot2)

# SEI-SIR with seasonal mosquito births
dengue_model_seasonal <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Seasonal birth rate of mosquitoes
    b_v <- b0 * (1 + a * sin(2 * pi * t / T))
    
    # Human dynamics
    dSh <- -beta_vh * Sh * Iv
    dIh <- beta_vh * Sh * Iv - gamma_h * Ih
    dRh <- gamma_h * Ih
    
    # Mosquito dynamics
    dSv <- b_v * (Sv + Ev + Iv) - beta_hv * Sv * Ih - mu_v * Sv
    dEv <- beta_hv * Sv * Ih - sigma_v * Ev - mu_v * Ev
    dIv <- sigma_v * Ev - mu_v * Iv
    
    list(c(dSh, dIh, dRh, dSv, dEv, dIv))
  })
}

# Initial state
initial_state <- c(
  Sh = 0.99,  Ih = 0.01,  Rh = 0.0,
  Sv = 0.95,  Ev = 0.04,  Iv = 0.01
)

# Parameters including seasonality
parameters <- c(
  beta_vh = 0.25,   # Mosquito → human transmission
  beta_hv = 0.3,    # Human → mosquito transmission
  gamma_h = 1 / 7,  # Recovery rate (humans)
  sigma_v = 1 / 10, # Incubation in mosquitoes
  mu_v = 1 / 14,    # Mosquito death rate (lifespan ~14 days)
  b0 = 1 / 14,      # Avg mosquito birth rate
  a = 0.3,          # Amplitude of seasonality (0 = none, 1 = max)
  T = 365           # Period (365 days = 1 year)
)

# Time vector
times <- seq(0, 365, by = 1)  # simulate one year

# Solve model
out <- ode(y = initial_state, times = times, func = dengue_model_seasonal, parms = parameters)
out_df <- as.data.frame(out)

# Plotting
ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = Sh, color = "Susceptible (Humans)")) +
  geom_line(aes(y = Ih, color = "Infected (Humans)")) +
  geom_line(aes(y = Rh, color = "Recovered (Humans)")) +
  geom_line(aes(y = Sv, color = "Susceptible (Mosquitoes)"), linetype = "dashed") +
  geom_line(aes(y = Ev, color = "Exposed (Mosquitoes)"), linetype = "dashed") +
  geom_line(aes(y = Iv, color = "Infected (Mosquitoes)"), linetype = "dashed") +
  labs(
    title = "SEI-SIR Dengue Model with Seasonal Mosquito Dynamics",
    x = "Time (days)",
    y = "Proportion of Population",
    color = "Compartments"
  ) +
  theme_minimal()
