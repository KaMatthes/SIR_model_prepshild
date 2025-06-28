# Load required libraries
library(deSolve)
library(tidyverse)

col_c  <- c( "#1A1F2BFF","#7030A0","#1EA7C4","#FF0066","#EB7746FF","#B9CA5DFF","#24719BFF")
col_sun <- c("#F6A1CAFF", "#1A1F2BFF","#EB7746FF", "#B9CA5DFF","#FBDB82FF", "#9F272EFF", "#D36B8DFF", "#EC9E58FF", "#24719BFF", "#6A1639FF")

lwd_size <- 2
size_plot <- 18
size_strip <- 18
axis_text_size <- 18
axis_title_size <- 18
plot_title_size <- 18
legend_text_size <- 25
# SEI-SIR Dengue Model with Death Compartment
dengue_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Human dynamics with deaths
    dSh <- -beta_vh * Sh * Iv
    dIh <- beta_vh * Sh * Iv - gamma_h * Ih - delta_h * Ih
    dRh <- gamma_h * Ih
    dDh <- delta_h * Ih  # Deaths from infection
    
    # Mosquito dynamics
    dSv <- -beta_hv * Sv * Ih
    dEv <- beta_hv * Sv * Ih - sigma_v * Ev
    dIv <- sigma_v * Ev
    
    list(c(dSh, dIh, dRh, dDh, dSv, dEv, dIv))
  })
}

# Initial state (proportions)
initial_state <- c(
  Sh = 0.99,  # Susceptible humans
  Ih = 0.01,  # Infected humans
  Rh = 0.0,   # Recovered humans
  Dh = 0.0,   # Dead humans
  Sv = 0.95,  # Susceptible mosquitoes
  Ev = 0.04,  # Exposed mosquitoes
  Iv = 0.01   # Infected mosquitoes
)

# Parameters
parameters <- c(
  beta_vh = 0.25,     # Transmission: mosquito → human
  beta_hv = 0.3,      # Transmission: human → mosquito
  gamma_h = 1 / 7,    # Recovery rate (humans)
  delta_h = 0.001,    # Death rate (humans)
  sigma_v = 1 / 10    # Incubation in mosquitoes
)

# Time vector
times <- seq(0, 120, by = 1)

# Run the model
out <- ode(y = initial_state, times = times, func = dengue_model, parms = parameters)
out_df <- as.data.frame(out) %>%
  gather(., comp, prop,Sh:Iv)

# Plot
ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "SEIRD Model Dengue",
       x = "Days",
       y = "Proportion")+
  scale_color_manual("",
                     breaks=c("Sh","Ih","Rh","Dh","Sv","Ev", "Iv"),
                     labels=c("Susceptible Human","Infected Human","Recovered Human","Deceased HUman",
                              "Susceptible Mosquito","Exposed Mosquito","Infected Mosquito"),
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
