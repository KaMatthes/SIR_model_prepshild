# Load required libraries
library(deSolve)
library(tidyverse)

options(scipen=999)
# "#FF0066"
col_c  <- c("#7030A0","#1EA7C4","#F6A1CAFF","#EB7746FF","#B9CA5DFF", "#1A1F2BFF")
col_sun <- c("#F6A1CAFF", "#1A1F2BFF","#EB7746FF", "#B9CA5DFF","#FBDB82FF", "#9F272EFF", "#D36B8DFF", "#EC9E58FF", "#24719BFF", "#6A1639FF")

lwd_size <- 3
size_plot <- 30
size_strip <- 30
axis_text_size <- 30
axis_title_size <- 30
plot_title_size <- 30
legend_text_size <- 30

# Define SEIRD Model
# Parameters mild
# R0 <- 1.5             # Basic reproduction number
# CFR <-0.0017      # Case fatality rate
# p_hosp <- 0.0023      # Proportion of infections hospitalized
# p_die_hosp <- 0.15     # Proportion of hospitalized who die

# Parameters severe
R0 <- 3             # Basic reproduction number
CFR <- 0.007      # Case fatality rate
p_hosp <- 0.03 # Proportion of infections hospitalized
p_die_hosp <- 0.25    # Proportion of hospitalized who die

# for both
incub_period <- 2
infection_period <- 5
Hr <-  10 # times in hospital

# calculate from given assumptions
sigma <- 1 / incub_period        # Incubation rate (2 days)
gamma <- 1 / infection_period        # Recovery rate (5 days)
mu <- p_die_hosp / Hr      # Death rate in hospital
rho <- (1 - p_die_hosp) / Hr  # Recovery rate in hospital
beta <- R0 * gamma   # Transmission rate


# Population and transmission rate
N <- 100000



# Initial state
initial_state <- c(
  S = N - 1, E = 0, I = 1, H = 0, R = 0, D = 0
)

# SEIRHD model function
seirhd_model <- function(time, state, parameters) {
  with(as.list(state), {
    
    dS <- -beta * S * I / N
    dE <-  beta * S * I / N - sigma * E
    dI <-  sigma * E - gamma * I
    dH <-  p_hosp * gamma * I - (rho + mu) * H
    dR <-  (1 - p_hosp - CFR) * gamma * I + rho * H
    dD <-  CFR * gamma * I + mu * H
    
    list(c(dS, dE, dI, dH, dR, dD))
  })
}


# Time vector
times <- seq(0, 300, by = 1)

# Run the model
out <- ode(y = initial_state, times = times, func = seirhd_model , parms = params)
out_df <- as.data.frame(out) %>%
  gather(., comp, prop,S:D) 

# Plot the results


ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "Influenza - Mild Scenario",
       x = "Days",
       y = "Individuals") +
  scale_color_manual("",
                     breaks=c("S","E","I","H","R","D"),
                     labels=c("Susceptible","Exposed","Infected","Hospitalisation","Recovered","Deceased"),
                     values = col_c) +
  scale_y_continuous(breaks = seq(0, 100000, by = 10000))+
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

ggsave("figures/influenza/mild.png",h=12,w=15)


ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "Influenza - Severe Scenario",
    x = "Days",
    y = "Individuals")+
  scale_color_manual("",
                     breaks=c("S","E","I","H","R","D"),
                     labels=c("Susceptible","Exposed","Infected","Hospitalisation","Recovered","Deceased"),
                     values = col_c) +
  scale_y_continuous(breaks = seq(0, 100000, by = 10000))+
      theme_bw()+
      theme(
        axis.text = element_text(size=axis_text_size),
        axis.title  = element_text(size=axis_title_size),
        legend.position = "bottom",
        legend.text=element_text(size=legend_text_size),
        plot.title = element_text(size=plot_title_size),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave("figures/influenza/severe.png",h=12,w=15)

# Zooom in

out_df <- as.data.frame(out) %>%
  gather(., comp, prop,S:D) %>%
  filter(comp %in% c("H", "D"))

col_c  <- c("#EB7746FF", "#1A1F2BFF")

ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "Influenza - Mild Scenario",
       x = "Days",
       y = "Individuals") +
  scale_color_manual("",
                     breaks=c("H","D"),
                     labels=c("Hospitalisation","Deceased"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000))+
  theme_bw()+
  theme(
    axis.text = element_text(size=axis_text_size+10),
    axis.title  = element_text(size=axis_title_size+10),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size+10),
    plot.title = element_text(size=plot_title_size+10),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggsave("figures/influenza/mild_zoom.png",h=10,w=10)


ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  labs(title = "Influenza - Severe Scenario",
       x = "Days",
       y = "Individuals")+
  scale_color_manual("",
                     breaks=c("H","D"),
                     labels=c("Hospitalisation","Deceased"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000))+
  theme_bw()+
    theme(
      axis.text = element_text(size=axis_text_size+10),
      axis.title  = element_text(size=axis_title_size+10),
      legend.position = "bottom",
      legend.text=element_text(size=legend_text_size+10),
      plot.title = element_text(size=plot_title_size+10),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank())

ggsave("figures/influenza/severe_zoom.png",h=10,w=10)


