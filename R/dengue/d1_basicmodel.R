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


# Parameters mild
beta_mh <- 0.3      # Mosquito -> Human
beta_hm <- 0.3     # Human -> Mosquito
CFR <-  0.00001 # cfr human
hmr <- 8 # human mosquito factor
p_hosp = 0.01     # Proportion hospitalized
p_die_hosp = 0.01   # Death in hospital
incub_period_m <- 15 # incubation mosquito

# Parameters sever
# beta_mh <- 0.5      # Mosquito -> Human
# beta_hm <- 0.5     # Human -> Mosquito
# CFR <-  0.005  # cfr human
# hmr <- 15 # human mosquito factor
# p_hosp = 0.15    # Proportion hospitalized
# p_die_hosp = 0.1  # Death in hospital
# incub_period_m <- 8 # incubation mosquito

# for both
incub_period_h <- 5 # incubation human
infection_period <- 5 # infection human
Hr <-  7 # times in hospital
mu_m = 1 / 14      # Mosquito death
N_h <- 100000
# Tc <- 365   # yearly cycle (days)


# calculate from given assumptions

sigma_h = 1 / incub_period_h   # Incubation (Human)
sigma_m = 1 / incub_period_m  # Incubation (Mosquito)
gamma_h = 1 / infection_period    # Recovery (Human)
N_m = N_h * hmr # number of mosquitos

rho <- (1 - p_die_hosp) / Hr  # Recovery from H
mu_hosp <- p_die_hosp / Hr     # Death in H

# Initial state
initial_state <- c(
  S_h =N_h, E_h = 0, I_h = 0, H_h = 0, R_h = 0, D_h = 0,
  S_m = N_m-1, E_m = 0, I_m = 1
)


# Model
dengue_model <- function(t, state, parameters) {
  with(as.list(state), {
    
    # 
    # beta_mh_t <- beta_mh * (1 + a * sin(2 * pi * t / Tc))
    # lambda_h <- beta_mh_t * I_m / N_m
    
    
    # Forces of infection
    lambda_h <- beta_mh * I_m / N_m
    lambda_m <- beta_hm * I_h / N_h  # Only I_h transmits
    
    # Human dynamics
    dS_h <- -lambda_h * S_h
    dE_h <- lambda_h * S_h - sigma_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h
    dH_h <- p_hosp * gamma_h * I_h - (rho + mu_hosp) * H_h
    dR_h <- (1 - p_hosp) * gamma_h * I_h + rho * H_h
    dD_h <- CFR * gamma_h * I_h + mu_hosp * H_h
    
    # Mosquito dynamics
  
    
    # Mosquito population and seasonal birth rate
    total_mosquitoes <- S_m + E_m + I_m
    birth_rate <- mu_m *  total_mosquitoes
    
    # Mosquito dynamics
    dS_m <- birth_rate - lambda_m * S_m - mu_m * S_m
    dE_m <- lambda_m * S_m - sigma_m * E_m - mu_m * E_m
    dI_m <- sigma_m * E_m - mu_m * I_m
    
    
    
    list(c(dS_h, dE_h, dI_h, dH_h, dR_h, dD_h, dS_m, dE_m, dI_m))
    
  })
}

# Run simulation

times <- seq(0, 300, by = 1)
out <- ode(y = initial_state, times = times, func = dengue_model, parms = NULL)
out_df <- as.data.frame(out)

dt_h <- out_df[, c(1,grep("_h", names(out_df)))] %>%
  gather(., comp, prop,2:7)  %>%
  mutate(
    comp = case_when(comp == "D_h" ~ "Deceased",
                    comp== "I_h" ~ "Infected",
                    comp== "H_h" ~ "Hospitalisation",
                    comp== "S_h" ~ "Susceptible",
                    comp== "R_h" ~ "Recovered",
                    comp== "E_h" ~ "Exposed"),
    comp  = factor(comp , levels=c("Susceptible","Exposed","Infected","Hospitalisation","Recovered","Deceased")),
    fac = "human"
  )

dt_m <- out_df[, c(1,grep("_m", names(out_df)))] %>%
  gather(., comp, prop,2:4)  %>%
  mutate(
    comp = case_when(comp== "I_m" ~ "Infected",
                     comp== "S_m" ~ "Susceptible",
                     comp== "E_m" ~ "Exposed"),
    comp  = factor(comp , levels=c("Susceptible","Exposed","Infected")),
    fac = "mosquito"
  )

dt <- rbind(dt_h, dt_m)

ggplot(dt, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  facet_wrap(~fac, ncol=2, scales = "free_y") +
  labs(title = "Dengue - Mild Scenario",
       x = "Days",
       y = "Individuals/Mosquitos")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
  theme_bw() +
  theme(
    strip.text = element_text(size=size_plot),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_title_size),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggsave("figures/dengue/mild.png",h=10,w=20)


ggplot(dt, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  facet_wrap(~fac, ncol=2, scales = "free_y") +
  labs(title = "Dengue - Severe Scenario",
       x = "Days",
       y = "Individuals/Mosquitos")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
  theme_bw() +
  theme(
    strip.text = element_text(size=size_plot),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_title_size),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggsave("figures/dengue/severe.png",h=10,w=20)



# Zooom in

dt_z <-dt %>%
  filter(comp %in% c("Hospitalisation","Deceased"))

col_c  <- c("#EB7746FF", "#1A1F2BFF")

ggplot(dt_z, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  facet_wrap(~fac, ncol=2, scales = "free_y") +
  labs(title = "Dengue - Mild Scenario",
       x = "Days",
       y = "Individuals/Mosquitos")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
  theme_bw() +
  theme(
    strip.text = element_text(size=size_plot),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_title_size),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggsave("figures/dengue/mild_zoom.png",h=10,w=10)

ggplot(dt_z, aes(x = time)) +
  geom_line(aes(y = prop, color = comp), lwd=lwd_size) +
  facet_wrap(~fac, ncol=2, scales = "free_y") +
  labs(title = "Dengue - Severe Scenario",
       x = "Days",
       y = "Individuals/Mosquitos")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_c) +
  # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
  theme_bw() +
  theme(
    strip.text = element_text(size=size_plot),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_title_size),
    legend.position = "bottom",
    legend.text=element_text(size=legend_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())


ggsave("figures/dengue/severe_zoom.png",h=10,w=10)



