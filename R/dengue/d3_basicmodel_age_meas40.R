# Load required libraries
library(deSolve)
library(tidyverse)

options(scipen=999)
# "#FF0066"
col_c  <- c("#7030A0","#1EA7C4","#F6A1CAFF","#EB7746FF","#B9CA5DFF", "#1A1F2BFF")
col_a  <- c("#FBDB82FF","#24719BFF","#6A1639FF")
col_sun <- c("#F6A1CAFF", "#1A1F2BFF","#EB7746FF", "#B9CA5DFF","#FBDB82FF", "#9F272EFF", "#D36B8DFF", "#EC9E58FF", "#24719BFF", "#6A1639FF")

lwd_size <- 3
size_plot <- 22
size_strip <- 22
axis_text_size <- 22
axis_title_size <- 22
plot_title_size <- 22
legend_text_size <- 25
ann_text_size <- 5

# Age groups
age_groups <- c("young", "adult", "elderly")
age_dist <- c(0.2, 0.6, 0.2)  # proportions of N_h

# Parameters mild
# beta_mh <- c(0.35, 0.3, 0.25)      # Mosquito -> Human
# beta_hm <- c(0.35, 0.3, 0.25)     # Human -> Mosquito
# CFR <- c(0.0000004, 0.000006, 0.00002)# cfr human
# p_hosp <- c(0.005, 0.01, 0.015) # Proportion hospitalized
# p_die_hosp <- c(0.005, 0.01, 0.015) # Death in hospital
# incub_period_m <- 15 # incubation mosquito
# hmr <- 8 # human mosquito factor

# Parameters sever
# # 
beta_mh <- c(0.55, 0.5, 0.45)      # Mosquito -> Human
beta_hm <- c(0.55, 0.5, 0.45)     # Human -> Mosquito
CFR <-  c(0.002, 0.005, 0.008) #  cfr human
p_hosp <- c(0.05, 0.15, 0.25) # Proportion hospitalized
p_die_hosp <- c(0.01, 0.1, 0.2) # Death in hospital
hmr <- 15 # human mosquito factor
incub_period_m <- 8 # incubation mosquito

# for both
incub_period_h <- 5 # incubation human
infection_period <- 5 # infection human
Hr <-  7 # times in hospital
mu_m = 1 / 14      # Mosquito death
N_h <- 100000

# Tc <- 365   # yearly cycle (days)

N_age <- N_h * age_dist
dt_pop <- data.frame(N_age) %>%
  mutate(
    age =row.names(.)) %>%
  rename(pop = N_age)

# calculate from given assumptions

sigma_h = 1 / incub_period_h   # Incubation (Human)
sigma_m = 1 / incub_period_m  # Incubation (Mosquito)
gamma_h = 1 / infection_period    # Recovery (Human)
N_m = N_h * hmr # number of mosquitos

rho <- (1 - p_die_hosp) / Hr  # Recovery from H
mu_hosp <- p_die_hosp / Hr     # Death in H

# Initial state
initial_state <- c(
  S_h_1 = N_age[1], E_h_1 = 0, I_h_1 = 0, H_h_1 = 0, R_h_1 = 0, D_h_1 = 0,
  S_h_2 = N_age[2], E_h_2 = 0, I_h_2 = 0, H_h_2 = 0, R_h_2 = 0, D_h_2 = 0,
  S_h_3 = N_age[3], E_h_3 = 0, I_h_3 = 0, H_h_3 = 0, R_h_3 = 0, D_h_3 = 0,
  S_m = N_m - 1, E_m = 0, I_m = 1
)


# Gradual behavioral/NPI effect
get_intervention_factor <- function(time) {
  if (time < 40) return(1.0)
  else return(1 - 0.3* (1 - exp(-0.1 * (time - 40))))
}

# Gradual mosquito birth reduction (source reduction)
get_mosquito_reduction_factor <- function(time) {
  if (time < 40) return(1.0)
  else return(1 - 0.3* (1 - exp(-0.1 * (time - 40))))
}


# Model
dengue_model <- function(t, state, parameters) {
  with(as.list(state), {
    # Total population (needed for λ_m)
    I_h_total <- I_h_1 + I_h_2 + I_h_3
    N_h_total <- sum(N_age)
    
    # Forces of infection
    # lambda_h <- beta_mh * I_m / N_m
    # lambda_m <- beta_hm * I_h_total / N_h_total
    
    # Human dynamics per age group
    beta_mh_eff <- beta_mh * get_intervention_factor(t)
    beta_hm_eff <- beta_hm * get_intervention_factor(t)
    
    lambda_h1 <-beta_mh_eff[1] * I_m / N_m
    dS_h_1 <- -lambda_h1 * S_h_1
    dE_h_1 <- lambda_h1 * S_h_1 - sigma_h * E_h_1
    dI_h_1 <- sigma_h * E_h_1 - gamma_h * I_h_1
    dH_h_1 <- p_hosp[1] * gamma_h * I_h_1 - (rho[1] + mu_hosp[1]) * H_h_1
    dR_h_1 <- (1 - p_hosp[1]) * gamma_h * I_h_1 + rho[1] * H_h_1
    dD_h_1 <- CFR[1] * gamma_h * I_h_1 + mu_hosp[1] * H_h_1
    
    lambda_h2 <- beta_mh_eff[2] * I_m / N_m
    dS_h_2 <- -lambda_h2 * S_h_2
    dE_h_2 <- lambda_h2 * S_h_2 - sigma_h * E_h_2
    dI_h_2 <- sigma_h * E_h_2 - gamma_h * I_h_2
    dH_h_2 <- p_hosp[2] * gamma_h * I_h_2 - (rho[2] + mu_hosp[2]) * H_h_2
    dR_h_2 <- (1 - p_hosp[2]) * gamma_h * I_h_2 + rho[2] * H_h_2
    dD_h_2 <- CFR[2] * gamma_h * I_h_2 + mu_hosp[2] * H_h_2
    
    lambda_h3 <- beta_mh_eff[3] * I_m / N_m
    dS_h_3 <- -lambda_h3 * S_h_3
    dE_h_3 <- lambda_h3 * S_h_3 - sigma_h * E_h_3
    dI_h_3 <- sigma_h * E_h_3 - gamma_h * I_h_3
    dH_h_3 <- p_hosp[3] * gamma_h * I_h_3 - (rho[3] + mu_hosp[3]) * H_h_3
    dR_h_3 <- (1 - p_hosp[3]) * gamma_h * I_h_3 + rho[3] * H_h_3
    dD_h_3 <- CFR[3] * gamma_h * I_h_3 + mu_hosp[3] * H_h_3
    
    # Mosquito dynamics
    # Calculate force of infection for mosquitoes:
    lambda_m <- ( beta_hm_eff[1] * I_h_1 +  beta_hm_eff[2] * I_h_2 +  beta_hm_eff[3] * I_h_3) / N_h_total
    
    # Mosquito dynamics
    total_mosquitoes <- S_m + E_m + I_m
    birth_rate <- mu_m * total_mosquitoes* get_mosquito_reduction_factor(t)
    # birth_rate <- mu_m * total_mosquitoes* get_mosquito_reduction_factor(t)
    
    dS_m <- birth_rate - lambda_m * S_m - mu_m * S_m
    dE_m <- lambda_m * S_m - sigma_m * E_m - mu_m * E_m
    dI_m <- sigma_m * E_m - mu_m * I_m
    
    # 40% of the mosquitos are directly killed by thermal fogging
    if (abs(t - 40) < 0.5) {
      dS_m <- dS_m - 0.4 * S_m
      dE_m <- dE_m - 0.4 * E_m
      dI_m <- dI_m - 0.4 * I_m
    }
    
    list(c(
      dS_h_1, dE_h_1, dI_h_1, dH_h_1, dR_h_1, dD_h_1,
      dS_h_2, dE_h_2, dI_h_2, dH_h_2, dR_h_2, dD_h_2,
      dS_h_3, dE_h_3, dI_h_3, dH_h_3, dR_h_3, dD_h_3,
      dS_m, dE_m, dI_m
    ))
  })
}


# Run simulation

times <- seq(0, 400, by = 1)
out <- ode(y = initial_state, times = times, func = dengue_model, parms = NULL)
out_df <- as.data.frame(out)

dt_h <- out_df[, c(1,grep("_h_", names(out_df)))] 

dt_hd <- dt_h[, c(1,grep("D_", names(dt_h)))] %>%
  gather(., comp, prop,2:4) 

dt_hi <- dt_h[, c(1,grep("I_", names(dt_h)))] %>%
  gather(., comp, prop,2:4) 

dt_hh <- dt_h[, c(1,grep("H_", names(dt_h)))] %>%
  gather(., comp, prop,2:4) 


dt <- rbind(dt_hd, dt_hh, dt_hi) %>%
  
  mutate(
    age = str_sub(comp, 5)
  ) %>%
  left_join(dt_pop) %>%
  mutate(
    age = as.factor(age),
    age = case_when(age == "1" ~ "0-18",
                    age == "2" ~ "19-64",
                    age == "3" ~ ">=65"),
    age  = factor(age , levels=c("0-18","19-64",">=65")),
    fac = str_sub(comp,1, 1),
    fac = as.factor(fac),
    fac = case_when(fac == "D" ~ "Deceased",
                    fac == "I" ~ "Infected",
                    fac == "H" ~ "Hospitalisation"),
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased")),
    mx = prop/pop *100000
  )

df1 <- data.frame(time= 10,
                  prop = c(0.1,0.0012, 0.00018), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "no NPI") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


df2 <- data.frame(time= 250,
                  prop = c(0.1,0.0012, 0.00018), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "integrated vector management") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


ggplot(dt, aes(x = time)) +
  facet_wrap(~fac, ncol=3, scales = "free_y") +
  annotate("rect",xmin=0,xmax=40,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
  geom_text(data = df1,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=40,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df2,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  geom_line(aes(y = prop, color = age), lwd=lwd_size) +
  labs(title = "Dengue - Mild Scenario - 40 days",
       x = "Days",
       y = "Individuals")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_a) +
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

ggsave("figures/dengue/mild_age_meas40.png",h=8,w=20)


df1 <- data.frame(time= 10,
                  prop = c(15,3, 5.5), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "no NPI") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


df2 <- data.frame(time= 250,
                  prop = c(15,3, 5.5), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "integrated vector management") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


ggplot(dt, aes(x = time)) +
  facet_wrap(~fac, ncol=3, scales = "free_y") +
  annotate("rect",xmin=0,xmax=40,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
  geom_text(data = df1,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=40,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df2,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  geom_line(aes(y = prop, color = age), lwd=lwd_size) +
  labs(title = "Dengue - Severe Scenario - 40 days",
       x = "Days",
       y = "Individuals")+
  scale_color_manual("",
                     # breaks=c("young","adult","elderly"),
                     # labels=c("0-18","19-64",">=65"),
                     values = col_a) +
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

ggsave("figures/dengue/severe_age_meas40.png",h=8,w=20)

