rm(list=ls())
library(deSolve)
library(tidyverse)
library(openxlsx)
library(DT)

options(scipen=999)



# Age groups
age_groups <- c("young", "adult", "elderly")

# Contact matrix
contact_matrix <- matrix(c(
  8, 4, 1,
  3, 6, 2,
  1, 2, 3
), nrow = 3, byrow = TRUE)
rownames(contact_matrix) <- age_groups
colnames(contact_matrix) <- age_groups

# Parameters
R0_vec <- c(young = 2.5, adult = 3, elderly = 3.5)  # Basic reproduction number
CFR_vec <- c(young = 0.0001, adult = 0.0003, elderly = 0.02) # Case fatality rate
p_hosp_vec <- c(young = 0.01, adult = 0.03, elderly = 0.05) # Proportion of infections hospitalized
p_die_hosp_vec <- c(young = 0.02, adult = 0.25, elderly = 0.5)  # Proportion of hospitalized who die
p_die_icu_vec <- c(young = 0.10, adult = 0.25, elderly = 0.5) # proportion of icu who die

Ht <- 5    # avg hospital stay (days)
Ct <- 10   # avg ICU stay (days)
p_icu <- 0.15  # fraction of hospitalized who go to ICU
Hosp_day <- 5 # how many days in hospital (no ICU)
ICU_day <- 10 # how many days in ICU
ICU_prop <- 0.15 # proportion of icu
sigma <- 1 / 2     # Incubation rate
gamma <- 1 / 5     # Recovery rate 
out_c_rate <- 1 / Ct # ICU recovery rate

# Total population by age group
N_age <- c(young = 20001, adult = 60000, elderly = 20000)

# Beta values per age group
contact_row_sums <- rowSums(contact_matrix) 
beta_vec <- R0_vec * gamma / contact_row_sums # transmission probability per contact per unit time for age group i

# period invervention, examples
period1 <- 25
period2 <- 50
period3 <- 100

# Intervention reduction factors
int1 <- 0.2
int2 <- 0.5
int3 <- 0.7

# period vaccination, examples
pervac1 <- 100
pervac2 <- 150
pervac3 <- 200

# Vaccination rates per day, needs to adjusted
vacc_young <- c(0, 0, 1000)
vacc_adult <- c(0, 1000, 1000)
vacc_elderly <- c(1000, 1000, 1000)

v1 <- as.integer(c(vacc_young[1],vacc_adult[1],vacc_elderly[1]))
v2 <- as.integer(c(vacc_young[2],vacc_adult[2],vacc_elderly[2]))
v3 <- as.integer(c(vacc_young[3],vacc_adult[3],vacc_elderly[3]))


# Initial state
init_state <- c(
  S_young=as.integer(N_age["young"])-1, E_young=0, I_young=1, H_young=0, C_young=0, R_young=0, D_young=0,
  S_adult=as.integer(N_age["adult"]), E_adult=0, I_adult=0, H_adult=0, C_adult=0, R_adult=0, D_adult=0,
  S_elderly=as.integer(N_age["elderly"]), E_elderly=0, I_elderly=0, H_elderly=0, C_elderly=0, R_elderly=0, D_elderly=0
)

get_intervention_factor <- function(time) {
  if (time < period1) 1.0
  else if (time < period2) int1
  else if (time < period3) int2
  else int3
}

get_vaccination_rate <- function(time) {
  if (time >= pervac1 & time < pervac2) c(young=v1[1], adult=v1[2], elderly=v1[3])
  else if (time >= pervac2 & time < pervac3) c(young=v2[1], adult=v2[2], elderly=v2[3])
  else if (time >= pervac3) c(young=v3[1], adult=v3[2], elderly=v3[3])
  else c(young=0, adult=0, elderly=0)
}


# SEIR model

seird_model_age <- function(time, state, parameters) {
  with(as.list(state), {
    d_state <- numeric(length(state)); names(d_state) <- names(state)
    
    I_vec <- sapply(age_groups, function(age) state[paste0("I_", age)])
    contact_matrix_eff <- contact_matrix * get_intervention_factor(time)
    lambda <- sapply(1:3, function(j) sum(contact_matrix_eff[j, ] * I_vec / N_age))
    vac_rates <- get_vaccination_rate(time)
    
    for (i in 1:3) {
      age <- age_groups[i]
      S <- state[paste0("S_", age)] 
      E <- state[paste0("E_", age)] 
      I <- state[paste0("I_", age)] 
      H <- state[paste0("H_", age)] 
      C <- state[paste0("C_", age)] 
      
      CFR <- CFR_vec[age]
      beta <- beta_vec[age]
      p_hosp <- p_hosp_vec[age]
      p_die_hosp <- p_die_hosp_vec[age]
      p_die_icu <- p_die_icu_vec[age]
      
      vac <- min(S, vac_rates[age])
      
      
      # 
      mu_h <- p_die_hosp /Hosp_day
      rho_h <- (1 - p_die_hosp) /Hosp_day
      
      mu_c <- p_die_icu  / ICU_day
      rho_c <- (1 - p_die_icu)/ ICU_day
      
      p_icu <- ICU_prop/100
      p_icu_safe <- ifelse(p_icu >= 1, 0.9999, p_icu)
      phi <- p_icu_safe / (1 - p_icu_safe) * (rho_h + mu_h)
      
      # phi <- p_icu / (1 - p_icu) * (rho_h + mu_h) 
      #     
      # ODEs
      dS <- -beta * S * lambda[i] - vac
      dE <-  beta * S * lambda[i] - sigma * E
      dI <-  sigma * E - gamma * I
      # dH <- p_hosp * gamma * I - (rho_h + mu_h) * H
      # dC <- p_icu * p_hosp * gamma * I -  (rho_c + mu_c)* C
      
      dH <- p_hosp * gamma * I - (rho_h + mu_h + phi) * H  # all hospitalizations
      dC <- phi * H - (rho_c + mu_c) * C                         # ICU subset of H
      
      # phi = rate at which H patients are transferred to ICU to achieve exactly p_icu fraction
      
      
      # Recovered and Deaths: accumulate from compartment outflows and vaccination
      dR <- (1 - p_hosp) * gamma * I + rho_h * H + rho_c * C + vac
      dD <- CFR * gamma * I + mu_h * H + mu_c * C
      
      d_state[paste0("S_", age)] <- dS # susceptible
      d_state[paste0("E_", age)] <- dE # exposed
      d_state[paste0("I_", age)] <- dI # exposed
      d_state[paste0("H_", age)] <- dH # in hospital with out ICU
      d_state[paste0("C_", age)] <- dC # in ICU, total number hospital would be H + C
      d_state[paste0("R_", age)] <- dR # recovered
      d_state[paste0("D_", age)] <- dD # dead
    }
    return(list(d_state))
  })
}

# Run simulation
time <- seq(0, 300, by = 1)
out <- ode(y = init_state, times = time, func = seird_model_age, parms = NULL)
out_df <- as.data.frame(out)

