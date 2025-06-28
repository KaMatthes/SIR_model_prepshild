# Age-structured SEIHRD model with contact matrix and time-varying vaccination

library(deSolve)

# Age groups
age_groups <- c("young", "adult", "elderly")

# Contact matrix
contact_matrix <- matrix(c(
  10, 4, 1,
  4, 8, 2,
  1, 2, 5
), nrow = 3, byrow = TRUE)
rownames(contact_matrix) <- age_groups
colnames(contact_matrix) <- age_groups

# Parameters by age group
R0_vec <- c(young = 2.5, adult = 2.0, elderly = 1.5)
CFR_vec <- c(young = 0.0005, adult = 0.001, elderly = 0.02)
p_hosp_vec <- c(young = 0.005, adult = 0.01, elderly = 0.2)
p_die_hosp_vec <- c(young = 0.01, adult = 0.05, elderly = 0.3)

# Vaccination rate function (same for all groups, could be extended)
vaccination_rate <- function(t) {
  max_rate <- 0.01
  mid_time <- 90
  steepness <- 0.1
  rate <- max_rate / (1 + exp(-steepness * (t - mid_time)))
  return(rate)
}

# Time-dependent functions
hospi_rate <- function(t) {
  if (t < 15) return(0.01)
  else if (t <= 30) return(0.005)
  else if (t <= 45) return(0.009)
  else if (t <= 90) return(0.003)
  else return(0.004)
}

p_die_hosp_f <- function(t) if (t < 40) 0.1 else 0.05

cfr <- function(t) if (t < 40) 0.001 else 0.0015

R0_f <- function(t) {
  if (t < 15) return(2.5)
  else if (t <= 30) return(1.0)
  else if (t <= 45) return(1.8)
  else if (t <= 90) return(0.9)
  else return(1.5)
}

# Initial population per age group
init_state <- c(
  S_young = 50000, I_young = 1, H_young = 0, R_young = 0, D_young = 0,
  S_adult = 40000, I_adult = 0, H_adult = 0, R_adult = 0, D_adult = 0,
  S_elderly = 10000, I_elderly = 0, H_elderly = 0, R_elderly = 0, D_elderly = 0
)

N_age <- c(young = 50000 + 1, adult = 40000, elderly = 10000)

# Model function
seihrd_model_age <- function(time, state, parameters) {
  with(as.list(state), {
    d_state <- numeric(length(state))
    names(d_state) <- names(state)
    
    I_vec <- sapply(age_groups, function(age) state[paste0("I_", age)])
    S_vec <- sapply(age_groups, function(age) state[paste0("S_", age)])
    
    lambda <- numeric(3)
    names(lambda) <- age_groups
    
    for (j in 1:3) {
      lambda[j] <- sum(contact_matrix[j, ] * I_vec / N_age)
    }
    
    for (i in 1:3) {
      age <- age_groups[i]
      S <- state[paste0("S_", age)]
      I <- state[paste0("I_", age)]
      H <- state[paste0("H_", age)]
      
      R0 <- R0_f(time)
      gamma <- 1 / 5
      beta <- R0 * gamma
      CFR <- CFR_vec[age]
      p_hosp <- p_hosp_vec[age]
      p_die_hosp <- p_die_hosp_vec[age]
      mu <- p_die_hosp / 10
      rho <- (1 - p_die_hosp) / 10
      v <- vaccination_rate(time)
      
      dS <- -beta * S * lambda[i] - min(v * S, S)
      dI <- beta * S * lambda[i] - (1 - p_hosp - CFR) * gamma * I - p_hosp * gamma * I - CFR * gamma * I
      dH <- p_hosp * gamma * I - (rho + mu) * H
      dR <- (1 - p_hosp - CFR) * gamma * I + rho * H + min(v * S, S)
      dD <- CFR * gamma * I + mu * H
      
      d_state[paste0("S_", age)] <- dS
      d_state[paste0("I_", age)] <- dI
      d_state[paste0("H_", age)] <- dH
      d_state[paste0("R_", age)] <- dR
      d_state[paste0("D_", age)] <- dD
    }
    
    return(list(d_state))
  })
}

# Solve
times <- seq(0, 180, by = 1)
out <- ode(y = init_state, times = times, func = seihrd_model_age, parms = NULL)
out <- as.data.frame(out)


# Quick plot
matplot(out$time, out[, grep("D_", names(out))], type = "l", lty = 1, col = 2:4,
        ylab = "Death", xlab = "Time", main = "Death per Age Group")
legend("topright", legend = age_groups, col = 2:4, lty = 1)


# Quick plot
matplot(out$time, out[, grep("H_", names(out))], type = "l", lty = 1, col = 2:4,
        ylab = "Hospital", xlab = "Time", main = "Hospital per Age Group")
legend("topright", legend = age_groups, col = 2:4, lty = 1)

# Quick plot
matplot(out$time, out[, grep("I_", names(out))], type = "l", lty = 1, col = 2:4,
        ylab = "Infected", xlab = "Time", main = "Infected per Age Group")
legend("topright", legend = age_groups, col = 2:4, lty = 1)
