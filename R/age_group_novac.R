#library(deSolve)

# Age groups
age_groups <- c("young", "adult", "elderly")

# Contact matrix (contacts per day between groups)
contact_matrix <- matrix(c(
  8.0, 2.5, 0.5,
  2.5, 6.0, 1.5,
  0.5, 1.5, 3.0
), nrow = 3, byrow = TRUE)
rownames(contact_matrix) <- age_groups
colnames(contact_matrix) <- age_groups

# Parameters by age group
R0_vec <- c(young = 2.5, adult = 2.0, elderly = 1.5)
CFR_vec <- c(young = 0.0005, adult = 0.001, elderly = 0.02)

# Initial population per age group
init_state <- c(
  S_young = 50000, E_young = 0, I_young = 1, R_young = 0, D_young = 0,
  S_adult = 40000, E_adult = 0, I_adult = 0, R_adult = 0, D_adult = 0,
  S_elderly = 10000, E_elderly = 0, I_elderly = 0, R_elderly = 0, D_elderly = 0
)

N_age <- c(young = 50000 + 1, adult = 40000, elderly = 10000)

# Model function without hospitalization or vaccination
seird_model_age <- function(time, state, parameters) {
  with(as.list(state), {
    d_state <- numeric(length(state))
    names(d_state) <- names(state)
    
    # Extract infectious and susceptible vectors
    I_vec <- sapply(age_groups, function(age) state[paste0("I_", age)])
    S_vec <- sapply(age_groups, function(age) state[paste0("S_", age)])
    E_vec <- sapply(age_groups, function(age) state[paste0("E_", age)])
    
    D <- 5
    gamma <- 1 / D  # Recovery rate (fixed for all)
    sigma <- 1 / 2  # Incubation rate (fixed for all)
    
    for (i in 1:3) {
      age <- age_groups[i]
      S <- state[paste0("S_", age)]
      E <- state[paste0("E_", age)]
      I <- state[paste0("I_", age)]
      
      R0 <- R0_vec[age]
      beta <- R0 * gamma
      CFR <- CFR_vec[age]
      
      # SEIRD equations without hospitalization or vaccination
      dS <- -beta * S * sigma
      dE <- beta * S * sigma - sigma * E
      dI <- sigma * E - gamma * I
      dR <- (1 - CFR) * gamma * I
      dD <- CFR * gamma * I
      
      d_state[paste0("S_", age)] <- dS
      d_state[paste0("E_", age)] <- dE
      d_state[paste0("I_", age)] <- dI
      d_state[paste0("R_", age)] <- dR
      d_state[paste0("D_", age)] <- dD
    }
    
    return(list(d_state))
  })
}

# Time vector
times <- seq(0, 180, by = 1)

# Run simulation
out <- ode(y = init_state, times = times, func = seird_model_age, parms = NULL)
out_df <- as.data.frame(out)

times <- seq(0, 180, by = 1)
out <- ode(y = init_state, times = times, func = seird_model_age, parms = NULL)
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
