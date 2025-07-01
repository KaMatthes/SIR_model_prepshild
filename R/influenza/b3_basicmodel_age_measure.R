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

library(deSolve)

# Age groups
age_groups <- c("young", "adult", "elderly")

# Base contact matrix
contact_matrix <- matrix(c(
  8, 4, 1,
  3, 6, 2,
  1, 2, 3
), nrow = 3, byrow = TRUE)
rownames(contact_matrix) <- age_groups
colnames(contact_matrix) <- age_groups

# Parameters (mild)
# R0_vec <- c(young = 1.8, adult = 1.5, elderly = 1.2)
# CFR_vec <- c(young = 0.00002, adult = 0.0002, elderly = 0.005)
# p_hosp_vec <- c(young = 0.001, adult = 0.002, elderly = 0.004)
# p_die_hosp_vec <- c(young = 0.002, adult = 0.15, elderly = 0.3)

# Parameters (severe)

R0_vec <- c(young = 3.5, adult = 3, elderly = 2.5)
CFR_vec <- c(young = 0.0001, adult = 0.0003, elderly = 0.02)
p_hosp_vec <- c(young = 0.01, adult = 0.03, elderly = 0.05)
p_die_hosp_vec <- c(young = 0.02, adult = 0.25, elderly = 0.5)

sigma <- 1 / 2     # Incubation rate
gamma <- 1 / 5     # Recovery rate

# Total population by age group
N_age <- c(young = 20001, adult = 60000, elderly = 20000)

# Beta values per age group
contact_row_sums <- rowSums(contact_matrix)
beta_vec <- R0_vec * gamma / contact_row_sums

# Initial state
init_state <- c(
  S_young = 20000, E_young = 0, I_young = 1, H_young = 0, R_young = 0, D_young = 0,
  S_adult = 64000, E_adult = 0, I_adult = 0, H_adult = 0, R_adult = 0, D_adult = 0,
  S_elderly = 20000, E_elderly = 0, I_elderly = 0, H_elderly = 0, R_elderly = 0, D_elderly = 0
)

# NPI function
# get_intervention_factor <- function(time) {
#   if (time < 30) {
#     return(1.0)    # No intervention
#   } else if (time %in% 30:60) {
#     return(0.2)    # Strong lockdown
#   } else if (time %in% 61:120) {
#     return(0.6)    # Partial reopening
#   } else {
#     return(1.0)    # Restrictions lifted
#   }
# }

# get_intervention_factor <- function(time) {
#   if (time < 30) {
#     return(1.0)    # No intervention
#   } else if (time %in% 30:60) {
#     return(0.2)    # Strong lockdown
#   } else if (time %in% 61:90) {
#     return(0.5)    # Partial reopening
#   } else if (time %in% 91:120) {
#     return(0.2)    # Strong lockdown
#   } else if (time %in% 121:150) {
#     return(0.2)   # Partial reopening
#   } else if (time %in% 151:250) {
#     return(0.7)    # wearing mak
#   } else {
#     return(1)    # Restrictions lifted
#   }
# }

get_intervention_factor <- function(time) {
  if (time < 30) {
    return(1.0)    # No intervention
  } else if (time >=30 & time < 90) {
    return(0.2)    # Strong lockdown
  } else if (time >=90 & time < 150) {
    return(0.5)    # Partial reopening
  } else if (time >=150 & time < 180) {
    return(0.2)    # Strong lockdown
  } else if (time >=180 &  time < 210) {
    return(0.5)   # Partial reopening
  # } else if (time >=210 &  time < 240) {
  #   return(0.7)    # wearing mask
  } else {
    return(0.7)     # wearing mask
  }
}


# SEIR-HDR model with NPIs
seird_model_age <- function(time, state, parameters) {
  with(as.list(state), {
    d_state <- numeric(length(state))
    names(d_state) <- names(state)
    
    # Infectious individuals
    I_vec <- sapply(age_groups, function(age) state[paste0("I_", age)])
    
    # Adjusted contact matrix with NPI effect
    contact_matrix_eff <- contact_matrix * get_intervention_factor(time)
    
    # Force of infection λ_i = sum_j β_i * C_ij * I_j / N_j
    lambda <- sapply(1:3, function(j) sum(contact_matrix_eff[j, ] * I_vec / N_age))
    
    for (i in 1:3) {
      age <- age_groups[i]
      S <- state[paste0("S_", age)]
      E <- state[paste0("E_", age)]
      I <- state[paste0("I_", age)]
      H <- state[paste0("H_", age)]
      
      CFR <- CFR_vec[age]
      beta <- beta_vec[age]
      p_hosp <- p_hosp_vec[age]
      p_die_hosp <- p_die_hosp_vec[age]
      
      mu <- p_die_hosp / 10       # death rate in hospital
      rho <- (1 - p_die_hosp) / 10  # recovery rate in hospital
      
      dS <- -beta * S * lambda[i]
      dE <-  beta * S * lambda[i] - sigma * E
      dI <-  sigma * E - gamma * I
      dH <- p_hosp * gamma * I - (rho + mu) * H
      dR <- (1 - p_hosp) * gamma * I + rho * H
      dD <- CFR * gamma * I + mu * H
      
      d_state[paste0("S_", age)] <- dS
      d_state[paste0("E_", age)] <- dE
      d_state[paste0("I_", age)] <- dI
      d_state[paste0("H_", age)] <- dH
      d_state[paste0("R_", age)] <- dR
      d_state[paste0("D_", age)] <- dD
    }
    
    return(list(d_state))
  })
}

# Run simulation
times <- seq(0, 300, by = 1)
out <- ode(y = init_state, times = times, func = seird_model_age, parms = NULL)
out_df <- as.data.frame(out)


dt_d <- out_df[, c(1,grep("D_", names(out_df)))] %>%
  gather(., comp, prop,2:4) 

dt_i <- out_df[, c(1,grep("I_", names(out_df)))] %>%
  gather(., comp, prop,2:4) 

dt_h <- out_df[, c(1,grep("H_", names(out_df)))] %>%
  gather(., comp, prop,2:4) 

dt_pop <- data.frame(N_age) %>%
  mutate(
    age =row.names(.)) %>%
  rename(pop = N_age)

dt <- rbind(dt_d, dt_h, dt_i) %>%
  mutate(
    age = str_sub(comp, 3),
         age = as.factor(age),
         fac = str_sub(comp,1, 1),
         fac = as.factor(fac)
    ) %>%
  mutate(
    fac = case_when(fac == "D" ~ "Deceased",
                         fac == "I" ~ "Infected",
                         fac == "H" ~ "Hospitalisation"),
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased")),
    ) %>%
  left_join(dt_pop) %>%
  mutate(
    mx = prop/pop *100000
  )

df1 <- data.frame(time= 10,
           prop = c(5.9, 0.016, 0.0325), 
           fac=c("Infected","Hospitalisation","Deceased"),
           text = "no NPI") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


df2 <- data.frame(time= 60,
                  prop = c(5.9, 0.016, 0.0325), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "strong lockdown") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df3 <- data.frame(time= 120,
                  prop = c(5.9, 0.016, 0.0325), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "partial reopening") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df4 <- data.frame(time= 165,
                  prop = c(5, 0.015, 0.032), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "strong lockdown") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df5 <- data.frame(time= 195,
                  prop = c(5, 0.015, 0.032), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "partial reopening") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df6 <- data.frame(time= 255,
                  prop = c(5, 0.015, 0.032), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "all open, but still wearing mask") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )
  
ggplot(dt, aes(x = time)) +
  facet_wrap(~fac, ncol=3, scales = "free_y") +
  annotate("rect",xmin=0,xmax=30,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
  geom_text(data = df1,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=30,xmax=90,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey5") +
  geom_text(data = df2,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=90,xmax=150,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df3,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=150,xmax=180,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey5") +
  geom_text(data = df4,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=180,xmax=210,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df5,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=210,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey70") +
  geom_text(data = df6,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  geom_line(aes(y = prop, color = age), lwd=lwd_size) +
  labs(title = "Influenza - Mild Scenario - Individuals",
       x = "Days",
       y = "Individuals")+
  scale_color_manual("",
                     breaks=c("young","adult","elderly"),
                     labels=c("0-18","19-64",">=65"),
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

ggsave("figures/influenza/mild_age_meas.png",h=8,w=20)

# 
# ggplot(dt, aes(x = time)) +
#   facet_wrap(~fac, ncol=3, scales = "free_y") +
#   annotate("rect",xmin=0,xmax=30,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
#   geom_text(data = df1,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=30,xmax=90,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey5") +
#   geom_text(data = df2,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=90,xmax=150,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey40") +
#   geom_text(data = df3,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=150,xmax=180,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey5") +
#   geom_text(data = df4,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=180,xmax=210,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey40") +
#   geom_text(data = df5,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=210,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey70") +
#   geom_text(data = df6,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   geom_line(aes(y = prop, color = age), lwd=lwd_size) +
#   labs(title = "Influenza - Mild Scenario - Rate",
#        x = "Days",
#        y = "Rate per 100'000")+
#   scale_color_manual("",
#                      breaks=c("young","adult","elderly"),
#                      labels=c("0-18","19-64",">=65"),
#                      values = col_a) +
#   # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
#   theme_bw() +
#   theme(
#     strip.text = element_text(size=size_plot),
#     axis.text = element_text(size=axis_text_size),
#     axis.title  = element_text(size=axis_title_size),
#     legend.position = "bottom",
#     legend.text=element_text(size=legend_text_size),
#     plot.title = element_text(size=plot_title_size),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank())
# 
# ggsave("figures/influenza/mild_age_meas_mx.png",h=8,w=20)

df1 <- data.frame(time= 15,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "no NPI") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


df2 <- data.frame(time= 60,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "strong lockdown") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df3 <- data.frame(time= 120,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "partial reopening") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df4 <- data.frame(time= 165,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "strong lockdown") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df5 <- data.frame(time= 195,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "partial reopening") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )

df6 <- data.frame(time= 220,
                  prop = c(3000, 155, 400), 
                  fac=c("Infected","Hospitalisation","Deceased"),
                  text = "all open, but still wearing mask") %>%
  mutate(
    fac  = factor(fac , levels=c("Infected","Hospitalisation","Deceased"))
  )


ggplot(dt, aes(x = time)) +
  facet_wrap(~fac, ncol=3, scales = "free_y") +
  annotate("rect",xmin=0,xmax=30,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
  geom_text(data = df1,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=30,xmax=90,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey5") +
  geom_text(data = df2,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=90,xmax=150,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df3,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=150,xmax=180,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey5") +
  geom_text(data = df4,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=180,xmax=210,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey40") +
  geom_text(data = df5,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  annotate("rect",xmin=210,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey70") +
  geom_text(data = df6,
            aes(x=time, y=prop,label = text),
            angle = 90,
            size=ann_text_size) +
  geom_line(aes(y = prop, color = age), lwd=lwd_size) +
  labs(title = "Influenza - Severe Scenario - Individuals",
       x = "Days",
       y = "Individuals")+
  scale_color_manual("",
                     breaks=c("young","adult","elderly"),
                     labels=c("0-18","19-64",">=65"),
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

ggsave("figures/influenza/severe_meas_age.png",h=8,w=20)


# ggplot(dt, aes(x = time)) +
#   facet_wrap(~fac, ncol=3, scales = "free_y") +
#   annotate("rect",xmin=0,xmax=30,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey90") +
#   geom_text(data = df1,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=30,xmax=90,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey5") +
#   geom_text(data = df2,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=90,xmax=150,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey40") +
#   geom_text(data = df3,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=150,xmax=180,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey5") +
#   geom_text(data = df4,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=180,xmax=210,ymin=-Inf,ymax=Inf,alpha=0.3,fill="grey40") +
#   geom_text(data = df5,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   annotate("rect",xmin=210,xmax=Inf,ymin=-Inf,ymax=Inf,alpha=0.2,fill="grey70") +
#   geom_text(data = df6,
#             aes(x=time, y=prop,label = text),
#             angle = 90) +
#   geom_line(aes(y = prop, color = age), lwd=lwd_size) +
#   labs(title = "Influenza - Severe Scenario - Rate",
#        x = "Days",
#        y = "Rate per 100'000")+
#   scale_color_manual("",
#                      breaks=c("young","adult","elderly"),
#                      labels=c("0-18","19-64",">=65"),
#                      values = col_a) +
#   # scale_y_continuous(breaks = seq(0, 100000, by = 10000)) +
#   theme_bw() +
#   theme(
#     strip.text = element_text(size=size_plot),
#     axis.text = element_text(size=axis_text_size),
#     axis.title  = element_text(size=axis_title_size),
#     legend.position = "bottom",
#     legend.text=element_text(size=legend_text_size),
#     plot.title = element_text(size=plot_title_size),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank())
# 
# ggsave("figures/influenza/sever_age_meas_mx.png",h=8,w=20)
