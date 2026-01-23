###############################
# LIBRARIES
###############################
pacman::p_load(
  GillespieSSA2, dplyr, ggplot2, zoo, plotly
)

###############################
# CLIMATE DATA
###############################
source("mombasa_climate_data.R")

days <- 1825
days_per_month <- 30

temp_scaled <- (mombasa_climate$mean_temp - min(mombasa_climate$mean_temp)) /
  (max(mombasa_climate$mean_temp) - min(mombasa_climate$mean_temp))

temp_daily <- rep(temp_scaled, each = days_per_month)[1:days]

###############################
# CLIMATE-DRIVEN PARAMETERS
###############################

Lambda_m0 <- 20000
mu_m0     <- 1/10
sigma_m0  <- 1/5

Lambda_m_daily <- Lambda_m0 * (0.3 + temp_daily)
mu_m_daily     <- mu_m0     * (1.2 - temp_daily)
sigma_m_daily  <- sigma_m0 * (0.5 + temp_daily)

###############################
# INITIAL STATE
###############################
initial_state <- c(
  S1 = 1400000 - 532,
  E1 = 0,
  S2 = 50000 - 35,
  E2 = 0,
  Im = 503,
  Is = 64,
  C  = 0,
  R  = 0,
  
  Sm   = 10000,
  Em   = 500,
  Im_m = 500,
  
  CumulativeCases = 567,
  ImportedCases   = 0
)

###############################
# PARAMETERS (STATIC)
###############################
params <- c(
  mu = 1/(62*365),
  sigma = 1/3,
  gamma_m = 1/3,
  gamma_s = 1/7,
  gamma_c = 1/30,
  
  severe_prob_1 = 0.23,
  severe_prob_2 = 0.49,
  chronic_prob_s = 0.44,
  
  contact_rate_1 = 0.7,
  contact_rate_2 = 0.3,
  
  mu_m = mu_m0,
  sigma_m = sigma_m0,
  
  b = 0.8,
  beta_hm = 0.442,
  
  alpha = 1/(59*365),
  Lambda_m = Lambda_m0,
  
  imported_rate = 1,
  beta_mh = 0.333
)

###############################
# REACTIONS
###############################
reactions <- list(
  
  ## ---------------- HUMAN INFECTION ----------------
  reaction(
    ~ (b * beta_mh * Im_m / (S1 + E1 + S2 + E2 + Im + Is + R + C)) *
      contact_rate_1 * S1,
    c(S1 = -1, E1 = +1, CumulativeCases = +1),
    "S1_to_E1"
  ),
  
  reaction(
    ~ imported_rate * 0.7 * S1 / (S1 + S2),
    c(S1 = -1, E1 = +1, ImportedCases = +1, CumulativeCases = +1),
    "imported_S1"
  ),
  
  reaction(~ sigma * (1 - severe_prob_1) * E1,
           c(E1 = -1, Im = +1),
           "E1_to_Im"),
  
  reaction(~ sigma * severe_prob_1 * E1,
           c(E1 = -1, Is = +1),
           "E1_to_Is"),
  
  reaction(~ alpha * S1,
           c(S1 = -1, S2 = +1),
           "aging"),
  
  reaction(
    ~ (b * beta_mh * Im_m / (S1 + E1 + S2 + E2 + Im + Is + R + C)) *
      contact_rate_2 * S2,
    c(S2 = -1, E2 = +1, CumulativeCases = +1),
    "S2_to_E2"
  ),
  
  reaction(
    ~ imported_rate * 0.3 * S2 / (S1 + S2),
    c(S2 = -1, E2 = +1, ImportedCases = +1, CumulativeCases = +1),
    "imported_S2"
  ),
  
  reaction(~ sigma * (1 - severe_prob_2) * E2,
           c(E2 = -1, Im = +1),
           "E2_to_Im"),
  
  reaction(~ sigma * severe_prob_2 * E2,
           c(E2 = -1, Is = +1),
           "E2_to_Is"),
  
  reaction(~ gamma_m * Im, c(Im = -1, R = +1), "Im_to_R"),
  reaction(~ gamma_s * Is * (1 - chronic_prob_s), c(Is = -1, R = +1), "Is_to_R"),
  reaction(~ gamma_s * Is * chronic_prob_s, c(Is = -1, C = +1), "Is_to_C"),
  reaction(~ gamma_c * C, c(C = -1, R = +1), "C_to_R"),
  
  ## ---------------- MOSQUITOES ----------------
  reaction(~ Lambda_m, c(Sm = +1), "mosq_birth"),
  
  reaction(
    ~ (b * beta_hm * (Im + Is) / (S1 + E1 + S2 + E2 + Im + Is + R + C)) * Sm,
    c(Sm = -1, Em = +1),
    "Sm_to_Em"
  ),
  
  reaction(~ sigma_m * Em, c(Em = -1, Im_m = +1), "Em_to_Im_m"),
  
  reaction(~ mu_m * Sm, c(Sm = -1), "Sm_death"),
  reaction(~ mu_m * Em, c(Em = -1), "Em_death"),
  reaction(~ mu_m * Im_m, c(Im_m = -1), "Im_m_death")
)

###############################
# COMPILE
###############################
compiled_reactions <- compile_reactions(
  reactions = reactions,
  state_ids = names(initial_state),
  params = params
)

###############################
# HYBRID SSA (DAILY STEPPING)
###############################
set.seed(1000)

state <- initial_state
output <- list()

for (d in 1:days) {
  
  
  params["Lambda_m"] <- Lambda_m_daily[d]
  params["mu_m"]     <- mu_m_daily[d]
  params["sigma_m"]  <- sigma_m_daily[d]
  
  sim <- ssa(
    initial_state = state,
    reactions = compiled_reactions,
    params = params,
    method = ssa_exact(),
    final_time = 1,
    census_interval = 1
  )
  
  state <- as.numeric(tail(sim$state, 1))
  names(state) <- colnames(sim$state)
  
  output[[d]] <- data.frame(time = d, t(state))
}

sim_df <- bind_rows(output)

###############################
# PLOTS
###############################
#Plot the infectious humans
ggplot(sim_df, aes(x = time)) +
  geom_line(aes(y = Im+Is, color = "Infectious"), size = 1)+
  labs(
    x = "Time (days)",
    y = "Number of Infectious Humans",
    title = "Infectious Human Population Over Time"
  ) 

#Plot the human population dynamics
ggplot(sim_df, aes(x=time))+
  geom_line(aes(y=E1+E2, color="Exposed"), size=1)+
  geom_line(aes(y=C, color="Chronic"), size=1)+
  geom_line(aes(y=R, color="Recovered"), size=1)+
  labs(
    x="Time(days)",
    y="Population",
    title="Human Population Dynamics")

#Plot the s1
ggplot(sim_df, aes(x = time)) +
  geom_line(aes(y = S1, color = "S1"), size = 1)+
  labs(
    x = "Time (days)",
    y = "Number of Susceptible Humans",
    title = "Susceptible Human Population Over Time"
  ) +
  scale_color_manual(values = c("S1" = "blue"), 
                     name = "Compartments") +
  theme_minimal(base_size = 14)

#Plot the s2 group over time
ggplot(sim_df, aes(x = time)) +
  geom_line(aes(y = S2, color = "S2"), size = 1) +
  labs(
    x = "Time (days)",
    y = "Number of Susceptible Humans",
    title = "Susceptible Human Population Over Time"
  ) +
  scale_color_manual(values = c("S2" = "green"), 
                     name = "Compartments") +
  theme_minimal(base_size = 14)

#Plot the mosquito population ( Sm, Em, Im_m) over time separately
ggplot(sim_df, aes(x = time)) +
  geom_line(aes(y = Sm, color = "Sm"), size = 1) +
  geom_line(aes(y = Em, color = "Em"), size = 1) +
  geom_line(aes(y = Im_m, color = "Im_m"), size = 1) +
  labs( 
    x = "Time (days)",
    y = "Number of Mosquitoes",
    title = "Mosquito Population Over Time"
  ) +
  scale_color_manual(values = c("Sm" = "orange", "Em" = "purple", "Im_m" = "red"), 
                     name = "Compartments") +
  theme_minimal(base_size = 14)


