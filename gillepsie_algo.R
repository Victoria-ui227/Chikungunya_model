library(pacman)
pacman::p_load(deSolve, ggplot2, reshape2, dplyr, gridExtra, GillespieSSA2)

# -------------------------------
# Initial state
# -------------------------------
initial_state <- c(
  S1 = 1400000-532,    
  E1 = 0,                  
  S2 = 50000-35,    
  E2 = 0,                  
  Im = 503,     
  Is = 64,      
  C = 0,        
  R = 0,         
  Sm = 2418,    
  Em = 1413,    
  Im_m = 502,
  CumulativeCases = 567      # NEW: track cumulative infections
)

# -------------------------------
# Parameters
# -------------------------------
params <- c(
  mu = 1/(62*365),           # Human natural death rate
  sigma = 1/3,               # Incubation period
  gamma_m = 1/3,             # Recovery from mild
  gamma_s = 1/7,             # Recovery from severe
  gamma_c = 1/30,            # Recovery from chronic
  severe_prob_1 = 0.23,      
  severe_prob_2 = 0.49,      
  chronic_prob_s = 0.44,     
  contact_rate_1 = 0.7,      
  contact_rate_2 = 0.3,      
  prop_age1_Im = 0.765,      
  prop_age1_Is = 0.003,      
  prop_age1_C = 0.015,       
  prop_age1_R = 0.85,        
  mu_m = 1/10,               
  sigma_m = 1/5,             
  b = 0.41,                 
  beta_hm = 0.442,           
  beta_mh = 0.333,           
  alpha = 1/365,        # Aging S1 -> S2
  Lambda_m = 35.5            # Mosquito recruitment rate
)

# -------------------------------
# Define reactions with explicit formulas
# -------------------------------
reactions <- list(
  # Human births into S1 (using total human population N_h)
  reaction(
    ~ mu * (S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_R * R + 
              S2 + E2 + (1 - prop_age1_Im) * Im + (1 - prop_age1_Is) * Is + 
              prop_age1_C * C + (1 - prop_age1_R) * R),
    c(S1 = +1), 
    "birth_S1"
  ),
  
  # Infection S1 -> E1
  reaction(
    ~ (b * beta_mh * Im_m / (S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_R * R + 
                               S2 + E2 + (1 - prop_age1_Im) * Im + (1 - prop_age1_Is) * Is + 
                               prop_age1_C * C + (1 - prop_age1_R) * R)) * contact_rate_1 * S1,
    c(S1 = -1, E1 = +1, CumulativeCases = +1),
    "S1_to_E1"
  ),
  
  # E1 -> Im (mild)
  reaction(~ sigma * (1 - severe_prob_1) * E1, 
           c(E1 = -1, Im = +1), 
           "E1_to_Im"),
  
  # E1 -> Is (severe)
  reaction(~ sigma * severe_prob_1 * E1, 
           c(E1 = -1, Is = +1), 
           "E1_to_Is"),
  
  # Aging S1 -> S2
  reaction(~ alpha * S1,
           c(S1 = -1, S2 = +1),
           "S1_to_S2"),
  
  # Infection S2 -> E2
  reaction(
    ~ (b * beta_mh * Im_m / (S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_R * R + 
                               S2 + E2 + (1 - prop_age1_Im) * Im + (1 - prop_age1_Is) * Is + 
                               prop_age1_C * C + (1 - prop_age1_R) * R)) * contact_rate_2 * S2,
    c(S2 = -1, E2 = +1, CumulativeCases = +1),
    "S2_to_E2"
  ),
  
  # E2 -> Im
  reaction(~ sigma * (1 - severe_prob_2) * E2, 
           c(E2 = -1, Im = +1), 
           "E2_to_Im"),
  
  # E2 -> Is
  reaction(~ sigma * severe_prob_2 * E2, 
           c(E2 = -1, Is = +1), 
           "E2_to_Is"),
  
  # Im -> R
  reaction(~ gamma_m * Im, 
           c(Im = -1, R = +1), 
           "Im_to_R"),
  
  # Is -> R
  reaction(~ gamma_s * Is * (1 - chronic_prob_s), 
           c(Is = -1, R = +1), 
           "Is_to_R"),
  
  # Is -> C
  reaction(~ gamma_s * Is * chronic_prob_s, 
           c(Is = -1, C = +1), 
           "Is_to_C"),
  
  # C -> R
  reaction(~ gamma_c * C, 
           c(C = -1, R = +1), 
           "C_to_R"),
  
  # Deaths humans
  reaction(~ mu * S1, c(S1 = -1), "death_S1"),
  reaction(~ mu * E1, c(E1 = -1), "death_E1"),
  reaction(~ mu * S2, c(S2 = -1), "death_S2"),
  reaction(~ mu * E2, c(E2 = -1), "death_E2"),
  reaction(~ mu * Im, c(Im = -1), "death_Im"),
  reaction(~ mu * Is, c(Is = -1), "death_Is"),
  reaction(~ mu * R, c(R = -1), "death_R"),
  reaction(~ mu * C, c(C = -1), "death_C"),
  
  # Mosquitoes
  # Birth of susceptible mosquitoes (constant recruitment)
  reaction(~ Lambda_m, c(Sm = +1), "mosq_birth"),
  
  # Infection of susceptible mosquitoes by infectious humans
  reaction(
    ~ (b * beta_hm * (Im + Is) / (S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_R * R + 
                                    S2 + E2 + (1 - prop_age1_Im) * Im + (1 - prop_age1_Is) * Is + 
                                    prop_age1_C * C + (1 - prop_age1_R) * R)) * Sm,
    c(Sm = -1, Em = +1),
    "Sm_to_Em"
  ),
  
  # Progression from exposed to infectious mosquitoes
  reaction(~ sigma_m * Em, c(Em = -1, Im_m = +1), "Em_to_Im_m"),
  
  # Death of mosquitoes
  reaction(~ mu_m * Sm,   c(Sm = -1),   "Sm_death"),
  reaction(~ mu_m * Em,   c(Em = -1),   "Em_death"),
  reaction(~ mu_m * Im_m, c(Im_m = -1), "Im_m_death")
)

# -------------------------------
# Compile reactions
# -------------------------------
compiled_reactions <- compile_reactions(
  reactions = reactions,
  state_ids = names(initial_state),
  params = params
)

# -------------------------------
# Run Gillespie SSA
# -------------------------------
set.seed(1000)
gillespie_output <- ssa(
  initial_state = initial_state,
  reactions = compiled_reactions,
  params = params,
  method = ssa_exact(),
  final_time = 1825,
  census_interval = 1,
  verbose = TRUE
)
# 
# # -------------------------------
# # Extract simulation results
# # -------------------------------
sim_df <- as.data.frame(gillespie_output$state)
sim_df$time <- gillespie_output$time

# Remove x. prefix from column names if present
colnames(sim_df) <- gsub("^x\\.", "", colnames(sim_df))

# -------------------------------
# Plot cumulative cases
# -------------------------------
ggplot(sim_df, aes(x = time, y = CumulativeCases)) +
  geom_line(color = "red", size = 1) +
  labs(
    x = "Time (days)",
    y = "Cumulative autochthonous cases",
    title = "Stochastic simulation of cumulative human infections"
  ) +
  theme_minimal(base_size = 14)
# 
# # Plot daily new infections
sim_df <- sim_df %>%
  mutate(
    daily_new_infections = c(0, diff(CumulativeCases))
  )

ggplot(sim_df, aes(x = time, y = daily_new_infections)) +
  geom_line(color = "blue", size = 1) +
  labs(
    x = "Time (days)",
    y = "Daily new infections",
    title = "Daily new human infections"
  ) +
  theme_minimal(base_size = 14)
# Check if total cases = cumulative infections
sim_df <- sim_df %>%
  mutate(
    total_ever_infected = E1 + E2 + Im + Is + C + R,
    model_cumulative = CumulativeCases
  )
# 
# # These should be roughly equal (allowing for deaths)
ggplot(sim_df) +
  geom_line(aes(x = time, y = total_ever_infected), color = "blue") +
  geom_line(aes(x = time, y = model_cumulative), color = "red", linetype = "dashed") +
  labs(x = "Time", y = "Count", title = "Verification: Total ever infected vs CumulativeCases")

# Number of simulations
n_sims <- 1000

# Prepare a list to store results
sim_results <- vector("list", n_sims)

set.seed(123)  # For reproducibility

for (i in 1:n_sims) {
  sim <- ssa(
    initial_state = initial_state,
    reactions = compiled_reactions,
    params = params,
    method = ssa_exact(),
    final_time = 90,
    census_interval = 1
  )
  
  df <- as.data.frame(sim$state)
  df$time <- sim$time
  colnames(df) <- gsub("^x\\.", "", colnames(df))
  
  # Compute cumulative cases from E1 + E2
  df <- df %>%
    mutate(
      new_infections = (E1 + E2) - lag(E1 + E2, default = 0),
      new_infections = ifelse(new_infections < 0, 0, new_infections),
      cumulative_cases = cumsum(new_infections)
    )
  
  sim_results[[i]] <- df$cumulative_cases
}

# Convert list to matrix: rows = time points, cols = simulations
sim_matrix <- do.call(cbind, sim_results)
time_points <- df$time

# Compute median, mean, and 95th percentile
summary_df <- data.frame(
  time = time_points,
  median = apply(sim_matrix, 1, median),
  mean = apply(sim_matrix, 1, mean),
  p95 = apply(sim_matrix, 1, quantile, probs = 0.975),
  p05 = apply(sim_matrix, 1, quantile, probs = 0.025)
)
ggplot(summary_df, aes(x = time)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = median), color = "red", size = 1) +
  geom_line(aes(y = mean), color = "blue", linetype = "dashed", size = 1) +
  labs(
    x = "Time (days)",
    y = "Cumulative cases",
    title = "Stochastic simulation of outbreak (median, mean, 95% interval)"
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(expand = c(0, 0))

#Plot the susceptible human population over time separately
ggplot(sim_df, aes(x = time)) +
  geom_line(aes(y = S1, color = "S1"), size = 1)+
  labs(
    x = "Time (days)",
    y = "Number of Susceptible Humans",
    title = "Susceptible Human Population Over Time"
  ) +
  scale_color_manual(values = c("S1" = "blue", "S2" = "green"), 
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

#Plot the human population dynamics(Exposed/ Infectious/recovered)
ggplot(sim_df, aes(x=time))+
  geom_line(aes(y=E1+E2, color="Exposed"), size=1)+
  geom_line(aes(y=Im+Is, color="Infectious"), size=1)+
  geom_line(aes(y=C, color="Chronic"), size=1)+
  geom_line(aes(y=R, color="Recovered"), size=1)+
  labs(
    x="Time(days)",
    y="Population",
    title="Human Population Dynamics")