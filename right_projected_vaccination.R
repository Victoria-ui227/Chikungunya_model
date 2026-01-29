###############################
# LIBRARIES
###############################
pacman::p_load(
  GillespieSSA2,
  dplyr,
  ggplot2
)

###############################
# PROJECTED TEMPERATURE DATA
###############################
projected_days <- 3650
days_per_month <- 30

mombasa_projected_temp <- read.csv("mombasa_forecasted_temperature.csv")

# Scale projected temperature (0–1)
temp_projected_scaled <-
  (mombasa_projected_temp$mean_lst_celsius -
     min(mombasa_projected_temp$mean_lst_celsius)) /
  (max(mombasa_projected_temp$mean_lst_celsius) -
     min(mombasa_projected_temp$mean_lst_celsius))

temp_projected_daily <-
  rep(temp_projected_scaled, each = days_per_month)[1:projected_days]
###############################
# CLIMATE-DRIVEN PARAMETERS
###############################
Lambda_m0 <- 20000
mu_m0     <- 1/10
sigma_m0  <- 1/5

# How projected temperature affects mosquito population
projected_Lambda_m_daily <- Lambda_m0 * (0.3 + temp_projected_daily)
projected_mu_m_daily     <- mu_m0 * (1.2 - temp_projected_daily)
projected_sigma_m_daily  <- sigma_m0 * (0.5 + temp_projected_daily)

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
  V1=0,
  V2=0,
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
  
  prop_age1_Im = 0.765,      
  prop_age1_Is = 0.003,      
  prop_age1_C = 0.015,       
  prop_age1_R = 0.85,
  prop_age1_V=0.8,
  
  vaccination_rate_1=0.03,
  vaccination_rate_2=0.01,
  waning_rate = 1/365,
  
  contact_rate_1 = 0.7,
  contact_rate_2 = 0.3,
  
  amp = 0.5,
  phi = 0,
  pi=3.141592653589793,
  
  mu_m = mu_m0,
  sigma_m = sigma_m0,
  
  b = 0.6,
  beta_hm = 0.442,
  beta_mh = 0.333,
  
  # alpha = 1/(59*365),
  Lambda_m = Lambda_m0,
  
  imported_rate = 1
)


###############################
# REACTIONS
###############################
reactions <- list(
  
  ## ---------------- HUMAN Reactions ----------------
  
  # Human births
  reaction(
    ~ mu * (S1 + E1 + prop_age1_Im*Im + prop_age1_Is*Is + prop_age1_C*C + prop_age1_R*R),
    c(S1 = +1),
    "birth_S1"
  ),
  
  # Infection S1 -> E1
  reaction(
    ~ (1 + amp*cos(2*3.141592653589793*(time/(365*5)-phi))) *
      (b * beta_mh * Im_m / (S1 + E1 + S2 + E2 + Im + Is + R + C)) * contact_rate_1 * S1,
    c(S1 = -1, E1 = +1, CumulativeCases = +1),
    "S1_to_E1"
  ),
  
  # Imported cases S1
  reaction(
    ~ imported_rate * 0.7 * S1 / (S1 + S2),
    c(S1 = -1, E1 = +1, ImportedCases = +1, CumulativeCases = +1),
    "imported_S1"
  ),
  
  # E1 -> Im / Is
  reaction(~ sigma * (1 - severe_prob_1) * E1, c(E1 = -1, Im = +1), "E1_to_Im"),
  reaction(~ sigma * severe_prob_1 * E1, c(E1 = -1, Is = +1), "E1_to_Is"),
  
  # Aging S1 -> S2
  # reaction(~ alpha * S1, c(S1 = -1, S2 = +1), "aging"),
  
  # S2 -> E2
  reaction(
    ~ (1 + amp*cos(2*3.141592653589793*(time/(365*5)-phi))) *
      (b * beta_mh * Im_m / (S1 + E1 + S2 + E2 + Im + Is + R + C)) * contact_rate_2 * S2,
    c(S2 = -1, E2 = +1, CumulativeCases = +1),
    "S2_to_E2"
  ),
  
  # Imported cases S2
  reaction(
    ~ imported_rate * 0.3 * S2 / (S1 + S2),
    c(S2 = -1, E2 = +1, ImportedCases = +1, CumulativeCases = +1),
    "imported_S2"
  ),
  
  # E2 -> Im / Is
  reaction(~ sigma * (1 - severe_prob_2) * E2, c(E2 = -1, Im = +1), "E2_to_Im"),
  reaction(~ sigma * severe_prob_2 * E2, c(E2 = -1, Is = +1), "E2_to_Is"),
  
  # Im -> Recovered
  reaction(~ gamma_m * Im, c(Im = -1, R = +1), "Im_to_R"),
  
  # Is -> Recovered or chronic
  reaction(~ gamma_s * Is * (1 - chronic_prob_s), c(Is = -1, R = +1), "Is_to_R"),
  reaction(~ gamma_s * Is * chronic_prob_s, c(Is = -1, C = +1), "Is_to_C"),
  
  # Chronic -> Recovered
  reaction(~ gamma_c * C, c(C = -1, R = +1), "C_to_R"),
  
  # Human deaths
  reaction(~ mu * S1, c(S1 = -1), "S1_death"),
  reaction(~ mu * S2, c(S2 = -1), "S2_death"),
  reaction(~ mu * E1, c(E1 = -1), "E1_death"),
  reaction(~ mu * E2, c(E2 = -1), "E2_death"),
  reaction(~ mu * Im, c(Im = -1), "Im_death"),
  reaction(~ mu * Is, c(Is = -1), "Is_death"),
  reaction(~ mu * C, c(C = -1), "C_death"),
  reaction(~ mu * R, c(R = -1), "R_death"),
  
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
# SIMULATION 1: BASELINE
###############################
set.seed(123)
state <- initial_state
out1 <- vector("list", projected_days)

for (d in 1:projected_days) {
  
  # Update daily mosquito parameters
  params["Lambda_m"] <- projected_Lambda_m_daily[d]
  params["mu_m"]     <- projected_mu_m_daily[d]
  params["sigma_m"]  <- projected_sigma_m_daily[d]
  
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
  
  out1[[d]] <- data.frame(time = d, t(state))
}

baseline_df <- bind_rows(out1)
baseline_df$infectious <- baseline_df$Im + baseline_df$Is
baseline_df$type <- "Baseline"

###############################
# VACCINATION
###############################

####Reactions$###############################
vac_reactions <- list(
  
  ## ---------------- HUMAN Reactions ----------------
  
  # Human births
  reaction(
    ~ mu * (S1 + E1 + S2 + E2 + Im + Is + C + R + V1 + V2),
    c(S1 = +1),
    "birth_S1"
  ),
  
  # Infection S1 -> E1
  reaction(
    ~ (1 + amp*cos(2*pi*(time/(365*5)-phi))) *
      (b * beta_mh * Im_m /
         (S1 + E1 + S2 + E2 + Im + Is + R + C + V1 + V2)) *
      contact_rate_1 * S1,
    c(S1 = -1, E1 = +1, CumulativeCases = +1),
    "S1_to_E1"
  ),
  
  # Infection S2 -> E2
  reaction(
    ~ (1 + amp*cos(2*pi*(time/(365*5)-phi))) *
      (b * beta_mh * Im_m /
         (S1 + E1 + S2 + E2 + Im + Is + R + C + V1 + V2)) *
      contact_rate_2 * S2,
    c(S2 = -1, E2 = +1, CumulativeCases = +1),
    "S2_to_E2"
  ),
  
  # E → I
  reaction(~ sigma * (1 - severe_prob_1) * E1, c(E1 = -1, Im = +1)),
  reaction(~ sigma * severe_prob_1 * E1, c(E1 = -1, Is = +1)),
  reaction(~ sigma * (1 - severe_prob_2) * E2, c(E2 = -1, Im = +1)),
  reaction(~ sigma * severe_prob_2 * E2, c(E2 = -1, Is = +1)),
  
  # Recovery
  reaction(~ gamma_m * Im, c(Im = -1, R = +1)),
  reaction(~ gamma_s * Is * (1 - chronic_prob_s), c(Is = -1, R = +1)),
  reaction(~ gamma_s * Is * chronic_prob_s, c(Is = -1, C = +1)),
  reaction(~ gamma_c * C, c(C = -1, R = +1)),
  
  # Waning immunity
  reaction(~ waning_rate * V1, c(V1 = -1, S1 = +1)),
  reaction(~ waning_rate * V2, c(V2 = -1, S2 = +1)),
  
  # Deaths
  reaction(~ mu * S1, c(S1 = -1)),
  reaction(~ mu * S2, c(S2 = -1)),
  reaction(~ mu * E1, c(E1 = -1)),
  reaction(~ mu * E2, c(E2 = -1)),
  reaction(~ mu * Im, c(Im = -1)),
  reaction(~ mu * Is, c(Is = -1)),
  reaction(~ mu * C, c(C = -1)),
  reaction(~ mu * R, c(R = -1)),
  reaction(~ mu * V1, c(V1 = -1)),
  reaction(~ mu * V2, c(V2 = -1)),
  
  ## ---------------- MOSQUITOES ----------------
  reaction(~ Lambda_m, c(Sm = +1)),
  reaction(
    ~ (b * beta_hm * (Im + Is) /
         (S1 + E1 + S2 + E2 + Im + Is + R + C + V1 + V2)) * Sm,
    c(Sm = -1, Em = +1)
  ),
  reaction(~ sigma_m * Em, c(Em = -1, Im_m = +1)),
  reaction(~ mu_m * Sm, c(Sm = -1)),
  reaction(~ mu_m * Em, c(Em = -1)),
  reaction(~ mu_m * Im_m, c(Im_m = -1))
)

###############################
# COMPILE
###############################
vac_compiled_reactions <- compile_reactions(
  reactions = vac_reactions,
  state_ids = names(initial_state),
  params = params
)

###############################
# SIMULATION 2: VACCINATION PULSE
###############################
set.seed(456)
state <- initial_state
out2 <- vector("list", projected_days)

for (d in 1:projected_days) {
  
  # Update daily mosquito parameters
  params["Lambda_m"] <- projected_Lambda_m_daily[d]
  params["mu_m"]     <- projected_mu_m_daily[d]
  params["sigma_m"]  <- projected_sigma_m_daily[d]
  
  # Vaccination pulse on day 100 of year 1
  if (d == 100) {
    n_vax_S1 <- floor(0.05 * state["S1"])
    n_vax_S2 <- floor(0.05 * state["S2"])
    
    state["S1"] <- state["S1"] - n_vax_S1
    state["S2"] <- state["S2"] - n_vax_S2
    state["V1"] <- state["V1"] + n_vax_S1
    state["V2"] <- state["V2"] + n_vax_S2
  }
  
  sim <- ssa(
    initial_state = state,
    reactions = vac_compiled_reactions,
    params = params,
    method = ssa_exact(),
    final_time = 1,
    census_interval = 1
  )
  
  state <- as.numeric(tail(sim$state, 1))
  names(state) <- colnames(sim$state)
  
  out2[[d]] <- data.frame(time = d, t(state))
}

vax_df <- bind_rows(out2)
vax_df$infectious <- vax_df$Im + vax_df$Is
vax_df$type <- "Vaccination"

###############################
#vaccination at 10% each
#####################

set.seed(789)
state <- initial_state
out3 <- vector("list", projected_days)

for (d in 1:projected_days) {
  
  # Update daily mosquito parameters
  params["Lambda_m"] <- projected_Lambda_m_daily[d]
  params["mu_m"]     <- projected_mu_m_daily[d]
  params["sigma_m"]  <- projected_sigma_m_daily[d]
  
  # Vaccination pulse on day 100 of year 1
  if (d == 100) {
    n_vax_S1 <- floor(0.1 * state["S1"])
    n_vax_S2 <- floor(0.1 * state["S2"])
    
    state["S1"] <- state["S1"] - n_vax_S1
    state["S2"] <- state["S2"] - n_vax_S2
    state["V1"] <- state["V1"] + n_vax_S1
    state["V2"] <- state["V2"] + n_vax_S2
  }
  
  sim <- ssa(
    initial_state = state,
    reactions = vac_compiled_reactions,
    params = params,
    method = ssa_exact(),
    final_time = 1,
    census_interval = 1
  )
  
  state <- as.numeric(tail(sim$state, 1))
  names(state) <- colnames(sim$state)
  
  out3[[d]] <- data.frame(time = d, t(state))
}

ten_percent_vax_df <- bind_rows(out3)
ten_percent_vax_df$infectious_ten_percent <- ten_percent_vax_df$Im + ten_percent_vax_df$Is
ten_percent_vax_df$type <- "Vaccination"


############################
#Vaccinating 20%
###########################
set.seed(789)
state <- initial_state
out4 <- vector("list", projected_days)

for (d in 1:projected_days) {
  
  # Update daily mosquito parameters
  params["Lambda_m"] <- projected_Lambda_m_daily[d]
  params["mu_m"]     <- projected_mu_m_daily[d]
  params["sigma_m"]  <- projected_sigma_m_daily[d]
  
  # Vaccination pulse on day 100 of year 1
  if (d == 100) {
    n_vax_S1 <- floor(0.2 * state["S1"])
    n_vax_S2 <- floor(0.2 * state["S2"])
    
    state["S1"] <- state["S1"] - n_vax_S1
    state["S2"] <- state["S2"] - n_vax_S2
    state["V1"] <- state["V1"] + n_vax_S1
    state["V2"] <- state["V2"] + n_vax_S2
  }
  
  sim <- ssa(
    initial_state = state,
    reactions = vac_compiled_reactions,
    params = params,
    method = ssa_exact(),
    final_time = 1,
    census_interval = 1
  )
  
  state <- as.numeric(tail(sim$state, 1))
  names(state) <- colnames(sim$state)
  
  out4[[d]] <- data.frame(time = d, t(state))
}

twenty_percent_vax_df <- bind_rows(out4)
twenty_percent_vax_df$infectious_twenty_percent <- twenty_percent_vax_df$Im + twenty_percent_vax_df$Is
twenty_percent_vax_df$type <- "Vaccination"

combined_df <- combined_df %>%
  left_join(
    twenty_percent_vax_df %>% select(time, infectious_twenty_percent),
    by = "time"
  )

# Combine all scenarios into one dataframe
combined_df <- baseline_df %>%
  select(time, infectious) %>%
  rename(infectious_baseline = infectious) %>%
  left_join(
    vax_df %>% select(time, infectious) %>% rename(infectious_vax = infectious),
    by = "time"
  ) %>%
  left_join(
    ten_percent_vax_df %>% select(time, infectious_ten_percent),
    by = "time"
  ) %>%
  left_join(
    twenty_percent_vax_df %>% select(time, infectious_twenty_percent),
    by = "time"
  )

# Add actual dates
combined_df$date <- as.Date("2025-01-01") + (combined_df$time - 1)

ggplot(combined_df, aes(x = date)) +
  geom_line(aes(y = infectious_baseline, color = "Baseline"), size = 1) +
  geom_line(aes(y = infectious_vax, color = "Vaccination 5%"), size = 1) +
  geom_line(aes(y = infectious_ten_percent, color = "Vaccination 10%"), size = 1) +
  geom_line(aes(y = infectious_twenty_percent, color = "Vaccination 20%"), size = 1) +
  labs(
    x = "Date",
    y = "Infectious Humans",
    title = "Baseline vs Vaccination (Projected)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Baseline" = "black",
    "Vaccination 5%" = "blue",
    "Vaccination 10%" = "green",
    "Vaccination 20%" = "red"
  ))

