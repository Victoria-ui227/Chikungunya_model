#================================================
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GillespieSSA2)

#================================================
# 1. TEMPERATURE INTERPOLATION (raw °C, daily)
# Outbreak start : 14 March 2025 = day 0
# Outbreak end   : 29 June  2025 = day 107
# Drop Jan & Feb — data starts from March
#================================================

mombasa_temp <- read.csv("mombasa_forecasted_temperature.csv")

# Assign each monthly value to the 15th of its month
mombasa_temp$date <- as.Date(paste(mombasa_temp$year,
                                   mombasa_temp$month,
                                   15, sep = "-"))
mombasa_temp      <- mombasa_temp[order(mombasa_temp$date), ]

# STEP 1: Drop January and February
mombasa_temp <- mombasa_temp[mombasa_temp$month >= 3, ]

# STEP 2: Expand each monthly value to daily rows (2025 = non-leap year)
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

daily_temp_df <- data.frame(
  date = as.Date(unlist(lapply(seq_len(nrow(mombasa_temp)), function(i) {
    m <- mombasa_temp$month[i]
    y <- mombasa_temp$year[i]
    seq(as.Date(paste(y, m, 1, sep = "-")),
        by         = "day",
        length.out = days_in_month[m])
  })), origin = "1970-01-01"),
  temp_celsius = rep(mombasa_temp$mean_lst_celsius,
                     times = days_in_month[mombasa_temp$month])
)

# STEP 3: Subset to simulation window
sim_start_date <- as.Date("2025-03-14")   # day 0
sim_end_date   <- as.Date("2025-06-29")   # day 107
projected_days <- as.integer(sim_end_date - sim_start_date) + 1L  # 108

daily_dates <- seq(sim_start_date, sim_end_date, by = "day")

temp_projected_daily <- daily_temp_df$temp_celsius[
  daily_temp_df$date %in% daily_dates
]

# Guard: confirm length matches simulation window
stopifnot(length(temp_projected_daily) == projected_days)
cat("Temperature days:", length(temp_projected_daily),
    "| Window:", format(sim_start_date), "to", format(sim_end_date), "\n")

# Fast integer-day lookup — day 0 → index 1, day 107 → index 108
temp_at <- function(day) {
  i <- floor(day) + 1L
  i <- max(1L, min(i, length(temp_projected_daily)))
  temp_projected_daily[i]
}

#================================================
# 2. TEMPERATURE-DEPENDENT MOSQUITO PARAMETERS
# Raw Celsius input — original polynomial fits
#================================================

# Extrinsic incubation completion rate
# Higher temp → faster virus maturation inside the mosquito
sigma_m_T <- function(T) (4 + exp(5.15 - 0.123 * T))

# Per-capita mosquito offspring rate
# Unimodal — peaks near 25–28°C, zero at thermal extremes
offspring_T <- function(T) {
  (-341.8585 + 108.7373*T - 12.3447*T^2 +
         0.6367451*T^3 - 0.0149469*T^4 + 0.0001294166*T^5)
}

# Mosquito daily mortality rate
# U-shaped — lowest near 25°C, rises at cold and hot extremes
mu_m_T <- function(T) {
  (0.8692 - 0.1590*T + 0.01116*T^2 -
         0.0003408*T^3 + 0.000003809*T^4)
}

#================================================
# 3. PARAMETERS & INITIAL CONDITIONS
#================================================

params <- list(
  mu             = 1 / (62 * 365),
  sigma          = 1 / 5,
  gamma_m        = 1 / 7,
  gamma_s        = 1 / 14,
  gamma_c        = 1 / 30,
  severe_prob_1  = 0.23,
  severe_prob_2  = 0.49,
  chronic_prob_s = 0.17,
  contact_rate_1 = 0.6,
  contact_rate_2 = 0.8,
  b              = 0.41,
  beta_hm        = 0.442,
  beta_mh        = 0.333,
  # Temperature-driven mosquito params —
  # initialised at day 0 temp, updated daily in Gillespie runner
  mu_m    = mu_m_T(temp_at(0)),
  sigma_m = sigma_m_T(temp_at(0)),
  Lambda_m = mu_m_T(temp_at(0)) * offspring_T(temp_at(0))
)

init <- c(
  S1      = 751956,
  E1      = 1,
  S2      = 46300,
  E2      = 0,
  Im      = 0,
  Is      = 0,
  C       = 0,
  R       = 0,
  Sm      = 1000,
  Em      = 250,
  Im_mosq = 80
)

#================================================
# 4. DETERMINISTIC MODEL (temperature-driven)
#================================================

chikungunya_det <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    # Temperature and derived mosquito params at time t
    Tnow     <- temp_at(t)
    sigma_m  <- sigma_m_T(Tnow)
    mu_m     <- mu_m_T(Tnow)
    Lambda_m <- mu_m * offspring_T(Tnow)
    
    # Population totals
    N       <- S1 + E1 + S2 + E2 + Im + Is + C + R
    I_total <- Im + Is
    M       <- Sm + Em + Im_mosq
    
    # Force of infection
    lambda_h <- b * beta_mh * Im_mosq / N
    lambda_m <- b * beta_hm * I_total / N
    
    # Human ODEs
    dS1 <- -mu*S1 - lambda_h*contact_rate_1*S1 + mu*N
    dE1 <-  lambda_h*contact_rate_1*S1 - sigma*E1 - mu*E1
    dS2 <- -mu*S2 - lambda_h*contact_rate_2*S2
    dE2 <-  lambda_h*contact_rate_2*S2 - sigma*E2 - mu*E2
    dIm <-  sigma*(1 - severe_prob_1)*E1 + sigma*(1 - severe_prob_2)*E2 - gamma_m*Im - mu*Im
    dIs <-  sigma*severe_prob_1*E1 + sigma*severe_prob_2*E2 - gamma_s*Is - mu*Is - chronic_prob_s*gamma_s*Is
    dC  <-  chronic_prob_s*gamma_s*Is - gamma_c*C - mu*C
    dR  <-  gamma_m*Im + (1 - chronic_prob_s)*gamma_s*Is + gamma_c*C - mu*R
    
    # Mosquito ODEs
    dSm      <- Lambda_m*M - mu_m*Sm - lambda_m*Sm
    dEm      <- lambda_m*Sm - sigma_m*Em - mu_m*Em
    dIm_mosq <- sigma_m*Em - mu_m*Im_mosq
    
    list(c(dS1, dE1, dS2, dE2, dIm, dIs, dC, dR, dSm, dEm, dIm_mosq))
  })
}

# Solve from day 0 to day 107 in 0.2-day steps
time_points     <- seq(0, projected_days - 1, by = 0.2)   # 0 to 107
det_out         <- as.data.frame(ode(init, time_points, chikungunya_det, params))
det_out$I_total <- det_out$Im + det_out$Is

#================================================
# 5. GILLESPIE STOCHASTIC MODEL (temperature-driven)
# Uses GillespieSSA2 — re-compiles params each day
# to inject updated mu_m, sigma_m, Lambda_m
#================================================
params <- c(
  mu             = 1 / (62 * 365),
  sigma          = 1 / 5,
  gamma_m        = 1 / 7,
  gamma_s        = 1 / 14,
  gamma_c        = 1 / 30,
  severe_prob_1  = 0.23,
  severe_prob_2  = 0.49,
  chronic_prob_s = 0.17,
  contact_rate_1 = 1.0,
  contact_rate_2 = 0.8,
  b              = 0.41,
  beta_hm        = 0.442,
  beta_mh        = 0.333,
  # Temperature-driven mosquito params — initialised at day 0
  # overwritten daily inside run_ssa_temperature()
  mu_m     = mu_m_T(temp_at(0)),
  sigma_m  = sigma_m_T(temp_at(0)),
  Lambda_m = mu_m_T(temp_at(0)) * offspring_T(temp_at(0))
)

# Initial state must match init exactly, using Im_m for mosquito infectious
initial_state <- c(
  S1    = 751956,
  E1    = 1,
  S2    = 46300,
  E2    = 0,
  Im    = 0,
  Is    = 0,
  C     = 0,
  R     = 0,
  Sm    = 5000,
  Em    = 2500,
  Im_m  = 80       # Im_mosq renamed Im_m for GillespieSSA2
)

#-------------------------------
# REACTIONS
#-------------------------------
reactions <- list(
  
  ## Human births
  reaction(~ mu * S1, c(S1 = +1), "birth_S1"),
  reaction(~ mu * S2, c(S2 = +1), "birth_S2"),
  
  ## S1 → E1
  reaction(
    ~ b * beta_mh * (Im_m / (S1+E1+S2+E2+Im+Is+C+R)) * contact_rate_1 * S1,
    c(S1 = -1, E1 = +1), "S1_to_E1"
  ),
  
  ## S2 → E2
  reaction(
    ~ b * beta_mh * (Im_m / (S1+E1+S2+E2+Im+Is+C+R)) * contact_rate_2 * S2,
    c(S2 = -1, E2 = +1), "S2_to_E2"
  ),
  
  ## E1 → Im / Is
  reaction(~ sigma*(1 - severe_prob_1)*E1, c(E1 = -1, Im = +1), "E1_to_Im"),
  reaction(~ sigma*severe_prob_1*E1,       c(E1 = -1, Is = +1), "E1_to_Is"),
  
  ## E2 → Im / Is
  reaction(~ sigma*(1 - severe_prob_2)*E2, c(E2 = -1, Im = +1), "E2_to_Im"),
  reaction(~ sigma*severe_prob_2*E2,       c(E2 = -1, Is = +1), "E2_to_Is"),
  
  ## Im → R
  reaction(~ gamma_m*Im, c(Im = -1, R = +1), "Im_to_R"),
  
  ## Is → R / C
  reaction(~ gamma_s*(1 - chronic_prob_s)*Is, c(Is = -1, R = +1), "Is_to_R"),
  reaction(~ gamma_s*chronic_prob_s*Is,       c(Is = -1, C = +1), "Is_to_C"),
  
  ## C → R
  reaction(~ gamma_c*C, c(C = -1, R = +1), "C_to_R"),
  
  ## Human deaths
  reaction(~ mu*S1, c(S1 = -1), "S1_death"),
  reaction(~ mu*S2, c(S2 = -1), "S2_death"),
  reaction(~ mu*E1, c(E1 = -1), "E1_death"),
  reaction(~ mu*E2, c(E2 = -1), "E2_death"),
  reaction(~ mu*Im, c(Im = -1), "Im_death"),
  reaction(~ mu*Is, c(Is = -1), "Is_death"),
  reaction(~ mu*C,  c(C  = -1), "C_death"),
  reaction(~ mu*R,  c(R  = -1), "R_death"),
  
  ## Mosquito reactions — Lambda_m, mu_m, sigma_m updated daily
  reaction(~ Lambda_m*(Sm + Em + Im_m),                                         c(Sm = +1),          "mosq_birth"),
  reaction(~ b*beta_hm*((Im+Is)/(S1+E1+S2+E2+Im+Is+C+R))*Sm, c(Sm = -1, Em = +1), "Sm_to_Em"),
  reaction(~ sigma_m*Em,  c(Em  = -1, Im_m = +1), "Em_to_Im_m"),
  reaction(~ mu_m*Sm,     c(Sm  = -1),             "Sm_death"),
  reaction(~ mu_m*Em,     c(Em  = -1),             "Em_death"),
  reaction(~ mu_m*Im_m,   c(Im_m = -1),            "Im_m_death")
)

#-------------------------------
# TEMPERATURE-DRIVEN SSA RUNNER
# Steps one day at a time,
# updating mosquito params each day
#-------------------------------
#-------------------------------
# COMPILE ONCE — outside the loop
# GillespieSSA2 uses params by
# reference in the numeric vector
# so updating params["mu_m"] etc.
# before each ssa() call is enough
#-------------------------------
compiled_reactions <- compile_reactions(
  reactions  = reactions,
  state_ids  = names(initial_state),
  params     = params
)

run_ssa_temperature <- function(initial_state, params, compiled,
                                t_end = 107, n_runs = 50) {
  
  all_runs <- vector("list", n_runs)
  n_states <- length(initial_state)
  
  for (run_id in seq_len(n_runs)) {
    
    set.seed(run_id)
    current_state <- initial_state
    
    out           <- matrix(NA_real_, nrow = t_end + 1L,
                            ncol = n_states + 1L)
    colnames(out) <- c("time", names(initial_state))
    out[1, ]      <- c(0, current_state)
    
    for (day in seq_len(t_end)) {
      
      Tnow <- temp_at(day - 1L)
      
      # Update the three temperature-driven values in the numeric vector
      # No recompile needed — ssa() reads params fresh each call
      params["mu_m"]     <- mu_m_T(Tnow)
      params["sigma_m"]  <- sigma_m_T(Tnow)
      params["Lambda_m"] <- params["mu_m"] * offspring_T(Tnow)
      
      result <- ssa(
        initial_state   = current_state,
        reactions       = compiled,        # same compiled object every day
        params          = params,          # updated numeric vector
        final_time      = 1,
        method          = ssa_etl(tau = 0.01),
        log_propensity  = FALSE,
        log_firings     = FALSE,
        log_buffer      = FALSE,
        census_interval = 1,
        verbose         = FALSE
      )
      
      final_row     <- result$state[nrow(result$state), ]
      current_state <- final_row[names(initial_state)]
      current_state <- pmax(current_state, 0)
      
      out[day + 1L, ] <- c(day, current_state)
    }
    
    out_df         <- as.data.frame(out)
    out_df$I_total <- out_df$Im + out_df$Is
    out_df$run_id  <- run_id
    all_runs[[run_id]] <- out_df
    
    gc()
    cat("Run", run_id, "of", n_runs, "complete\n")
  }
  
  dplyr::bind_rows(all_runs)
}

stoch_out <- run_ssa_temperature(
  initial_state = initial_state,
  params        = params,
  compiled      = compiled_reactions,
  t_end         = projected_days - 1L,
  n_runs        = 50
)

#================================================
# 6. SUMMARY & PLOTS
#================================================

# Stochastic summary across 50 runs
stoch_summary <- stoch_out %>%
  group_by(time) %>%
  summarise(
    med = median(I_total),
    lo  = quantile(I_total, 0.025),
    hi  = quantile(I_total, 0.975),
    .groups = "drop"
  )

# Add calendar dates to both outputs
det_out$date   <- sim_start_date + det_out$time
stoch_summary$date <- sim_start_date + stoch_summary$time
stoch_out$date     <- sim_start_date + stoch_out$time

# Plot — total infectious humans
ggplot() +
  geom_line(data = stoch_out,
            aes(x = date, y = I_total, group = run_id),
            colour = "#7fbbe3", alpha = 0.2, linewidth = 0.3) +
  geom_ribbon(data = stoch_summary,
              aes(x = date, ymin = lo, ymax = hi),
              fill = "#2196F3", alpha = 0.25) +
  geom_line(data = stoch_summary,
            aes(x = date, y = med),
            colour = "#1565C0", linewidth = 1) +
  geom_line(data = det_out,
            aes(x = date, y = I_total),
            colour = "#E53935", linewidth = 1, linetype = "dashed") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  labs(
    title    = "Chikungunya — total infectious humans",
    subtitle = "Blue: 50 stochastic runs (median + 95% CI)  |  Red dashed: deterministic ODE",
    x = "Date", y = "Infectious individuals (Im + Is)"
  ) +
  theme_minimal(base_size = 13)

# Plot — infectious mosquitoes
mosq_summary <- stoch_out %>%
  group_by(time) %>%
  summarise(
    med = median(Im_m),
    lo  = quantile(Im_m, 0.025),
    hi  = quantile(Im_m, 0.975),
    .groups = "drop"
  )
mosq_summary$date <- sim_start_date + mosq_summary$time

ggplot() +
  geom_line(data = stoch_out,
            aes(x = date, y = Im_m, group = run_id),
            colour = "#a5d6a7", alpha = 0.2, linewidth = 0.3) +
  geom_ribbon(data = mosq_summary,
              aes(x = date, ymin = lo, ymax = hi),
              fill = "#4CAF50", alpha = 0.25) +
  geom_line(data = mosq_summary,
            aes(x = date, y = med),
            colour = "#2E7D32", linewidth = 1) +
  geom_line(data = det_out,
            aes(x = date, y = Im_mosq),
            colour = "#E53935", linewidth = 1, linetype = "dashed") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  labs(
    title    = "Chikungunya — infectious mosquitoes",
    subtitle = "Green: 50 stochastic runs (median + 95% CI)  |  Red dashed: deterministic ODE",
    x = "Date", y = "Infectious mosquitoes"
  ) +
  theme_minimal(base_size = 13)