

# Load Libraries ----------------------------------------------------------

pacman:: p_load(
  tidyverse,
  deSolve,
  plotly,
  gridExtra,
  reshape2,
  dplyr
)



# define equations --------------------------------------------------------

#D:Define the equations
#I:Initialize the compartments
#P: Define parameters
#T: Define time interval
#O : solve the equations 


chikungunya_model <- function(
    t,y,parms){
  with(as.list(c(parms,y)),{
    
    #Age-specific compartments
    S1 <- y[1]; E1 <- y[2]  #Age Group 1 (<60 years): Susceptible, Exposed
    S2 <- y[3]; E2 <- y[4]  #Age group 2 (60+ years): Susceptible, Exposed
    
    #Age Independent compartments
    Im <- y[5] #Infectious (mild)
    Is <- y[6] #Infectious(severe)
    C <-  y[7] #Chronic
    R <-  y[8] #Recovered
    #V <- y[9] #Vaccinated
    
    #Mosquito compartments
    Sm <- y[9]; Em <- y[10]; Im_mosq <- y[11]
    
    #Total human populations by age group 
    N1 = S1 + E1 + prop_age1_Im * Im + prop_age1_Is * Is + prop_age1_C* C + prop_age1_R *R 
    N2 = S2 + E2 + (1- prop_age1_Im)*Im + (1- prop_age1_Is)*Is + (1- prop_age1_C)* C + (1- prop_age1_R) *R
    N_total <- N1 + N2
    
    #Total infectious humans 
    I_total <- Im +Is
    
    # Force of infection from mosquitoes to humans
    lambda_h <- b * beta_mh * Im_mosq / N_total
    
    # Force of infection from humans to mosquitoes
    lambda_m <- b * beta_hm * I_total/ N_total
    
    # Age-specific forces of infection
    lambda_h1 <- lambda_h * contact_rate_1
    lambda_h2 <- lambda_h * contact_rate_2
    
    #Age group 1 (<60 years)
    dS1 = mu *N1 - mu * S1 - lambda_h1 *S1 #mu; birth/death rate
    dE1 = lambda_h1* S1 - (sigma+mu)* E1
    
    #Age group 2 (60+ years)
    dS2 = mu * N2 - mu* S2 - lambda_h2 *S2
    dE2 = lambda_h2 *S2 - (sigma + mu)*E2
    
    #Infectious mild~ receives from both age groups based on severity probabilities
    dIm = (sigma * (1- severe_prob_1)* E1) + (sigma * (1 - severe_prob_2)* E2) - (mu + gamma_m + mu_ms)*Im
    
    #Infectious severe~ receives from both age groups based on severity probabilities
    dIs = (sigma * severe_prob_1 * E1)+ (sigma * severe_prob_2 * E2)+ (mu_ms * Im)- (mu + gamma_s + chronic_prob_s* gamma_s)* Is
    
    #Chronic 
    dC = (chronic_prob_s * gamma_s * Is)- (mu + gamma_c)* C
    
    #Recovered - receives from mild, severe, and chronic compartments
    dR = (gamma_m *Im)+ ((1- chronic_prob_s)* gamma_s *Is)+ (gamma_c * C)- (mu * R)
    
    #Mosquito equations
    dSm = lambda_m - mu_m * Sm - lambda_m * Sm
    dEm = lambda_m * Sm - (sigma_m + mu_m)* Em
    dIm_mosq = sigma_m * Em - mu_m * Im_mosq
    
    # Return 11 derivatives for the ode solver
    list(c(dS1, dE1,dS2,dE2,dIm, dIs, dC, dR, dSm,dEm, dIm_mosq))
    
  })
}


# function to run simulation ----------------------------------------------

run_chikungunya_simulation <- function (params, initial_conditions, time_points){
  
  #Solve ODEs
  out<- ode(y= initial_conditions,
            times= time_points,
            func = chikungunya_model,
            parms = params,
            method = "lsoda")
  
  # Convert to dataframe
  out_df <- as.data.frame(out)
  
  #Add column names
  colnames(out_df) <- c ("time", "S1", "E1", "S2", "E2", "Im", "Is", "C", "R", "Sm", "Em", "Im_mosq")
  
    return(out_df)
}


# Function to plot results ------------------------------------------------

plot_chikungunya_results <- function(results){
  
  #Calculate derived quantities
  results$Total_Infectious <- results$Im + results$Is
  results$Total_Susceptible <- results$S1 + results$S2
  results$Total_Exposed <- results$E1 + results$E2
  results$Total_Mosquito <- results$Sm + results$Em + results$Im_mosq
  
  # Plot 1: Age-specific Susceptible and exposed
  age_specific_data <- results |>
    select(time, S1, E1,S2, E2) |>
    melt(id.vars = "time", variable.name = "Compartment", value.name = "Population")
  
  age_specific_data$Age_Group <- ifelse(grepl("1", age_specific_data$Compartment),
                                        "Age group 1 (<60)", "Age group 2 (60+)")
  age_specific_data$State <- ifelse(grepl("S", age_specific_data$Compartment),
                                    "Susceptible", "Exposed")
  
  p1 <- ggplot(age_specific_data, aes(x= time, y= Population, color = State))+
    geom_line(linewidth = 1.2)+
    facet_wrap(~Age_Group, scales ="free_y")+
    labs(title = "Age-Specific Susceptible and Exposed Populations",
         x= "Time(days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  
  #Plot 2: Age-independent disease states
  disease_data <- results |>
    select(time, Im, Is, C, R)|>
    melt(id.vars ="time",variable.name = "Compartment", value.name = "Population")
  
  disease_data$Compartment <- factor(disease_data$Compartment,
                                     levels = c("Im", "Is", "C","R"),
                                     labels = c("Infectious(mild)", "Infectious(severe)",
                                                "Chronic", "Recovered"))
  
  p2 <- ggplot(disease_data, aes(x= time, y = Population, color = Compartment))+
    geom_line(linewidth= 1.2)+
    labs(title = "Disease States (Age-Independent)",
         x = "Time (days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Plot 3: Mosquito population
  mosquito_data <- results |>
    select(time, Sm, Em, Im_mosq)|>
    melt(id.vars = "time", variable.name = "Compartment", value.name = "Population")
  
  mosquito_data$Compartment <- factor(mosquito_data$Compartment,
                                      levels = c("Sm","Em", "Im_mosq"),
                                      labels = c("Susceptible", "Exposed", "Infectious"))
  
  p3 <- ggplot(mosquito_data, aes(x= time, y = Population, color = Compartment)) +
    geom_line(linewidth = 1.2)+
    labs(title = "Mosquito Population",
         x= "Time(days)", y = "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Plot 4: Key Epidemiological Indicators
  epi_data <- results |>
    select(time, Total_Infectious, C, Total_Susceptible)|>
    melt(id.vars ='time', variable.name = "Indicator", value.name = "Population")
  
  
  epi_data$Indicator <- factor(epi_data$Indicator,
                               levels = c("Total_Infectious", "C", "Total_Susceptible"),
                               labels = c("Total Infectious", "Chronic", "Total Susceptible"))
  
  p4 <- ggplot(epi_data, aes(x= time, y = Population, color = Indicator))+
    geom_line(linewidth= 1.2)+
    labs(title = "Key Epidemiological Indicators",
         x= "Time(days)", y= "Population")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  #Arrange plots 
  grid.arrange(p1,p2,p3,p4,ncol=2)
  
} 

# Function to calculate RO ------------------------------------------------

calculate_RO <- function(params){
  with(params,{
    # Average infectious period considerinng chronic progression
    avg_infectious_period_mild <- 1 /(gamma_m + mu)
    avg_infectious_period_severe <- 1 /(gamma_s + mu)
    
    # Weighted average based on severity probabilities
    avg_severe_prob <- (severe_prob_1 + severe_prob_2)/2
    avg_infectious_period <- (1 - avg_severe_prob)* avg_infectious_period_mild + 
      avg_severe_prob * avg_infectious_period_severe
    
    avg_infectious_period_mosquito <- 1 / mu_m
    prob_transmission <- (sigma_m / ( sigma_m + mu_m))
    
    RO <- (b^2 * beta_hm * beta_mh * prob_transmission * avg_infectious_period * avg_infectious_period_mosquito)/mu_m
    
    return(RO)
  })
}

# Function to calculate cumulative incidence by age group -----------------

calculate_cumulative_incidence <- function(results, params) {
  
  # Calculate new infections over time (approximation)
  results$new_infections_1 <- c(0, diff(results$E1)) * params$sigma
  results$new_infections_2 <- c(0, diff(results$E2)) * params$sigma
  
  # Cumulative incidence
  results$cum_incidence_1 <- cumsum(pmax(0, results$new_infections_1))
  results$cum_incidence_2 <- cumsum(pmax(0, results$new_infections_2))
  
  return(results)
}


# Main execution function -------------------------------------------------

main_simulation <- function() {
  
  # Define parameters
  params <- list(
    # Human parameters
    mu = 1/(62*365),              # Human natural death rate
    sigma = 1/5,                  # 1/incubation period (days)
    gamma_m = 1/7,                # Recovery rate from mild infection
    gamma_s = 1/14,               # Recovery rate from severe infection
    gamma_c = 1/60,               # Recovery rate from chronic state (slower)
    severe_prob_1 = 0.1,          # Probability of severe disease in age group 1
    severe_prob_2 = 0.3,          # Probability of severe disease in age group 2
    #chronic_prob_m = 0.15,        # Probability of chronic from mild infection
    chronic_prob_s = 0.25,        # Probability of chronic from severe infection
    mu_ms = 1/10,          # Progression rate from Im to Is (1/10 days)
    #vaccination_rate_1 = 0.001,   # Vaccination rate for age group 1
    #vaccination_rate_2 = 0.002,   # Vaccination rate for age group 2
    contact_rate_1 = 1.0,         # Relative contact rate for age group 1
    contact_rate_2 = 0.8,         # Relative contact rate for age group 2
    
    # Proportions for age distribution in shared compartments (for demographic calculation)
    prop_age1_Im = 0.8,           # Proportion of mild infectious in age group 1
    prop_age1_Is = 0.7,           # Proportion of severe infectious in age group 1  
    prop_age1_C = 0.75,           # Proportion of chronic in age group 1
    prop_age1_R = 0.8,            # Proportion of recovered in age group 1
    #prop_age1_V = 0.8,            # Proportion of vaccinated in age group 1
    
    # Mosquito parameters
    Lambda_m = 1000,              # Mosquito recruitment rate
    mu_m = 1/14,                  # Mosquito death rate
    sigma_m = 1/7,                # 1/extrinsic incubation period
    
    # Transmission parameters
    b = 0.5,                      # Biting rate
    beta_hm = 0.1,                # Transmission probability human to mosquito
    beta_mh = 0.1                 # Transmission probability mosquito to human
  )
  
  # Initial conditions
  N1_initial <- 50000  # Initial population age group 1
  N2_initial <- 10000  # Initial population age group 2
  
  initial_conditions <- c(
    S1 = N1_initial - 10,    # Susceptible age group 1
    E1 = 0,                  # Exposed age group 1
    S2 = N2_initial - 5,     # Susceptible age group 2
    E2 = 0,                  # Exposed age group 2
    Im = 7,                  # Infectious mild (mixed ages)
    Is = 8,                  # Infectious severe (mixed ages)
    C = 0,                   # Chronic
    R = 0,                   # Recovered  
    #V = 0,                   # Vaccinated
    Sm = 9000,               # Susceptible mosquitoes
    Em = 500,                # Exposed mosquitoes
    Im_mosq = 500            # Infectious mosquitoes
  )
  
  # Time points
  time_points <- seq(0, 365*5, by = 1) 
  
  # Run simulation
  cat("Running modified chikungunya simulation...\n")
  results <- run_chikungunya_simulation(params, initial_conditions, time_points)
  
  # Calculate cumulative incidence
  results <- calculate_cumulative_incidence(results, params)
  
  # Calculate RO
  RO <- calculate_RO(params)
  cat(sprintf("Basic reproduction number (RO): %.3f\n", RO))
  
  # Plot results
  cat("Generating plots...\n")
  plot_chikungunya_results(results)
  
  # Print final population sizes
  final_results <- tail(results, 1)
  cat("\nFinal population sizes:\n")
  cat(sprintf("Age group 1 - S&E: %.0f\n", final_results$S1 + final_results$E1))
  cat(sprintf("Age group 2 - S&E: %.0f\n", final_results$S2 + final_results$E2))
  cat(sprintf("Total Infectious: %.0f\n", final_results$Im + final_results$Is))
  cat(sprintf("Chronic cases: %.0f\n", final_results$C))
  cat(sprintf("Recovered: %.0f\n", final_results$R))
  #cat(sprintf("Vaccinated: %.0f\n", final_results$V))
  cat(sprintf("Total mosquitoes: %.0f\n", final_results$Sm + final_results$Em + final_results$Im_mosq))
  
  # Print cumulative incidence
  cat(sprintf("\nCumulative incidence by end of simulation:\n"))
  cat(sprintf("Age group 1: %.0f cases\n", tail(results$cum_incidence_1, 1)))
  cat(sprintf("Age group 2: %.0f cases\n", tail(results$cum_incidence_2, 1)))
  
  return(list(results = results, params = params, RO = RO))
}


# Function for chronic disease burden disease analysis --------------------

analyze_chronic_burden <- function(results){
  
  # Peak Chronic cases
  peak_chronic <- max(results$C)
  peak_chronic_time <- results$time[which.max(results$C)]
  
  #Chronic prevalence over time 
  chronic_data <- data.frame(
    time = results$time,
    chronic_cases = results$C,
    chronic_prevalence = results$C / (results$S1 + results$E1 + results$S2 + results$E2 +
                                        results$Im + results$Is + results$C + results$R)
  )
  
  # Plot chronic burden 
  p_chronic <- ggplot(chronic_data, aes(x= time))+
    geom_line(aes(y= chronic_cases, color = "Chronic Cases"), linewidth = 1.2)+
    geom_line(aes(y = chronic_prevalence* 1000 , color = "Prevalence (per 1000)"),linewidth= 1.2)+
    scale_y_continuous(
      name = "Chronic Cases",
      sec.axis = sec_axis(~./1000, name = "Prevalence (per 1000)")
    )+ labs(title = "Chronic Chikungunya Burden Over Timme",
            x= "Time(days)", color = "Measure")+
    theme_minimal()+
    theme(legend.position = "bottom")
  
  print(p_chronic)
  
  cat(sprintf("Peak chronic cases: %.0f at day %.If\n", peak_chronic, peak_chronic_time))
  
  return(chronic_data)
    
}



# Run the main simulation -------------------------------------------------

if(interactive()){
  simulation_results <- main_simulation()
  cat("Simulation completed successfully!\n")
  
  #Analyze chronic disease burden 
  cat("\nAnalyzing chronic disease burden...\n")
  chronic_analysis <- analyze_chronic_burden(simulation_results$results)
}



# ADDING STOCHASTICITY TO OUR MODEL, YAYY!!! ------------------------------


# Load libraries ----------------------------------------------------------


pacman:: p_load(GillespieSSA2,
                tidyr)



# Define Parameters -------------------------------------------------------

params <- list(
  # Human parameters
  mu = 1/(62*365),              # Human natural death rate
  sigma = 1/5,                  # 1/incubation period (days)
  gamma_m = 1/7,                # Recovery rate from mild infection
  gamma_s = 1/14,               # Recovery rate from severe infection
  gamma_c = 1/60,               # Recovery rate from chronic state (slower)
  severe_prob_1 = 0.1,          # Probability of severe disease in age group 1
  severe_prob_2 = 0.3,          # Probability of severe disease in age group 2
  #chronic_prob_m = 0.15,        # Probability of chronic from mild infection
  chronic_prob_s = 0.25,        # Probability of chronic from severe infection
  mu_ms = 1/10,          # Progression rate from Im to Is (1/10 days)
  #vaccination_rate_1 = 0.001,   # Vaccination rate for age group 1
  #vaccination_rate_2 = 0.002,   # Vaccination rate for age group 2
  contact_rate_1 = 1.0,         # Relative contact rate for age group 1
  contact_rate_2 = 0.8,         # Relative contact rate for age group 2
  amp = 0.3,                    #Seasonality amplitude (0-1)
  phi = 0.2,                   #phase shift (0-1, fraction of year)
  
  # Proportions for age distribution in shared compartments (for demographic calculation)
  prop_age1_Im = 0.8,           # Proportion of mild infectious in age group 1
  prop_age1_Is = 0.7,           # Proportion of severe infectious in age group 1  
  prop_age1_C = 0.75,           # Proportion of chronic in age group 1
  prop_age1_R = 0.8,            # Proportion of recovered in age group 1
  #prop_age1_V = 0.8,            # Proportion of vaccinated in age group 1
  
  # Mosquito parameters
  Lambda_m = 1000,              # Mosquito recruitment rate
  mu_m = 1/14,                  # Mosquito death rate
  sigma_m = 1/7,                # 1/extrinsic incubation period
  
  # Transmission parameters
  b = 0.5,                      # Biting rate
  beta_hm = 0.1,                # Transmission probability human to mosquito
  beta_mh = 0.1                 # Transmission probability mosquito to human
)



# Define Initial State ----------------------------------------------------

# Initial conditions
#N1_initial <- 50000  # Initial population age group 1
#N2_initial <- 10000  # Initial population age group 2

initial_conditions <- c(
  S1 = N1_initial - 10,    # Susceptible age group 1
  E1 = 0,                  # Exposed age group 1
  S2 = N2_initial - 5,     # Susceptible age group 2
  E2 = 0,                  # Exposed age group 2
  Im = 7,                  # Infectious mild (mixed ages)
  Is = 8,                  # Infectious severe (mixed ages)
  C = 0,                   # Chronic
  R = 0,                   # Recovered  
  #V = 0,                   # Vaccinated
  Sm = 9000,               # Susceptible mosquitoes
  Em = 500,                # Exposed mosquitoes
  Im_mosq = 500            # Infectious mosquitoes
)


# Define the Reactions ----------------------------------------------------

#We translate every term in the d/dt equations into a reaction
params$amp <- 0.5  # Example amplitude
params$phi <- 0    # Example phase shift
params$PI <- pi

reactions<- list (
  #Age-group 1
  reaction(
    ~ (1 + amp * cos(2 * PI * (time / (365*5) - phi)))* (b * beta_mh * Im_mosq / (S1+E1+S2+E2+Im+Is+C+R)) *contact_rate_1 * S1, c(S1 = -1, E1 = +1)),
  reaction(~ mu * (S1+E1+Im+Is+C+R), c(S1 =+1)), # Birth (approx)
  reaction(~ mu * S1, c( S1= -1)),   #Natural Death S1
  
  #Progression Age Group 1
  reaction ( ~ sigma * (1- severe_prob_1)* E1, c(E1 = -1,Im = +1)),
  reaction(~sigma * severe_prob_1 *E1, c(E1 = -1, Is = +1)),
  
  #Age-Group 2 
  reaction(
     ~ (1 + amp * cos(2 * PI * (time / (365*5) - phi))) * (b * beta_mh * Im_mosq / (S1+E1+S2+E2+Im+Is+C+R)) * contact_rate_2 * S2, c(S2 = -1, E2 = +1)),
  reaction(~ mu * S2, c(S2 = -1)),
  
  # Progression Age Group 2
  reaction(~ sigma * (1 - severe_prob_2) * E2, c(E2 = -1, Im = +1)),
  reaction(~ sigma * severe_prob_2 * E2, c(E2 = -1, Is = +1)),
  
  # Recovery & Chronic Transitions
  reaction(~ gamma_m * Im, c(Im = -1, R = +1)),
  reaction(~ mu_ms * Im, c(Im = -1, Is = +1)), # Progression mild to severe
  reaction(~ (1 - chronic_prob_s) * gamma_s * Is, c(Is = -1, R = +1)),
  reaction(~ chronic_prob_s * gamma_s * Is, c(Is = -1, C = +1)),
  reaction(~ gamma_c * C, c(C = -1, R = +1)),
  
  # Mosquito Lifecycle
  reaction(
    ~ (1 + amp * cos(2 * PI * (time / (365*5) - phi))) *(b * beta_hm * (Im + Is) / (S1+E1+S2+E2+Im+Is+C+R)) * Sm, c(Sm = -1, Em = +1)),
  reaction(~ sigma_m * Em, c(Em = -1, Im_mosq = +1)),
  reaction(~ mu_m * Sm, c(Sm = -1)),
  reaction(~ mu_m * Em, c(Em = -1)),
  reaction(~ mu_m * Im_mosq, c(Im_mosq = -1))
)

# Run the simulation ------------------------------------------------------

set.seed(123) #For reproducibility

out <- ssa(
  initial_state = initial_conditions,
  reactions = reactions,
  params = unlist(params),
  final_time = 365 * 5, 
  method = ssa_exact()
)

# Plotting ----------------------------------------------------------------

plot_data <- as.data.frame(out$time)|>
  rename(time = 1)|>
  bind_cols(as.data.frame(out$state))|>
  pivot_longer(-time, names_to = 'Compartment', values_to = 'Count')

ggplot(plot_data, aes(x=time, y= Count, color= Compartment))+
  geom_step()+
  facet_wrap(~Compartment, scales = "free_y")+
  theme_minimal()+
  labs(title = "Stochastic Chikunginua Simulation", c= "Days", y = "Count")

out_df <- as.data.frame(out$state)
out_df$time <- as.numeric(rownames(out_df))

# plot human compartments -------------------------------------------------

human_df <- out_df %>%
  mutate(
    I_total = Im + Is
  ) %>%
  select(time, S1, S2, E1, E2, I_total, C, R) %>%
  pivot_longer(-time, names_to = "compartment", values_to = "count")

ggplot(human_df, aes(x = time, y = count, color = compartment)) +
  geom_step(linewidth = 1) +
  labs(
    title = "Stochastic Chikungunya Dynamics (Humans)",
    x = "Time (days)",
    y = "Population count"
  ) +
  theme_minimal()


inf_df <- out_df %>%
  select(time, Im, Is) %>%
  pivot_longer(-time, names_to = "infection_type", values_to = "count")

ggplot(inf_df, aes(x = time, y = count, color = infection_type)) +
  geom_step(size = 1.1) +
  labs(
    title = "Mild vs Severe Infections (Stochastic)",
    x = "Time (days)",
    y = "Number of individuals"
  ) +
  theme_minimal()


# Plot Mosquito dynamics --------------------------------------------------


mosq_df <- out_df %>%
  select(time, Sm, Em, Im_mosq) %>%
  pivot_longer(-time, names_to = "compartment", values_to = "count")

ggplot(mosq_df, aes(x = time, y = count, color = compartment)) +
  geom_step(size = 1) +
  labs(
    title = "Mosquito Population Dynamics",
    x = "Time (days)",
    y = "Mosquito count"
  ) +
  theme_minimal()


runs <- 20

multi_out <- lapply(1:runs, function(i){
  
  sim <- ssa(
    initial_state = initial_conditions,
    reactions = reactions,
    params = unlist(params),
    final_time = 365 * 5,
    method = ssa_exact()
  )
  
  df <- as.data.frame(sim$state)   
  df$time <- sim$time               
  df$run <- i
  
  df
})

multi_df <- bind_rows(multi_out)

multi_df <- multi_df %>%
  mutate(I_total = Im + Is)

ggplot(multi_df, aes(time, I_total, group = run)) +
  geom_step(alpha = 0.3) +
  labs(
    title = "Stochastic Variability in Total Infections",
    x = "Time (days)",
    y = "Total infectious (Im + Is)"
  ) +
  theme_minimal()


t_seq <- seq(5, 365, by = 1)
seas <- 1 + params$amp * cos(2 * phi * (t_seq / 365 - params$phi))

plot(t_seq, seas, type = "l",
     xlab = "Day",
     ylab = "Seasonal multiplier",
     main = "Seasonal forcing function")


















