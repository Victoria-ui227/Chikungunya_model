#The following model is focused on having a uniform SEIR-SEI model that does not divide the model into the population dynamics and also the infectious are in one group

#Required libraries
library(pacman)

pacman::p_load(deSolve, ggplot2, reshape2, dplyr, gridExtra)

#Function
seir_sei_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #Total human populations
    N_h <- S_h + E_h + I_h + R_h
    
    #Force of infection from mosquitoes to humans
    lambda_h <- b * beta_vh * I_v / N_h
    
    #Force of infection from humans to mosquitoes
    lambda_v <- b * beta_hv * I_h / N_h
    
    # Human population dynamics
    dS_h <- mu_h * N_h - lambda_h * S_h - mu_h * S_h
    dE_h <- lambda_h * S_h - sigma_h * E_h - mu_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h - mu_h * I_h
    dR_h <- gamma_h * I_h - mu_h * R_h
    
    # Vector population dynamics
    dS_v <- - lambda_v * S_v - mu_v * S_v
    dE_v <- lambda_v * S_v - sigma_v * E_v - mu_v * E_v
    dI_v <- sigma_v * E_v - mu_v * I_v
    
    list(c(dS_h, dE_h, dI_h, dR_h, dS_v, dE_v, dI_v))
  })
}

#Function to run simulation
run_seir_sei_simulation <- function(params, initial_state, time_points) {
  out <- ode(
    y = initial_state, 
    times = time_points, 
    func = seir_sei_model, 
    parms = params,
    method = "lsoda"
  )
  
  #convert to dataframe
  out <- as.data.frame(out)
  
  #Add column names
  colnames(out) <- c("time", "S_h", "E_h", "I_h", "R_h", "S_v", "E_v", "I_v")
  
  return(out)
}

#Function to plot results
plot_seir_sei_results <- function(results) {
  
  #Susceptible and Exposed humans
  data <- results |>
    dplyr::select(time, S_h, E_h, I_h, R_h) |>
    reshape2::melt(id.vars = "time",
         variable.name = "Compartment",
         value.name = "Population")
  
  data$Compartment <- factor(
    data$Compartment,
    levels = c("S_h", "E_h", "I_h", "R_h"),
    labels = c("Susceptible", "Exposed", "Infectious", "Recovered")
  )
  
  p1 <- ggplot(data, aes(x = time, y = Population, color = Compartment)) +
    geom_line(linewidth = 1.2) +
    labs(title = "Disease States",
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 3: Mosquito Population
  mosquito_data <- results |>
    dplyr::select(time, S_v, E_v, I_v) |>
    reshape2::melt(id.vars = "time",
         variable.name = "Compartment",
         value.name = "Population")
  
  mosquito_data$Compartment <- factor(
    mosquito_data$Compartment,
    levels = c("S_v", "E_v", "I_v"),
    labels = c("Susceptible", "Exposed", "Infectious")
  )
  
  p2 <- ggplot(mosquito_data, aes(x = time, y = Population, color = Compartment)) +
    geom_line(linewidth = 1.2) +
    labs(title = "Mosquito Population",
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  #Key Epidemiological Indicators
  epi_data <- results |>
    dplyr::select(time, S_h, I_h) |>
    reshape2::melt(id.vars = "time",
         variable.name = "Indicator",
         value.name = "Population")
  
  epi_data$Indicator <- factor(
    epi_data$Indicator,
    levels = c("S_h", "I_h"),
    labels = c("Susceptible", "Infectious")
  )
  
  p3 <- ggplot(epi_data, aes(x = time, y = Population, color = Indicator)) +
    geom_line(linewidth = 1.2) +
    labs(title = "Key Epidemiological Indicators",
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Arrange plots
  grid.arrange(p1, p2, p3, ncol = 2)
}

# calculate basic reproduction number
calculate_R0 <- function(parameters) {
  with(as.list(parameters), {
    R0_hv <- (b * beta_hv * sigma_v) / ((mu_v + sigma_v) * mu_v)
    R0_vh <- (b * beta_vh * sigma_h) / ((mu_h + gamma_h) * (mu_h + sigma_h))
    R0 <- sqrt(R0_hv * R0_vh)
    return(R0)
  })
}

#calculate the cumulative incidence
calculate_cumulative_incidence <- function(results) {
  initial_infected <- results$I_h[1]
  final_infected <- tail(results$I_h, n = 1)
  cumulative_incidence <- final_infected - initial_infected
  return(cumulative_incidence)
}

main_simulation <- function() {
  # Define parameters
  parameters <- c(
    #human parameters
    mu_h = 1/(62*365),   # Natural death rate of humans
    sigma_h = 1/5,       # 1/incubation period
    gamma_h = 1/7,       # Recovery rate in humans
    
    #vector parameters
    lambda_v = 0.1,        # Recruitment rate of vectors
    mu_v = 1/14,         # Natural death rate of vectors
    sigma_v = 1/5,       # 1/EIP
    
    #Transmission parameters
    beta_vh = 0.442,       # Transmission rate from vector to human
    beta_hv = 0.333,       # Transmission rate from human to vector
    b = 0.41            # Biting rate
 
  )
  #Initial conditions
  N_intial<-1450000 #initial human population
    
  # Initial state
  initial_state <- c(
    S_h = N_intial,   # Susceptible humans
    E_h = 0,     # Exposed humans
    I_h = 589,    # Infectious humans
    R_h = 0,     # Recovered humans
    S_v = 4327,  # Susceptible vectors
    E_v =1912,     # Exposed vectors
    I_v = 636     # Infectious vectors
  )
  
  # Time points
  time_points <- seq(0, 59, by = 1)
  
  # Run simulation
  cat("Running the chikungunya simulation---\n")
  results <- run_seir_sei_simulation(parameters, initial_state, time_points)
  
  #calculate R0
  R0 <- calculate_R0(parameters)
  cat(sprintf("Basic reproduction number (R0): %.2f\n", R0))
  
  # Plot results
  cat("Generating plots---\n")
  plot_seir_sei_results(results)
  
  #print the final population sizes
  final_results <- tail(results, n = 1)
  cat("Final population sizes:\n")
  cat(sprintf("Susceptible Humans: %.0f\n", final_results$S_h))
  cat(sprintf("Exposed Humans: %.0f\n", final_results$E_h))
  cat(sprintf("Infectious Humans: %.0f\n", final_results$I_h))
  cat(sprintf("Recovered Humans: %.0f\n", final_results$R_h))
  
  cat(sprintf("Susceptible Vectors: %.0f\n", final_results$S_v))
  cat(sprintf("Exposed Vectors: %.0f\n", final_results$E_v))
  cat(sprintf("Infectious Vectors: %.0f\n", final_results$I_v))
  
  return(list(results = results, params = parameters, R0 = R0, final_results = final_results))
}

analyze_burden<-function(results){
  
  # Peak infection cases
  peak_infection <- max(results$I_h)
  peak_infection_time <- results$time[which.max(results$I_h)]
  
  # Chronic prevalence over time
  infection_data <- data.frame(
    time = results$time,
    infection_cases = results$I_h,
    infectious_prevalence= results$I_h / (results$S_h + results$E_h + results$I_h + results$R_h))
  
  #plot the infectious prevalence over time
  p_infectious_prevalence <- ggplot(infection_data, aes(x = time, y = infectious_prevalence)) +
    geom_line(color = "blue", size = 1.2) +
    labs(title = "Infectious Prevalence Over Time",
         x = "Time (days)", y = "Infectious Prevalence") +
    theme_minimal()
  
  print(p_infectious_prevalence)
  
  cat(sprintf("Peak infectious cases: %.0f at day %.1f\n", peak_infection, peak_infection_time))

    
  
}

# Run the main simulation
if (interactive()) {
  simulation_results <- main_simulation()
  cat("Simulation completed successfully!\n")
  
  #Analyze infection disease burden
  cat("Analyzing disease burden---\n")
  analyze_burden(simulation_results$results)
}
