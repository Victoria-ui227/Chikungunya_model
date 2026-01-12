##############################
## SIS–SI vector–host model ##
##############################
library(deSolve)
library(tidyverse)
library(plotly)


# -------------------------
# Parameters 
# -------------------------
biting = 0.3
seeking = 0.1
infectivity = 0.3
alpha = biting * seeking * infectivity
k = 5 # mosquito to human ratio

params <- list(
  alpha     = alpha,      # transmission parameter from mosquitoes to humans
  beta = 0.5,      # P(human infection per bite from  an infectious mosquito)
  gamma = 1/25,    # human recovery rate
  mu_m  = 1/10,     # mosquito death rate (1/lifespan)
  Nh    = 1000,     # total number of humans and we assume that the system is closed
  Nm    = 1000*k       # mosquito:human ratio  (Nm = k * Nh)
)


# -------------------------
# Model state variables
# Sh, Ih : susceptible & infectious humans
# Sm, Im : susceptible & infectious mosquitoes
# -------------------------

state0 <- c(
  Sh = params$Nh - 10,  # 1% initially infectious
  Ih = 10,
  Sm = params$Nm - 0.02*params$Nm,  # 2% infectious mosquitoes
  Im = 0.02*params$Nm
)

# -------------------------
# ODE system
# -------------------------
sis_si <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih # should remain constant 
    Nm <- Sm + Im  # should remain constant 

    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + gamma * Ih
    dIh <-  (alpha * (Im / Nh)   * Sh) - gamma * Ih
    
    # Mosquitoes
    births  <- mu_m * Nm                   
    dSm <- births - (beta * (Ih / Nh)* Sm) - mu_m * Sm
    dIm <-  (beta * (Ih / Nh)* Sm) - mu_m * Im
    
    list(c(dSh, dIh, dSm, dIm))
  })
}

# -------------------------
# Solve
# -------------------------
times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params))



### plotting

subplot(
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human popn."),
  
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito popn."),
  
  nrows = 1
)

#####################################################################
#       interventions  ##############################################
#####################################################################
#1) put vector intervention - reduce biting rate by 40%
# -------------------------
# Parameters 
# -------------------------
#update the paramters to show the effect of rge vector intervention
biting = 0.3 * (1-0.4) # reduction by 40%
seeking = 0.1
infectivity = 0.3
alpha = biting * seeking * infectivity
k = 5 # mosquito to human ratio

params_vec_int <- list(
  alpha     = alpha,      
  beta = 0.5,      
  gamma = 1/25,   
  mu_m  = 1/10,     
  Nh    = 1000,     
  Nm    = 1000*k       
)
# solve the modle with the updated parameter values
out_vec_int <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params_vec_int))


### plot with the updated vector intervetions model outpus

subplot(
  plot_ly(out_vec_int, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human popn. with vector intervention"),
  
  plot_ly(out_vec_int, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito popn. with vector intervention"),
  
  nrows = 1
)

#calculating the infections averted of implmenting the vector intervention
infections_averted_humans_vec = out$Ih[181] - out_vec_int$Ih[181]
infections_averted_mosq_vec = out$Im[181] - out_vec_int$Im[181]

#2) put treatment intervention - reduces the average infectious period from 25 to 10 days
# -------------------------
# Parameters 
# -------------------------
#update the paramters to show the effect of rge vector intervention
biting = 0.3 
seeking = 0.1
infectivity = 0.3
alpha = biting * seeking * infectivity
k = 5 # mosquito to human ratio

params_treat_int <- list(
  alpha     = alpha,      
  beta = 0.5,      
  gamma = 1/10,    # human recovery rate (reduce the infectious period from 25 days to 10)
  mu_m  = 1/10,     
  Nh    = 1000,     
  Nm    = 1000*k      
)
# solve the model with the updated parameter values
out_treat_int <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params_treat_int))


### plot with the updated treat interventions model outputs

subplot(
  plot_ly(out_treat_int, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human popn. with treatment intervention"),
  
  plot_ly(out_treat_int, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito popn. with treatment intervention"),
  
  nrows = 1
)

#calculating the infections averted of implmenting the vector intervention
infections_averted_humans_treat = out$Ih[181] - out_treat_int$Ih[181]
infections_averted_mosq_treat = out$Im[181] - out_treat_int$Im[181]


#3) put both vector and treatmment interventions 
# -------------------------
# Parameters 
# -------------------------
#update the paramters to show the effect of both interventions
biting = 0.3 *(1-0.4) # reduction by 40%  
seeking = 0.1
infectivity = 0.3
alpha = biting * seeking * infectivity
k = 5 # mosquito to human ratio

params_treat_vec_int <- list(
  alpha     = alpha,      
  beta = 0.5,      
  gamma = 1/10,    # human recovery rate (reduce the infectious period from 25 days to 10)
  mu_m  = 1/10,   
  Nh    = 1000,     
  Nm    = 1000*k      
)
# solve the model with the updated parameter values
out_treat_vec_int <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params_treat_vec_int))


### plot with the updated both vector and treatment interventions model outputs

subplot(
  plot_ly(out_treat_vec_int, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human popn. with  both interventions"),
  
  plot_ly(out_treat_vec_int, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito popn. with both interventions"),
  
  nrows = 1
)

#calculating the infections averted of implementing both interventions
infections_averted_humans_both = out$Ih[181] - out_treat_vec_int$Ih[181]
infections_averted_mosq_both = out$Im[181] - out_treat_vec_int$Im[181]
