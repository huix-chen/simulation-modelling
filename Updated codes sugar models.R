#working directory
setwd("Desktop/Data/")

#clear variables
rm(list = ls())      # clear memory (removes all the variables from the workspace)

#####Load packages
if (!require('pacman')) install.packages('pacman'); library(pacman) # use this package to conveniently install other packages
# load (install if required) packages from CRAN
p_load("dplyr", "tidyr", "reshape2", "devtools", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "knitr", "stringr", "diagram", "dampack")                                               
# load (install if required) packages from GitHub
# install_github("DARTH-git/darthtools", force = TRUE) #Uncomment if there is a newer version
p_load_gh("DARTH-git/darthtools")

## General setup
v.age <- seq(5,99,1)
cycle_length    <- 1                             # cycle length equal to one year (use 1/12 for monthly)
v_names_states  <- c("Healthy", "Sick", "Dead")  # state names
n_states        <- length(v_names_states)        # number of health states 

### Discounting factors 
d_c <- 0.03 # annual discount rate for costs 
d_e <- 0.03 # annual discount rate for QALYs

### Strategies 
v_names_str     <- c("Standard of Care",         # store the strategy names
                     "Treatment A", 
                     "Treatment B")  
n_str           <- length(v_names_str)           # number of strategies

## Within-cycle correction (WCC) using Simpson's 1/3 rule 
v_wcc <- gen_wcc(n_cycles = n_cycles,  method = "Simpson1/3")


# Effect sizes
IR_HS_perGram  <- 1.006948
ES_HS_TrtA_Gram<- -17
ES_HS_TrtA     <- IR_HS_perGram^ES_HS_TrtA_Gram

ES_HS_TrtB_Gram<- -8.75
ES_HS_TrtB     <- IR_HS_perGram^ES_HS_TrtB_Gram

ES_HD_Trt_perGram <- 1.007339


# incidence,prevalence,disability and mort rates
library(readxl)
datainputs <- read_excel("Data Input Table.xlsx")
datainputs_male = filter(datainputs, sex=="Male")
datainputs_female = filter(datainputs, sex=="Females")

##To loop through both gender, age from 5 to 99:
inc.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, incidence))
prev.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, prevalence))
disab.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, dr))
mort.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, mortality))
morb.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, morbidity))
popcount.df <- subset(datainputs, age >= 5 & age <= 99, select = c(age, sex, popcount))

### State rewards
#### DALYs 
#c_H       <- 400   # cost of one cycle in healthy state
#c_S       <- 1000  # cost of one cycle in sick state
#c_D       <- 0     # cost of one cycle in dead state
#c_trtA    <- 800   # cost of treatment A (per cycle) in healthy state
#c_trtB    <- 1500  # cost of treatment B (per cycle) in healthy state
#### Utilities
u_H       <- 1     # utility when healthy 
u_S       <- 1-0.01   # utility when sick
u_D       <- 0     # utility when dead

### Discount weight for costs and effects 
v_dwc     <- 1 / ((1 + (d_e * cycle_length)) ^ (0:99))
v_dwe     <- 1 / ((1 + (d_c * cycle_length)) ^ (0:99))

# Initialize the list to store cohort traces and population count for all individuals
l_m_M <- list()
l_pop_M <- list()

for (i in 1:nrow(inc.df)){
  age <- as.numeric(inc.df[i, 1])
  n_cycles <- 99-age
  v_names_cycles  <- paste("cycle", 0:n_cycles) 
  #Transition probabilities
  # probability of becoming sick when healthy, conditional on surviving, under standard of care
  p_HS_SoC  <- 1-exp(-rowSums(inc.df[i,3,drop=FALSE])) 
  # probability of becoming sick when healthy, conditional on surviving, under treatment A
  p_HS_trtA <- 1-exp(-rowSums(inc.df[i,3,drop=FALSE])*ES_HS_TrtA) 
  # probability of becoming sick when healthy, conditional on surviving, under treatment B
  p_HS_trtB <- 1-exp(-rowSums(inc.df[i,3,drop=FALSE])*ES_HS_TrtB) 
  # probabilities of dying when healthy (age-dependent) - this is now a sequence of numbers from life table
  
  ##p_SD was not defined !!
  v_p_HD    <-p_SD<- 1-exp(-rowSums(mort.df[i,3]))
  
  # 04 Construct state-transition models
  
  m_P_diag <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(v_names_states, v_names_states))
  m_P_diag["Healthy", "Sick" ]     = "" 
  m_P_diag["Healthy", "Dead" ]     = ""
  m_P_diag["Healthy", "Healthy" ]  = ""
  m_P_diag["Sick"   , "Dead" ]     = ""
  m_P_diag["Sick"   , "Sick" ]     = ""
  m_P_diag["Dead"   , "Dead" ]     = ""
  layout.fig <- c(2, 1)
  plotmat(t(m_P_diag), t(layout.fig), self.cex = 0.5, curve = 0, arr.pos = 0.8,
          latex = T, arr.type = "curved", relsize = 0.85, box.prop = 0.8,
          cex = 0.8, box.cex = 0.7, lwd = 1)

  ## 04.1 Initial state vector
  start.prev       <- sum(prev.df[i,3])
  
  # Starting prevalence
  v_s_init <- c("Healthy" = 1-start.prev, "Sick" = start.prev, "Dead" = 0)  
  v_s_init
  
  ## 04.2 Initialize cohort traces
  ### Initialize cohort trace for SoC 
  m_M_SoC <- matrix(0, 
                    nrow = n_cycles + 1, ncol = n_states, 
                    dimnames = list(v_names_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_SoC[1, ] <- v_s_init
  
  ## Initialize cohort traces for treatments A and B
  # Structure and initial states are the same as for SoC
  m_M_trtA <- m_M_trtB <- m_M_SoC 
  
  ## Initialize population matrices for all three treatments
  m_pop_trtA <- m_pop_trtB <- m_pop_SoC<- matrix(0, nrow = (n_cycles + 1), ncol = n_states, 
                                                 dimnames = list(v_names_cycles, v_names_states))
  
  #population count at the beginning for each cohort
  m_pop_SoC[1, ] = m_pop_trtA[1, ]=  m_pop_trtB[1, ] = m_M_SoC[1, ] * rowSums(popcount.df[i, 3])
    
  ## Create transition probability matrices for strategy SoC 
  ### Initialize transition probability matrix for strategy SoC 
  # All transitions to a non-death state are assumed to be conditional on survival 
  m_P_SoC  <- matrix(0,
                     nrow = n_states, ncol = n_states,
                     dimnames = list(v_names_states, v_names_states)) 
  ### Fill in matrix 
  # from Healthy
  m_P_SoC["Healthy", "Healthy"] <- (1-v_p_HD) * (1-p_HS_SoC)
  m_P_SoC["Healthy", "Sick"]    <- (1 - v_p_HD) * p_HS_SoC
  m_P_SoC["Healthy", "Dead"]    <-      v_p_HD
  # from Sick
  m_P_SoC["Sick", "Sick"] <- 1 - p_SD
  m_P_SoC["Sick", "Dead"] <-     p_SD
  # from Dead
  m_P_SoC["Dead", "Dead"] <- 1
  
  ## Treatment A
  # Store the same matrix of SoC in Trt A, only replace different parameters
  m_P_trtA <- m_P_SoC 
  
  ## I changed p_HD to v_p_HD for treatment A and treatment B
  
  m_P_trtA["Healthy", "Healthy"] <- (1 - v_p_HD) * (1 - p_HS_trtA)
  m_P_trtA["Healthy", "Sick"]    <- (1 - v_p_HD) *      p_HS_trtA
  
  ## Treatment B
  # Store the same matrix of SoC in Trt B, only replace different parameters
  m_P_trtB <- m_P_SoC
  m_P_trtB["Healthy", "Healthy"] <- (1 - v_p_HD) * (1 - p_HS_trtB)
  m_P_trtB["Healthy", "Sick"]    <- (1 - v_p_HD) *      p_HS_trtB
  
  ## Check if transition probability matrices are valid
  ### Check that transition probabilities are in [0, 1]
  check_transition_probability(m_P_SoC,  verbose = TRUE)
  check_transition_probability(m_P_trtA, verbose = TRUE)
  check_transition_probability(m_P_trtB, verbose = TRUE)
  ### Check that all rows sum to 1
  check_sum_of_transition_array(m_P_SoC,  n_states = n_states, verbose = TRUE)
  check_sum_of_transition_array(m_P_trtA, n_states = n_states, verbose = TRUE)
  check_sum_of_transition_array(m_P_trtB, n_states = n_states, verbose = TRUE)
  
  # Initialize cumulative death counts
  cum_death_SoC <- 0
  cum_death_trtA <- 0
  cum_death_trtB <- 0

  # 05 Run Markov model
  # Iterative solution of time-independent cSTM
  if (n_cycles > 0){
    for (t in 1:n_cycles){
      #update transition matrices after each cycle
      # probability of becoming sick when healthy, conditional on surviving, under standard of care
      p_HS_SoC  <- 1-exp(-rowSums(inc.df[i + t,3,drop=FALSE])) 
      # probability of becoming sick when healthy, conditional on surviving, under treatment A
      p_HS_trtA <- 1-exp(-rowSums(inc.df[i + t,3,drop=FALSE])*ES_HS_TrtA) 
      # probability of becoming sick when healthy, conditional on surviving, under treatment B
      p_HS_trtB <- 1-exp(-rowSums(inc.df[i +t ,3,drop=FALSE])*ES_HS_TrtB) 
      # probability of death
      v_p_HD    <-p_SD<- 1-exp(-rowSums(mort.df[i + t,3]))
      m_P_SoC["Healthy", "Healthy"] <- (1-v_p_HD) * (1-p_HS_SoC)
      m_P_SoC["Healthy", "Sick"]    <- (1 - v_p_HD) * p_HS_SoC
      m_P_SoC["Healthy", "Dead"]    <-      v_p_HD
      # from Sick
      m_P_SoC["Sick", "Sick"] <- 1 - p_SD
      m_P_SoC["Sick", "Dead"] <-     p_SD
      # from Dead
      m_P_SoC["Dead", "Dead"] <- 1
      
      ## Treatment A
      m_P_trtA <- m_P_trtB <-m_P_SoC
      m_P_trtA["Healthy", "Healthy"] <- (1 - v_p_HD) * (1 - p_HS_trtA)
      m_P_trtA["Healthy", "Sick"]    <- (1 - v_p_HD) *      p_HS_trtA
      
      ## Treatment B
      m_P_trtB["Healthy", "Healthy"] <- (1 - v_p_HD) * (1 - p_HS_trtB)
      m_P_trtB["Healthy", "Sick"]    <- (1 - v_p_HD) *      p_HS_trtB
      
      # For SoC
      m_M_SoC [t + 1, ] <- m_M_SoC [t, ] %*% m_P_SoC
      # For treatment A
      m_M_trtA[t + 1, ] <- m_M_trtA[t, ] %*% m_P_trtA
      # For treatment B
      m_M_trtB[t + 1, ] <- m_M_trtB[t, ] %*% m_P_trtB
    
      # Calculate the population in each state by multiplying cohort traces with adjusted popcounts
      m_pop_SoC[t + 1, ] <- m_M_SoC[t + 1, ] * (rowSums(popcount.df[i + t, 3]) - cum_death_SoC)
      m_pop_trtA[t + 1, ] <- m_M_trtA[t + 1, ] * (rowSums(popcount.df[i + t, 3]) - cum_death_trtA)
      m_pop_trtB[t + 1, ] <- m_M_trtB[t + 1, ] * (rowSums(popcount.df[i + t, 3]) - cum_death_trtB)
      
      # Check if the sums of populations are less than 0, and if so, reset populations to 0
      if (sum(m_pop_SoC[t + 1, ]) < 0) {
        m_pop_SoC[t + 1, ] <- 0
      }
      if (sum(m_pop_trtA[t + 1, ]) < 0) {
        m_pop_trtA[t + 1, ] <- 0
      }
      if (sum(m_pop_trtB[t + 1, ]) < 0) {
        m_pop_trtB[t + 1, ] <- 0
      }
    
      # Update the cumulative death counts
      cum_death_SoC <- cum_death_SoC + m_pop_SoC[t + 1, "Dead"]
      cum_death_trtA <- cum_death_trtA + m_pop_trtA[t + 1, "Dead"]
      cum_death_trtB <- cum_death_trtB + m_pop_trtB[t + 1, "Dead"]
    }
  
  }
    # Store the cohort traces in the list for the current individual
    l_m_M[[i]] <- list(m_M_SoC = m_M_SoC, m_M_trtA = m_M_trtA, m_M_trtB = m_M_trtB)
    l_pop_M[[i]] <- list(m_pop_SoC = m_pop_SoC, m_pop_trtA = m_pop_trtA, m_pop_trtB = m_pop_trtB)
  }

l_m_M
l_pop_M

# Initialize lists to store DALYs for each cohort and each treatment strategy
l_dalys_SoC <- list()
l_dalys_trtA <- list()
l_dalys_trtB <- list()

# Iterate through the list of population matrices for each cohort
for (i in 1:nrow(disab.df)) {
  # Get the disability weight for the current cohort
  disab_rate <- disab.df[i, "dr"]
  age <- as.numeric(inc.df[i, 1])
  n_cycles <- 99-age

  # Initialize DALYs for the current cohort and each treatment strategy
  daly_values_SoC <- numeric(n_cycles)
  daly_values_trtA <- numeric(n_cycles)
  daly_values_trtB <- numeric(n_cycles)

  # Iterate through cycles
  for (j in 1:n_cycles ) { 
    # Calculate DALYs for the prevalent cases in the "Sick" state for each strategy at the current cycle
    daly_values_SoC[j] <- l_pop_M[[i]]$m_pop_SoC[j, "Sick"] * disab_rate
    daly_values_trtA[j] <- l_pop_M[[i]]$m_pop_trtA[j , "Sick"] * disab_rate
    daly_values_trtB[j] <- l_pop_M[[i]]$m_pop_trtB[j , "Sick"] * disab_rate
    disab_rate <- disab.df[i+j, "dr"]
  }

  # Store DALYs in the lists for each strategy for the current cohort
  l_dalys_SoC[[i]] <- daly_values_SoC
  l_dalys_trtA[[i]] <- daly_values_trtA
  l_dalys_trtB[[i]] <- daly_values_trtB
}

l_dalys_SoC
l_dalys_trtA
l_dalys_trtB


# Convert nested lists to matrices
mat_list_dalys_SoC <- lapply(l_dalys_SoC, function(x) matrix(unlist(x), ncol = 1))
mat_list_dalys_trtA <- lapply(l_dalys_trtA, function(x) matrix(unlist(x), ncol = 1))
mat_list_dalys_trtB <- lapply(l_dalys_trtB, function(x) matrix(unlist(x), ncol = 1))


calculate_sum_dalys <- function(mat_list_dalys_strategy, cycle_length) {
# Initialize a list to store the calculated sums for each sublist
  sums_list <- list()

  # Loop through the list of sublists
  for (i in 1:length(mat_list_dalys_strategy)) {
    sublist <- mat_list_dalys_strategy[[i]]  # Get the current sublist
    sums <- sapply(seq(1, length(sublist), by = cycle_length), function(start) {
     end <- min(start + cycle_length-1, length(sublist))
     sum(sublist[start:end])
    })
    sums_list[[i]] <- sums  # Store the calculated sums for the current sublist
  }
  return(sums_list)
}

# Sum DALYs for each cohort for 10 cycles, 25 cycles and for lifetime under SoC
calculate_sum_dalys(mat_list_dalys_SoC,10)
calculate_sum_dalys(mat_list_dalys_SoC,25)
calculate_sum_dalys(mat_list_dalys_SoC,length(mat_list_dalys_SoC[[1]]))
# Sum DALYs for each cohort for 10 cycles, 25 cycles and for lifetime under treatment A
calculate_sum_dalys(mat_list_dalys_trtA,10)
calculate_sum_dalys(mat_list_dalys_trtA,25)
calculate_sum_dalys(mat_list_dalys_trtA,length(mat_list_dalys_trtA[[1]]))
# Sum DALYs for each cohort for 10 cycles, 25 cycles and for lifetime under Treatment B
calculate_sum_dalys(mat_list_dalys_trtB,10)
calculate_sum_dalys(mat_list_dalys_trtB,25)
calculate_sum_dalys(mat_list_dalys_trtB, length(mat_list_dalys_trtB[[1]]))


# Initialize lists to store morbidity values for each cohort and each treatment strategy
l_morbidity_SoC <- list()
l_morbidity_trtA <- list()
l_morbidity_trtB <- list()

# Iterate through the list of population matrices for each cohort
for (i in 1:nrow(morb.df)) {
  # Get the morbidity rates for the current cohort
  morbidity_rate <- morb.df[i, "morbidity"]
  age <- as.numeric(inc.df[i, 1])
  n_cycles <- 99-age
  
  # Initialize morbidity values for the current cohort and each treatment strategy
  morbidity_values_SoC <- numeric(n_cycles)
  morbidity_values_trtA <- numeric(n_cycles)
  morbidity_values_trtB <- numeric(n_cycles)
  
  # Iterate through cycles
  for (t in 1:n_cycles) {
    # Calculate morbidity values for people alive (in "Healthy" or "Sick" state) for each strategy at the current cycle
    morbidity_values_SoC[t] <- (l_pop_M[[i]]$m_pop_SoC[t, "Healthy"] + l_pop_M[[i]]$m_pop_SoC[t, "Sick"]) * morbidity_rate
    morbidity_values_trtA[t] <- (l_pop_M[[i]]$m_pop_trtA[t, "Healthy"] + l_pop_M[[i]]$m_pop_trtA[t, "Sick"]) * morbidity_rate
    morbidity_values_trtB[t] <- (l_pop_M[[i]]$m_pop_trtB[t, "Healthy"] + l_pop_M[[i]]$m_pop_trtB[t, "Sick"]) * morbidity_rate
    morbidity_rate <- morb.df[i + t, "morbidity"]
  }
  
  # Store morbidity values in the lists for each strategy for the current individual
  l_morbidity_SoC[[i]] <- morbidity_values_SoC
  l_morbidity_trtA[[i]] <- morbidity_values_trtA
  l_morbidity_trtB[[i]] <- morbidity_values_trtB
}
l_morbidity_SoC
l_morbidity_trtA 
l_morbidity_trtB 

# Initialize a list to store the result of subtraction
haly_SoC <- list()
haly_trtA <- list()
haly_trtB <- list()

# Iterate through each element of the list and perform subtraction
for (i in 1:length(l_morbidity_SoC)) {
  # Perform subtraction element-wise
  haly_SoC[[i]] <- as.numeric(l_morbidity_SoC[[i]]) - as.numeric(l_dalys_SoC[[i]])
  haly_trtA[[i]] <- as.numeric(l_morbidity_trtA[[i]]) - as.numeric(l_dalys_trtA[[i]])
  haly_trtB[[i]] <- as.numeric(l_morbidity_trtB[[i]]) - as.numeric(l_dalys_trtB[[i]])
}

# Print the resulting list
print(haly_SoC)
print(haly_trtA)
print(haly_trtB)

calculate_sum_haly <- function(haly_list, cycle_length) {
  # Initialize a list to store the calculated sums for each cohort
  cohort_sums <- list()
  
  # Loop through the list of HALY values for each cohort
  for (i in 1:length(haly_list)) {
    haly_values <- haly_list[[i]]  # Get the HALY values for the current cohort
    
    # Calculate the number of cycles and initialize a vector to store sums
    num_cycles <- length(haly_values)
    num_sums <- ceiling(num_cycles / cycle_length)
    cohort_sum <- numeric(num_sums)
    
    # Calculate the sum for each set of cycles
    for (j in 0:(num_sums - 1)) {
      start_cycle <- j * cycle_length + 1
      end_cycle <- min(start_cycle + (cycle_length - 1), num_cycles)
      cohort_sum[j + 1] <- sum(haly_values[start_cycle:end_cycle])
    }
    
    cohort_sums[[i]] <- cohort_sum  # Store the calculated sums for the current cohort
  }
  
  return(cohort_sums)
}

# Sum HALYs for each cohort for 10 cycles, 25 cycles and for lifetime under SoC
calculate_sum_haly(haly_SoC, 10)
calculate_sum_haly(haly_SoC, 25)
calculate_sum_haly(haly_SoC, length(haly_SoC[[1]]))

# Sum HALYs for each cohort for 10 cycles, 25 cycles and for lifetime under treatment A
calculate_sum_haly(haly_trtA, 10)
calculate_sum_haly(haly_trtA, 25)
calculate_sum_haly(haly_trtA, length(haly_trtA[[1]]))

# Sum HALYs for each cohort for 10 cycles, 25 cycles and for lifetime under treatment B
calculate_sum_haly(haly_trtB, 10)
calculate_sum_haly(haly_trtB, 25)
calculate_sum_haly(haly_trtB, length(haly_trtB[[1]]))


# Extract the sublist at index 37
haly_male_SoC_42 <- unlist(haly_SoC[37])
haly_male_trtA_42 <- unlist(haly_trtA[37])
haly_male_trtB_42 <- unlist(haly_trtB[37])

haly_female_SoC_42 <- unlist(haly_SoC[132])
haly_female_trtA_42 <- unlist(haly_trtA[132])
haly_female_trtB_42 <- unlist(haly_trtB[132])

haly_diff <- function(haly_gender_strategy_age) {
  n <- length(haly_gender_strategy_age)
  haly_diff <- list()  # Initialize list to store HALY differences
  
  for (i in 1:(n - 1)) {
    if (haly_gender_strategy_age[i + 1] == 0) {
      break  # Exit the loop if the next value becomes 0
    }
    
    haly_diff <- append(haly_diff, haly_gender_strategy_age[i + 1] - haly_gender_strategy_age[i])
  }
  
  return(haly_diff)
}

library(ggplot2)

data_male <- data.frame(
  Cycle = 1:length(unlist(haly_diff(haly_male_SoC_42))),
  haly_male_SoC_42 = unlist(haly_diff(haly_male_SoC_42)),
  haly_male_trtA_42 = unlist(haly_diff(haly_male_trtA_42)),
  haly_male_trtB_42 = unlist(haly_diff(haly_male_trtB_42))
)

data_female <- data.frame(
  Cycle = 1:length(unlist(haly_diff(haly_female_SoC_42))),
  haly_female_SoC_42 = unlist(haly_diff(haly_female_SoC_42)),
  haly_female_trtA_42 = unlist(haly_diff(haly_female_trtA_42)),
  haly_female_trtB_42 = unlist(haly_diff(haly_female_trtB_42))
)

pdf("my_plot.pdf")  # Use pdf() as the graphics device
ggplot() +
  geom_point(data = data_male, aes(x = Cycle, y = haly_male_SoC_42, color = "haly_male_SoC_42")) +
  geom_point(data = data_male, aes(x = Cycle, y = haly_male_trtA_42, color = "haly_male_trtA_42")) +
  geom_point(data = data_male, aes(x = Cycle, y = haly_male_trtB_42, color = "haly_male_trtB_42")) +
  geom_point(data = data_female, aes(x = Cycle, y = haly_female_SoC_42, color = "haly_female_SoC_42")) +
  geom_point(data = data_female, aes(x = Cycle, y = haly_female_trtA_42, color = "haly_female_trtA_42")) +
  geom_point(data = data_female, aes(x = Cycle, y = haly_female_trtB_42, color = "haly_female_trtB_42")) +
  labs(title = "HALY Values for Cycle 42") +
  xlab("Cycle") + ylab("HALY Values") +
  scale_color_manual(values = c(
    "haly_male_SoC_42" = "blue", "haly_male_trtA_42" = "red", "haly_male_trtB_42" = "green",
    "haly_female_SoC_42" = "cyan", "haly_female_trtA_42" = "magenta", "haly_female_trtB_42" = "yellow"
  )) +
  theme_minimal()
dev.off()  # Close the graphics device









# 06 Plot Outputs

## 06.1 Plot the cohort trace for strategies SoC

plot_trace(m_M_SoC) 


## 06.2 Overall Survival (OS)

#Print the overall survival for the Standard of Care

v_os_SoC <- 1 - m_M_SoC[, "Dead"]    # calculate the overall survival (OS) probability
v_os_SoC <- rowSums(m_M_SoC[, 1:2])  # alternative way of calculating the OS probability   

plot(v_os_SoC, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")  # create a simple plot showing the OS

# add grid 
grid(nx = n_cycles, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), 
     equilogs = TRUE) 


## 06.2.1 Life Expectancy (LE)


le_SoC <- sum(v_os_SoC)  # summing probability of OS over time  (i.e. life expectancy)


## 06.2.2 Disease prevalence


v_prev <- m_M_SoC[, "Sick"]/v_os_SoC
plot(v_prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")


# 07 State Rewards 


## Scale by the cycle length 
# Vector of state utilities under strategy SoC
v_u_SoC    <- c(H  = u_H, 
                S  = u_S,
                D  = u_D) * cycle_length
# Vector of state costs under strategy SoC
v_c_SoC    <- c(H  = c_H, 
                S  = c_S,
                D  = c_D) * cycle_length
# Vector of state utilities under treatment A
v_u_trtA   <- c(H  = u_H, 
                S  = u_S, 
                D  = u_D) * cycle_length
# Vector of state costs under treatment A
v_c_trtA   <- c(H  = c_H + c_trtA, 
                S  = c_S, 
                D  = c_D) * cycle_length
# Vector of state utilities under treatment B
v_u_trtB   <- c(H  = u_H, 
                S  = u_S, 
                D  = u_D) * cycle_length
# Vector of state costs under treatment B
v_c_trtB   <- c(H  = c_H + c_trtB, 
                S  = c_S, 
                D  = c_D) * cycle_length

## Store state rewards 
# Store the vectors of state utilities for each strategy in a list 
l_u   <- list(SQ = v_u_SoC,
              A  = v_u_trtA,
              B  = v_u_trtB)
# Store the vectors of state cost for each strategy in a list 
l_c   <- list(SQ = v_c_SoC,
              A  = v_c_trtA,
              B  = v_c_trtB)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str


# 08 Compute expected outcomes 


# Create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

## Loop through each strategy and calculate total utilities and costs 
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  
  ### Expected QALYs and costs per cycle 
  ## Vector of QALYs and Costs
  # Apply state rewards 
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
  v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
  
  ### Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
  # QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  # Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}


# 09 Cost-effectiveness analysis (CEA) 


## Incremental cost-effectiveness ratios (ICERs) 
df_cea <- calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea



## CEA table in proper format 
table_cea <- format_table_cea(df_cea) 
table_cea



## CEA frontier 
plot(df_cea, label = "all", txtsize = 14) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.3))


untrace()