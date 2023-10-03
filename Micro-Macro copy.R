seed=1
coverage_1 = 0.5
coverage_2 = 0.8
coverage_3 = 1
age_end =100
gender_M = "male"
gender_F = "female"
v_names_states <- c("Healthy", "Dental caries", "Edentulousness", "Mortality")  # state names
n_states <- length(v_names_states)        # number of health states 
# read in data from csv files
incidence <- read.csv('updated incidence.csv',sep=',')
population <- read.csv('pop.csv',sep=',')
mortality <- read.csv('updated mortality.csv',sep=',') 

#Calculate transition probabilities
Probs <- function(M_it, p_HC, p_CD, p_mort, trt){ 
  
  v.p.it <- rep(NA, n_states )     # create vector of state transition probabilities
  names(v.p.it) <- v_names_states       # name the vector
  
  # update v.p.it with the appropriate probabilities   
  v.p.it[M_it == "Healthy"]  <- c(1 - p_HC*trt, p_HC*trt, 0, 0)  # transition probabilities when healthy
  v.p.it[M_it == "Dental caries"] <- c(0, 1 - p_CD *trt , p_CD*trt, 0)   # transition probabilities when sick
  v.p.it[M_it == "Edentulousness"]  <- c(0, 0, 1- p_mort , p_mort)  # transition probabilities when dead   
  v.p.it[M_it == "Mortality"]  <- c(0, 0, 0, 1) 
  ifelse(sum(v.p.it) == 1, return(v.p.it), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
} 

# Function to keep the first "Mortality" in each row and set the rest to NA
#once an individual has experienced mortality, they are no longer part of the population, 
keep_first_mortality <- function(matrix, n_individuals, age) {
  result_matrix <- matrix(NA, nrow(matrix), ncol(matrix), 
                          dimnames = list(paste("ind", 1:n_individuals, sep = " "), 
                              paste("cycle", age:age_end, sep = " ")))
  
  for (i in 1:nrow(matrix)) {
    found_mortality <- FALSE
    for (j in 1:ncol(matrix)) {
      if (!is.na(matrix[i, j]) && matrix[i, j] == "Mortality") {
        if (!found_mortality) {
          found_mortality <- TRUE
          result_matrix[i, j] <- "Mortality"
        } else {
          result_matrix[i, j] <- NA
        }
      } else {
        result_matrix[i, j] <- matrix[i, j]
      }
    }
  }
  return(result_matrix)
}

Micro_sim <- function(mortality, incidence, coverage, m.M, n_individuals, age, gender_offset){
  set.seed(1)
  #Randomly choose individuals to receive treatment
  treated_individuals <- sample(1:n_individuals, n_individuals * coverage)
  
  # Perform simulation for all the individuals within the same cohort
  for (i in 1: n_individuals){
    #number of cycles each individual experience
    n_cycles <- age_end - age 
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    #Assume all individuals start with "Healthy"
    v.M_1 <- rep("Healthy", n_individuals) 
    m.M[, 1] <- v.M_1 # indicate the initial health state 
    set.seed(seed + i)                
    
    if (n_cycles > 0){
      for (t in 1:n_cycles) {
        #Update all the transition probabilities after each cycle
        #Transition probability: Healthy -> Caries
        p_HC  <- 1 - exp(-incidence[age + t + gender_offset, 3])
        #Transition probability: Edentulouness -> Mortality
        p_mort <- 1 - exp(-mortality[age + t + gender_offset, 3])
        
        # Mortality rate for male and female with age 100
        if ((age + 1) == 101 ) {
          p_mort <- 1  # Set mortality to 1 for these individuals
        }
        
        #Transition probability: Caires -> Edentulousness
        if ((age + t < 15)) {
          p_CD <- 0
        } else if ((age + t  >= 15 && age + t <= 34)){
          p_CD <- 1 - exp(-rbeta(1, 0.0, 0.0003))
        } else if ((age + t >= 35 && age + t <= 54)) {
          p_CD <- 1 - exp(-rbeta(1, 0.017, 0.0023))
        } else if ((age + t  >= 55 && age + t <= 74)){
          p_CD <- 1 - exp(-rbeta(1, 0.14, 0.0064))
        } else {
          p_CD <- 1 - exp(-rbeta(1, 0.36, 0.016))
        }
        
        #Treatment effect
        if (i %in% treated_individuals) {
          trt <- rnorm(1, 0.154, 0.0237)
        } else {
          trt <- 1 #Do not receive fluoride
        }
        v.p <- Probs(m.M[i, t], p_HC, p_CD, p_mort, trt) # calculate the transition probabilities at cycle t 
        m.M[i, t + 1] <- sample(v_names_states, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      } 
    } 
  }
  cat("Dimensions of m.M:", dim(m.M), "\n")#To trace progress and debug
  m.M <- keep_first_mortality(m.M, n_individuals, age)
  return(m.M)
}

#for each cycle, calculate the proportion of people at each state
Trace <- function(m.M, age ){
  TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
  n_non_na <- apply(!is.na(m.M), 2, sum)  # Count non-NA values in each column
  TR <- TR / n_non_na  # Create a distribution trace based on non-NA counts
  colnames(TR) <- v_names_states
  return(TR)
}

Cohort_sim<- function(prop_each_state, mortality, incidence, population, gender,
                      age, gender_offset){
  m_M_list <- list()
  m_pop_list <- list()
  for (j in 1: nrow(prop_each_state)){
    cohort_age = age+j-1
    n_cycles <- age_end - cohort_age
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    ##print(paste()) is used for debugging
    #print(paste("j-th cohort:", j, "cycle number", n_cycles, "index", cohort_age + gender_offset, "current age", cohort_age))
    ### Initialize cohort matrix and population matrix
    m_M <-m_pop<- matrix(0, nrow = n_cycles + 1, ncol = n_states, 
                         dimnames = list(v_names_cycles, v_names_states))
    #Take the proportion of individuals at each state for each cycles as the initial distribution for the cohort simulation
    m_M[1, ] <- prop_each_state[j,]
    
    #population count at the beginning for each cohort
    m_pop[1, ] = m_M[1, ] * as.numeric(population[cohort_age + gender_offset+1, 3])
    # Initialize transition probability matrix 
    m_P <-matrix(0,nrow = n_states, ncol = n_states,
                 dimnames = list(v_names_states, v_names_states))
    
    cum_death <- 0 #Total number of death
    # Run Markov model
    if (n_cycles > 0){
      for (t in 1:n_cycles){
        #update transition matrices after each cycle
        #Transition probability: Healthy -> Caries
        p_HC  <- 1 - exp(-incidence[cohort_age + t + gender_offset, 3])
        #Transition probability: Edentulouness -> Mortality
        p_mort <- 1 - exp(-mortality[cohort_age + t + gender_offset, 3])
        
        # Mortality rate for male and female with age 100
        if (cohort_age + 1 == 101 ) {
          p_mort <- 1  # Set mortality to 1 for these individuals
        }
        
        #Transition probability: Caires -> Edentulousness
        if ((cohort_age + t < 15)) {
          p_CD <- 0
        } else if ((cohort_age + t >= 15 && cohort_age + t<= 34)){
          p_CD <- 1 - exp(-rbeta(1, 0.0, 0.0003))
        } else if ((cohort_age + t >= 35 && cohort_age + t <= 54)) {
          p_CD <- 1 - exp(-rbeta(1, 0.017, 0.0023))
        } else if ((cohort_age + t >= 55 && cohort_age + t <= 74)){
          p_CD <- 1 - exp(-rbeta(1, 0.14, 0.0064))
        } else {
          p_CD <- 1 - exp(-rbeta(1, 0.36, 0.016))
        }
        
        ### Fill in transition matrix without fluoride
        # from Healthy
        m_P["Healthy", "Healthy"] <- 1 - p_HC
        m_P["Healthy", "Dental caries"] <- p_HC
        # from Dental caries
        m_P["Dental caries", "Dental caries"] <- 1 - p_CD
        m_P["Dental caries", "Edentulousness"] <- p_CD
        # from Edentulousness
        m_P["Edentulousness", "Edentulousness"] <- 1- p_mort
        m_P["Edentulousness", "Mortality"] <- p_mort
        # from Mortality
        m_P["Mortality", "Mortality"] <- 1
        
        # For state matrix
        m_M[t + 1, ] <- m_M[t, ] %*% m_P
        # Calculate the population in each state by multiplying cohort traces with adjusted popcounts
        m_pop[t + 1, ] <- m_M[t + 1, ] * (as.numeric(population[cohort_age + 1 + t, 3])-cum_death)
        
        # Check if the sums of populations are less than 0, and if so, reset populations to 0
        if (sum(m_pop[t + 1, ]) < 0) {
          m_pop[t + 1, ] <- 0
        }
        cum_death <- cum_death + m_pop[t + 1, "Mortality"]
      }
      # Append matrices to the lists
      m_M_list[[j]] <- m_M
      m_pop_list[[j]] <- m_pop
    }else{
      m_M_list[[j]] <- m_M
      m_pop_list[[j]] <- m_pop
    }
  }
  
  # Return lists of matrices for all cohorts
  result <- list("m_M" = m_M_list, "m_pop" = m_pop_list)
  return(result)
}

Micro_cohort <- function(mortality, incidence, coverage, population,gender){
  
  # Determine start index based on gender
  if (gender == "male") {
    n_cohort <- 101
    gender_offset <- 0
  } else if (gender == "female") {
    n_cohort <- 101
    gender_offset <- 101
  } else {
    stop("Invalid gender specified. Use 'male' or 'female'.")
  }
  # Create empty lists to store results for each cohort
  state_matrix_list <- list()
  prop_each_state_list <- list()
  m_M_list_of_list <- list()
  m_pop_list_of_list <- list()
  for (i in 1:n_cohort){
    age <- as.numeric(mortality[i + gender_offset, 2])
    n_cycles <- age_end-age
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    
    n_states <- length(v_names_states)
    num_ind <- as.numeric(population[i, 6]) 
    state.M <- matrix(nrow = num_ind, ncol = age_end + 1 - age,
             dimnames = list(paste("ind", 1:num_ind, sep = " "), 
                             paste("cycle", age:age_end, sep = " ")))
    state_matrix <- Micro_sim(mortality, incidence, coverage, m.M=state.M, n_individuals= num_ind, age, gender_offset)
    
    prop_each_state <- Trace(state_matrix, age)
    
    #for each cohort, there is a state matrix for each individual and proportion of people at each state
    #Store them in a list 
    state_matrix_list[[i]] <- state_matrix
    prop_each_state_list[[i]] <- prop_each_state
    
    result_after_cohortSim <- Cohort_sim(prop_each_state, mortality, incidence, population, gender, age, gender_offset)
    m_M_list_of_list[[i]] <- result_after_cohortSim$m_M
    m_pop_list_of_list[[i]] <- result_after_cohortSim$m_pop
    }
    #return(list(state_matrix = state_matrix_list, prop_each_state = prop_each_state_list))
    return(list("state_matrix" <- state_matrix_list, "prop_each_state" <- prop_each_state_list
                , "m_M_list_of_list" <- m_M_list_of_list, "m_pop_list_of_list" <- m_pop_list_of_list))
}
    
male_50 <- Micro_cohort(mortality, incidence, coverage=coverage_1, population,gender=gender_M)
#male_50[[1]]: For 101 cohorts, there is a matrix to trace the health state for each individual at each cycle.
#male_50[[2]]: For 101 cohorts, there is a matrix showing the proportion of individuals at each state for each cycle.
#male_50[[3]]: For each cohort, there is a distribution showing the proportion of individuals at each state for 
#each cycle. Treat each row as the initial distribution for the cohort simulation model 
#and run a Markov process for it.
#male_50[[4]]: Multiply male_50[[3]] by the poupulation at rach cohort to get the population at each state
male_80 <- Micro_cohort(mortality, incidence, coverage=coverage_2, population,gender=gender_M)
male_100 <- Micro_cohort(mortality, incidence, coverage=coverage_3, population,gender=gender_M)

female_50 <- Micro_cohort(mortality, incidence, coverage=coverage_1, population,gender=gender_F)
female_80 <- Micro_cohort(mortality, incidence, coverage=coverage_2, population,gender=gender_F)
female_100 <- Micro_cohort(mortality, incidence, coverage=coverage_3, population,gender=gender_F)

#Multiply disability rate by column "Dental caries" to calculate DALYs
DALY <- function(mortality, gender, result) {
  # Initialize lists to store DALYs for each cohort 
  list_list_dalys <- list()
  
  # Determine start index based on gender
  if (gender == "male") {
    gender_offset <- 0
  } else if (gender == "female") {
    gender_offset <- 101
  } else {
    stop("Invalid gender specified. Use 'male' or 'female'.")
  }
  n_cohort <- 101
  
  # Iterate through the list of population matrices for each cohort
  for (i in 1:n_cohort) {
    l_dalys <- list()
    # Get the disability weight for the current cohort
    disab_rate <- 0.057
    
    # Iterate through matrices within the current cohort
    for (j in 1:(n_cohort - i + 1)) {
      daly_values <- numeric(nrow(result[[4]][[i]][[j]]))
      
      # Iterate through cycles within the current matrix
      for (k in 1:nrow(result[[4]][[i]][[j]])) {
        # Calculate DALYs for the prevalent cases in the "Dental caries" state for each cycle
        daly_values[k] <- result[[4]][[i]][[j]][k, "Dental caries"] * disab_rate
      }
      # Store DALYs in the lists for each matrix within the current cohort
      l_dalys[[j]] <- daly_values
    }
    # Store the DALYs for the current cohort in the overall list
    list_list_dalys[[i]] <- l_dalys
  }
  return(list_list_dalys)
}

DALY_male_50 = DALY(mortality, "male", male_50)
DALY_male_80 = DALY(mortality, "male", male_80)
DALY_male_100 = DALY(mortality, "male", male_100)

DALY_female_50 = DALY(mortality, "female", female_50)
DALY_female_80 = DALY(mortality, "female", female_80)
DALY_female_100 = DALY(mortality, "female", female_100)


