seed=1
coverage_1 = 0.5
coverage_2 = 0.8
coverage_3 = 1
age_end =100
# read in data from csv files
incidence <- read.csv('updated incidence.csv',sep=',')
population <- read.csv('updated population.csv',sep=',')
mortality <- read.csv('updated mortality.csv',sep=',') 

v_names_states <- c("Healthy", "Dental caries", "Edentulousness", "Mortality")  # state names
n_states <- length(v_names_states)        # number of health states 

#state matrix for each individual(male and female seperately) at each cycle
n_individuals_per_gender <- 101

shifted_m.M.cov1.male <- shifted_m.M.cov2.male <- shifted_m.M.cov3.male <- 
  shifted_m.M.cov1.female <- shifted_m.M.cov2.female <- shifted_m.M.cov3.female <- 
  m.M.cov1.male <- m.M.cov2.male <- m.M.cov3.male <-
  m.M.cov1.female <- m.M.cov2.female <- m.M.cov3.female <- 
  matrix(nrow = n_individuals_per_gender, ncol = age_end + 1,
         dimnames = list(paste("ind", 1:n_individuals_per_gender, sep = " "), 
                         paste("cycle", 0:age_end, sep = " ")))

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
keep_first_mortality <- function(row) {
  found_mortality <- FALSE
  for (i in 1:length(row)) {
    if (!is.na(row[i]) && row[i] == "Mortality") {
      if (!found_mortality) {
        found_mortality <- TRUE
      } else {
        row[i] <- NA
      }
    }
  }
  return(row)
}

Micro_sim <- function(mortality, incidence, coverage, seed, v_names_states, m.M, gender, shifted_m.M){
  set.seed(1)
  
  # Determine the start index based on gender
  if (gender == "male") {
    n_individuals <- 101
    gender_offset <- 0
  } else if (gender == "female") {
    n_individuals <- 101
    gender_offset <- 101
  } else {
    stop("Invalid gender specified. Use 'male' or 'female'.")
  }
  
  #Randomly choose individuals to receive treatment
  treated_individuals <- sample(1:n_individuals, n_individuals * coverage)
  
  # Perform simulation based on gender
  for (i in 1: n_individuals){
    age <- as.numeric(mortality[i + gender_offset, 2])
    #number of cucles each individual experience
    n_cycles <- age_end - age
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    
    #Assume all individuals start with "Healthy"
    v.M_1 <- rep("Healthy", n_individuals) 
    m.M[, 1] <- v.M_1 # indicate the initial health state 
    set.seed(seed + i)                
    
    if (n_cycles > 0){
      for (t in 1:n_cycles) {
        #Update all the transition probabilities based on their ages
        #Transition probability: Healthy -> Caries
        p_HC  <- 1 - exp(-incidence[i + t - 1 + gender_offset, 3])
        #Transition probability: Edentulouness -> Mortality
        p_mort <- 1 - exp(-mortality[i + t - 1 + gender_offset, 3])
        
        # Mortality rate for male and female with age 100
        if (i == 101 ) {
          p_mort <- 1  # Set mortality to 1 for these individuals
        }
        
        #Transition probability: Caires -> Edentulousness
        if ((i + t - 1 < 15)) {
          p_CD <- 0
        } else if ((i + t - 1 >= 15 && i + t - 1 <= 34)){
          p_CD <- 1 - exp(-rbeta(1, 0.0, 0.0003))
        } else if ((i + t - 1 >= 35 && i + t - 1 <= 54)) {
          p_CD <- 1 - exp(-rbeta(1, 0.017, 0.0023))
        } else if ((i + t - 1 >= 55 && i + t - 1 <= 74)){
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
  for (i in 1:nrow(m.M)) {
    non_na_indices <- which(!is.na(m.M[i, ]))
    shifted_indices <- non_na_indices + (i - 1)
    shifted_m.M[i, shifted_indices] <- m.M[i, non_na_indices]
  }
  # Apply the keep_first_mortality function to each row of shifted_m.M
  shifted_m.M <- t(apply(shifted_m.M, 1, keep_first_mortality))
  return(shifted_m.M)
}

m.M.cov1.male <- Micro_sim(mortality, incidence, coverage_1, seed, v_names_states, m.M=m.M.cov1.male, gender="male", shifted_m.M.cov1.male)
m.M.cov2.male <- Micro_sim(mortality, incidence, coverage_2, seed, v_names_states, m.M=m.M.cov2.male, gender="male", shifted_m.M.cov2.male)
m.M.cov3.male <- Micro_sim(mortality, incidence, coverage_3, seed, v_names_states, m.M=m.M.cov3.male, gender="male", shifted_m.M.cov3.male)

m.M.cov1.female <- Micro_sim(mortality, incidence, coverage_1, seed, v_names_states, m.M=m.M.cov1.female, gender="female", shifted_m.M.cov1.female)
m.M.cov2.female <- Micro_sim(mortality, incidence, coverage_2, seed, v_names_states, m.M=m.M.cov2.female, gender="female", shifted_m.M.cov2.female)
m.M.cov3.female <- Micro_sim(mortality, incidence, coverage_3, seed, v_names_states, m.M=m.M.cov3.female, gender="female", shifted_m.M.cov3.female)

# Save male matrices to CSV
write.csv(m.M.cov1.male, "m.M.cov1.male.csv")
write.csv(m.M.cov2.male, "m.M.cov2.male.csv")
write.csv(m.M.cov3.male, "m.M.cov3.male.csv")

# Save female matrices to CSV
write.csv(m.M.cov1.female, "m.M.cov1.female.csv")
write.csv(m.M.cov2.female, "m.M.cov2.female.csv")
write.csv(m.M.cov3.female, "m.M.cov3.female.csv")

#for each cycle, calculate the proportion of people at each state
Trace <- function(shifted_m.M,v_names_states ){
  TR <- t(apply(shifted_m.M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
  n_non_na <- apply(!is.na(shifted_m.M), 2, sum)  # Count non-NA values in each column
  TR <- TR / n_non_na  # Create a distribution trace based on non-NA counts
  rownames(TR) <- paste("Age", 0:age_end, sep = " ")     # Name the rows 
  colnames(TR) <- v_names_states
  return(TR)
}

TR.cov1.male <- Trace(m.M.cov1.male, v_names_states)
TR.cov2.male <- Trace(m.M.cov2.male, v_names_states)
TR.cov3.male <- Trace(m.M.cov3.male, v_names_states)

TR.cov1.female <- Trace(m.M.cov1.female, v_names_states)
TR.cov2.female <- Trace(m.M.cov2.female, v_names_states)
TR.cov3.female <- Trace(m.M.cov3.female, v_names_states)


cohort_sim <- function(TR, mortality, incidence, coverage, population,gender, v_names_states ){
  
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
  
  m_M_list <- list()
  m_pop_list <- list()
  for (i in 1:n_cohort){
    age <- as.numeric(mortality[i + gender_offset, 2])
    n_cycles <- age_end-age
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    
    n_states=length(v_names_states)
    ### Initialize cohort matrix and population matrix
    m_M <-m_pop<- matrix(0, nrow = n_cycles + 1, ncol = n_states, 
                      dimnames = list(v_names_cycles, v_names_states))
    #Each row in the 'TR' matrix represents the initial distribution of health states for a cohort
    m_M[1, ] <- TR[i,]
    
    #population count at the beginning for each cohort
    m_pop[1, ] = m_M[1, ] * as.numeric(population[i, 3])
    
    # Initialize transition probability matrix 
    m_P <-matrix(0,nrow = n_states, ncol = n_states,
                       dimnames = list(v_names_states, v_names_states)) 
    
    if (runif(1) <= coverage) {
      trt <- rnorm(1, 0.154, 0.0237)  # Cohort receives fluoride treatment
    } else {
      trt <- 1  # Cohort does not receive fluoride treatment
    }
    
    cum_death <- 0 #Total number of death
    # Run Markov model
    if (n_cycles > 0){
      for (t in 1:n_cycles){
        #update transition matrices after each cycle
        #Transition probability: Healthy -> Caries
        p_HC  <- 1 - exp(-incidence[i + t - 1 + gender_offset, 3])
        #Transition probability: Edentulouness -> Mortality
        p_mort <- 1 - exp(-mortality[i + t - 1 + gender_offset, 3])
        
        # Mortality rate for male and female with age 100
        if (i == 101 ) {
          p_mort <- 1  # Set mortality to 1 for these individuals
        }
        
        #Transition probability: Caires -> Edentulousness
        if ((i + t - 1 < 15)) {
          p_CD <- 0
        } else if ((i + t - 1 >= 15 && i + t - 1 <= 34)){
          p_CD <- 1 - exp(-rbeta(1, 0.0, 0.0003))
        } else if ((i + t - 1 >= 35 && i + t - 1 <= 54)) {
          p_CD <- 1 - exp(-rbeta(1, 0.017, 0.0023))
        } else if ((i + t - 1 >= 55 && i + t - 1 <= 74)){
          p_CD <- 1 - exp(-rbeta(1, 0.14, 0.0064))
        } else {
          p_CD <- 1 - exp(-rbeta(1, 0.36, 0.016))
        }
        
        ### Fill in transition matrix without fluoride
        # from Healthy
        m_P["Healthy", "Healthy"] <- 1 - p_HC*trt
        m_P["Healthy", "Dental caries"] <- p_HC*trt
        # from Dental caries
        m_P["Dental caries", "Dental caries"] <- 1 - p_CD*trt
        m_P["Dental caries", "Edentulousness"] <- p_CD*trt
        # from Edentulousness
        m_P["Edentulousness", "Edentulousness"] <- 1- p_mort
        m_P["Edentulousness", "Mortality"] <- p_mort
        # from Mortality
        m_P["Mortality", "Mortality"] <- 1

        # For state matrix
        m_M[t + 1, ] <- m_M[t, ] %*% m_P
        # Calculate the population in each state by multiplying cohort traces with adjusted popcounts
        m_pop[t + 1, ] <- m_M[t + 1, ] * (as.numeric(population[i + t, 3])-cum_death)
    
        # Check if the sums of populations are less than 0, and if so, reset populations to 0
        if (sum(m_pop[t + 1, ]) < 0) {
          m_pop[t + 1, ] <- 0
        }
        cum_death <- cum_death + m_pop[t + 1, "Mortality"]
      }
    }
    # Append matrices to the lists
    m_M_list[[i]] <- m_M
    m_pop_list[[i]] <- m_pop
  }
  # Return lists of matrices for all cohorts
  result <- list("m_M" = m_M_list, "m_pop" = m_pop_list)
  return(result)
}

#Run cohort simulation model for both gender and all 3 coverages
male_50= cohort_sim(TR.cov1.male, mortality, incidence, coverage=0.5, population, "male", v_names_states )
female_50= cohort_sim(TR.cov1.female, mortality, incidence, coverage=0.5, population, "female", v_names_states )

male_80= cohort_sim(TR.cov2.male, mortality, incidence, coverage=0.8, population, "male", v_names_states )
female_80= cohort_sim(TR.cov2.female, mortality, incidence, coverage=0.8, population, "female", v_names_states )

male_100= cohort_sim(TR.cov3.male, mortality, incidence, coverage=1, population, "male", v_names_states )
female_100= cohort_sim(TR.cov3.female, mortality, incidence, coverage=1, population, "female", v_names_states )

#For each cohort, calculate DALY values for each cycle
DALY<- function(mortality, gender, result){
  # Initialize lists to store DALYs for each cohort and each treatment strategy
  l_dalys <- list()
  
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
    # Get the disability weight for the current cohort
    disab_rate <- 0.057
    age <- as.numeric(mortality[i + gender_offset, 2])
    n_cycles <- age_end-age
    
    # Initialize DALYs for the currrent cohort
    daly_values <- numeric(n_cycles)
  
    # Iterate through cycles
    for (j in 1:n_cycles ) { 
      # Calculate DALYs for the prevalent cases in the "Dental caries" state for each strategy at the current cycle
      daly_values[j] <- result$m_pop[[i]][j, "Dental caries"] * disab_rate
    }
    # Store DALYs in the lists for each strategy for the current cohort
    l_dalys[[i]] <- daly_values
  }
  return(l_dalys)
}

DALY(mortality, "male", male_50)
DALY(mortality, "male", male_80)
DALY(mortality, "male", male_100)

DALY(mortality, "female", female_50)
DALY(mortality, "female", female_80)
DALY(mortality, "female", female_100)


