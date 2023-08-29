# #define variables
# dw <- 0.057
# coverage.current <- 0.76
# coverage.int <- 0.98
# int.effect.mean <- 0.1540
# int.effect.sd <- 0.0237
# age.end <- 100
# iterations <- 1000
seed=1
coverage_1 = 0.5
coverage_2 = 0.8
coverage_3 = 1
# # three scenarios to evaluate 
# sensitivity.scenario <- matrix(c(0.5, 0.03, 0, 0, 0, 0.03), nrow=3, byrow=TRUE)
# colnames(sensitivity.scenario) <- c('cost.wt', 'discount')

# read in data from csv files
incidence <- read.csv('incidence.csv',sep=',')
duration <- read.csv('duration.csv',sep=',')
#d.coverage <- read.csv('coverage.csv',sep=',')
#population.states <- read.csv('populationstates.csv',sep=',') 
#population.start <- read.csv('population5yr.csv',sep=',') 
#population <- read.csv('population.csv',sep=',') 
mortality <- read.csv('mortality.csv',sep=',') 

v_names_states <- c("Healthy", "Dental caries", "Edentulousness", "Mortality")  # state names
n_states <- length(v_names_states)        # number of health states 

#state matrix for each individual at each cycle
shifted_m.M.cov1 <- shifted_m.M.cov2 <- shifted_m.M.cov3 <- m.M.cov1 <- m.M.cov2 <-m.M.cov3 <-matrix(nrow = nrow(mortality), ncol = 102, 
              dimnames = list(paste("ind", 1:nrow(mortality), sep = " "), 
                              paste("cycle", 0:101, sep = " ")))  

Micro_sim <- function(mortality, incidence, coverage, seed, v_names_states){
  #Choose individuals to receive treatment
  treated_individuals <- sample(1:nrow(mortality), nrow(mortality) * coverage)
  #print(treated_individuals)

  for (i in 1: nrow(mortality)){
    age <- as.numeric(mortality[i, 2])
    n_cycles <- 101-age
    v_names_cycles  <- paste("cycle", 0:n_cycles) 
    #Transition probabilities
    p_HC  <- 1-exp(-incidence[i,3])

    v.M_1 <- rep("Healthy", nrow(mortality)) 
    m.M[, 1] <- v.M_1 # indicate the initial health state 
    set.seed(seed + i)                
  
    if (n_cycles > 0){
      for (t in 1:n_cycles) {
        p_HC  <- 1-exp(-incidence[i + t - 1,3])
        p_mort<- 1-exp(-mortality[i + t - 1,3])
        
        # Mortality rate for male and female with age 100
        if (i == 101 || i == 202) {
          p_mort <- 1  # Set mortality to 1 for these individuals
        }
        
        if (i + t - 1 < 15) {
          p_CD <- 0
        } else if (i + t - 1 >= 15 && i + t - 1 <= 34) {
          p_CD <-1-exp(- rbeta(1, 0.0, 0.0003))
        } else if (i + t - 1 >= 35 && i + t - 1 <= 54) {
          p_CD <-1-exp(- rbeta(1, 0.017, 0.0023))
        } else if (i + t - 1 >= 55 && i + t - 1 <= 74) {
          p_CD <-1-exp(-rbeta(1, 0.14, 0.0064))
        } else {
          p_CD <-1-exp(-rbeta(1, 0.36, 0.016))
        }

        #Treatment effect
        if (i %in% treated_individuals) {
          trt <- 0.5 #Receive fluoride (What is the effect of fluoride?)
        } else {
          trt <- 1 #Do not receive fluoride
        }
         #print(paste("individual", i, "at cycle", t, "p_HC:", p_HC, "p_mort:", p_mort, trt)) #DEBUG
          v.p <- Probs(m.M[i, t], p_HC, p_CD, p_mort, trt) # calculate the transition probabilities at cycle t 
        #  print(v.p) #DEBUG
         # print(m.M[i,t]) #DEBUG
          m.M[i, t + 1] <- sample(v_names_states, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      } 
    }
  }
  return(m.M) 
}

m.M.cov1 <- Micro_sim(mortality, incidence, coverage_1, seed, v_names_states)
m.M.cov2 <- Micro_sim(mortality, incidence, coverage_2, seed, v_names_states)
m.M.cov3 <- Micro_sim(mortality, incidence, coverage_3, seed, v_names_states)

Shifted_matrix <- function(m.M, shifted_m.M){
for (i in 1:nrow(m.M)) {
    non_na_indices <- which(!is.na(m.M[i, ]))
    if (i <= 101) {
      shifted_indices <- non_na_indices + (i - 1)
    } else {
      shifted_indices <- non_na_indices + (i - 102)
    }
    shifted_m.M[i, shifted_indices] <- m.M[i, non_na_indices]
    }
  return(shifted_m.M)
}

shifted_m.M.cov1 =Shifted_matrix(m.M.cov1, shifted_m.M.cov1)
shifted_m.M.cov2 = Shifted_matrix(m.M.cov2, shifted_m.M.cov2)
shifted_m.M.cov3 = Shifted_matrix(m.M.cov3, shifted_m.M.cov3)

write.csv(shifted_m.M.cov1, "shifted_m_M_cov1.csv")
write.csv(shifted_m.M.cov2, "shifted_m_M_cov2.csv")
write.csv(shifted_m.M.cov3, "shifted_m_M_cov3.csv")

Trace <- function(shifted_m.M,v_names_states ){
  #for each cycle, calculate the proportion of people at each state
  TR <- t(apply(shifted_m.M, 2, function(x) table(factor(x, levels = v_names_states, ordered = TRUE))))
  n_non_na <- apply(!is.na(shifted_m.M), 2, sum)  # Count non-NA values in each column
  TR <- TR / n_non_na  # Create a distribution trace based on non-NA counts
  rownames(TR) <- paste("Cycle", 0:101, sep = " ")     # Name the rows 
  colnames(TR) <- v_names_states
  return(TR)
}
TR.cov1 <- Trace(shifted_m.M.cov1, v_names_states)
TR.cov2 <- Trace(shifted_m.M.cov2, v_names_states)
TR.cov3 <- Trace(shifted_m.M.cov3, v_names_states)

write.csv(TR.cov1, "trace1.csv")
write.csv(TR.cov2, "trace2.csv")
write.csv(TR.cov3, "trace3.csv")


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














