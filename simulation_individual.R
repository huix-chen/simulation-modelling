##Test
v_n   <- c("H","S","D")  
n_s   <- length(v_n) 
u_H       <- 1     # utility when healthy 
u_S       <- 0.5   # utility when sick
u_D       <- 0     # utility when dead
u_Treatment <- 0.75  #utility with treatment

n_i <- 100 # number of simulated individuals
v_M_initial <- rep("H", n_i)  
start_age=20
end_age=30
n_t <- end_age - start_age
d_e <- 0.03
p.HS = 0.5
p.HD = 0.1
p.SD = 0.1

individual_sim(v_M_initial, n_i, n_t, v_n, d_e, coverage=0.9) 
  

individual_sim <- function(v_M_initial, n_i, n_t, v_n, d_e, coverage,TR.out = TRUE,seed = 1) {
  
  v_dwe <- 1 / (1 + d_e) ^ (0:n_t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m_M <- m_E <-  matrix(nrow = n_i, ncol = n_t + 1, 
                               dimnames = list(paste("ind", 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_t, sep = " ")))  
  
  m_M[, 1] <- v_M_initial # initial health state for all individuals  
 
  for (i in 1:n_i) {
    set.seed(seed + i)  # set the seed for every individual for the random number generator
    
    # New treatment assignment process
    if (i <= n_i * coverage) {
      Treatment_i <- TRUE  # Assign treatment to the first 90% of individuals
    } else {
      Treatment_i <- FALSE  # The remaining 10% do not receive treatment
    }
    
    m_E[i, 1] <- Quality(m_M[i, 1], Treatment_i)  # Estimate QALYs per individual for the initial health state conditional on treatment
    
    for (t in 1:n_t) {
      v_p <- Probs(m_M[i, t])  # calculate the transition probabilities at cycle t
      
      m_M[i, t + 1] <- sample(v_n, prob = v_p, size = 1)  # sample the next health state and store that state in matrix m.M
      m_E[i, t + 1] <- Quality(m_M[i, t + 1], Treatment_i)  # estimate QALYs per individual during cycle t + 1 conditional on treatment
      
    }
  } 
  
  te <- m_E %*% v_dwe       # total (discounted) QALYs per individual 
  te_hat <- mean(te)        # average (discounted) QALYs

  # create a trace from the individual trajectories
  if (TR.out == TRUE) { 
    TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE))))
    TR <- TR / n_i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n_t, sep = " ")     # name the rows 
    colnames(TR) <- v_n                                  # name the columns 
  } else {
    TR <- NULL
  }
  
  results <- list(m_M = m_M, m_E = m_E, te = te, te_hat = te_hat, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.
Probs <- function(M_it) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  v.p.it <- rep(NA, n_s)     # create vector of state transition probabilities
  names(v.p.it) <- v_n       # name the vector
  
  # update v.p.it with the appropriate probabilities   
  v.p.it[M_it == "H"]  <- c(1 - p.HS - p.HD, p.HS, p.HD) # transition probabilities when healthy
  v.p.it[M_it == "S"] <- c(0, 1- p.SD, p.SD)   # transition probabilities when sick
  v.p.it[M_it == "D"]  <- c(0, 0, 1) # transition probabilities when dead   
  ifelse(sum(v.p.it) == 1, return(v.p.it), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}       

#utility function
Quality <- function (M_it, Treatment=TRUE, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Treatment: take fluoride or not
  # cl:   cycle length (default is 1)
  u.it <- numeric(length(M_it)) 
  u.it[M_it == "H"]  <- u_H      # update the utility if healthy
  u.it[M_it == "S"] <- Treatment * u_Treatment + (1 - Treatment) * u_S  # update the utility if sick conditional on treatment
  u.it[M_it == "D"] <- u_D     # update the utility if sicker
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}

