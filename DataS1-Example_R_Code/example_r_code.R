# example_r_code_fastverse.R
# Version 1.2.0
# R code for examples demonstrating estimation of uncertain identification using multi-observer methods
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# These 5 examples provide diverse examples of multi-observer methods. Each include simulated survey data, input files specifying initial values and constraints for parameters, a function for conducting statistical analyses, and example statistical output. Models, parameter naming conventions, and methods for executing analyses are described in 'Metadata S1 : Example R code for estimating uncertain species identification using multi-observer methods.' Statistical methods are described in the companion article. Code developed and tested in R version 4.1.

###############################################################################
#               Load R packages
###############################################################################

# # Required packages
# library(dplyr)      # Data frame manipulation
# library(fastverse)  # High-performance statistical computing, data manipulation
# # library(Rfast)      High-performance data analysis
# 
# # Optional package for estimating variances and standard errors of back-transformed parameter estimates
# library(mrds)   # Delta method for arbitrary functions

library(fastverse)

fastverse_detach(data.table, kit, fst)

fastverse_extend(dplyr, mrds, Rfast)

###############################################################################
#               Load R workspace
###############################################################################

# This R workspace contains example survey data (data_1 to data_5), example output of statistical analyses (model_output_1 to model_output_5), functions (example.1.f to example.5.f, below) for conducting example statistical analyses, required supplemental functions (below), and  the 'likelihood_equations' list required for analyses with group size >=1. 

load("example_data.RData")

###############################################################################
#           Required supporting functions
###############################################################################

# These functions are included in the R workspace 'example_data.RData'

# Function: {ztpois.f} Solves for 'lambda' of zero-truncated Poisson distribution. Accepts inputs of mean group size 'size' and 'lambda' parameter of a zero-truncated Poisson distribution, returns the absolute difference between mean group size in 'size' and the mean group size computed from 'lambda'. Solving for difference of 0 gives 'lambda' for the mean group size in 'size'.

ztpois.f <- function(lambda, size){
  abs(size - lambda / (1 - exp(-lambda)))
}

# Function: {ztpois.probmass.f} Probability mass function for zero-truncated Poisson distribution. Accepts inputs of 'the lambda' parameter of a zero-truncated Poisson distribution and vector of observed group sizes 'size', returns vector with probability mass for group sizes 1 to maximum 'observed group size For 'lambda' = 0, returns probability = 1 for group size of 1. 

ztpois.probmass.f <- function(lambda, size) {
  if (lambda == 0) {
    return(as.numeric(size == 1))
  }else{
    lambda[1] ^ (1:max(size)) / ((exp(lambda[1]) - 1) * factorial(1:max(size)))
  }
}

# Function: {penalty.f} Penalty function to -log(likelihood) when sum of probabilities >1. Accepts vector of probabilities and returns vector with penalty term that scales with the sum. 
penalty.f <- function(x) {2000 * (sum(x - 1) + 0.0001) ^ 2}

# Function: {group.probability.constant.f}: Calculates group probabilities (pi) for 2 species with homogeneous or heterogeneous groups. For heterogeneous groups, uses the constant heterogeneous group probability model described in the companion article. Accepts inputs of a vector of true species probabilities 'p', heterogeneous group probability 'm', and a vector of unique sizes of observed groups 'g'. Returns a list with group probabilities 'p' and a penalty term 'pen' to the -log(likelihood) for violating model constraints.

group.probability.constant.f <- function(p, m, g) {
  
  # For all models, calculate group probabilities (pi) from true species probabilities 'p' (psi) accounting for mean group sizes.
  p <- (p / g) / sum(p / g)
  pen <- 0 # Penalty term for violating model constraints.
  
  ## Probabilities for heterogeneous groups. Here, 'p' is "delta" (group probabilities prior to heterogeneous groups forming).
  if (any(m > 0)) {
    
    # Transform constant heterogeneous group probability parameter (pi.12) to proportion of all groups (prior to heterogeneous groups forming) that each species contributes to heterogeneous groups. 'm' is equivalent to "epsilon" and 'mix_g' to "epsilon / (1 + epsilon)".  
    mix_g <- m / (1 + m) 
    
    # For each species, test if heterogeneous group probability 'mix_g' exceeds 'p' (i.e, heterogeneous group probabilities for a species exceed overall group probability for that species). If so, reduce heterogeneous groups and apply penalty term to -log(likelihood). 
    if (sum(p < mix_g) > 0) { 
      test <- mix_g - p
      p[which(test > 0)] <- p[which(test > 0)] + max(test) + 0.01
      p[which(test < 0)] <- 1 - p[which(test > 0)]
      test_diff <- test[which(test > 0)] + 1
      for (i in 1:length(test_diff)) {
        pen <- pen + penalty.f(test_diff[i])
      }
    }
    
    # Calculate group probability (pi) for each species (pi.1, pi.2) and for heterogeneous groups (pi.12).
    p <- c((p - mix_g) / (1 - mix_g), m)
  }
  return(list(p, pen))
}

# Function: {group.true.probability.key.f}: With heterogeneous groups for 2 species, computes probabilities for all possible true groups. Accepts inputs of group probabilities 'group_probability', probability mass of group sizes (1 to the maximum observed group size) for each species 'size_probability', and observed group sizes 'size'. Returns a list with named elements for each unique group size composed of a vector of probabilities for each possible true group. 

group.true.probability.key.f <- function(group_probability, size_probability, size){
  
  # Calculate probabilities first for homogeneous groups at the start/end of each vector, then add probabilities for heterogeneous groups with decreasing numbers of species 1 and increasing numbers of species 2.
  
  B <- dim(size_probability)[2]
  
  output <-
    vapply(1:B, function(x)
      group_probability[x] * size_probability[size, x], numeric(length(size))) %>%
    {
      lapply(size, function(x)
        c(
          .[which(size == x) , 1],
          group_probability[3] * size_probability[0:(x - 1), 2] * size_probability[(x - 1):0, 1],
          .[which(size == x) , 2]
        ))
    } %>%
    {
      lapply(1:length(.), function(x)
        setColnames(.[[x]], paste0(size[x]:0, 0:size[x])))
    } %>%
    setColnames(., paste(size))
  
}

# Function: {group.observed.cprobability.f} Computes conditional (on the true group) probabilities for observed groups with heterogeneous groups. Accepts inputs of an observed group 'd', a vector of classification probabilities 'p', the number of possible combinations (with repetition) of true groups 'states', and group size 's'. Returns vector with a conditional probability for each true group.
# This function relies on an external list object 'likelihood_equations' (provided in the file 'likelihood_equations_a2b2g25.RData') with pre-computed values

group.observed.cprobability.f <- function(d, p, states, s){  
  # Return probability = 1 for missing observed groups (i.e. no classifications).
  if (s < 1) return(rep(1, states))
  
  # Compute probabilities for each combination of true states (possible true groups) and each permutation of observation states in the observed group by raising classification probabilities to the exponent of the numbers of individuals in corresponding observation states and then multiplying by the multinomial coefficient for each permutation.

  # Compute overall probabilities for possible true states by summing across permutations of observation states in the observed group with identical index values.
  
  t1 <- qM(likelihood_equations[[paste0("observed.", paste0(d, collapse = "."))]])
  
  t2 <-
    dapply(t1[, -1], MARGIN = 1, function(z)
      prod(p ^ z[-1] , z[1])) %>%
    fsum(., t1[, 1]) 
  
}

# Function: {group.observed.cprobability.homogeneous.f} Computes conditional (on true group) probabilities for homogeneous observed groups. Accepts inputs of an observed group 'd', a vector of classification probabilities 'p', and the number of possible combinations (with repetition) of true groups 'states', and group size 's'. Returns vector with a conditional probability for each true group.
# This function relies on an external list object 'likelihood_equations' (provided in the file 'likelihood_equations_a2b2g25.RData' and in the R workspace 'example_data.RData') with pre-computed values.

group.observed.cprobability.homogeneous.f <- function(d, p, states, s){
  # Return probability = 1 for missing observed groups (i.e. no classifications).
  if (s < 1) return(rep(1, states))
  
  # Compute probabilities for each combination of true states (possible true groups) and each permutation of observation states for the observed group by raising classification probabilities to the exponent of the numbers of individuals in corresponding observation states and then multiplying by the multinomial coefficient for each permutation.
  
  t1 <- qM(likelihood_equations[[paste0("observed.", paste0(d, collapse = "."))]])
  
  # Row numbers in 'homogeneous' exclude heterogeneous groups 
  
  homogeneous <- 
    which(likelihood_equations[[paste0("true.permutations.count.g.", s)]] == 1) %>%
    vapply(., function(x) which(t1[, 1] == x), numeric(1)) 
  
  t2 <-
    dapply(t1[homogeneous, -1], MARGIN = 1, function(z)
      prod(p ^ z[-1] , z[1])) %>%
    fsum(., t1[homogeneous, 1])
  
}

# Function: {mlogit.regress.predict.f} Predicts group-level probabilities for multinomial logistic regression. Accepts input of vector of covariate values 'cov_mat' and a matrix of beta parameters 'beta_mat' (columns for intercept and slope, rows for each logit regression). Returns a matrix with species-specific (columns 1 to B) predicted probabilities for each covariate value (rows).

mlogit.regress.predict.f <- function(cov_mat, beta_mat) {

  # tidy up the inputs to matrices with 2 covariates, 3 beta parameters
  
  if (!is.matrix(cov_mat))
    cov_mat <- qM(cov_mat)
  if (is.vector(beta_mat))
    beta_mat <- matrix(beta_mat, byrow = TRUE, ncol = 2)
  else if (is.array(beta_mat))
    beta_mat <- qM(beta_mat)
  
  n_cat <- dim(beta_mat)[1] + 1
  
  # Compute distribution of predicted values by category
  
  distribution <-
    t(vapply(1:dim(beta_mat)[1], function(b)
      exp(beta_mat[b, 1] + beta_mat[b, 2] * cov_mat[, 1]), 
      numeric(dim(cov_mat)[1]) )) %>%
    {
      . %r/% (1 + fsum(.))
    } %>%
    {
      t(rbind(., 1 - colsums(.)))
    } 
  
  distribution
  
}

###############################################################################
#           Example 1: 3 species
###############################################################################

# Function: {example.1.f} Example code for a MOM (multi-observation method) model with 3 observation states, 3 true species states, groups of size 1, 1 primary observer, and 3 identical secondary observers. Accepts inputs of initial parameter values 'parameters' and survey observation data 'dat', returns the -log(likelihood) (plus any penalty term). The function includes an example of a penalty function that adds a penalty term to the -log(likelihood) if the sum of the true species probability parameters exceeds 1. 

# Required function(s) and object(s)
# data_1  (data frame)

example.1.f <- function(parameters, dat) {
  
  ## Parameter import, formatting, and constraints ---------
  
  # Back-transform probability parameters from logit scale to probability scale
  parameters <- plogis(parameters)
  
  # Array 'theta_arr' contains matrices with classification probabilities (theta) for each observer. 
  # Dimension 1 (column) = true species states 1 to 3, dimension 2 (row) = observation states 1 to 3, dimension 3 (matrix) = observers 1 to 4 (1 primary, 3 identical secondary)
  
  theta_arr <- array(0, dim = c(3, 3, 4), 
                     dimnames = list(c(paste0("spp_", 1:3)), c(paste0("class_", 1:3)), c("obs_p", paste0("obs_s", 1:3))))
  
  theta_arr[, , 1] <- 
    c(
      1 - sum(parameters[1:2]),
      parameters[c(3, 5, 1)],
      1 - sum(parameters[3:4]),
      parameters[c(6, 2, 4)],
      1 - sum(parameters[5:6])
    )
  theta_arr[, , 2:4] <-
    c(
      1 - sum(parameters[7:8]),
      parameters[c(9, 11, 7)],
      1 - sum(parameters[9:10]),
      parameters[c(12, 8, 10)],
      1 - sum(parameters[11:12])
    )
  
  # True species probabilities (psi) for species 1 and 2
  psi <- c(parameters[13:14])
  
  # Enforce constraint that probabilities for true species parameters sum to 1. If violated, the 'penalty' term increases with greater deviation from 1 according the the penalty function, and true species probabilities are reduced to sum to 1 to avoid an inadmissible (negative) value of true species probability for species 3 (1 - sum of species 1 and 2). 
  penalty <- 0
  
  if (sum(psi) > 1) {
    penalty <- 2000 * (sum(sum(psi) - 1) + 0.0001) ^ 2
    psi <- (psi / sum(psi)) * 0.999
  }
  
  # Vector 'psi' now contains true species probabilities for species 1 to 3 and sums to 1
  psi <- 
    c(parameters[13:14], 1 - sum(parameters[13:14])) %>%
    setColnames(., c(paste0("spp_", 1:3)))
  
  ## Compute likelihood for each unique observation history  ---------
  
  nrow_dat <- dim(dat)[1]
  
  # For each unique observation history, the array 'group_observation_probability' gives probabilities (conditional on possible true species states) for observations states of each observer. 
  # Dimension 1 = unique observation histories, dimension 2 = true species states 1 to 3, dimension 3 = observers 1 to 4
  
  group_observation_probability <- 
    array(0,
          dim = c(nrow_dat, 3, 4),
          dimnames = list(c(paste0(
            "obs_history_", 1:nrow_dat
          )), c(paste0("spp_", 1:3)), c("obs_p", paste0("obs_s", 1:3))))
  
  # Add conditional probabilities for each observer
  for (obs in 1:4) {
    # Matrix 'dat_obs' contains subset of observed groups for the current observer 'obs'
    
    dat_obs <- dat[, (((obs - 1) * 3) + 1):(obs * 3)]
    
    group_observation_probability[, , obs] <-
      t(apply(dat_obs, 1, function(x)
        apply(theta_arr[, , obs], 1, function(y)
          dmultinom(x, prob = y))))
  }
  
  # For each true species state 1 to 3 (columns), the matrix 'probability_mat' gives probabilities of each unique observation history (rows), computed as the product of probabilities for observation states of each observer and of true species probabilities (psi)
  
  probability_mat <-
    matrix(0, nrow_dat, 3,
           dimnames = (list(c(
             paste0("obs_history_", 1:nrow_dat)), 
             c(paste0("spp_", 1:3)) )))
  
  # Add probabilities for each true species state
  
  for (b in 1:3) {
    probability_mat[, b] <-
      apply(group_observation_probability[, b,], 1, function(x)
        psi[b] * prod(x))
  }
  
  # Likelihoods for each unique observation history are the summed probabilities across true species states
  likelihood <- rowsums(probability_mat)
  
  # Compute the -log(likelihood) as the sum of the product of the likelihood and count for each unique observation history. Add penalty term for violating constraint(s).
  sum(dat[, "count"] * -log(likelihood)) + penalty
}

# Specify initial parameter values. Entered probability values are logit-transformed before analyses.

parameters_ini_1 <-
  qlogis(c(seq(0.2, 0.1, length.out = 12), 0.3, 0.25))

# Specify lower box constraints for parameter values
constraints_low_1 <- c(rep(-9.21024, 14))

names(parameters_ini_1) <- names(constraints_low_1) <- 
  c(
    "theta_p_21", "theta_p_31", "theta_p_12", "theta_p_32", "theta_p_13", "theta_p_23",
    "theta_s1_21", "theta_s1_31", "theta_s1_12", "theta_s1_32", "theta_s1_13", "theta_s1_23",
    "psi_1", "psi_2"
  )

# The 'optim' function solves for parameter values minimizing the -log(likelihood) (plus any penalty terms). The vector 'constraints_low_1' specifies lower box constraints. 

model_1 <-
  optim(
    parameters_ini_1,
    example.1.f,
    gr = NULL,
    dat = data_1,
    hessian = TRUE,
    method = c("L-BFGS-B"),
    lower = constraints_low_1
  )

(model_1)

# Estimated variance-covariance matrix and standard errors of estimated parameters, computed from the Hessian matrix.

model_1_var_covar <- solve(model_1$hessian)

model_1_estimates_se <- 
  matrix(
    c(model_1$par,
      (diag(solve(model_1$hessian))) ^ 0.5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("estimate", "SE"), names(model_1$par))
  )

(model_1_estimates_se)

# The 'DeltaMethod' function in the optional R library 'mrds' conveniently estimates variances and standard errors when back-transforming estimated logit parameters to probabilities. 

model_1_estimates_probabilities_se <-
  matrix(
    c(plogis(model_1$par),
      (diag(
        DeltaMethod(
          as.vector(model_1$par),
          plogis,
          solve(model_1$hessian),
          0.000001
        )$variance) ^ 0.5)),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("estimate", "SE"), names(model_1$par))
  )

(model_1_estimates_probabilities_se)

###############################################################################
#           Example 2: 2 species, homogeneous groups
###############################################################################

# Function: {example.2.f} Example code for a SOM (single observation method) model with 2 observation states, 2 true species states, homogeneous groups, and 4 identical secondary observers. Accepts inputs of initial parameter values 'parameters' and survey observation data 'dat', returns the -log(likelihood). 

# Required function(s) and object(s)
# data_2                                    (data frame)
# likelihood_equations                      (list)
# ztpois.f                                  (function)
# ztpois.probmass.f                         (function)
# group.probability.constant.f              (function)
# group.observed.cprobability.homogeneous.f (function)
# penalty.f                                 (function)

example.2.f <- function(parameters, dat) {
  
  ## Data and parameter import, formatting ----------
  
  # Summarize group sizes from data and from group size parameters
  
  g <- unique(dat$group_size) # Observed group sizes
  g1 <- sum(dat$group_size == 1) # Count of unique observation histories with group size = 1
  g_spp <- parameters[4:5] # Mean group size parameters for species 1 and 2
  
  # Back-transform probability parameters from logit scale to probability scale
  
  parameters <- plogis(parameters[1:3])
  
  # Vector 'psi' contains true species probabilities for species 1 and 2
  
  psi <- 
    c(parameters[3], 1 - parameters[3]) %>%
    setColnames(c(paste0("spp_", 1:2)))
  
  # Matrix 'theta_mat' summarizes classification probabilities (theta) of each observer as a vector for efficient computation
  # Dimension 1 (column) = classification probabilities, dimension 2 (row) = observers 1 to 4 (4 secondary observers with identical classification probabilities)
  
  theta_mat <- 
    matrix(
      rep(c(1 - parameters[1], parameters[2], parameters[1], 1 - parameters[2]), 4),
      ncol = 4,
      byrow = TRUE,
      dimnames = list(c(paste0("obs_s", 1:4)), 
                      c("theta_11", "theta_12", "theta_21", "theta_22"))
    )
  
  # Vector 'lambda' contains lambda parameters (defining mean group size) of a zero-truncated Poisson distribution for species 1 and 2
  lambda <- rep(0, 2)
  
  lambda[1:2] <-
    vapply(g_spp[1:2], function(x)
      optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))
  
  names(lambda) <- c(paste0("spp_", 1:2))
  
  ## For each species, compute group probabilities and probability mass of group sizes ----------
  
  # Matrix 'group_size_probmass' contains probability mass (mu) of group sizes (rows, from 1 to the maximum observed group size) for each true species state 1 to 2 (column).
  
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = 2) %>%
    setDimnames(., list(c(paste0("group_size_", 1:dim(.)[1])), 
                        c(paste0("spp_", 1:2))))
  
  # Function 'group_probability.constant.f' calculates group probabilities (pi) for each true species state, with the '0' function argument specifying homogeneous groups.
  
  output <- 
    group.probability.constant.f(psi, 0, g_spp)
  
  group_probability <- output[[1]]
  
  ## Compute likelihood for each unique observation history ----------
  # Unique observation histories with group size = 1
  
  # Vector 'likelihood' for probabilities of each unique observation history
  
  likelihood <- numeric(dim(dat)[1])
  
  # For each unique observation history, array 'group_observation_probability' contains conditional (on possible true species states) probabilities for observed groups of each observer. 
  # Dimension 1 (row) = unique observation histories, dimension 2 (column) = true species state 1 to 2, dimension 3 (matrix) = observers 1 to 4 (4 identical secondary)
  
  group_observation_probability <- 
    array(0, dim = c(g1, 2, 4),
          dimnames = list(c(paste0("obs_history_", 1:g1)), 
                          c(paste0("spp_", 1:2)), c(paste0("obs_s", 1:4))))
  
  # Add conditional probabilities for observed groups of each observer.
  for (obs in 1:4) {
    
    # Matrix 'dat_obs' contains subset of observed groups for the current observer 'obs' and group size = 1
    
    dat_obs <- 
      fselect(dat[1:g1, ], (((obs - 1) * 2) + 1):(obs * 2))
    
    group_observation_probability[, , obs] <- 
      t(apply(dat_obs, 1, function(x)
        apply(matrix(theta_mat[obs, ], ncol = 2), 1, function(y)
          dmultinom(x, prob = y))))
  }
  
  # For each true species state 1 to 2 (columns), matrix 'probability_mat' gives probabilities of each unique observation history (rows), calculated as the product of probabilities for observed groups of each observer and of true species probabilities (psi)
  
  probability_mat <-
    matrix(0, g1, 2,
           dimnames = (list(c(paste0("obs_history_", 1:g1)), 
                            c(paste0("spp_", 1:2)))) )
  
  # Add probabilities for each true species state b
  
  for (b in 1:2) {
    probability_mat[, b] <-
      apply(group_observation_probability[, b,], 1, function(x)
        group_probability[b] * group_size_probmass[1, b] * prod(x))
  }
  
  # Likelihoods for each unique observation history are the summed probabilities across true species states
  
  likelihood[1:g1] <- rowsums(probability_mat)
  
  # Unique observation histories with group sizes > 1 ---------
  
  # Matrix 'n_group_size' gives the count of unique observation histories 'n' for each unique 'group_size'
  
  n_group_size <-  
    fgroup_by(dat, group_size) %>%
    fselect(., "count") %>%
    fnobs(.)
  
  # Compute likelihoods for each group size > 1
  
  for (i in g[g > 1]) {
    
    # For the current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
    
    n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]]
    rows_i <- which(dat$group_size == i) 
    
    # Columns in matrix 'probability_mat' contain conditional (on possible true groups) probabilities for observed groups of each observer and probability for one possible true group. Rows correspond to each possible true group for each unique observation history.
    
    probability_mat <-
      matrix(
        c(rep(numeric(1), n_i * 2 * 4),
          rep(group_size_probmass[i, ] * group_probability, n_i)), 
        ncol =  5, 
        dimnames = list(c(paste0("group_observed_", rep(rows_i, each = 2))), 
                        c(paste0("obs_s", 1:4), "psi")))
    
    # Add conditional probabilities for observed groups of each observer.
    
    for (obs in 1:4) {
      
      # Matrix 'dat_tmp' contains subset of observed groups for the current observer 'obs' and group size 'i'
      
      dat_tmp <-
        qM(fselect(dat[rows_i, ], c( (((obs - 1) * 2) + 1):(obs * 2) )))
      
      observation_sum <- rowsums(dat_tmp)
      
      probability_mat[, obs] <-
        as.vector(vapply(1:n_i, function(x)
          group.observed.cprobability.homogeneous.f(dat_tmp[x,],
                                                    theta_mat[obs,],
                                                    2,
                                                    observation_sum[x]),
          numeric(2))) 
    }
    
    # Integers in 'index' associate rows in 'probability_mat' with unique observation histories.  
    
    index <- rep(1:n_i, each = 2)
    
    # Add likelihoods for the current group size to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of probabilities in each row of 'probability_mat' summed across all possible combinations of true groups for each unique observation history specified by 'index'. 
    
    likelihood[rows_i] <-
      rowprods(probability_mat) %>%
      fsum(., index)
  }
  
  # Compute the -log(likelihood) as the sum of the product of the likelihood and count for each unique observation history. 
  sum(dat[, "count"] * -log(likelihood)) 
}

# Specify initial parameter values. Probabilities (thetas, psi) are logit-transformed before analyses, and group size parameters are mean group size.

parameters_ini_2 <- c(qlogis(c(0.2, 0.1, 0.4)), 1.2, 1.3)

# Specify lower and upper constraints for parameter values. Upper constraints of infinity (Inf) result in estimates without upper constraints. 

constraints_low_2 <- c(rep(-9.2, 3), rep(1.002, 2))
constraints_up_2 <- c(rep(-1.4, 2), rep(Inf, 3))
names(parameters_ini_2) <- names(constraints_low_2) <- names(constraints_up_2) <- 
  c("theta_s1_21",	"theta_s1_12", 	"psi_1",	"g_1",	"g_2")

# The 'optim' function solves for parameter values minimizing the -log(likelihood). Vectors 'constraints_low_2' and 'constraints_up_2' specify lower and upper box constraints. 

model_2 <- 
  optim(
    parameters_ini_2,
    example.2.f,
    gr = NULL,
    dat = data_2,
    hessian = TRUE,
    method = c("L-BFGS-B"),
    lower = constraints_low_2,
    upper = constraints_up_2
  )

(model_2)

###############################################################################
#           Example 3: Covariate predicting true species probabilities 
###############################################################################

# Function: {example.3.f} Example code for a MOM model with 3 observation states; 2 true species states; groups of size 1; a continuous, covariate predicting true species probabilities (psi); a primary observer; and 4 identical secondary observers. Accepts inputs of initial parameter values 'parameters' and survey observation data 'dat', returns the -log(likelihood). 
# A continuous, group-level, standard normal covariate is used to predict true species probabilities (psi) via multinomial logistic regression with an intercept coefficient and a slope coefficient. True species state 2 is the baseline category. 

# Required function(s) and object(s)
# data_3                      (data frame)
# mlogit.regress.predict.f    (function)

example.3.f <- function(parameters, dat) {
  
  ## Parameter import and formatting ---------
  
  # Extract regression coefficients predicting true species probability to vector 'psi_betas'
  
  psi_betas <- parameters[1:2]
  
  # Back-transform probability parameters from logit scale to probability scale
  
  parameters <- plogis(parameters[3:10])
  
  # Array 'theta_arr' contains matrices with classification probabilities (theta) for each observer. 
  # Dimension 1 (column) = true species states 1 to 2, dimension 2 (row) = observation states 1 to 3, dimension 3 (matrix) = observers 1 to 5 (1 primary, 4 secondary with identical classification probabilities)
  
  theta_arr <-
    array(0, c(2, 3, 5),
          dimnames = list(c(paste0("spp_", 1:2)), 
                          c(paste0("class_", 1:3)), 
                          c("obs_p", paste0("obs_s", 1:4))))
  
  theta_arr[, , 1] <-
    c(1 - sum(parameters[1:2]), 
      parameters[c(3, 1)],
      1 - sum(parameters[3:4]),
      parameters[c(2, 4)])
  
  theta_arr[, , 2:5] <-
    c(1 - sum(parameters[5:6]),
      parameters[c(7, 5)],
      1 - sum(parameters[7:8]),
      parameters[c(6, 8)])
  
  # Matrix 'psi_group_mat' gives group-level true species probabilities for each true species state predicted from multinomial logit regression coefficients and group-level covariate values
  
  nrow_dat <- dim(dat)[1]
  
  psi_group_mat <- 
    mlogit.regress.predict.f(dat$covariate_psi, psi_betas) %>%
    setDimnames(., list( c(paste0("obs_history_", 1:nrow_dat)), 
                         c("psi_1", "psi_2")))
  
  ## Compute likelihood for each unique observation history ---------
  
  # Matrix 'probability_mat' contains conditional (on possible true species states) probabilities of observation states for each observer (columns 1 to 5) and probabilities of true species states (last column). Rows correspond to each possible true species state for each unique observation history.
  
  probability_mat <-
    matrix(
      c(rep(numeric(1), nrow_dat * 2 * 5),
        t(psi_group_mat)), 
      ncol =  6,
      dimnames = list(c(paste0("obs_history_", rep(1:nrow_dat, each = 2))), 
                      c("obs_p", paste0("obs_s", 1:4), "psi"))
    )
  
  # Add conditional probabilities of observation states for each observer
  
  for (obs in 1:5) {
    
    # Matrix 'dat_obs' contains subset of observations for the current observer 'obs'
    
    dat_obs <-
      fselect(dat, (((obs - 1) * 3) + 1):(obs * 3))
    
    probability_mat[, obs] <-
      as.vector(apply(dat_obs, 1, function(x)
        apply(theta_arr[, , obs], 1, function(y)
          dmultinom(x, prob = y))))
  }
  
  # Integers in 'index' associate rows in 'probability_mat' with unique observation histories 
  
  index <- 
    rep(1:nrow_dat, each = 2)
  
  # Likelihoods for each unique observation history are calculated as the product of probabilities in each row of 'probability_mat' summed across all possible true species states for each unique observation history indexed by 'index'
  
  likelihood <-
    rowprods(probability_mat) %>%
    fsum(., index)
  
  # Compute the -log(likelihood) as the sum of likelihoods for each unique observation history. 
  
  sum(-log(likelihood))
}

# Specify initial parameter values. Probabilities (thetas) are logit-transformed before analyses.

parameters_ini_3 <-
  c(-1, 0.5, qlogis(seq(0.25, 0.1, length.out = 8)))

# Specify lower constraints for parameter values

constraints_low_3 <- c(-5, -5, rep(-9.2, 8))
names(parameters_ini_3) <- names(constraints_low_3) <- c(
  "b0_psi_1",	"b1_psi_1",	
  "theta_p_21",	"theta_p_31",	"theta_p_12",	"theta_p_32",	
  "theta_s1_21",	"theta_s1_31",	"theta_s1_12",	"theta_s1_32"
)

# The 'optim' function solves for parameter values minimizing the -log(likelihood). Vector 'constraints_low_3' specifies lower box constraints. 

model_3 <-
  optim(
    parameters_ini_3,
    example.3.f,
    gr = NULL,
    dat = data_3,
    hessian = TRUE,
    method = c("L-BFGS-B"),
    lower = constraints_low_3,
    control = list(trace = 3)
  )

(model_3)

###############################################################################
#           Example 4: Covariate predicting classification, heterogeneous groups
###############################################################################

# Function: {example.4.f} Example code for a MOM model with 2 observation states; 2 true species states; heterogeneous groups; a covariate predicting classification probabilities (theta); a primary observer; and 4 identical secondary observers. Accepts inputs of initial parameter values 'parameters' and survey observation data 'dat', returns the -log(likelihood). 
# A continuous, group-level, standard normal covariate predicts classification probabilities (theta) via separate multinomial logistic regressions for primary versus secondary observers. Each regression includes an intercept coefficient and a slope coefficient, and probability of correct identification is the baseline category. 

# Required function(s) and object(s)
# data_4                          (data frame)
# likelihood_equations            (list)
# mlogit.regress.predict.f        (function)
# ztpois.f                        (function)
# ztpois.probmass.f               (function)
# group.probability.constant.f    (function)
# group.true.probability.key.f    (function)
# group.observed.cprobability.f   (function)
# penalty.f                       (function)

example.4.f <- function(parameters, dat) {
  
  ## Data and parameter import, formatting ---------
  
  # Summarize group sizes
  
  g <- unique(dat$group_size) # Observed group sizes
  g_spp <- parameters[10:11] # Mean group size parameters for species 1 and 2
  
  # True species probabilities (psi) for species 1 and 2 and heterogeneous group probability (pi) are back-transformed from logit scale to probability scale
  
  psi <- 
    c(plogis(parameters[9]), 1 - plogis(parameters[9]))
  
  pi <- plogis(parameters[12])
  
  # Array 'theta_betas_arr' organizes regression coefficients (betas) predicting classification probabilities (theta), with a matrix for each observer
  # Dimension 1 (row) = true species states 1 to 2, dimension 2 (column) = beta.0 (intercept) and beta.1 (slope) coefficients, dimension 3 (matrix) = observers 1 to 5 (1 primary, 4 secondary with identical classification probabilities)
  
  theta_betas_arr <-
    array(
      c(parameters[1:2], parameters[5:6], parameters[3:4], parameters[7:8]),
      dim = c(2, 2, 2),
      dimnames = (list(
        c(paste0("spp_", 1:2)),  
        c(paste0("beta_", 0:1)) , 
        c("obs_p", "obs_s")
      ))
    )
  
  # Array 'theta_arr' contains matrices with classification probabilities (theta) for each observer and each unique observation history. 
  # Dimension 1 (column) = classification probabilities, dimension 2 (row) = unique observation histories, dimension 3 (matrix) = observers 1 to 5 (1 primary, 4 identical secondary)
  
  nrow_dat <- dim(dat)[1]
  
  theta_arr <- 
    array(0,
          dim = c(4, nrow_dat, 5),
          dimnames = list(
            c("theta_p_11",	"theta_p_12",	"theta_s1_21",	"theta_s1_22"),
            c(paste0("obs_history_", 1:nrow_dat)),
            c("obs_p", paste0("obs_s", 1:4))
          ))
  
  for (b in 1:2) {
    theta_arr[c(3, 1), , b] <- 
      t(mlogit.regress.predict.f(dat$covariate_theta, 
                                 theta_betas_arr[1, , b, drop = FALSE]))
    
    theta_arr[c(2, 4), , b] <- 
      t(mlogit.regress.predict.f(dat$covariate_theta, 
                                 theta_betas_arr[2, , b, drop = FALSE]))
  }
  
  # Secondary observers are assumed to have identical classification probabilities
  
  theta_arr[, , 5] <- theta_arr[, , 4] <- theta_arr[, , 3]  <- theta_arr[,  , 2]
  
  # Vector 'lambda' contains lambda parameters (defining mean group size) of a zero-truncated Poisson distribution for true species states 1 to 2.
  
  lambda <- rep(0, 2)
  
  lambda[1:2] <-
    vapply(g_spp[1:2], function(x)
      optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1)) %>%
    setColnames(., c(paste0("spp_", 1:2)))
  
  ## Compute group probabilities and group size probability mass for each species ----------
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each group size (row, from 1 to the maximum observed group size) and for each true species state 1 to 2 (column)
  
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = 2) %>%
    setDimnames(., list(c(paste0("group_size_", 1:dim(.)[1])), 
                        c(paste0("spp_", 1:2))))
  
  # Function 'group.probability.constant.f' computes group probabilities (pi) for each true species, assuming a constant heterogeneous group probability parameter (pi.12) as described in the companion article
  # Vector 'group_probability' gives group probabilities for homogeneous groups of species 1 to 2 (pi.1 and pi.2) and for heterogeneous groups (pi.12). Penalty term 'penalty' is >0 if the heterogeneous group probability takes an inadmissible value (i.e., greater than the available groups).
  
  output <- 
    group.probability.constant.f(psi, pi, g_spp)
  
  group_probability <- output[[1]] %>%
    setColnames(., c("pi_1", "pi_2", "pi_12"))
  
  penalty <- output[[2]]
  
  ## Compute likelihood for each unique observation history ----------
  
  # Vector 'likelihood' for the likelihood of each unique observation history
  
  likelihood <- numeric(nrow_dat)
  
  # Matrix 'n_group_size' gives the count of unique observation histories 'n' for each unique 'group_size'
  
  n_group_size <-  
    fgroup_by(dat, group_size) %>%
    fselect(., "covariate_theta") %>%
    fnobs(.)
  
  # List 'group.true.probability' contains probabilities for each possible true group
  # List elements (named with each group size) are vectors with probabilities for each possible true group for that group size (named by the counts of species 1 and 2). 
  
  group.true.probability <- 
    group.true.probability.key.f(group_probability, group_size_probmass, g)
  
  # Compute likelihoods for each observed group size
  
  for (i in g) {
    
    # For the current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
    
    n_i <- 
      n_group_size[[pmatch(i, n_group_size$group_size), 2]]
    
    rows_i <- which(dat$group_size == i) 
    
    B_states <- length(group.true.probability[[paste(i)]]) # Number of combinations of true species states for the current group size 
    
    # Matrix 'probability_mat' contains conditional probabilities for observed groups of each observer and probabilities for possible true groups. For each unique observation history, separate rows correspond to each possible true group.  
    # Dimension 1 (row) = unique observation histories (separate rows for each true group), dimension 2 (column) = probabilities for observers 1 to 4 and for the true group ('psi')
    
    probability_mat <-
      matrix(c(rep(numeric(1), n_i * B_states * 5),
               rep(group.true.probability[[paste(i)]], n_i)),
             ncol =  6,
             dimnames = list(c(paste0("obs_history_", rep(rows_i, each = B_states))), 
                             c("obs_p", paste0("obs_s", 1:4), "psi")))
    
    # Add condition probabilities for observed groups of each observer
    
    for (obs in 1:5) {
      
      # Matrix 'data.tmp' contains subset of observed groups for the current observer 'obs' and group size 'i'
      
      dat_tmp <- 
        qM(fselect(dat[rows_i, ], (((obs - 1) * 2) + 1):(obs * 2)) )
      
      observation_sum <- rowsums(dat_tmp)
      
      probability_mat[, obs] <- 
        as.vector(vapply(1:n_i, function(x) 
          group.observed.cprobability.f(dat_tmp[x, ], 
                                        theta_arr[, rows_i[x], obs], 
                                        B_states, 
                                        observation_sum[x]),
          numeric(B_states))) 
    }
    
    # Integers in 'index' associate rows in 'probability_mat' with unique observation histories
    
    index <- rep(1:n_i, each = B_states)
    
    # Add likelihoods for the current group sizes to 'likelihood'. Likelihoods for each unique observation history are computed as the product of probabilities in each row of 'probability_mat' summed across all possible combinations of true groups for each unique observation history indexed by 'index'.
    
    likelihood[rows_i] <-
      rowprods(probability_mat) %>%
      fsum(., index)
  }
  
  # Compute the -log(likelihood) as the sum of likelihoods for each unique observation history. Add penalty term for violating constraint(s).
  
  sum(-log(likelihood)) + penalty
}

# Specify initial parameter values. True species probability (psi) and heterogeneous group probability (pi) are logit-transformed before analyses, and group size parameters are mean group size.

parameters_ini_4 <-
  c(c(seq(-3,-5, length.out = 4), 
      seq(-1, 2, length.out = 4), 
      qlogis(0.4)),  
    1.2,  1.3, 
    qlogis(0.01))

# Specify lower constraints for parameter values

constraints_low_4 <- c(rep(-9.2, 4), rep(-5, 4), -9.2, rep(1.002, 2), -7)

names(parameters_ini_4) <- names(constraints_low_4) <- 
  c(
    "b0_theta_p_21", "b0_theta_p_12", "b0_theta_s1_21", "b0_theta_s1_12",
    "b1_theta_p_21", "b1_theta_p_12", "b1_theta_s1_21", "b1_theta_s1_12",
    "psi_1",
    "g_1","g_2",
    "pi_12"
  )

# The 'optim' function solves for parameter values minimizing the -log(likelihood). Vector 'constraints_low_4' specifies lower box constraints. 

model_4 <- 
  optim(
    parameters_ini_4,
    example.4.f,
    gr = NULL,
    dat = data_4,
    hessian = TRUE,
    method = c("L-BFGS-B"),
    lower = constraints_low_4,
    control = list(trace = 3)
  )

(model_4)

###############################################################################
#           Example 5: Distinct observers
###############################################################################

# Function: {example.5.f} Example code for a MOM model with 2 observation states, 2 true species states, a primary observer; and 3 distinct secondary observers. Accepts inputs of initial parameter values 'parameters' and survey observation data 'dat', returns the -log(likelihood). 

# Required function(s) and object(s)
# data_5        (data frame)

example.5.f <- function(parameters, dat) {
  
  ## Data and parameter import, formatting ----------
  
  # Back-transform probability parameters from logit scale to probability scale
  
  parameters <- plogis(parameters)
  
  # Array 'theta_arr' contains matrices with classification probabilities (theta) for each observer
  # Dimension 1 (column) = true species states 1 to 2, dimension 2 (row) = observation states 1 to 2, dimension 3 (matrix) = observers 1 to 4 (1 primary, 3 secondary with distinct classification probabilities)
  
  theta_arr <- array(0, dim = c(2, 2, 4), 
                     dimnames = list(c(paste0("spp_", 1:2)), 
                                     c(paste0("class_", 1:2)), 
                                     c("obs_p", paste0("obs_s", 1:3))))
  
  theta_arr[, , 1] <- c(1 - parameters[1], parameters[2:1], 1 - parameters[2])
  theta_arr[, , 2] <- c(1 - parameters[3], parameters[4:3], 1 - parameters[4])
  theta_arr[, , 3] <- c(1 - parameters[5], parameters[6:5], 1 - parameters[6])
  theta_arr[, , 4] <- c(1 - parameters[7], parameters[8:7], 1 - parameters[8])
  
  # True species probabilities (psi) for true species states 1 and 2
  
  psi <- 
    c(parameters[9], 1 - parameters[9]) %>%
    setColnames(c(paste0("spp_", 1:2)))
  
  ## Compute likelihoods for unique observation histories ----------
  
  # For each unique observation history, array 'group_observation_probability' gives probabilities (conditional on  true species states) for observed groups of each observer
  # Dimension 1 = unique observation histories, dimension 2 = true species states 1 to 2, dimension 3 = observers 1 to 4 (1 primary, 3 distinct secondary)
  
  nrow_dat <- dim(dat)[1]
  
  group_observation_probability <- 
    array(0,
          dim = c(nrow_dat, 2, 4),
          dimnames = list(c(paste0("obs_history_", 1:nrow_dat)), 
                          c(paste0("spp_", 1:2)), 
                          c("obs_p", paste0("obs_s", 1:3))) )
  
  # Add conditional probabilities for each observer
  
  for (obs in 1:4) {
    
    # Matrix 'dat_obs' contains subset of observed groups for the current observer 'obs'
    
    dat_obs <-
      qM( fselect(dat, (((obs - 1) * 2) + 1):(obs * 2)) )
    
    group_observation_probability[, , obs] <-
      t(apply(dat_obs, 1, function(x)
        apply(theta_arr[, , obs], 1, function(y)
          dmultinom(x, prob = y))))
  }
  
  # For each true species state 1 to 3 (columns), matrix 'probability_mat' gives probabilities for each unique observation history (rows), calculated as the product of probabilities for observed groups of each observer and of true species probabilities (psi)
  
  probability_mat <- 
    matrix(0, nrow_dat, 2,
           dimnames = (list(c(paste0("obs_history_", 1:nrow_dat)), 
                            c(paste0("spp_", 1:2)) )) )
  
  # Add probabilities for each true species state
  
  for (b in 1:2) {
    probability_mat[, b] <-
      apply(group_observation_probability[, b,], 1, function(x)
        psi[b] * prod(x))
  }
  
  # Likelihoods for each unique observation history are the summed probabilities across true species states
  
  likelihood <- rowsums(probability_mat)
  
  # Compute the -log(likelihood) as the sum of the product of the likelihood and count for each unique observation history. 
  
  sum(dat[, "count"] * -log(likelihood))
}

# Specify initial parameter values. Probabilities (theta, psi) are logit-transformed before analyses.

parameters_ini_5 <- 
  qlogis(c(seq(0.2, 0.1, length.out = 8), 0.4))

# Specify lower constraints for parameter values

constraints_low_5 <- c(rep(-9.21024, 9))

names(parameters_ini_5) <- names(constraints_low_5) <- 
  c(
    "theta_p_21", "theta_p_12",
    "theta_s1_21","theta_s1_12",
    "theta_s2_21","theta_s2_12",
    "theta_s3_21", "theta_s3_12",
    "psi_1"
  )

# The 'optim' function solves for parameter values minimizing the -log(likelihood). Vector 'constraints_low_5' specifies lower box constraints.

model_5 <-
  optim(
    parameters_ini_5,
    example.5.f,
    gr = NULL,
    dat = data_5,
    hessian = TRUE,
    method = c("L-BFGS-B"),
    lower = constraints_low_5
  )

(model_5)

