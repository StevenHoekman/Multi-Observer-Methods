# murrelet_r_code.R
# Version 1.1.0
# R code for estimating uncertain identification for 2 species of murrelets 
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# This code conducts statistical analyses using multi-observation method models to estimate uncertain identification for 2 species of murrelets during line transect surveys in Glacier Bay, Alaska, USA during July 2014. Details of conducting analyses and of parameter naming conventions are provided in Metadata S4: Input Files and R Computer Code for Analyses of Glacier Bay Murrelet Survey Data. Descriptions of the statistical model and its parameters are provided in the companion article in Appendix S5: Methods and Results for Application of MOM Models to Line Transect Surveys for Murrelets. Code was developed and tested using R version 4.1.

###############################################################################
#             Load R packages
###############################################################################

library(dplyr)  # Data frame manipulation

###############################################################################
#             Load R workspace
###############################################################################

# The R workspace 'murrelet_data.RData' contains R objects and functions for statistical analyses, including: supplemental functions for analyses (below), the function 'murrelet.model.f' (below), and the 'likelihood.equations' list provided in the file 'likelihood_equations_a3b2g25.RData'

load("murrelet_data.RData")

###############################################################################
#           Required supplemental functions
###############################################################################

# These functions are included in the R workspace 'murrelet_data.RData'

# Function: {mlogit.regress.predict.f} Predicts group-level probabilities for multinomial logistic regression. Accepts input of vector of covariate values 'r', a matrix of beta parameters 'beta' (columns for intercept and slope, rows for each logit regression), and vector with number of species 'n'. Returns a matrix with species-specific (columns 1 to B) predicted probabilities for each covariate value (rows).

mlogit.regress.predict.f <- function(r, beta, n) {
  # Make beta vectors into matrices
  if (is.vector(beta)) {
    beta <- matrix(beta, ncol = 2, byrow = T)
  }
  distribution <- matrix(0, length(r), n)
  denom <- matrix(0, n, length(r))
  
  beta <- t(apply(beta, 1, function(x) x))
  
  # 'denom' is summation term in denominator
  denom[1:(n - 1), ] <-
    t(apply(beta, 1, function(x)
      exp(x[1] + x[2] * r)))
  
  denom <- 1 + colSums(denom)
  
  # Compute probabilities for non-baseline categories
  distribution[, 1:(n - 1)] <-
    apply(beta, 1, function(x)
      exp(x[1] + x[2] * r) / denom)
  
  # Add probability for baseline category
  distribution[ , n] <- 1 - rowSums(distribution)
  distribution
}

# Function: {group.probability.psi.constant.f}: With heterogeneous groups for 2 species, computes group probabilities (pi) with a group-level covariate predicting true species probabilities (psi) for each unique observation history. Accepts inputs of a matrix of predicted true species probabilities 'p' for each observation history (rows), the heterogeneous group probability parameter 'm', and a vector of unique sizes for observed groups 'g'. Returns a list with a matrix of group probabilities for each observation history (rows), the heterogeneous group probability parameter, and a penalty term to be applied to the -log(likelihood) for violating model constraints.

group.probability.psi.constant.f <- function(p, m, g) {
  
  # For all models with group sizes >1, calculate group probabilities (pi) from true species probabilities 'p' (psi) accounting for mean group sizes.
  if (any(g > 1)) {
    p <- t(t(p) / g)
    p <- p / rowSums(p)
  }
  pen <- 0 # Penalty term for violating model constraints.
  
  # Calculations for heterogeneous groups. Here, 'p' is "delta" (group probabilities prior to group mixing).
  
  # Transform absolute mixing parameter to proportion of all groups (prior to mixing) that each species contributes to heterogeneous groups. 'm' is equivalent to "epsilon" and 'mix.g' to epsilon / (1 + epsilon)". 
  mix.g <- m / (1 + m) 
  
  # Calculate group probability (pi) for each species and heterogeneous groups.
  p <- cbind((p - mix.g), mix.g) / (1 - mix.g)
  
  # For each species, test if any group probabilities 'p' are <0. If so, adjust 'p' >0 and apply penalty to negative log-likelihood. 
  if (any(p < 0)) {
    negative_values <- which(p < 0)
    cat("Warning: Mixing probability of ", m," reduced because", length(negative_values), "probabilities of Psi for individual groups were <0", "\n")
    pen <- pen + (1 + -sum(p[negative_values]) * 10) ^ 2
    negative_rows <- negative_values %% nrow(p)
    negative_rows[which(negative_rows == 0)] <- nrow(p)
    p[negative_rows, ncol(p)] <- 
      p[negative_rows, ncol(p)] + p[negative_values]
    p[negative_values] <- 0
  }
  return(list(p, m, pen))
}

# Function: {group.true.probability.f}: With heterogeneous groups for 2 species, computes probabilities for possible true groups for each unique observation history. Accepts inputs of a matrix of group probabilities 'group_probability', a matrix 'size_probability' with the probability mass for group sizes (1 to the maximum observed group size) of each species, and the observed group size 'size'. Returns a vector with probabilities for each possible true group. 

group.true.probability.f <- function(group_probability, size_probability, size){
  
  # Calculate probabilities first for homogeneous groups at the start/end of each vector, then add probabilities for heterogeneous groups with decreasing numbers of species 1 and increasing numbers of species 2.
  homogeneous <-
    vapply(1:2, function(x)
      group_probability[x] * size_probability[size, x], numeric(length(size)))
  output <- c(homogeneous[1],
              group_probability[3] * size_probability[0:(size - 1), 2] * size_probability[(size - 1):0, 1],
              homogeneous[2])
  tmp <- size:0
  names(output) <- paste0(tmp, ".", size - tmp)
  output
}

# Function: {group.observed.cprobability.f} Computes conditional (on the true group) probabilities for group records with heterogeneous groups. Accepts inputs of vectors with a group record 'd' , classification probabilities 'p', the number of possible combinations (with repetition) of true groups 'states', and group size 's'. Returns vector with a conditional probability for each possible true group.

group.observed.cprobability.f <- function(d, p, states, s){  
  # Return probability = 1 for missing group records.
  if (s < 1) return(rep(1, states))
  
  # Compute probabilities for each combination of true states (possible true groups) and each permutation of observation states in the group record by raising classification probabilities to the exponent of the numbers of individuals in corresponding observation states and then multiplying by the multinomial coefficient for each permutation.
  t1 <- likelihood.equations[[paste0("observed.", paste0(d, collapse = "."))]]
  t2 <- apply(t1[-(1:2)], 1, function(z) prod(p ^ z))
  t2 <- t2 * t1[[2]]
  
  # Compute overall probabilities for possible true states by summing across permutations of observation states in the group record with identical index values.
  vapply(1:max(t1[[1]]), function(x) sum(t2[which(t1[[1]] == x)]), numeric(1))  
}

# Function: {group.probmass.f} Probability mass function for mixture distribution describing the distribution of group sizes. Accepts inputs of a data frame of survey data 'dat' and vector of parameter values 'parameters' for the mixture distribution. Returns probability mass for group size 1 to maximum observed group size.

group.probmass.f <- function(parameters, dat){
  names(parameters) <- NULL
  bernoulli_probability <- parameters[1] # P{group size = 2}, otherwise size = 1
  geometric_probability <- parameters[2] # Probability parameter for geometric distribution
  mixture_weight <- parameters[3] # Weight for Bernoulli relative to geometric distribution
  size <- unique(dat$group_size) - 2 # Group sizes of observations minus groups of sizes 1/2 
  size <- size[size > 0]
  
  # Weighted probability mass of mixture distribution
  c(
    c((1 - bernoulli_probability), bernoulli_probability) * mixture_weight,
    zero.truncated.geometric.probmass.f(geometric_probability, size) * (1 - mixture_weight)
  )
}

# Function: {zero.truncated.geometric.probmass.f} Probability mass function for a zero-truncated geometric distribution. Accepts inputs of vectors with the probability parameter 'prob' of a zero-truncated geometric distribution and with observed group sizes 'size'. Returns the probability mass for group sizes 1 to maximum observed group size.

zero.truncated.geometric.probmass.f <- function(prob, size) {
  if (prob == 0) {
    return(as.numeric(size == 1))
  } else{
    prob * (1 - prob) ^ ((1:max(size)) - 1)
  }
}

###############################################################################
#           Function for murrelet model
###############################################################################

# This function optimizes the model estimating uncertain identification for murrelets described in Appendix S5 in the companion article. It is included in the R workspace 'murrelet_data.RData'

# Function: {murrelet.model.f} Model for estimating uncertain identification for two species of murrelets using line transects employing multi-observer methods. Accepts inputs of a vector parameter values 'parameters' and a data frame of group records for survey data 'dat'. Returns the -log(likelihood) plus a penalty term for violating model constraints. 

murrelet.model.f <- function(parameters, dat){
  
  ## Import and summarize parameters for model structure and estimated parameters ----------
  A <- 3   # Number of observation states
  n_observers <- 4 # Number of observers
  
  # Array 'theta_betas' contains regression parameters for multinomial logit regression predicting uncertain identification probabilities. 
  # Dimension 1 (row) = true species state, dimension 2 (column) = b0 (intercept) and b1 (slope) regression coefficients, and dimension 3 (matrix) = observers (primary and secondary). 
  # Vector 'theta_betas_2' contains the 2nd order slope coefficient predicting partial identification by primary observers. 
  theta_betas <-
    array(
      c(parameters[1:2], parameters[5:6], parameters[3:4],  0, parameters[7]),
      dim = c(2, 2, 2),
      dimnames = list(c(paste0("true_state_", 1:2)),
                      c("b0", "b1"),
                      c("primary", "secondary"))
    )
  
  theta_betas_2 <- parameters[8]
  
  # Vector 'precipitation_betas' contains intercept adjustment parameters for predicting partial identification by primary/secondary observers with precipitation = 1 (rain, mist, or fog). 
  precipitation_betas <- parameters[9:10]
  
  # Vector 'psi_1' contains parameters for true species probabilities for species 1 (Kittlitz's murrelets) in density areas with low and high predicted densities of Kittlitz's murrelets
  psi_1 <- parameters[11:12]
  
  # Heterogeneous group probability (back-transformed from logit scale to probability scale)
  pi_12 <- plogis(parameters[13]) 
  
  # Matrix 'group_distribution_parameters' contains parameters for estimating the distribution of group sizes using a hurdle model. 
  # Dimension 1 (row) = true species states 1 to 2, dimension 2 (column) = Bernoulli probability (tau), geometric probability parameter (q), weight parameter (omega). The geometric parameter q is equal between species. 
  
  group_distribution_parameters <-
    matrix(
      parameters[c(14:17, 15, 18)] ,
      byrow = T,
      nrow = 2,
      dimnames = list(c(paste0("true_state_", 1:2)),
                      c("tau", "q", "omega"))
    )
  
  # Matrix 'distance' contains distance (perpendicular distance from the transect center line in decameters) and distance^2 covariate values for each observation history
  
  distance <- matrix(
    c(dat$perpendicular_distance, dat$perpendicular_distance ^ 2),
    ncol = 2,
    dimnames = list(c(
      paste0("obs_history_", 1:nrow(dat))
    ), c("distance", "distance^2"))
  )
  
  # Array 'theta_arr' contains classification probabilities (arranged as vectors for efficient computation). Multinomial logit link functions enforce that probabilities for each observer and true species state sum to 1. The reference category is correct classification y = z.
  # Dimension 1 (row) = classification probabilities (theta), dimension 2 (column) = observation histories, and dimension 3 (matrix) = observers (primary and 3 secondary). 
  
  theta_arr <- array(0, dim = c(6, nrow(dat), 4), 
                     dimnames = list(
                       c("theta_11",	"theta_12",	"theta_21",	"theta_22",	"theta_31",	"theta_32"),
                       c(paste0("obs_history_", 1:nrow(dat))),
                       c("primary", paste0("secondary_", 1:3))
                     ))
  
  # 'denom' is summation term in denominator of the multinomial logit regression equations
  
  denom <- 1 + exp(theta_betas[1, 1, 1] + theta_betas[1, 2, 1] * distance[, 1]) +
    exp(theta_betas[2, 1, 1] + theta_betas[2, 2, 1] * distance[, 1] + theta_betas_2[1] * distance[, 2] + precipitation_betas[1] * dat$precipitation)
  
  theta_arr[1,  , 1] <- theta_arr[4,  , 1] <- 1 / denom
  theta_arr[3,  , 1] <- theta_arr[2,  , 1] <- 
    exp(theta_betas[1, 1, 1] + theta_betas[1, 2, 1] * distance[, 1]) / denom
  theta_arr[5,  , 1] <- theta_arr[6,  , 1] <- 
    exp(theta_betas[2, 1, 1] + theta_betas[2, 2, 1] * distance[, 1] + theta_betas_2[1] * distance[, 2] + precipitation_betas[1] * dat$precipitation) / denom
  
  denom <- 1 + exp(theta_betas[1, 1, 2] + theta_betas[1, 2, 2] * distance[, 1]) +
    exp(theta_betas[2, 1, 2] + theta_betas[2, 2, 2] * distance[, 1] + precipitation_betas[2] * dat$precipitation)
  
  theta_arr[1,  , 2] <- theta_arr[4,  , 2] <- 1 / denom
  theta_arr[3,  , 2] <- theta_arr[2,  , 2] <- 
    exp(theta_betas[1, 1, 2] + theta_betas[1, 2, 2] * distance[, 1]) / denom
  theta_arr[5,  , 2] <- theta_arr[6,  , 2] <- 
    exp(theta_betas[2, 1, 2] + theta_betas[2, 2, 2] * distance[, 1] + precipitation_betas[2] * dat$precipitation) / denom
  
  # Secondary observers assumed to have identical classification probabilities
  theta_arr[, , 4] <- theta_arr[, , 3]  <- theta_arr[,  , 2]
  
  # Extract mean group size parameters for each species to vector 'group_size_mean'
  
  group_size_mean <- apply(group_distribution_parameters, 1, function(x)
    (1 + x[1]) * x[3] + (2 + 1 / x[2]) * (1 - x[3]))
  
  ## Calculate group probabilities and group size probability mass for each species ----------
  # The matrix 'group_size_probmass' contains probability mass (mu) for each species (column) and group size (row) up to the maximum observed group size.
  
  group_size_probmass <- apply(group_distribution_parameters, 1, group.probmass.f, dat = dat)
  
  rownames(group_size_probmass) <- c(paste0("group_size_", 1:nrow(group_size_probmass)))
  
  # Make matrix of group probabilities (prior to formation of heterogeneous groups) for each observation history
  group_probability <- mlogit.regress.predict.f(dat$density_area, psi_1, 2)
  
  output <- group.probability.psi.constant.f(group_probability, pi_12, group_size_mean)
  
  # The vector 'group_probability' gives the group probabilities (pi) for occurrence of groups of each species (pi.1, pi.2) and heterogeneous groups (pi.12).
  group_probability <- output[[1]]
  dimnames(group_probability) <- list(c(paste0("obs_history_", 1:nrow(dat))), c("pi.1", "pi.2", "pi.12"))
  
  # Enforce constraint the group probabilities for species 1 and 2 don't exceed 1. The 'penalty' term added to the -log(likelihood) scales with the magnitude of violation of the constraint. 
  penalty <- output[[3]] 
  
  ## Compute -log{Likelihood} (nLL) for data conditional on estimated parameter values ----------
  # Vector 'likelihood' contains probabilities for each unique observation history
  likelihood <- numeric(nrow(dat))
  n_group_size <-  dat  %>% count(group_size) # Sample of observation histories by group size
  
  # Loop for each observed group size
  for (i in unique(dat$group_size)) {
    n_histories <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] # Sample of observation histories for current group size
    rows_i <- which(dat$group_size == i) # Row numbers for current group size
    
    # Matrix 'group_true_probability' contains probabilities of each possible true group (rows, labeled with counts of species1.species2) for each unique observation history (columns).
    group_true_probability <-
      apply(group_probability[rows_i, , drop = F], 1, group.true.probability.f, group_size_probmass, i)
    
    B_states <- nrow(group_true_probability) # Number of possible true groups
    
    # For the current group size, matrix 'likelihood_i' contains probabilities of observed groups (conditional on possible true groups) for each observer (columns 1 to 4) and probabilities of true groups (column 5). 
    likelihood_i <-
      as_tibble(matrix(
        c(
          rep(numeric(1), n_histories * B_states * n_observers),
          group_true_probability
        ), 
        ncol =  n_observers + 1,
        dimnames = list(c(), c("obs_p", paste0("obs_s", 1:3), "psi"))))
    
    for (obs in 1:n_observers) {
      # Matrix 'dat.tmp contains observed groups for the current observer 'obs' and group size
      dat_tmp <- filter(dat, group_size == i) %>%
        select(any_of((((obs - 1) * A) + 1):(obs * A)))
      size_tmp <- rowSums(dat_tmp)
      likelihood_i[, obs] <- 
        as.vector(vapply(1:n_histories, function(x) 
          group.observed.cprobability.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B_states, size_tmp[x]),
          numeric(B_states))) 
    }
    
    # Integers in 'index' associate rows in 'probability.mat' with sequentially numbered unique observation histories with the current group size.  
    index <- rep(1:n_histories, each = B_states)
    
    # Add likelihoods for the current group sizes to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of probabilities in each row of 'likelihood_i' summed across all possible combinations of true groups for each unique observation history specified by 'index'.
    
    likelihood[rows_i] <-
      transmute(likelihood_i,
                index = index,
                b_product = apply(likelihood_i, 1, function(x) prod(x))) %>%
      group_by(index)  %>%
      summarise(likelihood = sum(b_product)) %>%
      select(likelihood) %>%
      unlist(.)
  } # End of loop for group size
  
  # Compute the -log{likelihood} as the sum of the probabilities for each observation history. Add penalty term for violating constraint(s).
  sum(-log(likelihood)) + penalty
}

###############################################################################
#           Code for executing statistical analyses
###############################################################################

# Parameter names follow conventions defined in Metadata S4. Parameters are described in Appendix S5 in the companion article

## Specify initial parameter values ----------

# Specify initial parameter values for each parameter by altering the values assigned below. All initial values must fall within the constraints specified below. 

# Initial values: Intercept coefficients (beta 0) of multinomial logit regression models predicting classification probabilities. Enter as probabilities. These will be transformed (below) to log odds ratios for multinomial logistic regression analyses. 
b0_distance_p <- matrix(c(0.03, 0.02),
               byrow = T, ncol = 2)
b0_distance_s <- matrix(c(0.01, 0.04),
               byrow = T, ncol = 2)

# Initial values: Slope coefficients (beta 1, beta 2) of multinomial logit regression models predicting classification probabilities 
b1_distance_p <- c(1, 0)
b1_distance_s <- c(0.5)
b2_distance_p <- c(0.5)
b1_precipitation <- c(0.2, 0.3)

# Initial values: True species probabilities for areas of low and high expected densities of Kittlitz's murrelets. The first value is true species probability (which will be logit-transformed), the second value is an intercept adjustment term to the first term.
psi_1 <- c(0.4, 0)

# Initial vales: Heterogeneous group probability. Enter as a probability.
pi_12 <- 0.05

# Initial values: Hurdle model of group size distribution
# Row 1 = Kittlitz's murrelet, row 2 = marbled murrelet
# Column 1 = weight parameter (omega), column 2 = Bernoulli probability (tau), column 3 = geometric distribution probability parameter (q)
# The geometric distribution probability parameter (q) is assumed equal between species

group_size_distribution <-  c(0.5, 0.2, 0.7,
                              0.5,      0.7)

# Transform classification probabilities to log odds ratios
tmp_p <-
  as.vector(t(log(b0_distance_p / (
    1 - rowSums(b0_distance_p)
  ))))
tmp_s <-
  as.vector(t(log(b0_distance_s / (
    1 - rowSums(b0_distance_s)
  ))))

# Combine initial parameter values into one file and add parameter names
parameters_ini <- c(tmp_p, tmp_s, b1_distance_p, b1_distance_s, b2_distance_p, b1_precipitation, qlogis(psi_1[1]), psi_1[2], qlogis(pi_12), group_size_distribution)
parameter_names <- 
  c("b0_distance_p_21", "b0_distance_p_31", "b0_distance_s_21", "b0_distance_s_31", "b1_distance_p_21", "b1_distance_p_31", "b1_distance_s_21", "b2_distance_p_31", 
    "b1_precipitation_p", "b1_precipitation_s", "psi_1_low", "psi_1_high",  "pi_12", "tau_1", "q", "omega_1", "tau_2", "omega_2")
names(parameters_ini) <- parameter_names

## Specify lower and upper box constraints for each parameter ----------

# Specify lower and upper constraints for values of each each parameter by altering the values assigned below. 

# Constraints: Heterogeneous group probabilities (pi_12) and true species probabilities (psi_1) as probabilities
constraint_pi_12 <- qlogis(c(0.000001, 0.15))
constraint_psi_1 <- qlogis(c(0.0001, 0.99))

# Constraints: Intercept and slope coefficients  of multinomial logistic regressions predicting classification probabilities
# Enter values for intercept parameters (beta 0) of multinomial logit regression models as probabilities. Enter other values (beta 1, beta 2) as slope coefficients.
constraint_b0_distance <- c(0.0005, 0.75)
constraint_b1_distance <- c(-10, 10)
constraint_b2_distance <- c(-15, 15)
constraint_b1_precipitation <- c(-5, 5)

# Constraints: Hurdle model of group size distribution
# Column 1 = weight parameter (omega), column 2 = Bernoulli probability (tau), column 3 = geometric distribution probability parameter (q). Rows give lower and upper constraints.
constraint_group_size_distribution <- matrix(c(0.01, 0.01, 0.01,
                                             0.999, 0.9, 0.99),
                                           byrow = T,
                                           nrow = 2)

# Combine lower and upper box constraints into separate files and add parameter names
constraint_low <-
  c(
    rep(qlogis(constraint_b0_distance[1]), 4),
    rep(constraint_b1_distance[1], 3),
    rep(constraint_b2_distance[1], 1),
    rep(constraint_b1_precipitation[1], 2),
    rep(constraint_psi_1[1], 2),
    constraint_pi_12[1],
    constraint_group_size_distribution[1, c(1:3, 1, 3)]
  )

constraint_up <-
  c(
    rep(constraint_b0_distance[2], 4),
    rep(constraint_b1_distance[2], 3),
    rep(constraint_b2_distance[2], 1),
    rep(constraint_b1_precipitation[2], 2),
    rep(constraint_psi_1[2], 2),
    constraint_pi_12[2],
    constraint_group_size_distribution[2, c(1:3, 1, 3)]
  )

names(constraint_low) <- names(constraint_up) <- parameter_names

## Execute statistical analyses ----------

# The 'optim' function solves for parameter values minimizing the -log{likelihood} (plus any penalty terms). Vectors 'constraint_low' and 'constraint_up' specify lower and upper box constraints. 
# Optimization will likely allow enough time to make a cup of tea

model <- 
  optim(parameters_ini, murrelet.model.f, gr = NULL, 
        dat = murrelet_data, 
        hessian = T, 
        method = c("L-BFGS-B"),
        lower = constraint_low,
        upper = constraint_up,
        control = list(trace = 3) 
  )
