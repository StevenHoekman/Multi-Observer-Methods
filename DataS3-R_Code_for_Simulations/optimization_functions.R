# optimization_functions.R
# Functions for likelihood optimization of models estimating uncertain identification using multi-observer methods, version 1.1.2
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# R computer code for optimizing multi-observation method (MOM) and single-observation method (SOM) models for estimating uncertain species identification by minimizing the -log(likelihood). These functions are designed to conduct simulation analyses described in the companion article and Appendices S1 to S3. Each function optimizes models with differing predictive covariates, as described in comments of each function and in MetadataS3.pdf. Code developed and tested in R version 4.1.

# This code should be executed prior to conducting simulation analyses in 'simulations_R_code.R'

###############################################################################
#               Function Inputs/Outputs 
###############################################################################

# These functions accept inputs of:
# (1) Initial parameter values 'param'
# (2) Formatted simulated data for a single survey 'dat'
# (3) Optional tables of key values for unique observed groups 'keys' or unique binned covariate values 'keys_psi'
# (4) A simulation profile defining simulation and model inputs 'sim_profile' 

# Likelihoods are computed using keyed tables (if present) or otherwise using unique observation histories 
# User-defined penalty functions constrain optimization to avoid unrealistic or inadmissible parameter values

# These functions produce output of:
# (1) Parameter estimates
# (2) The -log(likelihood) plus any penalty terms
# (3) Codes indicating success of optimization
# (4) The Hessian matrix

# These functions require an external input file containing a list object 'likelihood_equations' with pre-computed values for likelihood computations. See DataS2. 

###############################################################################
#          Functions for likelihood optimization
###############################################################################

# Function: {optimize.M.f} Optimization of MOM/SOM models with misidentification and partial identification. Supports <=4 true species states B and observation states A. Accepts inputs of initial parameter values 'param', simulated survey data 'dat', optional key values for unique observed groups 'keys', and the simulation profile 'sim_profile'. Returns the -log(likelihood) plus any penalty term(s). 

optimize.M.f <- function(param, dat, keys, sim_profile){

## ----- Import and summarize true and estimated parameters -----
  
  # Summarize group sizes
  g <- unique(dat$group_size) # Observed group sizes
  g1 <- sum(dat$group_size == 1) # Number of observation histories with group size = 1
  
  # Extract true simulation parameters from 'sim_profile'
  B <- sim_profile$B; A <- sim_profile$A # Number of true species states (B) and observation states (A)
  n_O_ps <- c(sim_profile$O_p, sim_profile$O_s) # Number of primary/secondary observers
  n_observers <- sum(n_O_ps) # Total # of observers
  het_true <- sim_profile[grep("mix", names(sim_profile))] # True heterogeneous group parameter(s)
  mx_model <- sim_profile[grep("mx_model", names(sim_profile))] # Heterogeneous group model

  # Summarize values of estimated parameters in 'param'
  
  # Extract classification probabilities (theta) to 'theta_tmp' array. Dimension 1 = mis- and partial ID parameters, dimension 2 = true species states 1 to B, dimension 3 = individual observers.
  theta_col <- grep("theta", names(sim_profile))
  
  # If any observers have certain identification (probability of uncertain ID = 0), remove from estimated parameters
  if (any(sim_profile[theta_col] == 0)) {
    theta_col <- theta_col[-which(sim_profile[theta_col] == 0)]
    theta_tmp <- array(c(plogis(param[1:length(theta_col)]), rep(0L, B * (A - 1)) ), c(A - 1, B, n_observers))
    
  # With distinct secondary observers, record classification probabilities of each in separate arrays
  }else if (length(theta_col) == B * (A - 1) * n_observers) {
    theta_tmp <- array(plogis(param[1:length(theta_col)]), c(A - 1, B, n_observers))
    
  # With identical secondary observers, apply classification probabilities for 1st to all
  }else{
    theta_mat <- matrix(plogis(param[1:length(theta_col)]), ncol = B * (A - 1), byrow = TRUE)
    theta_tmp <- array(0, c(A - 1, B, n_observers))
    for (i in 1:n_observers) {
      if (dim(theta_mat)[1] >= i) {
        theta_tmp[, , i] <- theta_mat[i, ]
      }else{
        theta_tmp[, , i] <- theta_tmp[, , (i - 1)]
      }
    }
  }
  
  # True species probabilities for species 1 to (B - 1)
  psi <- plogis(param[(length(theta_col) + 1):(length(theta_col) + B - 1)])
  
  # Extract value(s) for heterogeneous group parameters to vector 'mix'
  if (any(het_true > 0)) {
    mix <- plogis(param[(length(param) - length(het_true) + 1):length(param)])
  }else{
    mix <- rep(0, length(het_true))
  }
  
  # Enforce constraints for probabilities that must sum to 1. The 'penalty' term added to the -log(likelihood) scales with the magnitude of violations of constraints. 
  penalty <- 0
  
  if (sum(psi) > 1) {
    penalty <- penalty.f(sum(psi))
  psi <- sum1.f(psi)
  }
  
  if (any(colSums(theta_tmp) > 1)) {
    sums <- colSums(theta_tmp)
    penalty <- penalty + penalty.f(sums[sums > 1])
    col <- which(colSums(theta_tmp) > 1, arr.ind = TRUE)
    for (i in 1:dim(col)[1]) {
      theta_tmp[, col[i, 1], col[i, 2]] <- sum1.f(theta_tmp[, col[i, 1], col[i, 2]])
    }
  }
  
  # True species probabilities for species 1 to B
  psi <- c(psi, 1 - sum(psi))

  # With group sizes >=1, extract group size parameters for species 1 to B to vector 'g_spp'
  g_spp <- rep(1, B) 
  if (any(g > 1)) {
    g_spp <- param[(length(param) - n_parameters[5] - B + 1):(length(param) - n_parameters[5])]
  }
  
  # Matrix 'theta_mat' summarizes classification probabilities (theta) of each observer as a vector for efficient computation
  # row = observers, col = classification probabilities 1 to A in sequence for species 1 to B
  theta_mat <- t(vapply(1:n_observers, theta4.f, numeric(B * A), theta_tmp))
  
  # 'theta_arr' summarizes classification probabilities (theta)  as an array for computing probabilities for group size = 1.
  # dimension 1 (row) = true species state 1 to B, dimension 2 (col) = classification probabilities for observation states 1 to A, dimension 3 (matrix) = observer
  theta_arr <- array(t(theta_mat), c(B, A, n_observers))
  
  # Vector 'lambda' contains lambda parameters (defined by mean group size) of zero-truncated Poisson distribution for species 1 to B
  lambda <- rep(0, B)
  lambda[which(g_spp > 1)] <-
    vapply(g_spp[which(g_spp > 1)], function(x)
      optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))

## ----- Compute group probabilities and group size probability mass for each species -----
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each spp (column) and group size (row) up to maximum observed group size
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = B)
  
  # Compute group probabilities for the specified heterogeneous group model
  if (mx_model == "constant") {
    # Function for "constant" heterogeneous group model
    output <- group.probability.constant.f(psi, mix, g_spp)
  } else if (mx_model == "encounter") {
    # Function for "encounter" heterogeneous group model
    output <- group.probability.encounter.f(psi, mix, g_spp)
  }

  # Vector 'group_probability' gives group probabilities (pi) for homogeneous (true states 1 to B) and then heterogeneous groups (if present)
  group_probability <- output[[1]]
  # Penalty term added to the -log(likelihood) if heterogeneous group parameters take inadmissible values 
  penalty <- penalty + output[[3]]

## ----- Compute -log(likelihood) for data conditional on estimated parameter values -----
  
  # Vector 'likelihood' contains likelihoods for each unique observation history
  likelihood <- numeric(dim(dat)[1])

## ----- Compute likelihoods: group size = 1 -----
  if (g1 > 0) { 
    # Array 'group1_obs_cprobability' gives probabilities of classifications (conditional on B possible true species states) for each observer for each unique observation history
    # Dimension 1 = unique observation histories, dimension 2 = true species state 1 to B, dimension 3 = observer
    
    group1_obs_cprobability <- array(0, dim = c(g1, B, n_observers))
    for (obs in 1:n_observers) {
      group1_obs_cprobability[, , obs] <-
        apply(theta_arr[, , obs], 1, function(x)
          dmnom(dat[(1:g1), (((obs - 1) * A) + 1):(obs * A)], size = 1, prob = matrix(x, nrow = 1)))
    }
    
    group1_obs_cprobability[group1_obs_cprobability == 0] <- 1
    group1_obs_cprobability[is.na(group1_obs_cprobability)] <- 1
    
    # For each true species state 1 to B (columns), matrix 'probability_mat' contains probabilities of each unique observation history (row) computed as the product of probabilities of classifications for all observers and probabilities of true species states
    
    probability_mat <- matrix(0, g1, B)
    for (b in 1:B) {
      probability_mat[, b] <- apply(group1_obs_cprobability[, b, , drop = FALSE], 1, function(x) 
        group_probability[b] * group_size_probmass[1, b] * prod(x))
    }
    # Likelihoods for each unique observation history are the summed probabilities for each true species state
    likelihood[1:g1] <- rowSums(probability_mat)
  }

## ----- Compute likelihoods: group size >= 1 and homogeneous groups -----
    if (any(g > 1 & all(het_true == 0))) {
      # 'n_group_size' = count of unique observation histories by group size
      n_group_size <-  dat  %>% count(group_size)
      
      # If key table of unique observed groups is present, compute likelihoods from keyed table of probabilities
      if (!is.null(keys)) {
        key_col <- grep("key", names(dat))

        # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level). Matrices for each observer contain probabilities of the keyed observed group (columns) conditional on each possible true group (rows).
        
        group_observed_cprobability_key <- lapply(1:n_observers, function(x)
          apply(keys, 1, function(y)
            group.observed.cprobability.homogeneous.f(y[1:A] , theta_mat[x, ],  B, y[A + 3])
          )
        )
        
        # 'dat_tmp' = observation histories for group size >1
        dat_tmp <- filter(dat, group_size > 1)
        n_group_size_tmp <- which(n_group_size$group_size > 1)
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories.  
        index <- tibble(index = unlist(lapply(n_group_size_tmp, function(x) 
          rep(1:n_group_size[[x, 2]], each = B) + sum(n_group_size[c(n_group_size_tmp[1]:x) , 2]) - n_group_size[[x, 2]] )))
        
        # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Row correspond to unique combinations of possible true groups and unique observation histories.
        
          probability_mat <- matrix(
            c(
              sapply(1:n_observers, function(x)
                unlist(group_observed_cprobability_key[[x]][, unlist(dat_tmp[, key_col[[x]]])])),
              unlist(apply(n_group_size[n_group_size_tmp, ], 1, function(x)
                rep(group_size_probmass[x[1], ] * group_probability, x[2])))
            ), 
          ncol = n_observers + 1)
        
        # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are computed as the product of each row in 'probability_mat' summed across all possible combinations of true groups for that unique observation history. 
          
        likelihood[pmatch(0, likelihood):length(likelihood)] <-
          bind_cols(index, data.frame(product = vapply(1:dim(probability_mat)[1], function(x)
            prod(probability_mat[x, ]), numeric(1)))) %>%
          group_by(index) %>%
          summarise(likelihood = sum(product)) %>%
          select(likelihood) %>%
          unlist(.)
      }else{
    ## If key table of unique observed groups is NOT present, compute likelihoods from individual unique observation histories
        
        # Loop for observed group sizes > 1
        for (i in g[g > 1]) {
          # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
          n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]]
          rows_i <- which(dat$group_size == i) 
          
          # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
          
          probability_mat <-
            as_tibble(matrix(
              c(
                rep(numeric(1), n_i * 2 * n_observers),
                rep(group_size_probmass[i, ] * group_probability, n_i)
              ), 
              ncol =  n_observers + 1),
              .name_repair = NULL)
          
          # Add conditional probabilities of observed groups for each observer
          
          for (obs in 1:n_observers) {
            dat_tmp <- dat[rows_i, ] %>%
              select(all_of((((obs - 1) * A) + 1):(obs * A)))
            record_sum <- rowSums(dat_tmp)
            probability_mat[, obs] <- 
              as.vector(vapply(1:n_i, function(x) 
                group.observed.cprobability.homogeneous.f(dat_tmp[x, ], theta_mat[obs, ], 2, record_sum[x]),
                numeric(2))) 
          }
          
          # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
          index <- rep(1:n_i, each = 2)
          
          # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of each row in 'probability_mat' summed across all possible combinations of true groups for that unique observation history. 

          likelihood[rows_i] <-
            transmute(probability_mat,
                      index = index,
                      b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
            group_by(index)  %>%
            summarise(likelihood = sum(b.product)) %>%
            select(likelihood) %>%
            unlist(.)
        }
      }
    }
  
## ----- Compute likelihoods:group size >= 1 and heterogeneous groups -----
  
  if (any(g > 1 & any(het_true > 0))) {

    # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each spp 1 to B) giving probabilities for each possible true group
    group_true_probability <- group.true.probability.key.f(group_probability, group_size_probmass, g)
    
    # n_group_size = sample of observation histories for each group size
    n_group_size <-  dat  %>% count(group_size)
    
    # If key table of unique observed groups is present, compute likelihoods from keyed table of probabilities
    if (!is.null(keys)) {
      key_col <- grep("key", names(dat))
      n_group_size_tmp <- which(n_group_size$group_size > 1)
      
      # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level) and key value (2nd list level), with vectors for each key value giving probabilities of observed groups for each possible true group
      
      group_observed_cprobability_key <- lapply(1:n_observers, function(x)
        apply(keys, 1, function(y)
          group.observed.cprobability.f(y[1:A], theta_mat[x, ],  y[A + 2], y[A + 3])
        )
      )
      
      # 'B_states' is the # of possible combinations of true species states for true groups of each observed group size > 1

      B_states <-
        tibble::enframe(c(2, vapply(n_group_size[n_group_size_tmp, 1][[1]], function(x)
          length(group_true_probability[[paste0(x)]]), numeric(1))), name = NULL)
      names(B_states) <- "B_states"
      n_group_size <- bind_cols(n_group_size, B_states)
      
      # 'dat_tmp' is observation histories for group sizes >1
      dat_tmp <- filter(dat, group_size > 1)
      n_group_size_tmp <- which(n_group_size$group_size > 1)
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      
      index <- tibble(index = unlist(lapply(n_group_size_tmp, function(x)
        rep(1:n_group_size[[x, 2]], each = n_group_size[[x, 3]]) + sum(n_group_size[c(n_group_size_tmp[1]:x) , 2]) - n_group_size[[x, 2]])))
      
      # Calculate likelihoods if no observers are missing observations

      # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
      
      if (min(dat_tmp[, unlist(key_col)]) > 1) {
        probability_mat <- matrix(
          c(
            sapply(1:n_observers, function(x)
              unlist(group_observed_cprobability_key[[x]][unlist(dat_tmp[, key_col[[x]]])])),
            unlist(apply(n_group_size[n_group_size_tmp,], 1, function(x)
              rep(group_true_probability[[paste0(x[1])]], x[2])))
          ), 
          ncol = n_observers + 1)
      }else{
        # Calculate likelihoods if any secondary observers are missing observations

        # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories
        
        probability_mat <- matrix(
          c(
            sapply(1:n_observers, function(x)
              unlist(apply(bind_cols(dat_tmp['group_size'], dat_tmp[key_col[[x]]]), 1,
                           probability.key.f, x, n_group_size, group_observed_cprobability_key))),
            unlist(apply(n_group_size[n_group_size_tmp,], 1, function(x)
              rep(group_true_probability[[paste0(x[1])]], x[2])))
          ),
          ncol = n_observers + 1)
      }
      
      # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of each row in 'probability_mat' summed across all possible combinations of true groups for that unique observation history. 
      
      likelihood[pmatch(0, likelihood):length(likelihood)] <-
        bind_cols(index, as_tibble(vapply(1:dim(probability_mat)[1], function(x)
          prod(probability_mat[x,]), numeric(1)))) %>%
        group_by(index) %>%
        summarise(likelihood = sum(value)) %>%
        select(likelihood) %>%
        unlist(.)
    }else{
      ## If key table of unique observed groups is NOT present, compute likelihoods from individual unique observation histories
      
      # Loop for observed group sizes > 1
      for (i in g[g > 1]) {
        # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
        n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]]
        rows_i <- which(dat$group_size == i) 
        B_states <- length(group_true_probability[[paste(i)]]) # Number of combinations of true states for the current group size
        
        # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
        
        probability_mat <-
          as_tibble(matrix(
            c(
              rep(numeric(1), n_i * B_states * n_observers),
              rep(group_true_probability[[paste(i)]], n_i)
            ), 
            ncol =  n_observers + 1),
            .name_repair = NULL)
        
        # Add condition probabilities of observed groups for each observer
        
        for (obs in 1:n_observers) {
          dat_tmp <- dat[rows_i, ] %>%
            select(all_of((((obs - 1) * A) + 1):(obs * A)))
          record_sum <- rowSums(dat_tmp)
          probability_mat[, obs] <- 
            as.vector(vapply(1:n_i, function(x) 
              group.observed.cprobability.f(dat_tmp[x, ], theta_mat[obs, ], B_states, record_sum[x]),
              numeric(B_states))) 
        }
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
        index <- rep(1:n_i, each = B_states)
        
        # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of each row in 'probability_mat' summed across all possible combinations of true groups for that unique observation history. 
        
        likelihood[rows_i] <-
          transmute(probability_mat,
                    index = index,
                    b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
          group_by(index)  %>%
          summarise(likelihood = sum(b.product)) %>%
          select(likelihood) %>%
          unlist(.)
      }
    }
  }
  
  ## Compute -log(likelihood) as the product of likelihoods for each observation history and the count for each history. Add penalty term(s) accrued for violating model constraints.
  sum(dat[, "count"] * -log(likelihood)) + penalty
}

# Function: {optimize.M.theta.f} Optimization of MOM/SOM models with misidentification, partial identification, and a covariate predicting classification probabilities (theta). Supports <=3 true species states B and observation states A. Accepts inputs of initial parameter values 'param', simulated survey data 'dat', optional key values for unique observed groups 'keys', and the simulation profile 'sim_profile'. Returns the -log(likelihood) plus any penalty term(s).

optimize.M.theta.f <- function(param, dat, keys, sim_profile){
  
## ----- Import and summarize true and estimated parameters -----
  
  # Summarise group sizes
  g <- unique(dat$group_size) # Observed group sizes
  
  # Extract true simulation parameters from 'sim_profile'
  B <- sim_profile$B; A <- sim_profile$A # Number of true spp states (B) and observation states (A)
  n_O_ps <- c(sim_profile$O_p, sim_profile$O_s) # Number of primary/secondary observers
  n_observers <- sum(n_O_ps) # Total number of observers
  het_true <- sim_profile[grep("mix", names(sim_profile))] # True heterogeneous group parameter(s)
  mx_model <- sim_profile[grep("mx_model", names(sim_profile))] # Heterogeneous group  model
 
  # Summarize values of estimated parameters in 'param'
  
  # Extract value(s) for heterogeneous group parameters to vector 'mix'
  if (any(het_true > 0)) {
    mix <- plogis(param[(length(param) - length(het_true) + 1):length(param)])
  }else{
    mix <- rep(0, length(het_true))
  }
  
  # With key table use unique combinations of observed groups and covariate values for likelihood computations. Without, use unique combinations of observation histories and covariate values.
  if (is.null(keys)) {
    n_unique <- dim(dat)[1] # Number of unique observation histories
    theta_cov <- rbind(dat$covariate_theta, dat$covariate_theta) # Covariate values for unique observation histories
    if (sim_profile$Model == 'M.theta.ps') {
      theta_cov[2, ] <- dat$covariate_theta_s
    }
  }else{
    n_unique <- dim(keys)[1] # Number of keyed observed groups
    theta_cov <- rbind(keys$covariate_theta, keys$covariate_theta) # Covariate values for keyed groups

  }

# Array 'theta_arr' summarizes classification probabilities (theta) of observers as vectors for efficient computation. Vectors are stored in an array with Dimension 1 (row) = classification probabilities 1 to A for species 1 to B, dimension 2 (column) = unique observation histories or keyed observed groups, and dimension 3 (matrix) = observer. Multinomial logit link functions enforce that each column (dimension 2) of each array sums to 1. Baseline category is correct classification y = z. Regression coefficients predicting classification probabilities based on group-level covariates are summarized in the 'betas' arrays.

# Vector 'psi' contains estimated true species probabilities (psi) for species 1 to (B - 1)

  if (B == 2) {
    if (A == 2) {
      # True species states (B) = observation states (A) = 2

      # Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 (row) = true species states 1 to B, dimension 2 (col) = regression coefficients (intercepts in column 1, slopes in column 2), dimension 3 (matrix) = observer. 
      
      if (n_O_ps[1] > 0) {
        betas_arr <- array(c(param[1:2], param[5:6], param[3:4], param[7:8])
                           , dim = c(2, 2, 2))
        psi <- plogis(param[9])
      }else{
        betas_arr <- array(c(param[1:4])
                           , dim = c(2, 2, 2))
        psi <- plogis(param[5])
      }
      
      theta_arr <- array(0, dim = c(4, n_unique, 5))
      for (obs in 1:2) {
        theta_arr[c(3, 1), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[1, , obs, drop = FALSE], A))
        theta_arr[c(2, 4), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[2, , obs, drop = FALSE], A))
      }

    }else {
    # True species states (B) = 2, observation states (A) = 3 (partial identification)
      
    # Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 = observation states 1 to A, dimension 2 = regression coefficients (intercepts in column 1, slopes in column 2), dimension 3 = observer (primary = 1/2, secondary = 3/4) and true species state (state 1 = 1/3, state 2 = 2/4).
      
    theta_arr <- array(0, dim = c(6, n_unique, 5))

    if (n_O_ps[1] > 0) {
      betas_arr <- array(c(param[1:2], param[9:10], param[3:4], param[11:12],
                           param[5:6], param[13:14], param[7:8], param[15:16])
                         , dim = c(2, 2, 4))
      psi <- plogis(param[17])
    }else{
      betas_arr <- array(c(param[1:2], param[5:6], param[3:4], param[7:8])
                         , dim = c(2, 2, 4))
      psi <- plogis(param[9])
    }
    
    for (obs in 1:2) {
      theta_arr[c(3, 5, 1), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[, , obs + (obs - 1), drop = FALSE], A))
      theta_arr[c(2, 6, 4), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[, , obs * 2, drop = FALSE], A))
    }
  }
}else if (B == 3 & A == 3) {
  # True species states (B) = observation states (A) = 3
  
  # 'Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 = observation states 1 to A, dimension 2 = regression coefficients (intercepts in column 1, slopes in column 2), dimension 3 = true species state 1 to B, dimension 4 = observer.
  
  theta_arr <- array(0, dim = c(9, n_unique, 5))
  
  if (n_O_ps[1] > 0) {
    betas_arr <- array(c(param[1:2], param[13:14], param[3:4], param[15:16],
                         param[5:6], param[17:18], param[7:8], param[19:20],
                         param[9:10], param[21:22], param[11:12], param[23:24])
                       , dim = c(2, 2, 3, 2))
    psi <- plogis(param[25:26])
  }else{
    betas_arr <- array(c(param[1:2], param[7:8], param[3:4], param[9:10],
                         param[5:6], param[11:12])
                       , dim = c(2, 2, 3, 2))
    psi <- plogis(param[13:14])
  }
  
  for (obs in 1:2) {
    theta_arr[c(4, 7, 1), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[, , 1, obs, drop = FALSE], A))
    theta_arr[c(2, 8, 5), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[, , 2, obs, drop = FALSE], A))
    theta_arr[c(3, 6, 9), , obs] <- t(mlogit.regress.predict.f(theta_cov[obs, ], betas_arr[, , 3, obs, drop = FALSE], A))
  }
}
  
  # Secondary observers assumed to have identical classification probabilities
  theta_arr[, , 5] <- theta_arr[, , 4] <- theta_arr[, , 3]  <- theta_arr[,  , 2]
  
  # Enforce constraint that true species probability parameters sum to 1. The 'penalty' term added to the -log(likelihood) scales with the magnitude of violations of constraints.
  penalty <- 0

  if (sum(psi) > 1) {
    penalty <- penalty.f(sum(psi))
    psi <- sum1.f(psi)
  }
  
  # If any groups >1, extract mean group size parameters for species 1 to B to vector 'g_spp'
  g_spp <- rep(1, B) 
  if (any(g > 1)) {
    g_spp <- param[(sum(n_parameters[1:3]) + 1):sum(n_parameters[1:4])]
  }
  
  # True species probabilities for species 1 to B
  psi <- c(psi, 1 - sum(psi))
  
  # Vector 'lambda' contains lambda parameters (defined by mean group size) of zero-truncated Poisson distribution for species 1 to B
  lambda <- rep(0, B)
  lambda[which(g_spp > 1)] <- vapply(g_spp[which(g_spp > 1)], function(x)
    optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))
  
## ----- Compute group probabilities and group size probability mass for each species -----
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each spp (column) and group size (row) up to maximum observed group size
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = B)
  
  # Compute group probabilities
  if (mx_model == "constant") {
    # Function for "constant" heterogeneous group model and for homogeneous groups
    output <- group.probability.constant.f(psi, mix, g_spp)
  }
  
  # Vector 'group_probability' gives group probabilities (pi) for homogeneous (true states 1 to B) and then heterogeneous groups (if present)
  group_probability <- output[[1]]
  
  # Penalty term added to the -log(likelihood) if heterogeneous group parameters take inadmissible values 
  penalty <- penalty + output[[3]]
  
## ----- Compute -log(likelihood) for data conditional on estimated parameter values -----
  
  # Vector 'likelihood' contains likelihoods for each unique observation history
  likelihood <- numeric(dim(dat)[1])
  # 'n_group_size' = count of unique observation histories by group size
  n_group_size <-  dat  %>% count(group_size) 
  
  if (all(het_true == 0)) { ## ----- Compute likelihoods: homogeneous groups -----

    # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each spp 1 to B) giving probabilities for each possible true group

    group_true_probability <-
      lapply(1:dim(group_size_probmass)[1], function(x)
        group_size_probmass[x, ] * matrix(group_probability, nrow = 1))
    names(group_true_probability) <- paste(1:max(g))
    
    # If key table of unique observed groups is present, compute likelihoods from keyed table of probabilities
    if (!is.null(keys)) {
      key_col <- grep("key", names(dat))
      
      # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level). Matrices for each observer contain probabilities of the keyed observed group (columns) conditional on each possible true group (rows).
      
      group_observed_cprobability_key <- lapply(1:n_observers, function(x)
        apply(keys, 1, function(y)
          group.observed.cprobability.homogeneous.f(y[1:A] , theta_arr[ ,  y[A + 2], x],  B, y[A + 4])
        )
      )
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- tibble(index = unlist(lapply(1:dim(n_group_size)[1], function(x)
        rep(1:n_group_size[[x, 2]], each = B) + sum(n_group_size[n_group_size[[1, 1]]:n_group_size[[x, 1]] , 2]) - n_group_size[[x, 2]] ))
      )
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.

      probability_mat <- matrix(
        c(
          sapply(1:n_observers, function(x)
            unlist(group_observed_cprobability_key[[x]][, unlist(dat[, key_col[[x]]])])),
          unlist(apply(n_group_size, 1, function(x)
            rep(group_true_probability[[x[1]]], x[2])))
        ),
        ncol = n_observers + 1)
 
      # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups. 
      
      likelihood <-
        bind_cols(index, data.frame(product = vapply(1:dim(probability_mat)[1], function(x)
          prod(probability_mat[x, ]), numeric(1)))) %>%
        group_by(index) %>%
        summarise(likelihood = sum(product)) %>%
        select(likelihood) %>%
        unlist(.)
    }else{
      ## If key table of unique observed groups is NOT present, compute likelihoods from individual unique observation histories
      
      # Loop for observed group sizes
      for (i in g) {
        # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
        n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
        rows_i <- which(dat$group_size == i) 
        
        # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
        
        probability_mat <-
          as_tibble(matrix(
            c(
              rep(numeric(1), n_i * B * n_observers),
              rep(group_true_probability[[i]], n_i)
            ), 
            ncol =  n_observers + 1,
            dimnames = list(c(), c(paste0("obs", 1:n_observers), "group"))))

        # Add condition probabilities of observed groups for each observer
        
        for (obs in 1:n_observers) {
          dat_tmp <- dat[rows_i, ] %>%
            select(all_of((((obs - 1) * A) + 1):(obs * A)))
          record_sum <- rowSums(dat_tmp)
          probability_mat[, obs] <- 
            as.vector(vapply(1:n_i, function(x) 
              group.observed.cprobability.homogeneous.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B, record_sum[x]),
              numeric(B))) 
        }
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
        index <- rep(1:n_i, each = B)
        
        # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
        
        likelihood[rows_i] <-
          transmute(probability_mat,
                    index = index,
                    b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
          group_by(index)  %>%
          summarise(likelihood = sum(b.product)) %>%
          select(likelihood) %>%
          unlist(.)
      }
    }
    

  } else if (!is.null(keys)) { ## ----- Compute likelihoods:heterogeneous groups, WITH key table -----

    
    key_col <- grep("key", names(dat))
    
    # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each species 1 to B) giving probabilities for each possible true group
    group_true_probability <- group.true.probability.key.f(group_probability, group_size_probmass, g)
    
    # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level) and key value (2nd list level), with vectors for each key value giving probabilities of observed groups for each possible true group
    
    group_observed_cprobability_key <- lapply(1:n_observers, function(x)
      apply(keys, 1, function(y)
        group.observed.cprobability.f(y[1:A] , theta_arr[ ,  y[A + 2], x],  B, y[A + 4])
      )
    )
    key_0 <- which(keys$sum == 0)
    
    # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      dat_tmp <- dat[rows_i, ] 
      
      # 'B_states' is the number of possible true groups for the current group size
      B_states <- length(group_true_probability[[paste0(i)]])
      
      # For missing observations (observation sum = 0), set classification probability = 1 for each possible true group so that these records don't influence likelihoods.
      for (j in 1:n_observers) {
        for (k in key_0) {
          group_observed_cprobability_key[[j]][[k]] <- rep(1, B_states)
        }
      }
      
      # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories
      
      probability_mat <- matrix(
        c(
          sapply(1:n_observers, function(x)
            unlist(group_observed_cprobability_key[[x]][unlist(dat_tmp[, key_col[[x]]])])),
          vapply(unlist(dat_tmp$covariate_theta), function(x)
            group_true_probability[[paste0(i)]], numeric(B_states))
        ),
        ncol = n_observers + 1)
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- tibble(index = rep(1:n_i, each = B_states))
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      likelihood[rows_i] <-
        bind_cols(index, data.frame(product = vapply(1:dim(probability_mat)[1], function(x)
          prod(probability_mat[x, ]), numeric(1)))) %>%
        group_by(index) %>%
        summarise(likelihood = sum(product)) %>%
        select(likelihood) %>%
        unlist(.)
    }
  }else{ ##  ----- Compute likelihoods: heterogeneous groups, WITHOUT key table -----

  # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each species 1 to B) giving probabilities for each possible true group
  group_true_probability <- group.true.probability.key.f(group_probability, group_size_probmass, g)
  
  # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      B_states <- length(group_true_probability[[paste(i)]]) # Number of combinations of true species states for the current group size 
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.

      probability_mat <-
        as_tibble(matrix(
          c(
            rep(numeric(1), n_i * B_states * n_observers),
            rep(group_true_probability[[paste(i)]], n_i)
          ), 
          ncol =  n_observers + 1))
      
      # Add condition probabilities of observed groups for each observer
      
      for (obs in 1:n_observers) {
        dat_tmp <- dat[rows_i, ] %>%
          select(all_of((((obs - 1) * A) + 1):(obs * A)))
        record_sum <- rowSums(dat_tmp)
        probability_mat[, obs] <- 
          as.vector(vapply(1:n_i, function(x) 
            group.observed.cprobability.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B_states, record_sum[x]),
            numeric(B_states))) 
      }
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories 
      index <- rep(1:n_i, each = B_states)
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      likelihood[rows_i] <-
        transmute(probability_mat,
                  index = index,
                  b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
        group_by(index)  %>%
        summarise(likelihood = sum(b.product)) %>%
        select(likelihood) %>%
        unlist(.)
    }
  }
## Compute -log(likelihood) as the product of likelihoods for each observation history and the count for each history. Add penalty term(s) accrued for violating model constraints.
  sum(dat[, "count"] * -log(likelihood)) + penalty
}

# Function: {optimize.M.theta.p.f} Optimization of MOM models with misidentification, partial identification, and a covariate predicting classification probabilities (theta) of primary observers only. Supports <=3 true species states B and observation states A. Accepts inputs of initial parameter values 'param', simulated survey data 'dat', and the simulation profile 'sim_profile'. Returns the -log(likelihood) plus any penalty term(s).

optimize.M.theta.p.f <- function(param, dat, sim_profile){
## ----- Import and summarize true and estimated parameters -----
  
  # Summarize group sizes
  g <- unique(dat$group_size) # Observed group sizes
  
  # Extract true simulation parameters from 'sim_profile'
  B <- sim_profile$B; A <- sim_profile$A # Number of true species states and observation states
  n_O_ps <- c(sim_profile$O_p, sim_profile$O_s) # Number of primary/secondary observers
  n_observers <- sum(n_O_ps) # Total observers
  het_true <- sim_profile[grep("mix", names(sim_profile))] # True heterogeneous group parameter(s)
  mx_model <- sim_profile[grep("mx_model", names(sim_profile))] # Heterogeneous group  model
  
  # Summarize values of estimated parameters in 'param'
  
  # Extract value(s) for heterogeneous group parameters to vector 'mix'
  if (any(het_true > 0)) {
    mix <- plogis(param[(length(param) - length(het_true) + 1):length(param)])
  }else{
    mix <- rep(0, length(het_true))
  }
  
  # If secondary observers have certain ID (probability of correct ID = 1), set 'id_certain' = 1
  theta_col <- grep("theta", names(sim_profile))
  id_certain <- (any(sim_profile[theta_col] == 0)) 
  
  # Array 'theta_arr' summarizes classification probabilities (theta) of observers as vectors for efficient computation. Vectors are stored in an array with Dimension 1 (row) = classification probabilities 1 to A for species 1 to B, dimension 2 (column) = unique observation histories or keyed observed groups, and dimension 3 (matrix) = observer. Baseline category is correct classification y = z. Regression coefficients predicting classification probabilities based on group-level covariates are summarized in the 'betas' arrays.
  
  # True species probabilities (psi) for species 1 to (B - 1)
  psi <- plogis(param[(sum(n_parameters[1:3]) - B + 2):sum(n_parameters[1:3])])
  
  n_unique <- dim(dat)[1] # Number of unique observation histories

  if (B == 2) {
    if (A == 2) {
      # True species states (B) = observation states (A) = 2
      
      # Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 (row) = true species states 1 to B, dimension 2 (columns)  = regression coefficients (intercepts in column 1, slopes in column 2). 
      
      betas_arr <- array(param[1:4]
                     , dim = c(2, 2))
      
      # With certain ID by secondary observers, set probability of uncertain ID to 0
      if (id_certain) {
        param_s <- rep(0L, 2)
      }else{
        param_s <- plogis(param[5:6])
      }
      
      theta_arr <- array(0, dim = c(4, n_unique, 4))
      
      theta_arr[c(3, 1), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[1, , drop = FALSE], A))
      theta_arr[c(2, 4), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[2, , drop = FALSE], A))
      
      theta_arr[c(1, 3),  , 2] <- c(1 - param_s[1], param_s[1])
      theta_arr[c(2, 4),  , 2] <- c(param_s[2], 1 - param_s[2])
      
    }else{
      # True species states (B) = 2, observation states (A) = 3
      
      # Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 (row) = observation states 1 to A, dimension (column) 2 = regression coefficients (intercepts in column 1, slopes in column 2), dimension 3 (matrix) = true species state 1 to B.
      
      betas_arr <- array(c(param[1:2], param[5:6], param[3:4], param[7:8]), dim = c(2, 2, 2))
      
      # With certain ID by secondary observers, set probability of uncertain ID to 0
      if (id_certain) {
        param_s <- rep(0L, 4)
      } else{
        param_s <- plogis(param[9:12])
      }
      
      theta_arr <- array(0, dim = c(6, n_unique, 4))
      
      theta_arr[c(3, 5, 1), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[, , 1, drop = FALSE], A))
      theta_arr[c(2, 6, 4), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[, , 2, drop = FALSE], A))
      
      theta_arr[c(1, 3, 5),  , 2] <- c(1 - sum(param_s[1:2]), param_s[1:2])
      theta_arr[c(2, 4, 6),  , 2] <- c(param_s[3], 1 - sum(param_s[3:4]), param_s[4])
      
    }
  } else if (B == 3 & A == 3) {
    # True species states (B) = observation states (A) = 3
    
    # 'Array 'betas_arr' contains beta parameters for multinomial logit regression predicting classification probabilities (theta). Dimension 1 (rows) = observation states 1 to A, dimension (column) 2 = regression coefficients (intercepts in column 1, slopes in column 2), dimension 3 (matrix) = true species state 1 to B.
    
    betas_arr <- array(c(
      param[1:2], param[7:8], param[3:4], param[9:10], param[5:6], param[11:12]), dim = c(2, 2, 3))
    if (id_certain) {
      param_s <- rep(0L, 6)
    } else{
      param_s <- plogis(param[13:18])
    }
    
    theta_arr <- array(0, dim = c(9, n_unique, 4))
    
    theta_arr[c(4, 7, 1), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[, , 1, drop = FALSE], A))
    theta_arr[c(2, 8, 5), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[, , 2, drop = FALSE], A))
    theta_arr[c(3, 6, 9), , 1] <- t(mlogit.regress.predict.f(dat$covariate_theta, betas_arr[, , 3, drop = FALSE], A))
    
    theta_arr[c(1, 4, 7),  , 2] <- c(1 - sum(param_s[1:2]), param_s[1:2])
    theta_arr[c(2, 5, 8),  , 2] <- c(param_s[3], 1 - sum(param_s[3:4]), param_s[4])
    theta_arr[c(3, 6, 9),  , 2] <- c(param_s[5:6], 1 - sum(param_s[5:6]))
  }
  
  # Secondary observers assumed to have identical classification probabilities
  theta_arr[, , 4] <- theta_arr[, , 3]  <- theta_arr[,  , 2]
  
  # Enforce constraint that true species probability parameters sum to 1. The 'penalty' term added to the -log(likelihood) scales with the magnitude of violations of constraints.
  penalty <- 0
  
  if (sum(psi) > 1) {
    penalty <- penalty.f(sum(psi))
    psi <- sum1.f(psi)
  }
  
  # If any groups >1, extract mean group size parameters for species 1 to B to vector 'g_spp'
  g_spp <- rep(1, B) 
  if (any(g > 1)) {
    g_spp <- param[(sum(n_parameters[1:3]) + 1):sum(n_parameters[1:4])]
  }
  
  # True species probabilities (psi) for species 1 to B
  psi <- c(psi, 1 - sum(psi))
  
  # Vector 'lambda' contains lambda parameters of zero-truncated Poisson distribution defining group size by species from mean group size parameters
  lambda <- rep(0, B)
  lambda[which(g_spp > 1)] <- vapply(g_spp[which(g_spp > 1)], function(x)
    optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))
  
## ----- Compute group probabilities and group size probability mass for each species -----
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each species (column) and group size (row) up to maximum observed group size
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = B)
  
  # Compute group probabilities
  if (mx_model == "constant") {
    # Function for "constant" heterogeneous group model and for homogeneous groups
    output <- group.probability.constant.f(psi, mix, g_spp)
  }
  
  # Vector 'group_probability' gives group probabilities (pi) for homogeneous (true states 1 to B) and then heterogeneous groups (if present)
  group_probability <- output[[1]]
  
  # Penalty term added to the -log(likelihood) if heterogeneous group parameters take inadmissible values 
  penalty <- penalty + output[[3]]
  
## ----- Compute -log(likelihood) for data conditional on estimated parameter values -----
  
  # Vector 'likelihood' contains likelihoods for each unique observation history
  likelihood <- numeric(dim(dat)[1])
  
  # 'n_group_size' = count of unique observation histories by group size
  n_group_size <-  dat  %>% count(group_size) 
  
  if (all(het_true == 0)) { ## ----- Compute likelihoods: homogeneous groups -----
    
    # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each species 1 to B) giving probabilities for each possible true group
    
    group_true_probability <-
      lapply(1:dim(group_size_probmass)[1], function(x)
        group_size_probmass[x,] * matrix(group_probability, nrow = 1))
    names(group_true_probability) <- paste(1:max(g))
    
    # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
      
      probability_mat <-
        as_tibble(matrix(
          c(
            rep(numeric(1), n_i * B * n_observers),
            rep(group_true_probability[[i]], n_i)
          ), 
          ncol =  n_observers + 1))
      
      # Add condition probabilities of observed groups for each observer
      
      for (obs in 1:n_observers) {
        dat_tmp <- filter(dat, group_size == i) %>%
          select(all_of((((obs - 1) * A) + 1):(obs * A)))
        record_sum <- rowSums(dat_tmp)
        probability_mat[, obs] <- 
          as.vector(vapply(1:n_i, function(x) 
            group.observed.cprobability.homogeneous.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B, record_sum[x]),
            numeric(B))) 
      }

      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- rep(1:n_i, each = B)
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      
      likelihood[rows_i] <-
        transmute(probability_mat,
                  index = index,
                  b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
        group_by(index)  %>%
        summarise(likelihood = sum(b.product)) %>%
        select(likelihood) %>%
        unlist(.)
    }
  }else{ ## ----- Compute likelihoods: heterogeneous groups -----
    
    # List 'group_true_probability' contains elements named with each observed group size each containing a vector (named with counts of each species 1 to B) giving probabilities for each possible true group
    group_true_probability <- group.true.probability.key.f(group_probability, group_size_probmass, g)
    
    # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      B_states <- length(group_true_probability[[paste(i)]]) # Number of combinations of true groups for current group size
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
      
      probability_mat <-
        as_tibble(matrix(
          c(
            rep(numeric(1), n_i * B_states * n_observers),
            rep(group_true_probability[[paste(i)]], n_i)
          ), 
          ncol =  n_observers + 1))
      
      # Add condition probabilities of observed groups for each observer
      
      for (obs in 1:n_observers) {
        dat_tmp <- filter(dat, group_size == i) %>%
          select(all_of((((obs - 1) * A) + 1):(obs * A)))
        record_sum <- rowSums(dat_tmp)
        probability_mat[, obs] <- 
          as.vector(vapply(1:n_i, function(x) 
            group.observed.cprobability.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B_states, record_sum[x]),
            numeric(B_states))) 
      }

      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- rep(1:n_i, each = B_states)
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      
      likelihood[rows_i] <-
        transmute(probability_mat,
                  index = index,
                  b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
        group_by(index)  %>%
        summarise(likelihood = sum(b.product)) %>%
        select(likelihood) %>%
        unlist(.)
    }
  } 

## Compute -log(likelihood) as the product of likelihoods for each observation history and the count for each history. Add penalty term(s) accrued for violating model constraints.
  sum(dat[, "count"] * -log(likelihood)) + penalty
}

# Function: {optimize.M.psi.f} Optimization of MOM models with misidentification, partial identification, and a covariate predicting true species probabilities (psi). Supports <=3 true species states B and observation states A. Accepts inputs of initial parameter values 'param', simulated survey data 'dat', optional tables of key values for unique observed groups 'keys' and unique covariate values 'keys_psi', and the simulation profile 'sim_profile'. Returns the -log(likelihood) plus any penalty term(s).

optimize.M.psi.f <- function(param, dat, keys, keys_psi, sim_profile){

## ----- Import and summarize true and estimated parameters -----
  
  # Summarize group sizes
  g <- unique(dat$group_size) # Observed group sizes
  
  # Extract true simulation parameters from 'sim_profile'
  B <- sim_profile$B; A <- sim_profile$A # Number of true species states (B) and observation states (A)
  n_O_ps <- c(sim_profile$O_p, sim_profile$O_s) # Number of primary/secondary observers
  n_observers <- sum(n_O_ps) # Total observers
  het_true <- sim_profile[grep("mix", names(sim_profile))] # True heterogeneous group parameter(s)
  mx_model <- sim_profile[grep("mx_model", names(sim_profile))] # Heterogeneous group model
 
  # Summarize values of estimated parameters in 'param'
  
  # Extract value(s) for heterogeneous group parameters to vector 'mix'
  if (any(het_true > 0)) {
    mix <- plogis(param[(length(param) - length(het_true) + 1):length(param)])
  }else{
    mix <- rep(0, length(het_true))
  }
  
  # Extract classification probabilities (theta) to 'theta_tmp' array. Dimension 1 = mis- and partial ID parameters, dimension 2 = true species states 1 to B, dimension 3 = individual observers.
  theta_col <- grep("theta", names(sim_profile))
  
  # If any observers have certain identification (probability of uncertain ID = 0), remove from estimated parameters
  if (any(sim_profile[theta_col] == 0)) {
    theta_col <- theta_col[-which(sim_profile[theta_col] == 0)]
    theta_tmp <- 
      array(c(plogis(param[ (sum(n_parameters[1:2]) + 1):(sum(n_parameters[1:3])) ]), rep(0L, B * (A - 1)) ), c(A - 1, B, n_observers))
    
    # With distinct secondary observers, record classification probabilities of each in separate arrays
  }else if (length(theta_col) == B * (A - 1) * n_observers) {
    theta_tmp <- array(plogis(param[(sum(n_parameters[1:2]) + 1):(sum(n_parameters[1:3]))]), c(A - 1, B, n_observers))
    
    # With identical secondary observers, apply classification probabilities for 1st to all
  }else{
    theta_mat <- matrix(plogis(param[(sum(n_parameters[1:2]) + 1):(sum(n_parameters[1:3]))]),
                       ncol = B * (A - 1), byrow = TRUE)
    theta_tmp <- array(0, c(A - 1, B, n_observers))
    for (i in 1:n_observers) {
      if (dim(theta_mat)[1] >= i) {
        theta_tmp[, , i] <- theta_mat[i, ]
      }else{
        theta_tmp[, , i] <- theta_tmp[, , (i - 1)]
      }
    }
  }
  
  # Extract true species probability regression coefficients to matrix 'psi_betas_mat'
  psi_betas_mat <- matrix(param[1:(sum(n_parameters[1:2]))], 
                      nrow = (B - 1), byrow = TRUE)
  
  # Enforce constraints for probabilities that must sum to 1. The 'penalty' term added to the -log(likelihood) increases with the magnitude of violations of constraints. 
  penalty <- 0
  
  if (any(colSums(theta_tmp) > 1)) {
    sums <- colSums(theta_tmp)
    penalty <- penalty + penalty.f(sums[sums > 1])
    col <- which(colSums(theta_tmp) > 1, arr.ind = TRUE)
    for (i in 1:dim(col)[1]) {
      theta_tmp[, col[i, 1], col[i, 2]] <- sum1.f(theta_tmp[, col[i, 1], col[i, 2]])
    }
  }
  
  # With group sizes >=1, extract group size parameters for species 1 to B to vector 'g_spp'
  g_spp <- rep(1, B) 
  if (any(g > 1)) {
    g_spp <- param[(sum(n_parameters[1:3]) + 1):sum(n_parameters[1:4])]
  }
  
  # Matrix 'theta_mat' summarizes classification probabilities (theta) of each observer as a vector for efficient computation
  # row = observers, column = classification probabilities 1 to A in sequence for spp 1 to B

  n_unique <- dim(dat)[1]
  theta_mat <- t(vapply(1:n_observers, theta4.f, numeric(B * A), theta_tmp))
  
  # Vector 'lambda' contains lambda parameters (defined by mean group size) of zero-truncated Poisson distribution for species 1 to B
  lambda <- rep(0, B)
  lambda[which(g_spp > 1)] <- vapply(g_spp[which(g_spp > 1)], function(x)
    optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))
  
## ----- Compute group probabilities and group size probability mass for each species -----
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each species (column) and group size (row) up to maximum observed group size
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = B)
  
  # Compute group probabilities
  
  if (!is.null(keys) & (any(het_true > 0))) {
    # If key table for observed groups exists and heterogeneous groups present, compute group probabilities for keyed table
    group_probability_key <- mlogit.regress.predict.f(keys_psi$covariate_psi, psi_betas_mat, B)
    
    # Compute group probability for the specified heterogeneous group model
    if (mx_model == "constant") {
      # Function for "constant" heterogeneous group model
      output <- group.probability.psi.constant.f(group_probability_key, mix, g_spp)
    } else if (mx_model == "encounter") {
      # Function for "encounter" heterogeneous group model
      output <- group.probability.psi.encounter.f(group_probability_key, mix, g_spp)
    }

    group_probability_key <- output[[1]]
  }else{
    # Without key table, compute matrix of group probabilities for unique observation histories
    group_probability <- mlogit.regress.predict.f(dat$covariate_psi, psi_betas_mat, B)
    
    # Compute group probability for the specified heterogeneous group model
    if (mx_model == "constant") {
      # Function for "constant" heterogeneous group model
      output <- group.probability.psi.constant.f(group_probability, mix, g_spp)
    } else if (mx_model == "encounter") {
      # Function for encounter" heterogeneous group model
      output <- group.probability.psi.encounter.f(group_probability, mix, g_spp)
    }
    
    # Matrix 'group_probability' gives group probabilities (pi) for homogeneous (true states 1 to B) and then heterogeneous groups (if present) for each unique observation history (rows)
    group_probability <- output[[1]]
  }
  
  if (is.null(keys)) {
    ## If key table of unique observed groups is NOT present, compute likelihoods from individual unique observation histories
    
    # Array 'theta_arr' summarizes classification probabilities (theta)  as an array for computing probabilities for group size = 1.
    # dimension 1 (row) = true spp state 1 to B, dimension 2 (column) = classification probabilities for observation states 1 to A, dimension 3 (matrix) = observer. 
    # Secondary observers assumed to have identical classification probabilities
    theta_arr <- array(
      apply(theta_mat, 1, function(x) matrix(x, nrow = n_unique, ncol = B * A))
      , dim = c(B * A, n_unique, n_observers))
  }
  
  # Penalty term added to the -log(likelihood) if heterogeneous group parameters take inadmissible values 
  penalty <- penalty + output[[3]]
  
## ----- Compute -log(likelihood) for data conditional on estimated parameter values -----
  
  # Vector 'likelihood' contains likelihoods for each unique observation history
  likelihood <- numeric(dim(dat)[1])
  # 'n_group_size' = count of unique observation histories by group size
  n_group_size <-  dat  %>% count(group_size)
  
  if (all(het_true == 0)) { ## ----- Compute likelihoods: homogeneous groups -----
    
    # Matrix 'group_true_probability_mat' contains probabilities for possible true groups, with rows for each unique observation history and columns for each true species state 1 to B.
    group_true_probability_mat <- group_size_probmass[dat$group_size, ] * group_probability
    
    # If key table of unique observed groups is present, compute likelihoods from keyed table of probabilities
    if (!is.null(keys)) {
      key_col <- grep("key", names(dat))
      key_psi_col <- grep("psi_key", names(dat))
      key_col <- setdiff(key_col, key_psi_col)
      
      # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level). Matrices for each observer contain probabilities of the keyed observed group (columns) conditional on each possible true group (rows).
      
      group_observed_cprobability_key <- lapply(1:n_observers, function(x)
        apply(keys, 1, function(y)
          group.observed.cprobability.homogeneous.f(y[1:A] , theta_mat[x, ],  B, y[A + 3])
        )
      )

      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- tibble(index = unlist(lapply(1:dim(n_group_size)[1], function(x)
        rep(1:n_group_size[[x, 2]], each = B) + sum(n_group_size[n_group_size[[1, 1]]:n_group_size[[x, 1]] , 2]) - n_group_size[[x, 2]] ))
      )
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
      
      probability_mat <- matrix(
        c(
          sapply(1:n_observers, function(x)
            unlist(group_observed_cprobability_key[[x]][, unlist(dat[, key_col[[x]]])])),
          t(group_true_probability_mat)
        ),
        ncol = n_observers + 1)
      
      # Add likelihoods for group sizes >1 to 'likelihood'. Likelihoods for each unique observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups. 
      
      likelihood <-
        bind_cols(index, data.frame(product = vapply(1:dim(probability_mat)[1], function(x)
          prod(probability_mat[x, ]), numeric(1)))) %>%
        group_by(index) %>%
        summarise(likelihood = sum(product)) %>%
        select(likelihood) %>%
        unlist(.)
      
    }else{
      ## If key table of unique observed groups is NOT present, compute likelihoods from individual unique observation histories
      
      # Loop for observed group sizes
      for (i in g) {
        # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
        n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
        rows_i <- which(dat$group_size == i) 
        
        # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.

        probability_mat <-
          as_tibble(matrix(
            c(
              rep(numeric(1), n_i * B * n_observers),
              t(group_true_probability_mat[rows_i, ])
            )
            , ncol =  n_observers + 1))
        
        # Add conditional probabilities of observed groups for each observer
        
        for (obs in 1:n_observers) {
          dat_tmp <- filter(dat, group_size == i) %>%
            select(all_of((((obs - 1) * A) + 1):(obs * A)))
          record_sum <- rowSums(dat_tmp)
          probability_mat[, obs] <-
            as.vector(vapply(1:n_i, function(x)
              group.observed.cprobability.homogeneous.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B, record_sum[x]),
              numeric(B)))
        }
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
        index <- rep(1:n_i, each = B)
        
        # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
        
        likelihood[rows_i] <-
          transmute(probability_mat,
                    index = index,
                    b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
          group_by(index)  %>%
          summarise(likelihood = sum(b.product)) %>%
          select(likelihood) %>%
          unlist(.)
      }
    } # End of loop for homogeneous groups
  }else{ ## ----- Compute likelihoods: heterogeneous groups -----
    
    # WITH key table of unique observed groups, compute likelihoods from keyed table of probabilities
      if (!is.null(keys)) {
      key_col <- grep("key", names(dat))
      key_psi_col <- grep("psi_key", names(dat))
      key_col <- setdiff(key_col, key_psi_col)
      
      # List 'group_observed_cprobability_key' contains keyed probabilities for each observer (1st list level) and key value (2nd list level), with vectors for each key value giving probabilities of observed groups for each possible true group
      
      group_observed_cprobability_key <- lapply(1:n_observers, function(x)
        apply(keys, 1, function(y)
          group.observed.cprobability.f(y[1:A] , theta_mat[x, ],  B, y[A + 3])
        )
      )

      # List 'group_true_probability_key' contains probabilities of possible true groups, with first list level corresponding the key value of the unique covariate from 'keys_psi', the second list level corresponding to group size containing a vector with probabilities for possible true groups
      group_true_probability_key <- apply(group_probability_key, 1, group.true.probability.key.f, group_size_probmass, g)
      
      # Loop for observed group sizes
      for (i in g) {
        # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
        n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
        rows_i <- which(dat$group_size == i) 
        dat_tmp <- filter(dat, group_size == i) 
        
        # 'B_states' is the number of possible true groups for the current group size
        B_states <- length(group_true_probability_key[[1]][[paste0(i)]])
        
        # For missing observations (observation key # = 1), set classification probability = 1 for each possible true group so that these records don't influence likelihoods
        for (j in 1:n_observers) {
          group_observed_cprobability_key[[j]][[1]] <- rep(1, B_states)
        }
        
        # Matrix 'probability_mat' contains probabilities of observations by each observer (cols = 1 to n_observers) and probabilities of true groups (last col). Rows correspond to each possible true group for unique observation histories.
        
        probability_mat <- matrix(
          c(
            sapply(1:n_observers, function(x)
              unlist(group_observed_cprobability_key[[x]][unlist(dat_tmp[, key_col[[x]]])])),
            vapply(unlist(dat_tmp[, key_psi_col]), function(x)
              group_true_probability_key[[x]][[paste0(i)]], numeric(B_states)
            )
          ),
          ncol = n_observers + 1)
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
        index <- tibble(index = rep(1:n_i, each = B_states))
        
        # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
        
        likelihood[rows_i] <-
          bind_cols(index, data.frame(product = vapply(1:dim(probability_mat)[1], function(x)
            prod(probability_mat[x, ]), numeric(1)))) %>%
          group_by(index) %>%
          summarise(likelihood = sum(product)) %>%
          select(likelihood) %>%
          unlist(.)
      }
    }else{
      # WITHOUT key table of unique observed groups, compute likelihoods from individual unique observation histories
      
      ## Loop for observed group sizes
      for (i in g) {
        # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
        n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
        rows_i <- which(dat$group_size == i) 
        
        # Matrix 'group_true_probability' contains probabilities of possible true groups, with columns for each unique observation history of the current group size and rows giving probabilities for each possible true group.
        group_true_probability <-
          apply(group_probability[rows_i, , drop = FALSE], 1, group.true.probability.f, group_size_probmass, i)
        
        # 'B_states' is the number of possible true groups for the current group size
        B_states <- dim(group_true_probability)[1]
        
        # Matrix 'probability_mat' contains probabilities of observations by each observer (cols = 1 to n_observers) and probabilities of true groups (last col). Rows correspond to each possible true group for unique observation histories.
        
        probability_mat <-
          as_tibble(matrix(
            c(
              rep(numeric(1), n_i * B_states * n_observers),
              group_true_probability
            ),
            ncol =  n_observers + 1))
        
        # Add conditional probabilities of observed groups for each observer
        
        for (obs in 1:n_observers) {
          dat_tmp <- filter(dat, group_size == i) %>%
            select(all_of((((obs - 1) * A) + 1):(obs * A)))
          record_sum <- rowSums(dat_tmp)
          probability_mat[, obs] <-
            as.vector(vapply(1:n_i, function(x)
              group.observed.cprobability.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B_states, record_sum[x]),
              numeric(B_states)))
        }
        
        # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
        index <- rep(1:n_i, each = B_states)
        
        # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
        
        likelihood[rows_i] <-
          transmute(probability_mat,
                    index = index,
                    b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
          group_by(index)  %>%
          summarise(likelihood = sum(b.product)) %>%
          select(likelihood) %>%
          unlist(.)
      }
    }
  }
## Compute -log(likelihood) as the product of likelihoods for each observation history and the count for each history. Add penalty term(s) accrued for violating model constraints.
  sum(dat[, "count"] * -log(likelihood)) + penalty
}

# Function: {optimize.M.theta.psi.f} Optimization of MOM models with misidentification, partial identification, and independent covariates predicting classification probabilities (theta) and true species probabilities (psi). Supports 2 true species states B and observation states A. Accepts inputs of initial parameter values 'param', simulated survey data 'dat', and the simulation profile 'sim_profile'. Returns the -log(likelihood) plus any penalty term(s).

optimize.M.theta.psi.f <- function(param, dat, sim_profile){
  
## ----- Import and summarize true and estimated parameters -----
  
  # Summarize group sizes
  g <- unique(dat$group_size) # Group sizes of observations
  
  # Extract true simulation parameters from 'sim_profile'
  B <- sim_profile$B; A <- sim_profile$A # Number of true spp states and observation states
  n_O_ps <- c(sim_profile$O_p, sim_profile$O_s) # Number of primary/secondary observers
  n_observers <- sum(n_O_ps) # Total observers
  het_true <- sim_profile[grep("mix", names(sim_profile))] # True heterogeneous group parameter(s)
  mx_model <- sim_profile[grep("mx_model", names(sim_profile))] # Heterogeneous group  model
  
  # Summarize values of estimated parameters in 'param'
  
  # Extract value(s) for heterogeneous group parameters to vector 'mix'
  if (any(het_true > 0)) {
    mix <- plogis(param[(length(param) - length(het_true) + 1):length(param)])
  }else{
    mix <- rep(0, length(het_true))
  }

  # Extract multinomial regression coefficients for true species probabilities to matrix 'psi.betas,mat' (intercepts in column 1 and slopes in column 2) and true species state 1 to B (rows).
  # Extract multinomial regression coefficients for classification probabilities to matrices 'theta_b0' and 'theta_b1'  (b0 = intercept parameters and b1 = slope parameters). Columns give coefficients for each classification probability, and observers are on each row. 

  n_psi_param <- length(grep("psi", names(sim_profile)))
  psi_betas_mat <- matrix(param[1:n_psi_param], nrow = (B - 1), byrow = TRUE)
  theta_b0 <- length(grep("b0_theta", names(sim_profile)))
  theta_b1 <- length(grep("b1_theta", names(sim_profile)))
  theta_b1 <- matrix(param[(n_psi_param + theta_b0 + 1):(sum(n_psi_param, theta_b0, theta_b1))]
                   , nrow = 2, ncol = 2, byrow = TRUE)
  theta_b0 <- matrix(param[(n_psi_param + 1):(n_psi_param + theta_b0)]
                   , nrow = 2, ncol = 2, byrow = TRUE)

  # The 'penalty' term is added to the -log(likelihood) and increases with the magnitude of violations of constraints
  penalty <- 0
  
  # If any groups >1, extract mean group size parameters for species 1 to B to vector 'g_spp'
  if (any(g > 1)) {
    g_spp <- param[(sum(n_parameters[1:3]) + 1):sum(n_parameters[1:4])]
  }else {
    g_spp <- rep(1, B)
  }
  
  # Array 'theta_arr' summarizes classification probabilities (theta) of observers as vectors for efficient computation. Vectors are stored in an array with Dimension 1 (row) = classification probabilities 1 to A for species 1 to B, dimension 2 (column) = unique observation histories, and dimension 3 (matrix) = observer. Multinomial logit link functions enforce that each column (dimension 2) of each array sums to 1. Baseline category is correct classification y = z. Regression coefficients predicting classification probabilities based on group-level covariates are summarized in the 'theta_betas_arr' array.
  
  n_unique <- dim(dat)[1] # Number of unique observation histories

  # Code for B = 2 with A = 2 or A = 3 only
  if (B == 2) {
    if (A == 2) {
      # True species states (B) = observation states (A) = 2
      
      # Array 'theta_betas_arr' summarizes regression coefficients predicting classification probabilities. Dimension 1 (row) = true species states 1 to B, dimension 2 (column) = intercept (column 1) and slope (column 2) coefficients, dimension 3 (matrix) = observer.
      theta_arr <- array(0, dim = c(4, n_unique, 4))
      theta_betas_arr <- array(c(theta_b0[1, ], theta_b1[1, ], theta_b0[2, ], theta_b1[2, ])
                     , dim = c(2, 2, 2))    
      
      theta_arr <- array(0, dim = c(4, n_unique, 4))
      for (b in 1:B) {
        theta_arr[c(3, 1), , b] <- t(mlogit.regress.predict.f(dat$covariate_theta, theta_betas_arr[1, , b, drop = FALSE], A))
        theta_arr[c(2, 4), , b] <- t(mlogit.regress.predict.f(dat$covariate_theta, theta_betas_arr[2, , b, drop = FALSE], A))
      }
    }else{
      # True species states (B) = 2, observation states (A) = 3 
      
      # Array 'theta_betas_arr' summarizes regression coefficients predicting classification probabilities. Dimension 1 (row) = observation states 1 to A, dim 2 (col) = regression intercepts (col 1) and slope coefficients (col 2), dim (matrix) 3 = observer (primary = 1/2, secondary = 3/4) and true species state (state 1 = 1/3, state 2 = 2/4).
      theta_arr <- array(0, dim = c(6, n_unique, 4))
      theta_betas_arr <- array(c(theta_b0[1, ], theta_b1[1, ], theta_b0[2, ], theta_b1[2, ], 
                       theta_b0[3, ], theta_b1[3, ], theta_b0[4, ], theta_b1[4, ])
                     , dim = c(2, 2, 4))
      
      for (b in 1:B) {
        theta_arr[c(3, 5, 1), , b] <- t(mlogit.regress.predict.f(dat$covariate_theta, theta_betas_arr[, , b + (b - 1), drop = FALSE], A))
        theta_arr[c(2, 6, 4), , b] <- t(mlogit.regress.predict.f(dat$covariate_theta, theta_betas_arr[, , b * 2, drop = FALSE], A))
      }
    }
  }
  
  # Secondary observers assumed to have identical classification probabilities
  theta_arr[, , 4] <-   theta_arr[, , 3] <- theta_arr[, , 2]
  
  # Vector 'lambda' contains lambda parameters (defined by mean group size) of zero-truncated Poisson distribution for species 1 to B
  lambda <- rep(0, B)
  lambda[which(g_spp > 1)] <- vapply(g_spp[which(g_spp > 1)], function(x)
    optimize(ztpois.f, interval = c(0.000001, 100), size = x)$minimum, numeric(1))
  
## ----- Compute group probabilities and group size probability mass for each species -----
  
  # Matrix 'group_size_probmass' contains probability mass (mu) for each species (column) and group size (row) up to maximum observed group size
  group_size_probmass <-
    matrix(vapply(lambda, function(x)
      ztpois.probmass.f(x, g), numeric(max(g))),
      ncol = B)
  
  # Compute matrix of group probabilities for unique observation histories
  group_probability <- mlogit.regress.predict.f(dat$covariate_psi, psi_betas_mat, B)
  
  # Compute group probability for the specified heterogeneous group model
  if (mx_model == "constant") {
    # Function for "constant" heterogeneous group model and for homogeneous groups
    output <- group.probability.psi.constant.f(group_probability, mix, g_spp)
  } else if (mx_model == "encounter") {
    # Function for proportional encounter mixing model
    output <- group.probability.psi.encounter.f(group_probability, mix, g_spp)
  }
  
  # Matrix 'group_probability' gives group probabilities (pi) for homogeneous (true states 1 to B) and then heterogeneous groups (if present) for each unique observation history (rows)
  group_probability <- output[[1]]
  
  # Penalty term added to the -log(likelihood) if heterogeneous group parameters take inadmissible values 
  penalty <- penalty + output[[3]]
  
## ----- Compute -log(likelihood) for data conditional on estimated parameter values -----
  
  # Vector 'likelihood' contains likelihoods for each unique observation history
  likelihood <- numeric(dim(dat)[1])
  
  # 'n_group_size' = count of unique observation histories by group size
  n_group_size <-  dat  %>% count(group_size)

  if (het_true == 0) { ## ----- Compute likelihoods: homogeneous groups -----
    
    # Matrix 'group_true_probability_mat' contains probabilities for possible true groups, with rows for each unique observation history and columns for each true species state 1 to B.
    group_true_probability_mat <- group_size_probmass[dat$group_size, ] * group_probability
    
    # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      
      # Columns in matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups. Rows correspond to possible true groups for unique observation histories.
      
      probability_mat <-
        as_tibble(matrix(
          c(
            rep(numeric(1), n_i * B * n_observers),
            t(group_true_probability_mat[rows_i, ])
          )
          , ncol =  n_observers + 1,
          dimnames = list(c(), c(paste0("obs", 1:n_observers), "group"))))
      
      # Add conditional probabilities of observed groups for each observer
      
      for (obs in 1:n_observers) {
        dat_tmp <- filter(dat, group_size == i) %>%
          select(all_of((((obs - 1) * A) + 1):(obs * A)))
        record_sum <- rowSums(dat_tmp)
        probability_mat[, obs] <- 
          as.vector(vapply(1:n_i, function(x) 
            group.observed.cprobability.homogeneous.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B, record_sum[x]),
            numeric(B))) 
      }
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- rep(1:n_i, each = B)
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      
      likelihood[rows_i] <-
        transmute(probability_mat,
                  index = index,
                  b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
        group_by(index)  %>%
        summarise(likelihood = sum(b.product)) %>%
        select(likelihood) %>%
        unlist(.)
    } # End homogeneous groups loop
  }else{ ## ----- Compute likelihoods: heterogeneous groups -----
    
    # Loop for observed group sizes
    for (i in g) {
      # For current group size, 'n_i' and 'rows_i' give the count of unique observation histories and a vector of row numbers in 'dat'
      n_i <- n_group_size[[pmatch(i, n_group_size$group_size), 2]] 
      rows_i <- which(dat$group_size == i) 
      
      # Matrix 'group_true_probability_mat' contains probabilities for possible true groups, with rows for each unique observation history columns for each true species state 1 to B.
      group_true_probability_mat <-
        apply(group_probability[rows_i, , drop = FALSE], 1, group.true.probability.f, group_size_probmass, i)
      
      B_states <- dim(group_true_probability_mat)[1] # Number of combinations of true states for the current group size
      
      # Columns in the matrix 'probability_mat' contain probabilities of observed groups (conditional on possible true groups) for each observer and probabilities of true groups (last col). Rows correspond to possible true groups for unique observation histories.
      
      probability_mat <-
        as_tibble(matrix(c(
          rep(numeric(1), n_i * B_states * n_observers),
          group_true_probability_mat
        ),
        ncol =  n_observers + 1))
      
      # Add conditional probabilities of observed groups for each observer
      
      for (obs in 1:n_observers) {
        dat_tmp <- filter(dat, group_size == i) %>%
          select(all_of((((obs - 1) * A) + 1):(obs * A)))
        record_sum <- rowSums(dat_tmp)
        probability_mat[, obs] <- 
          as.vector(vapply(1:n_i, function(x) 
            group.observed.cprobability.f(dat_tmp[x, ], theta_arr[, rows_i[x], obs], B_states, record_sum[x]),
            numeric(B_states))) 
      }
      
      # Positive integers in 'index' associate rows in 'probability_mat' with sequentially numbered unique observation histories
      index <- rep(1:n_i, each = B_states)
      
      # Add likelihoods for current group size to 'likelihood'. Likelihoods for each observation history are calculated as the product of each row in 'probability_mat' summed across possible true groups.
      
      likelihood[rows_i] <-
        transmute(probability_mat,
                  index = index,
                  b.product = apply(probability_mat, 1, function(x) prod(x))) %>%
        group_by(index)  %>%
        summarise(likelihood = sum(b.product)) %>%
        select(likelihood) %>%
        unlist(.)
    }
  }
  
## Compute -log(likelihood) as the product of likelihoods for each observation history and the count for each history. Add penalty term(s) accrued for violating model constraints.
  sum(dat[, "count"] * -log(likelihood)) + penalty
}
