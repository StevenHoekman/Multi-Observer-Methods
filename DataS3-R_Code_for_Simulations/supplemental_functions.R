# Supplemental_functions.R
# Supplemental functions for analyses and summaries, version 1.2.3
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# R computer code with supplemental functions for simulation analyses for estimating uncertain identification using multi-observer methods. These functions provide supplemental services such as: drawing random samples from probability distributions; computing probabilities and other values used in likelihood computations; and generating, formatting, summarizing, and error-checking simulation data and statistical output. Functions are grouped according by purpose. Comments with each function describe its purpose, inputs and outputs, and functioning of the code. Code developed and tested in R version 4.1.

# This code should be executed prior to conducting simulation analyses in 'simulations_R_code.R'

###############################################################################
#             Utility functions
###############################################################################

# Function: {multinomial.f} Multinomial function. Accepts vector of multinomial probabilities (summing to <1), returns vector of multinomial beta parameters.

multinomial.f <- function(p) {
  vapply(p, function(x) log(x / (1 - sum(p))), numeric(1))}

# Function: {multinomial.inv.f} Inverse multinomial function. Accepts matrix of sets of multinomial beta parameters (rows), returns matrix of multinomial probabilities.

multinomial.inv.f <- function(betas) {
  denom <- sapply(betas, function(x) exp(x)) 
  denom <- 1 + sum(denom)
  est <- sapply(betas, function(x) exp(x) / denom)
  est <- c(est, 1 - sum(est))
}

###############################################################################
#             Functions for probability distributions
###############################################################################

# Function: {ztpois.f} Solves for 'lambda' parameter of zero-truncated Poisson distribution. Accepts inputs of mean group size 'size' and 'lambda' parameter of a zero-truncated Poisson distribution, returns the absolute difference between mean group size for these. Solving for difference of 0 gives 'lambda' for a given 'size'.

ztpois.f <- function(lambda, size){
  abs(size - lambda / (1 - exp(-lambda)))
}

# Function: {ztpois.probmass.f} Probability mass for zero-truncated Poisson distribution. Accepts inputs of 'lambda' parameter of a zero-truncated Poisson distribution and vector 'size' of observed group sizes, returns vector with probability mass for group sizes 1 to maximum 'size'. For 'lambda' = 0, returns probability = 1 for group size = 1.  

ztpois.probmass.f <- function(lambda, size) {
  if (lambda == 0) {
    return(as.numeric(size == 1))
  }else{
    lambda[1] ^ (1:max(size)) / ((exp(lambda[1]) - 1) * factorial(1:max(size)))
  }
}

# Function: {logistic.dist.f} Distribution for logistic model. Accepts input of vector of beta parameters (intercept and slope) for the logistic model, returns vector with estimated mean and standard deviation of the distribution.

logistic.dist.f <- function(beta) {
  dist <- plogis(beta[1] + beta[2] * rnorm(10^6))
  c(mean(dist), sd(dist))
}

# Function: {mlogit.dist.f} Distribution for multinomial logit model. Accepts input of matrix or vector beta parameters (intercept and slope) for separate multinomial logit models, returns matrix with estimated mean and standard deviation of each distribution.

mlogit.dist.f <- function(beta) {
  # Convert vectors containing a set of regression coefficients to matrix with 1 row
  if (is.vector(beta)) {
    beta <- matrix(beta, ncol = 2, byrow = TRUE)
  }
  r <- rnorm(10^7.3)
  
  # Linear predictors for each logit regression
  predict <-
    t(vapply(1:dim(beta)[1], function(x)
      exp(beta[x, 1] + beta[x, 2] * r), numeric(length(r)) ))
  
  # Compute distribution of predicted values by category
  distribution <- 
    predict %r/% (1 + fsum(predict)  ) %>%
    {
      t(rbind(., 1 - colsums(.) ))  
    } %>%
    {
      cbind(colmeans(.), fsd(.))
    }
  distribution
}

# Function: {mlogit.regress.predict.f} Predicts probabilities for multinomial logistic regression. Accepts input of vector of covariate values 'cov_mat' and a matrix of beta parameters 'beta_mat' (columns for intercept and slope, rows for each logit regression equation). Returns a matrix with category-specific (columns 1 to 'n') predicted probabilities for each covariate value (rows).

mlogit.regress.predict.f <- function(cov_mat, beta_mat) {
  
  # tidy up the inputs to matrices with 2 covariates, 3 beta parameters
  
  if (!is.matrix(cov_mat))
    cov_mat <- qM(cov_mat)
  if (is.vector(beta_mat))
    beta_mat <- matrix(beta_mat, byrow = TRUE, ncol = 2)
  else if (is.array(beta_mat))
    beta_mat <- qM(beta_mat)
  
  if (dim(cov_mat)[2] == 1)
    cov_mat <- cbind(cov_mat, 0)
  if (dim(beta_mat)[2] == 2)
    beta_mat <- cbind(beta_mat, 0)
  
  n_cat <- dim(beta_mat)[1] + 1
  
  # Compute distribution of predicted values by category
  
  distribution <-
    t(vapply(1:dim(beta_mat)[1], function(b)
      exp(beta_mat[b, 1] + beta_mat[b, 2] * cov_mat[, 1] + 
            beta_mat[b, 3] * cov_mat[, 2]), 
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
#             Functions for likelihood computations
###############################################################################

# Several functions in this section rely on an external 'likelihood.equations' list objects containing pre-calculated values for equations used to compute likelihoods. See DataS2 and DataS3 for additional details on these lists.

# Function: {group.observed.cprobability.f} Computes conditional (on the true group) probabilities for group records with heterogeneous groups. Accepts inputs of vectors with a group record 'd', classification probabilities 'p', the number of possible combinations (with repetition) of true groups 'states', and group size 's'. Returns vector with conditional probabilities for each possible true group. Relies on external 'likelihood.equations' list objects.

group.observed.cprobability.f <- function(d, p, states, s){  
  # Return probability = 1 for missing data
  if (s < 1) return(rep(1, states))
  
  # For each combination of true states (possible true groups) and each permutation of possible observation states, compute probabilities by 1) raising classification probabilities to the power of the number of classifications in corresponding observation states and 2) then multiplying by the multinomial coefficient for each permutation.
  # For observed group 'd', matrix 't1' contains an 'index' value corresponding to each possible true group, a multinomial coefficient term 'coefficient', and additional columns containing vectors representing permutations of classifications giving rise to each observed group.

  # Compute overall probabilities for possible true groups by summing across permutations of observation states in the observed group with identical index values
  
  t1 <- likelihood.equations[[paste0("observed.", paste0(d, collapse = "."))]]
  
  t2 <-
    dapply(t1[, -1], MARGIN = 1, function(z)
      prod(p ^ z[-1] , z[1])) %>%
    fsum(., t1[, 1]) 
}

# Function: {group.observed.cprobability.homogeneous.f} Computes conditional (on true group) probabilities for group records with homogeneous groups. Accepts inputs of a group record 'd', a vector of classification probabilities 'p', and the number of possible combinations (with repetition) of true groups 'states', and group size 's'. Returns vector with conditional probabilities for each possible true group. Relies on external 'likelihood.equations' list objects.

group.observed.cprobability.homogeneous.f <- function(d, p, states, s){
  # Return probability = 1 for missing data
  if (s < 1) return(rep(1, states))
  
  # For each of 2 possible true groups and each permutation of possible observation states, compute probabilities by 1) raising classification probabilities to the power of the number of classifications in corresponding observation states and 2) then multiplying by the multinomial coefficient for each permutation.
  # For observed group 'd', matrix 't1' contains an 'index' value corresponding to each possible true group, a multinomial coefficient term 'coefficient', and additional columns containing vectors representing permutations of classifications giving rise to each observed group.
  
  t1 <- likelihood.equations[[paste0("observed.", paste0(d, collapse = "."))]]
  
  # Row numbers in 'homogeneous' exclude heterogeneous groups 
  
  homogeneous <- 
    which(likelihood.equations[[paste0("true.permutations.count.g.", s)]] == 1) %>%
    vapply(., function(x) which(t1[, 1] == x), numeric(1)) 
  
  t2 <-
    dapply(t1[homogeneous, -1], MARGIN = 1, function(z)
      prod(p ^ z[-1] , z[1])) %>%
    fsum(., t1[homogeneous, 1])
}

# Function: {multinom.likelihood.f} Computes likelihood for multinomial data with no covariates. Used for the multinomial models  estimating observed species proportions assuming no misidentification. Accepts input of vector with counts for observation classes (1 to A), returns 'multinom' model object.  

multinom.likelihood.f <- function(dat) {
  dat <- matrix(dat, nrow = 1)
  multinom(dat ~ 1, Hess = TRUE, trace = FALSE)
}

# Function: {probability.key.f} Retrieves probabilities for heterogeneous observed groups from a keyed table. Accepts inputs of 'd' (group size and key value for an observed group), an index 'o' to the observer (1 to the total # of observers), the 'n.group' table, and a keyed table for unique observed groups 'prob.key'. Returns a vector (length of possible combinations of true groups) with corresponding probabilities of the observed group (or for missing data, probabilities = 1).

probability.key.f <- function(d, o, n.group, prob.key) {
  if (d[2] != 1) {
    return(prob.key[[o]][d[2]])
  }else{
    rep(1, n.group$B_states[which(n.group$group_size == d[1])])
  }
}

# Function: {probability.key.homogeneous.f} Retrieves probabilities for homogeneous observed groups from a keyed table. Accepts inputs of 'd' (group size and key value for an observed group), an index 'o' to the observer (1 to the total # of observers), and a keyed table for unique observed groups 'prob.key'.  Returns a vector (length of 2 possible true groups) with corresponding probabilities of the observed group (or for missing data, probabilities = 1).

probability.key.homogeneous.f <- function(d, o, prob.key) {
  if (d[2] != 0) {
    return(prob.key[[o]][d[2]])
  }else{
    rep(1, B)
  }
}

# Function: {sum1.f} Enforces constraint that probabilities must sum  to 1. Accepts vector of probabilities and returns vector with probabilities normalized to sum to slightly <1.
sum1.f <- function(x) {(x / sum(x)) * 0.999}

# Function: {penalty.f} Penalty function to -log(likelihood) when sum of probabilities >1. Accepts vector of probabilities and returns vector with penalty term that scales with the sum. 
penalty.f <- function(x) {2000 * (sum(x - 1) + 0.0001) ^ 2}

# Function: {group.true.probability.key.f}: With heterogeneous groups (composed of species pairs) for 2 or 3 species, computes probabilities for all possible true groups. Accepts inputs of group probabilities 'group_probability', probability mass of group sizes (1 to maximum observed group size) for each species 'size_probability', and observed group sizes 'size'. Returns list with named elements for each group size composed of a vector with probabilities for each possible true group. For true species states = 3, relies on external 'likelihood.equations' list to provide possible true groups.

group.true.probability.key.f <- function(group_probability, size_probability, size){

  # True species states B is 'BB' to avoid conflict with collapse::B 
  BB <- dim(size_probability)[2]
  
  if (BB == 2) {
    # True species states B = 2
    
    # Compute probabilities for homogeneous groups at the start/end of each vector, then add probabilities for heterogeneous groups with decreasing/increasing numbers of species 1/species 2
    # Name vector elements with count of species 1/species 2
    
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
    
  }else if (BB == 3) {
    # True species states B = 3
    
    # Add probability = 1 for group size = max(group size) + 1
    size_probability <- 
      rbind(size_probability, 1) %>%
      setRownames(., NULL)
    
    # Add extra vector elements with group probability = 0 for heterogeneous groups with 3 species (assumed to not occur)
    group_probability <- c(group_probability, 0, 0)
    
    # List 'counts' contains a matrix for each group size with the count of each species in possible true groups (from 'likelihood.equations') on each row
    counts <-
      lapply(size, function(x)
        likelihood.equations[[paste0("true.combinations.g.", x)]]) %>%
      lapply(., function(x)
        vapply(1:B, function(y)
          apply(x, 1, function(z)
            sum(z == y)), integer(dim(x)[1])))
    
    # Counts of each species in each true group
    labels.spp <- NULL
    for (i in 1:length(size)) {
      labels.spp[[i]] <-
        apply(counts[[i]], 1, function(x) paste0(x, collapse = ""))
    }
    
    # Each vector in 'index.group' corresponds to matrices in 'counts', with vector elements containing the position of the group probability in 'group_probability' for the true group on each row. Heterogeneous groups with 3 species receive index values corresponding to group probability = 0. 
    index.group <- 
      lapply(counts, function(x)
        apply(x, 1, function(y)
          sum((y > 0) * 1:B) + sum(y > 0) - 1))
    
    # Index locations in each matrix with value = 0, replace 0's with value = max(group size) + 1
    counts.z <- 
      lapply(counts, function(x) which(x == 0))
    
    for (i in 1:length(size)) {
      counts[[i]][counts.z[[i]]] <- max(size) + 1
    }
    
    # Insert into matrices probability of groups of specified sizes for each species. Groups of size  = max(size) + 1 (representing groups of 0) get probability = 1 so these don't contribute to probabilities. 
    prob.mat <-
      lapply(counts, function(y)
        vapply(1:B, function(z)
          size_probability[y[, z] , z], numeric(dim(y)[1])))

    # Use 'index.group' to append group probabilities for each possible true group
    for (i in 1:length(size)) {
      prob.mat[[i]] <- 
        cbind(prob.mat[[i]], group_probability[index.group[[i]]])
    }
    
    # Calculate probability of each true group as row products
    output <-
      lapply(prob.mat, function(x)
        rowprods(x)) %>%
      {
        lapply(1:length(.), function(x)
          setColnames(.[[x]], labels.spp[[x]] ))
      } %>%
      setColnames(., paste(size))
  }
  output
}

# Function: {group.true.probability.f}: With heterogeneous groups (composed of species pairs) for 2 or 3 species, computes probabilities for possible true groups for each unique observation history. Accepts inputs of group probabilities 'group_probability', probability mass of group sizes (1 to maximum observed group size) for each species 'size_probability', and observed group sizes 'size'. Returns a vector with probabilities for each possible true group. For true species states = 3, relies on external 'likelihood.equations' list to provide possible true groups.

group.true.probability.f <- function(group_probability, size_probability, size){
  if (B == 2) {
    # True species states B = 2
    
    # Compute probabilities for homogeneous groups at the start/end of each vector, then add probabilities for heterogeneous groups with decreasing/increasing numbers of species 1/species 2
    homogeneous <-
      vapply(1:B, function(x)
        group_probability[x] * size_probability[size, x], numeric(length(size)))
    
    output <- c(homogeneous[1],
                group_probability[3] * size_probability[0:(size - 1), 2] * size_probability[(size - 1):0, 1],
                homogeneous[2])

  } else if (B == 3) {
    # True species states B = 3
    
    # Add probability = 1 for group size = max(group size) + 1
    size_probability <- rbind(size_probability, 1)
    rownames(size_probability) <- NULL
    
    # Add extra vector elements with group probability = 0 for heterogeneous groups with 3 species (assumed to not occur)
    group_probability <- c(group_probability, 0, 0)
    
    # List 'counts' contains a matrix for each group size with the count of each species in possible true groups (from 'likelihood.equations') on each row
    true_states <- likelihood.equations[[paste0("true.combinations.g.", size)]]
    
    counts <- 
      vapply(1:B, function(y)
        apply(true_states, 1, function(z)
          sum(z == y)), integer(dim(true_states)[1]))
    
    # Counts of each species in each true group
    spp_labels <- apply(counts, 1, function(x) paste0(x, collapse = ""))

    # Each vector in 'index.group' corresponds to matrices in 'counts', with vector elements containing the position of the group probability in 'group_probability' for the true group on each row. Heterogeneous groups with 3 species receive index values corresponding to group probability = 0. 
    index_group <- 
      apply(counts, 1, function(y)
        sum((y > 0) * 1:B) + sum(y > 0) - 1)
      
    # Index locations in each matrix with value = 0, replace 0's with value = max(group size) + 1
    counts_z <- which(counts == 0)
    counts[counts_z] <- dim(size_probability)[1]

    # Insert into matrices probability of groups of specified sizes for each species. Groups of size  = max(size) + 1 (representing groups of 0) get probability = 1 so these don't contribute to probabilities.
    probability_mat <- 
      vapply(1:B, function(z) size_probability[counts[, z] , z], numeric(dim(counts)[1]))

    # Use 'index.group' to append group probabilities for each possible true group
      probability_mat <- 
        cbind(probability_mat, group_probability[index_group])
      
      # Calculate probability of each true group as row products
    output <- 
      apply(probability_mat, 1, prod)
      
    # Add group size and count of each species for true groups
    names(output) <- spp_labels
  }
  output
}

# Function: {group.probability.constant.f}: Calculates group probabilities (pi) for homogeneous or heterogeneous (limited to 3 true species states) groups. For heterogeneous groups, uses the "constant" group probability model described in equation 7 the companion article. Accepts inputs of a vector of true species probabilities 'p', heterogeneous group probabilities 'm', and a vector of unique observed group sizes 'g'. Returns a list with group probabilities 'pi', heterogeneous group probabilities 'm', and a penalty term 'pen' to the -log(likelihood) for violating model constraints.  

group.probability.constant.f <- function(p, m, g) {
  
  # For all models, compute group probabilities (pi) from true species probabilities 'p' (psi) accounting for mean group sizes
  pi <- (p / g) / sum(p / g)
  pen <- 0 # Penalty term for violating model constraints

  # True species states B is 'BB' to avoid conflict with collapse::B 
  BB <- length(p)

  if (any(m > 0)) {
    ## Heterogeneous groups. For clarity, terms here follow equation 7 in companion article
    # Here, 'pi' is "delta" term in companion article (group probabilities prior to heterogeneous groups forming)
    delta <- pi

    if (BB == 2) {
      # True species states B = 2
      
      # Here, heterogeneous group probability 'm' is equivalent to "epsilon"in the companion article
      # Transform heterogeneous group probability 'epsilon' to proportion of all groups (prior to heterogeneous groups forming) that each species contributes to heterogeneous groups, with the resulting value 'epsilon_fraction' equivalent to the "epsilon / (1 + epsilon)" fraction in equation for computing pi for heterogeneous groups in the companion article. 
      epsilon <- m
      epsilon_fraction <- epsilon / (1 + epsilon) 
      
      # For each species, test if heterogeneous group probability 'epsilon_fraction' exceeds 'delta' (i.e. heterogeneous group probabilities for a species exceed overall group probability for that species). If so, reduce heterogeneous group probabilities and compute penalty term for -log(likelihood).
      
      if (sum(delta < epsilon_fraction) > 0) { 
        test <- epsilon_fraction - delta
        delta[which(test > 0)] <- delta[which(test > 0)] + max(test) + 0.01
        delta[which(test < 0)] <- 1 - delta[which(test > 0)]
        test_diff <- test[which(test > 0)] + 1
        for (i in 1:length(test_diff)) {
          pen <- pen + penalty.f(test_diff[i])
        }
      }
      
      # Compute group probabilities (pi) for each species (pi.1, pi.2) and for heterogeneous groups (pi.12)
      pi <- c((delta - epsilon_fraction) / (1 - epsilon_fraction), epsilon)
      
    } else if (BB == 3) {
      # True species states B = 3
      
      # Compute sum of heterogeneous group probabilities 'epsilon' and sums of heterogeneous group probabilities for each species 'spp_sum', which is equivalent to the summation term for pi[bx] in the companion article
      epsilon <- sum(m)
      pairs <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE)
      spp_sum <- vapply(1:3, function(x)
        sum(m[pairs[x, ]]), numeric(1))
      
      # For each species, test if summed heterogeneous group probabilities for a species 'spp.sum' exceeds 'delta' (i.e, heterogeneous group probabilities for a species exceed overall group probability for that species). If so, reduce heterogeneous groups and apply penalty term to -log(likelihood).
      
      test <- delta - spp_sum / (1 + epsilon)
      if (any(test < 0)) {
        xit <- 1
        repeat {
          diff_max <- which.min(test)
          r <- (sum(m[pairs[diff_max, ]]) - (-test[diff_max] + 0.02)) / sum(m[pairs[diff_max, ]])
          m[pairs[diff_max, ]] <- m[pairs[diff_max, ]] * max(r, 0)
          
          epsilon <- sum(m)
          spp_sum <- vapply(1:3, function(x)
            sum(m[pairs[x, ]]), numeric(1))
          test <- delta - spp_sum / (1 + epsilon)
          if (all(test > 0)) break
          xit <- xit + 1
          if (xit > 10) break
        }
      }
      # Compute group probability (pi) for each species (pi.1, pi.2, pi.3) and heterogeneous groups (pi.12, pi.13, pi.23).
      pi <- c((delta - spp_sum / (1 + epsilon)) / (1 - epsilon / (1 + epsilon)), m)
    }
  }
  return(list(pi, m, pen))
}

# Function: {group.probability.psi.constant.f}: With homogeneous or heterogeneous groups (limited to 3 species) and a covariate predicting true species probabilities (psi), computes group probabilities (pi) for each unique observation history (assuming the "constant" model for group probabilities described in eq. 7 in the companion article). Accepts inputs of a matrix 'p' of true species probabilities 1 to B for each unique observation history (rows), heterogeneous group probabilities 'm', and a vector of unique sizes for observed groups 'g'. Returns a list with a matrix of group probabilities for each observation history (rows), heterogeneous group probabilities, and a penalty term to be applied to the -log(likelihood) for violating model constraints.

group.probability.psi.constant.f <- function(p, m, g) {
  
  # For all models, compute group probabilities (pi) for each unique observation history from true species probabilities 'p' (psi) accounting for mean group sizes
  if (any(g > 1)) {
    p <- t(t(p) / g)
    pi <- p / rowsums(p)
  }else{
    pi <- p
  }
  pen <- 0 # Penalty term for violating model constraints
  
  if (any(m > 0)) {
    ## Heterogeneous groups using "constant" model. For clarity, terms follow equation 7 in companion article.
    # Here, 'pi' is "delta" term in companion article (group probabilities prior to heterogeneous groups forming)
    delta <- pi
    
    if (B == 2) {
      # True species states = 2
      
      # Transform heterogeneous group probability 'epsilon' to proportion of all groups (prior to heterogeneous groups forming) that each species contributes to heterogeneous groups, with the resulting value 'epsilon_fraction' equivalent to the "epsilon / (1 + epsilon)" fraction in equation 7 in the companion article for computing pi for heterogeneous groups. 
      epsilon <- m
      epsilon_fraction <- epsilon / (1 + epsilon) 
      
      # For each unique observation history, compute group probability (pi) for each species (pi.1, pi.2) and for heterogeneous groups (pi.12)
      pi <- cbind((delta - epsilon_fraction), epsilon_fraction) / (1 - epsilon_fraction)
      
      # For each species, test if any group probabilities 'pi' are <0. If so, adjust 'pi' to >0 and compute penalty for -log(likelihood). 
      if (any(pi < 0)) {
        negative_values <- which(pi < 0)
        cat("Warning: Mixing probability of ", m," reduced because", length(negative_values), "probabilities of Psi for individual groups were <0", "\n")
        pen <- pen + (1 + -sum(pi[negative_values]) * 10) ^ 2
        negative_rows <- negative_values %% dim(pi)[1]
        negative_rows[which(negative_rows == 0)] <- dim(pi)[1]
        p[negative_rows, dim(pi)[2]] <- 
          p[negative_rows, dim(pi)[2]] + pi[negative_values]
        pi[negative_values] <- 0
      }
    } else if (B == 3) {
      # True species states = 3
      
      # Compute sum of heterogeneous group probabilities 'epsilon' and sums of heterogeneous group probabilities for each species 'spp_sum', which is equivalent to the summation term for pi[bx] in the companion article
      epsilon <- sum(m)
      pairs <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE)
      spp_sum <-
        matrix(vapply(1:3, function(x)
          sum(m[pairs[x,]]), numeric(1)),
          ncol = 3,
          nrow = dim(delta)[1],
          byrow = TRUE)
      
      # For each species, test if summed heterogeneous group probabilities for a species 'spp.sum' exceeds 'delta' (i.e, heterogeneous group probabilities for a species exceed overall group probability for that species). If so, reduce heterogeneous groups and compute penalty term for -log(likelihood).

      test <- delta - spp_sum / (1 + epsilon)
      if (any(test < 0)) {
        xit <- 1
        repeat {
          diff_max <- which.min(test)
          diff_max_col <- (diff_max %% dim(test)[1]) + 1
          r <- (sum(m[pairs[diff_max_col, ]]) - (-test[diff_max] + 0.02)) / sum(m[pairs[diff_max_col, ]])
          m[pairs[diff_max_col, ]] <- m[pairs[diff_max_col, ]] * max(r, 0)
          
          epsilon <- sum(m)
          spp_sum <- 
            matrix(vapply(1:3, function(x)
              sum(m[pairs[x,]]), numeric(1)),
              ncol = 3,
              nrow = dim(delta)[1],
              byrow = TRUE)
          test <- delta - spp_sum / (1 + epsilon)
          if (all(test > 0)) break
          xit <- xit + 1
          if (xit > 10) break
        }
      }
      # Compute group probabilities (pi) for each species (pi.1, pi.2, pi.3) and heterogeneous groups (pi.12, pi.13, pi.23)
      pi <- cbind(
        (delta - spp_sum / (1 + epsilon)) / (1 - epsilon / (1 + epsilon)),
        matrix(
          m,
          ncol = 3,
          nrow = dim(delta)[1],
          byrow = TRUE))
    }
  }
  return(list(pi, m, pen))
}

# Function: {group.probability.encounter.pro.f}: Calculates group probabilities (pi) for homogeneous or heterogeneous (limited to 3 true species states) groups. For heterogeneous groups, uses the "encounter" group probability model described in the companion article in Appendix S3: Eqs. S3, S4. Accepts inputs of a matrix of true species probabilities 'p' for individual observation records (rows), heterogeneous group parameters 'm', and a vector of unique observed group sizes 'g'. Returns a list with group probabilities, heterogeneous group parameters, and a penalty term to the -log(likelihood) for violating model constraints.
group.probability.psi.encounter.f <- function(p, m, g) {
  
  # For all models, compute group probabilities (pi) for each unique observation history from true species probabilities 'p' (psi) accounting for mean group sizes
  if (any(g > 1)) {
    p <- t(t(p) / g)
    pi <- p / rowsums(p)
  }else{
    pi <- p
  }
  pen <- 0 # Penalty term for violating model constraints
  
  if (any(m > 0)) {
    ## Group-level group probabilities with heterogeneous groups using "encounter" model. For clarity, terms follow Appendix S3: Eqs. S3, S4 in the companion article
    # Here, 'pi' is "delta" term in Appendix S3 (group probabilities prior to heterogeneous groups forming)
    delta <- pi
    
    if (B == 2) {
      # True species states = 2
      
      # For all groups (prior to mixing), compute matrix of probabilities 'zeta' of species forming heterogeneous groups given true species probabilities for unique observation histories (row)
      zeta <- m[1] * delta[, 1] * delta[, 2]
      
      # Compute group probabilities for each species (pi.1, pi.2) and heterogeneous groups (pi.12) for unique observation histories (row)
      pi <-  cbind(delta - zeta, zeta) / (1 - zeta)

      # Deprecated
      # For each species, test if any group probabilities 'pi' are <0. If so, adjust 'pi' >0 and compute penalty term for -log(likelihood)
      # if (any(pi < 0)) {
      #   negative.values <- which(pi < 0)
      #   cat("Warning: Mixing probability of ", m," reduced because", length(negative.values), "probabilities of Psi for individual groups were <0", "\n")
      #   pen <- pen + (1 + -sum(p[negative.values]) * 10) ^ 2
      #   negative.rows <- negative.values %% nrow(pi)
      #   negative.rows[which(negative.rows == 0)] <- nrow(pi)
      #   p[negative.rows, ncol(pi)] <-
      #     p[negative.rows, ncol(pi)] + pi[negative.values]
      #   pi[negative.values] <- 0
      # }
    } else if (B == 3) {
      # True species states = 3
      
      # For all groups (prior to mixing), compute matrix of group probabilities 'zeta' for homogeneous groups of each species (cols 1-3; pi.1, pi.2, pi.3) and for heterogeneous groups (cols 4-6; pi.12, pi.13, pi.23) given true species probabilities for each unique observation record (rows)
      zeta <- matrix(0, dim(delta)[1], B + choose(B, 2))
      pairs <- matrix(c(
        1, 2, 1, 3, 2, 3
      ), ncol = 2, byrow = TRUE)
      
      zeta[, 4:6] <- vapply(1:3, function(x)
        m[x] * delta[, pairs[x, 1]] * delta[, pairs[x, 2]]
        , numeric(dim(zeta)[1]))
      
      zeta[, 1:3] <- vapply(1:3, function(x)
        delta[, x] -
          rowsums(zeta[, B + pairs[x, ]])
        , numeric(dim(zeta)[1]))
      
      # Calculate group probabilities (pi) for each species (pi.1, pi.2, pi.3) and heterogeneous groups (p.12, pi.13, pi.23) for each unique observation record
      pi <- zeta / (1 - rowsums(zeta[, 4:6]))

      # epsilon <- sum(m)
      # pairs <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE)
      # spp.sum <- vapply(1:3, function(x)
      #   sum(m[pairs[x, ]]), numeric(1))
      # 
      # test <- p - spp.sum / (1 + epsilon)
      # if (any(test < 0)) {
      #   xit <- 1
      #   repeat {
      #     diff.max <- which.min(test)
      #     r <- (sum(m[pairs[diff.max, ]]) - (-test[diff.max] + 0.02)) / sum(m[pairs[diff.max, ]])
      #     m[pairs[diff.max, ]] <- m[pairs[diff.max, ]] * max(r, 0)
      # 
      #     epsilon <- sum(m)
      #     spp.sum <- vapply(1:3, function(x)
      #       sum(m[pairs[x, ]]), numeric(1))
      #     test <- p - spp.sum / (1 + epsilon)
      #     if (all(test > 0)) break
      #     xit <- xit + 1
      #     if (xit > 10) break
      #   }
      #   
      #   
      #   # cat("mixing greater than group proportion for spp", which(test < 0))
      #   
      # }
      # p <- c((p - spp.sum / (1 + epsilon)) / (1 - epsilon / (1 + epsilon)), m)
      
    }
  }
  return(list(pi, m, pen))
}


###############################################################################
#             Functions for generating simulation inputs
###############################################################################

# Function: {theta.calc.f} Closure function for outputting functions that produce integer matrices (integers specify matrix position from top left by columns) used to assign parameter values to matrices containing classification probabilities for an observer (e.g. 'theta_mat' arrays in generated by optimization functions in 'optimization_functions.R'). Accepts inputs of the maximum possible number of true species states 'B.max' and number of true species states 'B' and observation states 'A' for the current model. Returns a function that generates matrices containing classification probabilities for the specified 'A' and 'B' states. 

theta.calc.f <- function(B.max, B , A) {
  
  # Bring A and B into the function environment
  tmp <- c(A, B)
  
  int.mat <- list(NULL)
  int.mat[[1]] <- matrix(rep(1:B.max, 2), ncol = 2)
  for (i in 1:(B.max - 1)) {
    tmp <- as.matrix(rev(expand.grid(1:(i + 1), 1:(i + 1))))
    dimnames(tmp) <- NULL
    int.mat[[i + 1]] <- tmp
    tmp <- as.matrix(rev(expand.grid(1:(i + 2), 1:(i + 1))))
    dimnames(tmp) <- NULL
    int.mat[[i + B.max]] <- tmp
  }
  int.mat[-1] <-
    lapply(int.mat[-1], function(x) x[-which(x[, 1] == x[, 2]),])
  
  # Output function that accepts inputs of number of observers 'obs' and matrices of classification probability parameters (i.e. without probabilities of correct identification) 'mat', and returns matrices of classification probabilities for each observer with B rows and A columns. 
  function(obs, mat) {

    theta <- matrix(0, B, A)
    # Add correct classification probabilities
    theta[int.mat[[1]][1:B, ]] <- 1 - colSums(matrix(mat[, , obs], ncol = B))
    
    # Add uncertain identification probabilities
    theta[int.mat[[B + (B.max - 1) * (A - B)]]] <- mat[, , obs, drop = FALSE] 
    theta    
  }
}

# Function: {generate.sim.profiles.f} Combines simulation inputs into single data frame with one row for each distinct model. Accepts inputs of simulation parameters common across models 'const' and simulation parameters unique to each model 'var', returns formatted data frame. 

generate.sim.profiles.f <- function(const, var) {
  
  qTBL(const[1:length(const)]) %>%
    {add_vars(.[rep(1, dim(var)[1]), ], var)}
} 

# Function: {sim.const.f} Formats list of simulation parameters 'sim.const' common across models, removing empty and blank entries.

format.sim.constant.f <- function(sim.const) {
  for (i in 1:10) 
  {sim.const[[i]] <- sim.const[[i]][1]}
  # Group probability model defaults to "constant" model described in eq. 7 in the companion article
  if (sim.const$mx_model != "constant" & sim.const$mx_model != "encounter" & sim.const$mx_model != "proportional")
    sim.const$mx_model = "constant"
  return(sim.const)
}

# Function: {generate.ini.f} Generates alternate sets of initial parameter values for model optimization. Accepts inputs of simulation parameters for the current distinct model 'profiles' and vector 'par' with numbers of each parameter type (regression intercept, regression slope, logit, group size, heterogeneous group). Returns matrix 'par_ini' with alternative sets of initial parameter values on each of 5 rows. 

generate.ini.f <- function(profiles, par) {

  # Extract parameters from the simulation profile
  B <- profiles$B; A <- profiles$A
  g.spp <- profiles[grep("^g_", names(profiles))]
  mix <- profiles[grep("mix", names(profiles))]
  psi <- length(grep("psi_", names(profiles)))
  p <- par[3] - psi
  
  par_ini <- matrix(numeric(1), 5, sum(par),
                    dimnames = list(c(), parameters_names))
  
  # Multinomial regression coefficients
  if (length(grep("theta|psi", profiles$Model)) > 0) {
    par_ini[1, 1:par[1]] <- c(seq(-1.5, -2, length.out = par[1]))
    par_ini[2, 1:par[1]] <- c(seq(-2, -4, length.out = par[1]))
    par_ini[3, 1:par[1]] <- c(rep(-3, par[1]))
    par_ini[4, 1:par[1]] <- c(seq(-1, -5, length.out = par[1]))
    par_ini[5, 1:par[1]] <- c(seq(-5, -1, length.out = par[1]))
    
    par_ini[1, (par[1] + 1):sum(par[1:2])] <- c(rep(0, par[2]))
    par_ini[2, (par[1] + 1):sum(par[1:2])] <- c(seq(-0.5, 0.5, length.out = par[2]))
    par_ini[3, (par[1] + 1):sum(par[1:2])] <- c(rep(1, par[2]))
    par_ini[4, (par[1] + 1):sum(par[1:2])] <- c(seq(1, -1, length.out = par[2]))
    par_ini[5, (par[1] + 1):sum(par[1:2])] <- c(seq(0.2, -0.3, length.out = par[2]))   
  }
  
  # Multinomial logit parameters
  if (par[3] > 0) {
  par_ini[1, (sum(par[1:2]) + 1):sum(par[1:3])] <-
    c(seq(A / (A ^ 2 + 1), 0.1, length.out = p), 
      pmin(seq(B / (B ^ 2 - B + 1), B / (B ^ 3), length.out = psi), rep(0.4, psi)))
  par_ini[2, (sum(par[1:2]) + 1):sum(par[1:3])] <-
    c(seq(0.01, A / (A ^ 2 + 1), length.out = p),
      pmin(seq(B / (B ^ 2 - B + 1), B / (B ^ 3.5), length.out = psi), rep(0.2, psi)))
  par_ini[3, (sum(par[1:2]) + 1):sum(par[1:3])] <-
    c( rep(c(0.1, 0.02, 0.05), len = p), 
      seq(1 / (A * 3 / A), 1 / (A * 4 / A), length.out = psi))
  par_ini[4, (sum(par[1:2]) + 1):sum(par[1:3])] <-
    c(seq(A / (A ^ 2 + 1), 0.01, length.out = p),
      seq(1 / (A * 3), 1 / (A + 1), length.out = psi))
  par_ini[5, (sum(par[1:2]) + 1):sum(par[1:3])] <-
    c(seq(0.01, (A - 1) / (3 * A), length.out = p), seq(1 / 20, 1 / 10, length.out = psi))
  
  par_ini[, (sum(par[1:2]) + 1):sum(par[1:3])] <- 
    qlogis(par_ini[, (sum(par[1:2]) + 1):sum(par[1:3])])
  }
  
  # Group size parameters
  if (any(g.spp > 1)) {
    par_ini[, (sum(par[1:3]) + 1):sum(par[1:4])] <- 
      matrix(rep(c(1.2, 1.5, 1.3), B * 5), nrow = 5, ncol = B)
  }
  
  # Heterogeneous group parameters
  if (any(mix > 0)) {
    par_ini[, (sum(par[1:4]) + 1):sum(par[1:5])] <- 
      matrix(c(-4.5, -3.5, -4, -5, -3), nrow = 5, ncol = length(mix), byrow = TRUE)
  }
  par_ini
}

###############################################################################
#             Functions for formatting simulation data
###############################################################################

# Function: {format.MOM.data.f} Summarizes and formats simulated survey observations. If applicable, also adds 1) unique key values (to survey data) and a separate key table for unique observed groups for each observer in each simulation replicate and 2) unique key values (to survey data) and a separate key table for unique covariate values (predicting true species probabilities) for all simulation replicates.
# Accepts inputs of simulated survey observations 'data.obs', number of observation states 'A', number of primary/secondary survey observers 'O_ps', heterogeneous group parameters 'mix', and the number of bins for covariate values 'n_bins'. Returns list with 3 elements: 1) data frame with summarized and formatted simulated survey data, 2) data frame with key table for unique observed groups, and 3) data frame with key table for unique covariate values. If not applicable, keys values are not included in simulated survey data and corresponding list elements (2) and (3) for key tables are NULL.

format.MOM.data.f <- function(data.obs, A, O_ps, mix, n_bins) {

  # Group, count, and sort unique observation histories
  var.name <- as.list(names(data.obs))
  
    var.sort <- c("id",
                  "group_size",
                  grep("y_p", var.name, value = TRUE))
    var.sort <- syms(var.sort)
    var.sort.2 <- quo(desc(count))
  
    dat <- group_by_all(data.obs) %>%
      summarise(., count = n(), .groups = "drop") %>%
      arrange(.,!!! var.sort, !! var.sort.2)

    # Return formatted observation histories without keyed table(s) if: (A) group size = 1 and no covariate predicting true species probabilities is present, if (B) continuous predictive covariates are present, or if (C) predictive covariates for true species probabilities and classification probabilities are present.

  if (max(dat$group_size) == 1 & (length(grep("covariate_psi", colnames(dat))) == 0))
    return(list(dat, NULL))
  if (length(grep("covariate_psi", colnames(dat))) > 0) {
    if  (n_bins < 3 | n_bins > 50) {
      return(list(dat, NULL))
    }
  }
  if (length(grep("covariate_psi", colnames(dat))) > 0 & length(grep("covariate_theta", colnames(dat))) > 0) {
    return(list(dat, NA, NA))
  }
   
    if (length(grep("covariate_theta", colnames(dat))) > 1) {
      return(list(dat, NULL))
    }
    
  # Generate list 'tmp' with 2 elements containing data frames with key values for observed groups ('tmp[[1]]') and with a key table for unique observed groups ('tmp[[2]]') 

  if (length(grep("covariate_theta", colnames(dat))) == 1) {
    tmp <- generate.keys.theta.f(dat, sum(O_ps), A)
  }else{
    tmp <- generate.keys.f(dat, sum(O_ps), A)
  }
    
  # Add key values to formatted survey data 
  add_vars(tmp[[1]], "front") <- dat
  
  # Add "B_states" (number of possible combinations of true groups) and 'sum' total count of classifications for each keyed observed group to the key table
  if (any(mix > 0)) {
    
    add_vars(tmp[[2]]) <-
      ftransform(tmp[[2]], g = rowSums(tmp[[2]][1:A])) %>%
      ftransform(.,
                 B_states = sapply(.$g, function(x)
                   as.integer(
                     factorial(x + B - 1) / prod(factorial(x), factorial(B - 1))
                   )),
                 sum = as.integer(rowsums( qM(tmp[[2]][, 1:A]) ))) %>%
      fselect(., c("B_states", "sum"))
    
  }else{
    
    settransform(tmp[[2]], 
                 B_states = rep(B, dim(tmp[[2]])[1]),
                 sum = as.integer(rowsums( qM(tmp[[2]][, 1:A]) )))
  }
  
  # For models with covariate predicting true species probabilities, add key value for unique covariate values across all simulation replicates 
  if (length(grep("covariate_psi", colnames(dat))) > 0) {
    
    unique.cov.psi <-
      funique(dat, cols = "covariate_psi") %>%
      roworder(., "covariate_psi") %>%
      ftransform(., psi_key = 1:dim(.)[1]) %>%
      fselect(., "covariate_psi", "psi_key")
    
    tmp[[1]] <- left_join(tmp[[1]], unique.cov.psi, by = c("covariate_psi")) 

    # Add key table for unique covariate values
    tmp[[3]] <- unique.cov.psi
  }
  return(tmp)
} 

# Function: {generate.keys.f} Generates unique key values for unique observed groups in each simulation replicate and key tables for unique observed groups. Accepts inputs of simulated survey data 'dat', # of survey observers 'obs', and # of observation states 'A'. Returns list with 2 elements: 1) matrix 'out' containing key values for unique observed groups for each observer in each simulation replicate and 2) matrix 'unique.key' with a key table for unique observed groups. 

generate.keys.f <- function(dat, obs, A) {
  
  # Isolate classification data in 'd'
  d <- get_vars(dat, "^y" , regex = TRUE)
  col.n <- names(d)
  
  # Matrix 'col' gives column numbers for classifications of each observer (row)
  col <- t(sapply(1:obs, function(x)
    ((x - 1) * A + 1):(x * A)))
  names(d) <- c(rep(paste0("A_", 1:A), obs))
  
  # Combine all observed groups A columns, adding one observed group with missing data
  # 'unique.key' = Unique observed groups, with missing data getting key = 1
  
  unique.key <-
    bind_rows(lapply(1:obs, function(x)
      d[, col[x, ]]),
      qDF(matrix(
        0L, ncol = A, dimnames = list(c(), paste0("A_", 1:A))
      )) ) %>%
    funique(.) %>%
    arrange_all(.) %>%
    ftransform(., key = 1:dim(.)[1]) %>%
    setColnames(c(paste0('A_', 1:(dim(.)[2] - 1)), "key"))
  
  # Data frame 'out' contains unique key values for unique observed groups within each simulation replicate (indexed by 'id') with columns labeled for each observer
  
  out <-
    do.call(add_vars, lapply(1:obs, function(x)
      left_join(d[, c(col[x, ])], unique.key, by = c(paste0("A_", 1:A))) %>%
        fselect(., key))) 

  key.n <-
    sapply(strsplit(col.n[col[, 1]], "[.]"), function(x)
      paste0(substr(x[1], 1, nchar(x[1]) - 2), "_key"))
  names(out) <- key.n
  
  return(list(out, unique.key))
}

# Function: {generate.keys.theta.f} With a covariate predicting classification probabilities, generates unique key values for unique combinations of observed groups and covariate values and generates a corresponding key table of unique combinations of observed groups and covariate values. Accepts inputs of simulated survey data 'dat', # of survey observers 'obs', and # of observation states 'A'. Returns list with 2 elements: 1) matrix 'out' containing unique key values for unique combinations of observed groups and covariate values and 2) matrix 'unique.key' containing a corresponding key table.

generate.keys.theta.f <- function(dat, obs, A) {

  # Isolate observation data in 'd'.
  
  d <- get_vars(dat, c("^y", "^cov"), regex = TRUE ) %>%
    setColnames(., c(rep(paste0("A_", 1:A), obs), "covariate_theta"))
  
  # Matrix 'col' gives column numbers for classifications of each observer (row)
  col <- t(sapply(1:obs, function(x)
    ((x - 1) * A + 1):(x * A)))
  col.cov <- dim(d)[2]

  # Combine observed groups for all observers in A columns, with an additional column for covariate values. Add one group observation with missing data.
  # unique.key = Unique combinations of observed groups and covariate values. Missing observations get keys = 1 to (# of covariate bins).
  
  unique.key <-
    bind_rows(lapply(1:obs, function(x)
      d[, c(col[x,], col.cov)]),
      qDF(matrix(
        0L, ncol = A, dimnames = list(c(), paste0("A_", 1:A))
      )) ) %>%
    funique(.) %>%
    arrange_all(.) %>%
    ftransform(., key = 1:dim(.)[1]) %>%
    setColnames(c(paste0('A_', 1:(dim(.)[2] - 2)), "covariate_theta", "key"))

  # Data frame 'out' contains unique key values for unique combinations of observed groups and covariate values, with columns labeled for each observer
  
  out <-
    do.call(add_vars, lapply(1:obs, function(x)
      left_join(d[, c(col[x, ], col.cov)], unique.key, by = c(paste0("A_", 1:A), "covariate_theta")) %>% 
        fselect(., key)))

  col.n <- names(dat)[c(1:(obs * A))]
  
  key.n <-
    sapply(strsplit(col.n[col[, 1]], "[.]"), function(x)
      paste0(substr(x[1], 1, nchar(x[1]) - 2), "_key"))
  names(out) <- key.n
  
  return(list(out, unique.key))
}

###############################################################################
#             Functions for summarizing simulation output
###############################################################################

# Function: {model.results.f} For models from simulation analyses, adds standard errors and confidence limits to parameter estimates, and summarizes estimates, standard errors, and confidence interval coverage. Accepts inputs of list of model results 'results', vector of true parameter values 'par.t', count of successful simulation replicates 'r', and vector with number of parameters of each type 'n' (regression intercept, regression slope, logit, group size, heterogeneous group). Returns data frame with 3 matrices containing parameter estimates, estimated standard errors, and confidence interval coverage for parameter estimates. In each matrix, columns correspond to parameters and rows to simulation replicates. 

model.results.f <- function(results, par.t, r, n) {

  # Apply logit transform to multinomial logit true parameter values
  if (sum(n[c(3, 5)]) > 0) {
    logit.col <- NULL
    if ((n[3]) > 0) {
      logit.col <- c((sum(n[1:2]) + 1):sum(n[1:3]))
    }
    if (n[5] > 0) logit.col <- c(logit.col, (sum(n[1:4]) + 1):sum(n))
    par.t[logit.col] <- as.list(qlogis(unlist(par.t[logit.col])))
  }
  
  # Calculate variance and standard errors of estimated parameters from Hessian matrix, using Delta Method to apply inverse logit function to standard error of multinomial logit parameters
  out <-
    results  %>% rowwise() %>% 
    do(par = as.vector(.$par), 
       var = diag(solve(.$hessian)),
       se = (diag(solve(.$hessian))) ^ 0.5,
       se.logit = diag(DeltaMethod(as.vector(.$par), plogis, solve(.$hessian), 0.000001)$variance) ^ 0.5)
  
  # Insert standard errors of multinomial logit parameters 
  if (sum(n[c(3, 5)]) > 0) {
    for (i in 1:dim(out)[1]) {
      out$se[[i]][logit.col] <- out$se.logit[[i]][logit.col]
    }
  }
  
  # Build matrices with parameter estimates, SEs, 95% CI coverage (values of 0/1 if 95% confidence interval doesn't/does include true values)
  out <- 
    out %>% 
    rowwise() %>% 
    do(l = par.t > (.$par - 1.96 * sqrt(.$var)), 
       u = par.t < (.$par + 1.96 * sqrt(.$var))) %>% 
    do(ci.cov = .$l * .$u) %>% 
    bind_cols(out[, "par"], out[, "se"], .) %>% 
    ungroup() %>%
    do(par = (matrix(  unlist(.$par), nrow = r , byrow = TRUE  )),
       se = (matrix(  unlist(.$se), nrow = r , byrow = TRUE  )  ),
       ci.cov = (matrix(unlist(.$ci.cov), nrow = r, byrow = TRUE ) )
       )
  
  # Apply inverse logit transform to multinomial logit parameters
  if (sum(n[c(3, 5)]) > 0) {
    out[[1]][[1]][, logit.col] <- plogis(out[[1]][[1]][, logit.col])
  }
  out
}

# Function: {model.results.het.s.f} For models from simulation analyses with un-modeled heterogeneity among secondary observers, adds standard errors and confidence limits to parameter estimates, and summarizes estimates, standard errors, and confidence interval coverage. Accepts inputs of list of model results 'results', vector of true parameter values 'par.t', count of successful simulation replicates 'r', vector with number of true parameters of each type 'n' (regression intercept, regression slope, logit, group size, heterogeneous group), vector with number of estimated parameters of each type 'n.est', and vector 'par.d' indicating classification parameters in true model but not estimation model. Returns data frame with 3 matrices containing parameter estimates, estimated standard errors, and confidence interval coverage for parameter estimates. In each matrix, columns correspond to parameters and rows to simulation replicates. 

model.results.het.s.f <- function(results, par.t, r, n, n.est, par.d) {
  
  # Apply logit transform to multinomial logit true parameter values
  if (sum(n[c(3, 5)]) > 0) {
    logit.col <- NULL
    logit.col.est <- NULL
    if ((n[3]) > 0) {
      logit.col <- c((sum(n[1:2]) + 1):sum(n[1:3]))
      logit.col.est <- c((sum(n.est[1:2]) + 1):sum(n.est[1:3]))
    }
    if (n[5] > 0) {
      logit.col <- c(logit.col, (sum(n[1:4]) + 1):sum(n))
      logit.col.est <- c(logit.col.est, (sum(n.est[1:4]) + 1):sum(n.est))
      }
    par.t[logit.col] <- as.list(qlogis(unlist(par.t[logit.col])))
  }
  
  # Calculate variance and standard errors of estimated parameters from Hessian matrix, using Delta Method to apply inverse logit function to standard error of multinomial logit parameters
  out <-
    results  %>% rowwise() %>% 
    do(par = as.vector(.$par), 
       var = diag(solve(.$hessian)),
       se = (diag(solve(.$hessian))) ^ 0.5,
       se.logit = diag(DeltaMethod(as.vector(.$par), plogis, solve(.$hessian), 0.000001)$variance) ^ 0.5)
  
  # Insert standard errors of multinomial logit parameters 
  if (sum(n[c(3, 5)]) > 0) {
    for (i in 1:dim(out)[1]) {
      out$se[[i]][logit.col.est] <- out$se.logit[[i]][logit.col.est]
    }
  }
  
  # Add duplicate estimates for additional (un-modeled) secondary observers
  col <- c(1:par.d[3], 
           rep.int(par.d[2]:par.d[3], par.d[1]),
           (par.d[3] + 1):sum(n.est))
  for (i in 1:dim(out)[1]) {
    out$par[[i]] <- out$par[[i]][col]
    out$var[[i]] <- out$var[[i]][col]
    out$se[[i]] <- out$se[[i]][col]
  }
  
  
  # Build matrices with parameter estimates, SEs, 95% CI coverage (values of 0/1 if 95% confidence interval doesn't/does include true values)
  
  out <- 
    out %>% 
    rowwise() %>% 
    do(l = par.t > (.$par - 1.96 * sqrt(.$var)), 
       u = par.t < (.$par + 1.96 * sqrt(.$var))) %>% 
    do(ci.cov = .$l * .$u) %>% 
    bind_cols(out[, "par"], out[, "se"], .) %>% 
    ungroup() %>%
    do(par = (matrix(  unlist(.$par), nrow = r , byrow = TRUE  )),
       se = (matrix(  unlist(.$se), nrow = r , byrow = TRUE  )  ),
       ci.cov = (matrix(unlist(.$ci.cov), nrow = r, byrow = TRUE ) )
    )
  
  # Apply inverse logit transform to multinomial logit parameters
  if (sum(n[c(3, 5)]) > 0) {
    out[[1]][[1]][, logit.col] <- plogis(out[[1]][[1]][, logit.col])
  }
  out
}

# Function: {summarise.results.f} Summarizes model output for each distinct model across successful replicates. Accepts inputs of data frame parameter output for each simulation replicate 'parameters', a vector of true values of estimated parameters 'par.t', and the count of successful simulation replicates 'r'. Returns a vector with (for each estimated parameter): means, standard deviations, mean standard errors, mean CV (coefficient of variation = SE of estimate / mean of estimate), mean error (parameter estimate - true parameter value), root mean square error, and 95% confidence interval coverage. 

summarise.results.f <- function(parameters, par.t, r){
  
  # For all estimated parameters, compute means, standard deviations, mean standard errors, and 95% confidence interval coverage
  
  summary <-
    parameters %>%
    do(
      mean = colMeans(parameters$par[[1]], na.rm = TRUE),
      sd = aaply(parameters$par[[1]], 2, sd, na.rm = TRUE),
      se = colMeans(parameters$se[[1]], na.rm = TRUE),
      ci.cov = aaply(parameters$ci.cov[[1]], 2, function(x) 
        sum(x, na.rm = TRUE) / length(which(x < 2)))
    )
  
  # Add mean error, root mean square error 
  
  unlist(
    c(
      summary$mean,
      summary$sd,
      summary$se,
      summary$sd[[1]] / summary$mean[[1]],
      summary$mean[[1]] -  par.t,
      colSums((as.matrix(parameters$par[[1]]) -  
                 matrix(unlist(par.t), byrow = TRUE, nrow = r, ncol = length(par.t))) ^ 2) / (r - 1),
      summary$ci.cov,
      r
    )
  )
}

# Function: {psi.est.f} Estimates realized mean true species probabilities (psi) for models with a covariate predicting true species probabilities and up to 3 true species states and observation states. Accepts inputs of vector of estimated parameter values 'par', simulated survey data 'dat', and vector with number of parameters of each type 'n_parameters' (regression intercept, regression slope, logit, group size, heterogeneous group). Returns vector of estimated realized mean true species probabilities. 

psi.est.f <- function(par, dat, n_parameters) {
  # ----- True species states = 2 -----
  if (B == 2) {
   
    if (n_parameters[4] == 0) {
      # Group size = 1

      # Estimate mean realized true species probability (psi) across covariate values
      
      psi.mean <-
        dat %>%
        ftransform(., 
                   psi.1 = count * mlogit.regress.predict.f(covariate_psi , par[1:2])[, 1]) %>%
        fselect(., psi.1, count) %>%
        fsum(., , drop = FALSE) %>%
        fcompute(., 
                 psi.1 = psi.1 / count,
                 psi.2 = 1 - (psi.1 / count))
      
    } else if (n_parameters[5] == 0) {
      # Homogeneous groups
      
      # Estimated mean group size for species 1 to B
      par.g <- par[1:2 + sum(n_parameters[1:3])]
      
      # Estimate mean realized true species probability (psi) across covariate values
      
      psi.mean <- 
        dat %>%
        ftransform(., 
                   psi.1 = mlogit.regress.predict.f(covariate_psi , par[1:2])[, 1]) %>%
        ftransform(.,
                   psi.2 = 1 - psi.1) %>%
        ftransform(., 
                   delta.1 = psi.1 / par.g[1], 
                   delta.2 = psi.2 / par.g[2]) %>%
        fcompute(., 
                 prob.g.1 = (delta.1 / (delta.1 + delta.2)) * count, 
                 prob.g.2 = (delta.2 / (delta.1 + delta.2)) * count
        ) %>% 
        fmean(., , drop = FALSE) %>%
        fcompute(., 
                 psi.1 = prob.g.1 * par.g[1] / (prob.g.1 * par.g[1] + prob.g.2 * par.g[2] ),
                 psi.2 = 1 - (prob.g.1 * par.g[1] / (prob.g.1 * par.g[1] + prob.g.2 * par.g[2]) ))
         
    } else {
      # Heterogeneous groups 
      
      # Estimated mean group size for species 1 to B
      par.g <- par[1:2 + sum(n_parameters[1:3])]
      
      # Estimated heterogeneous group affinity parameter 'rho' for "encounter" heterogeneous group model
      par.mix <- par[length(par)]
      
      # Estimate mean realized true species probability (psi) across covariate values
      
      psi.mean <- 
        dat %>%
        ftransform(., 
                   psi.1 = mlogit.regress.predict.f(covariate_psi , par[1:2])[, 1]) %>%
        ftransform(., 
                   psi.2 = 1 - psi.1) %>%
        ftransform(., 
                   delta.1 = psi.1 / par.g[1], 
                   delta.2 = psi.2 / par.g[2]) %>%
        ftransform(., 
                   prob.g.1 = delta.1 / (delta.1 + delta.2), 
                   prob.g.2 = delta.2 / (delta.1 + delta.2)) %>%
        ftransform(., 
                   prob.g.3 = (par.mix * pmin(prob.g.1, prob.g.2)^2 * pmax(prob.g.1, prob.g.2)) / pmin(prob.g.1, prob.g.2)) %>%
        fcompute(., 
                 prob.g.1 = ((prob.g.1 - prob.g.3) / (1 - prob.g.3)) * count,
                 prob.g.2 = ((prob.g.2 - prob.g.3) / (1 - prob.g.3)) * count,
                 prob.g.3 = (prob.g.3 / (1 - prob.g.3)) * count
        ) %>%
        fmean(., , drop = FALSE) %>%
        fcompute(., 
                 psi.1 = ( prob.g.1 * par.g[1] + prob.g.3 * par.g[1] )  / 
                   ( prob.g.1 * par.g[1] + prob.g.2 * par.g[2] + prob.g.3 * sum(par.g) ),
                 psi.2 = 1 -  ( prob.g.1 * par.g[1] + prob.g.3 * par.g[1] )  / 
                   ( prob.g.1 * par.g[1] + prob.g.2 * par.g[2] + prob.g.3 * sum(par.g) ))
    }
  }
  
  # ----- True species states = 3 -----
  if (B == 3) {

    if (n_parameters[4] == 0) {
      # Group size = 1

      # Estimate realized true species probabilities (psi) for each covariate value 'tmp' and overall mean for the model 'psi.mean'
      tmp <- dat$count * mlogit.regress.predict.f(dat$covariate_psi , par[1:4])
      
      psi.mean <- 
        dat %>%
        ftransform(., 
               psi.1 = tmp[, 1],
               psi.2 = tmp[, 2]
        ) %>%
        summarise(.,
                  psi.1 = sum(psi.1) / sum(count), 
                  psi.2 = sum(psi.2) / sum(count),
                  psi.3 = 1 - psi.1 - psi.2
        ) 
    } else if (n_parameters[5] == 0) {
      # Homogeneous groups
      
      # Estimated mean group size for species 1 to B
      par.g <- par[1:3 + sum(n_parameters[1:3])]
      
      # Estimate realized true species probabilities (psi) for each covariate value 'tmp' and overall mean for the model 'psi.mean'
      tmp <- mlogit.regress.predict.f(dat$covariate_psi , par[1:4])
      
      psi.mean <- 
        dat %>%
        mutate(., 
               psi.1 = tmp[, 1],
               psi.2 = tmp[, 2],
               psi.3 = tmp[, 3],
               delta.1 = psi.1 / par.g[1], 
               delta.2 = psi.2 / par.g[2], 
               delta.3 = psi.3 / par.g[3], 
               prob.g.1 = delta.1 / (delta.1 + delta.2 + delta.3), 
               prob.g.2 = delta.2 / (delta.1 + delta.2 + delta.3),
               prob.g.3 = delta.3 / (delta.1 + delta.2 + delta.3)
        ) %>%
        summarise(.,
                  psi.1 = (
                    mean(prob.g.1 * count) * par.g[1]) / 
                    (mean(prob.g.1 * count) * par.g[1] + mean(prob.g.2 * count) * par.g[2] + mean(prob.g.3 * count) * par.g[3]), 
                  psi.2 = (
                    mean(prob.g.2 * count) * par.g[2]) / 
                    (mean(prob.g.1 * count) * par.g[1] + mean(prob.g.2 * count) * par.g[2] + mean(prob.g.3 * count) * par.g[3]), 
                  psi.3 = 1 - psi.1 - psi.2
        ) 
    } else {
      # Heterogeneous groups 
      
      # Estimated mean group size for species 1 to B
      par.g <- par[1:3 + sum(n_parameters[1:3])]
      
      # Estimated mixing parameter 'rho' for encounter mixing model
      par.mix <- par[1:3 + sum(n_parameters[1:4])]
      
      # Estimate realized true species probabilities (psi) for each covariate value 'tmp' and overall group probabilities before mixing 'delta'
      tmp <- mlogit.regress.predict.f(dat$covariate_psi , par[1:4])
      
      delta <- as.matrix(
        dat %>%
          mutate(., 
                 psi.1 = tmp[, 1],
                 psi.2 = tmp[, 2],
                 psi.3 = tmp[, 3]
          ) %>%
          transmute(., 
                    tmp.1 = psi.1 / par.g[1], 
                    tmp.2 = psi.2 / par.g[2], 
                    tmp.3 = psi.3 / par.g[3]
          ) %>%
          transmute(., 
                    delta.1 = tmp.1 / (tmp.1 + tmp.2 + tmp.3),
                    delta.2 = tmp.2 / (tmp.1 + tmp.2 + tmp.3),
                    delta.3 = tmp.3 / (tmp.1 + tmp.2 + tmp.3)
          )
      )
      
      # Calculate matrix of probabilities 'zeta' of homogeneous groups of each species (cols 1-3) and heterogeneous groups (cols 4-6) given delta for each unique observation history (rows)
      zeta <- matrix(0, dim(delta)[1], B + choose(B, 2))
      pairs <- matrix(c(1, 2, 1, 3, 2, 3), 
                      ncol = 2, 
                      byrow = TRUE)
      
      zeta[, 4:6] <- vapply(1:3, function(x)
        par[c(x + sum(n_parameters[1:4]))] *
          (pmin(delta[, pairs[x, 1]], delta[, pairs[x, 2]])) ^ 2 *
          pmax(delta[, pairs[x, 1]], delta[, pairs[x, 2]]) /
          pmin(delta[, pairs[x, 1]], delta[, pairs[x, 2]])
        , numeric(dim(zeta)[1]))
      
      zeta[, 1:3] <- vapply(1:3, function(x)
        delta[, x] -
          rowsums(zeta[, B + pairs[x, ]])
        , numeric(dim(zeta)[1]))
      
      # Calculate group probabilities (pi) for homogeneous and heterogeneous groups for each unique observation history
      psi.mean <- zeta / (1 - rowsums(zeta[, 4:6]))
      colnames(psi.mean) <- c(paste0("prob.g.", 1:B), "prob.g.12", "prob.g.13", "prob.g.23")
      dat <- bind_cols(dat, as_tibble(psi.mean))
      
      # Estimate mean realized true species probability (psi) across covariate values
      psi.mean <-
        summarise(dat,
                  psi.1 = (
                    (mean(prob.g.1 * count) * par.g[1] + mean(prob.g.12 * count) * par.g[1] + mean(prob.g.13 * count) * par.g[1])
                    / (mean(prob.g.1 * count) * par.g[1] + mean(prob.g.2 * count) * par.g[2] + mean(prob.g.3 * count) * par.g[3] +
                         mean(prob.g.12 * count) * sum(par.g[c(1:2)]) + 
                         mean(prob.g.13 * count) * sum(par.g[c(1, 3)]) +
                         mean(prob.g.23 * count) * sum(par.g[c(2:3)]))
                  ),
                  psi.2 = (
                    (mean(prob.g.2 * count) * par.g[2] + mean(prob.g.12 * count) * par.g[2] + mean(prob.g.23 * count) * par.g[2])
                    / (mean(prob.g.1 * count) * par.g[1] + mean(prob.g.2 * count) * par.g[2] + mean(prob.g.3 * count) * par.g[3] +
                         mean(prob.g.12 * count) * sum(par.g[c(1:2)]) + 
                         mean(prob.g.13 * count) * sum(par.g[c(1, 3)]) +
                         mean(prob.g.23 * count) * sum(par.g[c(2:3)]))
                  ),
                  psi.3 = 1 - psi.1 - psi.2
        )
    }
  }
  unlist(psi.mean)
}

# Function: {theta.est.f} Estimates realized mean classification probabilities (theta) for models with a covariate predicting classification probabilities and up to 3 true species states and observation states. Accepts inputs of vector of estimated parameter values 'par' and simulated survey data 'dat'. Returns vector of estimated realized mean classification probabilities. 

theta.est.f <- function(par, dat) {
  
  if (B == 2 & A == 2) {
    # ----- True species states = observation states = 2 -----
    
    if (O_ps[1] > 0 & sim_profiles[1, ]$Model != "M.theta.p") {
      # Multi-observation method (MOM) models 
      
      # Vector 'par.beta' contains estimated intercept and slope regression coefficients 
      if (length(grep("covariate_psi", colnames(dat))) > 0) {
        par.beta <- 
          t(sapply(0:3, function(x) 
            par[c((3 + x), (7 + x))]))
      } else {
        par.beta <- 
          t(sapply(0:3, function(x) 
            par[c((1 + x), (5 + x))]))
      }
      
      # Estimate overall mean for classification probabilities
      theta.mean <-
        dat %>%
        ftransform(.,
               theta.p.21 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[1, ])[, 1],
               theta.p.12 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[2, ])[, 1],
               theta.s.21 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[3, ])[, 1],
               theta.s.12 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[4, ])[, 1]
        ) %>%
        summarise(.,
                  theta.p.21 = sum(theta.p.21) / sum(count * group_size),
                  theta.p.12 = sum(theta.p.12) / sum(count * group_size),
                  theta.s.21 = sum(theta.s.21) / sum(count * group_size),
                  theta.s.12 = sum(theta.s.12) / sum(count * group_size)
        )
    } else {
      # Single-observation method (SOM) models 
      
      # Vector 'par.beta' contains estimated intercept and slope regression coefficients
      if (length(grep("covariate_psi", colnames(dat))) > 0) {
        par.beta <- 
          t(sapply(0:1, function(x) 
            par[c((3 + x), (5 + x))]))
      } else {
        par.beta <- 
          t(sapply(0:1, function(x) 
            par[c((1 + x), (3 + x))]))
      }
      
      # Estimate overall mean for classification probabilities
      theta.mean <-
        dat %>%
        ftransform(.,
               theta.s.21 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[1, ])[, 1],
               theta.s.12 = count * group_size * mlogit.regress.predict.f(covariate_theta , par.beta[2, ])[, 1]
        ) %>%
        summarise(.,
                  theta.s.21 = sum(theta.s.21) / sum(count * group_size),
                  theta.s.12 = sum(theta.s.12) / sum(count * group_size)
        )
      if (sim_profiles[1, ]$Model == "M.theta.p") {
        colnames(theta.mean) <- c("theta_p_21", "theta_p_12")
      }
    }
  } else if (B == 2 & A == 3) {
    # ----- True species states = 2, observation states = 3 -----
    
    if (O_ps[1] > 0 & sim_profiles[1, ]$Model != "M.theta.p") {
      # Multi-observation method (MOM) models 
      
      # Estimate realized classification probabilities (theta) for each covariate value 'tmp' and overall mean for the model 'theta.mean'
      tmp <- 
        matrix(
          sapply(seq(0, 6, by = 2), function(x)
            mlogit.regress.predict.f(
              dat$covariate_theta, 
              matrix(unlist(par[c((1:2 + x), (9:10 + x))]), ncol = 2) )[, c(1:2)]
          ), ncol = 8
        )
      
      theta.mean <- 
        dat %>%
        ftransform(., 
               theta.p.21 = count * tmp[, 1], 
               theta.p.31 = count * tmp[, 2],
               theta.p.12 = count * tmp[, 3], 
               theta.p.32 = count * tmp[, 4],
               theta.s.21 = count * tmp[, 5],
               theta.s.31 = count * tmp[, 6],
               theta.s.12 = count * tmp[, 7],
               theta.s.32 = count * tmp[, 8]
        ) %>%
        summarise(., 
                  theta.p.21 = sum(theta.p.21) / sum(count),
                  theta.p.31 = sum(theta.p.31) / sum(count),
                  theta.p.12 = sum(theta.p.12) / sum(count),
                  theta.p.32 = sum(theta.p.32) / sum(count),
                  theta.s.21 = sum(theta.s.21) / sum(count),
                  theta.s.31 = sum(theta.s.31) / sum(count),
                  theta.s.12 = sum(theta.s.12) / sum(count),
                  theta.s.32 = sum(theta.s.32) / sum(count)
        )
    }else{
      # Single-observation method (SOM) models 
      
      # Estimate realized classification probabilities (theta) for each covariate value 'tmp' and overall mean for the model 'theta.mean'
      tmp <- 
        matrix(
          sapply(seq(0, 2, by = 2), function(x)
            mlogit.regress.predict.f(
              dat$covariate_theta, 
              matrix(unlist(par[c((1:2 + x), (5:6 + x))]), ncol = 2) )[, c(1:2)]
          ), ncol = 4
        )
      
      theta.mean <- 
        dat %>%
        ftransform(., 
               theta.s.21 = count * tmp[, 1],
               theta.s.31 = count * tmp[, 2],
               theta.s.12 = count * tmp[, 3],
               theta.s.32 = count * tmp[, 4]
        ) %>%
        summarise(., 
                  theta.s.21 = sum(theta.s.21) / sum(count),
                  theta.s.31 = sum(theta.s.31) / sum(count),
                  theta.s.12 = sum(theta.s.12) / sum(count),
                  theta.s.32 = sum(theta.s.32) / sum(count)
        )
      if (sim_profiles[1, ]$Model == "M.theta.p") {
        colnames(theta.mean) <- c("theta_p_21", "theta_p_31", "theta_p_12", "theta_p_32")
      }
    }
   
  } else if (B == 3 & A == 3) {
    # ----- True species states = observation states = 3 -----
    
    if (O_ps[1] > 0 & sim_profiles[1, ]$Model != "M.theta.p") {
      # Multi-observation method (MOM) models 
      
      # Estimate realized classification probabilities (theta) for each covariate value 'tmp' and overall mean for the model 'theta.mean'
      tmp <- 
        matrix(
          sapply(seq(0, 10, by = 2), function(x)
            mlogit.regress.predict.f(
              dat$covariate_theta, 
              matrix(unlist(par[c((1:2 + x), (13:14 + x))]), ncol = 2) )[, c(1:2)]
          ), ncol = 12
        )
      
      theta.mean <- 
        dat %>%
        ftransform(., 
               theta.p.21 = count * tmp[, 1], 
               theta.p.31 = count * tmp[, 2],
               theta.p.12 = count * tmp[, 3], 
               theta.p.32 = count * tmp[, 4],
               theta.p.13 = count * tmp[, 5], 
               theta.p.23 = count * tmp[, 6],
               
               theta.s.21 = count * tmp[, 7],
               theta.s.31 = count * tmp[, 8],
               theta.s.12 = count * tmp[, 9],
               theta.s.32 = count * tmp[, 10],
               theta.s.13 = count * tmp[, 11],
               theta.s.23 = count * tmp[, 12]
        ) %>%
        summarise(., 
                  theta.p.21 = sum(theta.p.21) / sum(count),
                  theta.p.31 = sum(theta.p.31) / sum(count),
                  theta.p.12 = sum(theta.p.12) / sum(count),
                  theta.p.32 = sum(theta.p.32) / sum(count),
                  theta.p.13 = sum(theta.p.13) / sum(count),
                  theta.p.23 = sum(theta.p.23) / sum(count),
                  
                  theta.s.21 = sum(theta.s.21) / sum(count),
                  theta.s.31 = sum(theta.s.31) / sum(count),
                  theta.s.12 = sum(theta.s.12) / sum(count),
                  theta.s.32 = sum(theta.s.32) / sum(count),
                  theta.s.13 = sum(theta.s.13) / sum(count),
                  theta.s.23 = sum(theta.s.23) / sum(count)
        )

  }else{
    # Single-observation method (SOM) models 
    
    # Estimate realized classification probabilities (theta) for each covariate value 'tmp' and overall mean for the model 'theta.mean'
    tmp <- 
      matrix(
        sapply(seq(0, 4, by = 2), function(x)
          mlogit.regress.predict.f(
            dat$covariate_theta, 
            matrix(unlist(par[c((1:2 + x), (7:8 + x))]), ncol = 2)  )[, c(1:2)]
        ), ncol = 6
      )
    
    theta.mean <- 
      dat %>%
      ftransform(., 
             theta.s.21 = count * tmp[, 1],
             theta.s.31 = count * tmp[, 2],
             theta.s.12 = count * tmp[, 3],
             theta.s.32 = count * tmp[, 4],
             theta.s.13 = count * tmp[, 5],
             theta.s.23 = count * tmp[, 6]
      ) %>%
      summarise(., 
                theta.s.21 = sum(theta.s.21) / sum(count),
                theta.s.31 = sum(theta.s.31) / sum(count),
                theta.s.12 = sum(theta.s.12) / sum(count),
                theta.s.32 = sum(theta.s.32) / sum(count),
                theta.s.13 = sum(theta.s.13) / sum(count),
                theta.s.23 = sum(theta.s.23) / sum(count)
      )
    if (sim_profiles[1, ]$Model == "M.theta.p") {
      colnames(theta.mean) <- c("theta_p_21", "theta_p_31", "theta_p_12", "theta_p_32", "theta_p_13", "theta_p_23")
    }
  }
  }
  unlist(theta.mean)
}

# Function: {logit.ci.f} Estimates confidence intervals for probability parameters assuming the sampling distribution of the logit-transformed parameter estimates is normal. Accepts inputs of vectors of estimated parameter values 'p' and estimated standard errors 's', returns vector of estimated lower and upper 95% confidence limits for each parameter. 

logit.ci.f <- function(p, s) {
  
  # Apply logit transform to parameter estimates, estimate standard errors for logit(parameters) using Delta Method 
  p.logit <- sapply(1:length(p), function(x) qlogis(p[x]))
  p.logit.se <- 
    sapply(1:length(p), function(x)
      (DeltaMethod(
        p[x],
        qlogis,
        vcov = unlist(s[x] ^ 2),
        delta = 0.0000001
      )$variance) ^ 0.5
    )
  
  # Generate 95% confidence limits assuming asymptotic normality in logit scale
  out <- 
    as.vector(
      sapply(1:length(p), function(x)
        plogis(c(p.logit[x] - 1.96 * p.logit.se[x], p.logit[x] + 1.96 * p.logit.se[x]))
      )
    )
}

# Function: {psi.derived.f} Produces overall mean estimates of true species probabilities (psi) and classification probabilities (theta) derived from regression coefficients (betas) from multinomial logistic regression models along with associated standard errors, 95% confidence limits, and 95% confidence interval coverage. Supports models with regression coefficients predicting true species probabilities (true species states B = 2 or 3) or regression coefficients predicting true species probabilities (B = 2) AND classification probabilities (A = 2).

# Accepts inputs of:
# 'sim_profiles' containing simulation profile for the current distinct model 
# 'sim.data.tmp' containing simulated survey data 
# 'model.results' containing simulation output from model optimization 
# 'parameter.output' containing summarized statistics for estimated parameters 
# 'output.betas' containing summarized statistics for overall means for true species probabilities and for classification probabilities from prior distinct models 

# Returns list of output including;
# output.betas' containing  summarized statistics (mean, root mean square error, 95% CI coverage) across simulation replicates for overall mean true species probabilities and classification probabilities. Contains results for the current and all prior distinct models.
# 'betas' containing (for each simulation replicate) estimated regression coefficients and estimated means, SEs, confidence limits, and 95% confidence interval coverage for overall mean true species probabilities and classification probabilities

psi.derived.f <- function(sim_profiles, sim.data, sim.data.tmp, model.results, parameter.output, output.betas) {
  
  # Optional output for derived estimates of overall mean true species probabilities (psi) and classification probabilities (theta)
  
  # 'tmp.id' and 'n.sim' contain an index and the count of simulation replicates successfully optimizing
  # 'tmp.betas.psi' and 'tmp.betas.theta' contain estimated intercept and slope parameters for regression coefficients predicting true species probabilities (psi) and classification probabilities (theta)
  tmp.id <- unique(sim.data.tmp$id)
  n.sim <- rep(sim, reps)
  tmp.betas.psi <- cbind(parameter.output$par[[1]][, c(grep("psi", parameters_names))])
  colnames(tmp.betas.psi) <- c(grep("psi", parameters_names, value = TRUE))
  if (sim_profiles$Model[1] == "M.theta.psi" & B == 2 & A == 2 | sim_profiles$Model[1] == "M.theta+psi" & B == 2 & A == 2) {
    tmp.betas.theta <- cbind(parameter.output$par[[1]][, c(grep("b0...[^p]|b1...[^p]", parameters_names))])
    colnames(tmp.betas.theta) <- c(grep("b0...[^p]|b1...[^p]", parameters_names, value = TRUE))
  }

    if (B == 2) {
      ## ----- True species states (B) = 2 -----

      # Estimate overall mean true species probabilities (psi) and associated SEs (using Delta Method)
      tmp.psi <-
        t(sapply(1:reps, function(x)
          psi.est.f(parameter.output$par[[1]][x,],
                    sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]),],
                    n_parameters)))
      
      tmp.se <-
        t(sapply(1:reps, function(x)
          diag(
            DeltaMethod(
              parameter.output$par[[1]][x, ],
              psi.est.f,
              vcov = solve(model.results$hessian[[x]]),
              delta = 0.00001,
              dat = sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ],
              n_parameters = n_parameters
            )$variance
          ) ^ 0.5
        ))
      colnames(tmp.se) <- paste0("se.psi.", 1:B)
      
      # Apply logit transformation to estimates and standard errors, add 95% confidence limits assuming asymptotic normal distribution (in logit scale), and apply inverse logit transformation to CLs. Confidence limits from logit transform have '.lg.L' and '.lg.U' suffixes. 
      
      tmp.lg <-
        t(sapply(1:reps, function(x)
          logit.ci.f(tmp.psi[x, 1:B], tmp.se[x, 1:B])))
      colnames(tmp.lg) <- c("psi.1.lg.L", "psi.1.lg.U", "psi.2.lg.L", "psi.2.lg.U")
      
      # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
      
      tmp.psi <-
        qTBL(cbind(tmp.psi, tmp.se, tmp.lg)) %>%
        ftransform(.,
                   psi.1.L = psi.1 - 1.96 * se.psi.1,
                   psi.1.U = psi.1 + 1.96 * se.psi.1,
                   psi.2.L = psi.2 - 1.96 * se.psi.2,
                   psi.2.U = psi.2 + 1.96 * se.psi.2) %>%
        ftransform(.,
                   ci.psi.1 = (psi.1.L < sim.data[[4]][[1]][1]) * (psi.1.U > sim.data[[4]][[1]][1]),
                   ci.psi.2 = (psi.2.L < sim.data[[4]][[1]][2]) * (psi.2.U > sim.data[[4]][[1]][2]),
                   ci.lg.psi.1 = (psi.1.lg.L < sim.data[[4]][[1]][1]) * (psi.1.lg.U > sim.data[[4]][[1]][1]),
                   ci.lg.psi.2 = (psi.2.lg.L < sim.data[[4]][[1]][2]) * (psi.2.lg.U > sim.data[[4]][[1]][2])
        )
      
      if (sim_profiles$Model[1] == "M.theta.psi" & B == 2 & A == 2 | sim_profiles$Model[1] == "M.theta+psi" & B == 2 & A == 2) {
        ## ----- Models also includes covariate predicting classification probabilities (theta) ----
        
        # Estimate overall mean classification probabilities and associated SEs (using Delta Method)
        tmp.theta <- 
          t(sapply(1:reps, function(x)
            theta.est.f(parameter.output$par[[1]][x, ],  sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ])
          ))
        
        tmp.se <- 
          t(sapply(1:reps, function(x)
            diag(
              DeltaMethod(
                parameter.output$par[[1]][x, ],
                theta.est.f,
                vcov = solve(model.results$hessian[[x]]),
                delta = 0.00001,
                dat = sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ]
              )$variance
            )) ^ 0.5
          )
        colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.12", "se.theta.s.21", "se.theta.s.12")[c(O_ps[1]*1:2, 3:4)]
        
        # Apply logit transformation to estimates and standard errors, add 95% confidence limits assuming asymptotic normal distribution (in logit scale), and apply inverse logit transformation to CLs. Confidence limits from logit transform have '.lg.L' and '.lg.U' suffixes. 
        tmp.lg <- 
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, ]), unlist(tmp.se[x, ]))))
        
        colnames(tmp.lg) <- c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", 
                              "theta.s.21.lg.L", "theta.s.21.lg.U", "theta.s.12.lg.L", "theta.s.12.lg.U")[c(O_ps[1]*1:4, 5:8)]
        
        if (O_ps[1] > 0) {
          # ----- Multi-observation method (MOM) models (primary observes present) -----
          
          # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)

          tmp.theta <- 
            qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
            ftransform(., 
                       theta.p.21.L = theta.p.21 - 1.96 * se.theta.p.21,
                       theta.p.21.U = theta.p.21 + 1.96 * se.theta.p.21,
                       theta.p.12.L = theta.p.12 - 1.96 * se.theta.p.12,
                       theta.p.12.U = theta.p.12 + 1.96 * se.theta.p.12, 
                       theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                       theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                       theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                       theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12) %>%
            
            ftransform(., 
                       
                       ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1]) * (theta.p.21.U > sim.data[[4]][[2]][1]),
                       ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][2]) * (theta.p.12.U > sim.data[[4]][[2]][2]),
                       ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][3]) * (theta.s.21.U > sim.data[[4]][[2]][3]),
                       ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][4]) * (theta.s.12.U > sim.data[[4]][[2]][4]),
                       
                       ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1]),
                       ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][2]),
                       ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][3]) * (theta.s.21.lg.U > sim.data[[4]][[2]][3]),
                       ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][4]) * (theta.s.12.lg.U > sim.data[[4]][[2]][4])
            )
          
          # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean true species probabilities and classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
          
          betas <- rbind(betas, 
                         cbind(tmp.betas.psi, 
                               tmp.betas.theta, 
                               as.matrix(tmp.psi), 
                               psi.1.t = sim.data[[4]][[1]][1], 
                               psi.2.t = sim.data[[4]][[1]][2], 
                               as.matrix(tmp.theta), 
                               theta.p.21.t = sim.data[[4]][[2]][1], 
                               theta.p.12.t = sim.data[[4]][[2]][2], 
                               theta.s.21.t = sim.data[[4]][[2]][3], 
                               theta.s.12.t = sim.data[[4]][[2]][4],
                               n.sim = n.sim)
          )
          
          # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
          
          output.betas <- 
            add_vars(tmp.psi, tmp.theta) %>%
            ftransform(., 
                       psi.1.err = psi.1 - sim.data[[4]][[1]][[1]],
                       theta.p.21.err = theta.p.21 - sim.data[[4]][[2]][[1]],
                       theta.p.12.err = theta.p.12 - sim.data[[4]][[2]][[2]],
                       theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][[3]],
                       theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][[4]]
            ) %>%
            summarise(., 
                      me.psi.1 = mean(psi.1.err), 
                      me.theta.p.21 = mean(theta.p.21.err),
                      me.theta.p.12 = mean(theta.p.12.err),
                      me.theta.s.21 = mean(theta.s.21.err),
                      me.theta.s.12 = mean(theta.s.12.err),
                      
                      rm.psi.1 = (sum(psi.1.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                      rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      
                      ci.lg.psi.1 = mean(ci.lg.psi.1),
                      ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE),
                      ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE),  
                      ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE),
                      ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE)
            ) %>%
            bind_rows(output.betas, .)
          
          # Return list with output
          return(list(output.betas, betas))
          
        }else { 
          # ----- Single-observation method (SOM) models (primary observes NOT present) -----
          
          # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
          
          tmp.theta <- 
            qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
            ftransform(., 
                       theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                       theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                       theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                       theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12) %>%
            
            ftransform(., 
                       
                       ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1]) * (theta.s.21.U > sim.data[[4]][[2]][1]),
                       ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][2]) * (theta.s.12.U > sim.data[[4]][[2]][2]),
                       
                       ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1]),
                       ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][2]) * (theta.s.12.lg.U > sim.data[[4]][[2]][2])
            )
          
          # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean true species probabilities and classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
          
          betas <- rbind(betas, 
                         cbind(tmp.betas.psi, 
                               tmp.betas.theta, 
                               as.matrix(tmp.psi), 
                               psi.1.t = sim.data[[4]][[1]][1], 
                               psi.2.t = sim.data[[4]][[1]][2], 
                               as.matrix(tmp.theta), 
                               theta.s.21.t = sim.data[[4]][[2]][1], 
                               theta.s.12.t = sim.data[[4]][[2]][2], 
                               n.sim = n.sim)
          )
          
          # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
          
          output.betas <- 
            add_vars(tmp.psi, tmp.theta) %>%
            ftransform(., 
                       psi.1.err = psi.1 - sim.data[[4]][[1]][[1]],
                       theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][[1]],
                       theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][[2]]
            ) %>%
            summarise(., 
                      me.psi.1 = mean(psi.1.err), 
                      me.theta.s.21 = mean(theta.s.21.err),
                      me.theta.s.12 = mean(theta.s.12.err),
                      
                      rm.psi.1 = (sum(psi.1.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                      rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                      
                      ci.lg.psi.1 = mean(ci.lg.psi.1),
                      ci.lg.theta.s.21 = mean(ci.lg.theta.s.21),
                      ci.lg.theta.s.12 = mean(ci.lg.theta.s.12)
                      
            ) %>%
            bind_rows(output.betas, .)
          
          # Return list with output
          return(list(output.betas, betas))
        }  
      } else {
        ## ----- Model only includes covariate predicting true species probabilities -----
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean true species probabilities and classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- rbind(betas, 
                       cbind(tmp.betas.psi, 
                             as.matrix(tmp.psi), 
                             psi.1.t = sim.data[[4]][[1]][[1]], 
                             psi.2.t = sim.data[[4]][[1]][[2]], 
                             n.sim = n.sim)
        )
        
        # Optional: also un-comment code in next assignment
        # psi.group <-
        #   ldply(1:reps, function(x)
        #     psi.est.f(
        #                     parameter.output$par[[1]][x,],
        #                     sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]),],
        #                     n_parameters))
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          mutate(tmp.psi, 
                 psi.1.err = psi.1 - sim.data[[4]][[1]][[1]]) %>%
          summarise(., 
                    me.psi.1 = mean(psi.1.err), 
                    rm.psi.1 = (sum(psi.1.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                    ci.lg.psi.1 = mean(ci.lg.psi.1, na.rm = TRUE)
          ) %>%
          # bind_cols(., summarise_all(psi.group, mean)) %>%
          bind_rows(output.betas, .)
      }
      
      # Return list of estimates
      return(list(output.betas, betas))
      
    }else if (B == 3) {
      ## ----- True species states (B) = 3 -----
      # Models with B = 3 and regression coefficients predicting true species probabilities AND classification probabilities not supported
      
      # Estimate overall mean true species probabilities (psi) and associated SEs (using Delta Method)
      tmp.psi <-
        t(sapply(1:reps, function(x)
          psi.est.f(parameter.output$par[[1]][x,], 
                    sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]),], 
                    n_parameters)))
      
      tmp.se <- 
        t(sapply(1:reps, function(x)
          diag(
            DeltaMethod(
              parameter.output$par[[1]][x, ],
              psi.est.f,
              vcov = solve(model.results$hessian[[x]]),
              delta = 0.00001,
              dat = sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ],
              n_parameters = n_parameters
            )$variance
          )) ^ 0.5
        )
      colnames(tmp.se) <- paste0("se.psi.", 1:B)
      
      # Apply logit transformation to estimates and standard errors, add 95% confidence limits assuming asymptotic normal distribution (in logit scale), and apply inverse logit transformation to CLs. Confidence limits from logit transform have '.lg.L' and '.lg.U' suffixes. 
      
      tmp.lg <- 
        t(sapply(1:reps, function(x)
          logit.ci.f(tmp.psi[x, 1:B], tmp.se[x, 1:B])))
      colnames(tmp.lg) <- c("psi.1.lg.L", "psi.1.lg.U", "psi.2.lg.L", "psi.2.lg.U", "psi.3.lg.L", "psi.3.lg.U")
      
      # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
      
      tmp.psi <-
        qTBL(cbind(tmp.psi, tmp.se, tmp.lg)) %>%
        ftransform(.,
                   psi.1.L = psi.1 - 1.96 * se.psi.1,
                   psi.1.U = psi.1 + 1.96 * se.psi.1,
                   psi.2.L = psi.2 - 1.96 * se.psi.2,
                   psi.2.U = psi.2 + 1.96 * se.psi.2,
                   psi.3.L = psi.3 - 1.96 * se.psi.3,
                   psi.3.U = psi.3 + 1.96 * se.psi.3) %>%
        
        ftransform(.,
                   
                   ci.psi.1 = (psi.1.L < sim.data[[4]][[1]][1]) * (psi.1.U > sim.data[[4]][[1]][1]),
                   ci.psi.2 = (psi.2.L < sim.data[[4]][[1]][2]) * (psi.2.U > sim.data[[4]][[1]][2]),
                   ci.psi.3 = (psi.3.L < sim.data[[4]][[1]][3]) * (psi.3.U > sim.data[[4]][[1]][3]),
                   ci.lg.psi.1 = (psi.1.lg.L < sim.data[[4]][[1]][1]) * (psi.1.lg.U > sim.data[[4]][[1]][1]),
                   ci.lg.psi.2 = (psi.2.lg.L < sim.data[[4]][[1]][2]) * (psi.2.lg.U > sim.data[[4]][[1]][2]),
                   ci.lg.psi.3 = (psi.3.lg.L < sim.data[[4]][[1]][3]) * (psi.3.lg.U > sim.data[[4]][[1]][3])
        )
      
      # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean true species probabilities and classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
      
      betas <- rbind(betas, 
                     cbind(tmp.betas.psi, 
                           as.matrix(tmp.psi), 
                           psi.1.t = sim.data[[4]][[1]][[1]], 
                           psi.2.t = sim.data[[4]][[1]][[2]], 
                           psi.3.t = sim.data[[4]][[1]][[3]], 
                           n.sim = n.sim))
      
      # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
      
      output.betas <- 
        ftransform(tmp.psi, 
               psi.1.err = psi.1 - sim.data[[4]][[1]][[1]],
               psi.2.err = psi.2 - sim.data[[4]][[1]][[2]],
               psi.3.err = psi.3 - sim.data[[4]][[1]][[3]]
        ) %>%
        summarise(., 
                  me.psi.1 = mean(psi.1.err), 
                  me.psi.2 = mean(psi.2.err),
                  me.psi.3 = mean(psi.3.err),
                  rm.psi.1 = (sum(psi.1.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                  rm.psi.2 = (sum(psi.2.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                  rm.psi.3 = (sum(psi.3.err ^ 2) / (dim(tmp.psi)[1] - 1)) ^ 0.5,
                  ci.lg.psi.1 = mean(ci.lg.psi.1, na.rm = TRUE),
                  ci.lg.psi.2 = mean(ci.lg.psi.2, na.rm = TRUE),
                  ci.lg.psi.3 = mean(ci.lg.psi.3, na.rm = TRUE)
        ) %>%
        bind_rows(output.betas, .)
      
      # Return list with output
      return(list(output.betas, betas))
    } 
  if (B > 3 | B < 2) {
    cat(paste0("Derived parameters for ", B, " true species states not supported \n"))
  }else {
    cat("Derived parameter estimates not produced \n")
  }
}

# Function: {theta.derived.f} Produces overall mean estimates of classification probabilities (theta) derived from regression coefficients (betas) from multinomial logistic regression models along with associated standard errors, 95% confidence limits, and 95% confidence interval coverage. Supports models with regression coefficients predicting classification probabilities with B = 2 or 3 true species states and A = 2 or 3 observation states.

# Accepts inputs of:
# 'sim.data' and 'sim.data.tmp' containing simulated survey data 
# 'model.results' containing simulation output from model optimization 
# 'parameter.output' containing summarized statistics for estimated parameters 
# 'output.betas' containing summarized statistics for overall means for classification probabilities from prior distinct models 

# Returns list of output with:
# 'output.betas' containing  summarized statistics (mean, root mean square error, 95% CI coverage) across simulation replicates for overall mean classification probabilities. Contains results for the current and all prior distinct models.
# 'betas' containing (for each simulation replicate) estimated regression coefficients and estimated means, SEs, confidence limits, and 95% confidence interval coverage for overall mean true species probabilities and classification probabilities

theta.derived.f <- function(sim.data, sim.data.tmp, model.results, parameter.output, output.betas) {

  # Optional output for derived estimates of overall mean true species probabilities (psi) and classification probabilities (theta)-
  
  # 'tmp.id' and 'n.sim' contain an index and the count of simulation replicates successfully optimizing
  # 'tmp.betas' contains estimated intercept and slope parameters for regression coefficients predicting classification probabilities
  
  tmp.id <- unique(sim.data.tmp$id)
  n.sim <- rep(sim, reps)
  tmp.betas <- cbind(parameter.output$par[[1]][, c(grep("b0|b1", parameters_names))])
  colnames(tmp.betas) <- c(grep("b0|b1", parameters_names, value = TRUE))

  # Estimate overall mean classification probabilities (psi) and associated SEs (using Delta Method)
  tmp.theta <- 
    t(sapply(1:reps, function(x)
      theta.est.f(parameter.output$par[[1]][x, ],  sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ])
    ))
  
  tmp.se <- 
    t(sapply(1:reps, function(x)
      diag(
        DeltaMethod(
          parameter.output$par[[1]][x, ],
          theta.est.f,
          vcov = solve(model.results$hessian[[x]]),
          delta = 0.00001,
          dat = sim.data.tmp[which(sim.data.tmp$id == tmp.id[x]), ]
        )$variance
      )) ^ 0.5
    )

  # Apply logit transformation to estimates and standard errors, add 95% confidence limits assuming asymptotic normal distribution (in logit scale), and apply inverse logit transformation to CLs. Confidence limits from logit transform have '.lg.L' and '.lg.U' suffixes.  
  tmp.lg <- 
    t(sapply(1:reps, function(x)
      logit.ci.f(unlist(tmp.theta[x, ]), unlist(tmp.se[x, ]))))
  
  if (sim_profiles[1, ]$Model != "M.theta.p") {
    if (B == 2 & A == 2) {
      ## ----- True species states (B) = observation states (A) = 2 -----
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and standard errors ('se.' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.12", "se.theta.s.21", "se.theta.s.12")[c(O_ps[1]*1:2, 3:4)]
      
      colnames(tmp.lg) <- c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", 
                            "theta.s.21.lg.L", "theta.s.21.lg.U", "theta.s.12.lg.L", "theta.s.12.lg.U")[c(O_ps[1]*1:4, 5:8)]

      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     
                     theta.p.21.L = theta.p.21 - 1.96 * se.theta.p.21,
                     theta.p.21.U = theta.p.21 + 1.96 * se.theta.p.21,
                     theta.p.12.L = theta.p.12 - 1.96 * se.theta.p.12,
                     theta.p.12.U = theta.p.12 + 1.96 * se.theta.p.12, 
                     
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12) %>%
          
          ftransform(., 
                     
                     ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1]) * (theta.p.21.U > sim.data[[4]][[2]][1]),
                     ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][2]) * (theta.p.12.U > sim.data[[4]][[2]][2]),
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][3]) * (theta.s.21.U > sim.data[[4]][[2]][3]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][4]) * (theta.s.12.U > sim.data[[4]][[2]][4]),
                     
                     ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1]),
                     ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][2]),
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][3]) * (theta.s.21.lg.U > sim.data[[4]][[2]][3]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][4]) * (theta.s.12.lg.U > sim.data[[4]][[2]][4])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1], 
                      theta.p.12.t = sim.data[[4]][[2]][2], 
                      theta.s.21.t = sim.data[[4]][[2]][3], 
                      theta.s.12.t = sim.data[[4]][[2]][4],
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta.p.21 - sim.data[[4]][[2]][[1]],
                 theta.p.12.err = theta.p.12 - sim.data[[4]][[2]][[2]],
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][[3]],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][[4]]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err),
                    me.theta.p.12 = mean(theta.p.12.err),
                    me.theta.s.21 = mean(theta.s.21.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE),
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE),  
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE)
                    
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
      } else {
        # ----- Single-observation method (SOM) models (primary observers NOT present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval covarage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12) %>%
          
          ftransform(., 
                     
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1]) * (theta.s.21.U > sim.data[[4]][[2]][1]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][2]) * (theta.s.12.U > sim.data[[4]][[2]][2]),
                     
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][2]) * (theta.s.12.lg.U > sim.data[[4]][[2]][2])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.s.21.t = sim.data[[4]][[2]][1], 
                      theta.s.12.t = sim.data[[4]][[2]][2],
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][[1]],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][[2]]
          ) %>%
          summarise(., 
                    me.theta.s.21 = mean(theta.s.21.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE)
                    
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
      } 
    }else if (B == 2 & A == 3) {
      # ----- True species states = 2, observation states = 3 -----
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and associated standard errors ('se.theta.' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.31", "se.theta.p.12", "se.theta.p.32", 
                            "se.theta.s.21", "se.theta.s.31", "se.theta.s.12", "se.theta.s.32")[c(O_ps[1]*1:4, 5:8)]
      
      colnames(tmp.lg) <- 
        c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.31.lg.L", "theta.p.31.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", "theta.p.32.lg.L", "theta.p.32.lg.U", 
          "theta.s.21.lg.L", "theta.s.21.lg.U", "theta.s.31.lg.L", "theta.s.31.lg.U", "theta.s.12.lg.L", "theta.s.12.lg.U", "theta.s.32.lg.L", "theta.s.32.lg.U")[c(O_ps[1]*1:8, 9:16)]
      
      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     theta.p.21.L = theta.p.21 - 1.96 * se.theta.p.21,
                     theta.p.21.U = theta.p.21 + 1.96 * se.theta.p.21,
                     theta.p.31.L = theta.p.31 - 1.96 * se.theta.p.31,
                     theta.p.31.U = theta.p.31 + 1.96 * se.theta.p.31,
                     
                     theta.p.12.L = theta.p.12 - 1.96 * se.theta.p.12,
                     theta.p.12.U = theta.p.12 + 1.96 * se.theta.p.12, 
                     theta.p.32.L = theta.p.32 - 1.96 * se.theta.p.32,
                     theta.p.32.U = theta.p.32 + 1.96 * se.theta.p.32, 
                     
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.31.L = theta.s.31 - 1.96 * se.theta.s.31,
                     theta.s.31.U = theta.s.31 + 1.96 * se.theta.s.31,
                     
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12, 
                     theta.s.32.L = theta.s.32 - 1.96 * se.theta.s.32,
                     theta.s.32.U = theta.s.32 + 1.96 * se.theta.s.32) %>%
          
          ftransform(., 
                     ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.U > sim.data[[4]][[2]][1, 1]),
                     ci.theta.p.31 = (theta.p.31.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.U > sim.data[[4]][[2]][2, 1]),
                     ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.U > sim.data[[4]][[2]][1, 2]),
                     ci.theta.p.32 = (theta.p.32.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.U > sim.data[[4]][[2]][2, 2]),
                     
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1, 3]) * (theta.s.21.U > sim.data[[4]][[2]][1, 3]),
                     ci.theta.s.31 = (theta.s.31.L < sim.data[[4]][[2]][2, 3]) * (theta.s.31.U > sim.data[[4]][[2]][2, 3]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][1, 4]) * (theta.s.12.U > sim.data[[4]][[2]][1, 4]),
                     ci.theta.s.32 = (theta.s.32.L < sim.data[[4]][[2]][2, 4]) * (theta.s.32.U > sim.data[[4]][[2]][2, 4]),
                     
                     ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1, 1]),
                     ci.lg.theta.p.31 = (theta.p.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.lg.U > sim.data[[4]][[2]][2, 1]),
                     ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][1, 2]),
                     ci.lg.theta.p.32 = (theta.p.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.lg.U > sim.data[[4]][[2]][2, 2]),
                     
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1, 3]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1, 3]),
                     ci.lg.theta.s.31 = (theta.s.31.lg.L < sim.data[[4]][[2]][2, 3]) * (theta.s.31.lg.U > sim.data[[4]][[2]][2, 3]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][1, 4]) * (theta.s.12.lg.U > sim.data[[4]][[2]][1, 4]),
                     ci.lg.theta.s.32 = (theta.s.32.lg.L < sim.data[[4]][[2]][2, 4]) * (theta.s.32.lg.U > sim.data[[4]][[2]][2, 4])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1, 1], 
                      theta.p.31.t = sim.data[[4]][[2]][2, 1], 
                      theta.p.12.t = sim.data[[4]][[2]][1, 2], 
                      theta.p.32.t = sim.data[[4]][[2]][2, 2], 
                      theta.s.21.t = sim.data[[4]][[2]][1, 3], 
                      theta.s.31.t = sim.data[[4]][[2]][2, 3],
                      theta.s.12.t = sim.data[[4]][[2]][1, 4],
                      theta.s.32.t = sim.data[[4]][[2]][2, 4],
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta.p.21 - sim.data[[4]][[2]][1, 1], 
                 theta.p.31.err = theta.p.31 - sim.data[[4]][[2]][2, 1], 
                 theta.p.12.err = theta.p.12 - sim.data[[4]][[2]][1, 2], 
                 theta.p.32.err = theta.p.32 - sim.data[[4]][[2]][2, 2], 
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][1, 3], 
                 theta.s.31.err = theta.s.31 - sim.data[[4]][[2]][2, 3],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][1, 4],
                 theta.s.32.err = theta.s.32 - sim.data[[4]][[2]][2, 4]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err), 
                    me.theta.p.31 = mean(theta.p.31.err), 
                    me.theta.p.12 = mean(theta.p.12.err), 
                    me.theta.p.32 = mean(theta.p.32.err), 
                    me.theta.s.21 = mean(theta.s.21.err), 
                    me.theta.s.31 = mean(theta.s.31.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    me.theta.s.32 = mean(theta.s.32.err),
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.31 = (sum(theta.p.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.32 = (sum(theta.p.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.s.31 = (sum(theta.s.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.32 = (sum(theta.s.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE), 
                    ci.lg.theta.p.31 = mean(ci.lg.theta.p.31, na.rm = TRUE), 
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE), 
                    ci.lg.theta.p.32 = mean(ci.lg.theta.p.32, na.rm = TRUE), 
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE), 
                    ci.lg.theta.s.31 = mean(ci.lg.theta.s.31, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE),
                    ci.lg.theta.s.32 = mean(ci.lg.theta.s.32, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
      }else{
        # ----- Single-observation method (SOM) models (primary observers NOT present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.31.L = theta.s.31 - 1.96 * se.theta.s.31,
                     theta.s.31.U = theta.s.31 + 1.96 * se.theta.s.31,
                     
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12, 
                     theta.s.32.L = theta.s.32 - 1.96 * se.theta.s.32,
                     theta.s.32.U = theta.s.32 + 1.96 * se.theta.s.32) %>%
          
          ftransform(., 
                     
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.U > sim.data[[4]][[2]][1, 1]),
                     ci.theta.s.31 = (theta.s.31.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.U > sim.data[[4]][[2]][2, 1]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.U > sim.data[[4]][[2]][1, 2]),
                     ci.theta.s.32 = (theta.s.32.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.U > sim.data[[4]][[2]][2, 2]),
                     
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1, 1]),
                     ci.lg.theta.s.31 = (theta.s.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.lg.U > sim.data[[4]][[2]][2, 1]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.lg.U > sim.data[[4]][[2]][1, 2]),
                     ci.lg.theta.s.32 = (theta.s.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.lg.U > sim.data[[4]][[2]][2, 2])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.s.21.t = sim.data[[4]][[2]][1, 1], 
                      theta.s.31.t = sim.data[[4]][[2]][2, 1],
                      theta.s.12.t = sim.data[[4]][[2]][1, 2],
                      theta.s.32.t = sim.data[[4]][[2]][2, 2],
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][1, 1], 
                 theta.s.31.err = theta.s.31 - sim.data[[4]][[2]][2, 1],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][1, 2],
                 theta.s.32.err = theta.s.32 - sim.data[[4]][[2]][2, 2]
          ) %>%
          summarise(., 
                    me.theta.s.21 = mean(theta.s.21.err), 
                    me.theta.s.31 = mean(theta.s.31.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    me.theta.s.32 = mean(theta.s.32.err),
                    
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.s.31 = (sum(theta.s.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.32 = (sum(theta.s.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE), 
                    ci.lg.theta.s.31 = mean(ci.lg.theta.s.31, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE),
                    ci.lg.theta.s.32 = mean(ci.lg.theta.s.32, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
      }
    } else if (B == 3 & A == 3) {
      # ----- True species states B = observation states A = 3 -----
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and associated standard errors ('se.theta' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.31", "se.theta.p.12", "se.theta.p.32", "se.theta.p.13", "se.theta.p.23",
                            "se.theta.s.21", "se.theta.s.31", "se.theta.s.12", "se.theta.s.32", "se.theta.s.13", "se.theta.s.23")[c(O_ps[1]*1:6, 7:12)]
      
      colnames(tmp.lg) <- 
        c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.31.lg.L", "theta.p.31.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", "theta.p.32.lg.L", "theta.p.32.lg.U", "theta.p.13.lg.L", "theta.p.13.lg.U", "theta.p.23.lg.L", "theta.p.23.lg.U",
          "theta.s.21.lg.L", "theta.s.21.lg.U", "theta.s.31.lg.L", "theta.s.31.lg.U", "theta.s.12.lg.L", "theta.s.12.lg.U", "theta.s.32.lg.L", "theta.s.32.lg.U", "theta.s.13.lg.L", "theta.s.13.lg.U", "theta.s.23.lg.L", "theta.s.23.lg.U")[c(O_ps[1]*1:12, 13:24)]
      
      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     theta.p.21.L = theta.p.21 - 1.96 * se.theta.p.21,
                     theta.p.21.U = theta.p.21 + 1.96 * se.theta.p.21,
                     theta.p.31.L = theta.p.31 - 1.96 * se.theta.p.31,
                     theta.p.31.U = theta.p.31 + 1.96 * se.theta.p.31,
                     
                     theta.p.12.L = theta.p.12 - 1.96 * se.theta.p.12,
                     theta.p.12.U = theta.p.12 + 1.96 * se.theta.p.12, 
                     theta.p.32.L = theta.p.32 - 1.96 * se.theta.p.32,
                     theta.p.32.U = theta.p.32 + 1.96 * se.theta.p.32,
                     
                     theta.p.13.L = theta.p.13 - 1.96 * se.theta.p.13,
                     theta.p.13.U = theta.p.13 + 1.96 * se.theta.p.13, 
                     theta.p.23.L = theta.p.23 - 1.96 * se.theta.p.23,
                     theta.p.23.U = theta.p.23 + 1.96 * se.theta.p.23,
                     
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.31.L = theta.s.31 - 1.96 * se.theta.s.31,
                     theta.s.31.U = theta.s.31 + 1.96 * se.theta.s.31,
                     
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12, 
                     theta.s.32.L = theta.s.32 - 1.96 * se.theta.s.32,
                     theta.s.32.U = theta.s.32 + 1.96 * se.theta.s.32,
                     
                     theta.s.13.L = theta.s.13 - 1.96 * se.theta.s.13,
                     theta.s.13.U = theta.s.13 + 1.96 * se.theta.s.13, 
                     theta.s.23.L = theta.s.23 - 1.96 * se.theta.s.23,
                     theta.s.23.U = theta.s.23 + 1.96 * se.theta.s.23) %>%
          
          ftransform(., 
                     
                     ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.U > sim.data[[4]][[2]][1, 1]),
                     ci.theta.p.31 = (theta.p.31.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.U > sim.data[[4]][[2]][2, 1]),
                     ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.U > sim.data[[4]][[2]][1, 2]),
                     ci.theta.p.32 = (theta.p.32.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.U > sim.data[[4]][[2]][2, 2]),
                     ci.theta.p.13 = (theta.p.13.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.U > sim.data[[4]][[2]][1, 3]),
                     ci.theta.p.23 = (theta.p.23.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.U > sim.data[[4]][[2]][2, 3]),
                     
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1, 4]) * (theta.s.21.U > sim.data[[4]][[2]][1, 4]),
                     ci.theta.s.31 = (theta.s.31.L < sim.data[[4]][[2]][2, 4]) * (theta.s.31.U > sim.data[[4]][[2]][2, 4]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][1, 5]) * (theta.s.12.U > sim.data[[4]][[2]][1, 5]),
                     ci.theta.s.32 = (theta.s.32.L < sim.data[[4]][[2]][2, 5]) * (theta.s.32.U > sim.data[[4]][[2]][2, 5]),
                     ci.theta.s.13 = (theta.s.13.L < sim.data[[4]][[2]][1, 6]) * (theta.s.13.U > sim.data[[4]][[2]][1, 6]),
                     ci.theta.s.23 = (theta.s.23.L < sim.data[[4]][[2]][2, 6]) * (theta.s.23.U > sim.data[[4]][[2]][2, 6]),
                     
                     ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1, 1]),
                     ci.lg.theta.p.31 = (theta.p.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.lg.U > sim.data[[4]][[2]][2, 1]),
                     ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][1, 2]),
                     ci.lg.theta.p.32 = (theta.p.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.lg.U > sim.data[[4]][[2]][2, 2]),
                     ci.lg.theta.p.13 = (theta.p.13.lg.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.lg.U > sim.data[[4]][[2]][1, 3]),
                     ci.lg.theta.p.23 = (theta.p.23.lg.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.lg.U > sim.data[[4]][[2]][2, 3]),
                     
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1, 4]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1, 4]),
                     ci.lg.theta.s.31 = (theta.s.31.lg.L < sim.data[[4]][[2]][2, 4]) * (theta.s.31.lg.U > sim.data[[4]][[2]][2, 4]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][1, 5]) * (theta.s.12.lg.U > sim.data[[4]][[2]][1, 5]),
                     ci.lg.theta.s.32 = (theta.s.32.lg.L < sim.data[[4]][[2]][2, 5]) * (theta.s.32.lg.U > sim.data[[4]][[2]][2, 5]),
                     ci.lg.theta.s.13 = (theta.s.13.lg.L < sim.data[[4]][[2]][1, 6]) * (theta.s.13.lg.U > sim.data[[4]][[2]][1, 6]),
                     ci.lg.theta.s.23 = (theta.s.23.lg.L < sim.data[[4]][[2]][2, 6]) * (theta.s.23.lg.U > sim.data[[4]][[2]][2, 6])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <-  
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1, 1],
                      theta.p.31.t = sim.data[[4]][[2]][2, 1],
                      theta.p.12.t = sim.data[[4]][[2]][1, 2],
                      theta.p.32.t = sim.data[[4]][[2]][2, 2],
                      theta.p.13.t = sim.data[[4]][[2]][1, 3],
                      theta.p.23.t = sim.data[[4]][[2]][2, 3],
                      
                      theta.s.21.t = sim.data[[4]][[2]][1, 4],
                      theta.s.31.t = sim.data[[4]][[2]][2, 4],
                      theta.s.12.t = sim.data[[4]][[2]][1, 5],
                      theta.s.32.t = sim.data[[4]][[2]][2, 5],
                      theta.s.13.t = sim.data[[4]][[2]][1, 6],
                      theta.s.23.t = sim.data[[4]][[2]][2, 6],
                      
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta.p.21 - sim.data[[4]][[2]][1, 1],
                 theta.p.31.err = theta.p.31 - sim.data[[4]][[2]][2, 1],
                 theta.p.12.err = theta.p.12 - sim.data[[4]][[2]][1, 2],
                 theta.p.32.err = theta.p.32 - sim.data[[4]][[2]][2, 2],
                 theta.p.13.err = theta.p.13 - sim.data[[4]][[2]][1, 3],
                 theta.p.23.err = theta.p.23 - sim.data[[4]][[2]][2, 3],
                 
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][1, 4],
                 theta.s.31.err = theta.s.31 - sim.data[[4]][[2]][2, 4],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][1, 5],
                 theta.s.32.err = theta.s.32 - sim.data[[4]][[2]][2, 5],
                 theta.s.13.err = theta.s.13 - sim.data[[4]][[2]][1, 6],
                 theta.s.23.err = theta.s.23 - sim.data[[4]][[2]][2, 6]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err),
                    me.theta.p.31 = mean(theta.p.31.err),
                    me.theta.p.12 = mean(theta.p.12.err),
                    me.theta.p.32 = mean(theta.p.32.err),
                    me.theta.p.13 = mean(theta.p.13.err),
                    me.theta.p.23 = mean(theta.p.23.err),
                    
                    me.theta.s.21 = mean(theta.s.21.err),
                    me.theta.s.31 = mean(theta.s.31.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    me.theta.s.32 = mean(theta.s.32.err),
                    me.theta.s.13 = mean(theta.s.13.err),
                    me.theta.s.23 = mean(theta.s.23.err),
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.31 = (sum(theta.p.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.32 = (sum(theta.p.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.13 = (sum(theta.p.13.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.23 = (sum(theta.p.23.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.31 = (sum(theta.s.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.32 = (sum(theta.s.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.13 = (sum(theta.s.13.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.23 = (sum(theta.s.23.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE),
                    ci.lg.theta.p.31 = mean(ci.lg.theta.p.31, na.rm = TRUE),
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE),
                    ci.lg.theta.p.32 = mean(ci.lg.theta.p.32, na.rm = TRUE),
                    ci.lg.theta.p.13 = mean(ci.lg.theta.p.13, na.rm = TRUE),
                    ci.lg.theta.p.23 = mean(ci.lg.theta.p.23, na.rm = TRUE),
                    
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE),
                    ci.lg.theta.s.31 = mean(ci.lg.theta.s.31, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE),
                    ci.lg.theta.s.32 = mean(ci.lg.theta.s.32, na.rm = TRUE),
                    ci.lg.theta.s.13 = mean(ci.lg.theta.s.13, na.rm = TRUE),
                    ci.lg.theta.s.23 = mean(ci.lg.theta.s.23, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
      }else {
        # ----- Single-observation method (SOM) models (primary observers NOT present) -----
        
        # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval covarage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
        
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          ftransform(., 
                     theta.s.21.L = theta.s.21 - 1.96 * se.theta.s.21,
                     theta.s.21.U = theta.s.21 + 1.96 * se.theta.s.21,
                     theta.s.31.L = theta.s.31 - 1.96 * se.theta.s.31,
                     theta.s.31.U = theta.s.31 + 1.96 * se.theta.s.31,
                     
                     theta.s.12.L = theta.s.12 - 1.96 * se.theta.s.12,
                     theta.s.12.U = theta.s.12 + 1.96 * se.theta.s.12, 
                     theta.s.32.L = theta.s.32 - 1.96 * se.theta.s.32,
                     theta.s.32.U = theta.s.32 + 1.96 * se.theta.s.32,
                     
                     theta.s.13.L = theta.s.13 - 1.96 * se.theta.s.13,
                     theta.s.13.U = theta.s.13 + 1.96 * se.theta.s.13, 
                     theta.s.23.L = theta.s.23 - 1.96 * se.theta.s.23,
                     theta.s.23.U = theta.s.23 + 1.96 * se.theta.s.23) %>%
          
          ftransform(., 
                     
                     ci.theta.s.21 = (theta.s.21.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.U > sim.data[[4]][[2]][1, 1]),
                     ci.theta.s.31 = (theta.s.31.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.U > sim.data[[4]][[2]][2, 1]),
                     ci.theta.s.12 = (theta.s.12.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.U > sim.data[[4]][[2]][1, 2]),
                     ci.theta.s.32 = (theta.s.32.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.U > sim.data[[4]][[2]][2, 3]),
                     ci.theta.s.13 = (theta.s.13.L < sim.data[[4]][[2]][1, 3]) * (theta.s.13.U > sim.data[[4]][[2]][1, 3]),
                     ci.theta.s.23 = (theta.s.23.L < sim.data[[4]][[2]][2, 3]) * (theta.s.23.U > sim.data[[4]][[2]][2, 3]),
                     
                     ci.lg.theta.s.21 = (theta.s.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.lg.U > sim.data[[4]][[2]][1, 1]),
                     ci.lg.theta.s.31 = (theta.s.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.lg.U > sim.data[[4]][[2]][2, 1]),
                     ci.lg.theta.s.12 = (theta.s.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.lg.U > sim.data[[4]][[2]][1, 2]),
                     ci.lg.theta.s.32 = (theta.s.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.lg.U > sim.data[[4]][[2]][2, 2]),
                     ci.lg.theta.s.13 = (theta.s.13.lg.L < sim.data[[4]][[2]][1, 3]) * (theta.s.13.lg.U > sim.data[[4]][[2]][1, 3]),
                     ci.lg.theta.s.23 = (theta.s.23.lg.L < sim.data[[4]][[2]][2, 3]) * (theta.s.23.lg.U > sim.data[[4]][[2]][2, 3])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <-  
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.s.21.t = sim.data[[4]][[2]][1, 1],
                      theta.s.31.t = sim.data[[4]][[2]][2, 1],
                      theta.s.12.t = sim.data[[4]][[2]][1, 2],
                      theta.s.32.t = sim.data[[4]][[2]][2, 2],
                      theta.s.13.t = sim.data[[4]][[2]][1, 3],
                      theta.s.23.t = sim.data[[4]][[2]][2, 3],
                      
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.s.21.err = theta.s.21 - sim.data[[4]][[2]][1, 1],
                 theta.s.31.err = theta.s.31 - sim.data[[4]][[2]][2, 1],
                 theta.s.12.err = theta.s.12 - sim.data[[4]][[2]][1, 2],
                 theta.s.32.err = theta.s.32 - sim.data[[4]][[2]][2, 2],
                 theta.s.13.err = theta.s.13 - sim.data[[4]][[2]][1, 3],
                 theta.s.23.err = theta.s.23 - sim.data[[4]][[2]][2, 3]
          ) %>%
          summarise(., 
                    
                    me.theta.s.21 = mean(theta.s.21.err),
                    me.theta.s.31 = mean(theta.s.31.err),
                    me.theta.s.12 = mean(theta.s.12.err),
                    me.theta.s.32 = mean(theta.s.32.err),
                    me.theta.s.13 = mean(theta.s.13.err),
                    me.theta.s.23 = mean(theta.s.23.err),
                    
                    rm.theta.s.21 = (sum(theta.s.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.31 = (sum(theta.s.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.12 = (sum(theta.s.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.32 = (sum(theta.s.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.13 = (sum(theta.s.13.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.s.23 = (sum(theta.s.23.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.s.21 = mean(ci.lg.theta.s.21, na.rm = TRUE),
                    ci.lg.theta.s.31 = mean(ci.lg.theta.s.31, na.rm = TRUE),
                    ci.lg.theta.s.12 = mean(ci.lg.theta.s.12, na.rm = TRUE),
                    ci.lg.theta.s.32 = mean(ci.lg.theta.s.32, na.rm = TRUE),
                    ci.lg.theta.s.13 = mean(ci.lg.theta.s.13, na.rm = TRUE),
                    ci.lg.theta.s.23 = mean(ci.lg.theta.s.23, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
      }
    }
  }else {
    ## ----- For models with covariate predicting classification probabilities ONLY for primary observers -----
    # Specified as model = 'M.theta.p'
    
    if (B == 2 & A == 2) {
      ## True species states (B) = observation states (A) = 2
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and associated standard errors ('se.theta.' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.12")
      colnames(tmp.lg) <- c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U")
      
      # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
      
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          mutate(., 
                 theta.p.21.L = theta_p_21 - 1.96 * se.theta.p.21,
                 theta.p.21.U = theta_p_21 + 1.96 * se.theta.p.21,
                 theta.p.12.L = theta_p_12 - 1.96 * se.theta.p.12,
                 theta.p.12.U = theta_p_12 + 1.96 * se.theta.p.12, 
                 
                 ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1]) * (theta.p.21.U > sim.data[[4]][[2]][1]),
                 ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][2]) * (theta.p.12.U > sim.data[[4]][[2]][2]),
                 
                 ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1]),
                 ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][2])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1], 
                      theta.p.12.t = sim.data[[4]][[2]][2], 
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta_p_21 - sim.data[[4]][[2]][[1]],
                 theta.p.12.err = theta_p_12 - sim.data[[4]][[2]][[2]]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err),
                    me.theta.p.12 = mean(theta.p.12.err),
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE),
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE)
                    
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
    }else if (B == 2 & A == 3) {
      ## True species states (B) = 2, observation states (A) = 3
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and associated standard errors ('se.theta.' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.31", "se.theta.p.12", "se.theta.p.32")
      colnames(tmp.lg) <- 
        c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.31.lg.L", "theta.p.31.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", "theta.p.32.lg.L", "theta.p.32.lg.U")
      
      # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval covarage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
      
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          mutate(., 
                 theta.p.21.L = theta_p_21 - 1.96 * se.theta.p.21,
                 theta.p.21.U = theta_p_21 + 1.96 * se.theta.p.21,
                 theta.p.31.L = theta_p_31 - 1.96 * se.theta.p.31,
                 theta.p.31.U = theta_p_31 + 1.96 * se.theta.p.31,
                 
                 theta.p.12.L = theta_p_12 - 1.96 * se.theta.p.12,
                 theta.p.12.U = theta_p_12 + 1.96 * se.theta.p.12, 
                 theta.p.32.L = theta_p_32 - 1.96 * se.theta.p.32,
                 theta.p.32.U = theta_p_32 + 1.96 * se.theta.p.32,  
                 
                 ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.U > sim.data[[4]][[2]][1, 1]),
                 ci.theta.p.31 = (theta.p.31.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.U > sim.data[[4]][[2]][2, 1]),
                 ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.U > sim.data[[4]][[2]][1, 2]),
                 ci.theta.p.32 = (theta.p.32.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.U > sim.data[[4]][[2]][2, 2]),
                 
                 ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1, 1]),
                 ci.lg.theta.p.31 = (theta.p.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.lg.U > sim.data[[4]][[2]][2, 1]),
                 ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][1, 2]),
                 ci.lg.theta.p.32 = (theta.p.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.lg.U > sim.data[[4]][[2]][2, 2])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <- 
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1, 1], 
                      theta.p.31.t = sim.data[[4]][[2]][2, 1], 
                      theta.p.12.t = sim.data[[4]][[2]][1, 2], 
                      theta.p.32.t = sim.data[[4]][[2]][2, 2], 
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta_p_21 - sim.data[[4]][[2]][1, 1], 
                 theta.p.31.err = theta_p_31 - sim.data[[4]][[2]][2, 1], 
                 theta.p.12.err = theta_p_12 - sim.data[[4]][[2]][1, 2], 
                 theta.p.32.err = theta_p_32 - sim.data[[4]][[2]][2, 2]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err), 
                    me.theta.p.31 = mean(theta.p.31.err), 
                    me.theta.p.12 = mean(theta.p.12.err), 
                    me.theta.p.32 = mean(theta.p.32.err), 
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.31 = (sum(theta.p.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    rm.theta.p.32 = (sum(theta.p.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5, 
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE), 
                    ci.lg.theta.p.31 = mean(ci.lg.theta.p.31, na.rm = TRUE), 
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE), 
                    ci.lg.theta.p.32 = mean(ci.lg.theta.p.32, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
        
    } else if (B == 3 & A == 3) {
      ## True species states (B) = observation states (A) = 3
      
      # Add parameter names for estimated classification probabilities ('theta.' prefix) and associated standard errors ('se.theta.' prefix)
      colnames(tmp.se) <- c("se.theta.p.21", "se.theta.p.31", "se.theta.p.12", "se.theta.p.32", "se.theta.p.13", "se.theta.p.23")
      colnames(tmp.lg) <- 
        c("theta.p.21.lg.L", "theta.p.21.lg.U", "theta.p.31.lg.L", "theta.p.31.lg.U", "theta.p.12.lg.L", "theta.p.12.lg.U", "theta.p.32.lg.L", "theta.p.32.lg.U", "theta.p.13.lg.L", "theta.p.13.lg.U", "theta.p.23.lg.L", "theta.p.23.lg.U")
      
      # Add 95% confidence limits ('.L' and '.U' suffixes) constructed assuming asymptotic normal distribution (in natural scale) and associated 95% coverage interval coverage for natural ('ci.' prefix) and logit scale ('ci.lg.' prefix)
      
        tmp.theta <- 
          qTBL(cbind(tmp.theta, tmp.se, tmp.lg)) %>%
          mutate(., 
                 theta.p.21.L = theta_p_21 - 1.96 * se.theta.p.21,
                 theta.p.21.U = theta_p_21 + 1.96 * se.theta.p.21,
                 theta.p.31.L = theta_p_31 - 1.96 * se.theta.p.31,
                 theta.p.31.U = theta_p_31 + 1.96 * se.theta.p.31,
                 
                 theta.p.12.L = theta_p_12 - 1.96 * se.theta.p.12,
                 theta.p.12.U = theta_p_12 + 1.96 * se.theta.p.12, 
                 theta.p.32.L = theta_p_32 - 1.96 * se.theta.p.32,
                 theta.p.32.U = theta_p_32 + 1.96 * se.theta.p.32,
                 
                 theta.p.13.L = theta_p_13 - 1.96 * se.theta.p.13,
                 theta.p.13.U = theta_p_13 + 1.96 * se.theta.p.13, 
                 theta.p.23.L = theta_p_23 - 1.96 * se.theta.p.23,
                 theta.p.23.U = theta_p_23 + 1.96 * se.theta.p.23,
                 
                 ci.theta.p.21 = (theta.p.21.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.U > sim.data[[4]][[2]][1, 1]),
                 ci.theta.p.31 = (theta.p.31.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.U > sim.data[[4]][[2]][2, 1]),
                 ci.theta.p.12 = (theta.p.12.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.U > sim.data[[4]][[2]][1, 2]),
                 ci.theta.p.32 = (theta.p.32.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.U > sim.data[[4]][[2]][2, 2]),
                 ci.theta.p.13 = (theta.p.13.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.U > sim.data[[4]][[2]][1, 3]),
                 ci.theta.p.23 = (theta.p.23.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.U > sim.data[[4]][[2]][2, 3]),
                 
                 ci.lg.theta.p.21 = (theta.p.21.lg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.lg.U > sim.data[[4]][[2]][1, 1]),
                 ci.lg.theta.p.31 = (theta.p.31.lg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.lg.U > sim.data[[4]][[2]][2, 1]),
                 ci.lg.theta.p.12 = (theta.p.12.lg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.lg.U > sim.data[[4]][[2]][1, 2]),
                 ci.lg.theta.p.32 = (theta.p.32.lg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.lg.U > sim.data[[4]][[2]][2, 2]),
                 ci.lg.theta.p.13 = (theta.p.13.lg.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.lg.U > sim.data[[4]][[2]][1, 3]),
                 ci.lg.theta.p.23 = (theta.p.23.lg.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.lg.U > sim.data[[4]][[2]][2, 3])
          )
        
        # Output matrix 'betas' contains results for individual simulation replicates: estimated regression coefficients, estimated overall mean classification probabilities, 95% confidence limits and confidence interval coverage, true values, and an index the distinct model (simulation profile)
        
        betas <-  
          rbind(betas, 
                cbind(tmp.betas, 
                      as.matrix(tmp.theta), 
                      theta.p.21.t = sim.data[[4]][[2]][1, 1],
                      theta.p.31.t = sim.data[[4]][[2]][2, 1],
                      theta.p.12.t = sim.data[[4]][[2]][1, 2],
                      theta.p.32.t = sim.data[[4]][[2]][2, 2],
                      theta.p.13.t = sim.data[[4]][[2]][1, 3],
                      theta.p.23.t = sim.data[[4]][[2]][2, 3],
                      
                      n.sim = n.sim)
          )
        
        # Output data frame 'output.betas' contains results for each simulation summarized across replicates: mean error ('me.' prefix), root mean square error ('rm.' prefix), and 95% confidence interval coverage (with confidence limits constructed in logit scale) for estimated overall means ('ci.lg' prefix)
        
        output.betas <- 
          ftransform(tmp.theta, 
                 theta.p.21.err = theta_p_21 - sim.data[[4]][[2]][1, 1],
                 theta.p.31.err = theta_p_31 - sim.data[[4]][[2]][2, 1],
                 theta.p.12.err = theta_p_12 - sim.data[[4]][[2]][1, 2],
                 theta.p.32.err = theta_p_32 - sim.data[[4]][[2]][2, 2],
                 theta.p.13.err = theta_p_13 - sim.data[[4]][[2]][1, 3],
                 theta.p.23.err = theta_p_23 - sim.data[[4]][[2]][2, 3]
          ) %>%
          summarise(., 
                    me.theta.p.21 = mean(theta.p.21.err),
                    me.theta.p.31 = mean(theta.p.31.err),
                    me.theta.p.12 = mean(theta.p.12.err),
                    me.theta.p.32 = mean(theta.p.32.err),
                    me.theta.p.13 = mean(theta.p.13.err),
                    me.theta.p.23 = mean(theta.p.23.err),
                    
                    rm.theta.p.21 = (sum(theta.p.21.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.31 = (sum(theta.p.31.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.12 = (sum(theta.p.12.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.32 = (sum(theta.p.32.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.13 = (sum(theta.p.13.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    rm.theta.p.23 = (sum(theta.p.23.err ^ 2) / (dim(tmp.theta)[1] - 1)) ^ 0.5,
                    
                    ci.lg.theta.p.21 = mean(ci.lg.theta.p.21, na.rm = TRUE),
                    ci.lg.theta.p.31 = mean(ci.lg.theta.p.31, na.rm = TRUE),
                    ci.lg.theta.p.12 = mean(ci.lg.theta.p.12, na.rm = TRUE),
                    ci.lg.theta.p.32 = mean(ci.lg.theta.p.32, na.rm = TRUE),
                    ci.lg.theta.p.13 = mean(ci.lg.theta.p.13, na.rm = TRUE),
                    ci.lg.theta.p.23 = mean(ci.lg.theta.p.23, na.rm = TRUE)
          ) %>%
          bind_rows(output.betas, .)
        
        # Return list of output
        return(list(output.betas, betas))
    }
  }
}

# Function: {theta.bootstrap.f} Estimates bootstrap standard errors, 95% confidence limits, and 95% confidence interval coverage for estimates of overall mean classification probabilities (theta) derived from regression coefficients of multinomial logit regression models with up 3 true species states and 3 observation states. Results are summarized for each distinct model. Confidence intervals are estimated using the bootstrap estimate of the standard error for parameters to place confidence limits on logit-transformed estimates assuming the sampling distribution follows a normal distribution. Models with multinomial logit regressions predicting both true species probabilities and classification probabilities are not supported. 

# Accepts inputs of:
# 'sim.data' containing containing simulated survey data
# 'tmp.theta' containing  estimates of mean classification probabilities and associated standard errors and confidence intervals for each simulation replicate (rows) of the current distinct model
# 'tmp.theta.boot' containing estimated mean classification probabilities for each bootstrap re-sample of each simulation replicate (rows) of the current distinct model 
# 'output.boot' with summarized bootstrap estimates of standard errors and confidence interval coverage for estimates of mean classification probabilities for each distinct simulation model (rows)

# Returns data frame 'output.boot', appending summarized bootstrap estimates for the current distinct simulation model

theta.bootstrap.f <- function(sim.data, tmp.theta, tmp.theta.boot, output.boot) {
  
  # All models except "M.theta.p"  
  if (model != "M.theta.p") {
    if (B == 2 & A == 2) {
      ## ----- True species states (B) = observation states (A) = 2 -----
      
      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes).
        
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.p.21 = sd(theta.p.21, na.rm = TRUE),
                    se.bt.theta.p.12 = sd(theta.p.12, na.rm = TRUE),
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE)
          )
        
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:4]), unlist(tmp.boot.se[x, -1]))))
        colnames(tmp.boot.lg) <- c(
          "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.12.btlg.L", "theta.p.12.btlg.U",
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.12.btlg.L", "theta.s.12.btlg.U"
        )
  
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes) 
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1]),
                 ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][2]),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][3]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][3]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][4]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][4])
          ) %>%
          select(., 9:12) %>%
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarize(., 
        #                   se.bt.theta.p.21 = sd(theta.p.21),
        #                   se.bt.theta.p.12 = sd(theta.p.12),
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   theta.p.21.bt.L = quantile(theta.p.21, 0.025),
        #                   theta.p.12.bt.L = quantile(theta.p.12, 0.025),
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.p.21.bt.U = quantile(theta.p.21, 0.975),
        #                   theta.p.12.bt.U = quantile(theta.p.12, 0.975),
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.p.21 = (theta.p.21.bt.L < sim.data[[4]][[2]][1]) * (theta.p.21.bt.U > sim.data[[4]][[2]][1]),
        #                ci.bt.theta.p.12 = (theta.p.12.bt.L < sim.data[[4]][[2]][2]) * (theta.p.12.bt.U > sim.data[[4]][[2]][2]),
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][3]) * (theta.s.21.bt.U > sim.data[[4]][[2]][3]),
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][4]) * (theta.s.12.bt.U > sim.data[[4]][[2]][4])
        #         ) %>%
        #         select(., 2:5, 14:17) %>%
        #         summarize_all(., mean)
        #     )
        #   )
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        # return(output.boot)
      }else {
        
        # ----- Single-observation method (SOM) models (primary observers NOT present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
  
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE)
          )
  
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:2]), unlist(tmp.boot.se[x, -1]))))
        
        colnames(tmp.boot.lg) <- c(
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.12.btlg.L", "theta.s.12.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1]), 
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][2]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][2])  
          ) %>%
          select(., 5:6) %>%  
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Dataframe 'output.boot' summarizes confidence interval coverage for each method.
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarise(., 
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][1]) * (theta.s.21.bt.U > sim.data[[4]][[2]][1]), 
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][2]) * (theta.s.12.bt.U > sim.data[[4]][[2]][2])  
        #         ) %>%
        #         select(., 2:3, 8:9) %>% 
        #         summarise_all(., mean)
        #     )
        #   )
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        # return(output.boot)
      }
    } else if (B == 2 & A == 3) {
      ## ----- True species states (B) = 2, observation states (A) = 3 -----
      
      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes).
        
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.p.21 = sd(theta.p.21, na.rm = TRUE),
                    se.bt.theta.p.31 = sd(theta.p.31, na.rm = TRUE),
                    se.bt.theta.p.12 = sd(theta.p.12, na.rm = TRUE),
                    se.bt.theta.p.32 = sd(theta.p.32, na.rm = TRUE),
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.31 = sd(theta.s.31, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE),
                    se.bt.theta.s.32 = sd(theta.s.32, na.rm = TRUE)
          )
        
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:8]), unlist(tmp.boot.se[x, -1]))))
        colnames(tmp.boot.lg) <- c(
          "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.31.btlg.L", "theta.p.31.btlg.U",
          "theta.p.12.btlg.L", "theta.p.12.btlg.U", "theta.p.32.btlg.L", "theta.p.32.btlg.U",
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.31.btlg.L", "theta.s.31.btlg.U",
          "theta.s.12.btlg.L", "theta.s.12.btlg.U", "theta.s.32.btlg.L", "theta.s.32.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1, 1]),
                 ci.btlg.theta.p.31 = (theta.p.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.btlg.U > sim.data[[4]][[2]][2, 1]),
                 ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][1, 2]),
                 ci.btlg.theta.p.32 = (theta.p.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.btlg.U > sim.data[[4]][[2]][2, 2]),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1, 3]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1, 3]),
                 ci.btlg.theta.s.31 = (theta.s.31.btlg.L < sim.data[[4]][[2]][2, 3]) * (theta.s.31.btlg.U > sim.data[[4]][[2]][2, 3]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][1, 4]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][1, 4]),
                 ci.btlg.theta.s.32 = (theta.s.32.btlg.L < sim.data[[4]][[2]][2, 4]) * (theta.s.32.btlg.U > sim.data[[4]][[2]][2, 4])
          ) %>%
          select(., 17:24) %>%
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarize(., 
        #                   se.bt.theta.p.21 = sd(theta.p.21),
        #                   se.bt.theta.p.31 = sd(theta.p.31),
        #                   se.bt.theta.p.12 = sd(theta.p.12),
        #                   se.bt.theta.p.32 = sd(theta.p.32),
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.31 = sd(theta.s.31),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   se.bt.theta.s.32 = sd(theta.s.32),
        #                   
        #                   theta.p.21.bt.L = quantile(theta.p.21, 0.025),
        #                   theta.p.31.bt.L = quantile(theta.p.31, 0.025),
        #                   theta.p.12.bt.L = quantile(theta.p.12, 0.025),
        #                   theta.p.32.bt.L = quantile(theta.p.32, 0.025),
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.31.bt.L = quantile(theta.s.31, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.s.32.bt.L = quantile(theta.s.32, 0.025),
        #                   
        #                   theta.p.21.bt.U = quantile(theta.p.21, 0.975),
        #                   theta.p.31.bt.U = quantile(theta.p.31, 0.975),
        #                   theta.p.12.bt.U = quantile(theta.p.12, 0.975),
        #                   theta.p.32.bt.U = quantile(theta.p.32, 0.975),
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.31.bt.U = quantile(theta.s.31, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975),
        #                   theta.s.32.bt.U = quantile(theta.s.32, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.p.21 = (theta.p.21.bt.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.bt.U > sim.data[[4]][[2]][1, 1]),
        #                ci.bt.theta.p.31 = (theta.p.31.bt.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.bt.U > sim.data[[4]][[2]][2, 1]),
        #                ci.bt.theta.p.12 = (theta.p.12.bt.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.bt.U > sim.data[[4]][[2]][1, 2]),
        #                ci.bt.theta.p.32 = (theta.p.32.bt.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.bt.U > sim.data[[4]][[2]][2, 2]),
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][1, 3]) * (theta.s.21.bt.U > sim.data[[4]][[2]][1, 3]),
        #                ci.bt.theta.s.31 = (theta.s.31.bt.L < sim.data[[4]][[2]][2, 3]) * (theta.s.31.bt.U > sim.data[[4]][[2]][2, 3]),
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][1, 4]) * (theta.s.12.bt.U > sim.data[[4]][[2]][1, 4]),
        #                ci.bt.theta.s.32 = (theta.s.32.bt.L < sim.data[[4]][[2]][2, 4]) * (theta.s.32.bt.U > sim.data[[4]][[2]][2, 4])
        #         ) %>%
        #         select(., 2:9, 26:33) %>%
        #         summarize_all(., mean)
        #     )
        #   )
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        # return(output.boot)
      }else {
        
        # ----- Single-observation method (SOM) models (primary observers NOT present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
        
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.31 = sd(theta.s.31, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE),
                    se.bt.theta.s.32 = sd(theta.s.32, na.rm = TRUE)
          )
        
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:4]), unlist(tmp.boot.se[x, -1]))))
        colnames(tmp.boot.lg) <- c(
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.31.btlg.L", "theta.s.31.btlg.U",
          "theta.s.12.btlg.L", "theta.s.12.btlg.U", "theta.s.32.btlg.L", "theta.s.32.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1, 1]),
                 ci.btlg.theta.s.31 = (theta.s.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.btlg.U > sim.data[[4]][[2]][2, 1]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][1, 2]),
                 ci.btlg.theta.s.32 = (theta.s.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.btlg.U > sim.data[[4]][[2]][2, 2])
          ) %>%
          select(., 9:12) %>%
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Dataframe 'output.boot' summarizes confidence interval coverage for each method. 
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarise(., 
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.31 = sd(theta.s.31),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   se.bt.theta.s.32 = sd(theta.s.32),
        #                   
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.31.bt.L = quantile(theta.s.31, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.s.32.bt.L = quantile(theta.s.32, 0.025),
        #                   
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.31.bt.U = quantile(theta.s.31, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975),
        #                   theta.s.32.bt.U = quantile(theta.s.32, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.bt.U > sim.data[[4]][[2]][1, 1]),
        #                ci.bt.theta.s.31 = (theta.s.31.bt.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.bt.U > sim.data[[4]][[2]][2, 1]),
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.bt.U > sim.data[[4]][[2]][1, 2]),
        #                ci.bt.theta.s.32 = (theta.s.32.bt.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.bt.U > sim.data[[4]][[2]][2, 2])
        #         ) %>%
        #         select(., 2:5, 14:17) %>%
        #         summarize_all(., mean)
        #     )
        #   )
        # Return 'output.boot' with current simulation model results appended as a new row
        # return(output.boot)
      }
    } else if (B == 3) {
      # ----- True species states = observation states = 3 ----
  
      if (O_ps[1] > 0) {
        # ----- Multi-observation method (MOM) models (primary observers present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes).
        
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.p.21 = sd(theta.p.21),
                    se.bt.theta.p.31 = sd(theta.p.31),
                    se.bt.theta.p.12 = sd(theta.p.12),
                    se.bt.theta.p.32 = sd(theta.p.32),
                    se.bt.theta.p.13 = sd(theta.p.13),
                    se.bt.theta.p.23 = sd(theta.p.23),
                    se.bt.theta.s.21 = sd(theta.s.21),
                    se.bt.theta.s.31 = sd(theta.s.31),
                    se.bt.theta.s.12 = sd(theta.s.12),
                    se.bt.theta.s.32 = sd(theta.s.32),
                    se.bt.theta.s.13 = sd(theta.s.13),
                    se.bt.theta.s.23 = sd(theta.s.23)
          )
        
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:12]), unlist(tmp.boot.se[x, 2:13]) )))
        colnames(tmp.boot.lg) <- c(
          "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.31.btlg.L", "theta.p.31.btlg.U", 
          "theta.p.12.btlg.L", "theta.p.12.btlg.U", "theta.p.32.btlg.L", "theta.p.32.btlg.U",
          "theta.p.13.btlg.L", "theta.p.13.btlg.U", "theta.p.23.btlg.L", "theta.p.23.btlg.U",
          
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.31.btlg.L", "theta.s.31.btlg.U", 
          "theta.s.12.btlg.L", "theta.s.12.btlg.U", "theta.s.32.btlg.L", "theta.s.32.btlg.U",
          "theta.s.13.btlg.L", "theta.s.13.btlg.U", "theta.s.23.btlg.L", "theta.s.23.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1, 1]),
                 ci.btlg.theta.p.31 = (theta.p.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.btlg.U > sim.data[[4]][[2]][2, 1]),
                 ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][1, 2]),
                 ci.btlg.theta.p.32 = (theta.p.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.btlg.U > sim.data[[4]][[2]][2, 2]),
                 ci.btlg.theta.p.13 = (theta.p.13.btlg.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.btlg.U > sim.data[[4]][[2]][1, 3]),
                 ci.btlg.theta.p.23 = (theta.p.23.btlg.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.btlg.U > sim.data[[4]][[2]][2, 3]),
                 
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1, 4]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1, 4]),
                 ci.btlg.theta.s.31 = (theta.s.31.btlg.L < sim.data[[4]][[2]][2, 4]) * (theta.s.31.btlg.U > sim.data[[4]][[2]][2, 4]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][1, 5]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][1, 5]),
                 ci.btlg.theta.s.32 = (theta.s.32.btlg.L < sim.data[[4]][[2]][2, 5]) * (theta.s.32.btlg.U > sim.data[[4]][[2]][2, 5]),
                 ci.btlg.theta.s.13 = (theta.s.13.btlg.L < sim.data[[4]][[2]][1, 6]) * (theta.s.13.btlg.U > sim.data[[4]][[2]][1, 6]),
                 ci.btlg.theta.s.23 = (theta.s.23.btlg.L < sim.data[[4]][[2]][2, 6]) * (theta.s.23.btlg.U > sim.data[[4]][[2]][2, 6])
          ) %>%
          select(., 25:36) %>%
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarise(., 
        #                   se.bt.theta.p.21 = sd(theta.p.21),
        #                   se.bt.theta.p.31 = sd(theta.p.31),
        #                   se.bt.theta.p.12 = sd(theta.p.12),
        #                   se.bt.theta.p.32 = sd(theta.p.32),
        #                   se.bt.theta.p.13 = sd(theta.p.13),
        #                   se.bt.theta.p.23 = sd(theta.p.23),
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.31 = sd(theta.s.31),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   se.bt.theta.s.32 = sd(theta.s.32),
        #                   se.bt.theta.s.13 = sd(theta.s.13),
        #                   se.bt.theta.s.23 = sd(theta.s.23),
        #                   theta.p.21.bt.L = quantile(theta.p.21, 0.025),
        #                   theta.p.31.bt.L = quantile(theta.p.31, 0.025),
        #                   theta.p.12.bt.L = quantile(theta.p.12, 0.025),
        #                   theta.p.32.bt.L = quantile(theta.p.32, 0.025),
        #                   theta.p.13.bt.L = quantile(theta.p.13, 0.025),
        #                   theta.p.23.bt.L = quantile(theta.p.23, 0.025),
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.31.bt.L = quantile(theta.s.31, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.s.32.bt.L = quantile(theta.s.32, 0.025),
        #                   theta.s.13.bt.L = quantile(theta.s.13, 0.025),
        #                   theta.s.23.bt.L = quantile(theta.s.23, 0.025),
        #                   
        #                   theta.p.21.bt.U = quantile(theta.p.21, 0.975),
        #                   theta.p.31.bt.U = quantile(theta.p.31, 0.975),
        #                   theta.p.12.bt.U = quantile(theta.p.12, 0.975),
        #                   theta.p.32.bt.U = quantile(theta.p.32, 0.975),
        #                   theta.p.13.bt.U = quantile(theta.p.13, 0.975),
        #                   theta.p.23.bt.U = quantile(theta.p.23, 0.975),
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.31.bt.U = quantile(theta.s.31, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975),
        #                   theta.s.32.bt.U = quantile(theta.s.32, 0.975),
        #                   theta.s.13.bt.U = quantile(theta.s.13, 0.975),
        #                   theta.s.23.bt.U = quantile(theta.s.23, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.p.21 = (theta.p.21.bt.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.bt.U > sim.data[[4]][[2]][1, 1]),
        #                ci.bt.theta.p.31 = (theta.p.31.bt.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.bt.U > sim.data[[4]][[2]][2, 1]),
        #                ci.bt.theta.p.12 = (theta.p.12.bt.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.bt.U > sim.data[[4]][[2]][1, 2]),
        #                ci.bt.theta.p.32 = (theta.p.32.bt.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.bt.U > sim.data[[4]][[2]][2, 2]),
        #                ci.bt.theta.p.13 = (theta.p.13.bt.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.bt.U > sim.data[[4]][[2]][1, 3]),
        #                ci.bt.theta.p.23 = (theta.p.23.bt.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.bt.U > sim.data[[4]][[2]][2, 3]),
        #                
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][1, 4]) * (theta.s.21.bt.U > sim.data[[4]][[2]][1, 4]),
        #                ci.bt.theta.s.31 = (theta.s.31.bt.L < sim.data[[4]][[2]][2, 4]) * (theta.s.31.bt.U > sim.data[[4]][[2]][2, 4]),
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][1, 5]) * (theta.s.12.bt.U > sim.data[[4]][[2]][1, 5]),
        #                ci.bt.theta.s.32 = (theta.s.32.bt.L < sim.data[[4]][[2]][2, 5]) * (theta.s.32.bt.U > sim.data[[4]][[2]][2, 5]),
        #                ci.bt.theta.s.13 = (theta.s.13.bt.L < sim.data[[4]][[2]][1, 6]) * (theta.s.13.bt.U > sim.data[[4]][[2]][1, 6]),
        #                ci.bt.theta.s.23 = (theta.s.23.bt.L < sim.data[[4]][[2]][2, 6]) * (theta.s.23.bt.U > sim.data[[4]][[2]][2, 6])
        #         ) %>%
        #         select(., 2:13, 38:49) %>%
        #         summarise_all(., mean)
        #     )
        #   )
        # Return 'output.boot' with current simulation model results appended as a new row
        # return(output.boot)
      }else{
        
        # ----- Single-observation method (SOM) models (primary observers NOT present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
        
        tmp.boot.se <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.s.21 = sd(theta.s.21),
                    se.bt.theta.s.31 = sd(theta.s.31),
                    se.bt.theta.s.12 = sd(theta.s.12),
                    se.bt.theta.s.32 = sd(theta.s.32),
                    se.bt.theta.s.13 = sd(theta.s.13),
                    se.bt.theta.s.23 = sd(theta.s.23)
          )
        
        tmp.boot.lg <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.theta[x, 1:6]), unlist(tmp.boot.se[x, -1]))))
        colnames(tmp.boot.lg) <- c(
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.31.btlg.L", "theta.s.31.btlg.U", 
          "theta.s.12.btlg.L", "theta.s.12.btlg.U", "theta.s.32.btlg.L", "theta.s.32.btlg.U",
          "theta.s.13.btlg.L", "theta.s.13.btlg.U", "theta.s.23.btlg.L", "theta.s.23.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov <- 
          mutate(as_tibble(tmp.boot.lg),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1, 1]),
                 ci.btlg.theta.s.31 = (theta.s.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.btlg.U > sim.data[[4]][[2]][2, 1]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][1, 2]),
                 ci.btlg.theta.s.32 = (theta.s.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.btlg.U > sim.data[[4]][[2]][2, 2]),
                 ci.btlg.theta.s.13 = (theta.s.13.btlg.L < sim.data[[4]][[2]][1, 3]) * (theta.s.13.btlg.U > sim.data[[4]][[2]][1, 3]),
                 ci.btlg.theta.s.23 = (theta.s.23.btlg.L < sim.data[[4]][[2]][2, 3]) * (theta.s.23.btlg.U > sim.data[[4]][[2]][2, 3])
          ) %>%
          select(., 13:18) %>%
          summarise_all(., mean)
        
        # Return data frame 'output.boot' with results for current distinct model appended as a new row
        tmp.boot.se <- select(tmp.boot.se, -1) %>%
          summarise_all(., mean)
        output.boot <-
          bind_rows(output.boot,
                    bind_cols(tmp.boot.lg.cov, tmp.boot.se))
        
        # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 
        
        # output.boot <- 
        #   bind_rows(
        #     output.boot, 
        #     bind_cols(
        #       tmp.boot.lg.cov,
        #       group_by(as_tibble(tmp.theta.boot), rep) %>%
        #         summarise(., 
        #                   se.bt.theta.s.21 = sd(theta.s.21),
        #                   se.bt.theta.s.31 = sd(theta.s.31),
        #                   se.bt.theta.s.12 = sd(theta.s.12),
        #                   se.bt.theta.s.32 = sd(theta.s.32),
        #                   se.bt.theta.s.13 = sd(theta.s.13),
        #                   se.bt.theta.s.23 = sd(theta.s.23),
        #                   theta.s.21.bt.L = quantile(theta.s.21, 0.025),
        #                   theta.s.31.bt.L = quantile(theta.s.31, 0.025),
        #                   theta.s.12.bt.L = quantile(theta.s.12, 0.025),
        #                   theta.s.32.bt.L = quantile(theta.s.32, 0.025),
        #                   theta.s.13.bt.L = quantile(theta.s.13, 0.025),
        #                   theta.s.23.bt.L = quantile(theta.s.23, 0.025),
        #                   
        #                   theta.s.21.bt.U = quantile(theta.s.21, 0.975),
        #                   theta.s.31.bt.U = quantile(theta.s.31, 0.975),
        #                   theta.s.12.bt.U = quantile(theta.s.12, 0.975),
        #                   theta.s.32.bt.U = quantile(theta.s.32, 0.975),
        #                   theta.s.13.bt.U = quantile(theta.s.13, 0.975),
        #                   theta.s.23.bt.U = quantile(theta.s.23, 0.975)
        #         ) %>%
        #         mutate(.,
        #                ci.bt.theta.s.21 = (theta.s.21.bt.L < sim.data[[4]][[2]][1, 1]) * (theta.s.21.bt.U > sim.data[[4]][[2]][1, 1]),
        #                ci.bt.theta.s.31 = (theta.s.31.bt.L < sim.data[[4]][[2]][2, 1]) * (theta.s.31.bt.U > sim.data[[4]][[2]][2, 1]),
        #                ci.bt.theta.s.12 = (theta.s.12.bt.L < sim.data[[4]][[2]][1, 2]) * (theta.s.12.bt.U > sim.data[[4]][[2]][1, 2]),
        #                ci.bt.theta.s.32 = (theta.s.32.bt.L < sim.data[[4]][[2]][2, 2]) * (theta.s.32.bt.U > sim.data[[4]][[2]][2, 2]),
        #                ci.bt.theta.s.13 = (theta.s.13.bt.L < sim.data[[4]][[2]][1, 3]) * (theta.s.13.bt.U > sim.data[[4]][[2]][1, 3]),
        #                ci.bt.theta.s.23 = (theta.s.23.bt.L < sim.data[[4]][[2]][2, 3]) * (theta.s.23.bt.U > sim.data[[4]][[2]][2, 3])
        #         ) %>%
        #         select(., 2:7, 20:25) %>%
        #         summarise_all(., mean)
        #     )
        #   )
        # Return 'output.boot' with current simulation model results appended as a new row
        # return(output.boot)
      }
    }
  return(output.boot)  
  } else if (model == "M.theta.p") {
    # Only model "M.theta.p"
    
    if (B == 2 & A == 2) {
      ## ----- True species states (B) = observation states (A) = 2 -----
      
      # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
      
      tmp.boot.se <- 
        group_by(as_tibble(tmp.theta.boot), rep) %>%
        summarise(., 
                  se.bt.theta.p.21 = sd(theta_p_21, na.rm = TRUE),
                  se.bt.theta.p.12 = sd(theta_p_12, na.rm = TRUE)
        )
      
      tmp.boot.lg <-
        t(sapply(1:reps, function(x)
          logit.ci.f(unlist(tmp.theta[x, 1:2]), unlist(tmp.boot.se[x, -1]))))
      
      colnames(tmp.boot.lg) <- c(
        "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.12.btlg.L", "theta.p.12.btlg.U"
      )
      
      # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
      
      tmp.boot.lg.cov <- 
        mutate(as_tibble(tmp.boot.lg),
               ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1]), 
               ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][2])  
        ) %>%
        select(., 5:6) %>%  
        summarise_all(., mean)
      
      # Return data frame 'output.boot' with results for current distinct model appended as a new row
      tmp.boot.se <- select(tmp.boot.se, -1) %>%
        summarise_all(., mean)
      output.boot <-
        bind_rows(output.boot,
                  bind_cols(tmp.boot.lg.cov, tmp.boot.se))
    }else if (B == 2 & A == 3) {
      ## ----- True species states (B) = 2, observation states (A) = 3 -----
      
      # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
      
      tmp.boot.se <- 
        group_by(as_tibble(tmp.theta.boot), rep) %>%
        summarise(., 
                  se.bt.theta.p.21 = sd(theta_p_21, na.rm = TRUE),
                  se.bt.theta.p.31 = sd(theta_p_31, na.rm = TRUE),
                  se.bt.theta.p.12 = sd(theta_p_12, na.rm = TRUE),
                  se.bt.theta.p.32 = sd(theta_p_32, na.rm = TRUE)
        )
      
      tmp.boot.lg <-
        t(sapply(1:reps, function(x)
          logit.ci.f(unlist(tmp.theta[x, 1:4]), unlist(tmp.boot.se[x, -1]))))
      colnames(tmp.boot.lg) <- c(
        "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.31.btlg.L", "theta.p.31.btlg.U",
        "theta.p.12.btlg.L", "theta.p.12.btlg.U", "theta.p.32.btlg.L", "theta.p.32.btlg.U"
      )
      
      # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
      
      tmp.boot.lg.cov <- 
        mutate(as_tibble(tmp.boot.lg),
               ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1, 1]),
               ci.btlg.theta.p.31 = (theta.p.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.btlg.U > sim.data[[4]][[2]][2, 1]),
               ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][1, 2]),
               ci.btlg.theta.p.32 = (theta.p.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.btlg.U > sim.data[[4]][[2]][2, 2])
        ) %>%
        select(., 9:12) %>%
        summarise_all(., mean)
      
      # Return data frame 'output.boot' with results for current distinct model appended as a new row
      tmp.boot.se <- select(tmp.boot.se, -1) %>%
        summarise_all(., mean)
      output.boot <-
        bind_rows(output.boot,
                  bind_cols(tmp.boot.lg.cov, tmp.boot.se))
    }else if (B == 3 & A == 3) {
      # ----- True species states = observation states = 3 ----
      
      # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
      
      tmp.boot.se <- 
        group_by(as_tibble(tmp.theta.boot), rep) %>%
        summarise(., 
                  se.bt.theta.p.21 = sd(theta_p_21),
                  se.bt.theta.p.31 = sd(theta_p_31),
                  se.bt.theta.p.12 = sd(theta_p_12),
                  se.bt.theta.p.32 = sd(theta_p_32),
                  se.bt.theta.p.13 = sd(theta_p_13),
                  se.bt.theta.p.23 = sd(theta_p_23)
        )
      
      tmp.boot.lg <-
        t(sapply(1:reps, function(x)
          logit.ci.f(unlist(tmp.theta[x, 1:6]), unlist(tmp.boot.se[x, -1]))))
      colnames(tmp.boot.lg) <- c(
        "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.31.btlg.L", "theta.p.31.btlg.U", 
        "theta.p.12.btlg.L", "theta.p.12.btlg.U", "theta.p.32.btlg.L", "theta.p.32.btlg.U",
        "theta.p.13.btlg.L", "theta.p.13.btlg.U", "theta.p.23.btlg.L", "theta.p.23.btlg.U"
      )
      
      # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
      
      tmp.boot.lg.cov <- 
        mutate(as_tibble(tmp.boot.lg),
               ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1, 1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1, 1]),
               ci.btlg.theta.p.31 = (theta.p.31.btlg.L < sim.data[[4]][[2]][2, 1]) * (theta.p.31.btlg.U > sim.data[[4]][[2]][2, 1]),
               ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][1, 2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][1, 2]),
               ci.btlg.theta.p.32 = (theta.p.32.btlg.L < sim.data[[4]][[2]][2, 2]) * (theta.p.32.btlg.U > sim.data[[4]][[2]][2, 2]),
               ci.btlg.theta.p.13 = (theta.p.13.btlg.L < sim.data[[4]][[2]][1, 3]) * (theta.p.13.btlg.U > sim.data[[4]][[2]][1, 3]),
               ci.btlg.theta.p.23 = (theta.p.23.btlg.L < sim.data[[4]][[2]][2, 3]) * (theta.p.23.btlg.U > sim.data[[4]][[2]][2, 3])
        ) %>%
        select(., 13:18) %>%
        summarise_all(., mean)
      
      # Return data frame 'output.boot' with results for current distinct model appended as a new row
      tmp.boot.se <- select(tmp.boot.se, -1) %>%
        summarise_all(., mean)
      output.boot <-
        bind_rows(output.boot,
                  bind_cols(tmp.boot.lg.cov, tmp.boot.se))
      
    }
      
    return(output.boot)
  }
  }

# Function: {psi.bootstrap.f} Estimates bootstrap standard errors, 95% confidence limits, and 95% confidence interval coverage for estimates of overall mean true species probabilities (psi) derived from regression coefficients of multinomial logit regression models with up 3 true species states. Results are summarized for each distinct model. Confidence intervals are estimated using the bootstrap estimate of the standard error for parameters to place confidence limits on logit-transformed estimates assuming the sampling distribution follows a normal distribution. Models with multinomial logit regressions predicting both true species probabilities and classification probabilities are not supported. 

# Accepts inputs of:
# 'sim.data' and 'sim.data.tmp' containing simulated survey data 
# 'tmp.psi' containing  estimates of mean true species probabilities and associated standard errors and confidence intervals for each simulation replicate (rows) of the current distinct model
# 'tmp.psi.boot' containing estimated mean true species probabilities for each bootstrap re-sample of each simulation replicate (rows) of the current distinct model
# 'output.boot' with summarized bootstrap estimates of standard errors and confidence interval coverage for estimates of mean true species probabilities for each distinct simulation model (rows)

# Returns data frame 'output.boot', appending summarized bootstrap estimates for the current distinct simulation model

psi.bootstrap.f <- function(sim.data, tmp.psi, tmp.psi.boot, tmp.theta.boot, output.boot) {
  
  if (B == 2) {
    ## ----- True species states (B) = 2 -----
    
    # Estimate SEs of estimated mean true species probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
    
    tmp.boot.se <- 
      group_by(as_tibble(tmp.psi.boot), rep) %>%
      summarise(., 
                se.bt.psi.1 = sd(psi.1, na.rm = TRUE)
      )
    
    tmp.boot.lg <-
      t(sapply(tmp.boot.se$rep, function(x)
        logit.ci.f(
          unlist(tmp.psi[x, 1:B]), 
          rep(tmp.boot.se[[which(tmp.boot.se$rep == x), 2]], 2)  ))
      )
    colnames(tmp.boot.lg) <- c("psi.1.btlg.L", "psi.1.btlg.U", "psi.2.btlg.L", "psi.2.btlg.U")
    
    # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes) 
    
    tmp.boot.lg.cov <- 
      mutate(as_tibble(tmp.boot.lg),
             ci.btlg.psi.1 = (psi.1.btlg.L < sim.data[[4]][[1]][1]) * (psi.1.btlg.U > sim.data[[4]][[1]][1])
      ) %>%
      select(., 5) %>%
      summarise_all(., mean)

    if (sim_profiles[1, ]$Model == "M.theta.psi" | sim_profiles$Model[1] == "M.theta+psi") {
      if (O_ps[1] > 0) {
        
        # ----- Multi-observation method (MOM) models (primary observers present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes).
        
        tmp.boot.se.theta <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.p.21 = sd(theta.p.21, na.rm = TRUE),
                    se.bt.theta.p.12 = sd(theta.p.12, na.rm = TRUE),
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE)
          )
        
        theta_col_tmp <- grep("theta", colnames(tmp.psi))[1:4]
        
        tmp.boot.lg.theta <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.psi[x, theta_col_tmp]), unlist(tmp.boot.se.theta[x, -1]))))
        colnames(tmp.boot.lg.theta) <- c(
          "theta.p.21.btlg.L", "theta.p.21.btlg.U", "theta.p.12.btlg.L", "theta.p.12.btlg.U",
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.12.btlg.L", "theta.s.12.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes) 
        
        tmp.boot.lg.cov.theta <- 
          mutate(as_tibble(tmp.boot.lg.theta),
                 ci.btlg.theta.p.21 = (theta.p.21.btlg.L < sim.data[[4]][[2]][1]) * (theta.p.21.btlg.U > sim.data[[4]][[2]][1]),
                 ci.btlg.theta.p.12 = (theta.p.12.btlg.L < sim.data[[4]][[2]][2]) * (theta.p.12.btlg.U > sim.data[[4]][[2]][2]),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][3]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][3]),
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][4]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][4])
          ) %>%
          select(., 9:12) %>%
          summarise_all(., mean)
        
      }else{
        # ----- Single-observation method (SOM) models (primary observers NOT present)
        
        # Estimate SEs of estimated mean classification probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes). 
        
        tmp.boot.se.theta <- 
          group_by(as_tibble(tmp.theta.boot), rep) %>%
          summarise(., 
                    se.bt.theta.s.21 = sd(theta.s.21, na.rm = TRUE),
                    se.bt.theta.s.12 = sd(theta.s.12, na.rm = TRUE)
          )
        
        theta_col_tmp <- grep("theta", colnames(tmp.psi))[1:2]
        
        tmp.boot.lg.theta <-
          t(sapply(1:reps, function(x)
            logit.ci.f(unlist(tmp.psi[x, theta_col_tmp]), unlist(tmp.boot.se.theta[x, -1]))))
        
        colnames(tmp.boot.lg.theta) <- c(
          "theta.s.21.btlg.L", "theta.s.21.btlg.U", "theta.s.12.btlg.L", "theta.s.12.btlg.U"
        )
        
        # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes)
        
        tmp.boot.lg.cov.theta <- 
          mutate(as_tibble(tmp.boot.lg.theta),
                 ci.btlg.theta.s.21 = (theta.s.21.btlg.L < sim.data[[4]][[2]][1]) * (theta.s.21.btlg.U > sim.data[[4]][[2]][1]), 
                 ci.btlg.theta.s.12 = (theta.s.12.btlg.L < sim.data[[4]][[2]][2]) * (theta.s.12.btlg.U > sim.data[[4]][[2]][2])  
          ) %>%
          select(., 5:6) %>%  
          summarise_all(., mean)
      }
      # Add estimates for classification probabilities (theta) to estimates for true species probabilities
      
      tmp.boot.se <- bind_cols(tmp.boot.se, tmp.boot.se.theta[, -1])
      tmp.boot.lg.cov <- bind_cols(tmp.boot.lg.cov, tmp.boot.lg.cov.theta)
    }
    
    # Return data frame 'output.boot' with results for current distinct model appended as a new row
    tmp.boot.se <- select(tmp.boot.se, -1) %>%
      summarise_all(., mean)
    output.boot <-
      bind_rows(output.boot,
                bind_cols(tmp.boot.lg.cov, tmp.boot.se))
    
    # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 

    # output.boot <- 
    #   bind_rows(
    #     output.boot, 
    #     bind_cols(
    #       tmp.boot.lg.cov[, 3],
    #       group_by(as_tibble(tmp.psi.boot), rep) %>%
    #         summarise(., 
    #                   se.bt.psi.1 = sd(psi.1, na.rm = TRUE),
    #                   psi.1.bt.L = quantile(psi.1, 0.025),
    #                   psi.1.bt.U = quantile(psi.1, 0.975)
    #         ) %>%
    #         mutate(.,
    #                ci.bt.psi.1 = (psi.1.bt.L < sim.data[[4]][[1]][1]) * (psi.1.bt.U > sim.data[[4]][[1]][1])
    #         ) %>%
    #         select(., 2, 5) %>%
    #         summarise_all(., mean)
    #     )
    #   )
    
  } else if (B == 3) {
    ## ----- True species states (B) = 3 -----
    
    # Estimate SEs of estimated mean true species probabilities (psi) as the SD of the bootstrap estimates, then construct 95% confidence limits by: applying logit transformation to parameters and SEs, constructing confidence limits assuming sampling distributions follow normal distributions, and applying inverse logit transformation to CLs (lower and upper limits denoted by 'btlg.L' and 'btlg.U' suffixes).  
    
    tmp.boot.se <- 
      group_by(as_tibble(tmp.psi.boot), rep) %>%
      summarise(., 
                se.bt.psi.1 = sd(psi.1, na.rm = TRUE),
                se.bt.psi.2 = sd(psi.2, na.rm = TRUE),
                se.bt.psi.3 = sd(psi.3, na.rm = TRUE)
      )

    tmp.boot.lg <-
      t(sapply(tmp.boot.se$rep, function(x)
        logit.ci.f(
          unlist(tmp.psi[x, 1:B]), 
          tmp.boot.se[which(tmp.boot.se$rep == x), 2:(B + 1)]  )))
    
    colnames(tmp.boot.lg) <- c("psi.1.btlg.L", "psi.1.btlg.U", 
                               "psi.2.btlg.L", "psi.2.btlg.U",
                               "psi.3.btlg.L", "psi.3.btlg.U")
    
    # Compute 95% confidence interval coverage for each estimate from confidence limits constructed from logit-transformed values ('ci.btlg.' prefixes) 
    
    tmp.boot.lg.cov <- 
      mutate(as_tibble(tmp.boot.lg),
             ci.btlg.psi.1 = (psi.1.btlg.L < sim.data[[4]][[1]][1]) * (psi.1.btlg.U > sim.data[[4]][[1]][1]),
             ci.btlg.psi.2 = (psi.2.btlg.L < sim.data[[4]][[1]][2]) * (psi.2.btlg.U > sim.data[[4]][[1]][2]),
             ci.btlg.psi.3 = (psi.3.btlg.L < sim.data[[4]][[1]][3]) * (psi.3.btlg.U > sim.data[[4]][[1]][3])
      ) %>%
      select(., 7:9) %>%
      summarise_all(., mean)
    
    # Return data frame 'output.boot' with results for current distinct model appended as a new row
    tmp.boot.se <- select(tmp.boot.se, -1) %>%
      summarise_all(., mean)
    output.boot <-
      bind_rows(output.boot,
                bind_cols(tmp.boot.lg.cov, tmp.boot.se))
    
    # Deprecated :: Also compute 95% confidence limits and CI coverage using the 2.5 and 97.5 percentiles of the bootstrap distributions (lower and upper limits denoted by 'bt.L' and 'bt.U' suffixes) and associated confidence interval coverage ('ci.bt.' prefix). Data frame 'output.boot' summarizes confidence interval coverage for each method. 
    
    # output.boot <- 
    #   bind_rows(
    #     output.boot,
    #     bind_cols(
    #       tmp.boot.lg.cov,
    #       group_by(as_tibble(tmp.psi.boot), rep) %>%
    #         summarise(
    #           .,
    #           se.bt.psi.1 = sd(psi.1, na.rm = TRUE),
    #           se.bt.psi.2 = sd(psi.2, na.rm = TRUE),
    #           se.bt.psi.3 = sd(psi.3, na.rm = TRUE),
    #           psi.1.bt.L = quantile(psi.1, 0.025),
    #           psi.2.bt.L = quantile(psi.2, 0.025),
    #           psi.3.bt.L = quantile(psi.3, 0.025),
    #           psi.1.bt.U = quantile(psi.1, 0.975),
    #           psi.2.bt.U = quantile(psi.2, 0.975),
    #           psi.3.bt.U = quantile(psi.3, 0.975)
    #         )
    #       %>%
    #         mutate(.,
    #                ci.bt.psi.1 = (psi.1.bt.L < sim.data[[4]][[1]][1]) * (psi.1.bt.U > sim.data[[4]][[1]][1]),
    #                ci.bt.psi.2 = (psi.2.bt.L < sim.data[[4]][[1]][2]) * (psi.2.bt.U > sim.data[[4]][[1]][2]),
    #                ci.bt.psi.3 = (psi.3.bt.L < sim.data[[4]][[1]][3]) * (psi.1.bt.U > sim.data[[4]][[1]][3])
    #         ) %>%
    #         select(., 2:4, 11:13) %>%
    #         summarise_all(., mean)
    #     ) 
    #   )
  } # End of bootstrap output loop
}

# Function: {mlogit.group.predict.f} Predicts group-level probabilities for true species probabilities (psi) or classification probabilities (theta) and computes bias in predicted versus true values. Accepts inputs of data frames with simulated survey data 'dat', simulation profile for the current distinct model 'profile', and summarized statistical output from model estimates for the current distinct model 'output'. Returns a vector with the overall mean and overall standard deviation of predicted group-level probabilities from estimated parameters relative to true parameters. 

mlogit.group.predict.f <- function(dat, profile, output) {
  
  B <- profile$B
  A <- profile$A
  
  if (profile$Model == "M.psi" | profile$Model == "M.theta.psi") {
    # Predict true species probabilities 
    
    # Extract covariate values
    covariate <- unlist(dat[, grep("covariate_psi", names(dat))])
    
    
    # Extract estimated 'beta.m' and true 'beta.t' regression coefficients (column 1/2 = intercepts/slopes)
    beta.m <- 
      matrix(c(
        unlist(output[grep("^m.b[01]_psi", names(output))  ]) 
      ),
      ncol = 2,
      byrow = TRUE)
    
    beta.t <- 
      matrix(c(
        unlist(profile[grep("^b[01]_psi", names(profile))]) 
      ),
      ncol = 2,
      byrow = TRUE)
    
    # Predict group-level probabilities from estimated 'psi.m' and true 'psi.t' parameters, find difference 'diff'
    psi.m <- mlogit.regress.predict.f(covariate, beta.m)
    psi.t <- mlogit.regress.predict.f(covariate, beta.t)
    diff <- psi.m - psi.t 
    if (B == 2) {diff <- diff[, 1, drop = FALSE]}
    
    # Compute mean and SD accounting for counts of each unique observation history
    n.diff <- (ncol(diff) * sum(dat$count))
    n.count <-
      matrix(dat$count, 
             nrow = length(dat$count), 
             ncol = ncol(diff))
    diff.mean <- sum(diff * n.count) / n.diff
    diff.sd <- ((sum(diff ^ 2 * n.count) - sum(diff * n.count) ^ 2 / n.diff) / (n.diff - 1)) ^ 0.5
    
    out_psi <- c(diff.mean, diff.sd)
    names(out_psi) <- c("m.psi.diff", "sd.psi.diff")
    
    if (profile$Model == "M.psi") {
      return(out_psi)
    }
    
  }
  
  if (profile$Model == "M.theta" | profile$Model == "M.theta.psi") {
    # Predict classification probabilities 
    
    # Extract covariate values
    covariate <- unlist(dat[, grep("covariate_theta", names(dat))])
    
    
    if (B == 2 & A == 2) {
      # True species states B = observation states A = 2
      
      # Extract estimated 'beta.m' and true 'beta.t' regression coefficients (column 1/2 = intercepts/slopes)
      beta.m <- 
        matrix(c(
          unlist(output[grep("^m.b[01]_theta", names(output))]) 
        ),
        ncol = 2)
      
      beta.t <- 
        matrix(c(
          unlist(profile[grep("^b[01]_theta", names(profile))]) 
        ),
        ncol = 2)
      
      # Predict group-level probabilities from estimated 'theta.m' and true 'theta.t' parameters, find difference 'diff'
      theta.m <- apply(beta.m, 1, function(x)
        mlogit.regress.predict.f(covariate, x)[, 1])
      
      theta.t <- apply(beta.t, 1, function(x)
        mlogit.regress.predict.f(covariate, x)[, 1])
      
      diff <- theta.m - theta.t 
      
      # Compute mean and SD accounting for counts of each unique observation history
      n.diff <- (ncol(diff) * sum(dat$count))
      n.count <-
        matrix(dat$count, 
               nrow = dim(diff)[1], 
               ncol = dim(diff)[2])
      
      diff.mean <- sum(diff * n.count) / n.diff
      diff.sd <- ((sum(diff ^ 2 * n.count) - sum(diff * n.count) ^ 2 / n.diff) / (n.diff - 1)) ^ 0.5
      
      out <- c(diff.mean, diff.sd)
      
      names(out) <- c("m.theta.diff", "sd.theta.diff")
      
      if (profile$Model == "M.theta.psi") {
        return(c(out_psi, out))
      } else 
        return(out)
      
    }else if (A == 3) {
      # Observation states A = 3
      
      # Extract estimated 'beta.m' and true 'beta.t' regression coefficients (column 1/2 = intercepts/slopes)
      beta.m <-
        cbind(matrix(c(unlist(output[grep("^m.b0_theta", names(output))])),
                     ncol = 2,
                     byrow = TRUE),
              matrix(c(unlist(output[grep("^m.b1_theta", names(output))])),
                     ncol = 2,
                     byrow = TRUE))
      
      beta.t <- 
        cbind(matrix(c(unlist(profile[grep("^b0_theta", names(profile))])),
                     ncol = 2,
                     byrow = TRUE),
              matrix(c(unlist(profile[grep("^b1_theta", names(profile))])),
                     ncol = 2,
                     byrow = TRUE))
      
      # Predict group-level probabilities from estimated 'theta.m' and true 'theta.t' parameters, find difference 'diff'
      theta.m <- apply(beta.m, 1, function(x)
        mlogit.regress.predict.f(
          covariate, 
          matrix(x, ncol = 2, byrow = FALSE) ))
      
      theta.t <- apply(beta.t, 1, function(x)
        mlogit.regress.predict.f(
          covariate, 
          matrix(x, ncol = 2, byrow = FALSE) ))
      
      diff <- theta.m - theta.t 
      
      # Compute mean and SD accounting for counts of each unique observation history
      n.diff <- ((dim(diff)[1] / length(dat$count)) * ncol(diff) * sum(dat$count))
      n.count <-
        matrix(dat$count, 
               nrow = dim(diff)[1], 
               ncol = dim(diff)[2])
      diff.mean <- sum(diff * n.count) / n.diff
      diff.sd <- ((sum(diff ^ 2 * n.count) - sum(diff * n.count) ^ 2 / n.diff) / (n.diff - 1)) ^ 0.5
      
      out <- c(diff.mean, diff.sd)
      
      names(out) <- c("m.theta.diff", "sd.theta.diff")
      
      if (profile$Model == "M.theta.psi") {
        return(cbind(out_psi, out))
      } else 
        return(out)
    }
  }
}

###############################################################################
#             Functions for error-checking simulation output
###############################################################################

# Function: {failure.check.f} Checks model output from simulation analyses for errors (convergence failure, missing parameter estimates, estimates outside user-specified boundary conditions). Accepts inputs of list of model optimization output from simulation analyses 'results' and simulation profile for the current distinct model 'sim', and returns list identifying failed models by simulation and id #s with explanatory error messages and codes for each (codes defined for 'optim' function and 98 = missing estimates, 99 = estimates outside user-specified boundary conditions). 

failure.check.f <- function(results, sim){
  # Check for missing parameter estimates, assign error code = 98

  if (is.null(results$par)) {
    fail <- 1:length(results$message)
    msg <-
      as_tibble(
        list(
          sim = rep(sim, length(fail)), 
          id = fail, 
          msg = results$message[fail], 
          code = 98
        )
      )
    return(list(fail, msg))
  }
  
  err <- which(!laply(results$par, function(x)
    is.numeric(x)))
  msg <- NULL
  
  # Check for convergence errors generated by 'optim' function
  
  if (length(err) > 0) {
    msg <-
      as_tibble(
        list(
          sim = rep(sim, length(err)), 
          id = err, 
          msg = results$par[err], 
          code = results$value[err]
        )
      )
    results[err, "convergence"] <- 1
  }
  
  # Check for parameter estimates outside user-defined boundary conditions, assign error code = 99
  boundary.err <- boundary.check.f(as.data.frame(results$par))
  
  if (length(boundary.err) > length(err)) {
    results[boundary.err, "convergence"] <- 1
    boundary.err <- setdiff(boundary.err, err)
    msg <- bind_rows(msg, 
                     as_tibble(list(
                       sim = rep(sim, length(boundary.err)),
                       id = boundary.err,
                       msg = results$par[boundary.err],
                       code = as.list(99)
                     ))
    )
  }
  
  # 'fail' contains ID number of all failed 'model.results' 
  # 'msg' contains relevant explanatory error messages
  fail <- which(results$convergence > 0)
  return(list(fail, msg))
}

# Function: {boundary.check.f} Check if estimated parameters fall outside user-specified boundary conditions. Accepts inputs of estimated parameters 'par' from simulation results and returns which models have parameters at or near user-specified boundary conditions. 

boundary.check.f <- function(par) {
  # Check multinomial regression coefficients (b0 = intercept, b1 = slope)
  if (sum(n_parameters[1:2]) > 0) {
    b0.fail <- which(!between(unlist(slice(par, 1:n_parameters[1])),
                          constraints_b0[1] + 0.1,
                          constraints_b0[2] - 0.1)) + n_parameters[1] - 1
    
    b1.fail <- which(!between(unlist(slice(par, (n_parameters[1] + 1):sum(n_parameters[1:2]))),
                             constraints_b1[1] + 0.1,
                             constraints_b1[2] - 0.1)) + n_parameters[2] - 1
  }else{
    b0.fail <- b1.fail <- numeric(0)
  }

  # Check multinomial logit parameters
  logit.fail <- which(!between(unlist(slice(par, (sum(n_parameters[1:2]) + 1):sum(n_parameters[1:3]))),
                               constraints_logit[1] + 0.5,
                               constraints_logit[2] - 0.5)) + n_parameters[3] - 1
  
  # Check mean group size parameters
  if (n_parameters[4] > 0) {
  g.fail <- which(!between(unlist(slice(par, (sum(n_parameters[1:3]) + 1):sum(n_parameters[1:4]))),
                           constraints_g[1] + 0.0000001,
                           constraints_g[2] - 0.5)) + n_parameters[4] - 1
  }else{
    g.fail <- numeric(0)
  }
  
  # Check heterogeneous group parameters
  if (n_parameters[5] > 0) {
    m.fail <- which(!between(unlist(slice(par, (sum(n_parameters[1:4]) + 1):sum(n_parameters[1:5]))),
                             constraints_mix[1] + 0.001,
                             constraints_mix[2])) + n_parameters[5] - 1
  }else{
    m.fail <- numeric(0)
  }
  # Return ID number of model results with parameters outside boundary conditions
  unique(c(b0.fail %/% n_parameters[1], b1.fail %/% n_parameters[2], logit.fail %/% n_parameters[3], g.fail %/% n_parameters[4], m.fail %/% n_parameters[5]))
}

# Function: {refit.check.f} Checks if models fit using alternate sets of initial parameter values (because optimization with initial set failed) successfully optimize. Accepts input of output from model optimization 'results' for simulation models that previously failed to optimize correctly, and returns list of results with error codes for failed results (convergence = 2 denotes missing parameter estimates, convergence = 3 denotes parameter values outside boundary conditions). 

refit.check.f <- function(results){
  # Check for missing parameter estimates, label with "convergence" = 2
  err.again <-
    which(!laply(results[[1]], function(x) is.numeric(x)))
  results[err.again, "convergence"] <- 2
  
  # Check for parameter estimates outside user-specified boundary conditions, label with "convergence" = 3
  if (any(results$convergence == 0)) {
    boundary.err <-
      boundary.check.f(as_tibble(results$par[which(results$convergence == 0)]))
    results[boundary.err, "convergence"] <- 3
  }
  return(results)
}

# Function: {print.error.f} Prints any error messages and convergence failure reports during profiles simulations to user console. Accepts inputs of error messages 'err.msg' and 'converge.fail', prints output messages to user console. 

print.error.f <- function(err.msg, converge.fail){
  if (length(err.msg) > 0) {
    cat(dim(err.msg)[1], "model(s) returned error messages: see 'error_msg' \n")
    cat("sim  error message \n")
    for (i in 1:min(dim(err.msg)[1], 20)) {
      cat(err.msg$sim[i], unlist(err.msg$msg[i]), "\n")
      if (i == 20) {
        cat("Only first 20 error messages displayed \n")
      }
    }
  }
  
  if (!is.null(converge.fail)) {
    cat(
      dim(converge.fail)[1], "model(s) failed to converge and were removed from analyses: see 'converge_fail' \n"
    )
    converge.fail  %>% dplyr::count(., sim = sim)  %>% print(.)
  }
}

###############################################################################
#             Functions for generating simulation data
###############################################################################

# Function: {names.obs.f} Generates column labels for simulated survey observations. Accepts inputs of # of true species states 'B', observation states 'A', and # of primary/secondary survey observers 'O_ps'. Returns vector of names with format y_[observer]_[A], where [observer takes values "p" for primary observer(s) and values "s1", "s2", "s3" and "s4" for secondary observers and [A] takes values 1 to A for observation states.

names.obs.f <- function(B, A, O_ps){
  
  if (O_ps[1] > 0) {
    p.names <- c(paste0("y_p_", c(1:B, "p")))
  } else {
    p.names <- NULL
  }

  s.names <- c(paste0(c(1:B, "p")))
  s.names <- sapply(1:O_ps[2], function(x) paste("y_s", x, "_", s.names[1:A], sep = ""))
  c(p.names[1:A], s.names)
}
