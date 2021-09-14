# simulations_R_code.R
# Simulation engine for estimating uncertain identification using multi-observer methods, version 1.2.0
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# R computer code used to conduct "omnibus", "covariate", and "distinct observer" simulation analyses using multi-observer methods. Simulations rely on input in .csv file format (described in DataS2), which can be altered to conduct user-specified simulations. For multiple-observation method (MOM) and single-observation (SOM) method models and 2 to 4 species, code generates and analyzes simulated survey data, with summary output written to .csv format files. Descriptions of statistical methods are provided in the companion article and Appendices S1 to S3. Descriptions of R objects containing simulated survey data ('sim_data'), statistical output ('model_results'), and summarized simulation output ('output_global') are provided in DataS3. Code developed and tested in R version 4.1.
 
###############################################################################
#           Required functions
###############################################################################

# Execute R code in these files to load required functions into the R environment
# generate_simulation_data.R		Function for generating survey data
# optimization_functions.R	    Functions for likelihood optimization
# supplemental_functions.R      Provides supplemental functions for analyses  

###############################################################################
#           Required R packages
###############################################################################

## Install and load prior to executing code below
# library(fastverse)  High-performance statistical computing and data manipulation
#   Utilized fastverse packages are 'collapse', 'matrixStats', & 'magrittr'
# library(plyr)       Programming tools, data manipulation
# library(dplyr)      Data frame manipulation
# library(doSNOW)     Back end for parallel code execution
# library(mrds)       Delta method for arbitrary functions
# library(extraDistr) Probability distributions
# library(nnet)       Log-linear multinomial models (only required for observed proportions models) 
# library(Rfast)      High-performance data analysis

library(fastverse)

fastverse_detach(data.table, kit)

fastverse_extend(plyr, dplyr, extraDistr, nnet, doSNOW, mrds, Rfast)

###############################################################################
#     Register the 'doSNOW' parallel execution back end
###############################################################################

# Make cluster 'cl' with user-specified # of CPU workers (in parentheses). At least 1 worker is required for the 'foreach' construct, while >1 allows parallel processing of simulations. Typically, the number of workers should be ≤ the # of available CPU threads.

cl <- makeCluster(4)
registerDoSNOW(cl)

# AFTER completing all simulation analyses, un-register the cluster 'cl' to remove CPU workers from RAM.
stopCluster(cl) 

###############################################################################
#          Import simulation profile(s)
###############################################################################
  
## Run this code in this section before each simulation
# Exception: the 'likelihood.equations' lists only need to be loaded at start or between simulation analyses when altering true species states B or classification states A

# User-specified .csv format input/output files for simulations
# Do not include .csv file extension
in_filename <- "omnibus_m_22_example"
out_filename <- paste0(in_filename, "_output.csv")

# Vectors 'sim.const'/'sim_variables' contain simulation inputs that are constant/variable across distinct models
sim_constants <- as.list(read.csv(paste0(in_filename, ".csv"), header = TRUE, skip = 4, nrows = 1, as.is = TRUE))
sim_variables <- read.csv(paste0(in_filename, ".csv"), header = TRUE, skip = 9, as.is = TRUE)
sim_constants <- sim_constants[!is.na(sim_constants)]
sim_variables <- sim_variables[, !is.na(sim_variables[1, ])]
sim_constants <- format.sim.constant.f(sim_constants)

# Data frame 'sim_profiles' defines simulations, with 1 distinct model per row
sim_profiles <- generate.sim.profiles.f(sim_constants, sim_variables) 
profiles_names <- names(sim_profiles) # Field names

# Vector 'parameter_names_global' contains parameters names across all models
parameters_names_global <-
  c(grep("b0|b1|psi|theta", profiles_names, value = TRUE))
if (any(sim_profiles[grep("^g_[0123456789]", profiles_names)] > 1)) {
  parameters_names_global <-
    c(parameters_names_global,
      grep("^g_[0123456789]", profiles_names, value = TRUE))
}
mix_col <- grep("mix", profiles_names) # Column #s of parameters for estimating heterogeneous groups
if (any(sim_profiles[, mix_col] > 0)) {
  parameters_names_global <- c(parameters_names_global, profiles_names[mix_col])
}

# Set random number seed for repeatable results
if (any(colnames(sim_profiles) == "seed"))
  {set.seed(eval(parse(text = sim_profiles$seed[1]))) ; sim_data_list <- NULL}

model <- sim_profiles$Model[1] # Model specification
A <- sim_profiles$A[1]; B <- sim_profiles$B[1] # # of observation (A) and true species (B) states 
reps <- sim_profiles$reps[1] # Simulation replicates for each distinct model
n_bootstrap <- sim_profiles$n_bootstrap[1] # Number of bootstrap resamples (optional for models with covariates)

# Optional user-specified expression for upper constraints on classification probabilities (theta). See "MetadataS2.pdf" in DataS2-Input-Files-for-Simulations and "General Simulation Methods" in Appendix S3 for details.
theta_limit <- parse(text = sim_profiles[1, ]$t_limit)

# Function: {theta.limit.f} Defines upper box constraints for values of classification probabilities (theta) from user-specified expression
theta.limit.f <- function(theta) {
  min(eval(theta_limit), 0.999)
}

# Function: {theta4.f} Closure function producing functions that build matrices of classification probabilities (theta). See function 'theta.calc.f' in 'supplemental.functions.R' for details.
theta4.f <- theta.calc.f(4, B, A)

# Import appropriate 'likelihood.equations' R list object containing pre-computed values for likelihood computations
lookup_path <-
  switch(
    paste0(B, A),
    '22' = 'likelihood_equations_a2b2g25.RData',
    '23' = 'likelihood_equations_a3b2g25.RData',
    '33' = 'likelihood_equations_a3b3g18.RData',
    '44' = 'likelihood_equations_a4b4g10.RData'
  )
load(lookup_path)

###############################################################################
#            A few technical details for MOM/SOM simulations
###############################################################################

## Error files
# All MOM/SOM simulation routines create files ‘error_msg’ and ‘converge_fail’ containing error codes and output generated by the ‘optim’ function and also by internal error-checking functions (98 = “missing estimates”, 99 = “values outside user-specified boundaries”). 

## Limits for input parameters
# Code for models currently support:
# 0-1 primary observer and 1-4 secondary observers
# 1-2 observation methods
# 2-4 observation states and true species states
# For models including one covariate predicting psi or theta, 2-3 observation states and true species states
# For models including independent covariates predicting psi and theta, 2 observation states and true species states
# Heterogeneous groups include 2 species

## User-specified control parameters for 'optim' function are specified in the list 'cont' 
# 'factr' and 'maxit' set the required precision for likelihood optimization and max # of iterations
# 'lmm' (usually between 5-17) influences memory allocation and possibly optimization speed. See Zhu et al. 1997. Algorithm 778: L-BFGS-B: FORTRAN Subroutines for Large-Scale Bound-Constrained Optimization. ACM Transactions on Mathematical Software 23(4):550–560

## Constraints for estimated parameters
# The 'L-BFGS-B' optimization method allows box constraints specifying constraints on values for individual estimated parameters. Lower and upper constraints are summarized in the vectors 'constraints_low' and 'constraints_up'. Vectors (below) define user-specified values for lower and upper constraints by parameter type. Lower constraints often help optimization and are included by default. Upper constraints should only be included when necessary. Vector 'constraints_up' defaults to values of 'Inf' (e.g., no constraints). If users specify 'theta limits' in .csv input files (see MetadataS2 in DataS2) upper constraints will be used for model optimization, with upper constraints for classification probabilities (theta) defined by the 'theta limits' and upper constraints for other parameters taken from input vectors below. If estimated parameters fall very near or outside constraints, model optimization will fail and generate error code '99'.
# Constraints typically can be left unaltered, but may need to be changed to avoid model optimization problems occurring when true values are outside constraints or when models optimize on (extremely) incorrect values (e.g. multi-modality). 

# 'constraints_logit' Constraints for true species probabilities and classification probabilities (psi/theta)
# 'constraints_g' Constraints for mean group size
# 'constraints_mix' Constraints for heterogeneous group parameters (pi and rho)
# 'constraints_b0' and 'constraints_b1' Constraints for multinomial logit regression coefficients (betas, b0 = intercepts, b1 = slopes) predicting true species probabilities and classification probabilities (psi/theta)

###############################################################################
#           MOM/SOM model simulations 
###############################################################################

# These routines conduct "omnibus", "covariate", and "distinct observer" simulations for MOM/SOM models, assuming no un-modeled heterogeneity in the data (i.e., models for data generation and data analyses are identical). 
# Conduct analyses by selecting and executing all code in this section. 

if (model == "M") {
  # ------ MOM/SOM models without covariates ----------
  
  ## Define global objects (constant across all simulations)
  
  # For model optimization function 'optim', define user-specified control parameters and lower/upper constraints for parameter values 
  cont <- list(factr = 1e+07, maxit = 200, lmm = 12) # List of control parameters for 'optim' function
  constraints_logit <- qlogis(c(0.0001, 0.9999)) # Constraints for psi/theta parameters, transformed to logit scale
  constraints_g <- c(1.002, Inf) # Constraints for estimated mean group size
  constraints_mix <- c(qlogis(0.001), Inf) # Constraints for group mixing parameters
  
  t_start <- t_loop <- unclass(Sys.time()) # For recording simulation duration
  error_msg <- converge_fail <- NULL # For error messages and convergence warnings
  
  # Build empty data frame 'output_global' containing columns for global estimated parameters (across all simulations) and associated output 
  output_global <- matrix(numeric(1), 0, length(parameters_names_global) * 7 + 1) 
  output_names_global <-
    c(
      paste0("m.", parameters_names_global),
      paste0("sd.", parameters_names_global),
      paste0("se.", parameters_names_global),
      paste0("cv.", parameters_names_global),
      paste0("me.", parameters_names_global),
      paste0("rm.", parameters_names_global),
      paste0("ci.", parameters_names_global),
      "rep.tru"
    )
  colnames(output_global) <- output_names_global
  output_global <- data.frame(output_global)
  
  # Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
  for (sim in 1:dim(sim_profiles)[1]) { 
    
    ## Define objects for simulation of the current distinct model
    
    O_ps <- c(sim_profiles$O_p[sim], sim_profiles$O_s[sim]) # Count of primary and secondary observers
    
    # Vector 'n_parameters' contains the number of estimated parameters for the current simulation in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), un-transformed parameters (4), and heterogeneous group probabilities (pi) parameters (5)
    # Vectors 'parameters_names' and 'parameters_col' contain names and column #s of true parameter values in 'sim_profiles'
    
    # Add logit parameters for uncertain identification probabilities (theta) and true species probabilities (psi)
    n_parameters <- c(0, 0, length(grep("b0|b1|psi|theta", parameters_names_global)), 0, 0)
    parameters_col <- grep("b0|b1|psi|theta", profiles_names) # Column #s of parameters for estimating psi & theta
    parameters_names <- c( grep("b0|b1|psi|theta", parameters_names_global, value = TRUE))
    
    # Add parameters for mean group size g
    if (any(sim_profiles[sim, grep("^g_[0123456789]", profiles_names)] > 1)) {
      parameters_names <- c(parameters_names, grep("^g_[0123456789]", parameters_names_global, value = TRUE))
      parameters_col <- c(parameters_col, grep("^g_[0123456789]", profiles_names))
      n_parameters[4] <- length(parameters_names) - n_parameters[3]
    }
    
    # Add parameters for mixed groups: heterogeneous group probabilities (pi)
    if (any(sim_profiles[sim, mix_col] > 0)) {
      parameters_names <- c(parameters_names, profiles_names[mix_col])
      parameters_col <- c(parameters_col, grep("mix", profiles_names))
      n_parameters[5] <- length(mix_col)
    }
    
    # With certain ID by secondary observers (probability of uncertain ID = 0), remove these classification probabilities from estimated parameters
    if (any(sim_profiles[sim, c(grep("theta_s", profiles_names))] == 0)) {
      id_certain_col <- grep("theta_s", parameters_names)
      parameters_names <- parameters_names[-id_certain_col]
      parameters_col <- parameters_col[-id_certain_col]
      n_parameters[3] <- n_parameters[3] - length(id_certain_col)
      sim_profiles$O_s[sim] <- O_ps[2] <- 1   # Set # of secondary observers = 1
    }
    
    ## Generate initial objects 
    parameters_ini <- generate.ini.f(sim_profiles[sim, ], n_parameters) # Sets of alternative initial parameter values (rows)
    output_names <-
      c(
        paste0("m.", parameters_names),
        paste0("sd.", parameters_names),
        paste0("se.", parameters_names),
        paste0("cv.", parameters_names),
        paste0("me.", parameters_names),
        paste0("rm.", parameters_names),
        paste0("ci.", parameters_names),
        "rep.tru"
      )
    
    # Vector 'constraints_low' contains lower boundaries for each parameter
    constraints_low <- c(rep(constraints_logit[1], n_parameters[3]), 
                         rep(constraints_g[1], n_parameters[4]), 
                         rep(constraints_mix[1], n_parameters[5]))
    
    # Optional upper parameter boundaries default to infinity (no constraints) 
    constraints_up <- Inf
    
    # If present, apply user-specified upper constraints for classification probabilities (theta) exist in 'theta_limit'
    if (exists("theta_limit")) {
      if (eval(theta.limit.f(0.1)) > 0) {
        constraints_up <-
          c(rep(constraints_logit[2], n_parameters[3]),
            rep(constraints_g[2], n_parameters[4]),
            rep(constraints_mix[2], n_parameters[5]))
        
        theta_col <- c(grep("theta", parameters_names, value = TRUE))
        constraints_up[1:(length(theta_col))] <- 
          qlogis(
            unlist(lapply(sim_profiles[sim, theta_col], function(x)
              theta.limit.f(x)))
          )
      }
    }
    
    # List 'sim_data' contains formatted simulated MOM survey observations and (optionally) a keyed table of unique group records
    sim_data <- generate.simulation.data.f(sim_profiles[sim,]) 
    parameters_index <- 1 # Index for alternative sets of initial parameter values in 'parameters_ini'

  # ----- MOM/SOM Model Optimization ----- 
    
  # Function for simulation progress bar
  id_max <- max(sim_data[[1]]$id)
  pb <- txtProgressBar(max = id_max, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Fit models via the 'optimize.M.f' function utilizing parallel execution of the 'optim' routine using the "L-BFGS-B" method, with statistical output to list 'model_results'
  model_results <-
    as.data.frame(foreach(
      i = 1:id_max, 
      .combine = rbind, 
      .packages = c("plyr", "dplyr", "extraDistr"), 
      .errorhandling = "pass", 
      .options.snow = opts
    ) %dopar% {
      optim(
        parameters_ini[parameters_index, ],
        optimize.M.f,
        gr = NULL,
        dat = sim_data[[1]][which(sim_data[[1]]$id == i), ],
        keys = sim_data[[2]],
        sim_profile = sim_profiles[sim,],
        hessian = TRUE,
        method = c("L-BFGS-B"),
        lower = constraints_low,
        upper = constraints_up,
        control = cont
      )
    })

  # ----- Model Output: Error-checking -------  
  
  # Check for model optimization errors and convergence failures
  fail <- failure.check.f(model_results, sim)
  
  # Attempt to fit non-optimized models using alternate initial parameter values
  if (length(fail[[1]]) > 0) {
    error_msg <- bind_rows(error_msg, fail[[2]])
    fail_tmp <- fail <- fail[[1]]

    # Try alternate sets of initial parameter values until model(s) optimize successfully
    while (parameters_index < 5 && length(fail) > 0) {
      parameters_index <- parameters_index + 1 # Increment index for alternate sets of initial values
      # >1 models necessary to use 'foreach' function
      if (length(fail) == 1) {
        fail_tmp <- c(fail, fail) 
      }
      # Attempt optimization with alternate initial parameter values
      model_results_refit <-
        as.data.frame(foreach(
          i = 1:length(fail_tmp), 
          .combine = rbind, 
          .packages = c("plyr", "dplyr", "extraDistr"), 
          .errorhandling = "pass"
        ) %dopar% {
          optim(
            parameters_ini[parameters_index,], 
            optimize.M.f, gr = NULL, 
            dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]),], 
            keys = sim_data[[2]], 
            sim_profile = sim_profiles[sim, ],
            hessian = TRUE, 
            method = c("L-BFGS-B"), 
            lower = constraints_low, 
            upper = constraints_up, 
            control = cont
          )
        })
      
      # Add model(s) successfully optimizing to 'model_results'
      model_results_refit <- refit.check.f(model_results_refit)
      success <- which(model_results_refit[1:length(fail), "convergence"] == 0)
      if (length(success) > 0) {
        l_ply(success, function(x) {
          for (i in 1:6)
            model_results[[i]][fail[x]] <<- model_results_refit[[i]][x]
        })
        fail <- fail[-success] # Remove successful models from failed models
      }
    }
    
    # If any model(s) never optimize adequately: generate warnings, remove failed model(s) from 'model_results', and update # of (successful) simulation replicates
    if (length(fail) > 0) {
      model_results <- model_results  %>% do(model_results[-fail, 1:6])
      
      converge_fail <- rbind(converge_fail, as_tibble(list(
        sim = rep(sim, length(fail)), 
        id = fail, 
        code = model_results_refit$convergence[which(model_results_refit$convergence > 0)][1:length(fail)], 
        msg = model_results_refit$message[which(model_results_refit$convergence > 0)][1:length(fail)]
      )))
      
      reps <- reps - length(fail)
    }  
  }
  
  ## ----- Model Output: Summarize ----- 
  
  ## Back-transform estimates for logit parameters to probability scale and calculate summary statistics. Prefixes on variable names indicate summary statistics: m =  mean, sd = standard deviation, cv = coefficient of variation (sd / mean), me = mean error, rm = residual mean squared error, ci = 95% confidence interval coverage
  
  # Vector 'parameters_true' contains true parameter values
  parameters_true <- sim_profiles[sim, parameters_col]
  
  # Append summary statistics for this distinct model (vector 'output') as new row in data frame 'output_global'
  parameters_output <- model.results.f(model_results, parameters_true, reps, n_parameters)
  output <- summarise.results.f(parameters_output, parameters_true, reps)
  names(output) <- output_names
  output_global <- bind_rows(output_global, data.frame(t(output)))
  
  # For long simulation runs, prevent data loss by periodically writing results to a temporary file on hard drive 
  if (unclass(Sys.time()) - t_loop > 900) {
    try(
      write.csv(output_global, file = paste0("tmp_", out_filename))
    )
    t_loop <- unclass(Sys.time())
  }
  
  # Complete simulation loop for the current distinct model 
  reps <- sim_profiles$reps[1]
  cat("\n Completed Simulation ", sim, "\n")
  
  # With random seed assigned, append current simulation data to list
  if (any(profiles_names == "seed"))
    sim_data_list <- c(sim_data_list, list(sim_data))
  
  } # End of profile simulation loop
  
  ## ----- Output & Reports ---
  
  # Combine simulation inputs and outputs, and write 'output_global' as .csv format file to user-specified path
  output_global <- bind_cols(sim_profiles, output_global)
  try(
    write.csv(output_global, file = out_filename)
  )
  
  # With random seed assigned, export simulation data
  if (any(profiles_names == "seed"))
    if (is.integer(sim_profiles$seed[1]) & sim_profiles$seed[1] > 0)
      try(saveRDS(sim_data_list, file = paste0(in_filename, "_output", "_sim_data.Rdata"))
      )
  
  # Print reports on simulation speed and errors/warnings to the console
  elapsed <- unclass(Sys.time()) - t_start
  cat(round(elapsed / 60, 2),
      "m ",
      round(elapsed, 0),
      "s ",
      round(elapsed / (reps * dim(sim_profiles)[1]), 2),
      "s/rep \n")
  print.error.f(error_msg, converge_fail)
  
  # End of model M analysis loop
  
} else if (model == "M.theta.p" | model == "M.theta" | model == "M.theta.ps") {
  # ----- MOM/SOM models with covariates predicting classification probabilities (theta) -----
  # True states B = (2 or 3), classification states A = (2 or 3)
  
  ## Define global objects (constant across all simulations)
  
  # For model optimization function 'optim', define user-specified control parameters and lower/upper constraints for parameter values 
  cont <- list(factr = 1e+07, maxit = 200, lmm = 12) # List of control parameters for 'optim' function
  constraints_logit <- qlogis(c(0.0001, 0.9999)) # Constraints for psi/theta parameters, transformed to logit scale
  constraints_g <- c(1.002, 25) # Constraints for mean group size
  constraints_mix <- qlogis(c(0.0001, 0.8)) # Constraints for heterogeneous group parameters
  constraints_b0 <- c(-10, 10); constraints_b1 <- c(-5, 5) # Constraints for regression coefficients 
  
  t_start <- t_loop <- unclass(Sys.time()) # For recording simulation duration
  error_msg <- converge_fail <- NULL # For error messages and convergence warnings
  
  # Define objects which will have output appended after each simulation replicate
  betas <- NULL
  output_betas <- NULL
  output_bootstrap <- NULL
  theta_diff <- NULL

  # Build empty data frame 'output_global' containing columns for global estimated parameters (across all simulations) and associated output 
  output_global <- matrix(numeric(1), 0, length(parameters_names_global) * 7 + 1) 
  output_names_global <-
    c(
      paste0("m.", parameters_names_global),
      paste0("sd.", parameters_names_global),
      paste0("se.", parameters_names_global),
      paste0("cv.", parameters_names_global),
      paste0("me.", parameters_names_global),
      paste0("rm.", parameters_names_global),
      paste0("ci.", parameters_names_global),
      "rep.tru"
    )
  colnames(output_global) <- output_names_global
  output_global <- data.frame(output_global)
  
  # Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
  for (sim in 1:dim(sim_profiles)[1]) { 
    
    ## Define objects for simulation of the current distinct model
    
    O_ps <- c(sim_profiles$O_p[sim], sim_profiles$O_s[sim]) # Count of primary and secondary observers
    
    # Vector 'n_parameters' contains the number of estimated parameters for the current simulation in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), un-transformed parameters (4), and heterogeneous group probabilities (pi) parameters (5)
    # Vectors 'parameters_names' and 'parameters_col' contain names and column #s of true parameter values in 'sim_profiles'
    
    parameters_col <- grep("b0|b1|psi|theta", profiles_names) 
    parameters_names <- c( grep("b0|b1|psi|theta", parameters_names_global, value = TRUE))
    
    # Add regression coefficients (b0 = intercept, b1 = slope) and logit parameters for classification probabilities (theta)
    n_parameters <-
      c(length(grep("b0", parameters_names)), 
        length(grep("b1", parameters_names)), 
        length(grep("theta|psi", parameters_names)) - length(grep("b0|b1", parameters_names)), 
        0, 0)
    
    # Add parameters for mean group size g
    if (any(sim_profiles[sim, grep("^g_[0123456789]", profiles_names)] > 1)) {
      parameters_names <- c(parameters_names, grep("^g_[0123456789]", profiles_names, value = TRUE))
      parameters_col <- c(parameters_col, grep("^g_[0123456789]", profiles_names))
      n_parameters[4] <- length(parameters_names) - sum(n_parameters[1:3])
    }
    
    # Add parameters for mixed spp groups: heterogeneous group probabilities (pi) or heterogeneous group affinities (rho)
    if (any(sim_profiles[sim, mix_col] > 0)) {
      parameters_names <- c(parameters_names, profiles_names[mix_col])
      parameters_col <- c(parameters_col, grep("mix", profiles_names))
      n_parameters[5] <- length(mix_col)
    }
    
    # With certain ID by secondary observers (probability of uncertain ID = 0), remove these classification probabilities from estimated parameters
    if (any(sim_profiles[sim, c(grep("theta_s", profiles_names), 3)] == 0)) {
      id_certain_col <- grep("theta_s", parameters_names)
      parameters_names <- parameters_names[-id_certain_col]
      parameters_col <- parameters_col[-id_certain_col]
      n_parameters[3] <- n_parameters[3] - length(id_certain_col)
      sim_profiles$O_s[sim] <- O_ps[2] <- 1  # Set secondary observers = 1 
      theta_limit <- NULL  # Theta limits unneeded with certain ID
    }
    
    # Generate initial objects 
    parameters_ini <- generate.ini.f(sim_profiles[sim, ], n_parameters) # Sets of alternative initial parameter values (rows)
    output_names <-
      c(
        paste0("m.", parameters_names),
        paste0("sd.", parameters_names),
        paste0("se.", parameters_names),
        paste0("cv.", parameters_names),
        paste0("me.", parameters_names),
        paste0("rm.", parameters_names),
        paste0("ci.", parameters_names),
        "rep.tru"
      )
    
    # Vector 'constraints_low' contains lower boundary conditions for each parameter
    constraints_low <- c(rep(constraints_b0[1], n_parameters[1]),  
                         rep(constraints_b1[1], n_parameters[2]),  
                         rep(constraints_logit[1], n_parameters[3]), 
                         rep(constraints_g[1], n_parameters[4]), 
                         rep(constraints_mix[1], n_parameters[5]))
    
    # Optional upper parameter boundaries default to infinity (no constraints) 
    constraints_up <- Inf

    # If present, apply user-specified upper constraints for classification probabilities (theta) exist in 'theta_limit'
    if (!is.null(theta_limit) & sim_profiles[1, ]$Model == "M.theta.p") {
      if (eval(theta.limit.f(0.1)) > 0) {
        constraints_up <- c(
          rep(constraints_b0[2], (B - 1)),
          rep(constraints_b1[2], (B - 1)),
          rep(constraints_b0[2], (n_parameters[1] - (B - 1))),
          rep(constraints_b1[2], (n_parameters[2] - (B - 1))),
          rep(constraints_logit[2], n_parameters[3]),
          rep(constraints_g[2], n_parameters[4]),
          rep(constraints_mix[2], n_parameters[5])
        )
        theta_col <- c(grep("^theta", profiles_names))
        constraints_up[(sum(n_parameters[1:2]) + 1):(sum(n_parameters[1:2]) + length(theta_col))] <- 
          qlogis(
            unlist(lapply(sim_profiles[sim, theta_col], function(x)
              theta.limit.f(x)))
          )
      }
    }
    
    # List 'sim_data' contains formatted simulated survey observations and (optionally) a keyed table of unique groups records
    sim_data <- generate.simulation.data.f(sim_profiles[sim,]) 
    parameters_index <- 1 # Index for alternative sets of initial parameter values in 'parameters_ini'
    
    # Used for optional output for beta parameters
    sim_data_tmp <- sim_data[[1]]
    
    # ----- MOM/SOM Model Optimization ----- 
    
    # Function for simulation progress bar
    id_max <- max(sim_data[[1]]$id)
    pb <- txtProgressBar(max = id_max, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Fit models via function 'optimize.M.theta.p.f' or 'optimize.M.theta.f'  utilizing parallel execution of the 'optim' routine using the "L-BFGS-B" method, with statistical output to list 'model_results'
 
    model_results <-
      as.data.frame(foreach(
        i = 1:id_max, 
        .combine = rbind, 
        .packages = c("plyr", "dplyr", "collapse", "Rfast"), 
        .errorhandling = "pass", 
        .options.snow = opts
      ) %dopar% {
        if (sim_profiles$Model[1] == "M.theta.p") {
          # Optimization for models with logit regression predicting classification probabilities only for primary observers
          optim(
            parameters_ini[parameters_index, ], 
            optimize.M.theta.p.f, gr = NULL, 
            dat = sim_data[[1]][which(sim_data[[1]]$id == i),], 
            sim_profile = sim_profiles[sim, ],
            hessian = TRUE, method = c("L-BFGS-B"), 
            lower = constraints_low, 
            upper = constraints_up, 
            control = cont
          )
        }else if (sim_profiles$Model[1] == "M.theta" | model == "M.theta.ps") {
          # Optimization for models with logit regression predicting classification probabilities for all observers
          optim(
            parameters_ini[parameters_index, ], 
            optimize.M.theta.f, gr = NULL, 
            dat = sim_data[[1]][which(sim_data[[1]]$id == i),], 
            keys = sim_data[[2]],
            sim_profile = sim_profiles[sim, ],
            hessian = TRUE, 
            method = c("L-BFGS-B"), 
            lower = constraints_low, 
            control = cont
          )
        }
      }) 
    
    # ----- Model Output: Error-checking -------  
    
    # Check for model optimization errors/convergence failures
    fail <- failure.check.f(model_results, sim)
    
    # Attempt to fit non-optimized models using alternate initial parameter values
    if (length(fail[[1]]) > 0) {
      error_msg <- bind_rows(error_msg, fail[[2]])
      fail_tmp <- fail <- fail[[1]]
      
      # Try alternate sets of initial parameter values until model(s) optimize successfully
      while (parameters_index < 5 && length(fail) > 0) {
        parameters_index <- parameters_index + 1 # Increment index for alternate sets of initial values
        # >1 models necessary to use 'foreach' function
        if (length(fail) == 1) {
          fail_tmp <- c(fail, fail)
        }
        
        # Attempt optimization with alternate initial parameter values
        model_results_refit <-
          as.data.frame(foreach(
            i = 1:length(fail_tmp), 
            .combine = rbind, 
            .packages = c("plyr", "dplyr"), 
            .errorhandling = "pass"
          ) %dopar% {
            if (sim_profiles$Model[1] == "M.theta.p") {
              optim(
                parameters_ini[parameters_index ,], 
                optimize.M.theta.p.f, 
                gr = NULL, 
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]),], 
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                upper = constraints_up, 
                control = cont
              )
            }else if (sim_profiles$Model[1] == "M.theta" | model == "M.theta.ps") {
              optim(
                parameters_ini[parameters_index ,], 
                optimize.M.theta.f, 
                gr = NULL, 
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]),], 
                keys = sim_data[[2]],
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                upper = constraints_up, 
                control = cont
              )
            }
          })
        
        # Add model(s) successfully optimizing to 'model_results'
        model_results_refit <- refit.check.f(model_results_refit)
        success <- which(model_results_refit[1:length(fail), "convergence"] == 0)
        if (length(success) > 0) {
          l_ply(success, function(x) {
            for (i in 1:6)
              model_results[[i]][fail[x]] <<- model_results_refit[[i]][x]
          })
          fail <- fail[-success] # Remove successful models from failed models
        }
      }
      
      # If any model(s) never optimize adequately: generate warnings, remove failed model(s) from 'model_results', and update # of (successful) simulation replicates
      if (length(fail) > 0) {
        model_results <- model_results  %>% do(model_results[-fail, 1:6])
        converge_fail <- rbind(converge_fail, as_tibble(list(
          sim = rep(sim, length(fail)), 
          id = fail, 
          code = model_results_refit$convergence[which(model_results_refit$convergence > 0)][1:length(fail)], 
          msg = model_results_refit$message[which(model_results_refit$convergence > 0)][1:length(fail)]
        )))
        
        # For optional for beta parameters, remove data for non-optimized models
        for (i in fail) {
          sim_data_tmp <- filter(sim_data_tmp, id != i)
        }
        reps <- reps - length(fail)
      }  
    } # End of model re-fitting loop
    
    ## ----- Model Output: Summarize ----- 
    
    ## Back-transform estimates for logit parameters to probability scale and calculate summary statistics. Prefixes on variable names indicate summary statistics: m =  mean, sd = standard deviation, cv = coefficient of variation (sd / mean), me = mean error, rm = residual mean squared error, ci = 95% confidence interval coverage
    
    # Vector 'parameters_true' contains true parameter values
    parameters_true <- sim_profiles[sim, parameters_col]
    
    # Append summary statistics for this distinct model (vector 'output') as new row in data frame 'output_global'
    parameters_output <- model.results.f(model_results, parameters_true, reps, n_parameters)
    output <- summarise.results.f(parameters_output, parameters_true, reps)
    names(output) <- output_names
    output_global <- bind_rows(output_global, data.frame(t(output)))
    
    # Optional output summarizing estimated versus true group-level predictions for classification probabilities. Un-comment the following 2 lines to add results to data frame 'output'. 
    # tmp <- mlogit.group.predict.f(sim_data[[1]], sim_profiles[sim, ], output)
    # theta_diff <- bind_rows(theta_diff, as_tibble(t(as.matrix(tmp))))
    
    ## Output related to estimates of overall mean classification probabilities (theta) derived from regression coefficients of multinomial logistic regression and observed covariate values. Alternative 95% confidence limits for estimated means were constructed assuming 1) asymptotic normal distribution in natural scale (i.e, probabilities) and 2) asymptotic normal distribution of parameters in logit scale.
    # For each distinct model (rows), 'output_betas' contains summarized statistics (mean error, root mean square error, and 95% CI coverage) for overall mean estimates of classification probabilities
    # For each simulation replicate (rows), 'betas' contains estimated regression coefficients, estimated overall mean probabilities with associated standard errors, 95% confidence intervals ('lg' denotes limits estimated in logit scale), and 95% confidence interval coverage
    
      tmp <- theta.derived.f(sim_data, sim_data_tmp, model_results, parameters_output, output_betas)
      output_betas <- tmp[[1]]
      betas <- tmp[[2]]
      
    ## Loop for bootstrap replications
    # For each bootstrap re-sample, data are re-sampled with replacement from simulated survey data, models are fit to re-sampled data, and estimated parameters from re-sampled data are used to estimate mean values for classification probabilities (theta) across observed covariate values
      
    # Bootstrap only if more than 1 bootstrap re-samples (n_bootstrap > 1) are specified in simulation inputs  
    if (n_bootstrap > 1) {
      remove_col <- grep("id|count|key", colnames(sim_data[[1]]))
      tmp_theta_bootstrap <- NULL
      tmp_id <- unique(sim_data_tmp$id)
      parameters_index <- 1
      cat("\n Completed model optimization for simulation", sim, "\n")

      # Loop for bootstrap re-samples for each replicate (indexed by 'q') of simulated data
      for (q in 1:reps) {
        
        # Simulated survey data for current simulation replicate 
        data.rep <- filter(sim_data_tmp, id == tmp_id[q])
        
        # Sample groups (with replacement) from simulated data to generate n (= n_bootstrap) re-samples of the simulated survey data
        tmp <-
          ldply(1:n_bootstrap, function(x)
            sample_n(data.rep, size = sum(data.rep$count), replace = TRUE, weight = count)) %>%
          select(., !any_of(remove_col)) %>%
          bind_cols(., id = rep(1:n_bootstrap, each = sum(data.rep$count)))

        # List 'data.bootstrap' contains re-sampled data for the current simulation replicate, with formatting is identical to 'sim_data'
        data.bootstrap <-
          format.MOM.data.f(tmp, A, O_ps, sim_profiles[sim, mix_col], sim_profiles[sim, ]$n_bins)

        # Fit models to re-sampled data using the appropriate function and parallel optimization of the 'optim' routine
        pb <- txtProgressBar(max = n_bootstrap, style = 3)
        model_results_bootstrap <-
          as.data.frame(foreach(
            i = 1:n_bootstrap,
            .combine = rbind,
            .packages = c("plyr", "dplyr"),
            .errorhandling = "pass", 
            .options.snow = opts
          ) %dopar% {
            if (sim_profiles$Model[1] == "M.theta.p") {
              # Optimization for models with logit regression predicting classification probabilities only for primary observers
              optim(
                parameters_ini[parameters_index, ], 
                optimize.M.theta.p.f, 
                gr = NULL,
                dat = data.bootstrap[[1]][which(data.bootstrap[[1]]$id == i),],
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                upper = constraints_up, 
                control = cont
              )
            }else if (sim_profiles$Model[1] == "M.theta" | model == "M.theta.ps") {
              # Optimization for models with logit regression predicting classification probabilities for all observers
              optim(
                parameters_ini[parameters_index, ], 
                optimize.M.theta.f, 
                gr = NULL,
                dat = data.bootstrap[[1]][which(data.bootstrap[[1]]$id == i),],
                keys = data.bootstrap[[2]],
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                upper = constraints_up, 
                control = cont
              )
            }
          })
        
        # Check for non-optimized models and remove from output
        
        fail_bootstrap <- failure.check.f(model_results_bootstrap, sim)
        tmp_id_bootstrap <- 1:n_bootstrap
        
        if (length(fail_bootstrap[[1]]) > 0) {
          model_results_bootstrap <- model_results_bootstrap  %>% do(model_results_bootstrap[-fail_bootstrap[[1]], 1:6])
          tmp_id_bootstrap <- tmp_id_bootstrap[-fail_bootstrap[[1]]]
          if (length(tmp_id_bootstrap) == 0) {
            cat("\n All re-samples for for bootstrap", q, "failed to optimize \n" )
            tmp_theta_bootstrap <-
              rbind(
                tmp_theta_bootstrap,
                cbind(
                  t(rep(NA, 4)),
                  rep = q
                )
              )
            next
          }
        }  
        
        # Matrix 'par_mat' contains estimated parameter (columns) for each re-sample (rows, indexed by 'id')
        
        par_mat <-
          cbind(
            matrix(
              unlist(model_results_bootstrap$par),
              ncol = sum(n_parameters),
              byrow = TRUE
            ),
            id = tmp_id_bootstrap
          )
        
        # For each re-sample, estimate mean classification probabilities and associated SEs from estimated parameters and observed covariate values

        tmp_theta_bootstrap <-
          rbind(
            tmp_theta_bootstrap,
            cbind(
              t(sapply(tmp_id_bootstrap, function(x)
                theta.est.f(par_mat[which(par_mat[, "id"] == x),],
                          data.bootstrap[[1]][which(data.bootstrap[[1]]$id == x),]
                          ))),
              rep = q
            )
          )
        
        cat("\n Completed bootstrap", q, "\n")
      } # End of bootstrap estimation loop 
    
    # Remove any 'NA' results for bootstrap re-samples
    tmp_theta_bootstrap <- tmp_theta_bootstrap[(which(!is.na(tmp_theta_bootstrap[, 1]))), ]
      
    ## Estimate standard errors, 95% confidence limits, and 95% confidence interval coverage for classification probabilities (theta) from bootstrapped estimates. Confidence limits are estimated using the bootstrap estimated SE and assuming the sampling distribution of logit scale parameters follows a normal distribution.
    tmp.theta <- betas[which(betas[, dim(betas)[2]] == sim), -grep("^b[0123456789]_", colnames(betas))]
    output_bootstrap <- theta.bootstrap.f(sim_data, tmp.theta, tmp_theta_bootstrap, output_bootstrap)
    } # End of bootstrap loop

    # Complete simulation loop for the current distinct model 
    reps <- sim_profiles$reps[1]
    cat("\n Completed Simulation ", sim, "\n")
    
    # For long simulation runs, prevent data loss by periodically writing results to a temporary file on hard drive 
    if (unclass(Sys.time()) - t_loop > 900) {
      try(tmp <- bind_cols(sim_profiles[1:sim, ], output_global, output_betas))
      try(
        write.csv(tmp, file = paste0("tmp_", out_filename))
      )
      t_loop <- unclass(Sys.time())
    }
    
    # With random seed assigned, append current simulation data to list
    if (any(profiles_names == "seed"))
      sim_data_list <- c(sim_data_list, list(sim_data))

  } # End of simulation loop
  
  # ----- Output & Reports ---
  
  # Combine simulation inputs and outputs, and write 'output_global' as .csv format file to user-specified path
  output_global <- bind_cols(sim_profiles, output_global, output_betas, output_bootstrap, theta_diff)
  try(
    write.csv(output_global, file = out_filename)
  )
  
  #  Optional export of estimated regression coefficients for each simulation replicate
  # try(
  #   write.csv(betas, file = paste0(in_filename, "_out.betas.csv"))
  # )
    
  # With random seed assigned, export simulation data
  if (any(profiles_names == "seed"))
    if (is.integer(sim_profiles$seed[1]) & sim_profiles$seed[1] > 0)
      try(saveRDS(sim_data_list, file = paste0(in_filename, "_output", "_sim_data.Rdata"))
      )
  
  # Print reports on simulation speed and errors/warnings to the console
  elapsed <- unclass(Sys.time()) - t_start
  cat(round(elapsed / 60, 2),
      "m ",
      round(elapsed, 0),
      "s ",
      round(elapsed / (reps * dim(sim_profiles)[1]), 2),
      "s/rep \n")
  print.error.f(error_msg, converge_fail)
  
  # End of model M.phi analysis loop
} else if (model == "M.psi" | model == "M.theta.psi" | model == "M.theta+psi") {
  # ----- MOM/SOM models with covariates predicting true species probabilities (psi) -----
  # true states B = (2 or 3), classification states A = (2 or 3)
  
  ## Define global objects (constant across all simulations)
  
  # For model optimization function 'optim', define user-specified control parameters and lower/upper constraints for parameter values 
  cont <- list(factr = 1e+07, maxit = 100, lmm = 12) # List of control parameters for 'optim' function
  constraints_logit <- qlogis(c(0.0001, 0.9999)) # Constraints for psi/theta parameters, transformed to logit scale
  constraints_g <- c(1.002, Inf) # Constraints for estimated mean group size
  constraints_mix <- c(qlogis(0.0001), Inf)  # Constraints for heterogeneous group parameters
  constraints_b0 <- c(-10, 10); constraints_b1 <- c(-6, 6) # Constraints for regression coefficients (b0 = intercept, b1 = slope) predicting true species probabilities (psi)
  
  # Define objects which will have output appended after each simulation replicate
  betas <- NULL
  output_betas <- NULL
  output_bootstrap <- NULL
  psi_diff <- NULL
  
  # Generate initial objects 
  t_start <- t_loop <- unclass(Sys.time()) # For recording simulation duration
  error_msg <- converge_fail <- NULL # For error messages and convergence warnings
  
  # Build empty data frame 'output_global' containing columns for global estimated parameters (across all simulations) and associated output 
  output_global <- matrix(numeric(1), 0, length(parameters_names_global) * 7 + 1) 
  output_names_global <-
    c(
      paste0("m.", parameters_names_global),
      paste0("sd.", parameters_names_global),
      paste0("se.", parameters_names_global),
      paste0("cv.", parameters_names_global),
      paste0("me.", parameters_names_global),
      paste0("rm.", parameters_names_global),
      paste0("ci.", parameters_names_global),
      "rep.tru"
    )
  colnames(output_global) <- output_names_global
  output_global <- data.frame(output_global)

  # Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
  for (sim in 1:dim(sim_profiles)[1]) { 
    
    ## Define objects for simulation of the current distinct model
    
    O_ps <- c(sim_profiles$O_p[sim], sim_profiles$O_s[sim]) # Count of primary and secondary observers
    
    # Vector 'n_parameters' contains the number of estimated parameters for the current simulation in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), un-transformed parameters (4), and heterogeneous group probabilities (pi) parameters (5)
    # Vectors 'parameters_names' and 'parameters_col' contain names and column #s of true parameter values in 'sim_profiles'
    
    parameters_col <- grep("b0|b1|psi|theta", profiles_names) # Column #s of parameters for estimating psi & theta
    parameters_names <- c( grep("b0|b1|psi|theta", parameters_names_global, value = TRUE))
    
    # Add regression coefficients (b0 = intercept, b1 = slope) and logit parameters for classification probabilities (theta) and true species probabilities (psi)
    n_parameters <-
      c(length(grep("b0", profiles_names)), 
        length(grep("b1", profiles_names)), 
        length(grep("theta|psi", parameters_names)) - length(grep("b0|b1", profiles_names)), 
        0, 0)
    
    # Add parameters for mean group size g
    if (any(sim_profiles[sim, grep("^g_[0123456789]", profiles_names)] > 1)) {
      parameters_names <- c(parameters_names, grep("^g_[0123456789]", profiles_names, value = TRUE))
      parameters_col <- c(parameters_col, grep("^g_[0123456789]", profiles_names))
      n_parameters[4] <- length(parameters_names) - sum(n_parameters[1:3])
    }
    
    # Add parameters for heterogeneous group probabilities: heterogeneous group probabilities (pi) or heterogeneous group affinities (rho)
    if (any(sim_profiles[sim, mix_col] > 0)) {
      parameters_names <- c(parameters_names, profiles_names[mix_col])
      parameters_col <- c(parameters_col, grep("mix", profiles_names))
      n_parameters[5] <- length(mix_col)
    }
    
    # With certain ID by secondary observers (probability of uncertain ID = 0), remove these classification probabilities from estimated parameters
    if (length(sim_profiles[sim, grep("theta_s", profiles_names)]) > 0) {
      if (any(sim_profiles[sim, grep("theta_s", profiles_names)] == 0)) {
        id_certain_col <- grep("theta_s", parameters_names)
        parameters_names <- parameters_names[-id_certain_col]
        parameters_col <- parameters_col[-id_certain_col]
        n_parameters[3] <- n_parameters[3] - length(id_certain_col)
        # Make secondary observers = 1 (in case input file is incorrect)
        sim_profiles$O_s[sim] <- O_ps[2] <- 1
      } 
    }
    
    # Generate initial objects 
    parameters_ini <- generate.ini.f(sim_profiles[sim, ], n_parameters) # Sets of alternative initial parameter values on each row
    output_names <-
      c(
        paste0("m.", parameters_names),
        paste0("sd.", parameters_names),
        paste0("se.", parameters_names),
        paste0("cv.", parameters_names),
        paste0("me.", parameters_names),
        paste0("rm.", parameters_names),
        paste0("ci.", parameters_names),
        "rep.tru"
      )
    
    # Vector 'constraints_low' contains lower constraints for each parameter
    constraints_low <- c(rep(constraints_b0[1], (B - 1)),  rep(constraints_b1[1], (B - 1)),  
                         rep(constraints_b0[1], (n_parameters[1] - (B - 1))),  
                         rep(constraints_b1[1], (n_parameters[2] - (B - 1))),  
                         rep(constraints_logit[1], n_parameters[3]), 
                         rep(constraints_g[1], n_parameters[4]), 
                         rep(constraints_mix[1], n_parameters[5]))
    
    # Optional upper constraints default to infinity (no constraints)
    constraints_up <- Inf
    
    # List 'sim_data' contains formatted simulated survey observations and (optionally) a keyed table of unique groups records and a keyed table of unique values for a covariate predicting true species probabilities (psi) 
    sim_data <- generate.simulation.data.f(sim_profiles[sim,]) 
    parameters_index <- 1 # Index for alternative sets of initial parameter values in 'parameters_ini'
    
    # For optional output of regression coefficients
    sim_data_tmp <- sim_data[[1]]
      
    # If present, apply user-specified upper constraints for classification probabilities (theta) exist in 'theta_limit'
    if (exists("theta_limit") & sim_profiles[1, ]$Model == "M.psi") {
      if (eval(theta.limit.f(0.1)) > 0) {
        constraints_up <- c(
          rep(constraints_b0[2], (B - 1)),
          rep(constraints_b1[2], (B - 1)),
          rep(constraints_b0[2], (n_parameters[1] - (B - 1))),
          rep(constraints_b1[2], (n_parameters[2] - (B - 1))),
          rep(constraints_logit[2], n_parameters[3]),
          rep(constraints_g[2], n_parameters[4]),
          rep(constraints_mix[2], n_parameters[5])
        )
        theta_col <- c(grep("theta", profiles_names))
        constraints_up[(sum(n_parameters[1:2]) + 1):(sum(n_parameters[1:2]) + length(theta_col))] <- 
          qlogis(
            unlist(lapply(sim_profiles[sim, theta_col], function(x)
              theta.limit.f(x)))
          )
      }
    }
    
    ## ----- Model Optimization ----- 
    
    # Function for simulation progress bar
    id_max <- max(sim_data[[1]]$id)
    pb <- txtProgressBar(max = id_max, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Fit models via function 'optimize.M.psi.f', ''optimize.M.theta.psi.f' utilizing parallel execution of the 'optim' routine using the "L-BFGS-B" method, with statistical output to list 'model_results'
    
    model_results <-
      as.data.frame(foreach(
        i = 1:id_max, 
        .combine = rbind, 
        .packages = c("plyr", "dplyr"), 
        .errorhandling = "pass", 
        .options.snow = opts
      ) %dopar% {
        if (sim_profiles$Model[1] == "M.psi") {
          # Optimization for model with covariate predicting true species probabilities (psi)
          optim(
            parameters_ini[parameters_index, ], 
            optimize.M.psi.f, 
            gr = NULL,
            dat = sim_data[[1]][which(sim_data[[1]]$id == i), ],
            keys = sim_data[[2]],
            keys_psi = sim_data[[3]],
            sim_profile = sim_profiles[sim, ],
            hessian = TRUE, 
            method = c("L-BFGS-B"), 
            lower = constraints_low, 
            upper = constraints_up, 
            control = cont
          )
        }else if (sim_profiles$Model[1] == "M.theta.psi" | sim_profiles$Model[1] == "M.theta+psi") {
          # Optimization for model with covariates predicting classification probabilities (theta) and true species probabilities (psi)
            optim(
              parameters_ini[parameters_index, ], 
              optimize.M.theta.psi.f, 
              gr = NULL,
              dat = sim_data[[1]][which(sim_data[[1]]$id == i), ],
              sim_profile = sim_profiles[sim, ],
              hessian = TRUE, 
              method = c("L-BFGS-B"), 
              lower = constraints_low, 
              control = cont
            )
        }
      }) 
      
    ## ----- Model Output: Error-checking -------  
    
    # Check for model optimization errors/convergence failures
    fail <- failure.check.f(model_results, sim)
    
    # Attempt to fit non-optimized models using alternate initial parameter values
    if (length(fail[[1]]) > 0) {
      error_msg <- bind_rows(error_msg, fail[[2]])
      fail_tmp <- fail <- fail[[1]]
      
      # Try alternate sets of initial parameter values until model(s) optimize successfully
      while (parameters_index < 5 && length(fail) > 0) {
        parameters_index <- parameters_index + 1 # Increment index for alternate sets of initial values
        # >1 models necessary to use 'foreach' function
        if (length(fail) == 1) {
          fail_tmp <- c(fail, fail) 
        }
        
        # Attempt optimization with alternate initial parameter values
        model_results_refit <-
          as.data.frame(foreach(
            i = 1:length(fail_tmp), 
            .combine = rbind, 
            .packages = c("plyr", "dplyr"), 
            .errorhandling = "pass"
          ) %dopar% {
            if (sim_profiles$Model[1] == "M.psi") {
              optim(
                parameters_ini[parameters_index ,], 
                optimize.M.psi.f, 
                gr = NULL, 
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]),], 
                keys = sim_data[[2]],
                keys_psi = sim_data[[3]],
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                upper = constraints_up, 
                control = cont
              )
            }else if (sim_profiles$Model[1] == "M.theta.psi" | model == "M.theta+psi") {
              optim(
                parameters_ini[parameters_index ,], 
                optimize.M.theta.psi.f, 
                gr = NULL, 
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]),], 
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                control = cont
              )
            }
          })
        
        # Add model(s) successfully optimizing to 'model_results'
        model_results_refit <- refit.check.f(model_results_refit)
        success <- which(model_results_refit[1:length(fail), "convergence"] == 0)
        if (length(success) > 0) {
          l_ply(success, function(x) {
            for (i in 1:6)
              model_results[[i]][fail[x]] <<- model_results_refit[[i]][x]
          })
          fail <- fail[-success] # Remove successful models from failed models
        }
      }
      
      # If any model(s) never optimize adequately: generate warnings, remove failed model(s) from 'model_results', and update # of (successful) simulation replicates
      if (length(fail) > 0) {
        model_results <- model_results  %>% do(model_results[-fail, 1:6])
        converge_fail <- rbind(converge_fail, as_tibble(list(
          sim = rep(sim, length(fail)), 
          id = fail, 
          code = model_results_refit$convergence[which(model_results_refit$convergence > 0)][1:length(fail)], 
          msg = model_results_refit$message[which(model_results_refit$convergence > 0)][1:length(fail)]
        )))
        
        # For optional output for regression coefficients, remove data for non-optimized models
        for (i in fail) {
          sim_data_tmp <- filter(sim_data_tmp, id != i)
        }
        reps <- reps - length(fail)
      }  
    } # End of model re-fitting

    ## ----- Model results summary and output -----
    
    ## Back-transform estimates for logit parameters to probability scale and calculate summary statistics. Prefixes on variable names indicate summary statistics: m =  mean, sd = standard deviation, cv = coefficient of variation (sd / mean), me = mean error, rm = residual mean squared error, ci = 95% confidence interval coverage
    
    # Vector 'parameters_true' contains true parameter values
    parameters_true <- sim_profiles[sim, parameters_col]
    
    # Append summary statistics for this distinct model (vector 'output') as new row in data frame 'output_global'
    parameters_output <- model.results.f(model_results, parameters_true, reps, n_parameters)
    output <- summarise.results.f(parameters_output, parameters_true, reps)
    names(output) <- output_names
    output_global <- bind_rows(output_global, data.frame(t(output)))
    
    # Optional output summarizing estimated versus true group-level predictions for true species probabilities. Un-comment the following 'if' clause to append results to data frame 'output'.
    # if (model == "M.psi") {
    #   tmp <- mlogit.group.predict.f(sim_data[[1]], sim_profiles[sim, ], output)
    #   psi_diff <- bind_rows(psi_diff, as_tibble(t(as.matrix(tmp))))
    # }
    
    ## Output related to estimates of overall mean true species probabilities (psi) and mean classification probabilities (theta) derived from regression coefficients of multinomial logistic regression and observed covariate values. Alternative 95% confidence limits for estimated means were constructed assuming 1) asymptotic normal distribution in natural scale (i.e, probabilities) and 2) asymptotic normal distribution of parameters in logit scale.
    # For each distinct model (rows), 'output_betas' contains summarized statistics (mean error, root mean square error, and 95% CI coverage) for overall mean estimates of true species probabilities (psi) and classification probabilities
    # For each simulation replicate (rows), 'betas' contains estimated regression coefficients, estimated overall mean probabilities with associated standard errors, 95% confidence intervals ('lg' denotes limits estimated in logit scale), and 95% confidence interval coverage
    
    tmp <- psi.derived.f(sim_profiles, sim_data, sim_data_tmp, model_results, parameters_output, output_betas)
    output_betas <- tmp[[1]]
    betas <- tmp[[2]]
    
    ## Loop for bootstrap replications
    # For each bootstrap re-sample, data are re-sampled with replacement from simulated survey data, models are fit to re-sampled data, and estimated parameters from re-sampled data are used to estimate mean values for classification probabilities (theta) across observed covariate values
    
    # Bootstrap only if more than 1 bootstrap re-samples (n_bootstrap > 1) are specified in simulation inputs  
    if (n_bootstrap > 1 & B < 4) {
      if (sim_profiles$Model[1] == "M.theta.psi" & A > 2 | sim_profiles$Model[1] == "M.theta+psi" & A > 2) {
        cat("For models 'M.theta.psi' and 'M.theta+psi', bootstrap estimates only supported for A = 2 \n")
      }else{
      remove_col <- grep("id|count|key", colnames(sim_data[[1]]))
      tmp_psi_bootstrap <- NULL
      tmp_theta_bootstrap <- NULL
      tmp_id <- unique(sim_data_tmp$id)
      parameters_index <- 1
      cat("\n Completed model optimization for simulation", sim, "\n")

      # Loop for bootstrap re-samples for each replicate (indexed by 'q') of simulated data
      for (q in 1:reps) {

        # Simulated survey data for current simulation replicate 
        data.rep <- filter(sim_data_tmp, id == tmp_id[q])
        
        # Sample groups (with replacement) from simulated data to generate n (= n_bootstrap) re-samples of the simulated survey data
        tmp <-
          ldply(1:n_bootstrap, function(x)
            sample_n(
              data.rep,
              size = sum(data.rep$count),
              replace = TRUE,
              weight = count
            )) %>%
          select(., !any_of(remove_col)) %>%
          bind_cols(., id = rep(1:n_bootstrap, each = sum(data.rep$count)))
        
        # List 'data.bootstrap' contains re-sampled data for the current simulation replicate, with formatting is identical to 'sim_data'
        data.bootstrap <-
          format.MOM.data.f(tmp, A, O_ps, sim_profiles[sim, mix_col], sim_profiles[sim, ]$n_bins)
        
        # Fit models to re-sampled data using the appropriate function and parallel optimization of the 'optim' routine
        pb <- txtProgressBar(max = n_bootstrap, style = 3)
        
        model_results_bootstrap <-
          as.data.frame(foreach(
            i = 1:n_bootstrap,
            .combine = rbind,
            .packages = c("plyr", "dplyr"),
            .errorhandling = "pass", 
            .options.snow = opts
          ) %dopar% {
            if (sim_profiles$Model[1] == "M.psi") {
              # Optimization for model with covariate predicting true species probabilities (psi)
              optim(
                parameters_ini[parameters_index,],
                optimize.M.psi.f,
                gr = NULL,
                dat = data.bootstrap[[1]][which(data.bootstrap[[1]]$id == i), ],
                keys = data.bootstrap[[2]],
                keys_psi = data.bootstrap[[3]],
                sim_profile = sim_profiles[sim,],
                hessian = TRUE,
                method = c("L-BFGS-B"),
                lower = constraints_low,
                upper = constraints_up,
                control = cont
              )
            }else if (sim_profiles$Model[1] == "M.theta.psi" | sim_profiles$Model[1] == "M.theta+psi") {
              # Optimization for model with separate covariates predicting classification probabilities (theta) and true species probabilities (psi)
              optim(
                parameters_ini[parameters_index, ], 
                optimize.M.theta.psi.f, 
                gr = NULL,
                dat = data.bootstrap[[1]][which(data.bootstrap[[1]]$id == i), ],
                sim_profile = sim_profiles[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                control = cont
              )
            }
          })
        
        # Check for non-optimized models and remove from output
        fail_bootstrap <- failure.check.f(model_results_bootstrap, sim)
        tmp_id_bootstrap <- 1:n_bootstrap

        if (length(fail_bootstrap[[1]]) > 0) {
          model_results_bootstrap <- model_results_bootstrap  %>% do(model_results_bootstrap[-fail_bootstrap[[1]], 1:6])
          tmp_id_bootstrap <- tmp_id_bootstrap[-fail_bootstrap[[1]]]
          if (length(tmp_id_bootstrap) == 0) {
            cat("\n All re-samples for for bootstrap", q, "failed to optimize \n" )
            
            tmp_psi_bootstrap <-
              rbind(
                tmp_psi_bootstrap,
                cbind(
                  t(rep(NA, B)),
                  rep = q
                )
              )
            
            if (sim_profiles[1, ]$Model == "M.theta.psi") {
              tmp_theta_bootstrap <-
                rbind(
                  tmp_theta_bootstrap,
                  cbind(
                    t(rep(NA, B)),
                    rep = q
                  )
                )
            }
            next
          }
        }  
        
        # Matrix of parameter estimates (columns) for each re-sample (rows, indexed by 'id')
        
        par_mat <-
          cbind(
            matrix(
              unlist(model_results_bootstrap$par),
              ncol = sum(n_parameters),
              byrow = TRUE
            ),
            id = tmp_id_bootstrap
          )
          
        # For each re-sample, estimate mean true species probabilities and associated SEs from estimated parameters and observed covariate values
        
        tmp_psi_bootstrap <-
          rbind(
            tmp_psi_bootstrap,
            cbind(
              t(sapply(tmp_id_bootstrap, function(x)
                psi.est.f(par_mat[which(par_mat[, "id"] == x),][1:sum(n_parameters)],
                          data.bootstrap[[1]][which(data.bootstrap[[1]]$id == x),],
                          n_parameters))),
              rep = q
            )
          )
        
        if (sim_profiles[1, ]$Model == "M.theta.psi" | sim_profiles$Model[1] == "M.theta+psi") {
          
          # For each re-sample, estimate mean classification probabilities and associated SEs from estimated parameters and observed covariate values
          
          tmp_theta_bootstrap <-
            rbind(
              tmp_theta_bootstrap,
              cbind(
                t(sapply(tmp_id_bootstrap, function(x)
                  theta.est.f(par_mat[which(par_mat[, "id"] == x),],
                              data.bootstrap[[1]][which(data.bootstrap[[1]]$id == x),]
                  ))),
                rep = q
              )
            )
        }

        cat("\n Completed bootstrap", q, "\n")
      }
      }# End of bootstrap q loop 
      
      # Remove any 'NA' results for bootstrap re-samples
      tmp_psi_bootstrap <- tmp_psi_bootstrap[(which(!is.na(tmp_psi_bootstrap[, 1]))), ]
      tmp_theta_bootstrap <- tmp_theta_bootstrap[(which(!is.na(tmp_theta_bootstrap[, 1]))), ]
        
      ## Estimate standard errors, 95% confidence limits, and 95% confidence interval coverage for true species probabilities (psi) from bootstrapped estimates. Confidence limits are estimated using the bootstrap estimated SE and assuming the sampling distribution of logit scale parameters follows a normal distribution.
      tmp.psi <- betas[which(betas[, dim(betas)[2]] == sim)  ,  -grep("^b[0123456789]_", colnames(betas))]
      
      output_bootstrap <- psi.bootstrap.f(sim_data, tmp.psi, tmp_psi_bootstrap, tmp_theta_bootstrap, output_bootstrap)
      
    } # End of bootstrap code
    
    # Complete simulation loop for the current distinct model
    reps <- sim_profiles$reps[1]
    cat("\n Completed Simulation ", sim, "\n")
    
    # For long simulation runs, prevent data loss by periodically writing results to a temporary file on hard drive  
    if (unclass(Sys.time()) - t_loop > 900) {
      try(tmp <- bind_cols(sim_profiles[1:sim, ], output_global, output_betas))
      try(
        write.csv(tmp, file = paste0("tmp_", out_filename))
      )
      t_loop <- unclass(Sys.time())
    }
    
    # With random seed assigned, append current simulation data to list
    if (any(profiles_names == "seed"))
      sim_data_list <- c(sim_data_list, list(sim_data))

  } # End of loop for each profile simulation 
  
  # ----- Output & Reports ---
  
  # Combine simulation inputs and outputs, and write 'output_global' as .csv format file to user-specified path
  output_global <- bind_cols(sim_profiles, output_global, output_betas, output_bootstrap, psi_diff)
  try(
    write.csv(output_global, file = out_filename)
  )
  
  #  Optional export of estimated regression coefficients for each simulation replicate
  # try(
  #   write.csv(betas, file = paste0(in_filename, "_out.betas.csv"))
  # )

  # With random seed assigned, export simulation data
  if (any(profiles_names == "seed"))
    if (is.integer(sim_profiles$seed[1]) & sim_profiles$seed[1] > 0)
      try(saveRDS(sim_data_list, file = paste0(in_filename, "_output", "_sim_data.Rdata")))
  
  # Print reports on simulation speed and errors/warnings to the console
  elapsed <- unclass(Sys.time()) - t_start
  cat(round(elapsed / 60, 2),
      "m ",
      round(elapsed, 0),
      "s ",
      round(elapsed / (reps * dim(sim_profiles)[1]), 2),
      "s/rep \n")
  print.error.f(error_msg, converge_fail)
  
  # End of model M.psi and M.psi.theta analysis loop
} 
        
###############################################################################
#    MOM/SOM model simulations: Un-modeled heterogeneity (covariates)
###############################################################################

# These routines conduct "covariate" simulations for MOM/SOM models, assuming un-modeled heterogeneity in the data arising from covariates predicting true species probabilities (psi) or classification probabilities (theta). Data are generated with true covariate structure, but are analyzed without covariates. 
# Conduct analyses by selecting and executing all code in this section. 

# Output files include true parameters used to generate data, a column named "heterogeneity" taking value 1, and summary statistics for estimated parameters. For parameters with un-modeled heterogeneity, estimates of mean error and confidence interval coverage are found by comparing estimates to overall means predicted from true parameter values across all covariate values. 

if (model == "M.psi" | model == "M.theta" | model == "M.theta.p" | model == "M.theta.ps") {
  # True states B = (2 or 3), classification states A = (2 or 3)
  
  ## Define global objects (constant across all simulations)
  
  # Generate initial objects 
  t_start <- t_loop <- unclass(Sys.time()) # For recording simulation duration
  error_msg <- converge_fail <- NULL # For error messages and convergence warnings
  
  # For model optimization function 'optim', define user-specified control parameters and lower/upper constraints for parameter values
  cont <- list(factr = 1e+07, maxit = 200, lmm = 12) # List of control parameters for 'optim' function
  constraints_logit <- qlogis(c(0.0001, 0.9999)) # Constraints for psi/theta parameters, transformed to logit scale
  constraints_g <- c(1.002, Inf) # Constraints for estimated mean group size
  constraints_mix <- c(qlogis(0.0001), Inf)  # Constraints for group mixing parameters
  
  # Field 'heterogeneity' = 1 denotes un-modeled heterogeneity arising from covariates
  sim_profiles <- mutate(sim_profiles, heterogeneity = 1)
  
  # Suffix "_het' denotes objects used for analyses and output (i.e. with un-modeled heterogeneity)
  sim_profiles_het <- sim_profiles[, c(1:max(mix_col), grep("^g_[0123456789]", profiles_names))]
  sim_profiles_het$Model <- "M"
  sim_profiles_het$mx_model <- "constant"
  group_size_col <- grep("^g_[0123456789]", profiles_names) # Column #s of parameters for estimating group sizes
  
  # Add names of new parameters for heterogeneity analyses 
  if (B == 2) {
    if (A == 2) {
        sim_profiles_het <- 
          mutate(sim_profiles_het, 
                 theta_p_12 = 1, theta_p_21 = 1, theta_s1_12 = 1, theta_s1_21 = 1, psi_1 = 1)
    }else if (A == 3) {
        sim_profiles_het <- 
          mutate(sim_profiles_het,
                 theta_p_21 = 1, theta_p_31 = 1, theta_p_12 = 1, theta_p_13 = 1,  
                 theta_s1_21 = 1, theta_s1_31 = 1, theta_s1_12 = 1, theta_s1_13 = 1, psi_1 = 1)
    }
  }else if (B == 3) {
      sim_profiles_het <- 
        mutate(sim_profiles_het,
               theta_p_21 = 1, theta_p_31 = 1, theta_p_12 = 1, theta_p_32 = 1, theta_p_13 = 1, theta_p_23 = 1,
               theta_s1_21 = 1, theta_s1_31 = 1, theta_s1_12 = 1, theta_s1_32 = 1, theta_s1_13 = 1, theta_s1_23 = 1,
               psi_1 = 1, psi_2 = 1)
  }
  
  sim_profiles_het_names <- names(sim_profiles_het)
  if (any(sim_profiles_het$O_p == 0)) {
    sim_profiles_het <- sim_profiles_het[, -grep("theta_p", sim_profiles_het_names)]
    sim_profiles_het_names <- sim_profiles_het_names[-grep("theta_p", sim_profiles_het_names)]
  }
  
  # Vector 'parameter_names_het_global' contains parameters names across all simulations 
  parameters_names_het_global <- c(grep("theta|psi", sim_profiles_het_names, value = TRUE)) 
  
  # Add parameters for group size (g) and heterogeneous group parameters (pi or rho) 
  if (any(sim_profiles_het[grep("^g_[0123456789]", sim_profiles_het_names)] > 1)) {
    parameters_names_het_global <- c(parameters_names_het_global, grep("^g_[0123456789]", sim_profiles_het_names, value = TRUE))}
  if (any(sim_profiles_het[, mix_col] > 0)) {
    parameters_names_het_global <- c(parameters_names_het_global, sim_profiles_het_names[mix_col])}
  
  # Build empty data frame 'output_global' containing columns for global estimated parameters (across all simulations) and associated output 
  output_global <- matrix(numeric(1), 0, length(parameters_names_het_global) * 7 + 1) 
  output_names_global <-
    c(
      paste0("m.", parameters_names_het_global),
      paste0("sd.", parameters_names_het_global),
      paste0("se.", parameters_names_het_global),
      paste0("cv.", parameters_names_het_global),
      paste0("me.", parameters_names_het_global),
      paste0("rm.", parameters_names_het_global),
      paste0("ci.", parameters_names_het_global),
      "rep.tru"
    )
  colnames(output_global) <- output_names_global
  output_global <- data.frame(output_global)
  
  # Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
  for (sim in 1:dim(sim_profiles)[1]) { 
    
    ## Define objects for simulation of the current distinct model
    
    O_ps <- c(sim_profiles$O_p[sim], sim_profiles$O_s[sim]) # Count of primary and secondary observers
    
    # Vector 'parameters_names' contains names of estimated parameters for the current simulation  
    parameters_names <- c(grep("theta|psi", parameters_names_het_global, value = TRUE))

    # Vector 'n_parameters_het' contains the number of estimated parameters for the current simulation in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), un-transformed parameters (4), and heterogeneous group probabilities (pi) parameters (5)
    n_parameters_het <- c(0, 0, length(grep("theta|psi", parameters_names)), 0, 0)
    
    # Add parameters for mean group size g
    if (any(sim_profiles[sim, group_size_col] > 1)) {
      parameters_names <- c(parameters_names, grep("^g_[0123456789]", sim_profiles_het_names, value = TRUE))
      n_parameters_het[4] <- length(parameters_names) - sum(n_parameters_het[3])
    }
    
    # Add parameters for heterogeneous group probabilities: heterogeneous group probabilities (pi) or heterogeneous group affinities (rho)
    if (any(sim_profiles_het[sim, mix_col] > 0)) {
      parameters_names <- c(parameters_names, sim_profiles_het_names[mix_col])
      n_parameters_het[5] <- length(mix_col)
    }
    n_parameters <- n_parameters_het
    
    output_names <-
      c(
        paste0("m.", parameters_names),
        paste0("sd.", parameters_names),
        paste0("se.", parameters_names),
        paste0("cv.", parameters_names),
        paste0("me.", parameters_names),
        paste0("rm.", parameters_names),
        paste0("ci.", parameters_names),
        "rep.tru"
      )
    
    # Vector 'constraints_low' contains lower boundaries for all parameters
    constraints_low <- c(rep(constraints_logit[1], n_parameters_het[3]), 
                         rep(constraints_g[1], n_parameters_het[4]), 
                         rep(constraints_mix[1], n_parameters_het[5]))
    
    # Upper constraints default to infinity (no constraints) 
    constraints_up <- Inf
    
    parameters_ini <- generate.ini.f(sim_profiles_het[sim, ], n_parameters_het) # Sets of alternative initial parameter values on each row
    
    # List 'sim_data' contains formatted simulated survey observations generated from the true model (i.e. with predictive covariates) and (optionally) a keyed table of unique groups records
    sim_data <- generate.simulation.data.f(sim_profiles[sim,]) 
    parameters_index <- 1 # Index for set of values used as initial parameters for optimization
    
    # Vector 'parameters_true' contains true parameters values from the data-generating model for the current simulation
    # Matrix 'psi_betas_mat' gives beta intercepts (column 1) and slopes (column 2) for regression parameters for predicting true species probabilities (psi)

    if (sim_profiles[sim, ]$Model == "M.psi") {
      parameters_true <- c(
        unlist(sim_profiles[sim, c(grep("theta", profiles_names))]))
      
      psi_betas_mat <- matrix(unlist(sim_profiles[sim, c(grep("_psi", profiles_names))]), ncol = 2, byrow = TRUE) 
      
      if (B == 2) {
        parameters_true <- c(parameters_true,
                      sim_data[[4]][[1]][1])
      }else{
        psi_betas_mat <- matrix(unlist(sim_profiles[sim, c(grep("_psi", profiles_names))]), ncol = 2, byrow = TRUE)
        parameters_true <- c(parameters_true,
                      sim_data[[4]][[1]][1:2])
      }
    }
    
    # Matrix 'theta_betas_mat' gives beta intercepts (column 1) and slopes (column 2) for regression coefficients predicting classification probabilities (theta)
    
    # if (sim_profiles[sim, ]$Model == "M.theta") {
    #   theta_betas_mat <- matrix(c(
    #     unlist(sim_profiles[sim, c(grep("b0|b1", profiles_names))])),
    #     ncol = 2)
    #   parameters_true <- c(
    #     apply(theta_betas_mat, 1, logistic.dist.f)[1, ],
    #     unlist(sim_profiles[sim, c(grep("psi_", profiles_names))])
    #   )
    # }
    
    if (sim_profiles[sim, ]$Model == "M.theta" | sim_profiles[sim, ]$Model == "M.theta.p" | sim_profiles[sim, ]$Model == "M.theta.ps") {
      parameters_true <- c(
        as.vector(sim_data[[4]][[2]]), 
        unlist(sim_profiles[sim, c(grep("^theta_", profiles_names))]),
        unlist(sim_profiles[sim, c(grep("psi_", profiles_names))])
      )
    }
    
    if (any(sim_profiles[sim, group_size_col] > 1)) {
      parameters_true <- c(
        parameters_true, unlist(sim_profiles[sim,group_size_col])
      )
    }
    
    if (any(sim_profiles[sim, mix_col] > 0)) {
      if (sim_profiles$mx_model[1] == "encounter") {
        tmp_mix <- unlist(sim_profiles[sim, mix_col])
        if (any(sim_profiles[sim, group_size_col] > 1)) {
          g <- unlist(sim_profiles[sim, group_size_col])
        }
        
        dist_psi_tru <- mlogit.regress.predict.f(rnorm(10^6), psi_betas_mat, B)
        dist_psi_tru <- t(t(dist_psi_tru) / g)
        dist_psi_tru <- dist_psi_tru / rowSums(dist_psi_tru)
        
        # In vector 'mix_abs', compute mean values of heterogeneous group probability parameters across predicted values of true species probabilities (psi) and add to 'parameters_true'
        if (B == 2) {
          mix_abs <- aaply(dist_psi_tru, 1, function(x)
            tmp_mix * min(x) ^ 2 * max(x) / min(x))
          mix_abs <- matrix(mix_abs, ncol = 1)
          colnames(mix_abs) <- "pi.12"
        } else if (B == 3) {
          pairs <- matrix(c(
            1, 2, 1, 3, 2, 3
          ), ncol = 2, byrow = TRUE)
          mix_abs <- vapply(1:3, function(x)
            tmp_mix[x] *
              (pmin(dist_psi_tru[, pairs[x, 1]], dist_psi_tru[, pairs[x, 2]])) ^ 2 *
              pmax(dist_psi_tru[, pairs[x, 1]], dist_psi_tru[, pairs[x, 2]]) /
              pmin(dist_psi_tru[, pairs[x, 1]], dist_psi_tru[, pairs[x, 2]])
            , numeric(dim(dist_psi_tru)[1]))
          colnames(mix_abs) <- c("pi.12", "pi.13", "pi.23")
        }
        parameters_true <- c(
          parameters_true, colMeans(mix_abs))
      } else {
        parameters_true <- c(
          parameters_true, unlist(sim_profiles[sim, mix_col])
        )
      }
    }
    
    # If present, apply user-specified upper constraints for classification probabilities (theta) exist in 'theta_limit'
    
    if (exists("theta_limit")) {
      if (eval(theta.limit.f(0.1)) > 0) {
        constraints_up <- c(
          rep(constraints_logit[2], n_parameters_het[3]),
          rep(constraints_g[2], n_parameters_het[4]),
          rep(constraints_mix[2], n_parameters_het[5])
        )
        
        theta_col <- c(grep("theta", colnames(sim_profiles_het)))
        
        constraints_up[1:length(theta_col)] <- 
          qlogis(
            unlist(lapply(parameters_true[1:length(theta_col)], function(x)
              theta.limit.f(x)))
          )
      }
    }
    
    # ----- Model Optimization ----- 
    
    # Function for simulation progress bar
    id_max <- max(sim_data[[1]]$id)
    pb <- txtProgressBar(max = id_max, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Fit models via function 'optimize.M.f' utilizing parallel execution of the 'optim' routine using the "L-BFGS-B" method, with statistical output to list 'model_results'
    # Estimation models do not include predictive covariates
    
    model_results <-
      as.data.frame(foreach(
        i = 1:id_max, 
        .combine = rbind, 
        .packages = c("plyr", "dplyr", "extraDistr"), 
        .errorhandling = "pass", 
        .options.snow = opts
      ) %dopar% {
        optim(
          parameters_ini[parameters_index ,], 
          optimize.M.f, gr = NULL, 
          dat = sim_data[[1]][which(sim_data[[1]]$id == i),], 
          keys = sim_data[[2]], 
          sim_profile = sim_profiles_het[sim, ],
          hessian = TRUE, 
          method = c("L-BFGS-B"), 
          lower = constraints_low, 
          upper = constraints_up, 
          control = cont
        )
      })
    
    # ----- Model Output: Error-checking -------  
    
    # Check for model optimization errors/convergence failures
    fail <- failure.check.f(model_results, sim)
    
    # Attempt to fit non-optimized models using alternate initial parameter values
    if (length(fail[[1]]) > 0) {
      error_msg <- bind_rows(error_msg, fail[[2]])
      fail_tmp <- fail <- fail[[1]]
      while (parameters_index < 5 && length(fail) > 0) {
        parameters_index <- parameters_index + 1 # Increment index for alternate sets of initial values
        # >1 models necessary to use 'foreach' function
        if (length(fail) == 1) {
          fail_tmp <- c(fail, fail) 
        }
        
        # Attempt optimization with alternate initial parameter values
        
        model_results_refit <-
          as.data.frame(foreach(
            i = 1:length(fail_tmp), 
            .combine = rbind, 
            .packages = c("plyr", "dplyr", "extraDistr"), 
            .errorhandling = "pass"
          ) %dopar% {
              optim(
                parameters_ini[parameters_index ,], 
                optimize.M.f, gr = NULL, 
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]), ], 
                keys = sim_data[[2]], 
                sim_profile = sim_profiles_het[sim, ],
                hessian = TRUE, 
                method = c("L-BFGS-B"), 
                lower = constraints_low, 
                control = cont
              )
          })
        
        # Add model(s) successfully optimizing to 'model_results'
        model_results_refit <- refit.check.f(model_results_refit)
        success <- which(model_results_refit[1:length(fail), "convergence"] == 0)
        if (length(success) > 0) {
          l_ply(success, function(x) {
            for (i in 1:6)
              model_results[[i]][fail[x]] <<- model_results_refit[[i]][x]
          })
          fail <- fail[-success] # Remove successful models from failed models
        }
      }
      
      # If any model(s) never optimize adequately: generate warnings, remove failed model(s) from 'model_results', and update # of (successful) simulation replicates
      if (length(fail) > 0) {
        model_results <- model_results  %>% do(model_results[-fail, 1:6])
        converge_fail <- rbind(converge_fail, as_tibble(list(
          sim = rep(sim, length(fail)), 
          id = fail, 
          code = model_results_refit$convergence[which(model_results_refit$convergence > 0)][1:length(fail)], 
          msg = model_results_refit$message[which(model_results_refit$convergence > 0)][1:length(fail)]
        )))
        reps <- reps - length(fail)
      }  
    } # End of model refit loop
    
    # Back-transform estimates for logit parameters to probability scale and calculate summary statistics
    # Prefixes on variable names indicate summary statistics: m =  mean, sd = standard deviation, cv = coefficient of variation (sd / mean), me = mean error, rm = residual mean squared error, ci = 95% confidence interval coverage
    
    # Append summary statistics for this distinct model (vector 'output') as new row in data frame 'output_global'
    parameters_output <- model.results.f(model_results, parameters_true, reps, n_parameters_het)
    output <- summarise.results.f(parameters_output, parameters_true, reps)
    names(output) <- output_names
    output_global <- bind_rows(output_global, data.frame(t(output)))
    
    # Complete simulation loop for the current distinct model
    reps <- sim_profiles$reps[1]
    cat("\n Completed Simulation ", sim, "\n")
    
    # For long simulation runs, prevent data loss by periodically writing results to a temporary file on hard drive  
    if (unclass(Sys.time()) - t_loop > 900) {
      try(
        write.csv(output_global, file = paste0("tmp_", out_filename))
      )
      t_loop <- unclass(Sys.time())
    }
  } # End of profile simulation loop
  
  ## ----- Output & Reports -----
  
  # Combine simulation inputs and outputs, and write 'output_global' as .csv format file to user-specified path
  output_global <- bind_cols(sim_profiles, output_global)
  try(
    write.csv(output_global, file = out_filename)
  )
  
  # Print reports on simulation speed and errors/warnings to the console
  elapsed <- unclass(Sys.time()) - t_start
  cat(round(elapsed / 60, 2),
      "m ",
      round(elapsed, 0),
      "s ",
      round(elapsed / (reps * dim(sim_profiles)[1]), 2),
      "s/rep \n")
  print.error.f(error_msg, converge_fail)
} 
  
###############################################################################
#  MOM/SOM model simulations: Un-modeled heterogeneity (observers)
###############################################################################
  
# This routine conducts "distinct observer" simulations for MOM/SOM models without covariates and with un-modeled heterogeneity in the data arising from un-modeled differences in classification probabilities (theta) between secondary observers. Data are generated with distinct secondary observers, but are analyzed assuming identical classification probabilities among secondary observers. 
# Conduct analyses by selecting and executing all code in this section. 

# Data frame 'output' includes true parameters used to generate data, an indicator column named "heterogeneity" with value of 2, and statistics for estimated parameters. For classification probabilities, mean error and 95% confidence interval coverage are estimated by comparing overall estimates to true values for each distinct secondary observer.

if (model == "M") {

  # For model optimization function 'optim', define user-specified control parameters and lower/upper constraints for parameter values
  cont <- list(factr = 1e+07, maxit = 200, lmm = 10) # List of control parameters for 'optim' function
  constraints_logit <- qlogis(c(0.0001, 0.9999)) # Constraints for psi/theta parameters, transformed to logit scale
  constraints_g <- c(1.002, Inf) # Constraints for mean group size
  constraints_mix <- c(qlogis(0.001), Inf) # Constraints for group mixing parameters
  
  # Generate initial objects
  t_start <- t_loop <- unclass(Sys.time()) # For recording simulation duration
  error_msg <- converge_fail <- NULL # For error messages and convergence warnings
  
  parameters_names_global_true <- parameters_names_global
  
  # Build empty data frame 'output' containing columns for all estimated parameters and associated output 
  output_global <- matrix(numeric(1), 0, length(parameters_names_global) * 7 + 1) # Summarized simulation output
  output_names_global <-
    c(
      paste0("m.", parameters_names_global_true),
      paste0("sd.", parameters_names_global_true),
      paste0("se.", parameters_names_global_true),
      paste0("cv.", parameters_names_global_true),
      paste0("me.", parameters_names_global_true),
      paste0("rm.", parameters_names_global_true),
      paste0("ci.", parameters_names_global_true),
      "rep.tru"
    )
  colnames(output_global) <- output_names_global
  output_global <- data.frame(output_global)
  
  # Data frame 'sim_profiles_estimation' alters simulation profiles to remove classification parameters not included in the estimation model(s)
  sim_profiles <- mutate(sim_profiles, heterogeneity = 2)
  sim_profiles_estimation <- select(sim_profiles, -grep("theta_s2|theta_s3|theta_s4", profiles_names))
  parameters_names_global <- parameters_names_global_true[-c(grep("s2|s3|s4", parameters_names_global_true))]
  
  # Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
  for (sim in 1:dim(sim_profiles)[1]) {

    ## Define parameters for the current simulation
    O_ps <- c(sim_profiles$O_p[sim], sim_profiles$O_s[sim]) # Count of primary and secondary observers
    
    
    # 'parameters_names_true' and 'parameters_col' contain names and column #s of true parameter values in simulation profiles
    
    # Vector 'n_parameters' contains the number of true parameters in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), untransformed parameters (4), and heterogeneous group probabilities (pi) parameters (5). Vector 'n_parameters_estimated' contains corresponding numbers for estimated parameters. 
    
    parameters_names_true <- c(grep("psi|theta", parameters_names_global_true, value = TRUE))
    parameters_col_true <- c(grep("psi|theta", names(sim_profiles[sim, ])))
    
    n_parameters <- c(0, 0, length(parameters_names_true), 0, 0)
    
    # Add parameters for mean group size g
    if (any(sim_profiles[sim, grep("^g_", profiles_names)] > 1)) {
      parameters_names_true <-
        c(parameters_names_true, grep("^g_", profiles_names, value = TRUE))
      parameters_col_true <-
        c(parameters_col_true, grep("^g_", profiles_names))
      n_parameters[4] <- length(parameters_names_true) - n_parameters[3]
    }
    
    # Add parameters for heterogeneous group probabilities: heterogeneous group probabilities (pi) or heterogeneous group affinities (rho)
    if (any(sim_profiles[sim, mix_col] > 0)) {
      parameters_names_true <- c(parameters_names_true, profiles_names[mix_col])
      parameters_col_true <- c(parameters_col_true, grep("mix", profiles_names))
      n_parameters[5] <- length(mix_col)
    }
    
    n_parameters_estimated <- n_parameters
    n_parameters_estimated[3] <- n_parameters[3] - 
      (length(grep("psi|theta", parameters_names_true)) - length(grep("psi|theta_p|theta_s1", parameters_names_true)))
    parameters_names <- parameters_names_true[-c(grep("s2|s3|s4", parameters_names_true))]
    
    # Generate initial objects 
    parameters_ini <-
      generate.ini.f(sim_profiles_estimation[sim, ], n_parameters_estimated) # Sets of alternative initial parameter values on each row
    output_names <-
      c(
        paste0("m.", parameters_names_true),
        paste0("sd.", parameters_names_true),
        paste0("se.", parameters_names_true),
        paste0("cv.", parameters_names_true),
        paste0("me.", parameters_names_true),
        paste0("rm.", parameters_names_true),
        paste0("ci.", parameters_names_true),
        "rep.tru"
      )
    
    # Vector 'constraints_low' contains lower constraints for each parameter
    constraints_low <- c(rep(constraints_logit[1], n_parameters_estimated[3]), 
                         rep(constraints_g[1], n_parameters_estimated[4]), 
                         rep(constraints_mix[1], n_parameters_estimated[5]))
    
    # Upper constraints default to infinity (no constraints) 
    constraints_up <- Inf
    
    # List 'sim_data' contains formatted simulated survey observations generated using the true model (i.e. with distinct secondary observers) and (if needed) a keyed table of unique groups records
    sim_data <- generate.simulation.data.f(sim_profiles[sim, ]) 
    parameters_index <- 1 # Index for set of values used as initial parameters for optimization
    
    # ----- Model Optimization -----
    
    # Function for simulation progress bar
    id_max <- max(sim_data[[1]]$id)
    pb <- txtProgressBar(max = id_max, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Fit models using parallel execution of 'optim' function using "L-BFGS-B" method, with statistical output to list 'model_results'
    # Estimation models assume identical classification probabilities among secondary observers
    
    model_results <-
      as.data.frame(
        foreach(
          i = 1:id_max,
          .combine = rbind,
          .packages = c("plyr", "dplyr"),
          .errorhandling = "pass",
          .options.snow = opts
        ) %dopar% {
          optim(
            parameters_ini[parameters_index , ],
            optimize.M.f,
            gr = NULL,
            dat = sim_data[[1]][which(sim_data[[1]]$id == i), ],
            keys = sim_data[[2]],
            sim_profile = sim_profiles_estimation[sim,],
            hessian = TRUE,
            method = c("L-BFGS-B"),
            lower = constraints_low,
            upper = constraints_up,
            control = cont
          )
        }
      )
    
    # ----- Model Output: Error-checking -------
    
    # Check for model optimization errors/convergence failures
    fail <- failure.check.f(model_results, sim)
    
    # Attempt to re-fit non-optimized models
    if (length(fail[[1]]) > 0) {
      error_msg <- bind_rows(error_msg, fail[[2]])
      fail_tmp <- fail <- fail[[1]]
      while (parameters_index < 5 && length(fail) > 0) {
        parameters_index <-
          parameters_index + 1 # Increment index for alternate sets of initial values
        if (length(fail) == 1) {
          fail_tmp <-
            c(fail, fail) # >1 models needed for appropriate output formatting
        }
        
        # Try alternate sets of initial parameter values until model(s) optimize successfully
        
        model_results_refit <-
          as.data.frame(
            foreach(
              i = 1:length(fail_tmp),
              .combine = rbind,
              .packages = c("plyr", "dplyr"),
              .errorhandling = "pass"
            ) %dopar% {
              optim(
                parameters_ini[parameters_index , ],
                optimize.M.f,
                gr = NULL,
                dat = sim_data[[1]][which(sim_data[[1]]$id == fail_tmp[i]), ],
                keys = sim_data[[2]],
                sim_profile = sim_profiles_estimation[sim,],
                hessian = TRUE,
                method = c("L-BFGS-B"),
                lower = constraints_low,
                control = cont
              )
            }
          )
        
        # Add model(s) that optimize to 'model_results'
        model_results_refit <- refit.check.f(model_results_refit)
        success <-
          which(model_results_refit[1:length(fail), "convergence"] == 0)
        if (length(success) > 0) {
          l_ply(success, function(x) {
            for (i in 1:6)
              model_results[[i]][fail[x]] <<-
                model_results_refit[[i]][x]
          })
          fail <-
            fail[-success] # Remove successful models from failed models
        }
      }
      
      # If any models never optimize correctly, generate error/warning messages and remove model from output
      
      if (length(fail) > 0) {
        model_results <- model_results  %>% do(model_results[-fail, 1:6])
        converge_fail <-
          rbind(converge_fail, as_tibble(
            list(
              sim = rep(sim, length(fail)),
              id = fail,
              code = model_results_refit$convergence[which(model_results_refit$convergence > 0)][1:length(fail)],
              msg = model_results_refit$message[which(model_results_refit$convergence > 0)][1:length(fail)]
            )
          ))
        reps <- reps - length(fail)
      }
    }
    
    # ----- Model Output: Summarize -----
    
    # Back-transform estimates for logit parameters to probability scale and calculate summary statistics
    # Prefixes on variable names indicate summary statistics: m =  mean, sd = standard deviation, cv = coefficient of variation (sd / mean), me = mean error, rm = residual mean squared error, ci = 95% confidence interval coverage
    
    # True parameter values
    parameters_true <- sim_profiles[sim, parameters_col_true]
    
    # # of duplicates, start/end col 
    parameters_duplicate <- c(
      O_ps[2] - 1,
      sum(n_parameters_estimated[1:2]) + length(grep("theta_p", parameters_names_global_true)) + 1,
      sum(n_parameters_estimated[1:2]) + length(grep("theta_p|theta_s1", parameters_names_global_true))
    )
    
    # Append summary statistics for this distinct model as new row in data frame 'output'
    parameters_output <-
      model.results.het.s.f(model_results, parameters_true, reps, n_parameters, n_parameters_estimated, parameters_duplicate)
    
    output <-
      summarise.results.f(parameters_output, parameters_true, reps)
    names(output) <- output_names
    output_global <- bind_rows(output_global, data.frame(t(output)))
    
    # For long simulation runs, guard against data loss by periodically writing results to temporary file on hard drive 
    if (unclass(Sys.time()) - t_loop > 900) {
      try(
        write.csv(output_global, file = paste0("tmp_", out_filename))
      )
      t_loop <- unclass(Sys.time())
    }
    
    # Complete simulation loop
    reps <- sim_profiles$reps[1]
    cat("\n Completed Simulation ", sim, "\n")
    
    # With random seed assigned, append current simulation data to list
    if (any(profiles_names == "seed"))
      sim_data_list <- c(sim_data_list, list(sim_data))

  } # End of profile simulation loop
    
  # ----- Output & Reports -----
  
  # Upon completing simulations, combine simulation inputs and outputs, and write 'output' as .csv file format text file to user-specified path
  output_global <- bind_cols(sim_profiles, output_global)
  try(write.csv(output_global, file = out_filename)
      
  )
  
  # With random seed assigned, export simulation data
  if (any(profiles_names == "seed"))
    if (is.integer(sim_profiles$seed[1]) & sim_profiles$seed[1] > 0)
      try(saveRDS(sim_data_list, file = paste0(in_filename, "_output", "_sim_data.Rdata")))
  
  # Print reports on simulation speed and errors/warnings to the console
  elapsed <- unclass(Sys.time()) - t_start
  cat(round(elapsed / 60, 2),
      "m ",
      round(elapsed, 0),
      "s ",
      round(elapsed / (reps * dim(sim_profiles)[1]), 2),
      "s/rep \n")
  print.error.f(error_msg, converge_fail)
    
  } 

###############################################################################
#       Multinomial model simulations: Estimating observed proportions
###############################################################################

# This routine conducts "omnibus" simulations for multinomial logit models estimating observed species proportions (proportion of each species among individuals classified to species by primary observers), assuming no species misidentification. Data are generated using true parameter values for MOM models, but are analyzed using models that estimates observed species proportions from counts of individuals classified as each species by primary observers. 

## Define global objects (constant across all simulations) 

# Add field 'observed_p' with value 1 to sim_profiles to denote in output that estimates are for observed proportions
sim_profiles <- mutate(sim_profiles, observed_p = 1L)

# Vector 'n_parameters_output_global' contains the number of estimated parameters (for the multinomial logit model) across all simulation in each of 5 categories: regression coefficient intercept and slope parameters (1 & 2), logit parameters (3), un-transformed parameters (4), and heterogeneous group probabilities (pi) parameters (5). Only catergories 3 (observed species proportions) and 4 (mean group sizes) are used here. 

parameters_psi <- grep("psi", parameters_names_global)
parameters_group_size <- grep("^g_[0123456789]", parameters_names_global)
n_parameters_output_global <- c(0, 0, (B - 1), length(parameters_group_size), 0)
parameters_names_output_global <- c(parameters_names_global[c(parameters_psi, parameters_group_size)])

# Build empty data frame 'output_global' containing columns for global estimated parameters (across all simulations) and associated output 
output_global <- matrix(0, 0, sum(n_parameters_output_global) * 7 + 1)
output_names_global <-
  c(
    paste0("m.", parameters_names_output_global),
    paste0("sd.", parameters_names_output_global),
    paste0("se.", parameters_names_output_global),
    paste0("cv.", parameters_names_output_global),
    paste0("me.", parameters_names_output_global),
    paste0("rm.", parameters_names_output_global),
    paste0("ci.", parameters_names_output_global),
    "rep.tru"
  )
colnames(output_global) <- output_names_global
output_global <- data.frame(output_global)

# Loop for sequentially executing simulations for distinct models on each row of 'sim_profiles'
for (sim in 1:dim(sim_profiles)[1]) {
  
  # List 'sim_data' contains formatted simulated survey observations generated using the true model (including uncertain identification, heterogeneous groups)
  sim_data <- generate.simulation.data.f(sim_profiles[sim, ]) 
  
  # Matrix 'count.mat' summarizes observations as counts of each species
  count.mat <- 
    data.frame(llply(1:A, function(x) sim_data[[1]][x] * sim_data[[1]]$count)) %>%
    bind_cols(., data.frame(id = sim_data[[1]]$id, count = sim_data[[1]]$count)) %>%
    group_by(., id) %>% 
    summarise_all(., sum) %>%
    select(., all_of(2:(A + 1)))
  
  # Calculate mean group sizes if any groups sizes can exceed one
  size.groups <- size.groups.SE <- NULL
  if (any(sim_profiles[sim, grep("^g_[0123456789]", profiles_names)] > 1)) {
    if (B < A) {
      # No estimates of mean group size with partial identification
      size.groups <- size.groups.SE <- matrix(0, reps, B)
    }else{
      # For B = A, add estimated mean group sizes 
      spp.groups <- data.frame(llply(1:B, function(x) (sim_data[[1]][x] > 0) * sim_data[[1]]$count)) %>%
        bind_cols(., data.frame("id" = sim_data[[1]]$id)) %>%
        group_by(., id) %>%
        summarise_all(., sum) %>%
        select(., all_of(2:(A + 1)))
      
      size.groups <- count.mat/spp.groups
      
      squared.sum <- 
        data.frame(sim_data[[1]][1:A]^2 * sim_data[[1]]$count) %>%
        bind_cols(., data.frame(id = sim_data[[1]]$id)) %>%
        group_by(., id) %>%
        summarise_all(., sum) %>%
        select(., all_of(2:(A + 1)))
      
      size.groups.SE <- ((squared.sum - count.mat ^ 2 / spp.groups) / (spp.groups - 1)) ^ 0.5 / spp.groups ^ 0.5
    }
  } # End of group size estimation loop
  
  # Check for partial identification
  if (B < A) {
    # If present, pro-rate partially-identified individuals to species-specific counts based on observed species proportions
    for (i in 1:B) {
      count.mat[, i] <- 
        count.mat[, i] + count.mat[, A] * (count.mat[, i] / rowSums(count.mat[, 1:B]))
    }
    count.mat <- select(count.mat, all_of(1:B))
  }
  
  count.mat <- as.matrix(count.mat)
  
  # Rearrange 'count.mat' so last species becomes reference category for multinomial estimation of observed species proportions
  model_results <- apply(count.mat[, c(B, 1:(B - 1))], 1, multinom.likelihood.f)
  
  # Model error-checking
  converge <- fail <- NULL
  for (i in 1:reps) {
    converge <- c(converge, model_results[[i]]$convergence)
  }
  if (any(converge != 0)) {
    fail <- which(converge != 0)
    model_results <- model_results[-fail]
    size.groups <- size.groups[-fail, ]
    size.groups.SE <- size.groups.SE[-fail, ]
  }
  
  # True parameter values for estimated parameters (natural scale) from the data generating model
  parameters_group_size <- grep("^g_[0123456789]", names(sim_profiles))
  parameters_psi <- grep("psi", names(sim_profiles))
  if (any(sim_profiles[sim, parameters_group_size] > 1)) {
    parameters_true <- sim_profiles[sim, c(parameters_psi, parameters_group_size)]
    parameters_names_output <- names(sim_profiles)[c(parameters_psi, parameters_group_size)]
    n_parameters_output <- c(0, 0, (B - 1), length(parameters_group_size), 0)
  }else{
    parameters_true <- sim_profiles[sim, c(parameters_psi)]
    parameters_names_output <- names(sim_profiles)[c(parameters_psi)]
    n_parameters_output <- c(0, 0, (B - 1), 0, 0)
  }
  
  # Extract parameter estimates (natural scale) and standard errors from 'model_results', append estimated group sizes
  # Array 'estimates.array' contains parameter estimates in array 1, standard errors in array 2
  rep.tru <- length(model_results)
  estimates.array <- 
    array(0, dim = c(rep.tru, sum(n_parameters_output), 2))
  
  estimates.array[, , 1] <- 
    c(
      aaply(
        matrix(
          vapply(1:rep.tru, function(x) summary(model_results[[x]])$coefficients, numeric(B - 1)),
          ncol = (B - 1), byrow = TRUE), 
        1, 
        multinomial.inv.f
      )[, 1:(B - 1)], 
      unlist(size.groups)
    )
  
  Hessian.mat <-
    llply(1:rep.tru, function(x) solve(summary(model_results[[x]])$Hessian))
  
  estimates.array[, , 2] <- 
    as.matrix(
      cbind(
        t(vapply(
          1:rep.tru, 
          function(x)
            diag(DeltaMethod(estimates.array[x, 1:(B - 1), 1], multinomial.inv.f, Hessian.mat[[x]], 0.00001)$variance) ^ 0.5, 
          numeric((B))
        ))[ ,1:(B - 1), drop = FALSE],
        size.groups.SE
      )
    )
  
  # Compute 95% confidence interval coverage
  ci.cov <- 
    aaply(estimates.array[, , 1] - 1.96 * estimates.array[, , 2], 1, function(x) 
      parameters_true > x) *
    aaply(estimates.array[, , 1] + 1.96 * estimates.array[, , 2], 1, function(x) 
      parameters_true < x)
  
  # Place summary statistics for this distinct model into data frame 'output'
  # Special case for # of parameters = 1 to avoid inferno of dropped array dimensions
  if (B == 2 & !any(sim_profiles[sim, grep("^g_[0123456789]", profiles_names)] > 1)) {
    # With groups = 1
    output <- 
      unlist(c(
        colMeans(estimates.array[, , 1, drop = FALSE]),
        aaply(estimates.array[, , 1, drop = FALSE], 2, sd),
        colMeans(estimates.array[, , 2, drop = FALSE]),
        aaply(estimates.array[, , 1, drop = FALSE], 2, sd) / colMeans(estimates.array[, , 1, drop = FALSE]),
        colMeans(estimates.array[, , 1, drop = FALSE]) - parameters_true,
        sum((estimates.array[, 1, 1] -
               unlist(parameters_true)) ^ 2) / (rep.tru - 1),
        sum(ci.cov, na.rm = TRUE) / rep.tru,
        rep.tru))
  }else{
    # Group size >= 1
    output <- 
      unlist(c(
        colMeans(estimates.array[, , 1]),
        aaply(estimates.array[, , 1], 2, sd),
        colMeans(estimates.array[, , 2]),
        aaply(estimates.array[, , 1], 2, sd) / colMeans(estimates.array[, , 1]),
        colMeans(estimates.array[, , 1]) - parameters_true,
        colSums((estimates.array[, , 1] -
                   matrix(
                     unlist(parameters_true),
                     byrow = TRUE,
                     nrow = rep.tru,
                     ncol = length(parameters_true))) ^ 2) / (rep.tru - 1), 
        aaply(ci.cov, 2, function(x) sum(x, na.rm = TRUE) / rep.tru),
        rep.tru))
  } # End of summary loop
  
  output_names <-
    c(
      paste0("m.", parameters_names_output),
      paste0("sd.", parameters_names_output),
      paste0("se.", parameters_names_output),
      paste0("cv.", parameters_names_output),
      paste0("me.", parameters_names_output),
      paste0("rm.", parameters_names_output),
      paste0("ci.", parameters_names_output),
      "rep.tru"
    )
  
  names(output) <- output_names
  
  # Append output for this distinct model to the global output
  output_global <- bind_rows(output_global, data.frame(t(output)))
  
  reps <- sim_profiles$reps[1]
  cat("\n Completed Simulation ", sim, "\n")
} # End of profile simulation loop

# Upon completing simulations, combine simulation inputs and outputs, and write 'output' as .csv file format text file to user-specified file name
output_global <- bind_cols(sim_profiles, output_global)
write.csv(output_global, file = paste0("observed_proportions_", out_filename))     
