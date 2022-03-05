# likelihood_equations.R
# Generates R list objects with pre-computed values for likelihood equations, version 1.2.4
# Steven T. Hoekman, Wild Ginger Consulting, PO Box 182 Langley, WA 98260, steven.hoekman@protonmail.com

# This code generates list objects 'likelihood_equations' composed of three types of elements: 1) 'true.combinations.g.[size]' containing a matrix with possible combinations of true groups for groups of each [size], 2) 'true.permutations.count.g.[size]' containing a vector with counts of possible permutations for each combination for a group true group of each [size], and 3) 'observed.[observed group]' containing a matrix summarizing equations for computing probabilities of the [observed group] specified by counts of individuals classified to each true species state. These lists substantially increase speed of model optimization. Required lists for conducting simulation study described in the companion article and MetadataS3 are provided in DataS2. Code below produces these lists to user specifications for the number of observation states (A), number of true species states (B), and the range of group sizes (g). See MetadataS3.pdf in for addition details. Code developed and tested in R version 4.1.

# This code does NOT need to be executed prior to conducting simulation analyses in 'simulations_R_code.R'

###############################################################################
#             Required R libraries
###############################################################################

# These libraries must be installed and loaded prior to executing code

# library(fastverse)  High-performance statistical computing and data manipulation
#   Utilized fastverse packages are 'collapse', 'matrixStats', & 'magrittr'
# library(doSNOW)     # Back end for parallel execution
# library(gtools)     # Numeric tools

library(fastverse)

fastverse_detach(data.table, kit, fst)
fastverse_extend(doSNOW, gtools)

###############################################################################
#     Register the 'doSNOW' parallel execution back end
###############################################################################

# Make cluster 'cl' with user-specified # (>0) of CPU workers. Greater than 1 allows parallel processing, but # of workers typically shouldn't exceed the number of available CPU threads.

cl <- makeCluster(4)
registerDoSNOW(cl)

# AFTER completing all simulation analyses, un-register the cluster 'cl' to remove the workers from RAM.
stopCluster(cl) 

###############################################################################
#             Required supporting functions
###############################################################################

# These functions must be loaded prior to executing code

# Function: {terms.f} For a given group size g, accepts combinations of possible true groups 'B.cmbn', corresponding counts of each true species state in each combination 'B.ct', the row index 's' specifying a true species state, and combinations of possible observation states 'A.cmbn'. 
# Returns a matrix with integers 'index', 'coefficient', and A * B columns of counts summarizing possible true and observed states for an observed group. See 'Generate new likelihood.equations list' section for details.  

terms.f <- function(s , A.cmbn, B.cmbn, B.ct){ 
  
  B.num <- 
    B.ct[s, which(B.ct[s,] > 0)]
  
  X <-
    unique(combinations(
      length(A.cmbn),
      B.num[1],
      A.cmbn,
      set = FALSE,
      repeats.allowed = FALSE
    ))
  
  if (length(B.num) > 1) {
    
    for (b in 2:length(B.num)) {
      
      X <-
        do.call(rbind, lapply(1:nrow(X), function(x)
          B.byte.f(X[x,], B.num[b])))
      
    }
  }
  
  X <- 
    t((X - 1L) * B) + B.cmbn[s,]
  
  X <-
    t(apply(X, 2, function(y)
      vapply(1:(B * A), function(z)
        sum(y == z), integer(1))))
  
  matrix(c(rep(s, dim(X)[1]),
             apply(X, 1, function(y)
               factorial(g) / prod(factorial(y)))
             / list.new[[paste0("true.permutations.count.g.", g, collapse = ".")]][s]
           ,
           X),
         ncol = dim(X)[2] + 2) %>%
    setColnames(c("index", "coefficient", names.vec))
  
}     

# Function: {B.byte.f} Accepts a partial row of matrix terms from terms.f and adds all possible unique combinations of b additional terms. 

B.byte.f <- function(x, b){
  
  v <-
    observed.combinations[[ndex]][os,-(pmatch(x, observed.combinations[[ndex]][os,], duplicates.ok = FALSE))]
  
  new.col <-
    unique(combinations(length(v), b, v, set = FALSE, repeats.allowed = FALSE))
  
  cbind(matrix(
    x,
    ncol = length(x),
    nrow = dim(new.col)[1],
    byrow = TRUE
  ) , new.col)
}


###############################################################################
#             User-specified inputs
###############################################################################

# User-specified inputs: 
# B = # of true species states
# A = # of observation states
# g.limits = lower, upper limits for group size g 

B <- 2L
A <- 2L
g.limits <- c(1L, 5L)

###############################################################################
#               Generate lists of combinations
###############################################################################

# Generate combinations of true species states and observation states
# Lists 'true.combinations' and 'observed.combinations' contain a list element for each group size g containing a matrix with each possible combination of true species states (B) or observation states (A) on separate row, with columns corresponding to each group individual.
# Lists 'true.combinations.count' and 'observed.combinations.count'  contain a list element for each group size g containing a matrix with each possible combination of true species states (B) or observation states (A) on separate row, with columns enumerating the count of individuals in each true species state 1 to B or each observation state 1 to A. 

true.combinations <-
  lapply(g.limits[1]:g.limits[2], function(x)
    combinations(B, x, repeats.allowed = TRUE))

observed.combinations <-
  lapply(g.limits[1]:g.limits[2], function(x)
    combinations(A, x, repeats.allowed = TRUE))

true.combinations.count <-
  lapply(true.combinations, rowTabulates)

observed.combinations.count <-
  lapply(observed.combinations, rowTabulates)

# For each group size, list 'true.permutations.count' enumerates (via the multinomial coefficient) the count of possible permutations for each true group.

true.permutations.count <- 
  lapply(true.combinations.count, function(x)
    apply(x, 1, function(y)
      as.integer(factorial(sum(y)) / prod(factorial(y)))))

names(true.combinations) <- 
  paste("true.combinations.g.", g.limits[1]:g.limits[2], sep = "")

names(true.permutations.count) <- 
  paste("true.permutations.count.g.", g.limits[1]:g.limits[2], sep = "")


###############################################################################
#           Generate a list
###############################################################################

# Generate a list with elements for each possible observed group, named by the count of individuals in each observed state. Each element contains a matrix where each row corresponds to a possible permutation of observation states for a given true group. Matrix columns are:
# index', an integer corresponding to possible combinations of true groups (in 'true.combinations'). Separate rows with the same index value correspond to different possible permutations of classification states potentially giving rise to the observed group for that true group. 
# 'coeff', an integer giving possible permutations of classification states (the multinomial coefficient).
# 'V1, V2,  ...', integers giving a vector of counts for true species states versus observed states. 
# These vectors are one-dimensional representations of matrices having B rows for true species states and A columns for observation states, with the cell on row b and column a representing the count of individuals of true species state b classified to observation state a. (This matrix arrangement can be visualized using Tables 1 and 2 in the companion article) Storing these matrices in a vector (starting from the top-left cell and moving left to right along each row to finish at the bottom-right cell) increases computational efficiency. 

t.start <- unclass(Sys.time()) # For recording duration of list generation
names.vec <- unlist(lapply(1:A, function(x)
  paste0(x, "." , 1:B)))

observed.0 <-
  matrix(c(1L, 1L, rep(0L, A * B)), nrow = 1) %>%
  setColnames(c("index", "coefficient", names.vec)) %>%
  list(.) %>%
  setColnames(paste0("observed.", paste0(rep(0, A), collapse = ".")))

true.permutations.count.g.0 <- 1L %>%
  list(.) %>%
  setColnames("true.permutations.count.g.0")

list.new <- 
  c(true.combinations, true.permutations.count.g.0, true.permutations.count, observed.0) # list.new will contain the generated list

ndex <- 1 # 'ndex' indexes the group size g from 1 to the (upper limit - lower limit)

for (g in g.limits[1]:g.limits[2]) { # Loop for group size
  t.loop <- unclass(Sys.time()) # Duration of each loop
  o.states <- dim(observed.combinations[[ndex]])[1] # Number of possible observation states
  terms <- list(NULL)
  
  # Loop through possible observation states
  output <-
    foreach(
      os = 1:o.states,
      .packages = c("magrittr", "gtools", "collapse"),
      .inorder = TRUE,
      .errorhandling = "pass"
    ) %dopar% {
      terms <-
        do.call(
          rbind,
          lapply(
            1:dim(true.combinations[[ndex]])[1],
            terms.f,
            observed.combinations[[ndex]][os, ],
            true.combinations[[ndex]],
            true.combinations.count[[ndex]]
          )
        )
    }
  
  # Set storage mode to "integer" if coefficient values don't exceed numeric limit of 2e10^9
  
  for (q in seq_along(output)) {
    if (max(output[[q]]) < 2.1e9)
      storage.mode(output[[q]]) <- "integer"
  }
  
  # Name each list element by the count of individuals in each observation state 1 to A (separated by a decimal), append 'output' to the new list
  names(output) <- sapply(1:o.states, function(x)
    paste0("observed.", paste0(observed.combinations.count[[ndex]][x,], collapse = ".")))
  
  list.new <- c(list.new, output)
  cat("completed group size:", g, "duration (m):", round((unclass(Sys.time()) - t.loop)/60, 2), "\n")
  ndex <- ndex + 1 
}

cat("Total duration (m): ", round((unclass(Sys.time()) - t.start)/60, 2), "\n")


###############################################################################
#                  List management
###############################################################################

# Rename the new list and save. Lists can also be joined. 
likelihood.equations <- list.new

# Separate lists with the same number of true species states and observation states may be merged. 
# 'duplicate_entries' removes identical entries from the new list
# duplicate_entries <- duplicated(c(likelihood.equations, list.new))
# likelihood.equations <-c(likelihood.equations, list.new)[!duplicate_entries]

# Save the list for future use
save(likelihood.equations, file = "likelihood_equations_a2b2g5.RData")
