#' @name hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Hierarchical Embedded Bayesian Additive Regression Trees
#' @description This function runs a hebart model and returns the tree and
#' other results obtained in the last iteration of the MCMC
#' @param formula The model formula
#' @param data The data to be used in the modeling
#' @param group_variable The name of the grouping variable 
#' @param num_trees The number of trees (P)
#' @param control A list with control settings
#' @param priors A list with prior hyperparameters as defined by the model
#' @param inits A list with initial values for parameters
#' @param MCMC A list with MCMC parameters
#' @param k_1_pars A list to decide whether to sample k1 or not and if yes, from
#' which range
#' @return A list containing:
#'  Everything 
#' @details
#' Priors used ----------------------------------------------------------
#' y_{ij} ~ Normal(m_j, tau^-{1})
#' tau    ~ Gamma(alpha, beta)
#' mu     ~ Normal(0, tau_mu = k2*tau^-{1})
#' mu_j   ~ Normal(mu, k1*tau^-{1})
#' ----------------------------------------------------------------------

hebart <- function(
  formula, 
  data, 
  group_variable, # groups is the group number of each obs
  num_trees = 1, # Number of trees
  control = list(node_min_size = 5), # Size of smallest nodes
  priors = list(
    alpha = 0.95, # Prior control list
    beta = 2,
    k_1 = 1e-10,
    k_2 = 0.2,
    nu = 3,
    lambda = 0.1
  ),
  inits = list(tau = 1), # Initial values list
  MCMC = list(
    iter = 500, # Number of iterations
    burn = 150, # Size of burn in
    thin = 1
  ), # Amount of thinning
  k_1_pars = list(sample_k1 = TRUE, 
                  min_u     = 0, 
                  max_u     = 15, 
                  k1_prior  = TRUE)
) {
  
  if(MCMC$iter <= MCMC$burn){
    stop("Number of burn-in iterations is not smaller
than the total number of iterations")
  }
  
  # Handling formula interface
  formula_int <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula_int)[1]
  names_x <- all.vars(formula_int[[3]])
  
  #data    <- dplyr::select(data, c(!!response_name, !!names_x, !!group_variable))
  data    <- dplyr::select(data, c(!!response_name, !!names_x, "group"))
  names(data)[names(data) == group_variable] <- "group"
  groups  <- data$group
  mf      <- stats::model.frame(formula_int, data = data)
  X       <- as.matrix(stats::model.matrix(formula_int, mf))
  y       <- stats::model.extract(mf, "response")
  
  # Extract control parameters
  node_min_size <- control$node_min_size
  
  # Extract hyper-parameters
  alpha   <- priors$alpha # Tree shape parameter 1
  beta    <- priors$beta # Tree shape parameter 2
  k_1     <- priors$k_1 # Group standard deviation multiplier
  k_2     <- priors$k_2 # Overall mean precision multiplier
  nu      <- priors$nu # Parameter 1 for precision
  lambda  <- priors$lambda # Parameter 2 for precision
  
  # Extract initial values
  tau <- inits$tau
  sigma <- 1 / sqrt(tau)
  log_lik <- 0
  
  # Extract MCMC details
  iter <- MCMC$iter # Number of iterations
  burn <- MCMC$burn # Size of burn in
  thin <- MCMC$thin # Amount of thinning
  
  # k_1 sampling parameters
  sample_k1 <- k_1_pars$sample_k1
  min_u     <- k_1_pars$min_u
  max_u     <- k_1_pars$max_u
  k1_prior  <- k_1_pars$k1_prior
  
  # Storage containers
  store_size      <- (iter - burn) / thin
  tree_store      <- vector("list", store_size)
  sigma_store     <- rep(NA, store_size)
  y_hat_store     <- matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store   <- rep(NA, store_size)
  full_cond_store <- matrix(NA, ncol = num_trees, nrow = store_size)
  samples_k1      <- rep(NA, store_size)
  
  # Scale the response target variable
  y_mean   <- mean(y)
  y_sd     <- stats::sd(y)
  y_scale  <- (y - y_mean) / y_sd
  n        <- length(y_scale)
  
  # Get the group matrix M
  M           <- stats::model.matrix(~ factor(groups) - 1)
  group_sizes <- table(groups)
  num_groups  <- length(group_sizes)
  
  # Create a list of trees for the initial stump
  curr_trees <- create_stump(
    num_trees  = num_trees,
    num_groups = num_groups,
    y          = y_scale,
    X          = X
  )
  # predictions <- get_predictions(curr_trees, X, single_tree = num_trees == 1)
  predictions <- get_group_predictions(curr_trees, X, groups, single_tree = num_trees == 1)
  
  # Set up a progress bar
  pb <- utils::txtProgressBar(
    min = 1, max = iter,
    style = 3, width = 60,
    title = "Running rBART..."
  )
  
  # Start the iterations loop
  for (i in 1:iter) {
    utils::setTxtProgressBar(pb, i)
    
    # If at the right place store everything
    if ((i > burn) & ((i %% thin) == 0)) {
      curr <- (i - burn) / thin
      tree_store[[curr]] <- curr_trees
      sigma_store[curr] <- sigma
      y_hat_store[curr, ] <- predictions
      log_lik_store[curr] <- log_lik
    }
    
    # Start looping through trees
    for (j in 1:num_trees) {
      
      # Calculate partial residuals for current tree
      if (num_trees > 1) {
        partial_trees <- curr_trees
        partial_trees[[j]] <- NULL # Blank out that element of the list
        current_partial_residuals <- y_scale -
          get_group_predictions(partial_trees, X, groups,
                                single_tree = num_trees == 2
          )
        # current_partial_residuals <- y_scale -
        # get_predictions(partial_trees, X,
        #                 single_tree = num_trees == 2
        # )
      } else {
        current_partial_residuals <- y_scale
      }
      
      # Propose a new tree via grow/change/prune/swap
      #type <- sample(c("grow", "prune", "change", "swap"), 1)
      type <- sample(c("grow", "prune"), 1)
      if (i < max(floor(0.1 * burn), 10)) type <- "grow" # Grow for the first few iterations
      
      # Get a new tree!
      new_trees <- curr_trees
      new_trees[[j]] <- update_tree(
        y = y_scale,
        X = X,
        num_groups = num_groups,
        type = type,
        curr_tree = curr_trees[[j]],
        node_min_size = node_min_size
      )
      
      # Calculate the complete conditional and acceptance probability
      l_new <- tree_full_conditional_hebart(
        tree = new_trees[[j]],
        R = current_partial_residuals,
        k_1,
        k_2,
        M,
        nu,
        lambda
      ) +
        get_tree_prior(new_trees[[j]], alpha, beta)
      
      l_old <- tree_full_conditional_hebart(
        curr_trees[[j]],
        current_partial_residuals,
        k_1,
        k_2,
        M,
        nu,
        lambda
      ) +
        get_tree_prior(curr_trees[[j]], alpha, beta)
      
      # cat('tree', j,'\n')
      # cat('l_new = ',l_new,'; l_old = ',l_old,'\n')
      
      # If accepting a new tree update all relevant parts
      a <- exp(l_new - l_old)
      
      if ((i > burn) & ((i %% thin) == 0)) {
        full_cond_store[curr, j] <- l_old
      }
      if (a > stats::runif(1)) {
        # Make changes if accept
        curr_trees <- new_trees
      } # End of accept if statement
      
      # Update mu whether tree accepted or not
      curr_trees[[j]] <- simulate_mu_hebart(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        tau, k_1, k_2
      )
      
      # Finally update the group means:
      curr_trees[[j]] <- simulate_mu_groups_hebart(
        tree = curr_trees[[j]],
        current_partial_residuals,
        groups,
        tau, k_1, k_2
      )
      
      # Check the trees
      if (any(curr_trees$tree_matrix[, "node_size"] < node_min_size)) browser()
    } # End loop through trees
    
    # Calculate full set of predictions
    predictions <- get_group_predictions(curr_trees, X, groups, single_tree = num_trees == 1)
    # predictions <- get_predictions(curr_trees, X, single_tree = num_trees == 1)
    S <- sum((y_scale - predictions)^2)
    
    # Update tau and sigma
    tau <- update_tau(S, nu, lambda,
                      n = length(y_scale)
    )
    sigma <- 1 / sqrt(tau)
    
    # Sample k1
    if(sample_k1){
      # We can set these parameters more smartly
      sampled_k1 <- update_k1(y, min_u, max_u, k_1, k_2, M, nu, lambda, prior = k1_prior)
      
      samples_k1[i] <- k_1
      if(sampled_k1 != k_1){ samples_k1[i] <- k_1 <- sampled_k1 }
    }
    
    # Get the overall log likelihood
    log_lik <- sum(stats::dnorm(y_scale, mean = predictions, sd = sigma, log = TRUE))
  } # End iterations loop
  cat("\n") # Make sure progress bar ends on a new line
  
  
  result <- list(
    formula   = formula, 
    trees     = tree_store,
    sigma     = sigma_store,
    y_hat     = y_hat_store * y_sd + y_mean,
    log_lik   = log_lik_store,
    full_cond = full_cond_store,
    y = y,
    X = X,
    iter = iter,
    burn = burn,
    thin = thin,
    store_size = store_size,
    num_trees  = num_trees
  )
  
  # Saving the k_1 samples 
  if(sample_k1) result$samples_k1 <-  samples_k1

  
  # RMSE calculation
  pred <- predict_hebart(X, groups, result,
                        type = "mean")
  mse                 <- mean((pred - y)^2)
  r.squared           <- 1 - mse / stats::var(y)
  
  result$mse           <- mse
  result$r.squared     <- r.squared
  result$num_variables <- length(names_x)
  
  class(result) <- "hebart"
  
  
  return(result = result)
} # End main function
