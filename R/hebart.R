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
#' @return A list containing:
#'  Everything
#' @details
#' Priors used ----------------------------------------------------------
#' y_{ij} ~ Normal(m_j, tau^-1)
#' tau    ~ Gamma(nu/2, nu*lambda/2)
#' mu     ~ Normal(0, tau_mu^-1)
#' phi    ~ Normal(mu, sigma_phi^2 / T)
#' ----------------------------------------------------------------------

hebart <- function(formula,
                   data,
                   group_variable,
                   # X is the feature matrix, y is the target,
                   # groups, # groups is the group number of each obs
                   num_trees = 2, # Number of trees
                   control = list(node_min_size = 5), # Size of smallest nodes
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     nu = 2,
                     lambda = 0.1,
                     tau_mu = 16 * num_trees,
                     shape_sigma_phi = 0.5,
                     scale_sigma_phi = 1
                   ),
                   inits = list(
                     tau = 1,
                     tau_phi = 1
                   ), # Initial values list
                   MCMC = list(
                     iter = 250, # Number of iterations
                     burn = 50, # Size of burn in
                     thin = 1,
                     sigma_phi_sd = 2
                   )) {

  # Handling formula interface
  formula_int <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula_int)[1]
  names_x <- all.vars(formula_int[[3]])

  # Used in create_S to avoid error with stumps
  mod.mat <<- function(f) {
    if(nlevels(f) == 1) {
      m <- matrix(1, nrow = length(f), ncol = 1)
    } else {
      m <- stats::model.matrix(~ f - 1)  
    }
    return(m)
  }
  
  #-------------------------------------------------------
  data <- dplyr::select(data, c(!!response_name, !!names_x, !!group_variable))
  # data    <- dplyr::select(data, c(!!response_name, !!names_x, "group"))
  names(data)[names(data) == group_variable] <- "group"
  groups <- data$group
  mf <- stats::model.frame(formula_int, data = data)
  X <- as.matrix(stats::model.matrix(formula_int, mf))
  y <- stats::model.extract(mf, "response")

  # Extract control parameters
  node_min_size <- control$node_min_size

  # Extract hyper-parameters
  alpha <- priors$alpha # Tree shape parameter 1
  beta <- priors$beta # Tree shape parameter 2
  nu <- priors$nu # Parameter 1 for precision
  lambda <- priors$lambda # Parameter 2 for precision
  tau_mu <- priors$tau_mu # Overall mean precision
  shape_sigma_phi <- priors$shape_sigma_phi # Weibull prior parameters
  scale_sigma_phi <- priors$scale_sigma_phi

  # Extract initial values
  tau <- inits$tau
  sigma <- 1 / sqrt(tau)
  sigma_phi <- inits$sigma_phi
  tau_phi <- 1 / (sigma_phi^2)
  log_lik <- 0

  # Extract MCMC details
  iter <- MCMC$iter # Number of iterations
  burn <- MCMC$burn # Size of burn in
  thin <- MCMC$thin # Amount of thinning
  sigma_phi_sd <- MCMC$sigma_phi_sd # SD parameter for sigma_phi MH update
  
  # Storage containers
  store_size <- (iter - burn) / thin
  tree_store <- vector("list", store_size)
  sigma_store <- rep(NA, store_size)
  y_hat_store <- matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store <- rep(NA, store_size)
  sigma_phi_store <- rep(NA, store_size)

  # Scale the response target variable
  y_min <- min(y)
  y_max <- max(y)
  y_scale <- (y - y_min)/(y_max - y_min) - 0.5
  n <- length(y_scale)

  # Get the group matrix M
  M <- stats::model.matrix(~ factor(groups) - 1)
  group_sizes <- table(groups)
  num_groups <- length(group_sizes)

  # Create a list of trees for the initial stump
  curr_trees <- create_stump(
    num_trees = num_trees,
    groups = groups,
    y = y_scale,
    X = X
  )
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
      sigma_phi_store[curr] <- sigma_phi
    }

    # Start looping through trees
    for (j in 1:num_trees) {

      # Calculate partial residuals for current tree
      if (num_trees > 1) {
        partial_trees <- curr_trees
        partial_trees[[j]] <- NULL # Blank out that element of the list
        current_partial_residuals <- y_scale -
          get_group_predictions(
            trees = partial_trees, X, groups,
            single_tree = num_trees == 2
          )
      } else {
        current_partial_residuals <- y_scale
      }

      # Propose a new tree via grow/change/prune/swap
      type <- sample(c("grow", "prune", "change", "swap"), 1)
      if (i < max(floor(0.1 * burn), 10)) type <- "grow" # Grow for the first few iterations

      # Get a new tree!
      new_trees <- curr_trees

      new_trees[[j]] <- update_tree(
        y = y_scale,
        X = X,
        groups = groups,
        type = type,
        curr_tree = curr_trees[[j]],
        node_min_size = node_min_size
      )
      # Calculate the complete conditional and acceptance probability
      l_new <- full_conditional_hebart(
        tree = new_trees[[j]],
        R = current_partial_residuals,
        num_trees,
        tau,
        tau_phi,
        tau_mu,
        M
      ) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      l_old <- full_conditional_hebart(
        curr_trees[[j]],
        current_partial_residuals,
        num_trees,
        tau,
        tau_phi,
        tau_mu,
        M
      ) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # If accepting a new tree update all relevant parts
      log.alpha <- (l_new - l_old)
      accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
      if (accept) curr_trees <- new_trees
      
      # Update mu whether tree accepted or not
      curr_trees[[j]] <- simulate_mu_hebart(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        tau,
        tau_phi,
        tau_mu,
        M,
        num_trees
      )

      # Update phi as well
      curr_trees[[j]] <- simulate_phi_hebart(
        tree = curr_trees[[j]],
        R = current_partial_residuals,
        groups,
        tau,
        tau_phi,
        M,
        num_trees
      )

      # Check the trees
      if (any(curr_trees$tree_matrix[, "node_size"] < node_min_size)) browser()
    } # End loop through trees

    # Calculate full set of predictions
    predictions <- get_group_predictions(curr_trees, 
                                         X, 
                                         groups, 
                                         single_tree = num_trees == 1
    )

    # Update tau
    tau <- update_tau(
      y_scale,
      predictions,
      nu,
      lambda
    )
    sigma <- 1 / sqrt(tau)

    # Update tau_phi
    S1 <- create_S(curr_trees, groups)
    S2 <- create_S(curr_trees)
    
    sigma_phi <- update_sigma_phi(
      y_scale, S1, S2, sigma_phi, tau_mu, tau, shape_sigma_phi, scale_sigma_phi, 
      num_trees, sigma_phi_sd
    )
    tau_phi <- 1 / (sigma_phi^2)

    # Get the overall log likelihood
    Omega_y <- diag(n)/tau + tcrossprod(S1)/(num_trees * tau_phi) + 
      tcrossprod(S2)/tau_mu
    log_lik <- mvnfast::dmvn(
      y, rep(0, n), Omega_y, log = TRUE)
      
  } # End iterations loop
  cat("\n") # Make sure progress bar ends on a new line

  result <- list(
    trees = tree_store,
    sigma = sigma_store,
    y_hat = (y_hat_store + 0.5) * (max(y) - min(y)) + min(y),
    log_lik = log_lik_store,
    sigma_phi = sigma_phi_store,
    y = y,
    X = X,
    iter = iter,
    burn = burn,
    thin = thin,
    store_size = store_size,
    num_trees = num_trees,
    formula = formula,
    y_min = y_min,
    y_max = y_max
  )

  # RMSE calculation
  pred <- predict_hebart(newX = data, new_groups = groups, 
                         hebart_posterior = result, type = "mean")
  mse <- mean((pred - y)^2)
  rmse <- sqrt(mse)
  r.squared <- 1 - mse / stats::var(y)

  result$rmse <- rmse
  result$r.squared <- r.squared
  result$num_variables <- length(names_x)

  class(result) <- "hebart"

  return(result = result)
} # End main function