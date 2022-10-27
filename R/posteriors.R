# Tree conditional for HEBART
#' @name full_conditional_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Tree full conditional for hebart
#' @description A function that returns current tree full conditional
#' distribution value
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param num_trees The number of trees
#' @param tau The current value of tau
#' @param tau_phi The current value of tau_phi
#' @param tau_mu The current value of tau_mu
#' @param M The group matrix
full_conditional_hebart <- function(tree, R, num_trees, tau, tau_phi, tau_mu, M) {
  # Function to compute log full conditional distribution for an individual tree
  # R is a vector of partial residuals
  
  # hbeart version is
  # log_cond = sum (dmv(R, 0, Omega_R, log = TRUE) )
  # where M is the group allocation matrix.
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  log_cond <- 0
  for (i in 1:length(nj)) {
    M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
    R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
    Omega_R <- diag(nj[i])/tau + tcrossprod(M_j)/(num_trees * tau_phi) + 
      matrix(1/tau_mu, ncol = nj[i], nrow = nj[i])
    log_cond <- log_cond + mvnfast::dmvn(
      R_j, rep(0, nj[i]), Omega_R, log = TRUE
    )
  }
  
  return(log_cond)
}


#' @name get_tree_prior
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Get tree prior
#' @description Returns prior distribution value for a tree
#' @param tree The current tree
#' @param alpha The prior value for alpha (grow)
#' @param beta The prior value for beta (grow)

get_tree_prior <- function(tree, alpha, beta) {
  # Returns the tree log prior score
  
  # Need to work out the depth of the tree
  # First find the level of each node, then the depth is the maximum of the level
  level <- rep(NA, nrow(tree$tree_matrix))
  level[1] <- 0 # First row always level 0
  
  # Escpae quickly if tree is just a stump
  if (nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0
  }
  
  
  for (i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent <- tree$tree_matrix[i, "parent"]
    # This child must have a level one greater than it's current parent
    level[i] <- level[curr_parent] + 1
  }
  
  # Only compute for the internal nodes
  internal_nodes <- which(tree$tree_matrix[, "terminal"] == 0)
  log_prior <- 0
  for (i in 1:length(internal_nodes)) {
    log_prior <- log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]])
  }
  # Now add on terminal nodes
  terminal_nodes <- which(tree$tree_matrix[, "terminal"] == 1)
  for (i in 1:length(terminal_nodes)) {
    log_prior <- log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  }
  
  
  return(log_prior)
}




#' @name simulate_mu_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Simulate mu
#' @description Simulates mu for each terminal node
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param tau The  current value of tau
#' @param tau_phi The current value of tau_phi
#' @param tau_mu The value of tau_mu
#' @param M The value of M
#' @param num_trees The number of trees
simulate_mu_hebart <- function(tree, R, tau, tau_phi, tau_mu, M, num_trees) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal     <- which(tree$tree_matrix[, "terminal"] == 1)
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  for (i in 1:length(nj)) {
    M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
    R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
    Psi_R <- diag(nj[i])/tau + tcrossprod(M_j)/(num_trees*tau_phi)
    ones <- rep(1, nj[i])
    Prec_bit <- t(ones)%*%solve(Psi_R, ones) + tau_mu
    mean <- t(ones)%*%solve(Psi_R, R_j) / Prec_bit
    tree$tree_matrix[which_terminal[i], "mu"] <- stats::rnorm(1,
                                                              mean,
                                                              sd = 1/sqrt(Prec_bit))
  }
  tree$tree_matrix[which_non_terminal, "mu"] <- NA
  
  return(tree)
}

#' @name simulate_phi_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Simulate mu_js
#' @description Simulates mu_j for each group
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param groups The groups specification
#' @param tau The current value of tau
#' @param tau_phi The current value of tau_phi
#' @param M The group allocation matrix
#' @param num_trees The  number of trees
simulate_phi_hebart <- function(tree, R, groups, tau, tau_phi, M, num_trees) {
  
  # Simulate the group mu values for a given tree
  group_names <- unique(groups)
  num_groups <- length(unique(groups))
  group_col_names <- paste0("phi", group_names)

  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  #num_groups <- length(unique(groups))
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    curr_R <- R[tree$node_indices == which_terminal[i]]
    curr_mu <- tree$tree_matrix[which_terminal[i], "mu"]
    curr_M <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
    Prec_bit <- tau*crossprod(curr_M) + num_trees * tau_phi * diag(num_groups)
    mean <- solve(Prec_bit, tau * t(curr_M)%*%curr_R + rep(num_trees * tau_phi * curr_mu,
                                                           num_groups))
    tree$tree_matrix[which_terminal[i], 
                     sort(group_col_names)] <- mvnfast::rmvn(1,
                                                             mu = mean,
                                                             sigma = solve(Prec_bit))
    
  }
  
  tree$tree_matrix[which_non_terminal, sort(group_col_names)] <- NA
  
  return(tree)
}


#' # Simulate tau -----------------------------------------------
#' @name update_tau
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update tau
#' @description Samples values from the posterior distribution of tau
#' @param y The y vector
#' @param predictions The current set of predictions
#' @param nu The value of nu
#' @param lambda The value of lambda
update_tau <- function(y, predictions, nu, lambda) {
  
  # Sum of squared predictions
  S <- sum((y - predictions)^2)

  # Update
  tau <- stats::rgamma(1,
                       shape = (nu + length(y)) / 2,
                       rate = (S + nu * lambda) / 2
  )
  
  return(tau)
}

#' @name create_S
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Create S
#' @description Creates the S matrix 
#' @param curr_trees The current set of trees
#' @param groups if NULL just creates the terminal node matrix. If present 
#' allocates it into groups too
create_S <- function(curr_trees, groups = NULL){
  
  
  if(is.null(groups)) {
    node_fac <- factor(curr_trees[[1]]$node_indices)
    S <- mod.mat(node_fac)  
    if(length(curr_trees) > 1) {
      for(i in 2:length(curr_trees)) {
        node_fac <- factor(curr_trees[[i]]$node_indices)
        S <- cbind(S, mod.mat(node_fac))
      }
    }
  } else {
    # This is a hacky way of combining the terminal node and the groups
    # There must be a more elegant way
    node_groups <- factor(paste0(curr_trees[[1]]$node_indices, groups))
    S <- mod.mat(node_groups)
    if(length(curr_trees) > 1) {
      for(i in 2:length(curr_trees)) {
        node_groups <- factor(paste0(curr_trees[[i]]$node_indices, groups))
        S <- cbind(S, mod.mat(node_groups))
      }
    }
  }

  return(S)
}

#' @name update_sigma_phi
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update sigma_phi
#' @description Samples values from the posterior distribution of sigma_phi
#' @param y The response vector
#' @param S1 The group level node allocation matrix
#' @param S2 The terminal node allocation
#' @param tau The current value of tau
#' @param sigma_phi The current value of sigma_phi
#' @param tau_mu The current value of tau_mu
#' @param num_trees The number of trees
#' @param shape_sigma_phi Weibull shape parameter
#' @param scale_sigma_phi Weibull scale parameter
#' @param num_trees Number of trees
#' @param sigma_phi_sd Standard deviation of proposal distribution
update_sigma_phi <- function(y, S1, S2, sigma_phi, tau_mu, tau, 
                             shape_sigma_phi, scale_sigma_phi, num_trees, 
                             sigma_phi_sd = 0.1){
  
  repeat {
    # Proposal distribution
    new_sigma_phi <- sigma_phi + stats::rnorm(1, sd = sigma_phi_sd) 
    if (new_sigma_phi > 0)
      break
  }
  log_rat <- stats::pnorm(sigma_phi, sd = sigma_phi_sd, log = TRUE) - stats::pnorm(new_sigma_phi, sd = sigma_phi_sd, log = TRUE)
  new_tau_phi <- 1/(new_sigma_phi^2)
  tau_phi <- 1/(sigma_phi^2)
  
  n <- length(y)
  Omega_y_current <- diag(n)/tau + tcrossprod(S1)/(num_trees*tau_phi) + tcrossprod(S2)/tau_mu
  Omega_y_candidate <- diag(n)/tau + tcrossprod(S1)/(num_trees*new_tau_phi) + tcrossprod(S2)/tau_mu
  
  post_current <- mvnfast::dmvn(y, rep(0, n), sigma = Omega_y_current, log = TRUE)
  post_candidate <- mvnfast::dmvn(y, rep(0, n), sigma = Omega_y_candidate, log = TRUE)

  # Switching to gamma; 
  prior_current   <- stats::dgamma(sigma_phi,
                                    shape = shape_sigma_phi,
                                    scale = scale_sigma_phi, log = TRUE)
  prior_candidate <- stats::dgamma(new_sigma_phi,
                                    shape = shape_sigma_phi,
                                    scale = scale_sigma_phi, log = TRUE)
  
  log.alpha <- (post_candidate - post_current) + (prior_candidate - prior_current) + log_rat

  accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
  sigma_phi <- ifelse(accept, new_sigma_phi, sigma_phi)
  return(sigma_phi)
}
