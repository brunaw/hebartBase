#' @name tree_full_conditional
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Tree full conditional
#' @description A function that returns current tree full conditional 
#' distribution value
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param tau The current value of tau
#' @param tau_mu The current value of tau_mu
#' 
# Get complete conditions -------------------------------------------------

tree_full_conditional <- function(tree, R, tau, tau_mu) {
  # Function to compute log full conditional distirbution for an individual tree
  # R is a vector of partial residuals
  
  # Need to calculate log complete conditional, involves a sum over terminal nodes
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get sum of residuals and sum of residuals squared within each terminal node
  sumRsq_j <- stats::aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[, 2]
  S_j      <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
  
  # Now calculate the log posterior
  # log_post = 0.5 * length(R) * log(tau) + sum(0.5 * log( tau_mu / (tau_mu + nj * tau)) -
  #   0.5 * tau * (sumRsq_j - tau * S_j^2 / (tau_mu + nj * tau) ) )
  log_post <- 0.5 * length(R) * log(tau) +
    0.5 * (sum(log(tau_mu / (tau_mu + nj * tau))) -
             tau * sum(sumRsq_j) +
             tau^2 * sum(S_j^2 / (tau_mu + nj * tau)))
  return(log_post)
  #
  # New Mahdi version - slower
  # P1 = 0.5 * length(R) * log(tau)
  # P2 = 0.5 * sum( log( tau_mu / (tau_mu + nj * tau)))
  # P3 = -0.5 * tau * sum( sumRsq_j )
  # P4 = 0.5 * (tau^2) * sum ( (S_j^2) / (tau_mu + nj * tau) )
  #
  # return(P1 + P2 + P3 + P4)
}

# Tree conditional for HEBART
#' @name tree_full_conditional_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Tree full conditional for hebart
#' @description A function that returns current tree full conditional 
#' distribution value
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param k_1 The current value of k_1
#' @param k_2 The current value of k_2
#' @param M The group matrix
#' @param nu The current value of nu
#' @param lambda The current value of lambda
#' 
tree_full_conditional_hebart <- function(tree, R, k_1, k_2, M, nu, lambda) {
  # Function to compute log full conditional distribution for an individual tree
  # R is a vector of partial residuals
  
  # HEBART version is
  # log_cond = sum( log(gamma(n_j/2 + alpha)) - 0.5 * logdet(W) - (n_j/2 + alpha) *
  #    log(beta + 0.5 * t(R_j)%*%solve(W, R_j) )
  # where now W = k_2 * ones %*% t(ones) + k_1 * M %*% t(M) + diag(n_j)
  # where M is the group allocation matrix.
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  log_cond <- 0
  for (i in 1:length(nj)) {
    M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
    R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
    W <- k_2 * matrix(1, nrow = nj[i], ncol = nj[i]) + k_1 * M_j %*% t(M_j) + diag(nj[i])
    log_cond <- log_cond - 0.5 * logdet(W) + lgamma(nj[i] / 2 + nu / 2) - (nj[i] / 2 + nu / 2) *
      log(lambda / 2 + 0.5 * t(R_j) %*% solve(W, R_j))
    # There's also this term in the maths which I don't think is necessary
    # - 0.5 * nj[i] * log(2 * pi)
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
#' @param k_1 The  current value of k_1
#' @param k_2 The  current value of k_2 
#' 
# Simulate_mu -------------------------------------------------------------

simulate_mu_hebart <- function(tree, R, tau, k_1, k_2) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get sum of residuals in each terminal node
  sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
  
  # Now calculate mu values
  mu <- stats::rnorm(length(nj),
                     mean = (sumR / k_1) / (nj / k_1 + 1 / k_2),
                     sd = sqrt(1 / (tau * nj / k_1 + 1 / k_2))
  )
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, "mu"] <- NA
  
  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal, "mu"] <- mu
  
  return(tree)
}


#' @name simulate_mu_groups_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Simulate mu_js
#' @description Simulates mu_j for each group
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param groups The groups specification
#' @param tau The  current value of tau
#' @param k_1 The  current value of k_1
#' @param k_2 The  current value of k_2 
#' 
# Simulate mu groups HEBART -----------------------------------------------
simulate_mu_groups_hebart <- function(tree, R, groups, tau, k_1, k_2) {
  
  # Simulate the group mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  num_groups <- length(unique(groups))
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    curr_R           <- R[tree$node_indices == which_terminal[i]]
    curr_groups      <- groups[tree$node_indices == which_terminal[i]]
    curr_group_sizes <- table(curr_groups)
    group_R_means    <- stats::aggregate(curr_R, by = list(curr_groups), "sum")[, 2]
    curr_mu          <- tree$tree_matrix[which_terminal[i], "mu"]
    curr_group_mu     <- stats::rnorm(num_groups,
                                      mean = (curr_mu / k_1 + group_R_means) / (curr_group_sizes + 1 / k_1),
                                      sd = sqrt(1 / (curr_group_sizes + 1 / k_1))
    )
    tree$tree_matrix[which_terminal[i], paste0("mu", 1:num_groups)] <- curr_group_mu
    
  }
  
  # Wipe all the old mu groups out for other nodes
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  tree$tree_matrix[which_non_terminal, paste0("mu", 1:num_groups)] <- NA
  
  return(tree)
}



#' @name update_tau
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update tau
#' @description Samples values from the posterior distribution of tau
#' @param S The sum of squared residuals
#' @param nu .
#' @param lambda .
#' @param n The number of observations
#' 
# Update tau --------------------------------------------------------------
update_tau <- function(S, nu, lambda, n) {
  # Update from maths in Github folder
  tau <- stats::rgamma(1,
                       shape = (nu + n) / 2,
                       rate = (S + nu * lambda) / 2
  )
  # Alternative
  # tau = rgamma(1, shape = (nu + n) / 2 - 1, scale = 2 / (S + nu * lambda))
  
  return(tau)
}