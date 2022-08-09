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
tree_full_conditional_hebart <- function(tree, R, k_1, k_2, M, nu, lambda) {
  # Function to compute log full conditional distribution for an individual tree
  # R is a vector of partial residuals
  
  # hbeart version is
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
    W <- k_2 + k_1 * tcrossprod(M_j) + diag(nj[i])
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
#' @param num_trees The  number of trees
#'
# Simulate_mu -------------------------------------------------------------
simulate_mu_hebart <- function(tree, R, tau, k_1, k_2, num_trees) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get sum of residuals in each terminal node
  sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
  
  # Now calculate mu values
  mu <- stats::rnorm(length(nj),
                     mean = (sumR / (k_1/num_trees)) / (nj / (k_1/num_trees) + 1 / k_2),
                     sd = sqrt(1 / (tau * (nj / (k_1/num_trees) + 1 / k_2)))
  )
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, "mu"] <- NA
  
  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal, "mu"] <- mu
  
  return(tree)
}

simulate_mu_hebart2 <- function(tree, R, M, tau, k_1, k_2) {
  
  # Simulate mu values for a given tree
  
  # this is the marginalised mu version from
  # https://bookdown.org/connect/#/apps/56e67516-559f-4d69-91df-54702fbc2206/access
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, "mu"] <- NA
  
  # Loop through terminal nodes to get values
  for (i in 1:length(nj)) {
    M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
    R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
    Psi <- k_1 * tcrossprod(M_j) + diag(nj[i])
    #Psi <- k_1 * tcrossprod(M_j) 
    ones <- rep(1, nj[i])
    Prec_bit <- t(ones)%*%solve(Psi, ones) + 1/k_2
    mean <- t(ones)%*%Psi%*%R_j / Prec_bit
    tree$tree_matrix[which_terminal[i], "mu"] <- stats::rnorm(1,
                                                              mean,
                                                              sd = 1/sqrt(tau*Prec_bit))
  }
  
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
#' @param num_trees The  number of trees
simulate_mu_groups_hebart <- function(tree, R, groups, tau, k_1, k_2, num_trees) {
  
  # Simulate the group mu values for a given tree
  
    group_names     <- unique(groups)
    num_groups      <- length(unique(groups))
    group_col_names <- paste0("mu", group_names)

  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  #num_groups <- length(unique(groups))
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    
    curr_R           <- R[tree$node_indices == which_terminal[i]]
    curr_groups      <- groups[tree$node_indices == which_terminal[i]]
    curr_group_sizes <- table(curr_groups)
    group_R_means    <- stats::aggregate(curr_R, by = list(curr_groups), "sum")[, 2]
    curr_mu          <- tree$tree_matrix[which_terminal[i], "mu"]
    
    # For the missing mus
    if(length(group_R_means) < num_groups | sum(curr_group_sizes == 0) > 0){
      group_names <- as.character(group_names)
      df_groups <- data.frame(groups = group_names)
      curr_groups <- as.character(curr_groups)
      
      curr_group_sizes <- curr_group_sizes[curr_group_sizes > 0]
      df_actual <- data.frame(groups = unique(curr_groups), 
                              mus = group_R_means, 
                              n = as.vector(curr_group_sizes))
     df <- df_groups |> 
       dplyr::left_join(df_actual, by = "groups") |> 
       dplyr::mutate(mus = ifelse(is.na(mus), mean(curr_R), mus),
                     n = ifelse(is.na(n), 0, n))
     group_R_means <- df$mus
     curr_group_sizes <- df$n
    }
    
    curr_group_mu <- stats::rnorm(num_groups,
                                  mean = (curr_mu / (k_1/num_trees) + group_R_means) / (curr_group_sizes + 1 / (k_1/num_trees)),
                                  sd = sqrt(1 / (curr_group_sizes + 1 / (k_1/num_trees)))
    )

    #tree$tree_matrix[which_terminal[i], paste0("mu", 1:num_groups)] <- curr_group_mu
    tree$tree_matrix[which_terminal[i], sort(group_col_names)] <- curr_group_mu
    
  }
  
  # Wipe all the old mu groups out for other nodes
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  tree$tree_matrix[which_non_terminal,  sort(group_col_names)] <- NA
  
  return(tree)
}


#' # Simulate mu groups HEBART -----------------------------------------------
#' @name update_tau
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update tau
#' @description Samples values from the posterior distribution of tau
#' @param y The y vector
#' @param nu The current value of nu
#' @param M The current group model matrix
#' @param lambda The current value of lambda
#' @param num_groups The number of groups
#' @param k_1 The current value of k1
#' @param k_2 The current value of k2
#' # Update tau --------------------------------------------------------------
update_tau <- function(y, M, nu, lambda, num_groups, k_1, k_2, num_trees, last_trees) {
  
  W_tilde <- create_S(k_1, k_2, last_trees, num_trees)$W_tilde

  n <- length(y)
  #W_1 <- (k_2 * matrix(1, nrow = n, ncol = n)) + (k_1 * M %*% t(M)) + diag(n)
  S <- t(y) %*% solve(W_tilde, y)

  # New update
  tau <- stats::rgamma(1,
                       shape = (nu + n) / 2,
                       rate = (S + nu * lambda) / 2
  )
  
  return(tau)
}

#' @name create_S
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Create S
#' @description Creates the S matrix 
#' @param k_1 The current value of k1
#' @param k_2 The current value of k2
#' @param last_trees The last trees 
#' @param num_trees The number of trees

create_S <- function(k_1, k_2, last_trees, num_trees){
  # First tree
  tree <-   last_trees[[1]]
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  if(length(unique(tree$node_indices)) > 1){
    S <- stats::model.matrix(~factor(tree$node_indices) - 1)
  } else {
    S <- matrix(rep(1, length(tree$node_indices)))
  }
  n_nodes <- length(unique(tree$node_indices))
  
  if(num_trees > 1){
    for(i in 2:num_trees){
      tree <-   last_trees[[i]]
      which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
      if(length(unique(tree$node_indices)) > 1){
        Ss <- stats::model.matrix(~factor(tree$node_indices) - 1)
      } else {
        Ss <- matrix(rep(1, length(tree$node_indices)))
      }
  
      S <- cbind(S, Ss)
      n_nodes <- length(unique(tree$node_indices)) + n_nodes
    } 
  }
  
  W_tilde <- diag(nrow(S)) + ((k_1)/num_trees)*(S %*% t(S)) + ((k_2)/num_trees)*(S %*% t(S))
  
  return(list(S = S, n_nodes = n_nodes, W_tilde = W_tilde))
}

 
#' @name full_conditional_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Full conditional of all the data
#' @description A function that returns the full conditional value
#' @param y The y vector
#' @param k_1 The current value of k_1
#' @param k_2 The current value of k_2
#' @param last_trees The last set of trees
#' @param num_trees The number of trees
#' @param tau The current value of tau
#'
full_conditional_hebart <- function(y, k_1, k_2, last_trees, num_trees, tau) {
  # n = length(y)
  # W_1 <- (k_2 * matrix(1, nrow = n, ncol = n)) + (k_1 * M %*% t(M)) + diag(n)
  # # lgamma(n / 2 + nu / 2) = will be the same for either k_1
  # log_cond <- - 0.5 * logdet(W_1) - ((n / 2 + nu / 2) *
  #   log(lambda / 2 + 0.5 * t(y) %*% solve(W_1, y)))
  
  W_tilde <- create_S(k_1, k_2, last_trees, num_trees)$W_tilde
  sig_y    <- (1/tau ) * W_tilde
  #log_cond <-  sum(stats::dnorm(y, mean = 0, sd = diag(sd_y), log = TRUE))
  # log_cond <- mvtnorm::dmvnorm(y, mean = rep(0, length(y)), 
  #                              sigma = sig_y, log = TRUE)
  log_cond <- mvnfast::dmvn(y, mu = rep(0, length(y)), 
                               sigma = sig_y, log = TRUE)
  #log_cond <- log(log_cond + 0.00000001)
  return(log_cond)
}

#' @name update_k1
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update k1
#' @description Samples values from the posterior distribution of k1
#' @param y The y vector
#' @param min_u The lower bound of the sampled values of k1
#' @param max_u The upper bound of the sampled values of k1
#' @param k_1 The current value of k_1
#' @param k_2 The current value of k_2
#' @param last_trees The last set of trees
#' @param num_trees The number of trees
#' @param tau The current value of tau
#' @param prior Logical to decide whether to use a prior (set as
#' dweibull(1, 5) for now)
#'
# Update k1 --------------------------------------------------------------

update_k1 <- function(y, min_u, max_u, k_1, k_2, last_trees, num_trees, tau, prior = TRUE){
  # new_k1    <- stats::runif(1, min = min_u, max = max_u)
  k1_sd <- 0.1
  repeat {
    # Proposal distribution
    new_k1 <- k_1 + stats::rnorm(1, sd = k1_sd) 
    if (new_k1 > 0)
      break
  }
  log_rat <- stats::pnorm(k_1, sd = k1_sd, log = TRUE) - stats::pnorm(new_k1, sd = k1_sd, log = TRUE)
  
  current   <- full_conditional_hebart(y, k_1, k_2, last_trees, num_trees, tau)
  candidate <- full_conditional_hebart(y, new_k1, k_2, last_trees, num_trees, tau)

  if(prior){
    prior_current   <- stats::dweibull(k_1, shape = 1, 5, log = TRUE)
    prior_candidate <- stats::dweibull(new_k1, shape = 1, 5, log = TRUE)
    log.alpha <- (candidate - current) + (prior_candidate - prior_current) + log_rat # ll is the log likelihood, lprior the log prior
  } else {
    log.alpha <- candidate - current + log_rat
  }

  accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
  theta <- ifelse(accept, new_k1, k_1)
  return(theta)
}
