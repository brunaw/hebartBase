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
#' @param num_groups The number of groups 
#' 
# Simulate_mu -------------------------------------------------------------

simulate_mu_hebart <- function(tree, R, tau, k_1, k_2, num_groups) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get sum of residuals in each terminal node
  sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
  
  
  mu <- stats::rnorm(length(nj),
              mean = (sumR / k_1) / (nj / k_1 + 1 / k_2),
              sd = sqrt(1 / (tau * nj / k_1 + 1 / k_2))
  )
  
  # Now calculate mu values
  # mu <- stats::rnorm(length(nj),
  #                    mean = (sumR / k_1) / (nj / k_1 + 1 / k_2),
  #                    sd = sqrt(1 / (tau * (nj / k_1 + 1 / k_2)))
  # )
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, "mu"] <- NA
  
  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal, "mu"] <- mu
  
  return(tree)
}

#' @name simulate_mu_hebart_2
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Simulate mu
#' @description Simulates mu for each terminal node
#' @param tree The current tree
#' @param R The corresponding residuals for the tree
#' @param tau The  current value of tau
#' @param k_1 The  current value of k_1
#' @param k_2 The  current value of k_2
#' @param groups The groups 
#' @param type Action type
#' @param acc Acceptance or not of proposal
#' 
simulate_mu_hebart_2 <- function(tree, R, tau, k_1, k_2, groups, type,
                                 acc) {
  
  # psi  <- k_1 * M %*% t(M) + diag(n)
  # mean <- (rep(1, n) %*% solve(psi, R)) / (rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))
  # var  <- 1/((rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))*tau)
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  group_names     <- unique(groups)
  num_groups      <- length(unique(groups))
  group_col_names <- paste0("mu", group_names)
  #df_groups       <- data.frame(groups = group_names)
  
  # mu_js 
  mu_js <- tree$tree_matrix[, c("terminal", sort(group_col_names))]
  if(!"matrix" %in% class(mu_js)){
    mu_js <- matrix(mu_js, byrow = FALSE, ncol = length(c("terminal", sort(group_col_names))))
    colnames(mu_js) <- c("terminal", sort(group_col_names))
  }
  n_init <- nrow(stats::na.omit(mu_js))
  
  #if(nrow(mu_js) == 1){type = "same"}
  
  # dealing with the first iteration
  # if(nrow(tree$tree_matrix)==3){
  #   curr_sum_mu <- c(0, 0)
  #   correct_inds <- c(1, 2)
  # } else {
  # existing node indices:
  if(type == "grow"){
    if(nrow(mu_js) > 1){
      inds            <- unique(tree$node_indices)
      mu_js           <- cbind(mu_js, node_index = 1:max(inds))
      non_na_mu_js    <- stats::na.omit(mu_js)
      which_ind       <- non_na_mu_js[, "node_index"]
      which_to_dup    <- which_ind[which(!which_ind %in% inds)]
      non_dup         <- which_ind[which(which_ind %in% inds)]
      new_inds        <- inds[which(!inds %in% non_dup)]
      
      mu_js          <- rbind(mu_js, mu_js[which_to_dup, ])
      mu_js          <- stats::na.omit(mu_js)
      mu_js[mu_js[, "node_index"] == which_to_dup, "node_index"] <- new_inds
      correct_inds   <- mu_js[, "node_index"]
      curr_sum_mu    <- rowSums(mu_js[, sort(group_col_names)])   
    
      
    } else{
      curr_sum_mu <- sum(stats::na.omit(mu_js[, sort(group_col_names)]))
      correct_inds <- 2
    }
  } else if(type == "prune"){
    if(nrow(mu_js) > 1){
      # For prune, 
      inds             <- unique(tree$node_indices)
      mu_js            <- cbind(mu_js, node_index = 1:max(inds))
      mu_js            <- stats::na.omit(mu_js)
      correct_inds     <- mu_js[, "node_index"]
      curr_sum_mu      <- rowSums(mu_js[, sort(group_col_names)])              
      
    } else {
      curr_sum_mu <- sum(stats::na.omit(mu_js[, sort(group_col_names)]))
      correct_inds <- 1
    }
  }
  #}
  
  if(nrow(tree$tree_matrix)==3 & n_init ==1){
    curr_sum_mu <- c(0, 0)
    correct_inds <- c(3, 2)
  }

  mu <- stats::rnorm(length(nj),
                     mean = (curr_sum_mu / k_1) / (num_groups / k_1 + 1 / k_2),
                     sd = sqrt(1 / (tau * (num_groups / k_1 + 1 / k_2)))
  )
  
  # Get sum of residuals in each terminal node
  #sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[, "mu"] <- NA
  
  # Put in just the ones that are useful
  #tree$tree_matrix[which_terminal, "mu"] <- mu
  # a bug correction for when a grow was proposed but not accepted 
  if(is.na(acc)){ acc <- FALSE }
  if(length(curr_sum_mu) == 1 & !acc){
    correct_inds <- 1
  }
  tree$tree_matrix[correct_inds, "mu"] <- mu

  
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
  
  group_names     <- unique(groups)
  num_groups      <- length(unique(groups))
  group_col_names <- paste0("mu", group_names)
  df_groups       <- data.frame(groups = group_names)
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    # curr_R           <- R[tree$node_indices == which_terminal[i]]
    # curr_groups      <- groups[tree$node_indices == which_terminal[i]]
    # curr_group_sizes <- data.frame(groups = curr_groups) |> 
    #   dplyr::count(groups)
    # means            <- data.frame(curr_R = curr_R, 
    #                                groups = curr_groups)
    # # Correcting for missing groups in the node
    # means            <- dplyr::group_by(means, groups) |>
    #   # sum, not mean
    #   dplyr::summarise(mean_curr_R = sum(curr_R)) |> 
    #   dplyr::left_join(curr_group_sizes, by = "groups") |> 
    #   dplyr::right_join(df_groups, by = "groups") 
    # means <- stats::na.omit(means)
    # group_R_means     <- means$mean_curr_R
    # curr_mu           <- tree$tree_matrix[which_terminal[i], "mu"]
    # curr_group_sizes  <- means$n
    # curr_group_mu     <- stats::rnorm(
    #   num_groups,
    #   mean = (curr_mu / k_1 + group_R_means) / (curr_group_sizes + 1 / k_1),
    #   sd = sqrt(1 / (tau * (curr_group_sizes + 1 / k_1)))
    # )
    
    curr_R           <- R[tree$node_indices == which_terminal[i]]
    curr_groups      <- groups[tree$node_indices == which_terminal[i]]
    curr_group_sizes <- table(curr_groups)
    curr_group_sizes <- curr_group_sizes[curr_group_sizes>0]
    group_R_means    <- stats::aggregate(curr_R, by = list(curr_groups), "sum")[, 2]
    curr_mu          <- tree$tree_matrix[which_terminal[i], "mu"]
    curr_group_mu    <- stats::rnorm(num_groups,
                           mean = (curr_mu / k_1 + group_R_means) / (curr_group_sizes + 1 / k_1),
                           sd = sqrt(1 / (curr_group_sizes + 1 / k_1))
    )
    tree$tree_matrix[which_terminal[i], sort(group_col_names)] <- curr_group_mu
    
  }
  
  # Wipe all the old mu groups out for other nodes
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  tree$tree_matrix[which_non_terminal, sort(group_col_names)] <- NA
  
  return(tree)
}



#' @name update_tau
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Update tau
#' @description Samples values from the posterior distribution of tau
#' @param S The sum of squared residuals
#' @param res_mu_j The mu_j residuals
#' @param res_mu The mu residuals
#' @param nu The current value of nu
#' @param lambda The current value of lambda
#' @param n The number of observations
#' @param groups The grouping variable
#' @param k_1 The current value of k1
#' @param k_2 The current value of k2
# Update tau --------------------------------------------------------------
update_tau <- function(S, res_mu_j, res_mu, nu, lambda, n, groups, k_1, k_2) {
  
  num_groups      <- length(unique(groups))

  # Simple version: 
  tau <- stats::rgamma(1,
                       shape = (nu + n) / 2,
                       rate = (S + nu * lambda) / 2
  )
  
  # Update from maths in Github folder
  # tau <- stats::rgamma(1,
  #                      shape = (nu + n + num_groups + 1) / 2,
  #                      rate = (S + (nu * lambda) + res_mu_j/k_1 + res_mu/k_2) / 2
  # )
  #Alternative
  #tau = rgamma(1, shape = (nu + n) / 2 - 1, scale = 2 / (S + nu * lambda))
  
  return(tau)
}


#' @name full_conditional_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Full conditional of all the data 
#' @description A function that returns the full conditional value
#' @param y The y vector
#' @param k_1 The current value of k_1
#' @param k_2 The current value of k_2
#' @param M The group matrix
#' @param nu The current value of nu
#' @param lambda The current value of lambda
#' 
full_conditional_hebart <- function(y, k_1, k_2, M, nu, lambda) {
  n = length(y)
  W_1 <- (k_2 * matrix(1, nrow = n, ncol = n)) + (k_1 * M %*% t(M)) + diag(n)
  # lgamma(n / 2 + nu / 2) = will be the same for either k_1
  log_cond <- - 0.5 * logdet(W_1) - ((n / 2 + nu / 2) *
    log(lambda / 2 + 0.5 * t(y) %*% solve(W_1, y)))
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
#' @param M The group matrix
#' @param nu The current value of nu
#' @param lambda The current value of lambda
#' @param prior Logical to decide whether to use a prior (set as 
#' dweibull(1, 5) for now)
#' 
# Update k1 --------------------------------------------------------------

update_k1 <- function(y, min_u, max_u, k_1, k_2, M, nu, lambda, prior = TRUE){
  new_k1    <- stats::runif(1, min = min_u, max = max_u)
  current   <- full_conditional_hebart(y, k_1, k_2, M, nu, lambda)
  candidate <- full_conditional_hebart(y, new_k1, k_2, M, nu, lambda)
  
  if(prior){
    prior_current   <- stats::dweibull(k_1, shape = 1, 5, log = TRUE)
    prior_candidate <- stats::dweibull(new_k1, shape = 1, 5, log = TRUE)
    log.alpha <- (candidate - current) + (prior_candidate - prior_current) # ll is the log likelihood, lprior the log prior
  } else {
    log.alpha <- candidate - current
  }
  
  accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
  theta <- ifelse(accept, new_k1, k_1)
  return(theta)
}
