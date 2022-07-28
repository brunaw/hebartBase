#' #' @name tree_full_conditional
#' #' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' #' @export
#' #' @title Tree full conditional
#' #' @description A function that returns current tree full conditional 
#' #' distribution value
#' #' @param tree The current tree
#' #' @param R The corresponding residuals for the tree
#' #' @param tau The current value of tau
#' #' @param tau_mu The current value of tau_mu
#' #' 
#' # Get complete conditions -------------------------------------------------
#' 
#' tree_full_conditional <- function(tree, R, tau, tau_mu) {
#'   # Function to compute log full conditional distirbution for an individual tree
#'   # R is a vector of partial residuals
#'   
#'   # Need to calculate log complete conditional, involves a sum over terminal nodes
#'   
#'   # First find which rows are terminal nodes
#'   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#'   
#'   # Get node sizes for each terminal node
#'   nj <- tree$tree_matrix[which_terminal, "node_size"]
#'   
#'   # Get sum of residuals and sum of residuals squared within each terminal node
#'   sumRsq_j <- stats::aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[, 2]
#'   S_j      <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
#'   
#'   # Now calculate the log posterior
#'   # log_post = 0.5 * length(R) * log(tau) + sum(0.5 * log( tau_mu / (tau_mu + nj * tau)) -
#'   #   0.5 * tau * (sumRsq_j - tau * S_j^2 / (tau_mu + nj * tau) ) )
#'   log_post <- 0.5 * length(R) * log(tau) +
#'     0.5 * (sum(log(tau_mu / (tau_mu + nj * tau))) -
#'              tau * sum(sumRsq_j) +
#'              tau^2 * sum(S_j^2 / (tau_mu + nj * tau)))
#'   return(log_post)
#'   #
#'   # New Mahdi version - slower
#'   # P1 = 0.5 * length(R) * log(tau)
#'   # P2 = 0.5 * sum( log( tau_mu / (tau_mu + nj * tau)))
#'   # P3 = -0.5 * tau * sum( sumRsq_j )
#'   # P4 = 0.5 * (tau^2) * sum ( (S_j^2) / (tau_mu + nj * tau) )
#'   #
#'   # return(P1 + P2 + P3 + P4)
#' }

#' #' @name simulate_mu_hebart_2
#' #' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' #' @export
#' #' @title Simulate mu
#' #' @description Simulates mu for each terminal node
#' #' @param tree The current tree
#' #' @param R The corresponding residuals for the tree
#' #' @param tau The  current value of tau
#' #' @param k_1 The  current value of k_1
#' #' @param k_2 The  current value of k_2
#' #' @param groups The groups 
#' #' @param type Action type
#' #' @param acc Acceptance or not of proposal
#' #' 
#' simulate_mu_hebart_2 <- function(tree, R, tau, k_1, k_2, groups, type,
#'                                  acc) {
#'   
#'   # psi  <- k_1 * M %*% t(M) + diag(n)
#'   # mean <- (rep(1, n) %*% solve(psi, R)) / (rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))
#'   # var  <- 1/((rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))*tau)
#'   # Simulate mu values for a given tree
#'   
#'   # First find which rows are terminal nodes
#'   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#'   
#'   # Get node sizes for each terminal node
#'   nj <- tree$tree_matrix[which_terminal, "node_size"]
#'   
#'   group_names     <- unique(groups)
#'   num_groups      <- length(unique(groups))
#'   group_col_names <- paste0("mu", group_names)
#'   #df_groups       <- data.frame(groups = group_names)
#'   
#'   # mu_js 
#'   mu_js <- tree$tree_matrix[, c("terminal", sort(group_col_names))]
#'   if(!"matrix" %in% class(mu_js)){
#'     mu_js <- matrix(mu_js, byrow = FALSE, ncol = length(c("terminal", sort(group_col_names))))
#'     colnames(mu_js) <- c("terminal", sort(group_col_names))
#'   }
#'   n_init <- nrow(stats::na.omit(mu_js))
#'   
#'   #if(nrow(mu_js) == 1){type = "same"}
#'   
#'   # dealing with the first iteration
#'   # if(nrow(tree$tree_matrix)==3){
#'   #   curr_sum_mu <- c(0, 0)
#'   #   correct_inds <- c(1, 2)
#'   # } else {
#'   # existing node indices:
#'   if(type == "grow"){
#'     if(nrow(mu_js) > 1){
#'       inds            <- unique(tree$node_indices)
#'       mu_js           <- cbind(mu_js, node_index = 1:max(inds))
#'       non_na_mu_js    <- stats::na.omit(mu_js)
#'       which_ind       <- non_na_mu_js[, "node_index"]
#'       which_to_dup    <- which_ind[which(!which_ind %in% inds)]
#'       non_dup         <- which_ind[which(which_ind %in% inds)]
#'       new_inds        <- inds[which(!inds %in% non_dup)]
#'       
#'       mu_js          <- rbind(mu_js, mu_js[which_to_dup, ])
#'       mu_js          <- stats::na.omit(mu_js)
#'       mu_js[mu_js[, "node_index"] == which_to_dup, "node_index"] <- new_inds
#'       correct_inds   <- mu_js[, "node_index"]
#'       curr_sum_mu    <- rowSums(mu_js[, sort(group_col_names)])   
#'     
#'       
#'     } else{
#'       curr_sum_mu <- sum(stats::na.omit(mu_js[, sort(group_col_names)]))
#'       correct_inds <- 2
#'     }
#'   } else if(type == "prune"){
#'     if(nrow(mu_js) > 1){
#'       # For prune, 
#'       inds             <- unique(tree$node_indices)
#'       mu_js            <- cbind(mu_js, node_index = 1:max(inds))
#'       mu_js            <- stats::na.omit(mu_js)
#'       correct_inds     <- mu_js[, "node_index"]
#'       curr_sum_mu      <- rowSums(mu_js[, sort(group_col_names)])              
#'       
#'     } else {
#'       curr_sum_mu <- sum(stats::na.omit(mu_js[, sort(group_col_names)]))
#'       correct_inds <- 1
#'     }
#'   }
#'   #}
#'   
#'   if(nrow(tree$tree_matrix)==3 & n_init ==1){
#'     curr_sum_mu <- c(0, 0)
#'     correct_inds <- c(3, 2)
#'   }
#' 
#'   mu <- stats::rnorm(length(nj),
#'                      mean = (curr_sum_mu / k_1) / (num_groups / k_1 + 1 / k_2),
#'                      sd = sqrt(1 / (tau * (num_groups / k_1 + 1 / k_2)))
#'   )
#'   
#'   # Get sum of residuals in each terminal node
#'   #sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
#'   
#'   # Wipe all the old mus out for other nodes
#'   tree$tree_matrix[, "mu"] <- NA
#'   
#'   # Put in just the ones that are useful
#'   #tree$tree_matrix[which_terminal, "mu"] <- mu
#'   # a bug correction for when a grow was proposed but not accepted 
#'   if(is.na(acc)){ acc <- FALSE }
#'   if(length(curr_sum_mu) == 1 & !acc){
#'     correct_inds <- 1
#'   }
#'   tree$tree_matrix[correct_inds, "mu"] <- mu
#' 
#'   
#'   return(tree)
#' }
#' 
#' #' @name simulate_mu_hebart_2
#' #' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' #' @export
#' #' @title Simulate mu
#' #' @description Simulates mu for each terminal node
#' #' @param tree The current tree
#' #' @param R The corresponding residuals for the tree
#' #' @param tau The  current value of tau
#' #' @param k_1 The  current value of k_1
#' #' @param k_2 The  current value of k_2
#' #' @param groups The groups 
#' #' @param type Action type
#' #' @param acc Acceptance or not of proposal
#' #' 
#' simulate_mu_hebart_3 <- function(tree, R, tau, k_1, k_2, groups, type,
#'                                  acc) {
#'   
#'   # psi  <- k_1 * M %*% t(M) + diag(n)
#'   # mean <- (rep(1, n) %*% solve(psi, R)) / (rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))
#'   # var  <- 1/((rep(1, n) %*% solve(psi, rep(1, n)) + (1/k_2))*tau)
#'   # Simulate mu values for a given tree
#'   
#'   # First find which rows are terminal nodes
#'   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#'   
#'   # Get node sizes for each terminal node
#'   nj <- tree$tree_matrix[which_terminal, "node_size"]
#'   
#'   group_names     <- unique(groups)
#'   num_groups      <- length(unique(groups))
#'   group_col_names <- paste0("mu", group_names)
#'   #df_groups       <- data.frame(groups = group_names)
#'   
#'   curr_mus <- curr_var <- c()
#'   for(i in 1:length(nj)){
#'     n   <-  nj[i]
#'     node_groups <- groups[tree$node_indices == which_terminal[i]]
#'     node_R      <- R[tree$node_indices == which_terminal[i]]
#'     M_node <- stats::model.matrix(~ factor(node_groups) - 1)
#'     PSI <- (k_1 * M_node %*% t(M_node)) + diag(n)
#'     curr_mus[i] <- rep(1, n) %*% solve(PSI, node_R) / (rep(1, n) %*% solve(PSI, rep(1, n)) + 1/k_2)
#'     curr_var[i] <- tau * (rep(1, n) %*% solve(PSI, rep(1, n)) + 1/k_2)
#'     }
#'   
#'   
#'   
#'   #mu_js <- tree$tree_matrix[, c("terminal", sort(group_col_names))]
#'   
#'   mu <- stats::rnorm(length(nj),
#'                      mean = curr_mus,
#'                      sd = sqrt(1 / curr_var)
#'   )
#'   
#'   # Get sum of residuals in each terminal node
#'   #sumR <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
#'   
#'   # Wipe all the old mus out for other nodes
#'   tree$tree_matrix[, "mu"] <- NA
#'   tree$tree_matrix[which_terminal, "mu"] <- mu
#'   
#'   return(tree)
#' }
#' 
