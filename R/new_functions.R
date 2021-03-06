# Function to create stump ------------------------------------------------

# create_stump <- function(num_trees,
#                          groups,
#                          y,
#                          X) {
#   
#   num_groups      <- length(unique(groups))
#   # Each tree has 8+num_groups columns and 2 elements
#   # The 3 elements are the tree matrix, and the node indices
#   # The tree matrix has columns:
#   # Terminal (0 = no, 1 = yes)
#   # Child left
#   # Child right
#   # Node parents
#   # Split variable
#   # Split value
#   # mu values
#   # mu values for each group
#   # Node size
#   
#   # Create holder for trees
#   all_trees <- vector("list", length = num_trees)
#   # Loop through trees
#   for (j in 1:num_trees) {
#     # Set up each tree to have two elements in the list as described above
#     all_trees[[j]] <- vector("list", length = 2)
#     # Give the elements names
#     names(all_trees[[j]]) <- c(
#       "tree_matrix",
#       "node_indices"
#     )
#     # Create the two elements: first is a matrix
#     all_trees[[j]][[1]] <- matrix(NA, ncol = 8 + num_groups, nrow = 1)
#     
#     # Second is the assignment to node indices
#     all_trees[[j]][[2]] <- rep(1, length(y))
#     
#     # Create column names
#     colnames(all_trees[[j]][[1]]) <- c(
#       "terminal",
#       "child_left",
#       "child_right",
#       "parent",
#       "split_variable",
#       "split_value",
#       "mu",
#       paste0("mu", 1:num_groups),
#       "node_size"
#     )
#     
#     # Set values for stump
#     all_trees[[j]][[1]][1, ] <- c(1, 1, NA, NA, NA, NA, rep(0, num_groups + 1), length(y))
#   } # End of loop through trees
#   
#   return(all_trees)
# } # End of function

# Function to update trees ------------------------------------------------
# i.e. grow prune change swap

# update_tree <- function(y, # Target variable
#                         X, # Feature matrix
#                         groups, # Number of groups
#                         type = c(
#                           "grow", # Grow existing tree
#                           "prune", # Prune existing tree
#                           "change", # Change existing tree - change split variable and value for an internal node
#                           "swap"
#                         ), # Swap existing tree - swap a parent/child combo where both are internal
#                         curr_tree, # The current set of trees (not required if type is stump)
#                         node_min_size) { # The minimum size of a node to grow
#   num_groups      <- length(unique(groups))
#   # Each tree has 8 + num_groups columns and 2 elements
#   # The 3 elements are the tree matrix, and the node indices
#   # The tree matrix has columns:
#   # Terminal (0 = no, 1 = yes)
#   # Child left
#   # Child right
#   # Node parents
#   # Split variable
#   # Split value
#   # mu values
#   # mu values for each group
#   # Node size
#   
#   # Call the appropriate function to get the new tree
#   new_tree <- switch(type,
#                      grow = grow_tree(X, y, num_groups, curr_tree, node_min_size),
#                      prune = prune_tree(X, y, curr_tree),
#                      change = change_tree(X, y, curr_tree, node_min_size),
#                      swap = swap_tree(X, y, curr_tree, node_min_size)
#   )
#   
#   # Return the new tree
#   return(new_tree)
# } # End of update_tree function


# Grow_tree function ------------------------------------------------------

# grow_tree <- function(X, y, num_groups, curr_tree, node_min_size) {
#   
#   # Set up holder for new tree
#   new_tree <- curr_tree
#   
#   # Get the list of terminal nodes
#   terminal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 1) # Create the list of terminal nodes
#   
#   # Find terminal node sizes
#   terminal_node_size <- new_tree$tree_matrix[terminal_nodes, "node_size"]
#   
#   # Add two extra rows to the tree in question
#   new_tree$tree_matrix <- rbind(
#     new_tree$tree_matrix,
#     c(1, NA, NA, NA, NA, NA, rep(NA, num_groups + 1), NA), # Make sure they're both terminal
#     c(1, NA, NA, NA, NA, NA, rep(NA, num_groups + 1), NA)
#   )
#   
#   # Choose a random terminal node to split
#   node_to_split <- sample(terminal_nodes, 1,
#                           prob = as.integer(terminal_node_size > node_min_size)
#   ) # Choose which node to split, set prob to zero for any nodes that are too small
#   
#   # Choose a split variable uniformly from all columns
#   split_variable <- sample(1:ncol(X), 1)
#   # Choose a split value from the range of the current node but stop it from choosing empty nodex
#   # low_bound = min(X[new_tree$node_indices == node_to_split,
#   #                   split_variable]) + .Machine$double.eps
#   # high_bound = max(X[new_tree$node_indices == node_to_split,
#   #                   split_variable]) - .Machine$double.eps
#   # split_value = runif(1, low_bound, high_bound)
#   
#   # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
#   available_values <- sort(unique(X[
#     new_tree$node_indices == node_to_split,
#     split_variable
#   ]))
#   split_value <- sample(available_values[-c(1, length(available_values))], 1)
#   
#   curr_parent <- new_tree$tree_matrix[node_to_split, "parent"] # Make sure to keep the current parent in there. Will be NA if at the root node
#   new_tree$tree_matrix[node_to_split, 1:6] <- c(
#     0, # Now not temrinal
#     nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
#     nrow(new_tree$tree_matrix), # child_right is penultimate row
#     curr_parent,
#     split_variable,
#     split_value
#   )
#   
#   #  Fill in the parents of these two nodes
#   new_tree$tree_matrix[nrow(new_tree$tree_matrix), "parent"] <- node_to_split
#   new_tree$tree_matrix[nrow(new_tree$tree_matrix) - 1, "parent"] <- node_to_split
#   
#   # Now call the fill function on this tree
#   new_tree <- fill_tree_details(new_tree, X)
#   
#   # Reject this tree is any values are smaller the node_min_size
#   if (any(new_tree$tree_matrix[, "node_size"] < node_min_size)) new_tree <- curr_tree
#   
#   # Return new_tree
#   return(new_tree)
# } # End of grow_tree function


# Prune_tree function -----------------------------------------------------

# prune_tree <- function(X, y, curr_tree) {
#   
#   # Create placeholder for new tree
#   new_tree <- curr_tree
#   
#   if (nrow(new_tree$tree_matrix) == 1) {
#     return(new_tree)
#   } # No point in pruning a stump!
#   
#   # Get the list of terminal nodes
#   terminal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 1) # Create the list of terminal nodes
#   
#   # Pick a random termianl node to prune
#   # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
#   bad_node_to_prune <- TRUE # Assume a bad node pick
#   while (bad_node_to_prune) {
#     
#     # Choose a random terminal node
#     node_to_prune <- sample(terminal_nodes, 1)
#     
#     # Find the parent of this terminal node
#     parent_pick <- new_tree$tree_matrix[node_to_prune, "parent"]
#     
#     # Get the two children of this parent
#     child_left <- new_tree$tree_matrix[parent_pick, "child_left"]
#     child_right <- new_tree$tree_matrix[parent_pick, "child_right"]
#     
#     # See whether either are terminal
#     child_left_terminal <- new_tree$tree_matrix[child_left, "terminal"]
#     child_right_terminal <- new_tree$tree_matrix[child_right, "terminal"]
#     
#     # If both are terminal then great
#     if ((child_left_terminal == 1) & (child_right_terminal == 1)) {
#       bad_node_to_prune <- FALSE # Have chosen a pair of terminal nodes so exist while loop
#     }
#   } # End of bad node to prune while loop
#   
#   # Delete these two rows from the tree matrix
#   new_tree$tree_matrix <- new_tree$tree_matrix[-c(child_left, child_right), ,
#                                                drop = FALSE
#   ]
#   # Make this node terminal again with no children or split values
#   new_tree$tree_matrix[parent_pick, c(
#     "terminal",
#     "child_left",
#     "child_right",
#     "split_variable",
#     "split_value"
#   )] <- c(1, NA, NA, NA, NA)
#   
#   # If we're back to a stump no need to call fill_tree_details
#   if (nrow(new_tree$tree_matrix) == 1) {
#     new_tree$node_indices <- rep(1, length(y))
#   } else {
#     # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
#     if (node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
#       # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
#       bad_parents <- which(new_tree$tree_matrix[, "parent"] >= node_to_prune)
#       # Shift them back because you have removed two rows
#       new_tree$tree_matrix[bad_parents, "parent"] <- new_tree$tree_matrix[bad_parents, "parent"] - 2
#       
#       for (j in node_to_prune:nrow(new_tree$tree_matrix)) {
#         # Find the current parent
#         curr_parent <- new_tree$tree_matrix[j, "parent"]
#         # Find both the children of this node
#         curr_children <- which(new_tree$tree_matrix[, "parent"] == curr_parent)
#         # Input these children back into the parent
#         new_tree$tree_matrix[curr_parent, c("child_left", "child_right")] <- sort(curr_children)
#       } # End for loop of correcting parents and children
#     } # End if statement to fill in tree details
#     
#     # Call the fill function on this tree
#     new_tree <- fill_tree_details(new_tree, X)
#   }
#   
#   # Return new_tree
#   return(new_tree)
# } # End of prune_tree function


# change_tree function ----------------------------------------------------

# change_tree <- function(X, y, curr_tree, node_min_size) {
#   
#   # Change a node means change out the split value and split variable of an internal node. Need to make sure that this does now produce a bad tree (i.e. zero terminal nodes)
#   
#   # If current tree is a stump nothing to change
#   if (nrow(curr_tree$tree_matrix) == 1) {
#     return(curr_tree)
#   }
#   
#   # Create a holder for the new tree
#   new_tree <- curr_tree
#   
#   # Need to get the internal nodes
#   internal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 0)
#   terminal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 1)
#   
#   # Create a while loop to get good trees
#   # Create a counter to stop after a certain number of bad trees
#   max_bad_trees <- 2
#   count_bad_trees <- 0
#   bad_trees <- TRUE
#   while (bad_trees) {
#     # Re-set the tree
#     new_tree <- curr_tree
#     
#     # choose an internal node to change
#     node_to_change <- sample(internal_nodes, 1)
#     
#     # Use the get_children function to get all the children of this node
#     all_children <- get_children(new_tree$tree_matrix, node_to_change)
#     
#     # Now find all the nodes which match these children
#     use_node_indices <- !is.na(match(new_tree$node_indices, all_children))
#     
#     # Create new split variable and value based on ignorance
#     # then check this doesn't give a bad tree
#     new_split_variable <- sample(1:ncol(X), 1)
#     available_values <- sort(unique(X[
#       use_node_indices,
#       new_split_variable
#     ]))
#     new_split_value <- sample(available_values[-c(1, length(available_values))], 1)
#     
#     
#     # Update the tree details
#     new_tree$tree_matrix[
#       node_to_change,
#       c(
#         "split_variable",
#         "split_value"
#       )
#     ] <- c(
#       new_split_variable,
#       new_split_value
#     )
#     
#     # Update the tree node indices
#     new_tree <- fill_tree_details(new_tree, X)
#     
#     # Check for bad tree
#     if (any(new_tree$tree_matrix[terminal_nodes, "node_size"] == 0)) {
#       count_bad_trees <- count_bad_trees + 1
#     } else {
#       bad_trees <- FALSE
#     }
#     if (count_bad_trees == max_bad_trees) {
#       return(curr_tree)
#     }
#   } # end of while loop
#   
#   # Revert if the new tree has too small terminal node sizes
#   if (any(new_tree$tree_matrix[, "node_size"] < node_min_size)) new_tree <- curr_tree
#   
#   # Return new_tree
#   return(new_tree)
# } # End of change_tree function

# swap_tree function ------------------------------------------------------

# swap_tree <- function(X, y, curr_tree, node_min_size) {
#   
#   # Swap takes two neighbouring internal nodes and swaps around their split values and variables
#   
#   # If current tree is a stump nothing to change
#   if (nrow(curr_tree$tree_matrix) == 1) {
#     return(curr_tree)
#   }
#   
#   # Create a holder for the new tree
#   new_tree <- curr_tree
#   
#   # Need to get the internal nodes
#   internal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 0)
#   terminal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 1)
#   
#   # If less than 3 internal nodes return curr_tree
#   if (length(internal_nodes) < 3) {
#     return(curr_tree)
#   }
#   
#   # Find pairs of neighbouring internal nodes
#   parent_of_internal <- new_tree$tree_matrix[internal_nodes, "parent"]
#   pairs_of_internal <- cbind(internal_nodes, parent_of_internal)[-1, ]
#   
#   # Create a while loop to get good trees
#   # Create a counter to stop after a certain number of bad trees
#   max_bad_trees <- 2
#   count_bad_trees <- 0
#   bad_trees <- TRUE
#   while (bad_trees) {
#     # Re-set the tree
#     new_tree <- curr_tree
#     
#     # Pick a random pair
#     nodes_to_swap <- sample(1:nrow(pairs_of_internal), 1)
#     
#     # Get the split variables and values for this pair
#     swap_1_parts <- new_tree$tree_matrix[
#       pairs_of_internal[nodes_to_swap, 1],
#       c("split_variable", "split_value")
#     ]
#     swap_2_parts <- new_tree$tree_matrix[
#       pairs_of_internal[nodes_to_swap, 2],
#       c("split_variable", "split_value")
#     ]
#     
#     # Update the tree details - swap them over
#     new_tree$tree_matrix[
#       pairs_of_internal[nodes_to_swap, 1],
#       c(
#         "split_variable",
#         "split_value"
#       )
#     ] <- swap_2_parts
#     new_tree$tree_matrix[
#       pairs_of_internal[nodes_to_swap, 2],
#       c(
#         "split_variable",
#         "split_value"
#       )
#     ] <- swap_1_parts
#     
#     # Update the tree node indices
#     new_tree <- fill_tree_details(new_tree, X)
#     
#     # Check for bad tree
#     if (any(new_tree$tree_matrix[terminal_nodes, "node_size"] == 0)) {
#       count_bad_trees <- count_bad_trees + 1
#     } else {
#       bad_trees <- FALSE
#     }
#     if (count_bad_trees == max_bad_trees) {
#       return(curr_tree)
#     }
#   } # end of while loop
#   
#   # Reject if new tree has any too small terminal node sizes
#   if (any(new_tree$tree_matrix[, "node_size"] < node_min_size)) new_tree <- curr_tree
#   
#   # Return new_tree
#   return(new_tree)
# } # End of swap_tree function

# Fill_tree_details -------------------------------------------------------

# The fill tree details function takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in ecah terminal node
# fill_tree_details <- function(curr_tree, X) {
#   
#   # This code should essentially start from ignorance - no indices just a tree
#   # Fill in the number of observations and the node indices
#   
#   # Collect right bits of tree
#   tree_matrix <- curr_tree$tree_matrix
#   
#   # Create a new tree matrix to overwrite
#   new_tree_matrix <- tree_matrix
#   
#   # Start with dummy node indices
#   node_indices <- rep(1, nrow(X))
#   
#   # For all but the top row, find the number of observations falling into each one
#   for (i in 2:nrow(tree_matrix)) {
#     # Get the parent
#     curr_parent <- tree_matrix[i, "parent"]
#     
#     # Find the split variable and value of the parent
#     split_var <- tree_matrix[curr_parent, "split_variable"]
#     split_val <- tree_matrix[curr_parent, "split_value"]
#     
#     # Find whether it's a left or right terminal node
#     left_or_right <- ifelse(tree_matrix[curr_parent, "child_left"] == i,
#                             "left", "right"
#     )
#     if (left_or_right == "left") {
#       # If left use less than condition
#       new_tree_matrix[i, "node_size"] <- sum(X[node_indices == curr_parent, split_var] < split_val)
#       node_indices[node_indices == curr_parent][X[node_indices == curr_parent, split_var] < split_val] <- i
#     } else {
#       # If right use greater than condition
#       new_tree_matrix[i, "node_size"] <- sum(X[node_indices == curr_parent, split_var] >= split_val)
#       node_indices[node_indices == curr_parent][X[node_indices == curr_parent, split_var] >= split_val] <- i
#     }
#   } # End of loop through table
#   
#   return(list(
#     tree_matrix = new_tree_matrix,
#     node_indices = node_indices
#   ))
# } # End of function


# Get complete conditions -------------------------------------------------

# tree_full_conditional <- function(tree, R, tau, tau_mu) {
#   # Function to compute log full conditional distirbution for an individual tree
#   # R is a vector of partial residuals
#   
#   # Need to calculate log complete conditional, involves a sum over terminal nodes
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   # Get sum of residuals and sum of residuals squared within each terminal node
#   sumRsq_j <- stats::aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[, 2]
#   S_j <- stats::aggregate(R, by = list(tree$node_indices), sum)[, 2]
#   
#   # Now calculate the log posterior
#   # log_post = 0.5 * length(R) * log(tau) + sum(0.5 * log( tau_mu / (tau_mu + nj * tau)) -
#   #   0.5 * tau * (sumRsq_j - tau * S_j^2 / (tau_mu + nj * tau) ) )
#   log_post <- 0.5 * length(R) * log(tau) +
#     0.5 * (sum(log(tau_mu / (tau_mu + nj * tau))) -
#              tau * sum(sumRsq_j) +
#              tau^2 * sum(S_j^2 / (tau_mu + nj * tau)))
#   return(log_post)
#   #
#   # New Mahdi version - slower
#   # P1 = 0.5 * length(R) * log(tau)
#   # P2 = 0.5 * sum( log( tau_mu / (tau_mu + nj * tau)))
#   # P3 = -0.5 * tau * sum( sumRsq_j )
#   # P4 = 0.5 * (tau^2) * sum ( (S_j^2) / (tau_mu + nj * tau) )
#   #
#   # return(P1 + P2 + P3 + P4)
# }

# Tree conditional for HEBART
# logdet fuction useful for simple return of log determinant
# logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
# tree_full_conditional_hebart <- function(tree, R, k_1, k_2, M, nu, lambda) {
#   # Function to compute log full conditional distribution for an individual tree
#   # R is a vector of partial residuals
#   
#   # hbeart version is
#   # log_cond = sum( log(gamma(n_j/2 + alpha)) - 0.5 * logdet(W) - (n_j/2 + alpha) *
#   #    log(beta + 0.5 * t(R_j)%*%solve(W, R_j) )
#   # where now W = k_2 * ones %*% t(ones) + k_1 * M %*% t(M) + diag(n_j)
#   # where M is the group allocation matrix.
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   log_cond <- 0
#   for (i in 1:length(nj)) {
#     M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
#     R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
#     W <- k_2 * matrix(1, nrow = nj[i], ncol = nj[i]) + k_1 * tcrossprod(M_j) + diag(nj[i])
#     log_cond <- log_cond - 0.5 * logdet(W) + lgamma(nj[i] / 2 + nu / 2) - (nj[i] / 2 + nu / 2) *
#       log(lambda / 2 + 0.5 * t(R_j) %*% solve(W, R_j))
#     # There's also this term in the maths which I don't think is necessary
#     # - 0.5 * nj[i] * log(2 * pi)
#   }
#   
#   return(log_cond)
# }


# Get predictions ---------------------------------------------------------

# Gets the predicted values from a current set of trees
# get_predictions <- function(trees, X, single_tree = FALSE) {
#   
#   # Stop nesting problems in case of multiple trees
#   if (is.null(names(trees)) & (length(trees) == 1)) trees <- trees[[1]]
#   
#   # Normally trees will be a list of lists but just in case
#   if (single_tree) {
#     # Deal with just a single tree
#     if (nrow(trees$tree_matrix) == 1) {
#       predictions <- rep(trees$tree_matrix[1, "mu"], nrow(X))
#     } else {
#       # Loop through the node indices to get predictions
#       predictions <- rep(NA, nrow(X))
#       unique_node_indices <- unique(trees$node_indices)
#       # Get the node indices for the current X matrix
#       curr_X_node_indices <- fill_tree_details(trees, X)$node_indices
#       # Now loop through all node indices to fill in details
#       for (i in 1:length(unique_node_indices)) {
#         predictions[curr_X_node_indices == unique_node_indices[i]] <-
#           trees$tree_matrix[unique_node_indices[i], "mu"]
#       }
#     }
#     # More here to deal with more complicated trees - i.e. multiple trees
#   } else {
#     # Do a recursive call to the function
#     partial_trees <- trees
#     partial_trees[[1]] <- NULL # Blank out that element of the list
#     predictions <- get_predictions(trees[[1]], X, single_tree = TRUE) +
#       get_predictions(partial_trees, X,
#                       single_tree = length(partial_trees) == 1
#       )
#     # single_tree = !is.null(names(partial_trees)))
#     # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
#   }
#   
#   return(predictions)
# }


# Get group predictions ---------------------------------------------------

# Gets the predicted values from a current set of trees
# get_group_predictions <- function(trees, X, groups, single_tree = FALSE) {
#   
#   # Stop nesting problems in case of multiple trees
#   if (is.null(names(trees)) & (length(trees) == 1)) trees <- trees[[1]]
#   
#   num_groups <- length(unique(groups))
#   group_col_names <- unique(paste0("mu", groups))
#   
#   # Normally trees will be a list of lists but just in case
#   if (single_tree) {
#     # Deal with just a single tree
#     if (nrow(trees$tree_matrix) == 1) {
#       predictions <- trees$tree_matrix[1, group_col_names][groups]
#     } else {
#       # Loop through the node indices to get predictions
#       predictions <- rep(NA, nrow(X))
#       unique_node_indices <- unique(trees$node_indices)
#       # Get the node indices for the current X matrix
#       curr_X_node_indices <- fill_tree_details(trees, X)$node_indices
#       
#       # Now loop through all node indices to fill in details
#       for (i in 1:length(unique_node_indices)) {
#         curr_groups <- groups[curr_X_node_indices == unique_node_indices[i]]
#         predictions[curr_X_node_indices == unique_node_indices[i]] <-
#           trees$tree_matrix[unique_node_indices[i], paste0("mu", curr_groups)]
#       }
#     }
#     # More here to deal with more complicated trees - i.e. multiple trees
#   } else {
#     # Do a recursive call to the function
#     partial_trees <- trees
#     partial_trees[[1]] <- NULL # Blank out that element of the list
#     predictions <- get_group_predictions(trees[[1]], X, groups, single_tree = TRUE) +
#       get_group_predictions(partial_trees, X, groups,
#                             single_tree = length(partial_trees) == 1
#       )
#   }
#   
#   return(predictions)
# }


# Get tree priors ---------------------------------------------------------

# get_tree_prior <- function(tree, alpha, beta) {
#   # Returns the tree log prior score
#   
#   # Need to work out the depth of the tree
#   # First find the level of each node, then the depth is the maximum of the level
#   level <- rep(NA, nrow(tree$tree_matrix))
#   level[1] <- 0 # First row always level 0
#   
#   # Escpae quickly if tree is just a stump
#   if (nrow(tree$tree_matrix) == 1) {
#     return(log(1 - alpha)) # Tree depth is 0
#   }
#   
#   
#   for (i in 2:nrow(tree$tree_matrix)) {
#     # Find the current parent
#     curr_parent <- tree$tree_matrix[i, "parent"]
#     # This child must have a level one greater than it's current parent
#     level[i] <- level[curr_parent] + 1
#   }
#   
#   # Only compute for the internal nodes
#   internal_nodes <- which(tree$tree_matrix[, "terminal"] == 0)
#   log_prior <- 0
#   for (i in 1:length(internal_nodes)) {
#     log_prior <- log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]])
#   }
#   # Now add on terminal nodes
#   terminal_nodes <- which(tree$tree_matrix[, "terminal"] == 1)
#   for (i in 1:length(terminal_nodes)) {
#     log_prior <- log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
#   }
#   
#   
#   return(log_prior)
# }

# Simulate_mu -------------------------------------------------------------

# simulate_mu <- function(tree, R, tau, tau_mu) {
#   
#   # Simulate mu values for a given tree
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   # Get sum of residuals in each terminal node
#   sumR <- aggregate(R, by = list(tree$node_indices), sum)[, 2]
#   
#   # Now calculate mu values
#   mu <- rnorm(length(nj),
#               mean = tau * sumR / (nj * tau + tau_mu),
#               sd = sqrt(1 / (nj * tau + tau_mu))
#   )
#   
#   # Wipe all the old mus out for other nodes
#   tree$tree_matrix[, "mu"] <- NA
#   
#   # Put in just the ones that are useful
#   tree$tree_matrix[which_terminal, "mu"] <- mu
#   
#   return(tree)
# }

# simulate_mu_hebart <- function(tree, R, tau, k_1, k_2) {
#   
#   # Simulate mu values for a given tree
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   # Get sum of residuals in each terminal node
#   sumR <- aggregate(R, by = list(tree$node_indices), sum)[, 2]
#   
#   # Now calculate mu values
#   mu <- rnorm(length(nj),
#               mean = (sumR / k_1) / (nj / k_1 + 1 / k_2),
#               sd = sqrt(1 / (tau * nj / k_1 + 1 / k_2))
#   )
#   
#   # Wipe all the old mus out for other nodes
#   tree$tree_matrix[, "mu"] <- NA
#   
#   # Put in just the ones that are useful
#   tree$tree_matrix[which_terminal, "mu"] <- mu
#   
#   return(tree)
# }
# 
# simulate_mu_hebart2 <- function(tree, R, M, tau, k_1, k_2) {
#   
#   # Simulate mu values for a given tree
#   
#   # this is the marginalised mu version from
#   # https://bookdown.org/connect/#/apps/56e67516-559f-4d69-91df-54702fbc2206/access
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   # Wipe all the old mus out for other nodes
#   tree$tree_matrix[, "mu"] <- NA
#   
#   # Loop through terminal nodes to get values
#   for (i in 1:length(nj)) {
#     M_j <- M[tree$node_indices == which_terminal[i], , drop = FALSE]
#     R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
#     Psi <- k_1 * tcrossprod(M_j) + diag(nj[i])
#     ones <- rep(1, nj[i])
#     Prec_bit <- t(ones)%*%solve(Psi, ones) + 1/k_2
#     mean <- t(ones)%*%Psi%*%R_j / Prec_bit
#     tree$tree_matrix[which_terminal[i], "mu"] <- rnorm(1,
#                                                        mean,
#                                                        sd = 1/sqrt(tau*Prec_bit))
#   }
#   
#   return(tree)
# }


# Simulate mu groups hebart -----------------------------------------------

# simulate_mu_groups_hebart <- function(tree, R, groups, tau, k_1, k_2) {
#   
#   # Simulate the group mu values for a given tree
#   
#   # First find which rows are terminal nodes
#   which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
#   
#   # Get node sizes for each terminal node
#   nj <- tree$tree_matrix[which_terminal, "node_size"]
#   
#   num_groups <- length(unique(groups))
#   
#   # Get the group means in each terminal node
#   # Doing this with loops but probably can be faster
#   for (i in 1:length(nj)) {
#     curr_R <- R[tree$node_indices == which_terminal[i]]
#     curr_groups <- groups[tree$node_indices == which_terminal[i]]
#     curr_group_sizes <- table(curr_groups)
#     group_R_means <- aggregate(curr_R, by = list(curr_groups), "sum")[, 2]
#     curr_mu <- tree$tree_matrix[which_terminal[i], "mu"]
#     curr_group_mu <- rnorm(num_groups,
#                            mean = (curr_mu / k_1 + group_R_means) / (curr_group_sizes + 1 / k_1),
#                            sd = sqrt(1 / (curr_group_sizes + 1 / k_1))
#     )
#     tree$tree_matrix[which_terminal[i], paste0("mu", 1:num_groups)] <- curr_group_mu
#     
#   }
#   
#   # Wipe all the old mu groups out for other nodes
#   which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
#   tree$tree_matrix[which_non_terminal, paste0("mu", 1:num_groups)] <- NA
#   
#   return(tree)
# }


# Update tau --------------------------------------------------------------

# update_tau <- function(S, nu, lambda, n) {
#   # Update from maths in Github folder
#   tau <- rgamma(1,
#     shape = (nu + n) / 2,
#     rate = (S + nu * lambda) / 2
#   )
#   # Alternative
#   # tau = rgamma(1, shape = (nu + n) / 2 - 1, scale = 2 / (S + nu * lambda))
#
#   return(tau)
# }

# update_tau <- function(y, M, nu, lambda, num_groups, k_1, k_2) {
#   
#   n <- length(y)
#   W_1 <- (k_2 * matrix(1, nrow = n, ncol = n)) + (k_1 * M %*% t(M)) + diag(n)
#   S <- t(y) %*% solve(W_1, y)
#   
#   # Update from maths in Github folder
#   tau <- stats::rgamma(1,
#                        shape = (nu + n) / 2,
#                        rate = (S + nu * lambda) / 2
#   )
#   
#   return(tau)
# }

# get_children ------------------------------------------------------------

# A function which if, the current node is terminal, returns the node, or if not returns the children and calls the function again on the children
# get_children <- function(tree_mat, parent) {
#   # Create a holder for the children
#   all_children <- NULL
#   if (tree_mat[parent, "terminal"] == 1) {
#     # If the node is terminal return the list so far
#     return(c(all_children, parent))
#   } else {
#     # If not get the current children
#     curr_child_left <- tree_mat[parent, "child_left"]
#     curr_child_right <- tree_mat[parent, "child_right"]
#     # Return the children and also the children of the children recursively
#     return(c(
#       all_children,
#       get_children(tree_mat, curr_child_left),
#       get_children(tree_mat, curr_child_right)
#     ))
#   }
# }


# Predict function --------------------------------------------------------

# predict_hebart <- function(newX, new_groups, hebart_posterior,
#                            type = c("all", "median", "mean")) {
#   # Create predictions based on a new feature matrix
#   # Note that there is minimal error checking in this - newX needs to be right!
#   
#   # Create holder for predicted values
#   n_its <- length(hebart_posterior$sigma)
#   y_hat_mat <- matrix(NA,
#                       nrow = n_its,
#                       ncol = nrow(newX)
#   )
#   
#   # Now loop through iterations and get predictions
#   for (i in 1:n_its) {
#     # Get current set of trees
#     curr_trees <- hebart_posterior$trees[[i]]
#     # Use get_predictions function to get predictions
#     # y_hat_mat[i, ] <- get_predictions(curr_trees,
#     #                                         newX,
#     #                                         single_tree = length(curr_trees) == 1
#     # )
#     y_hat_mat[i, ] <- get_group_predictions(curr_trees,
#                                             newX,
#                                             new_groups,
#                                             single_tree = length(curr_trees) == 1
#     )
#   }
#   
#   # Sort out what to return
#   out <- switch(type,
#                 all = y_hat_mat,
#                 mean = apply(y_hat_mat, 2, "mean"),
#                 median = apply(y_hat_mat, 2, "median")
#   )
#   
#   return(out)
# } # end of predict function


# Simulate friedman -------------------------------------------------------

# sim_friedman <- function(n,
#                          pars = c(10, 20, 10, 5),
#                          p = 0, scale_err = 1) {
#   # Simulate some data using Friedman's example
#   # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
#   X <- matrix(runif(n * (5 + p)), nrow = n, ncol = 5 + p)
#   mean <- pars[1] * sin(pi * X[, 1] * X[, 2]) + pars[2] * (X[, 3] - 0.5)^2 +
#     pars[3] * X[, 4] + pars[4] * X[, 5]
#   y <- rnorm(n, mean, scale_err)
#   return(list(
#     y = y, X = X, df = cbind(y, X),
#     true_mean = mean, true_scale = scale_err
#   ))
# }
