#' @name get_predictions
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title HEBART Group predictions
#' @description A function that returns the group predictions for a hebart model
#' @param trees The current trees
#' @param X The set of covariates
#' @param single_tree Logical to indicate whether we only have one tree
#'
get_predictions <- function(trees, X, single_tree = FALSE) {
  
  # Stop nesting problems in case of multiple trees
  if (is.null(names(trees)) & (length(trees) == 1)) trees <- trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if (single_tree) {
    # Deal with just a single tree
    if (nrow(trees$tree_matrix) == 1) {
      predictions <- rep(trees$tree_matrix[1, "mu"], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions <- rep(NA, nrow(X))
      unique_node_indices <- unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices <- fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for (i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] <-
          trees$tree_matrix[unique_node_indices[i], "mu"]
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees <- trees
    partial_trees[[1]] <- NULL # Blank out that element of the list
    predictions <- get_predictions(trees[[1]], X, single_tree = TRUE) +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1
      )
    # single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }
  
  return(predictions)
}


#' @name get_group_predictions
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title HEBART Group predictions
#' @description A function that returns the group predictions for a hebart model
#' @param trees The current trees
#' @param X The set of covariates
#' @param groups The groups specification
#' @param single_tree Logical to indicate whether we only have one tree
#' @param old_groups The groups used during training
# Get group predictions ---------------------------------------------------

get_group_predictions <- function(trees, X, groups, single_tree = FALSE, 
                                  old_groups) {
  
  
  train_groups <- unique(old_groups)
  new_groups <- unique(groups)
  # Are those new groups?
  which_new <- new_groups[!(new_groups %in% train_groups)]
  if(length(which_new) > 1){
    X_new <- X[groups %in% which_new, ]
    inds_new <- which(groups %in% which_new)
    inds_old <- which(!groups %in% which_new)
    pred_new <- get_predictions(trees, X_new, single_tree = single_tree)
    X_old <- X[!(groups %in% which_new), ]
    if(nrow(X_new) == nrow(X)){
      return(pred_new)
    }
  } else{
    X_old <- X
  }
  
  # Stop nesting problems in case of multiple trees
  if (is.null(names(trees)) & (length(trees) == 1)) trees <- trees[[1]]
  
  group_names     <- unique(groups)
  num_groups      <- length(unique(groups))
  group_col_names <- paste0("phi", group_names)
  group_col_names_all <- paste0("phi", groups)
  
  # Normally trees will be a list of lists but just in case
  if (single_tree) {
    # Deal with just a single tree
    if (nrow(trees$tree_matrix) == 1) {
      predictions <- trees$tree_matrix[1, group_col_names][group_col_names_all]
    } else {
      # Loop through the node indices to get predictions
      predictions <- rep(NA, nrow(X_old))
      unique_node_indices <- unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices <- fill_tree_details(trees, X_old)$node_indices
      actual_node_indices  <- unique(curr_X_node_indices)
      # Now loop through all node indices to fill in details
      for (i in 1:length(actual_node_indices)) {
        curr_groups <- groups[curr_X_node_indices == actual_node_indices[i]]
        predictions[curr_X_node_indices == actual_node_indices[i]] <-
          trees$tree_matrix[actual_node_indices[i], paste0("phi", curr_groups)]
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees <- trees
    partial_trees[[1]] <- NULL # Blank out that element of the list
    predictions <- get_group_predictions(trees[[1]], X_old, groups, single_tree = TRUE, old_groups = old_groups) +
      get_group_predictions(partial_trees, X_old, groups,
                            single_tree = length(partial_trees) == 1,
                            old_groups = old_groups
      )
  }
  
  if(exists("pred_new")){
    final <- data.frame(ind = c(inds_new, inds_old), 
                        predictions = pred_new, predictions)
    predictions <- final$predictions
  }
  
  return(predictions)
}


#' @name predict_hebart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title HEBART Predictions
#' @description A function that returns the predictions for a hebart model
#' @param newX The new set of covariates
#' @param new_groups The new groups specification
#' @param hebart_posterior The posterior values from the model
#' @param type The prediction type ("all", "median" or "mean")
# Predict function --------------------------------------------------------
predict_hebart <- function(newX, new_groups, hebart_posterior,
                           type = c("all", "median", "mean")) {
  # Create predictions based on a new feature matrix
  # Note that there is minimal error checking in this - newX needs to be right!
  if(!is.data.frame(newX)){
    stop("Please use a data.frame")
  } 
  formula     <- hebart_posterior$formula
  response_name <- all.vars(formula)[1]
  names_x <- all.vars(formula[[3]])
  new_formula <- paste0("~", paste0(names_x, collapse = "+"))
  new_formula <- stats::as.formula(new_formula)
  formula_int <- stats::as.formula(paste(c(new_formula), "- 1"))
  mf   <- stats::model.frame(formula_int,  data = newX)
  newX <- as.matrix(stats::model.matrix(formula_int, mf))
  old_groups <- hebart_posterior$groups
#   if(is.data.frame(newX)){
#     newX <- matrix(newX[, names_x], ncol = length(names_x))
#   }
#   
#   if(ncol(newX) != length(names_x)){
#     stop("Please use a matrix with the same \
# number of covariates used during training")
#   }
  

  # Create holder for predicted values
  n_its <- length(hebart_posterior$sigma)
  y_hat_mat <- matrix(NA,
                      nrow = n_its,
                      ncol = nrow(newX)
  )
  
  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees <- hebart_posterior$trees[[i]]
    # Use get_predictions function to get predictions
    # y_hat_mat[i, ] <- get_predictions(curr_trees,
    #                                         newX,
    #                                         single_tree = length(curr_trees) == 1
    # )
    y_hat_mat[i, ] <- get_group_predictions(trees = curr_trees,
                                            X = newX,
                                            groups = new_groups,
                                            single_tree = length(curr_trees) == 1,
                                            old_groups = old_groups
    )
  }
  
  # Sort out what to return
  inv_scale <- function(x) (x + 0.5) * (hebart_posterior$y_max - hebart_posterior$y_min) + hebart_posterior$y_min
  out <- switch(type,
                all = inv_scale(y_hat_mat),
                mean = apply(inv_scale(y_hat_mat), 2, "mean"),
                median = apply(inv_scale(y_hat_mat), 2, "median")
  )
  
  return(out)
} # end of predict function
