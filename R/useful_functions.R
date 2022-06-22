#' @name logdet
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @title logdet fuction useful for simple return of log determinant
#' @description Gets tree children
#' @param A The matrix
#' @export
logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)


#' @name get_children
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Get children
#' @description Gets tree children
#' @param tree_mat The tree matrix
#' @param parent The corresponding parent
#' 
# A function which if, the current node is terminal, returns the node, or if not returns the children and calls the function again on the children
get_children <- function(tree_mat, parent) {
  # Create a holder for the children
  all_children <- NULL
  if (tree_mat[parent, "terminal"] == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left <- tree_mat[parent, "child_left"]
    curr_child_right <- tree_mat[parent, "child_right"]
    # Return the children and also the children of the children recursively
    return(c(
      all_children,
      get_children(tree_mat, curr_child_left),
      get_children(tree_mat, curr_child_right)
    ))
  }
}

#' @name sim_friedman
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Simulate Friedman
#' @description Simulates data accordingly to the Friedman equation
#' @param n The number of simulated points
#' @param pars The model parameters
#' @param p The number of noise variables
#' @param scale_err The scale parameter of the simulation

# Simulate friedman -------------------------------------------------------
sim_friedman <- function(n,
                         pars = c(10, 20, 10, 5),
                         p = 0, scale_err = 1) {
  # Simulate some data using Friedman's example
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  X <- matrix(stats::runif(n * (5 + p)), nrow = n, ncol = 5 + p)
  mean <- pars[1] * sin(pi * X[, 1] * X[, 2]) + pars[2] * (X[, 3] - 0.5)^2 +
    pars[3] * X[, 4] + pars[4] * X[, 5]
  y <- stats::rnorm(n, mean, scale_err)
  return(list(
    y = y, X = X, df = cbind(y, X),
    true_mean = mean, true_scale = scale_err
  ))
}
