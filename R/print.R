##' @title Print hebart
##' @param x Object of class 'hebart'
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{hebart}}
##' @author Bruna Wundervald
##' @export
print.hebart <- function(x, ...) {
  cat("# ------------------------------------------- #\n")
  cat("# HEBART result\n")
  cat("# ------------------------------------------- #\n")
  cat("Formula:\n", deparse(x$formula), "\n\n")
  cat("Number of trees:        ", x$num_trees, "\n")
  cat("Number of covariates:   ", x$num_variables, "\n")
  cat("Training error (RMSE):   ", x$rmse, "\n")
  cat("R Squared:              ", x$r.squared, "\n")
}
