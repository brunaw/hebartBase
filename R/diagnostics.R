#' @name get_avg_nodes
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Gets node averages
#' @description Get node averages per tree
#' @param model The hebart model object
#' @return The correspondent averages
#' 
get_avg_nodes <- function(model){
  all_trees <- model$trees
  n_iters   <- length(all_trees)
  avg_node <- vector(length = n_iters)
  
  for(i in 1:n_iters){
    trees_i   <- all_trees[[i]]
    trees_i   <- purrr::map(trees_i, "tree_matrix")
    avg_node[i]  <- base::mean(purrr::map_dbl(
      trees_i, 
      ~{which_terminal <- which(.x[, "terminal"] == 1)
      round(length(which_terminal), 3)
      }))
  }
  return(avg_node)
}


#' @name plot_avg_nodes
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Plots node averages
#' @description Gets node averages per tree and plots it
#' @param model The hebart model object
#' @return The correspondent averages

plot_avg_nodes <- function(model){
  avg_nodes <- dplyr::tibble(
    `Node Average` = get_avg_nodes(model)
  ) |> 
    dplyr::mutate(`Iteration` = 1:n())
  
  ggplot2::ggplot(data = avg_nodes, 
                  ggplot2::aes(y = `Node Average`, 
                               x = `Iteration`)) +
    ggplot2::geom_hline(yintercept = mean(avg_nodes$`Node Average`),
                        colour = '#c95a49', size = 0.5, linetype = 'dotted') +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::theme_bw(12) +
    ggplot2::labs(y = "Average Number of Nodes")
}


#' @name density_plot
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Density plots for tau and sigma_phi
#' @description Plots the sampled of values as a density
#' @param model The hebart model object
#' @param type Type of plot: tau or sigma_phi 
#' @param sqrt Logical to decide whether to plot 1/sqrt(tau) instead
#' (residual precision)
#' @return The correspondent plot
density_plot <- function(model, type = 'tau', sqrt = FALSE){
  if(type == 'tau'){
    df_tau <- data.frame(tau = model$tau)
    if(sqrt){
      df_tau$tau <- 1/sqrt(df_tau$tau)
      label_x <- expression('Samples of 1/'~sigma[tau])
    } else {
      label_x <- expression('Samples of '~tau)
    }
    
    ggplot2::ggplot(df_tau, ggplot2::aes(x = tau)) +
      ggplot2::geom_vline(xintercept = mean(df_tau$tau),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::labs(x = label_x, y = "Density") +
      ggplot2::theme_bw(12)
  } else if( type == 'sigma_phi'){
    
    df_sigma_phi <- data.frame(sigma_phi = model$sigma_phi)
    label_x <- expression('Samples of '~sigma[phi])
    
    ggplot2::ggplot(df_sigma_phi, ggplot2::aes(x = sigma_phi)) +
      ggplot2::geom_vline(xintercept = mean(df_sigma_phi$sigma_phi),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
      ggplot2::labs(x = label_x, y = "Density") +
      ggplot2::theme_bw(12)
  } else{
    stop("Type of plot not available")
  }
}

#' @name traceplot
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Traceplot for tau or sigma_phi 
#' @description Plots the sampled of values of tau, k1 or sqrt(k1/tau) per iteration
#' @param model The hebart model object
#' @param type Type of plot: tau or sigma_phi 
#' @param sqrt Logical to decide whether to plot 1/sqrt(tau) instead
#' (residual precision)
#' @return The correspondent plot

traceplot <- function(model, type = 'tau', sqrt = FALSE){
  if(type == 'tau'){
    df_tau <- data.frame(tau = model$tau, iter = 1:length(model$tau))
    if(sqrt){
      df_tau$tau <- 1/sqrt(df_tau$tau)
      label_y <- expression('Samples of 1/'~sigma[tau])
    } else {
      label_y <- expression('Samples of '~tau)
    }
    
    ggplot2::ggplot(df_tau, ggplot2::aes(y = tau, x = iter)) +
      ggplot2::geom_hline(yintercept = mean(df_tau$tau),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_line(alpha = 0.4) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(y = label_y, x = "Iteration") +
      ggplot2::theme_bw(12)
    
  } else if( type == 'sigma_phi'){
    
    df_sigma_phi <- data.frame(sigma_phi = model$sigma_phi,
                               iter = 1:length(model$sigma_phi))
    label_y <- expression('Samples of '~sigma[phi])
    
    ggplot2::ggplot(df_sigma_phi, ggplot2::aes(y = sigma_phi, x = iter)) +
      ggplot2::geom_hline(yintercept = mean(df_sigma_phi$sigma_phi),
                          colour = '#c95a49', size = 0.5, linetype = 'dotted') +
      ggplot2::geom_line(alpha = 0.4) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(y = label_y, x = "Iteration") +
      ggplot2::theme_bw(12)
    
  } else{
    stop("Type of plot not available")
  }
}


#' @name plot_mse_iter
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title MSE per iteration
#' @description Plots training MSE values per iteration
#' @param model The HEBART model object
#' @return The correspondent plot

plot_mse_iter <- function(model){
  
  df_avg <- data.frame(mse = model$mse, iter = 1:length(model$mse)) 
  # Plotting -----
  label_y <-  expression('MSE per iteration')
  ggplot2::ggplot(df_avg, ggplot2::aes(y = mse, x = iter)) +
    ggplot2::geom_hline(yintercept = mean(df_avg$mse),
                        colour = '#c95a49', size = 0.5, linetype = 'dotted') +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::labs(y = label_y, x = "Iteration") +
    ggplot2::theme_bw(12)
  
}



#' @name diagnostics
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Plots all diagnostics
#' @description Plots all diagnostics plots
#' @param model The HEBART model object
#' @return The correspondent plot

diagnostics <- function(model){
  library(patchwork)
  # Traceplots --------------------
  p1t <- traceplot(model, type = "tau", sqrt = TRUE)
  p2t <- traceplot(model, type = "sigma_phi")
  p1t <- p1t +
    ggplot2::ggtitle("Traceplots") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "italic"))
  
  # Density plots --------------------
  p1d <- density_plot(model, type = "tau", sqrt = TRUE)
  p2d <- density_plot(model, type = "sigma_phi")
  p1d <- p1d +
    ggplot2::ggtitle("Density plots") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "italic"))
  
  # Number of nodes ------------------
  p_node <- plot_avg_nodes(model) +
    ggplot2::ggtitle("Nodes per tree (average)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "italic"))
  
  # MSE per iteration ------------------
  p_mse <- plot_mse_iter(model) +
    ggplot2::ggtitle("MSE per iteration") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "italic"))
  
  (p1t + p2t) / (p1d + p2d)  / (p_node + p_mse)
  
}




