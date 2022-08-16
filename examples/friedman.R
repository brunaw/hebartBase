# In case you're installing, building, or removing the package:
# remove.packages("hebartBase")
# devtools::document()
# devtools::check()
# devtools::install()

# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(tidymodels)
library(hebartBase)

# Dataset split  ------------------------------------
set.seed(1)
data <- sim_friedman(n = 500)
test <- sim_friedman(n = 200)
group_test <- sample(1:3, size = length(test$y), replace = TRUE)
X <- data$X
y <- scale(data$y)[, 1]
group <- sample(1:3, size = length(y), replace = TRUE)
train <- data.frame(y, X, group)

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1 + X2 + X3 + X4 + X5

num_trees <- 3
pars   <- list(
  alpha = 0.95, beta = 2,
  nu = 3, lambda = 0.1,
  tau_mu = 16 * num_trees,
  shape_tau_phi = 0.5,
  scale_tau_phi = 1 # These give mean ~2 and sd ~4
)

# Running the model ----------------------------------

hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = num_trees,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     nu = 2,
                     lambda = 0.1,
                     tau_mu = 16 * num_trees,
                     shape_tau_phi = 0.5,
                     scale_tau_phi = 1
                   ), 
                   inits = list(tau = 1,
                                tau_phi = 1000),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1)
                   )
pp <- predict_hebart(test$X, group_test, hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.021
cor(pp, scale(test$y)) # 0.021
hb_model$sigma 
qplot(scale(test$y), pp) + geom_abline()
stop()

#----------------------------------------------------
# When we change k_1 --------------------------------
k_1_pars        <-  list(sample_k1 = TRUE,
                        min_u     = 0.1,
                        max_u     = 15,
                        k1_prior  = TRUE)

# when num_trees = 5
hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 10,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 1e-10,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 250, # Number of iterations
                     burn = 100, # Size of burn in
                     thin = 1
                   ))
hb_model
mean(hb_model$samples_k1)

pp <- predict_hebart(test$X, group_test, hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.05077673
cor(pp, scale(test$y)) # 0.021
hb_model$sigma 
qplot(test$y, pp)
# The means get really weird
#----------------------------------------------------
# A grouped version
data1 <- sim_friedman(n = 250,
                      scale_err = 1.5)
data2 <- sim_friedman(n = 250,
                      pars = c(10, 20, 10, 5)*5, 
                      scale_err = 1.5)
data3 <- sim_friedman(n = 250, 
                      pars = c(10, 20, 10, 5)*10,
                      scale_err = 1.5)
df <- bind_rows(as.data.frame(data1$X), 
                as.data.frame(data2$X), 
                as.data.frame(data3$X)) |> 
  mutate(y = c(data1$y, data2$y, data3$y),
         group = rep(c("1", "2", "3"), each = 250))
#test <- sim_friedman(n = 200)
#train <- data.frame(y, X, group)

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ V1 + V2 + V3 + V4 + V5
pars           <- list(
  alpha = 0.95, beta = 2, k_1 = 1,
  k_2 = 1, nu = 3, lambda = 0.1
)

# 1. Without varying k1: 
k1_pars        <-  list(sample_k1 = FALSE)

hb_model <- hebart(formula,
                   data           = df,
                   group_variable = "group", 
                   num_trees = 2,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.0001,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ),
                   MCMC = list(
                     iter = 1000, # Number of iterations
                     burn = 100, # Size of burn in
                     thin = 1
                   ))
pp <- predict_hebart(df, df$group, hb_model, type = "mean")
sqrt(mean(pp - scale(df$y))^2) # 0.021
cor(pp, scale(df$y)) # 0.7584812
# hb_model$sigma 
qplot(df$y, pp)
#----------------------------------------------------
# When we change k_1 --------------------------------
k_1_pars        <-  list(sample_k1 = TRUE,
                         min_u     = 0.1,
                         max_u     = 15,
                         k1_prior  = FALSE)

# when num_trees = 5
hb_model <- hebart(formula,
                   data           = df,
                   group_variable = "group", 
                   num_trees = 2,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 1e-10,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 1000, # Number of iterations
                     burn = 100, # Size of burn in
                     thin = 1
                   ))

hb_model
pp <- predict_hebart(df, df$group, hb_model, type = "mean")
sqrt(mean(pp - scale(df$y))^2) # 0.021
cor(pp, scale(df$y)) # 0.79
mean(hb_model$samples_k1)
#----------------------------------------------------