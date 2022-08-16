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
y <- data$y
group <- sample(1:3, size = length(y), replace = TRUE)
train <- data.frame(y, X, group)

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1 + X2 + X3 + X4 + X5

num_trees <- 4
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
                                tau_phi = 1),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1,
                               tau_phi_sd = 2)
                   )
pp <- predict_hebart(test$X, group_test, hb_model, type = "mean")
sqrt(mean(pp - test$y)^2)
cor(pp, scale(test$y))
qplot(1:length(hb_model$sigma), hb_model$sigma)
qplot(1:length(hb_model$sigma), hb_model$sigma_phi)
qplot(test$y, pp) + geom_abline()
