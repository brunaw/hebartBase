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
pars           <- list(
  alpha = 0.95, beta = 2, k_1 = 1,
  k_2 = 1, nu = 3, lambda = 0.1
)

# 1. Without varying k1: 
k1_pars        <-  list(sample_k1 = FALSE)

# Running the model ----------------------------------
# when num_trees = 1
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 5,
  k_1_pars       = k1_pars)

hb_model
pp <- predict_hebart(test$X, group_test, hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.021
cor(pp, scale(test$y)) # 0.021
hb_model$sigma 
qplot(test$y, pp)
# Formula:
#   y ~ X1 + X2 + X3 + X4 + X5 
# 
# Number of trees:         1 
# Number of covariates:    5 
# Training error (MSE):    0.3658532 
# R Squared:               0.6341468 


# when num_trees = 15
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 5,
  k_1_pars       = k1_pars)
hb_model
pp <- predict_hebart(test$X, group_test, hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.02792329
# Formula:
#   y ~ X1 + X2 + X3 + X4 + X5 
# 
# Number of trees:         5 
# Number of covariates:    5 
# Training error (MSE):    0.18127 
# R Squared:               0.81873 
#----------------------------------------------------
# When we change k_1 --------------------------------
k1_pars        <-  list(sample_k1 = TRUE,
                        min_u     = 0.1,
                        max_u     = 5,
                        k1_prior  = TRUE)

# when num_trees = 1
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 1,
  k_1_pars       = k1_pars)
hb_model
mean(hb_model$samples_k1)
# Results are OK, and even improved:
# Formula:
#   y ~ X1 + X2 + X3 + X4 + X5 
# 
# Number of trees:         1 
# Number of covariates:    5 
# Training error (MSE):    0.3347007 
# R Squared:               0.6652993 


# when num_trees = 3
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 3,
  k_1_pars       = k1_pars)

hb_model
# Formula:
# y ~ X1 + X2 + X3 + X4 + X5 
# 
# Number of trees:         3 
# Number of covariates:    5 
# Training error (MSE):    0.2022467 
# R Squared:               0.7977533 

# when num_trees = 5
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 5,
  k_1_pars       = k1_pars)

hb_model

# Formula:
#   y ~ X1 + X2 + X3 + X4 + X5 
# 
# Number of trees:         5 
# Number of covariates:    5 
# Training error (MSE):    0.1694411 
# R Squared:               0.8305589

# when num_trees = 12
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 10,
  k_1_pars       = k1_pars)

hb_model$samples_k1
#----------------------------------------------------