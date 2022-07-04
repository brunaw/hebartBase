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
set.seed(2022)
df_real     <- lme4::sleepstudy %>% set_names(c('y', 'X1', 'group'))
df_real$y   <- c(scale(df_real$y))
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <-  train$group
# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1
pars           <- list(
  alpha = 0.95, beta = 2, k_1 = 0.1,
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
  num_trees      = 1,
  k_1_pars       = k1_pars)

hb_model

# Results are OK:
# Formula:
#   y ~ X1 
# 
# Number of trees:         1 
# Number of covariates:    1 
# Training error (MSE):    0.6886846 
# R Squared:               0.3113154 


# when num_trees = 5
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 15,
  k_1_pars       = k1_pars)
hb_model
# It breaks: 
# Formula:
#   y ~ X1 
# 
# Number of trees:         5 
# Number of covariates:    1 
# Training error (MSE):    9627.78 
# R Squared:               -9626.78
#----------------------------------------------------
# When we change k_1 --------------------------------
k1_pars        <-  list(sample_k1 = TRUE,
                        min_u     = 0,
                        max_u     = 2,
                        k1_prior  = TRUE)

# when num_trees = 1
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 1,
  k_1_pars       = k1_pars)

# Results are OK, and even improved:
# Formula:
#   y ~ X1 
# 
# Number of trees:         1 
# Number of covariates:    1 
# Training error (MSE):    0.2816423 
# R Squared:               0.7183577 


# when num_trees = 5
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 5,
  k_1_pars       = k1_pars)


# It breaks: 
# Formula:
#   y ~ X1 
# 
# Number of trees:         5 
# Number of covariates:    1 
# Training error (MSE):    30237.64 
# R Squared:               -30236.64 
#----------------------------------------------------