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
table(df_real$group)
df_real$y   <- c(scale(df_real$y))
df_real$X1 <- df_real$X1 + rnorm(nrow(df_real))
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <-  train$group
# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1
pars           <- list(
  alpha = 0.95, beta = 2, k_1 = 0.05,
  k_2 = 1, nu = 3, lambda = 0.1
)

# 1. Without varying k1: 
k1_pars        <-  list(sample_k1 = FALSE)

# Running the model ----------------------------------
# when num_trees = 1
train$group <- as.numeric(train$group)
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 5,
  k_1_pars       = k1_pars)

hb_model
hebartBase::grow_tree
# Results are OK:
# Formula:
#   y ~ X1 
# 
# Number of trees:         1 
# Number of covariates:    1 
# Training error (MSE):    0.6927752 
# R Squared:               0.3072248 


# when num_trees = 15
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
# y ~ X1 
# 
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.4352189 
# R Squared:               0.5647811 
#----------------------------------------------------
# When we change k_1 --------------------------------
k1_pars        <-  list(sample_k1 = TRUE,
                        min_u     = 0.1,
                        max_u     = 1,
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
# Results are OK, and even improved:
# Formula:
#   y ~ X1 
# 
# Number of trees:         1 
# Number of covariates:    1 
# Training error (MSE):    0.481118 
# R Squared:               0.518882 


# when num_trees = 15
hb_model <- hebart(
  formula,
  group_variable = "group",
  data           = train,
  control        = pars,
  num_trees      = 3,
  k_1_pars       = k1_pars)

hb_model
# ------------------------------------------- #
# HEBART result
# ------------------------------------------- #
# Formula:
#   y ~ X1 
# 
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.4217759 
# R Squared:               0.5782241
#----------------------------------------------------