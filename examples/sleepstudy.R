# In case you're installing, building, or removing the package:
# remove.packages("hebartBase")
# devtools::document()
# devtools::check()
# devtools::install()

# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(lme4)
library(tidymodels)
library(hebartBase)

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- lme4::sleepstudy %>% set_names(c('y', 'X1', 'group'))
df_real$y   <- c(scale(df_real$y))
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <- train$group
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
# when num_trees = 15

hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 3,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 1e-10,
                     k_2 = 3,
                     nu = 2,
                     lambda = 0.1
                   ), 
                   MCMC = list(iter = 1500, burn = 250, thin = 10)
                   )
hb_model
# ------------------------------------------- #
# HEBART result
# ------------------------------------------- #
# Formula:
#   NULL 
# 
# Number of trees:         3 
# Number of covariates:    1 
# Training error (MSE):    0.7093291 
# R Squared:               0.2906709 

pp <- predict_hebart(newX = test, new_groups = test$group,
                     hebart_posterior  = hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.6845038
cor(pp, scale(test$y)) # 0.7556438

# qplot(test$y, pp)
# Comparison to BART --------------------------
bart_0 = dbarts::bart2(formula, 
                       #n.trees = 15,
                       data = train,
                       test = test, 
                       keepTrees = TRUE)
pp <- bart_0$yhat.test.mean
sqrt(mean(pp - scale(test$y))^2) # 0.08045881
cor(pp, scale(test$y)) #  0.4984845
qplot(test$y, pp)

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + (1|group), train)
pp <- predict(lme_ss, test)
sqrt(mean(pp - scale(test$y))^2) # 0.009018843
cor(pp, scale(test$y)) # 0.8426536
qplot(test$y, pp)


# When we change k_1 --------------------------------
k_1_pars        <-  list(sample_k1 = TRUE,
                         min_u     = 0.1,
                         max_u     = 1,
                         k1_prior  = FALSE)

# when num_trees = 5
hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 5,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.5,
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
# Number of trees:         5 
# Number of covariates:    1 
# Training error (MSE):    0.2563552 
# R Squared:               0.7436448 
plot(hb_model$samples_k1)
pp <- predict_hebart(newX = matrix(test$X1, ncol = 1), new_groups = test$group,
                     hebart_posterior  = hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.01975726 (much better)
cor(pp, scale(test$y)) # 0.8391904
qplot(test$y, pp) + geom_abline()


# when num_trees = 15 and iter  = 250
k_1_pars        <-  list(sample_k1 = TRUE,
                         min_u     = 0.1,
                         max_u     = 1,
                         k1_prior  = FALSE)

hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 15,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.5,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 250, # Number of iterations
                     burn = 50, # Size of burn in
                     thin = 1
                   ))
hb_model
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.279567 
# R Squared:               0.720433
plot(hb_model$samples_k1)
pp <- predict_hebart(newX = matrix(test$X1, ncol = 1), new_groups = test$group,
                     hebart_posterior  = hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.01617428 (much better)
cor(pp, scale(test$y)) # 0.7890124
qplot(test$y, pp) + geom_abline()

# when num_trees = 15 and iter  = 750
hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 15,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.5,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 750, # Number of iterations
                     burn = 250, # Size of burn in
                     thin = 1
                   ))
hb_model
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.266179 
# R Squared:               0.733821
plot(hb_model$samples_k1)
pp <- predict_hebart(newX = matrix(test$X1, ncol = 1), new_groups = test$group,
                     hebart_posterior  = hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.01617428 (much better)
cor(pp, scale(test$y)) # 0.7890124
qplot(test$y, pp) + geom_abline()


# when num_trees = 15 and iter  = 250, and k_1 can vary more
k_1_pars        <-  list(sample_k1 = TRUE,
                         min_u     = 0.1,
                         max_u     = 5,
                         k1_prior  = FALSE)

hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 15,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.5,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 250, # Number of iterations
                     burn = 50, # Size of burn in
                     thin = 1
                   ))
hb_model
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.2967165 
# R Squared:               0.7032835
# plot(hb_model$samples_k1)
# pp <- predict_hebart(newX = matrix(test$X1, ncol = 1), new_groups = test$group,
#                      hebart_posterior  = hb_model, type = "mean")
# sqrt(mean(pp - scale(test$y))^2) # 0.01617428 (much better)
# cor(pp, scale(test$y)) # 0.7890124
# qplot(test$y, pp) + geom_abline()

hb_model <- hebart(formula,
                   data           = train,
                   group_variable = "group", 
                   num_trees = 15,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     k_1 = 0.5,
                     k_2 = 0.2,
                     nu = 3,
                     lambda = 0.1
                   ), 
                   k_1_pars       = k_1_pars,
                   MCMC = list(
                     iter = 1500, # Number of iterations
                     burn = 250, # Size of burn in
                     thin = 1
                   ))
hb_model
# Number of trees:         15 
# Number of covariates:    1 
# Training error (MSE):    0.2818765 
# R Squared:               0.7181235 
plot(hb_model$samples_k1)
pp <- predict_hebart(newX = matrix(test$X1, ncol = 1), new_groups = test$group,
                     hebart_posterior  = hb_model, type = "mean")
sqrt(mean(pp - scale(test$y))^2) # 0.01633317
cor(pp, scale(test$y)) # 0.6996391
qplot(test$y, pp) + geom_abline()

#-------------------------------------------------
# Check one particular observation
which.max(pp) # 27
test[27,]
# y X1 group
# 100 2.847725  9   337

pp_27 <- predict_hebart(newX = matrix(test$X1[27], ncol = 1), new_groups = test$group[27],
                        hebart_posterior  = hb_model, type = "all")
pp_27[150] # 1.302005

# Now compare with the trees
hb_model$trees[[150]][[1]]$tree_matrix
hb_model$trees[[150]][[2]]$tree_matrix
hb_model$trees[[150]][[3]]$tree_matrix

# Observation 99 should be more or less identical in prediction
# and is in the training set
which(rownames(train) == "99") # 43
train[43,]
# 2.793536  8   337 
pp_43 <- predict_hebart(newX = matrix(train$X1[43], ncol = 1), new_groups = train$group[43],
                        hebart_posterior  = hb_model, type = "all")
pp_43[150]




