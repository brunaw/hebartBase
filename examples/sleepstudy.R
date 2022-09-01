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
library(dbarts)
#library(hebartBase)
devtools::load_all(".")

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- lme4::sleepstudy %>% set_names(c('y', 'X1', 'group'))
df_real$y   <- df_real$y
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <- train$group
num_trees   <- 15

# Running the model ----------------------------------

hb_model <- hebart(formula = y ~ X1,
                   data = train,
                   group_variable = "group", 
                   num_trees = num_trees,
                   priors = list(
                     alpha = 0.95, # Prior control list
                     beta = 2,
                     nu = 2,
                     lambda = 0.1,
                     tau_mu = 16 * num_trees,
                     shape_sigma_phi = 0.5,
                     scale_sigma_phi = 1,
                     sample_sigma_phi = TRUE
                   ), 
                   inits = list(tau = 1,
                                sigma_phi = 1),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1,
                               sigma_phi_sd = 0.5)
                   )
# Let's not use matrices 
pp <- predict_hebart(newX = test, new_groups = test$group,
                     hebart_posterior  = hb_model, 
                     type = "mean")
median(hb_model$sigma_phi)  # 0.1287825 
sqrt(mean((pp - test$y)^2)) # 31.35084
cor(pp, test$y)

diagnostics(hb_model)
qplot(test$y, pp)

# Comparison to BART --------------------------
bart_0 = dbarts::bart2(y ~ X1 + group, 
                       #n.trees = 15,
                       data = train,
                       test = test,
                       keepTrees = TRUE)
pp <- bart_0$yhat.test.mean
sqrt(mean((pp - test$y)^2)) # 53.28167 - 100 trees
cor(pp, test$y) #   0.4983258
qplot(test$y, pp) + geom_abline()

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + (1|group), train)
<<<<<<< HEAD
pp_lme <- predict(lme_ss, test)
sqrt(mean((pp_lme - test$y)^2))
cor(pp_lme, scale(test$y)) 
qplot(test$y, pp_lme)
#----------------------------------------------------------------------
||||||| c09f7be
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
#plot(hb_model$samples_k1)
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
||||||| 494c788
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
=======
sqrt(mean((pp - test$y)^2)) # 33.19528
cor(pp, test$y) # 0.8426536
>>>>>>> AP_Aug22
qplot(test$y, pp) + geom_abline()
=======
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
#plot(hb_model$samples_k1)
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
sqrt(mean((pp - test$y)^2)) # 33.19528
cor(pp, test$y) # 0.8426536

qplot(test$y, pp) + geom_abline()
>>>>>>> 2cf33eb5e63f1022110309ced37a4419ecb626d5




