# In case you're installing, building, or removing the package:
# remove.packages("hebartBase")
# devtools::document()
# devtools::check()
# devtools::install()

#----------------------------------------------------------------------
# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(lme4)
library(tidymodels)
library(dbarts)
library(hebartBase)
library(brms)

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
median(hb_model$sigma_phi)  # 0.1252214
sqrt(mean((pp - test$y)^2)) # 31.685
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
sqrt(mean((pp - test$y)^2)) # 34.8079 
cor(pp, test$y)
qplot(test$y, pp) + geom_abline()

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + (1|group), train)
pp_lme <- predict(lme_ss, test)
sqrt(mean((pp_lme - test$y)^2))  # 33.19528
cor(pp_lme, scale(test$y)) 
qplot(test$y, pp_lme)

# Comparison to BLME --------------------------
pr = prior(normal(0, 1), class = 'b')
blme <-  brm(
  y ~ X1 + (1 |group), data = train, prior = pr, cores = 4)
blme
pp_blme <- predict(blme, test)
sqrt(mean((pp_blme - test$y)^2))  # 33.19528
cor(pp_blme, scale(test$y)) 
qplot(test$y, pp_blme[, "Estimate"])
#----------------------------------------------------------------------
#----------------------------------------------------------------------
