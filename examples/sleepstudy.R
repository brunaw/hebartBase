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
# library(hebartBase)

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- lme4::sleepstudy %>% set_names(c('y', 'X1', 'group'))
df_real$y   <- df_real$y
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <- train$group

# Model parameters -----------------------------------
group_variable <-  "group"
formula <- y ~ X1
num_trees <- 1
pars   <- list(
  alpha = 0.95, beta = 2,
  nu = 3, lambda = 0.1,
  tau_mu = 16 * num_trees,
  shape_tau_phi = 0.5,
  scale_tau_phi = 1 # These give mean ~2 and sd ~4
)

# Running the model ----------------------------------

hb_model <- hebart(formula,
                   data = train,
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
                                tau_phi = 10000),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1,
                               tau_phi_sd = 0)
                   )
hb_model

pp <- predict_hebart(newX = test, new_groups = test$group,
                     hebart_posterior  = hb_model, 
                     type = "mean")

sqrt(mean(pp - test$y)^2) 
cor(pp, test$y)
qplot(test$y, pp) + geom_abline()
qplot(1:length(hb_model$sigma), hb_model$sigma)
qplot(1:length(hb_model$sigma), hb_model$sigma_phi)

stop()

# Comparison to BART --------------------------
bart_0 = dbarts::bart2(formula, 
                       #n.trees = 15,
                       data = train,
                       test = test, 
                       keepTrees = TRUE)
pp <- bart_0$yhat.test.mean
sqrt(mean(pp - test$y)^2) # 5.512617 - 100 trees
cor(pp, test$y) #   0.4983258
qplot(test$y, pp) + geom_abline()

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + (1|group), train)
pp <- predict(lme_ss, test)
sqrt(mean(pp - test$y)^2) # 1.527163
cor(pp, test$y) # 0.8426536
qplot(test$y, pp)




