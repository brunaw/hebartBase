library(hebartBase)
library(magrittr)
library(ggplot2)
library(tidymodels)
library(lme4)

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- lme4::sleepstudy %>% set_names(c('y', 'X1', 'group'))
df_real$y   <- c(scale(df_real$y))
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1
pars           <- list(
  k1 = 8, k2 = 5, alpha = 0.5, beta = 1, mu_mu = 0
)

# Running the model ----------------------------------
hb_model <- hebart(formula,
                   data = train,
                   group_variable = group_variable,
                   num_trees = 5,
                   MCMC = list(iter = 500, burn = 150, thin = 1), 
                   k_1_pars = list(sample_k1 = FALSE, min_u = 0, max_u = 20,
                                   k1_prior = TRUE))
hb_model

# Hebart result
# -----------------------------------
#   Formula:
#   y ~ X1 
# 
# Number of trees:         50 
# Number of covariates:    1 
# Prediction error (MSE):  0.150179 
# R squared:               0.8402177 


# Making predictions ----------------------------------
pred_test <- predict_hebart(
  hb_model, newdata = test, 
  group_variable = group_variable, formula = formula
)

# Predictions in the test set --------------------------
data.frame(
  y = test$y, pred = pred_test$pred, 
  group = pred_test$group 
)  %>% 
  ggplot(aes(x = y, y = pred)) +
  geom_point(aes(colour = factor(group)), size = 2) +
  labs(x = "True y", y = "Predictions") +
  theme_light(16) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  guides(colour = "none")
