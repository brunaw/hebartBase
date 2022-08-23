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
num_trees   <- 10

# Running the model ----------------------------------

hb_model <- hebart(y ~ X1,
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
                     scale_sigma_phi = 1
                   ), 
                   inits = list(tau = 1,
                                sigma_phi = 1),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1,
                               sigma_phi_sd = 0.5)
                   )
hb_model

pp <- predict_hebart(newX = test, new_groups = test$group,
                     hebart_posterior  = hb_model, 
                     type = "mean")

sqrt(mean((pp - test$y)^2)) # 27.6762
cor(pp, test$y)
qplot(test$y, pp) + geom_abline()
qplot(1:length(hb_model$sigma), hb_model$sigma)
qplot(1:length(hb_model$sigma), hb_model$sigma_phi)
qplot(1:length(hb_model$sigma), hb_model$log_lik)

stop()

# Comparison to BART --------------------------
bart_0 = dbarts::bart2(y ~ X1, 
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
pp <- predict(lme_ss, test)
sqrt(mean((pp - test$y)^2)) # 33.19528
cor(pp, test$y) # 0.8426536
qplot(test$y, pp) + geom_abline()




