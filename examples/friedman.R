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

X <- data$X
y <- data$y
X_test <- test$X
y_test <- test$y
group <- sample(1:3, size = length(y), replace = TRUE)
group_test <- sample(1:3, size = length(y_test), replace = TRUE)


train <- data.frame(y, X, group)
test  <- data.frame(y = y_test, X_test, group = group_test)

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1 + X2 + X3 + X4 + X5
num_trees <- 4

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
                     shape_sigma_phi = 0.5,
                     scale_sigma_phi = 1,
                     sample_sigma_phi = TRUE
                   ), 
                   inits = list(tau = 1,
                                sigma_phi = 1),
                   MCMC = list(iter = 1500, 
                               burn = 250, 
                               thin = 1,
                               sigma_phi_sd = 2)
                   )
pp <- predict_hebart(test, test$group, hb_model, type = "mean")
sqrt(mean(pp - test$y)^2)
cor(pp, scale(test$y))
qplot(test$y, pp) + geom_abline()
#-----------------------------------------------------------------------
# There is no need to make comparisons here; 
# The group is not actually relevant to the modelling
#-----------------------------------------------------------------------