# Create a really simple example that HEBART should do really well on

# Just one tree with one covariate and 3 groups
# Terminal node 1 comes from x <= 0.5, has mu = -5, phi 1 2 3 = -7 -5 -3
# Node 2 comes from x > 0.5 and x <= 0.75, has mu = 10, phi = 8 10 12
# Node 3 comes from x > 0.75 has mu  0 and phi = -2 0 2

# Clear workspace and load in packages
library(tidyverse)
devtools::load_all(".")

# Create this data set
set.seed(123)
n <- 400
x <- runif(n)
groups <- factor(rep(1:3, length = n))
mu <- rep(0, n)
mu[x<0.75] <- 10
mu[x<0.5] <- -5
# qplot(x, mu)
phi <- rep(NA, length = n)
phi_tn1 <- c(-7, -5, -3)
phi[mu == -5] <- phi_tn1[groups[mu == -5]]
phi_tn2 <- c(8, 10, 12)
phi[mu == 10] <- phi_tn2[groups[mu == 10]]
phi_tn3 <- c(-2, 0, 2)
phi[mu == 0] <- phi_tn3[groups[mu == 0]]
# qplot(x, phi, colour = groups)
res_sd <- 1
y <- rnorm(n, phi, res_sd)
qplot(x, y, colour = groups)

train <- data.frame(
  y = y[1:200], 
  x = x[1:200],
  groups = groups[1:200]
)
test <- data.frame(
  y = y[201:400], 
  x = x[201:400],
  groups = groups[201:400]
)


# Now run HEBART on it ----------------------------------------------------
num_trees <- 10
hb_model <- hebart(y ~ x,
                   data = train,
                   group_variable = "groups", 
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
                                sigma_phi = 0.01),
                   MCMC = list(iter = 1500, 
                               burn = 500, 
                               thin = 1,
                               sigma_phi_sd = 1)
)
hb_model

pp <- predict_hebart(newX = data.frame(x = test$x), 
                     new_groups = test$group,
                     hebart_posterior  = hb_model, 
                     type = "mean")

sqrt(mean((pp - test$y)^2))
cor(pp, test$y)
qplot(test$y, pp) + geom_abline()
qplot(1:length(hb_model$sigma), hb_model$sigma)
qplot(1:length(hb_model$sigma), hb_model$sigma_phi)
