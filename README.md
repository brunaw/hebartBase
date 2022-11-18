
**Author:** [Bruna Wundervald](http://brunaw.com/) **License:** [MIT](https://opensource.org/licenses/MIT)<br/>

  [![CRAN status](http://www.r-pkg.org/badges/version/hebart)](https://cran.r-project.org/package=hebart) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/hebart)](https://cran.r-project.org/package=hebart)  [![Travis-CI Build Status](https://travis-ci.org/brunaw/hebart.svg?branch=master)](https://travis-ci.org/brunaw/hebart)


`hebart`: A package for fitting of Hierachical Embedded Bayesian Additive Regression Trees, which is an extension of BART to grouped data. 
========================================

hebart description goes here

  Installation
------------------------


You can currently install the latest version of `hebart` from the github repository with:
  ``` r
# install.packages("devtools")
devtools::install_github("brunaw/hebartBase")
```

This package wil soon be on CRAN; 

``` r
install.packages("hebart")
```

Functionalities
------------------------


All of the functions and documentation can be found with:

``` r
library(hebart)
packageVersion("hebart")
ls("package:hebart")
help(package = "hebart")
```

Examples 
------------------------

``` r
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
                   MCMC = list(iter = 500, 
                               burn = 250, 
                               thin = 1,
                               sigma_phi_sd = 0.5)
)

hb_model

Hebart result
-----------------------------------
Formula:
 y ~ X1 

Number of trees:         10 
Number of covariates:    1 
Training error (RMSE):    0.3851619 
R Squared:                0.842164 


# Making predictions ----------------------------------
pred_test <- predict_hebart(newX = test, new_groups = test$group,
                            hebart_posterior  = hb_model, type = "mean")

# Predictions in the test set --------------------------
data.frame(
  y = test$y, pred = pred_test, 
  group = test$group 
)  %>% 
  ggplot(aes(x = y, y = pred)) +
  geom_point(aes(colour = factor(group)), size = 2) +
  labs(x = "True y", y = "Predictions") +
  theme_light(16) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  guides(colour = "none")
 
```

<img src="img/predictions.png" class="center">

```{r}
# Diagnostics --------------------------
diagnostics(hb_model)
#----------------------------------------------------
```
<img src="img/diagnostics.png" class="center">


To be implemented
------------------------


- [ ] The official `pkgdown` for `hebart`
- [ ] Variable importances


Citation
------------------------

To cite this package in publications, please use:

```
Bruna Wundervald (2022). hebartBase: Hierachical Embedded Bayesian Additive Regression Trees. R package version 0.1.0.
https://CRAN.R-project.org/package=hebartBase 
```

A BibTeX entry for LaTeX users is

```
@Manual{,
  title = {hebartBase: Hierachical Embedded Bayesian Additive Regression Trees},
  author = {Bruna Wundervald},
  year = {2022},
  note = {R package version 0.1.0},
  url = {https://CRAN.R-project.org/package=hebartBase},
}
```

This citation format can be obtained at any moment in `R` with:

  ``` r
citation('hebart')
```

Contributing
------------------------

  Contributions to this project are always highly incentivized. To do
so, please be aware that `git` is our main tool for version control.
The minimal steps for a contribution are:

  1. Fork this repository into your `GitHub` account and clone it
the way you prefer.
2. Do the changes, making sure everything is well documented,
examples are provided and checking if the package still correctly builds.
3. Push your changes to `git` and create a new pull request in
`GitHub`, explaining why and what are the changes made.
4. Done! Wait for review & acceptance of the pull request :)

To contributors who are new to writing R packages, we recommend
the ['R Packages' book](http://r-pkgs.had.co.nz/), by Hadley
Wickham. To those who are new to `git`/`GitHub`, we recommend
[this tutorial](http://brunaw.com/talk/git/). Many contributing
resources to open source projects can be found at
[this repository](https://github.com/freeCodeCamp/how-to-contribute-to-open-source).



