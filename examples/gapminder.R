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
library(firatheme)
library(hebartBase)
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0) |> 
  set_names(c('X1', 'group', 'y', "X2", "X3"))


years       <- unique(df_real$X1)
to_remove   <- sample(years, 15)
train       <- df_real |> filter(!(X1 %in% to_remove))
test        <- df_real |> filter(X1 %in% to_remove)
groups      <- train$group
num_trees   <- 5

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1 + X2 + X3

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
                   MCMC = list(iter = 300, 
                               burn = 250, 
                               thin = 1,
                               sigma_phi_sd = 2)
                   )
pp <- predict_hebart(test, test$group, hb_model, type = "mean")
sqrt(mean((pp - test$y)^2)) # 1.907164
cor(pp, scale(test$y))  # 0.9881025
pp_train <- predict_hebart(train, train$group, hb_model, type = "mean")
sqrt(mean(pp_train - train$y)^2) # 0.0001694319
diagnostics(hb_model)

# Comparison to BART --------------------------
bart_0 = dbarts::bart2(y ~ X1 + X2 + X3, 
                       data = train,
                       test = test,
                       keepTrees = TRUE)
ppbart <- bart_0$yhat.test.mean
sqrt(mean((ppbart - test$y)^2)) # 7.944524
cor(ppbart, test$y) #    0.698455

ppbart <- bart_0$yhat.train.mean
sqrt(mean((ppbart - train$y)^2)) # 0.8950348- 100 trees

# BART+Group
bart_0 = dbarts::bart2(y ~ X1 + X2 + X3 + group, 
                       data = train,
                       test = test,
                       keepTrees = TRUE)
ppbart <- bart_0$yhat.test.mean
sqrt(mean((ppbart - test$y)^2)) # 0.9425852
cor(ppbart, test$y) #    0.99683

ppbart <- bart_0$yhat.train.mean
sqrt(mean((ppbart - train$y)^2)) # 0.3275252

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + X2 + X3 + (1|group), train)
pplme <- predict(lme_ss, test)
sqrt(mean((pplme - test$y)^2)) # 3.991818
cor(pplme, test$y) # 0.936175

# Average predictions 
preds_y <- data.frame(test, pred = pp, pred_lme = pplme)

preds_y |> 
  filter(group %in% c("China", "South Africa", "Russia", 
                      "United States", "Mexico", 
                      "Canada")) |> 
  ggplot(aes(x = X1, y = y)) +
  geom_point(colour = "gray") +
  geom_line(aes(y = pred, colour= "#75E6DA"), size = 1.3) +
  geom_line(aes(y = pred), colour= "#75E6DA", size = 1,
            alpha = 0.7) +
  geom_line(aes(y = pred_lme, colour= "#F96209"),  size = 1.3) +
  geom_line(aes(y = pred_lme), colour= "#F96209",  size = 1, 
            alpha = 0.7) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = "Life expectancy", 
       y = 'Predictions', 
       #title = paste0("RMSE\nHE-BART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  facet_wrap(~group, scales = "free", ncol = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(14) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="gray", 
             `HEBART Prediction`="#75E6DA", 
             `LME Prediction`= '#F96209'), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3, 3), shape = c(16, 16, 16)))) + 
  theme_fira()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------




