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
devtools::load_all(".")
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)
df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0) |> 
  set_names(c('X1', 'group', 'y', "X2", "X3"))

dim(df_real)
data_split  <- initial_split(df_real)
train       <- training(data_split)
test        <- testing(data_split)
groups      <- train$group
num_trees   <- 10

# Model parameters -----------------------------------
group_variable <-  "group"
formula        <- y ~ X1 + X2 + X3

pars   <- list(
  alpha = 0.95, beta = 2,
  nu = 3, lambda = 0.1,
  tau_mu = 16 * num_trees,
  shape_sigma_phi = 0.5,
  scale_sigma_phi = 1 # These give mean ~2 and sd ~4
)

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
                     scale_sigma_phi = 1
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
qplot(1:length(hb_model$sigma), hb_model$sigma)
qplot(1:length(hb_model$sigma), hb_model$sigma_phi)
qplot(test$y, pp) + geom_abline()

# Comparison to BART --------------------------
bart_0 = dbarts::bart2(y ~ X1 + X2 + X3 + group, 
                       #n.trees = 15,
                       data = train,
                       test = test,
                       keepTrees = TRUE)
ppbart <- bart_0$yhat.test.mean
sqrt(mean((ppbart - test$y)^2)) # 53.28167 - 100 trees
cor(ppbart, test$y) #   0.4983258
qplot(test$y, ppbart) + geom_abline()

# Comparison to LME --------------------------
lme_ss <- lme4::lmer(y ~ X1 + X2 + X3 + (1|group), train)
pplme <- predict(lme_ss, test)
sqrt(mean((pplme - test$y)^2)) # 33.19528
cor(pplme, test$y) # 0.8426536
qplot(test$y, pplme) + geom_abline()


# Average predictions 

preds_y <- data.frame(test, pred = pp, pred_lme = pplme)

# preds_y |> 
#   group_by(group) |> 
#   summarise(
#     m = mean((y - pred)^2)
#   )

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
# theme(panel.spacing.x = unit(0.5, "lines"), 
#       legend.position = "bottom")

ggsave("results/pred_gapminder.png", width = 7, height = 8)
# ------------------------------------------------------------------------------




