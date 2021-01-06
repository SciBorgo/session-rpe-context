

#### Session-RPE study
# Author: JB, DB
# Date: August 2020

# Libraries
library(readr)
library(dplyr)
library(brms)
library(stringr)
library(emmeans)
library(tidybayes)
library(modelr)
library(ggplot2)
library(cowplot)
library(tidyverse)

# Data munging and model formula
sdat <- read_csv("session-rpe-single-time-point-data.csv")

sdat <- sdat %>%
  mutate(
    ID = as.factor(ID)
  )

# Merge post-cycling lactate scores
df <- read_csv("session-rpe-pre-post-data.csv")

diffd <- df %>% select(ID, Condition, Time, Lactate) %>%
  pivot_wider(names_from = Time, values_from = Lactate, names_prefix = "Time") %>%
  mutate(lactate_diff = Time2 - Time1)
diffd

df_lactate <- select(diffd, -Time1, -Time2)
sdat <- merge(sdat, df_lactate, by = c("ID","Condition"))

sdat$Lactate_s <- scale(sdat$lactate_diff, center = T, scale = T)
sdat$DALDA_B_s <- scale(sdat$DALDA_B, center = T, scale = T)

sdat <- sdat %>% 
  mutate(
    Condition = relevel(
      factor(Condition, 
             levels = 1:3, 
             labels = c("Written","Verbal","Group")), 
      ref = "Written")
  )

rpe_formula <- Session_RPE ~ DALDA_B_s + Lactate_s + Condition + (1|ID)

# Prior investigation
# default prior with flat prior for covariates
default_prior <- 
  get_prior(
    formula = rpe_formula,  
    family = cumulative(link = "logit"),
    data = sdat
    )

# replace flat prior for all covariates (i.e. not REs, not Condition) with (*):
default_covar_prior <- set_prior("normal(0, 1)", class = "b", coef = "")

updated_prior <- default_prior %>% 
  filter(!(coef == "" & str_detect(class, "b"))) %>% # remove old priors for Condition
  bind_rows(default_covar_prior) # add new priors

# replace default prior (as * above) for Condition with 
condition_prior <- set_prior("normal(0, 1)", class = "b", coef = paste0("Condition",levels(sdat$Condition)[-1]))

updated_prior <- updated_prior %>% 
  filter(!str_detect(coef, "Condition")) %>% # remove old priors for Condition
  bind_rows(condition_prior) # add new priors

md_prior_only_update <-  brm(
  formula = rpe_formula,
  family = cumulative(link = "logit"),
  prior = updated_prior,
  chains = 4, cores = 4, iter = 2000, 
  sample_prior = "only",
  data = sdat) 

# Prior odds ratio check:
emmeans(md_prior_only_update, pairwise~Condition, type="response")

# Prior predictive checks:
pp_check(md_prior_only_update, re_formula = NULL, nsamples = 50)
pp_check(md_prior_only_update, type = "bars_grouped", re_formula = NULL, group = "Condition", nsamples = 100)

# check prior predictive values (i.e. from prior only, not data)
sdat %>% 
  data_grid(Condition, ID) %>%
  mutate(Lactate_s = mean(sdat$Lactate_s), DALDA_B_s = mean(sdat$DALDA_B_s)) %>%
  add_fitted_draws(md_prior_only_update, 
                   re_formula = NULL, 
                   value = "Probability",
                   category = "Session_RPE",
                   ) %>%
  ggplot(aes(y = Probability, x = Session_RPE, fill = Condition)) +
  stat_halfeye(alpha = 0.3)

# Run model
md_s_rpe <- brm(
  formula = rpe_formula,
  family = cumulative(link = "logit"),
  prior = updated_prior,
  chains = 8, cores = 8, iter = 20000,
  control = list(adapt_delta = 0.95),
  seed = 123,
  data = sdat) 

# Posterior predictive checks
pp_check(md_s_rpe, re_formula = NULL, nsamples = 100)
pp_check(md_s_rpe, type = "bars_grouped", re_formula = NULL, group = "Condition", nsamples = 200)

# Results and Plots
# Summary
summary(md_s_rpe)

# Random effects
ranef(md_s_rpe)

# Effect of Condition
pd_95 <- conditional_effects(md_s_rpe, effects = "Condition", categorical = TRUE, plot = FALSE, prob = .95)[[1]] %>%
  as_tibble() %>% mutate(ci = "95%")
pd_66 <- conditional_effects(md_s_rpe, effects = "Condition", categorical = TRUE, plot = FALSE, prob = .66)[[1]] %>%
  as_tibble() %>% mutate(ci = "66%")

df <- bind_rows(pd_95, pd_66)

df %>%
  ggplot(aes(x = effect1__, y = estimate__, colour = cats__, group = cats__)) +
  geom_pointinterval(data = pd_66, aes(ymin = lower__, ymax = upper__),size = 5, position = position_dodge(width = 0.475)) +
  geom_interval(data = pd_95, aes(ymin = lower__, ymax = upper__),size = 0.5, position = position_dodge(width = 0.475)) +
  ylim(0, 1) + 
  ylab("Posterior Probability") + 
  xlab("Session-RPE Collection Context") +
  scale_colour_discrete(name = "s-RPE") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_segment(aes(x = 1, xend = 3, y = 0.875, yend = 0.875), colour = "gray50", size = 0.35) +
  geom_segment(aes(x = 1, xend = 1, y = 0.875, yend = 0.850), colour = "gray50", size = 0.35) +
  geom_segment(aes(x = 3, xend = 3, y = 0.875, yend = 0.850), colour = "gray50", size = 0.35) +
  annotate(geom = "text", x = 2, y = 0.925, label="OR = 5.31, 95% CrI = 1.62 to 17.64", size = 3.25)
ggsave("session-rpe.png", units = "in", width = 6, height = 3.5, dpi = 300)
ggsave("session-rpe.tiff", units = "in", width = 6, height = 3.5, dpi = 300, compression = "lzw")



# Predicted probabilities for each subject (after accounting for random effect)
# if each observation had just one condition (i.e. all observations used "written")
prob_gr <- function(p_x, p_y, equal = F){
  stopifnot(length(p_x) == length(p_y))
  n <- length(p_x)
  if(equal){
    term_collect <- outer(1:n, 1:n, ">=")
  } else {
    term_collect <- outer(1:n, 1:n, ">")
  }
  sum(outer(p_x,p_y,"*")[term_collect])
  
}

compare_groups_prob_gr <- function(probs_by, gr1, gr2, ...){
  sapply(1:nrow(sdat), function(i) prob_gr(probs_by[[gr1]][i,], probs_by[[gr2]][i,], ...)) %>%
    {c(mean = mean(.), sd = sd(.), quantile(., 0.025), quantile(., 0.975))}
}

compare_groups_prob_eq <- function(probs_by, gr1, gr2){

  compare_groups_prob_gr(probs_by, gr1, gr2, equal = TRUE) - compare_groups_prob_gr(probs_by, gr1, gr2)
  
}

probs_by <-
  c("group", "written", "verbel") %>%
  map( ~ mutate(sdat, Condition = .x)) %>%
  map( ~ predict(md_s_rpe, re_formula = NA, newdata = .x) ) %>%
  setNames(nm = c("group", "written", "verbel"))

# Pr(Group > Written) across participants
round(compare_groups_prob_gr(probs_by, "group", "written"),2)
# Pr(Group >= Written) across participants
round(compare_groups_prob_gr(probs_by, "group", "written", equal = TRUE),2)

# Pr(Group > Verbel) across participants
round(compare_groups_prob_gr(probs_by, "group", "verbel"),2)
# Pr(Group >= Verbel) across participants
round(compare_groups_prob_gr(probs_by, "group", "verbel", equal = TRUE),2)

# Pr(Written > Verbel) across participants
round(compare_groups_prob_gr(probs_by, "written", "verbel"),2)
# Pr(Written >= Verbel) across participants
round(compare_groups_prob_gr(probs_by, "written", "verbel", equal = TRUE),2)

# Pr(Written == Verbel) across participants
compare_groups_prob_eq(probs_by, "written", "verbel")
