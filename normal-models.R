

#### Session-RPE study
# Author: JB, DB
# Date: 26 Nov 2020

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
library(tidyr)


#### Heart rate
# Data munging and model formula

sdat <- read_csv("session-rpe-cycling-data.csv") 

sdat <- sdat %>%
  mutate(
    Condition = as.factor(Condition),
    ID = as.factor(ID),
    Interval = as.factor(Interval)
  )

sdat <- sdat %>% 
  mutate(
    Condition = relevel(
      factor(Condition, 
             levels = 1:3, 
             labels = c("Written","Verbal","Group")), 
      ref = "Written")
  )

# remove missing and potential bad data:
sdat <- 
  sdat %>% filter(
  !is.na(HR),
  ! (ID == "7" & Condition == "Verbal")
  )

sdat %>% ggplot() + geom_line(aes(x=Time, y = HR, colour = Condition)) + facet_wrap(~ID)

hr_formula <- bf(HR ~ s(Time, bs = "cr", k = 5) + Condition + (1|Interval) + (1|ID))

# Prior investigation
# default prior with flat prior for covariates
default_prior <- 
  get_prior(
    formula = hr_formula,  
    family = gaussian(),
    data = sdat
  )

# replace flat prior for all covariates (i.e. not REs, not Condition) with (*):
default_covar_prior <- set_prior("normal(0, 20)", class = "b", coef = "")

updated_prior <- default_prior %>% 
  filter(!(coef == "" & str_detect(class, "b"))) %>% # remove old priors for Condition
  bind_rows(default_covar_prior) # add new priors

# replace default prior (as * above) for Condition with 
condition_prior <- set_prior("normal(0, 20)", class = "b", coef = paste0("Condition",levels(sdat$Condition)[-1]))

updated_prior <- updated_prior %>% 
  filter(!str_detect(coef, "Condition")) %>% # remove old priors for Condition
  bind_rows(condition_prior) # add new priors

sd_prior <- set_prior("student_t(3, 0, 10)", class = "sd")

updated_prior <- updated_prior %>%
  filter(!str_detect(class, "sd")) %>% # remove old priors for Condition
  bind_rows(sd_prior) # add new priors


# Prior Cohen'd d check:
tb_Condition2 <- 
  gather_draws(md_prior_only_update, b_Condition2) %>%
  mutate(cohen = .value/sd(.value))
tb_Condition2 %>% ggplot() + geom_histogram(aes(x = cohen))
tb_Condition2 %>% mean_qi(cohen, .width = c(0.5, .95))

tb_Condition3 <- 
  gather_draws(md_prior_only_update, b_Condition3) %>%
  mutate(cohen = .value/sd(.value))
tb_Condition3 %>% ggplot() + geom_histogram(aes(x = cohen))
tb_Condition3 %>% mean_qi(cohen, .width = c(0.5, .95))

tb_Condition23 <- 
  gather_draws(md_prior_only_update, b_Condition2, b_Condition3) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(cohen = (b_Condition3 - b_Condition2) / sd(b_Condition3 - b_Condition2) )
tb_Condition23 %>% ggplot() + geom_histogram(aes(x = cohen))
tb_Condition23 %>% mean_qi(cohen, .width = c(0.5, .95))

pp_check(md_prior_only_update, re_formula = NULL, nsamples = 50)

pp_check(md_prior_only_update, type = 'ribbon_grouped', x= "Time", group = "Condition", nsamples = 100)
# the "strange" lines come from having multiple individuals


# Run model
md_s_hr <- brm(
  formula = hr_formula,
  family = gaussian(),
  prior = updated_prior,
  chains = 4, cores = 8, iter = 10000, 
  data = sdat, 
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  seed = 123) 

# Posterior predictive checks
pp_check(md_s_hr, re_formula = NULL, nsamples = 50)
pp_check(md_s_hr, type = 'ribbon_grouped', x= "Time", group = "Condition", nsamples = 100)

# Look at priors
prior_summary(md_s_hr)

# Summary
summary(md_s_hr)

# Plot
time_by_interval <- sdat %>% select(Time, Interval) %>% unique()

levels(sdat$Condition)

p <- sdat %>%
  group_by(ID) %>%
  data_grid(Condition, Time) %>%
  left_join(time_by_interval) %>%
  add_fitted_draws(md_s_hr) %>%
  ggplot(aes(x = Time, y = HR, colour = Condition, group = Condition)) +
  stat_pointinterval(aes(y = .value), .width = c(0.50,0.95), position = position_dodge(width = 1)) +
  #scale_colour_grey() +
  scale_colour_manual(values = c("#90b8da","#10385a","#2171b5")) +
  scale_x_continuous(breaks = seq(0,24, by = 3)) +
  scale_y_continuous(limits = c(57,195), breaks = seq(60,195, by = 15)) +
  ylab(Heart~Rate~(b~min^{-1})) +
  xlab("Time (min)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p + geom_segment(aes(x = 0, xend = 24, y = 185, yend = 185), colour = "gray50", size = 0.25) +
  annotate(geom = "text", x = 12, y = 187, label = "*", size = 6)

ggsave(file = "session-rpe-hr.png", units = "in", width = 7, height = 4, dpi = 300)
ggsave(file = "session-rpe-hr.tiff", units = "in", width = 7, height = 4, dpi = 300, compression = "lzw")




#### USG
d <- read_csv("session-rpe-single-time-point-data.csv") 

d <- d %>%
  mutate(
    Condition = as.factor(Condition),
    ID = as.factor(ID)
    )

d %>% ggplot() + 
  geom_boxplot(aes(x=Condition, y = USG)) +
  geom_jitter(aes(x=Condition, y = USG), width = 0.1, size = 3)

fit <- brm(
  formula = log(USG) ~ Condition + (1|ID),
  chains = 4, cores = 4, iter = 20000, 
  data = d,
  seed = 123) 

# Predictive check
pp_check(fit, re_formula = NULL, nsamples = 100)

# Summary
summary(fit)

# Pairwise
fit_usg <- emmeans::emmeans(fit, pairwise ~ Condition)
summary(fit_usg, "Condition", transform = "response") 






#### Baseline nude body mass
d <- read_csv("session-rpe-pre-post-data.csv") 

d <- d %>%
  mutate(
    Condition = as.factor(Condition),
    ID = as.factor(ID),
    Time = as.factor(Time)
  )

d %>% ggplot() + 
  geom_boxplot(aes(x=Condition, y = log(Mass))) +
  geom_jitter(aes(x=Condition, y = log(Mass)), width = 0.1, size = 3, alpha = 0.25) +
  facet_grid(~Time)

d %>% group_by(Time, Condition) %>%
  summarise(mu = mean(Mass),
            median = median(Mass),
            sd = sd(Mass))

dsub <- d %>% filter(Time == "1")

condition_prior <- set_prior("normal(0, 20)", class = "b", coef = paste0("Condition",2:3))

# Model
fit <- brm(
  formula = Mass ~ Condition + (1|ID),
  chains = 4, cores = 4, iter = 20000,
  prior = condition_prior,
  data = dsub,
  seed = 123) 

# Priors
prior_summary(fit)

# Predictive check
pp_check(fit, re_formula = NULL, nsamples = 100)

# Summary
summary(fit)

# Pairwise
fit_bm <- emmeans::emmeans(fit, pairwise ~ Condition)
summary(fit_bm, "Condition") 

