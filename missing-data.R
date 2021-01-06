
library(naniar)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

d1 = read.csv("session-rpe-single-time-point-data.csv")
d2 = read.csv("session-rpe-cycling-data.csv")
d3 = read.csv("session-rpe-pre-post-data.csv")

miss_1 <- d1 %>%
  select(USG, Temp, Condition) %>%
  group_by(Condition) %>%
  as.data.frame()

miss_2 <- d2 %>%
  select(HR, Condition) %>%
  group_by(Condition) %>%
  as.data.frame()

# USG, Temp
gg_miss_fct(x = miss_1, fct = Condition) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 20, face = "bold", hjust = -0.2)) +
  ylab("Variable") + ggtitle("(A)") -> plot1
plot1

miss_1 %>% select(-Condition) %>%
  gg_miss_var() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", hjust = -0.2)) + 
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, by = 1)) +
  ggtitle("(B)") -> plot2
plot2

# Missing data plot HR
gg_miss_fct(x = miss_2, fct = Condition) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 20, face = "bold", hjust = -0.2)) +
  ylab("Variable") + ggtitle("(C)") -> plot3
plot3

miss_2 %>% select(-Condition) %>%
  gg_miss_var() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 20, face = "bold", hjust = -0.2)) + 
  scale_y_continuous(limits = c(0,2), breaks = seq(0, 2, by = 1)) +
  ggtitle("(D)") -> plot4
plot4


# Cowplot
plot_grid(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2, align = 'v', axis = "lr")
ggsave(file = "supplementary-missing.png", units = "in", width = 8, height = 5, dpi = 300)


