#setwd("C:/Users/jy33/OneDrive/Desktop/R/winner_loser")

library(lme4)
library(tidyverse)
library(metafor)
library(ggplot2); theme_set(theme_classic())
library(emmeans)
library(car)
library(ggsci)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text.y = element_text(size = 16))

all_data <- read.csv("win_meta_2023.csv", stringsAsFactors = TRUE) %>% 
            filter(exclude != "exclude") # excluding one study that we realized doesn't qualify 

# Use escalc function to get effect sizes (in log-odds) and their variances for each comparison
all_data <- escalc(data = all_data, xi = ifelse(contest_outcome == "winner", wins, losses), 
                   ni = sample_size, measure = "PLO", append = TRUE)

############################## MAIN MODEL #############################
# Generate rma model for winner_loser effects for all studies
met_mod_all <-rma.mv(data = all_data, yi, vi, mods = ~ contest_outcome + protocol, 
                     random = ~ 1 | study/experiment, method="ML")
summary(met_mod_all)

# Get QM values for each moderator
anova(met_mod_all, btt = 2) # contest outcome (winner or loser)
anova(met_mod_all, btt = 3) # protocol type (random or self)

# Obtain 95% CI for winner and loser effects for plotting
em_mod_all <- qdrg(object = met_mod_all, data = all_data)
win_loss_95_CI <- as.data.frame(emmeans(em_mod_all, specs = ~ contest_outcome))

emmeans(em_mod_all, specs = ~ contest_outcome)

# Plot of effect sizes for winner and loser effects
ggplot(data = all_data, aes(x = yi, y = contest_outcome, color = contest_outcome)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.4, size = all_data$sample_size/6) + 
  My_Theme + ylab("") + xlab("Effect size (Log odds)") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) + 
  geom_vline(xintercept = 0, linetype = 2) + scale_color_manual(values = c("#0077b6", "#d00000")) +
  geom_pointrange(data = win_loss_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                  lwd = 2.5, fatten = 2, shape = 23)

# Obtain 95% CI for random vs. self selection for plotting
prot_95_CI <- as.data.frame(emmeans(em_mod_all, specs = ~ protocol))

# Plot of effect sizes for self vs. random selection
ggplot(data = all_data, aes(x = yi, y = protocol, color = protocol)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.4, size = all_data$sample_size/6) + 
  My_Theme + ylab("") + xlab("Effect size (Log odds)") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) + 
  geom_vline(xintercept = 0, linetype = 2) + scale_color_manual(values = c("#4c956c", "#f8961e")) +
  geom_pointrange(data = prot_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                  lwd = 2.5, fatten = 2, shape = 23) 

####################### TAXANOMIC GROUP MODEL #####################
met_mod_animal <-rma.mv(data = all_data, yi, vi, mods = ~ 0 + class, random = ~ 1 | study/experiment, method="ML")
summary(met_mod_animal)

# Another way to get 95% CI for winner and loser effects
em_mod_animal <- qdrg(object = met_mod_animal, data = all_data)
animal_95_CI <- as.data.frame(emmeans(em_mod_animal, specs = ~ class))
emmeans(em_mod_animal, specs = ~ class)

ggplot(data = all_data, aes(x = yi, y = class, color = class)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.3, size = all_data$sample_size/6) + 
  My_Theme + ylab("") + xlab("Effect size (Log odds)") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) +
  geom_vline(xintercept = 0, linetype = 2) + scale_color_manual(values = c("#353535","#ff7900","#00a5cf","#73a942","#724cf9","red")) +
  geom_pointrange(data = animal_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                lwd = 2, fatten = 2, shape = 23) 

############################ PROTOCOL TYPE AS A FUNCTION OF STUDY PUB DATE ##################################
study_info <- read.csv("study_info_2023.csv", stringsAsFactors = TRUE)

year_mod <- glm(data = study_info, protocol ~ year, family = binomial())
summary(year_mod)
confint.default(year_mod)

study_info$protocol <- as.numeric(study_info$protocol) - 1

Anova(year_mod, test.statistic = "Wald")

study_info$protocol <- as.numeric(study_info$protocol) - 1 # Turn study type into numeric variable for plotting

ggplot(data = study_info, aes(y = protocol, x = year)) + 
  geom_smooth(method = "glm", method.args = list(family = binomial), size = 1.5, color = "orangered3") + 
  geom_jitter(height = 0, width = 0.2, alpha = 0.4, size = 3) + 
  ylab("Frequency of studies using \n random assignment") + xlab("Year") + 
  My_Theme

## Breakdown of studies by animal group
animal_dat <- read.csv("data/counts_by_animals.csv")

ggplot(data = animal_dat, aes(x = frequency, y = animal_group, fill = animal_group, color = animal_group)) + 
       geom_bar(stat = "identity", alpha = 0.4) + theme(legend.position = "none") +
       xlab("Number of studies") + scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10 ,12)) + 
       ylab("") + My_Theme + scale_fill_nejm() + scale_color_nejm()


############# PUBLICATION BIAS #######################
with(met_mod_all, funnel(yi, vi))


########  EXTRA ANALYSES   ########### 
######## ONLY RANDOM MODEL ###########
random_dat <- all_data %>% 
  filter(protocol == "random")

met_mod_random <-rma.mv(data = random_dat, yi, vi, mods = ~ contest_outcome, random = ~ (1|study), method="ML")
summary(met_mod_random)

# Another way to get 95% CI for winner and loser effects
em_mod_random <- qdrg(object = met_mod_random, data = random_dat)
random_95_CI <- as.data.frame(emmeans(em_mod_random, specs = ~ contest_outcome))

ggplot(data = random_dat, aes(x = yi, y = contest_outcome, color = contest_outcome)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.4, size = random_dat$sample_size/6) + 
  My_Theme + ylab("") + xlab("Effect size (Log odds)") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) + 
  geom_vline(xintercept = 0, linetype = 2) + scale_fill_nejm() + scale_color_nejm() + 
  geom_pointrange(data = random_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                  lwd = 2.5, fatten = 2, shape = 23) 

######## ONLY SELF MODEL ###########
self_dat <- all_data %>% 
  filter(protocol == "self")

met_mod_self <-rma.mv(data = self_dat, yi, vi, mods = ~contest_outcome, random = ~ (1|study), method="ML")
summary(met_mod_self)

# Another way to get 95% CI for winner and loser effects
em_mod_self <- qdrg(object = met_mod_self, data = self_dat)
self_95_CI <- as.data.frame(emmeans(em_mod_self, specs = ~ contest_outcome))

ggplot(data = self_dat, aes(x = yi, y = contest_outcome, color = contest_outcome)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.4, size = self_dat$sample_size/6) + 
  My_Theme + ylab("") + xlab("Effect size (Log odds ratio)") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4)) + 
  geom_vline(xintercept = 0, linetype = 2) + scale_fill_nejm() + scale_color_nejm() +
  geom_pointrange(data = self_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                  lwd = 2.5, fatten = 2, shape = 23) 


## TEST 
ggplot(data = all_data, aes(x = yi, y = contest_outcome, color = contest_outcome)) + 
  geom_jitter(height = 0.25, width = 0, alpha = 0.4, size = all_data$sample_size/6) + 
  My_Theme + ylab("") + xlab("Odds of repeating prior outcome") + theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = c(-2.303, 0, 2.303), 
                                        label = c("0.1", "1", "10")) + 
  geom_vline(xintercept = 0, linetype = 2) + scale_color_manual(values = c("#0077b6", "#d00000")) +
  geom_pointrange(data = win_loss_95_CI, aes(x = emmean, xmin = asymp.LCL, xmax = asymp.UCL), color = "grey20", 
                  lwd = 2.5, fatten = 2, shape = 23)


######### TESTING FOR PUBLICATION BIAS ##########
met_mod_uni <- rma.uni(data = all_data, yi, vi, method="ML")
summary(met_mod_uni)

# Funnel plots
funnel(met_mod_uni, xlab = "Log odds ratio") # basic funnel plot
funnel(met_mod_uni, yaxis="seinv", xlab="Effect size (log odds ratio)", 
       back="white", col=rgb(0,153,76, max=255, alpha=125), digits = 1) # nicer funnel plot

# Trim and fill test
trimfill_left <- trimfill(met_mod_uni, side="left", estimator = "L0")
summary(trimfill_left)
funnel(trimfill_left, yaxis="seinv", xlab="Effect size (Log odds ratio)", digits=1,
       back="white", col=rgb(0,153,76, max=255, alpha=100))
abline(v = 0.82, lty = 2, col = "forestgreen")


regtest(met_mod_uni)

year_bias_test <- lmer(data = all_data, yi ~ year + (1|study))
Anova(year_bias_test)
