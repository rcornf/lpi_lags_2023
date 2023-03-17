###%%%%%%%%%%
### SI Analysis of lag results
### Cross-validation stability
### Effect of removing short population time-series (DP_min thresholds)
### Ecological subsetting
### Alternative environmental datasets
###%%%%%%%%%%

# Clear env
rm(list = ls())
graphics.off()

# Load packages
library(reshape2)       #~ Data wrangling
library(plyr)           #~ Data wrangling
library(lme4)           #~ LMMs
library(MuMIn)
library(ggplot2)        #~ Plotting....
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)
library(ggpubr)
library(dplyr)


###%%%%%%%%%%
# Main Code
###%%%%%%%%%%

source("00_analysis_functions.R")

# Initially, just get main results data
# lapply outut dir, read rds
lag_res <- lapply(list.files("../Results/main_models",
                             full.names = T), 
                  FUN = readRDS)
# bind together
lag_res <- do.call(rbind.fill, lag_res)

# Run main analysis
fig_ls <- analyses_wrap_(lag_res)
rm(lag_res)


##%%%%%
# Figures 
##%%%%%

##%%%%%
## CV coef results
##%%%%%
# load cv data
# _test
loo_res <- lapply(list.files("../Results/cv_models",
                             full.names = T), 
                  FUN = readRDS)

loo_df <- bind_rows(lapply(loo_res,
                           function(x){x$loo}))

loobsr_df <- bind_rows(lapply(loo_res,
                              function(x){x$loobsr}))
rm(loo_res)

# Summarise - iqr+, 95 quants, entire range
loo_summ_df <- loo_df %>%
    dplyr::group_by(model, class, type, cc_col, luc_col, dp_min, dur, off, Res, delta, coef_nm) %>%
    dplyr::summarise(q0 = min(coef_val),
                     q2.5 = quantile(coef_val, 0.025),
                     q25 = quantile(coef_val, 0.25),
                     q50 = quantile(coef_val, 0.5),
                     q75 = quantile(coef_val, 0.75),
                     q97.5 = quantile(coef_val, 0.975),
                     q100 = max(coef_val)
    ) %>%
    ungroup() %>%
    mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                               coef_nm == "cc_s" ~ "CC",
                               coef_nm == "luc_s" ~ "LUC",
                               coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                               coef_nm == "bm_s" ~ "BM",
                               coef_nm == "paYes" ~ "PA",
                               coef_nm == "Managed1" ~ "Man",
                               coef_nm == "Managed2" ~ "Man2",
                               coef_nm == "Utilised1" ~ "Use",
                               coef_nm == "Utilised2" ~ "Use2")) %>%
    mutate(coef_nm = factor(coef_nm,
                            levels = c("Int.", "CC", "LUC", "CC:LUC",
                                       "BM", "PA", "Man", "Use", "Man2", "Use2"))) 
####

loobsr_summ_df <- loobsr_df %>%
    dplyr::group_by(model, class, type, cc_col, luc_col, dp_min, dur, off, Res, delta, coef_nm) %>%
    # summarise(mean_ = mean(coef_val))
    dplyr::summarise(q0 = min(coef_val),
                     q2.5 = quantile(coef_val, 0.025),
                     q25 = quantile(coef_val, 0.25),
                     q50 = quantile(coef_val, 0.5),
                     q75 = quantile(coef_val, 0.75),
                     q97.5 = quantile(coef_val, 0.975),
                     q100 = max(coef_val)
    ) %>%
    ungroup() %>%
    mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                               coef_nm == "cc_s" ~ "CC",
                               coef_nm == "luc_s" ~ "LUC",
                               coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                               coef_nm == "bm_s" ~ "BM",
                               coef_nm == "paYes" ~ "PA",
                               coef_nm == "Managed1" ~ "Man",
                               coef_nm == "Managed2" ~ "Man2",
                               coef_nm == "Utilised1" ~ "Use",
                               coef_nm == "Utilised2" ~ "Use2")) %>%
    mutate(coef_nm = factor(coef_nm,
                            levels = c("Int.", "CC", "LUC", "CC:LUC",
                                       "BM", "PA", "Man", "Use", "Man2", "Use2"))) 


# plot loo cv/loobsr cv - stability
ggplot(subset(loo_summ_df, model == "m1b" & !(coef_nm %in% c("Man2", "Use2")))) +
    geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
    geom_hline(aes(yintercept = -3.50389 , 
                   alpha = "yellow"), 
               colour = "yellow2") +
    geom_hline(aes(yintercept =-6.696701 , 
                   alpha = "orange"), 
               colour = "orange") +
    
    geom_errorbar(aes(x = coef_nm, 
                      ymin = (10^q0 - 1) * 100,
                      ymax = (10^q100 - 1) * 100,
                      colour = delta, group = delta),
                  width = 0,
                  size = 1,
                  position = position_dodge(0.5)) +
    geom_point(aes(x = coef_nm, y = coef_val.1,
                   group = delta),
               colour = "firebrick3",
               size = 1.5,
               position = position_dodge(0.5),
               data = subset(fig_ls$coef_df,
                             model == "m1b" & !(coef_nm %in% c("Managed2", "Utilised2")) &
                                 type == "year" & delta < 6) %>%
                   mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                                              coef_nm == "cc_s" ~ "CC",
                                              coef_nm == "luc_s" ~ "LUC",
                                              coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                                              coef_nm == "bm_s" ~ "BM",
                                              coef_nm == "paYes" ~ "PA",
                                              coef_nm == "Managed1" ~ "Man",
                                              coef_nm == "Utilised1" ~ "Use")) %>%
                   mutate(coef_nm = factor(coef_nm,
                                           levels = c("Int.", "CC", "LUC", "CC:LUC",
                                                      "BM", "PA", "Man", "Use")))
    ) +
    scale_colour_gradient(low = "black",
                          high = "gainsboro",
                          limits = c(0,6),
                          breaks = c(0,3,6)) +
    labs(x = element_blank(),
         y = "Annual pop. change (%)",
         colour = expression(Delta~"AICc"),
         alpha = "IUCN A2\nthreshold") +
    scale_alpha_manual(values = c(1,1,1),
                       breaks = c("yellow","orange"),
                       labels = c("VU","EN")) +
    guides(alpha = guide_legend(override.aes= list(colour = c("yellow2",
                                                              "orange"),
                                                   size = 4),
                                direction = "vertical",
                                # title.position = "left",
                                title.vjust = 1)) +
    facet_grid(class~.,
               labeller = labeller("class" = c("bird" = "Birds",
                                               "mammal" = "Mammals"))) +
    theme_bw() +
    basic_thm +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1))

ggsave("../Results/Figs/sm_fig_2.pdf",
       device = "pdf", dpi = 300,
       width = 5.5, height = 7.5)


ggplot(subset(loobsr_summ_df, model == "m1b" & !(coef_nm %in% c("Man2", "Use2")) )) +
    geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
    geom_hline(aes(yintercept = -3.50389 , 
                   alpha = "yellow"), 
               colour = "yellow2") +
    geom_hline(aes(yintercept =-6.696701 , 
                   alpha = "orange"), 
               colour = "orange") +
    
    geom_errorbar(aes(x = coef_nm, 
                      ymin = (10^q0 - 1) * 100,
                      ymax = (10^q100 - 1) * 100,
                      colour = delta, group = delta),
                  width = 0,
                  size = 1,
                  position = position_dodge(0.5)) +
    
    geom_point(aes(x = coef_nm, y = coef_val.1,
                   group = delta),
               colour = "firebrick3",
               size = 1.5,
               position = position_dodge(0.5),
               data = subset(fig_ls$coef_df,
                             model == "m1b" & !(coef_nm %in% c("Managed2", "Utilised2")) &
                                 type == "year" & delta < 6) %>%
                   mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                                              coef_nm == "cc_s" ~ "CC",
                                              coef_nm == "luc_s" ~ "LUC",
                                              coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                                              coef_nm == "bm_s" ~ "BM",
                                              coef_nm == "paYes" ~ "PA",
                                              coef_nm == "Managed1" ~ "Man",
                                              coef_nm == "Utilised1" ~ "Use")) %>%
                   mutate(coef_nm = factor(coef_nm,
                                           levels = c("Int.", "CC", "LUC", "CC:LUC",
                                                      "BM", "PA", "Man", "Use")))
    ) +
    scale_colour_gradient(low = "black",
                          high = "gainsboro",
                          limits = c(0,6),
                          breaks = c(0,3,6)) +
    labs(x = element_blank(),
         y = "Annual pop. change (%)",
         colour = expression(Delta~"AICc"),
         alpha = "IUCN A2\nthreshold") +
    scale_alpha_manual(values = c(1,1,1),
                       breaks = c("yellow","orange"),
                       labels = c("VU","EN")) +
    guides(alpha = guide_legend(override.aes= list(colour = c("yellow2",
                                                              "orange"),
                                                   size = 4),
                                direction = "vertical",
                                # title.position = "left",
                                title.vjust = 1)) +
    facet_grid(class~.,
               labeller = labeller("class" = c("bird" = "Birds",
                                               "mammal" = "Mammals"))) +
    theme_bw() +
    basic_thm +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1))


# In res2 m1b and 2b yr, there are no pops/spp with Cook's D > 0.5 - 2b has loc with D>0.5...


##%%%%%
# DP min effects - lags, coefs
##%%%%%
# load data
# lapply output dir, read rds
# _test
dp_res <- lapply(list.files("../Results/dp_models",
                             full.names = T), 
                  FUN = readRDS)
# bind together
dp_res <- do.call(rbind.fill, dp_res)

# Analyse
dp_figs <- dlply(dp_res,
                 .(dp_min),
                 analyses_wrap_)

# plot of lags
p_dpmin_lags_sqrt <- bind_rows(fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                                 mutate(dp_min = "DP min: 3"),
                               dp_figs$`4`$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                                 mutate(dp_min = "DP min: 4"),
                               dp_figs$`5`$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                                 mutate(dp_min = "DP min: 5")
                               ) %>% 
  mutate(dp_min = factor(dp_min,
                         levels = c("DP min: 3", "DP min: 4", "DP min: 5"))) %>%
  mutate(sqrt_w_sum = sqrt(abs(w_sum)),
         sqrt_w_sum = case_when(env_var == "LUC" ~ -1 * sqrt_w_sum,
                                T ~ sqrt_w_sum)) %>%
  ggplot() +
  geom_segment(aes(x = -0.5, xend = 49.5,
                   y = 0, yend = 0),
               lty = 1, colour = "grey") +
  geom_histogram(aes(x = lag, y = sqrt_w_sum, fill = env_var), colour = NA,
                 stat = "identity", position = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("firebrick3", "dodgerblue3"),
                    breaks = c("CC", "LUC"),
                    labels = c("CC", "LUC")) +
  labs(x = "Lag (years)",
       y = expression(Sigma*"(Akaike weights)"),
       fill = "") +
  scale_y_continuous(breaks = c(-1, -sqrt(0.5), -sqrt(0.1), 0,
                                sqrt(0.1), sqrt(0.5), 1),
                     labels = c(1.0, 0.5, 0.1, 0.0, 0.1, 0.5, 1.0),
                     limits = c(-1, 1)) +
  facet_grid(dp_min~class, labeller = labeller(class = c("bird" = "Birds",
                                                         "mammal" = "Mammals")),
             scales = "free_x") +
  theme_bw() +
  basic_thm +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1,
        legend.position = c(0.9,0.05))

p_dpmin_lags_sqrt

# %%
dpmin_avg_coefs_df <- bind_rows(fig_ls$coef_df %>% 
                                  filter(model == "m1b") %>%
                                  mutate(dp_min = 3),
                                dp_figs$`4`$coef_df,
                                dp_figs$`5`$coef_df) %>%
  filter(type == "year" & delta < 6 & model == "m1b") %>%
  group_by(class, model, ecol_sub, dp_min, coef_nm) %>%
  mutate(aic_w = exp(-1/2*delta)/sum(exp(-1/2*delta))) %>%
  ungroup() %>%
  group_by(class, model, ecol_sub, dp_min, coef_nm) %>%
  dplyr::summarise(mumin_out = list(MuMIn::par.avg(coef_val, coef_se, aic_w))
  ) %>%
  mutate(mumin_out = dplyr:::map(mumin_out,~ as.data.frame.list(.x))) %>% 
  tidyr::unnest_legacy() %>%
  filter(!(coef_nm %in% c("Managed2", "Utilised2"))) %>%
  mutate(coef_nm_ = case_when(coef_nm == "(Intercept)" ~ "Int.",
                              coef_nm == "cc_s" ~ "CC",
                              coef_nm == "luc_s" ~ "LUC",
                              coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                              coef_nm == "bm_s" ~ "BM",
                              coef_nm == "paYes" ~ "PA",
                              coef_nm == "Managed1" ~ "Man",
                              coef_nm == "Utilised1" ~ "Use")) %>%
  mutate(coef_nm_ = factor(coef_nm_,
                           levels = c("Int.", "CC", "LUC", "CC:LUC",
                                      "BM", "PA", "Man", "Use")))

# plot of model averaged coefficients
p_dpmin_coefs_av <- ggplot(dpmin_avg_coefs_df) +
  geom_hline(yintercept = 0, lty = 1, colour = "darkgrey", size = 0.75) +
  geom_hline(aes(yintercept = -3.50389 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "yellow2") +
  geom_hline(aes(yintercept = -6.696701 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "orange") +
  geom_hline(aes(yintercept = -14.86601 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "red") +
  geom_errorbar(aes(x = coef_nm_, 
                    ymin = (10^Lower.CI - 1)*100,
                    ymax = (10^Upper.CI - 1)*100 ),
                width = 0,
                size = 1,
                position = position_dodge(0.5),
                show.legend = F) +
  
  geom_point(aes(x = coef_nm_, y = (10^Coefficient - 1)*100),
             size = 2.25,
             position = position_dodge(0.5),
             show.legend = F) +
  
  labs(x = element_blank(),
       y = "Annual pop. change (%)") +
  
  scale_y_continuous(limits = c(-16,16)) +
  facet_grid(dp_min~class,
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "dp_min" = c("3" = "DP min: 3",
                                              "4" = "DP min: 4",
                                              "5" = "DP min: 5"))) +
  coord_flip() +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )

p_dpmin_coefs_av

plot_grid(p_dpmin_lags_sqrt,
          p_dpmin_coefs_av,
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18)

ggsave("../Results/Figs/sm_fig_4.pdf",
       device = "pdf", dpi = 300,
       width = 15, height = 10)

# %%%%%%%%%%%%%%%%%%%%




##%%%%%
## Ecol sub analysis
##%%%%%

# Load data
# lapply output dir, read rds
# _test
ecol_res <- lapply(list.files("../Results/ecol_models",
                            full.names = T), 
                 FUN = readRDS)
# bind together
ecol_res <- do.call(rbind.fill, ecol_res)


# Focus on 1b/2b models
ecol_res <-subset(ecol_res, model %in% c("m1b", "m2b"))
# Analyses wrap
ecol_figs <- analyses_wrap_(ecol_res)


# sqrt lags v main..
p_ecol_lags_sqrt <- bind_rows(fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                                mutate(ecol_sub = "All"),
                              ecol_figs$plt_ls$m1b.carn$w_yr$w_bar_plt$data %>%
                                mutate(ecol_sub = "Carnivore"),
                              ecol_figs$plt_ls$m1b.herb$w_yr$w_bar_plt$data %>%
                                mutate(ecol_sub = "Herbivore"),
                              ecol_figs$plt_ls$m1b.tmpr$w_yr$w_bar_plt$data %>%
                                mutate(ecol_sub = "Temperate"),
                              ecol_figs$plt_ls$m1b.trop$w_yr$w_bar_plt$data %>%
                                mutate(ecol_sub = "Tropical")
                              ) %>% 
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("All", 
                                      "Herbivore", "Carnivore",
                                      "Temperate", "Tropical"))) %>%
  mutate(sqrt_w_sum = sqrt(abs(w_sum)),
         sqrt_w_sum = case_when(env_var == "LUC" ~ -1 * sqrt_w_sum,
                                T ~ sqrt_w_sum)) %>%
  ggplot() +
  geom_segment(aes(x = -0.5, xend = 49.5,
                   y = 0, yend = 0),
               lty = 1, colour = "grey") +
  geom_histogram(aes(x = lag, y = sqrt_w_sum, fill = env_var), colour = NA,
                 stat = "identity", position = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("firebrick3", "dodgerblue3"),
                    breaks = c("CC", "LUC"),
                    labels = c("CC", "LUC")) +
  labs(x = "Lag (years)",
       y = expression(Sigma*"(Akaike weights)"),
       fill = "") +
  scale_y_continuous(breaks = c(-1, -sqrt(0.5), -sqrt(0.1), 0,
                                sqrt(0.1), sqrt(0.5), 1),
                     labels = c(1.0, 0.5, 0.1, 0.0, 0.1, 0.5, 1.0),
                     limits = c(-1, 1)) +

  facet_grid(ecol_sub~class, labeller = labeller(class = c("bird" = "Birds",
                                                           "mammal" = "Mammals")),
             scales = "free_x") +
  theme_bw() +
  basic_thm +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1,
        legend.position = c(0.9,0.05))

p_ecol_lags_sqrt

# model av coefs,
avg_coefs_df <- bind_rows(subset(fig_ls$coef_df, type == "year" & model == "m1b"),
                          subset(ecol_figs$coef_df, type == "year" & model == "m1b" & 
                                   ecol_sub %in% c("carn", "herb", "tmpr", "trop"))) %>%
  filter(delta < 6) %>%
  group_by(class, model, ecol_sub, coef_nm) %>%
  mutate(aic_w = exp(-1/2*delta)/sum(exp(-1/2*delta))) %>%
  ungroup() %>%
  group_by(class, model, ecol_sub, coef_nm) %>%
  dplyr::summarise(mumin_out = list(MuMIn::par.avg(coef_val, coef_se, aic_w))
  ) %>%
  mutate(mumin_out = dplyr:::map(mumin_out,~ as.data.frame.list(.x))) %>% 
  tidyr::unnest_legacy() %>%
  filter(!(coef_nm %in% c("Managed2", "Utilised2"))) %>%
  mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                             coef_nm == "cc_s" ~ "CC",
                             coef_nm == "luc_s" ~ "LUC",
                             coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                             coef_nm == "bm_s" ~ "BM",
                             coef_nm == "paYes" ~ "PA",
                             coef_nm == "Managed1" ~ "Man",
                             coef_nm == "Utilised1" ~ "Use")) %>%
  mutate(coef_nm = factor(coef_nm,
                          levels = c("Int.", "CC", "LUC", "CC:LUC",
                                     "BM", "PA", "Man", "Use"))) %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("main", 
                                      "bm1", "bm2", "bm3",
                                      "herb", "carn", "tmpr", "trop")))

p_ecol_coefs_av <- ggplot(avg_coefs_df %>% 
                            filter(ecol_sub %in% c("main", 
                                                   "carn", "herb", "tmpr", "trop"))) +
  geom_hline(yintercept = 0, lty = 1, colour = "darkgrey", size = 0.75) +
  geom_hline(aes(yintercept = -3.50389 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "yellow2") +
  geom_hline(aes(yintercept = -6.696701 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "orange") +
  geom_hline(aes(yintercept = -14.86601 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "red") +
  geom_errorbar(aes(x = coef_nm, 
                    ymin = (10^Lower.CI - 1)*100,
                    ymax = (10^Upper.CI - 1)*100 ),
                width = 0,
                size = 1,
                position = position_dodge(0.5),
                show.legend = F) +
  
  geom_point(aes(x = coef_nm, y = (10^Coefficient - 1)*100),
             size = 2.25,
             position = position_dodge(0.5),
             show.legend = F) +
  
  labs(x = element_blank(),
       y = "Annual pop. change (%)"
  ) +
  
  scale_y_continuous(limits = c(-19.1,19.1),
                     breaks = c(-15,0,15)) +
  facet_grid(ecol_sub~class,
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "ecol_sub" = c("main" = "All",
                                                "carn" = "Carnivore",
                                                "herb" = "Herbivore",
                                                "tmpr" = "Temperate",
                                                "trop" = "Tropical"))) +
  coord_flip() +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )

p_ecol_coefs_av


# Note: luc col is specified based on the luc type used in the "best" lag-based 
# model for each ecological subset/class
p_ecol_coefs_row <- bind_rows(fig_ls$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "bird" &
                                         luc_col == "luh2_rng"
                                ),
                              
                              fig_ls$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "mammal" &
                                         luc_col == "luh2_norng"
                                ),
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "bird" &
                                         luc_col == "luh2_rng" &
                                         ecol_sub == "carn"),
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "mammal" &
                                         luc_col == "luh2_rng"&
                                         ecol_sub == "carn"),
                              
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "bird" &
                                         luc_col == "luh2_norng" &
                                         ecol_sub == "herb"),
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "mammal" &
                                         luc_col == "luh2_rng" &
                                         ecol_sub == "herb"),
                              
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "bird" &
                                         luc_col == "luh2_rng" &
                                         ecol_sub == "tmpr"),
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "mammal" &
                                         luc_col == "luh2_rng" &
                                         ecol_sub == "tmpr"),
                              
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "bird" &
                                         luc_col == "luh2_rng" &
                                         ecol_sub == "trop"),
                              ecol_figs$coef_df %>%
                                filter(model == "m1b" &
                                         type == "year" &
                                         class == "mammal" &
                                         luc_col == "luh2_norng" &
                                         ecol_sub == "trop")) %>%
  filter(cc_lag == luc_lag) %>%
  mutate(lag = cc_lag) %>%
  mutate(coef_nm = case_when(coef_nm == "(Intercept)" ~ "Int.",
                             coef_nm == "cc_s" ~ "CC",
                             coef_nm == "luc_s" ~ "LUC",
                             coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                             coef_nm == "bm_s" ~ "BM",
                             coef_nm == "paYes" ~ "PA",
                             coef_nm == "Managed1" ~ "Man",
                             coef_nm == "Utilised1" ~ "Use")) %>%
  mutate(coef_nm = factor(coef_nm,
                          levels = c("Int.", "CC", "LUC", "CC:LUC",
                                     "BM", "PA", "Man", "Use"))) %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("main", 
                                      "herb", "carn", "tmpr", "trop"))) %>%
  plot_coef_row() +
  facet_grid(ecol_sub~class,
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "ecol_sub" = c("main" = "All",
                                                "carn" = "Carnivore",
                                                "herb" = "Herbivore",
                                                "tmpr" = "Temperate",
                                                "trop" = "Tropical"
                                 )))


plot_grid(p_ecol_lags_sqrt + theme(legend.position = c(0.875,0.04)),
          p_ecol_coefs_av,
          p_ecol_coefs_row,
          ncol = 3,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b.", "c."),
          label_size = 18)

ggsave("../Results/Figs/sm_fig_5.pdf",
       device = "pdf", dpi = 300,
       width = 21, height = 13)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prediction surface and interaction type plot
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env_pred_df_av <- bind_rows(env_response_pred_av_wrap(subset(fig_ls$coef_df, 
                                                             model == "m1b" & 
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(fig_ls$plt_ls$m1b.main$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "All"),
                            
                            env_response_pred_av_wrap(subset(fig_ls$coef_df, 
                                                             model == "m1b" & 
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(fig_ls$plt_ls$m1b.main$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "All") ,
                            
                            # herb
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "herb" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.herb$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Herbivore"),
                            
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "herb" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.herb$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Herbivore") ,
                            # carn
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "carn" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.carn$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Carnivore"),
                            
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "carn" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.carn$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Carnivore") ,
                            # tmpr
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "tmpr" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.tmpr$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Temperate"),
                            
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "tmpr" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.tmpr$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Temperate") ,
                            # trop
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "trop" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.trop$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Tropical"),
                            
                            env_response_pred_av_wrap(subset(ecol_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "trop" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(ecol_figs$plt_ls$m1b.trop$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Tropical")
)   %>%
  dplyr::mutate(cc_eff_sign = sign(av_pred_cc),
                luc_eff_sign = sign(av_pred_luc),
                add_eff_sign = sign(av_pred_nointer),
                inter_eff_sign = sign(av_pred),
                piggot_class = case_when(cc_eff_sign == luc_eff_sign & 
                                           inter_eff_sign != cc_eff_sign ~ "Synergy",
                                         cc_eff_sign == luc_eff_sign & 
                                           inter_eff_sign == cc_eff_sign & 
                                           cc_eff_sign == 1 &
                                           av_pred > av_pred_nointer ~ "Synergy",
                                         cc_eff_sign == luc_eff_sign & 
                                           inter_eff_sign == cc_eff_sign & 
                                           cc_eff_sign == -1 &
                                           av_pred < av_pred_nointer ~ "Synergy",
                                         cc_eff_sign == luc_eff_sign & 
                                           inter_eff_sign == cc_eff_sign &
                                           cc_eff_sign == 1 &
                                           av_pred < av_pred_nointer ~ "Antagonism",
                                         cc_eff_sign == luc_eff_sign & 
                                           inter_eff_sign == cc_eff_sign &
                                           cc_eff_sign == -1 &
                                           av_pred > av_pred_nointer ~ "Antagonism",
                                         cc_eff_sign != luc_eff_sign &
                                           cc_eff_sign == 1 &
                                           av_pred < av_pred_cc &
                                           av_pred > av_pred_luc ~ "Antagonism",
                                         cc_eff_sign != luc_eff_sign &
                                           cc_eff_sign == -1 &
                                           av_pred > av_pred_cc &
                                           av_pred < av_pred_luc ~ "Antagonism",
                                         cc_eff_sign != luc_eff_sign &
                                           inter_eff_sign == cc_eff_sign &
                                           cc_eff_sign == 1 &
                                           av_pred > av_pred_cc ~ "Synergy",
                                         cc_eff_sign != luc_eff_sign &
                                           inter_eff_sign == cc_eff_sign &
                                           cc_eff_sign == -1 &
                                           av_pred < av_pred_cc ~ "Synergy",
                                         
                                         cc_eff_sign != luc_eff_sign &
                                           inter_eff_sign == luc_eff_sign &
                                           luc_eff_sign == 1 &
                                           av_pred > av_pred_luc ~ "Synergy",
                                         cc_eff_sign != luc_eff_sign &
                                           inter_eff_sign == luc_eff_sign &
                                           luc_eff_sign == -1 &
                                           av_pred < av_pred_luc ~ "Synergy")
  )

p_ecol_pred_av <- env_pred_df_av %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("All", "Herbivore", "Carnivore", "Temperate", "Tropical"))) %>%
  ggplot() +
  geom_tile(aes(x = cc_s, y = luc_s, 
                colour = (10^av_pred_int - 1)*100, fill = (10^av_pred_int - 1)*100)) +
  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, colour = "darkgrey") +
  geom_contour(aes(x = cc_s, y = luc_s, z = (10^av_pred_int - 1)*100), 
               colour = "yellow2",
               breaks = c(-3.50389)) +
  geom_contour(aes(x = cc_s, y = luc_s, z = (10^av_pred_int - 1)*100), 
               colour = "orange",
               breaks = c(-6.696701)) +
  geom_contour(aes(x = cc_s, y = luc_s, z = (10^av_pred_int - 1)*100), 
               colour = "red",
               breaks = c(-14.86601)) +
  scale_color_gradient2(low = "firebrick3",
                        mid = "gainsboro",
                        high = "dodgerblue3",
                        midpoint = 0,
                        # limits = coef_lim,
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "left",
                                               title.vjust = 1,
                                               label.hjust = 0.5,
                                               label.vjust = 0.5)) +
  scale_fill_gradient2(low = "firebrick3",
                       mid = "gainsboro",
                       high = "dodgerblue3",
                       midpoint = 0,
                       # limits = coef_lim,
                       guide = guide_colorbar(direction = "horizontal",
                                              title.position = "left",
                                              title.vjust = 1,
                                              label.hjust = 0.5,
                                              label.vjust = 0.5)) +
  facet_grid(ecol_sub~class) +
  labs(x = "CC",
       y = "LUC",
       colour = "Annual pop.\nchange (%)",
       fill = "Annual pop.\nchange (%)") +
  coord_equal(expand = F) +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))

p_ecol_pred_av

p_ecol_pred_av_interaction <- env_pred_df_av %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("All", "Small", "Medium", "Large",
                                      "Herbivore", "Carnivore",
                                      "Temperate", "Tropical"))) %>%
  ggplot() +
  geom_tile(aes(x = cc_s, y = luc_s,
                colour = piggot_class, fill = piggot_class)) +
  geom_contour(aes(x = cc_s, y = luc_s, z = piggot_class),
               colour = "black",
               breaks = c(0.5),
               # alpha = 0.5,
               size = 1,
               data = env_pred_df_av %>%
                 mutate(ecol_sub = factor(ecol_sub,
                                          levels = c("All", "Small", "Medium", "Large",
                                                     "Herbivore", "Carnivore",
                                                     "Temperate", "Tropical"))) %>%
                 filter(ecol_sub %in% c("All", "Herbivore", "Carnivore",
                                        "Temperate", "Tropical")) %>%
                 mutate(piggot_class = case_when(piggot_class == "Synergy" ~ 1,
                                                 piggot_class == "Antagonism" ~ 0))) +
  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, colour = "darkgrey") +
  facet_grid(ecol_sub~class) +
  labs(x = "CC",
       y = "LUC",
       colour = "Interaction type",
       fill = "Interaction type") +
  coord_equal(expand = F) +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))

p_ecol_pred_av_interaction

plot_grid(p_ecol_pred_av  +
            geom_contour(aes(x = cc_s, y = luc_s, z = piggot_class),
                         colour = "black",
                         breaks = c(0.5),
                         alpha = 0.5,
                         size = 0.5,
                         data = env_pred_df_av %>%
                           mutate(ecol_sub = factor(ecol_sub,
                                                    levels = c("All", "Small", "Medium", "Large",
                                                               "Herbivore", "Carnivore",
                                                               "Temperate", "Tropical"))) %>%
                           filter(ecol_sub %in% c("All", "Herbivore", "Carnivore",
                                                  "Temperate", "Tropical")) %>%
                           mutate(piggot_class = case_when(piggot_class == "Synergy" ~ 1,
                                                           piggot_class == "Antagonism" ~ 0))),
          p_ecol_pred_av_interaction,
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18
)
ggsave("../Results/Figs/sm_fig_7.pdf",
       width = 14, height = 13)


##%%%%%
## Env datasets
##%%%%%

# Load new results
# lapply output dir, read rds
# _test
env_res <- lapply(list.files("../Results/env_models",
                              full.names = T), 
                   FUN = readRDS)
# bind together
env_res <- do.call(rbind.fill, env_res)

# Analyses wrap
env_res$luc_col_ <- env_res$luc_col
env_res$luc_col_[grepl("luh2", env_res$luc_col)] <- "luh2"
env_res$luc_col_[grepl("hyde", env_res$luc_col)] <- "hyde"

env_figs <- dlply(env_res,
                  .(cc_col, luc_col_),
                  analyses_wrap_)


# SM Fig 8

p_env_lags_sqrt <- bind_rows(fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                               mutate(env_dat = "IPSL, LUH2"),
                             env_figs$ipsl.hyde$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                               mutate(env_dat = "IPSL, HYDE"),
                             env_figs$cru4.04.luh2$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                               mutate(env_dat = "CRU, LUH2"),
                             env_figs$cru4.04.hyde$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                               mutate(env_dat = "CRU, HYDE")
) %>% 
  mutate(env_dat = factor(env_dat,
                          levels = c("IPSL, LUH2", "IPSL, HYDE", 
                                     "CRU, LUH2", "CRU, HYDE"))) %>%
  mutate(sqrt_w_sum = sqrt(abs(w_sum)),
         sqrt_w_sum = case_when(env_var == "LUC" ~ -1 * sqrt_w_sum,
                                T ~ sqrt_w_sum)) %>%
  ggplot() +
  geom_segment(aes(x = -0.5, xend = 49.5,
                   y = 0, yend = 0),
               lty = 1, colour = "grey") +
  geom_histogram(aes(x = lag, y = sqrt_w_sum, fill = env_var), colour = NA,
                 stat = "identity", position = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("firebrick3", "dodgerblue3"),
                    breaks = c("CC", "LUC"),
                    labels = c("CC", "LUC")) +
  labs(x = "Lag (years)",
       y = expression(Sigma*"(Akaike weights)"),
       fill = "") +
  scale_y_continuous(breaks = c(-1, -sqrt(0.5), -sqrt(0.1), 0,
                                sqrt(0.1), sqrt(0.5), 1),
                     labels = c(1.0, 0.5, 0.1, 0.0, 0.1, 0.5, 1.0),
                     limits = c(-1, 1)) +
  
  facet_grid(env_dat~class, labeller = labeller(class = c("bird" = "Birds",
                                                          "mammal" = "Mammals")),
             scales = "free_x") +
  theme_bw() +
  basic_thm +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1,
        legend.position = c(0.9,0.05))

p_env_lags_sqrt

# %%
env_avg_coefs_df <- bind_rows(fig_ls$coef_df %>%
                                filter(model == "m1b") %>%
                                mutate(env_dat = "IPSL, LUH2"),
                              
                              env_figs$ipsl.hyde$coef_df %>%
                                mutate(env_dat = "IPSL, HYDE"),
                              
                              env_figs$cru4.04.luh2$coef_df %>%
                                mutate(env_dat = "CRU, LUH2"),
                              
                              env_figs$cru4.04.hyde$coef_df %>%
                                mutate(env_dat = "CRU, HYDE")
) %>%
  filter(type == "year" & delta < 6 & model == "m1b") %>%
  group_by(class, model, ecol_sub, env_dat, coef_nm) %>%
  mutate(aic_w = exp(-1/2*delta)/sum(exp(-1/2*delta))) %>%
  ungroup() %>%
  group_by(class, model, ecol_sub, env_dat, coef_nm) %>%
  dplyr::summarise(mumin_out = list(MuMIn::par.avg(coef_val, coef_se, aic_w))
  ) %>%
  mutate(mumin_out = dplyr:::map(mumin_out,~ as.data.frame.list(.x))) %>% 
  tidyr::unnest_legacy() %>%
  filter(!(coef_nm %in% c("Managed2", "Utilised2"))) %>%
  mutate(coef_nm_ = case_when(coef_nm == "(Intercept)" ~ "Int.",
                              coef_nm == "cc_s" ~ "CC",
                              coef_nm == "luc_s" ~ "LUC",
                              coef_nm == "luc_s:cc_s" ~ "CC:LUC",
                              coef_nm == "bm_s" ~ "BM",
                              coef_nm == "paYes" ~ "PA",
                              coef_nm == "Managed1" ~ "Man",
                              coef_nm == "Utilised1" ~ "Use")) %>%
  mutate(coef_nm_ = factor(coef_nm_,
                           levels = c("Int.", "CC", "LUC", "CC:LUC",
                                      "BM", "PA", "Man", "Use"))) %>%
  mutate(env_dat = factor(env_dat,
                          levels = c("IPSL, LUH2", "IPSL, HYDE", 
                                     "CRU, LUH2", "CRU, HYDE")))

View(env_avg_coefs_df)


p_env_coefs_av <- ggplot(env_avg_coefs_df) +
  geom_hline(yintercept = 0, lty = 1, colour = "darkgrey", size = 0.75) +
  geom_hline(aes(yintercept = -3.50389 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "yellow2") +
  geom_hline(aes(yintercept = -6.696701 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "orange") +
  geom_hline(aes(yintercept = -14.86601 ), 
             alpha = 0.5,
             size = 0.65,
             colour = "red") +
  geom_errorbar(aes(x = coef_nm_, 
                    ymin = (10^Lower.CI - 1)*100,
                    ymax = (10^Upper.CI - 1)*100 ),
                width = 0,
                size = 1,
                position = position_dodge(0.5),
                show.legend = F) +
  
  geom_point(aes(x = coef_nm_, y = (10^Coefficient - 1)*100),
             size = 2.25,
             position = position_dodge(0.5),
             show.legend = F) +
  
  labs(x = element_blank(),
       y = "Annual pop. change (%)"
  ) +
  
  scale_y_continuous(limits = c(-16,16)) +
  facet_grid(env_dat~class,
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "dp_min" = c("3" = "DP min: 3",
                                              "4" = "DP min: 4",
                                              "5" = "DP min: 5"),
                                 "model" = c("m1a" = "Base",
                                             "m1b" = "+MU",
                                             "m2a" = "+R",
                                             "m2b" = "+MUR"))) +
  coord_flip() +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )

p_env_coefs_av

plot_grid(p_env_lags_sqrt,
          p_env_coefs_av,
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18)

ggsave("../Results/Figs/env_split_lags_coefs_av.pdf",
       device = "pdf", dpi = 300,
       width = 15, height = 11.5)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ESA comparisons
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%
# ESA data extraction
# %%%%%


dat <- readRDS("../Data/anon_dat.rds")
bird_df1 <- subset(dat$bird_df, DP_to_2014 >= 3 & lambda_rsq > 0)
mam_df1  <- subset(dat$mam_df, DP_to_2014 >= 3 & lambda_rsq > 0)

bird_df2 <- subset(bird_df1, !(spp_idx %in% c("Gyps_bengalensis", "Podiceps_nigricollis")))
mam_df2 <- subset(mam_df1, !(ID %in% c(23514)))

loc_df <- bind_rows(mam_df2, bird_df2) %>%
  dplyr::select(loc_idx) %>%
  unique()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTE: due to anonymised data, below extraction doesn't work - lat/long missing
# Instead, we provide pre-extracted data..
# library(raster)         #~ Spatial data
# library(ncdf4)
# esa_loc_df <- unique(rbind(orig_dat$mam_df[,c("loc_id", "Longitude", "Latitude")],
#                            orig_dat$bird_df[,c("loc_id", "Longitude", "Latitude")]))
# esa_points <- SpatialPoints(esa_loc_df[, c("Longitude", "Latitude")])
# 
# esa_f_ls <- list.files("../Data/esa_lc",
#                        pattern = "ESACCI",
#                        full.names = T)
# 
# 
# m <- c(-Inf, 9, 0,  
#        10, 40, 1,  
#        41, Inf, 0)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# 
# r <- brick("~/Documents/LUH2/states.nc", var = "c3ann")
# r <- r[[1]]
# # empty raster 
# 
# 
# esa_pop_ls <- vector("list", length(esa_f_ls))
# 
# cntr <- 1
# for (f in esa_f_ls){
#   yr <- strsplit(gsub("-v2.0.7cds.nc", "", f), "-")[[1]][8]
#   esa_tmp <- raster(f)
#   # 1 - reclassify to agriculture/not
#   esa_tmp <- reclassify(esa_tmp, rclmat, right = NA)
#   # 2 - resample to 0.25 deg
#   esa_aggr <- resample(esa_tmp, r)
#   
#   # 3 - save
#   writeRaster(esa_aggr,
#               paste0("../Data/esa_agric/", yr, "_raster.tif"))
#   # 4 - Extract pop locs... - store to df
#   esa_pop <- data.frame(extract(esa_aggr, points))
#   nrow(esa_)
#   esa_pop$loc_idx <- loc_df$loc_idx
#   esa_pop$Year <- as.numeric(yr)
#   
#   esa_pop_ls[[cntr]] <- esa_pop
#   
#   cntr <- cntr + 1
# }

# # Load esa data files
# esa_f_ls_ <- list.files("~/Documents/QMEE_CDT/lag_manuscript/Revisions/Data/esa_agric",
#                         full.names = T)
# 
# 
# esa_f_ls_ <- list.files("../Data/esa_agric",
#                         full.names = T)
# # drop aux... files
# esa_f_ls_ <- esa_f_ls_[!grepl("aux", esa_f_ls_)]
# 
# esa_pop_ls <- vector("list", length(esa_f_ls_))
# 
# for (i in 1:length(esa_f_ls_)){
#   esa_aggr <- raster(esa_f_ls_[i])
#   # 4 - Extract pop locs... - store to df
#   esa_pop <- data.frame(extract(esa_aggr, esa_points))
#   
#   yr <- strsplit(gsub("_raster.tif", "", esa_f_ls_[i]), "esa_agric/")[[1]][2]
#   esa_pop$loc_id <- esa_loc_df$loc_id
#   esa_pop$Year <- as.numeric(yr)
#   
#   esa_pop_ls[[i]] <- esa_pop
#   
# }
# 
# rm(esa_aggr, esa_pop, yr)
# 
# esa_pop_df <- bind_rows(esa_pop_ls)
# rm(esa_pop_ls)
# 
# colnames(esa_pop_df) <- c("esa", "loc_id", "Year")
# 
# esa_pop_df <- left_join(esa_pop_df, 
#                         loc_df_key[c("loc_id", "loc_idx")] %>% unique())
# saveRDS(esa_pop_df[c("Year", "loc_idx", "esa")],
#         "../Data/esa_data.rds")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Get luh2 data from env
luh2_pop_df <- subset(dat$env_df, loc_idx %in% loc_df$loc_idx &
                        Year >= 1992)

luh2_crop <- readRDS("../Data/env_hist.rds")
luh2_crop <- subset(env_hist, Year >= 1992)
luh2_crop$crop <- base::rowSums(luh2_crop[,c("c3ann", "c3nfx", "c3per",
                                             "c4ann", "c4per")])

luh2_hyde <- left_join(luh2_pop_df,
                       unique(luh2_crop[,c("Year", "loc_idx", "crop")]))

# Calculate some luh2 and esa specs...
# Plots/compare
esa_change_df <- esa_pop_df %>% 
  filter(Year <= 2014) %>%
  group_by(loc_idx) %>%
  arrange(Year, .by_group = T) %>%
  summarise(esa_av_chng = mean(diff(esa)))

luh2_change_df <- luh2_hyde %>% 
  group_by(loc_idx) %>%
  arrange(Year, .by_group = T) %>%
  summarise(luh2_rng_av_chng = mean(diff(luh2_rng)),
            luh2_norng_av_chng = mean(diff(luh2_norng)),
            luh2_crop_av_chng = mean(diff(crop)),
            hyde3.2_crng_av_chng = mean(diff(hyde3.2_crng)),
            hyde3.2_grz_av_chng = mean(diff(hyde3.2_grz)))

luc_change_df <- left_join(luh2_change_df, 
                           esa_change_df)


# hyde lu types - 
# hyde..._crng = cropland + conv_rangeland + pasture
# hyde..._grz = cropland + grazing
# Grazing is a broader set that includes conv_rangeland, pasture and natural rangeland
# https://essd.copernicus.org/articles/9/927/2017/essd-9-927-2017.pdf


# ESA doesn't have a "pasture"/"grazing land" type
# 10-40 is "Agriculture"

# 10, 11, 12    Rainfed cropland
# 20            Irrigated cropland
# 30            Mosaic cropland (>50%) / natural vegetation (tree, shrub,
#               herbaceous cover) (<50%)
# 40            Mosaic natural vegetation (tree, shrub, herbaceous cover)
#               (>50%) / cropland (< 50%)



# Make some summary figs...

# Simple stats - lm, corr...
cor(luc_change_df$luh2_rng_av_chng,
    luc_change_df$esa_av_chng)

cor(luc_change_df$luh2_norng_av_chng,
    luc_change_df$esa_av_chng)

cor(luc_change_df$luh2_crop_av_chng,
    luc_change_df$esa_av_chng)

summary(lm(esa_av_chng ~ luh2_rng_av_chng, 
           data = luc_change_df))

summary(lm(esa_av_chng ~ luh2_norng_av_chng, 
           data = luc_change_df))

summary(lm(esa_av_chng ~ luh2_crop_av_chng, 
           data = luc_change_df))

data.frame(luh2_type = c("Cropland + Pasture + Rangeland",
                         "Cropland + Pasture",
                         "Cropland only"),
           slope_ = c(0.0349, 0.0373, 0.0408),
           rsq_ = c(0.0015, 0.0023, 0.0012)) %>%
  mutate(luh2_type = factor(luh2_type,
                            levels = c("Cropland + Pasture + Rangeland",
                                       "Cropland + Pasture",
                                       "Cropland only")))

bind_rows(luc_change_df[,c("loc_idx", "esa_av_chng", "luh2_rng_av_chng")] %>%
            rename(luh2_av_chng = luh2_rng_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture + Rangeland"),
          luc_change_df[,c("loc_idx", "esa_av_chng", "luh2_norng_av_chng")] %>%
            rename(luh2_av_chng = luh2_norng_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture"),
          luc_change_df[,c("loc_idx", "esa_av_chng", "luh2_crop_av_chng")] %>%
            rename(luh2_av_chng = luh2_crop_av_chng) %>%
            mutate(luh2_type = "Cropland only")
) %>%
  mutate(luh2_type = factor(luh2_type,
                            levels = c("Cropland + Pasture + Rangeland",
                                       "Cropland + Pasture",
                                       "Cropland only"))) %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, colour = "darkgrey") +
  geom_point(aes(x = luh2_av_chng,
                 y = esa_av_chng)) +
  geom_smooth(aes(x = luh2_av_chng,
                  y = esa_av_chng),
              method = "lm") +
  
  geom_text(aes(x = -0.015, y = 0.015, label = paste0("Slope: ", round(slope_, 2) )),
            hjust = 0,
            vjust = 1,
            size = 5,
            data = data.frame(luh2_type = c("Cropland + Pasture + Rangeland",
                                            "Cropland + Pasture",
                                            "Cropland only"),
                              slope_ = c(0.0349, 0.0373, 0.0408),
                              rsq_ = c(0.0015, 0.0023, 0.0012)) %>%
              mutate(luh2_type = factor(luh2_type,
                                        levels = c("Cropland + Pasture + Rangeland",
                                                   "Cropland + Pasture",
                                                   "Cropland only")))
  ) +
  
  labs(x = "LUH2 land-use change",
       y = "ESA crop cover change") +
  facet_grid(~luh2_type) +
  theme_bw() +
  basic_thm +
  scale_y_continuous(limits = c(-0.016,0.016),
                     breaks = c(-0.01,0,0.01)) +
  scale_x_continuous(limits = c(-0.016,0.016),
                     breaks = c(-0.01,0,0.01)) +
  coord_equal()

ggsave("../Results/Figs/luh2_v_esa_rates.pdf",
       width = 13, height = 6)



# Add slope(sig) to plots

# luh2_range v hyde grazing ?
bind_rows(luc_change_df[,c("loc_idx", "hyde3.2_grz_av_chng", "luh2_rng_av_chng")] %>%
            rename(luh2_av_chng = luh2_rng_av_chng,
                   hyde_av_chng = hyde3.2_grz_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture +\nRangeland",
                   hyde_type = "Cropland + Pasture +\nAll rangeland"),
          luc_change_df[,c("loc_idx", "hyde3.2_crng_av_chng", "luh2_rng_av_chng")] %>%
            rename(luh2_av_chng = luh2_rng_av_chng,
                   hyde_av_chng = hyde3.2_crng_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture +\nRangeland",
                   hyde_type = "Cropland + Pasture +\nConverted rangeland"),
          
          luc_change_df[,c("loc_idx", "hyde3.2_grz_av_chng", "luh2_norng_av_chng")] %>%
            rename(luh2_av_chng = luh2_norng_av_chng,
                   hyde_av_chng = hyde3.2_grz_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture",
                   hyde_type = "Cropland + Pasture +\nAll rangeland"),
          luc_change_df[,c("loc_idx", "hyde3.2_crng_av_chng", "luh2_norng_av_chng")] %>%
            rename(luh2_av_chng = luh2_norng_av_chng,
                   hyde_av_chng = hyde3.2_crng_av_chng) %>%
            mutate(luh2_type = "Cropland + Pasture",
                   hyde_type = "Cropland + Pasture +\nConverted rangeland")
) %>%
  mutate(luh2_type = factor(luh2_type,
                            levels = c("Cropland + Pasture +\nRangeland",
                                       "Cropland + Pasture")),
         hyde_type = factor(hyde_type,
                            levels = c("Cropland + Pasture +\nAll rangeland",
                                       "Cropland + Pasture +\nConverted rangeland"))) %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "darkgrey") +
  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
  geom_vline(xintercept = 0, lty = 2, colour = "darkgrey") +
  geom_point(aes(x = luh2_av_chng,
                 y = hyde_av_chng)) +
  geom_smooth(aes(x = luh2_av_chng,
                  y = hyde_av_chng),
              method = "lm") +
  
  geom_text(aes(x = -0.015, y = 0.015, label = paste0("Slope: ", round(slope_, 2), sig_ )),
            hjust = 0,
            vjust = 1,
            size = 5,
            data = data.frame(luh2_type = c("Cropland + Pasture +\nRangeland",
                                            "Cropland + Pasture",
                                            "Cropland + Pasture +\nRangeland",
                                            "Cropland + Pasture"),
                              hyde_type = c("Cropland + Pasture +\nAll rangeland",
                                            "Cropland + Pasture +\nAll rangeland",
                                            "Cropland + Pasture +\nConverted rangeland",
                                            "Cropland + Pasture +\nConverted rangeland"),
                              slope_ = c(0.812, 0.618, 0.773, 0.708),
                              sig_ = c("***", "***", "***", "***"),
                              rsq_ = c(0.787, 0.504, 0.709, 0.659)) %>%
              mutate(luh2_type = factor(luh2_type,
                                        levels = c("Cropland + Pasture +\nRangeland",
                                                   "Cropland + Pasture")),
                     hyde_type = factor(hyde_type,
                                        levels = c("Cropland + Pasture +\nAll rangeland",
                                                   "Cropland + Pasture +\nConverted rangeland")))
  ) +
  
  labs(x = "LUH2 land-use change",
       y = "HYDE land-use change") +
  facet_grid(hyde_type~luh2_type) +
  theme_bw() +
  basic_thm +
  # scale_y_continuous(breaks = c(-0.01,0,0.01)) +
  scale_y_continuous(limits = c(-0.016,0.016),
                     breaks = c(-0.01,0,0.01)) +
  scale_x_continuous(limits = c(-0.016,0.016),
                     breaks = c(-0.01,0,0.01)) +
  coord_equal()

ggsave("../Results/Figs/luh2_v_hyde_rates.pdf",
       width = 8, height = 8)


# Simple stats
cor(luc_change_df$luh2_rng_av_chng,
    luc_change_df$hyde3.2_grz_av_chng)

cor(luc_change_df$luh2_norng_av_chng,
    luc_change_df$hyde3.2_grz_av_chng)

cor(luc_change_df$luh2_rng_av_chng,
    luc_change_df$hyde3.2_crng_av_chng)

cor(luc_change_df$luh2_norng_av_chng,
    luc_change_df$hyde3.2_crng_av_chng)

summary(lm(hyde3.2_grz_av_chng ~ luh2_rng_av_chng, 
           data = luc_change_df))

summary(lm(hyde3.2_grz_av_chng ~ luh2_norng_av_chng, 
           data = luc_change_df))

summary(lm(hyde3.2_crng_av_chng ~ luh2_rng_av_chng, 
           data = luc_change_df))

summary(lm(hyde3.2_crng_av_chng ~ luh2_norng_av_chng, 
           data = luc_change_df))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lambda graphic 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lpd <- read.csv("../Data/lpd.csv",
                stringsAsFactors = F, na.strings = c("NULL", ""))


lemur_df <- subset(lpd, ID == 5635)
lemur_df <- data.frame("Year" = 1950:2014,
                       "abund" = unname(t(lemur_df[1,c(which(colnames(lemur_df) == "X1950"):
                                                         which(colnames(lemur_df) == "X2014"))])))
lemur_df <- na.omit(lemur_df)
lemur_df$Species = "Lemur"

hippo_df <- subset(lpd, ID == 9931)
hippo_df <- data.frame("Year" = 1950:2014,
                       "abund" = unname(t(hippo_df[1,c(which(colnames(hippo_df) == "X1950"):
                                                         which(colnames(hippo_df) == "X2014"))])))
hippo_df <- na.omit(hippo_df)
hippo_df$Species = "Hippo"
# lambda: 0.05197079

lemur_pred = data.frame("Year" = seq(min(lemur_df$Year), max(lemur_df$Year)))
smth_par <- max(round(nrow(lemur_df)/2),3)
s <- mgcv::`s`
lemur_gam <- mgcv::gam(log10(abund)~s(Year, k = smth_par), fx = T, data = lemur_df)
# Predict from pop start to end_yr

lemur_pred$abund_int <- predict(lemur_gam, lemur_pred)
lemur_lambda <- mean(diff(lemur_pred$abund), na.rm = T)

lemur_pred$Species <- "Lemur"

hippo_pred <- data.frame("Year" = min(hippo_df$Year):max(hippo_df$Year))
hippo_pred <- merge(hippo_pred, hippo_df,
                    by = "Year", all.x = T)
library(zoo)
hippo_pred$abund_int <- na.approx(log10(hippo_pred$abund))
hippo_lambda <- mean(diff(hippo_pred$abund_int), na.rm = T)

hippo_pred$Species <- "Hippo"

ggplot(bind_rows(lemur_df, hippo_df)) +
  
  geom_point(aes(x = Year, y = abund_int),
             colour = "firebrick3",
             alpha = 0.7,
             data = bind_rows(lemur_pred, hippo_pred)) +
  geom_line(aes(x = Year, y = abund_int),
            colour = "firebrick3",
            size = 1,
            alpha = 0.7,
            data = bind_rows(lemur_pred, hippo_pred)) +
  geom_point(aes(x = Year, y = log10(abund))) +
  geom_text(aes(x = x, y = y, label = lab),
            parse = T,
            size = 6,
            hjust = 0,
            vjust = 0,
            data = data.frame(Species = c("Hippo", "Lemur"),
                              x = c(1975, 1975),
                              y = c(3, 3),
                              lab = c("over(,lambda)~atop(,'= 0.05')", "over(,lambda)~atop(,'= 0.01')"))) +
  labs(y = expression("log"[10]~"Abundance")) +
  facet_grid(~Species,
             labeller = labeller("Species" = c("Hippo" = "DP < 6",
                                               "Lemur" = "DP >= 6"))) +
  theme_bw() +
  basic_thm

ggsave("../Results/Figs/lambda_calc_graphic.pdf",
       width = 10, height = 5)



