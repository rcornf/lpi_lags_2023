###%%%%%%%%%%
### Rscript to run analysis of main models + bm models
### Data summary
### Optimal lags
### Sensitivity analysis
### Coefficients (model-averaged)
### Prediction surface and interaction type
### Projections for pops in model-fitting
### Model checks - standard checks, spatio-, phyologenetic- and temporal-autocorrelation 
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
library(dplyr)          #~ Data wrangling
library(influence.ME)   #~ Sensitivity tests
library(raster)         #~ Spatial data wrangling
library(terra)
library(fasterize)
library(sf)
library(geosphere)
library(ape)            #~ Phylogenetic data wrangling
library(phytools)
library(rotl)
library(ggtree)
library(ggtreeExtra)


###%%%%%%%%%%
# Main Code
###%%%%%%%%%%

source("00_analysis_functions.R")

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Fig 1.
# Data summary
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

world_map <- map_data("world")
world_map_sub <- subset(world_map, region != "Antarctica")

dat <- readRDS("../Data/anon_dat.rds")

bird_df1 <- subset(dat$bird_df, DP_to_2014 >= 3 & lambda_rsq > 0)
mam_df1  <- subset(dat$mam_df, DP_to_2014 >= 3 & lambda_rsq > 0)

bird_df2 <- subset(bird_df1, !(spp_idx %in% c("Gyps_bengalensis", "Podiceps_nigricollis")))
mam_df2 <- subset(mam_df1, !(ID %in% c(23514)))

# Note: dat_summary_plttr will not run due to lack of lat/long in anaonymised data
dat_summ_plt_ls <- dat_summary_plttr(rbind.fill(mam_df2,
                                                bird_df2))

dat_summ_plt_ls$summ_plt

ggsave("main_fig_1.pdf",
       plot = dat_summ_plt_ls$summ_plt,
       path = "../Results/Figs/",
       device = pdf, dpi = 300,
       width = 19, height = 12)


###~~~~~
# Data summary tables
###~~~~~
mam_df2 %>%
    summarise(n_pop = n(),
              n_spp = length(unique(spp_idx)),
              n_loc = length(unique(loc_idx)))
# n_pop n_spp n_loc
#   921   287   440

bird_df2 %>%
    summarise(n_pop = n(),
              n_spp = length(unique(spp_idx)),
              n_loc = length(unique(loc_idx)))
# n_pop n_spp n_loc
#   830   425   238



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main Fig 2
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load results
###
# models for all populations
lag_res <- lapply(list.files("../Results/main_models",
                             full.names = T), 
                  FUN = readRDS)
lag_res <- do.call(rbind.fill, lag_res)
fig_ls <- analyses_wrap_(lag_res)

# %%%
# Also get results based on bm splits
bm_res <- lapply(list.files("../Results/bm_models",
                            full.names = T), 
                 FUN = readRDS)
bm_res <- do.call(rbind.fill, bm_res)
bm_figs <- analyses_wrap_(bm_res)


# Main fig 2 
# For Fig 2a: Distr of AICc support over lags,
p_bm_lags_sqrt <- bind_rows(fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
                              mutate(ecol_sub = "All"),
                            bm_figs$plt_ls$m1b.bm1$w_yr$w_bar_plt$data %>%
                              mutate(ecol_sub = "Small"),
                            bm_figs$plt_ls$m1b.bm2$w_yr$w_bar_plt$data %>%
                              mutate(ecol_sub = "Medium"),
                            bm_figs$plt_ls$m1b.bm3$w_yr$w_bar_plt$data %>%
                              mutate(ecol_sub = "Large")
                            ) %>% 
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("All", 
                                      "Small", "Medium",
                                      "Large"))) %>%
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

p_bm_lags_sqrt


# For Fig 2b: model av coefs
avg_coefs_df <- bind_rows(subset(fig_ls$coef_df, type == "year" & model == "m1b"),
                          subset(bm_figs$coef_df, type == "year" & model == "m1b")
                          ) %>%
  filter(delta < 6) %>%
  group_by(class, model, ecol_sub) %>%
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

# Check plot limits
avg_coefs_df %>% 
  filter(ecol_sub %in% c("main", 
                         "bm1", "bm2", "bm3")) %>%
  summarise(min = min((10^Lower.CI - 1)*100),
            max = max((10^Upper.CI - 1)*100))


p_bm_coefs_av <- ggplot(avg_coefs_df %>% 
                          filter(ecol_sub %in% c("main", 
                                                 "bm1", "bm2", "bm3"))) +
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
  
  scale_y_continuous(limits = c(-16,16)) +
  facet_grid(ecol_sub~class,
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "ecol_sub" = c("main" = "All",
                                                "bm1" = "Small",
                                                "bm2" = "Medium",
                                                "bm3" = "Large"))) +
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

p_bm_coefs_av

# For Fig 2c: coefficient values over lags
# Note: luc col is specified based on the luc type used in the "best" lag-based 
# model for each ecological subset/class
p_bm_coefs_row <- bind_rows(fig_ls$coef_df %>%
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
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "bird" &
                                       luc_col == "luh2_rng" &
                                       ecol_sub == "bm1"),
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "mammal" &
                                       luc_col == "luh2_norng"&
                                       ecol_sub == "bm1"),
                            
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "bird" &
                                       luc_col == "luh2_norng" &
                                       ecol_sub == "bm2"),
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "mammal" &
                                       luc_col == "luh2_norng" &
                                       ecol_sub == "bm2"),
                            
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "bird" &
                                       luc_col == "luh2_norng" &
                                       ecol_sub == "bm3"),
                            bm_figs$coef_df %>%
                              filter(model == "m1b" &
                                       type == "year" &
                                       class == "mammal" &
                                       luc_col == "luh2_norng" &
                                       ecol_sub == "bm3")) %>%
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
                                      "bm1", "bm2", "bm3"))) %>%
  plot_coef_row()


plot_grid(p_bm_lags_sqrt,
          p_bm_coefs_av,
          p_bm_coefs_row,
          ncol = 3,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b.", "c."),
          label_size = 18)

ggsave("../Results/Figs/main_fig_2.pdf",
       device = "pdf", dpi = 300,
       width = 21, height = 11.5)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# SM fig 17, 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fit null models for comparison
m_fit_null(mam_df2)
# model R2m       R2c      AICc
# mnull   0 0.3824281 -3105.338
m_fit_null(bird_df2)
# model R2m       R2c      AICc
# mnull   0 0.4074627 -2484.389


aic_sub <- subset(fig_ls$coef_df[,c("model", "AICc", "cc_col", "luc_col", "cc_lag", "luc_lag", "class", "type")],
                  type == "year" & model %in% c("m1a", "m1b", "m2a", "m2b")) %>%
  unique() %>%
  mutate(model_ = gsub("m", "", model))

ggplot(aic_sub) +
  facet_wrap(~class, scales = "free_y",
             labeller = labeller(class = c("bird" = "Birds",
                                           "mammal" = "Mammals"))) +
  # Null models
  geom_point(aes(x = 0.5, y = AICc),
             size = 3,
             alpha = 0.5,
             data = data.frame(class = c("bird", "mammal"),
                               AICc  = c(-2484.398, -3105.338))
  ) +
  
  geom_text(aes(x = 0.5, y = AICc, label = lab),
            colour = "grey30",
            size = 5,
            vjust = -1,
            data = data.frame(class = c("bird", "mammal"),
                              AICc  = c(-2484.398, -3105.338),
                              lab = "Null")
  ) +
  #                        model        R2m       R2c      AICc
  # m_fit_no_env(bird_df2) # mnull 0.00000000 0.4074627 -2484.389
  # m_fit_no_env(mam_df2)  # mnull 0.00000000 0.3824281 -3105.338
  
  geom_flat_violin(aes(x = as.factor(model_), y = AICc, group = interaction(class, model_),
                       fill = model_),
                   width = 0.75,
                   position = position_nudge(x = 0.175, y = 0),
                   alpha = 0.25,
                   show.legend = F) +
  
  geom_boxplot(aes(x = as.factor(model_), y = AICc,
                   fill = model_),
               width = 0.2,
               outlier.size = 1,
               alpha = 0.25,
               show.legend = F) +
  
  geom_point(aes(x = as.factor(model_), y = AICc),
             shape = 18,
             size = 3.75,
             alpha = 0.85,
             colour = "grey40",
             data = subset(aic_sub, cc_lag == 0 & luc_lag == 0 ) %>%
               group_by(model_, class) %>%
               summarise(AICc = mean(AICc))) +
  
  geom_segment(aes(x = x, xend = x,
                   y = low, yend = high),
               colour = "grey30",
               size = 0.9,
               arrow = arrow(length = unit(0.07, "inch")),
               data = aic_sub %>%
                 group_by(class) %>%
                 summarise(cent = (min(AICc) + max(AICc))/2 ,
                           rng = (diff(range(AICc)))/3,
                           low = cent + rng,
                           high = cent - rng) %>%
                 mutate(x = 0.5)
  ) +
  
  geom_text(aes(x = x, y = y, label = label, angle = 90),
            size = 5,
            fontface = "bold",
            colour = "grey30",
            data = aic_sub %>%
              group_by(class) %>%
              summarise(y = (min(AICc) + max(AICc))/2) %>%
              mutate(x = 0.3,
                     label = "Better")
  ) +
  
  scale_fill_manual(values = c("1a" = NA,
                               "1b" = "purple",
                               "2a" = NA,
                               "2b" = NA)) +
  scale_y_continuous(breaks = c(-2475, -2500, -2525, -2550, 
                                -3100,  -3150, -3200),
                     expand = expansion(add = c(5,10))) +
  scale_x_discrete(expand = expansion(add = c(0.9,0.6)),
                   ## Edits to model naming
                   breaks = c("1a","1b","2a","2b"), 
                   labels = c("Base", "+MU", "+R", "+MUR")
                   ##
  ) +
  labs(x = "Model structure") + 
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")) 

ggsave(filename = "../Results/Figs/sm_fig_17.pdf",
       device = "pdf", dpi = 300,
       width = 10, height = 5)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SM fig 3 - +MU v +MUR 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_mspec_lags_sqrt <- bind_rows(
  fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt$data %>%
    mutate(mspec = "+MU"),
  fig_ls$plt_ls$m2b.main$w_yr$w_bar_plt$data %>%
    mutate(mspec = "+MUR")
) %>% 
  mutate(mspec = factor(mspec,
                        levels = c("+MU", "+MUR"))) %>%
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
  facet_grid(mspec~class, labeller = labeller(class = c("bird" = "Birds",
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
        legend.position = c(0.9,0.075))

p_mspec_lags_sqrt


mspec_avg_coefs_df <- subset(fig_ls$coef_df, type == "year") %>%
  filter(delta < 6) %>%
  group_by(class, model, ecol_sub) %>%
  mutate(aic_w = exp(-1/2*delta)/sum(exp(-1/2*delta))) %>%
  ungroup() %>%
  group_by(class, model, ecol_sub, coef_nm) %>%
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
                              coef_nm == "Utilised1" ~ "Use",
                              coef_nm == "LPI_RealmIndo-Pacific" ~ "Indo-Pacific",
                              coef_nm == "LPI_RealmNearctic" ~ "Nearctic",    
                              coef_nm == "LPI_RealmNeotropical" ~ "Neotropical",
                              coef_nm == "LPI_RealmPalearctic" ~ "Palearctic")) %>%
  mutate(coef_nm_ = factor(coef_nm_,
                           levels = c("Int.", "CC", "LUC", "CC:LUC",
                                      "BM", "PA", "Man", "Use",
                                      "Indo-Pacific", "Nearctic",
                                      "Neotropical", "Palearctic"))) %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("main", 
                                      "bm1", "bm2", "bm3",
                                      "herb", "carn", "tmpr", "trop")))

p_mspec_coefs_av <- ggplot(mspec_avg_coefs_df %>% filter(model %in% c("m1b", "m2b"))) +
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
  facet_grid(model~class,
             scales = "free_y",
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "ecol_sub" = c("main" = "All",
                                                "bm1" = "Small",
                                                "bm2" = "Medium",
                                                "bm3" = "Large"),
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

p_mspec_coefs_av


plot_grid(p_mspec_lags_sqrt,
          p_mspec_coefs_av,
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18)

ggsave("../Results/Figs/sm_fig_3.pdf",
       device = "pdf", dpi = 300,
       width = 16, height = 8)


###~~~~~
# Summary of top models
###~~~~~

# Birds, m1b (+MU)
subset(fig_ls$plt_ls$m1b.main$metr_df, 
       class == "bird" & delta < 6 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank class cc_lag luc_lag  luc_col      AICc    delta       R2m       R2c
#    1  bird     14       3 luh2_rng -2553.318 0.000000 0.1834291 0.3174285
#    2  bird     14       2 luh2_rng -2549.859 3.458849 0.1730877 0.3186594
#    3  bird     19       3 luh2_rng -2549.005 4.313248 0.1743792 0.3329938
#    4  bird     14       1 luh2_rng -2548.948 4.369771 0.1702377 0.3210351
#    5  bird     41       2 luh2_rng -2548.871 4.446890 0.1756743 0.3901477
#    6  bird     41       3 luh2_rng -2548.502 4.815785 0.1763175 0.3874378

subset(fig_ls$plt_ls$m1b.main$metr_df, 
       class == "bird" & cc_lag == 0 & luc_lag == 0 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank class cc_lag luc_lag    luc_col      AICc    delta        R2m       R2c
# 1990  bird      0       0   luh2_rng -2506.574 46.74404 0.07417949 0.4032305
# 2167  bird      0       0 luh2_norng -2505.810 47.50817 0.08081041 0.4055484


# Mammals, m1b (+MU)
subset(fig_ls$plt_ls$m1b.main$metr_df, 
       class == "mammal" & delta < 6 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank  class cc_lag luc_lag    luc_col      AICc    delta       R2m       R2c
#    1 mammal     45       9 luh2_norng -3199.891 0.000000 0.1369968 0.4214862
#    2 mammal     45       7   luh2_rng -3194.344 5.547881 0.1384221 0.4421451

subset(fig_ls$plt_ls$m1b.main$metr_df, 
       class == "mammal" & cc_lag == 0 & luc_lag == 0 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank  class cc_lag luc_lag    luc_col      AICc    delta        R2m       R2c
# 2356 mammal      0       0 luh2_norng -3130.171 69.72039 0.06529260 0.3946373
# 2774 mammal      0       0   luh2_rng -3127.805 72.08673 0.06387996 0.4028054
    

# Birds, m2b (+MUR)
subset(fig_ls$plt_ls$m2b.main$metr_df, 
       class == "bird" & delta < 6 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank class cc_lag luc_lag  luc_col      AICc    delta       R2m       R2c
#    1  bird     14       3 luh2_rng -2552.943 0.000000 0.1971594 0.3214737
#    2  bird     14       2 luh2_rng -2549.569 3.373506 0.1878301 0.3202958
#    3  bird     14       1 luh2_rng -2548.425 4.517548 0.1850055 0.3220259
#    4  bird     19       3 luh2_rng -2547.837 5.105251 0.1844564 0.3362990

subset(fig_ls$plt_ls$m2b.main$metr_df, 
       class == "bird" & cc_lag == 0 & luc_lag == 0 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank class cc_lag luc_lag    luc_col      AICc    delta        R2m       R2c
# 1830  bird      0       0   luh2_rng -2505.543 47.40004 0.09524658 0.3876529
# 2100  bird      0       0 luh2_norng -2503.914 49.02854 0.09686121 0.3942151

# Mammals, m2b (+MUR)
subset(fig_ls$plt_ls$m2b.main$metr_df, 
       class == "mammal" & delta < 6 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank  class cc_lag luc_lag    luc_col      AICc      delta       R2m       R2c
#    1 mammal     45       7   luh2_rng -3202.087 0.00000000 0.1512826 0.4335078
#    2 mammal     45       9 luh2_norng -3202.024 0.06269648 0.1462607 0.4227943
#    3 mammal     45       8   luh2_rng -3201.915 0.17250576 0.1516658 0.4349049

subset(fig_ls$plt_ls$m2b.main$metr_df, 
       class == "mammal" & cc_lag == 0 & luc_lag == 0 & 
           type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                              "AICc", "delta", "R2m", "R2c")]
# Rank  class cc_lag luc_lag    luc_col      AICc    delta        R2m       R2c
# 2421 mammal      0       0 luh2_norng -3127.650 74.43722 0.07113252 0.3922652
# 2618 mammal      0       0   luh2_rng -3126.418 75.66868 0.07095993 0.3996402


# Can use similar code to generate stats re top models for each body-mass subset
# e.g.
subset(bm_figs$plt_ls$m1b.bm1$metr_df, 
       class == "mammal" & delta < 6 & 
         type == "year")[,c("Rank", "class", "cc_lag", "luc_lag", "luc_col",
                            "AICc", "delta", "R2m", "R2c")]


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sensitivity analysis - top models
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

infl_assess_specs <- unique(fig_ls$coef_df[fig_ls$coef_df$delta == 0,
                                           c("model","class","cc_lag","luc_lag","type",
                                             "cc_col","luc_col","dp_min","R2m","R2c",
                                             "AICc", "delta","id","dur","off","ecol_sub")]) %>%
    filter(type == "year" & model %in% c("m1a", "m1b", "m2a", "m2b"))

infl_assess_specs$id <- 1:nrow(infl_assess_specs)

mam_infl <- dlply(subset(infl_assess_specs, class == "mammal"),
                  .(model),
                  function(x, y, z){
                      infl_assess_wrap(y, x, z)
                  },  
                  y = mam_df2,
                  z = dat$env_df)

lapply(mam_infl,
       function(x){
           print(x$mod_spec$model)
           print(as.character(x$mod_spec$class))
           print(as.character(x$mod_spec$type))
           
           print(x$pop_cook[x$pop_cook$V1 > 0.5,])
           print(x$spp_cook[x$spp_cook$V1 > 0.5,])
           print(x$loc_cook[x$loc_cook$V1 > 0.5,])
       })

# No influential pops/spp/locs

bird_infl <- dlply(subset(infl_assess_specs, class == "bird"),
                   .(model),
                   function(x, y, z){
                       infl_assess_wrap(y, x, z)
                   },  
                   y = bird_df2,
                   z = dat$env_df)

lapply(bird_infl,
       function(x){
           print(x$mod_spec$model)
           print(as.character(x$mod_spec$class))
           print(as.character(x$mod_spec$type))
           
           print(x$pop_cook[x$pop_cook$V1 > 0.5,])
           print(x$spp_cook[x$spp_cook$V1 > 0.5,])
           print(x$loc_cook[x$loc_cook$V1 > 0.5,])
       })
# [1] "m1a"
# [1] "bird"
# [1] "year"
# V1        drop_id idx
# 0.6486315     108   1
# 0.5403033     157   2
# V1                     drop_id idx
# 0.6486315 Anser_brachyrhynchus   1
# V1            drop_id idx
# 0.6099781     loc_657   1
# 0.5403033     loc_638   2
# [1] "m1b"
# [1] "bird"
# [1] "year"
# 
# [1] "m2a"
# [1] "bird"
# [1] "year"
# V1       drop_id idx
# 0.6402953 loc_657   1
# [1] "m2b"
# [1] "bird"
# [1] "year"
# V1       drop_id idx
# 0.5581252 loc_657   1


bird_df2[108,] # Anser_brachyrhynchus
bird_df2[157,] # Oxyura_leucocephala
subset(bird_df2, loc_idx == "loc_657") # Lake Ashenghe, Ethiopia, East Africa
subset(bird_df2, loc_idx == "loc_638") # Oxyura_leucocephala


# Cooks D > 0.5 figs...
cooks_plttr(c(mam_infl, bird_infl), "year", "pop_cook") + theme_bw() + basic_thm + ggtitle("pop_cook")
cooks_plttr(c(mam_infl, bird_infl), "year", "spp_cook") + theme_bw() + basic_thm + ggtitle("spp_cook")
cooks_plttr(c(mam_infl, bird_infl), "year", "loc_cook") + theme_bw() + basic_thm + ggtitle("loc_cook")


## Set up loo-cv specs
# Prepare spec df for running loo/loobsr - only look at top models and 1b
cv_specs <- unique(fig_ls$coef_df[fig_ls$coef_df$delta < 6,
                                           c("model","class","cc_lag","luc_lag","type",
                                             "cc_col","luc_col","dp_min","R2m","R2c",
                                             "AICc", "delta","id","dur","off","ecol_sub")])

cv_specs <- subset(cv_specs, model %in% c("m1b", "m2b") & type == "year")

cv_specs$id <- 1:nrow(cv_specs)
saveRDS(cv_specs,
        "../Data/cv_specs.rds")  


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main Fig 3
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###~~~~~~~~~
# 3a Prediction surface
###~~~~~~~~~

# Get model-averaged predictions
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
                            
                            # small
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm1" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm1$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Small"),
                            
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm1" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm1$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Small") ,
                            # medium
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm2" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm2$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Medium"),
                            
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm2" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm2$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Medium") ,
                            # large
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm3" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "mammal"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm3$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "mammal")) %>%
                              dplyr::mutate(class = "Mammals",
                                            ecol_sub = "Large"),
                            
                            env_response_pred_av_wrap(subset(bm_figs$coef_df, 
                                                             model == "m1b" & 
                                                               ecol_sub == "bm3" &
                                                               delta < 6 &
                                                               type == "year" &
                                                               class == "bird"),
                                                      
                                                      subset(bm_figs$plt_ls$m1b.bm3$metr_df,
                                                             delta < 6 &
                                                               type == "year" &
                                                               class == "bird")) %>%
                              dplyr::mutate(class = "Birds",
                                            ecol_sub = "Large") 
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

p_bm_pred_av <- env_pred_df_av %>%
  filter(ecol_sub %in% c("All", "Small", "Medium", "Large")) %>%
  mutate(ecol_sub = factor(ecol_sub,
                           levels = c("All", "Small", "Medium", "Large"))) %>%
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
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "left",
                                               title.vjust = 1,
                                               label.hjust = 0.5,
                                               label.vjust = 0.5)) +
  scale_fill_gradient2(low = "firebrick3",
                       mid = "gainsboro",
                       high = "dodgerblue3",
                       midpoint = 0,
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

p_bm_pred_av



###~~~~~~~~~
# 3b Future projections
###~~~~~~~~~

# Extract env covariates
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTE: below code can't run as Long/Lat are removed from anonymised LPI data

# # Extract env data for pop locs
# loc_df <- bind_rows(bird_df2[,c("loc_idx", "Longitude", "Latitude")], 
#                     mam_df2[,c("loc_idx", "Longitude", "Latitude")]) %>%
#   unique() %>%
#   rename(x = Longitude,
#          y = Latitude)
# 
# # luh2_extract
# luh2_hist <- luh2_extract(loc_df)
# luh2_hist_sub <- subset(luh2_hist, Year >= 1950 & Year <= 2014)
# 
# luh2_126 <- luh2_fut_extract(loc_df, "ssp1_rcp2.6.nc")
# luh2_370 <- luh2_fut_extract(loc_df, "ssp3_rcp7.0.nc")
# luh2_585 <- luh2_fut_extract(loc_df, "ssp5_rcp8.5.nc")
# 
# # ipsl_extract
# ipsl_hist <- ipsl_extract(loc_df)
# ipsl_hist_sub <- subset(ipsl_hist, Year >= 1950 & Year <= 2014)
# 
# ipsl_126 <- ipsl_fut_extract(loc_df, "IPSL_126")
# ipsl_370 <- ipsl_fut_extract(loc_df, "IPSL_370")
# ipsl_585 <- ipsl_fut_extract(loc_df, "IPSL_585")
# 
# 
# # Combine ipsl and luh2
# env_hist <- left_join(luh2_hist_sub,
#                       ipsl_hist_sub)
# 
# env_126 <- left_join(luh2_126,
#                       ipsl_126)
# env_370 <- left_join(luh2_370,
#                       ipsl_370)
# env_585 <- left_join(luh2_585,
#                       ipsl_585)
# 
# # Save... 
# saveRDS(env_hist, 
#         "../Data/env_hist.rds")
# 
# saveRDS(env_126, 
#         "../Data/env_126.rds")
# saveRDS(env_370, 
#         "../Data/env_370t.rds")
# saveRDS(env_585, 
#         "../Data/env_585.rds")
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fit models
mam_lmm_dat <- lmm_bm_fit_store(subset(fig_ls$plt_ls$m1b.main$metr_df,
                                       type == "year" &
                                         model == "m1b" &
                                         delta < 6 & 
                                         class == "mammal"),
                                mam_df2,
                                dat$env_df
                                )
bird_lmm_dat <- lmm_bm_fit_store(subset(fig_ls$plt_ls$m1b.main$metr_df,
                                       type == "year" &
                                         model == "m1b" &
                                         delta < 6 & 
                                         class == "bird"),
                                bird_df2,
                                dat$env_df)

mam_bm1_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm1$metr_df,
                                       type == "year" &
                                         model == "m1b" &
                                         delta < 6 & 
                                         class == "mammal"),
                                mam_df2,
                                dat$env_df)


bird_bm1_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm1$metr_df,
                                        type == "year" &
                                          model == "m1b" &
                                          delta < 6 & 
                                          class == "bird"),
                                 bird_df2,
                                 dat$env_df)

mam_bm2_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm2$metr_df,
                                       type == "year" &
                                         model == "m1b" &
                                         delta < 6 & 
                                         class == "mammal"),
                                mam_df2,
                                dat$env_df)


bird_bm2_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm2$metr_df,
                                        type == "year" &
                                          model == "m1b" &
                                          delta < 6 & 
                                          class == "bird"),
                                 bird_df2,
                                 dat$env_df)

mam_bm3_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm3$metr_df,
                                       type == "year" &
                                         model == "m1b" &
                                         delta < 6 & 
                                         class == "mammal"),
                                mam_df2,
                                dat$env_df)


bird_bm3_dat <- lmm_bm_fit_store(subset(bm_figs$plt_ls$m1b.bm3$metr_df,
                                        type == "year" &
                                          model == "m1b" &
                                          delta < 6 & 
                                          class == "bird"),
                                 bird_df2,
                                 dat$env_df)


# Load env data
# Bind hist to each sc - to account for lagd
pop_env_hist <- readRDS("../Data/env_hist.rds")
pop_env_126  <- readRDS("../Data/env_126.rds")
pop_env_370  <- readRDS("../Data/env_370.rds")
pop_env_585  <- readRDS("../Data/env_585.rds")

pop_env_126 <- bind_rows(pop_env_hist,
                         pop_env_126)
pop_env_370 <- bind_rows(pop_env_hist,
                         pop_env_370)
pop_env_585 <- bind_rows(pop_env_hist,
                         pop_env_585)


# Project models
# Run for birds and mammals,
mam_pop_proj <- mod_proj(mam_lmm_dat, 
                         list("sc_126" = pop_env_126,
                              "sc_370" = pop_env_370,
                              "sc_585" = pop_env_585))

bird_pop_proj <- mod_proj(bird_lmm_dat, 
                          list("sc_126" = pop_env_126,
                               "sc_370" = pop_env_370,
                               "sc_585" = pop_env_585))

# Birds - bm splits
bird_bm1_proj <- mod_proj(bird_bm1_dat, 
                          list("sc_126" = pop_env_126,
                               "sc_370" = pop_env_370,
                               "sc_585" = pop_env_585))

bird_bm2_proj <- mod_proj(bird_bm2_dat, 
                          list("sc_126" = pop_env_126,
                               "sc_370" = pop_env_370,
                               "sc_585" = pop_env_585))

bird_bm3_proj <- mod_proj(bird_bm3_dat, 
                          list("sc_126" = pop_env_126,
                               "sc_370" = pop_env_370,
                               "sc_585" = pop_env_585))

# saveRDS(list(bird_bm1_proj = bird_bm1_proj,
#              bird_bm2_proj = bird_bm2_proj,
#              bird_bm3_proj = bird_bm3_proj),
#         "../Results/bird_bm_pop_proj.rds")

# Mammals - bm splits
mam_bm1_proj <- mod_proj(mam_bm1_dat, 
                         list("sc_126" = pop_env_126,
                              "sc_370" = pop_env_370,
                              "sc_585" = pop_env_585))

mam_bm2_proj <- mod_proj(mam_bm2_dat, 
                         list("sc_126" = pop_env_126,
                              "sc_370" = pop_env_370,
                              "sc_585" = pop_env_585))

mam_bm3_proj <- mod_proj(mam_bm3_dat, 
                         list("sc_126" = pop_env_126,
                              "sc_370" = pop_env_370,
                              "sc_585" = pop_env_585))

# saveRDS(list(mam_bm1_proj = mam_bm1_proj,
#              mam_bm2_proj = mam_bm2_proj,
#              mam_bm3_proj = mam_bm3_proj),
#         "../Results/mam_bm_pop_proj.rds")


# set 2010 as baseline
w_av_sc_fut_proj_df <- bind_rows(expand.grid(class = c("bird", "mammal"),
                                                 sc = c("sc_126","sc_370","sc_585"),
                                                 ecol_sub = c("All", "bm1", 
                                                              "bm2", "bm3")) %>%
                                       data.frame() %>%
                                       mutate(TS_end = 2010,
                                              pred_incl_ranef_mult = 1,
                                              pred_excl_ranef_mult = 1),
                                     
                                     # Alls
                                     bird_pop_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = "All"),
                                     mam_pop_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = "All"),
                                     # BM1
                                     mam_bm1_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm1")),
                                     bird_bm1_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm1")),
                                     # BM2
                                     mam_bm2_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm2")),
                                     bird_bm2_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm2")),
                                     # BM3
                                     mam_bm3_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm3")),
                                     bird_bm3_proj$w_av_sc_proj_df %>%
                                       mutate(ecol_sub = c("bm3")))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make projection figures
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_bm_proj_fut_incl_ranef_2100 <- ggplot(w_av_sc_fut_proj_df %>%
                                          filter(TS_end <= 2100)) +
  geom_hline(yintercept = 1, lty = 1, colour = "darkgrey", size = 0.75) +
  geom_rect(aes(xmin = 2010, xmax = 2010+min_lag,
                ymin = 0, ymax = Inf),
            alpha = 0.1,
            data = bind_rows(subset(fig_ls$coef_df,
                                    delta == 0 &
                                      model == "m1b" &
                                      type == "year") %>%
                               dplyr::mutate(ecol_sub = "All"),
                             subset(bm_figs$coef_df,
                                    delta == 0 &
                                      type == "year")) %>%
              dplyr::select(model, ecol_sub, class, cc_lag, luc_lag, R2m, R2c, AICc,
                            cc_col, luc_col, type, delta) %>%
              unique() %>%
              dplyr::mutate(min_lag = case_when(cc_lag <= luc_lag ~ cc_lag,
                                                T ~ luc_lag)) %>%
              left_join(w_av_sc_fut_proj_df %>%
                          dplyr::group_by(ecol_sub) %>%
                          dplyr::summarise(min_Idx = min(pred_incl_ranef_Idx)))) +
  geom_rect(aes(xmin = 2010, xmax = 2010+max_lag,
                ymin = 0, ymax = Inf),
            alpha = 0.1,
            data = bind_rows(subset(fig_ls$coef_df,
                                    delta == 0 &
                                      model == "m1b" &
                                      type == "year") %>%
                               dplyr::mutate(ecol_sub = "All"),
                             subset(bm_figs$coef_df,
                                    delta == 0 &
                                      type == "year")) %>%
              select(model, ecol_sub, class, cc_lag, luc_lag, R2m, R2c, AICc,
                     cc_col, luc_col, type, delta) %>%
              unique() %>%
              dplyr::mutate(max_lag = case_when(cc_lag >= luc_lag ~ cc_lag,
                                                T ~ luc_lag)) %>%
              left_join(w_av_sc_fut_proj_df %>%
                          dplyr::group_by(ecol_sub) %>%
                          dplyr::summarise(min_Idx = min(pred_incl_ranef_Idx)))) +
  geom_line(aes(x = TS_end, y = pred_incl_ranef_Idx, colour = sc),
            size = 1.2,
            alpha = 0.9) +
  scale_color_manual(values = c(
    "sc_126" = "dodgerblue3",
    "sc_370" = "orange",
    "sc_585" = "firebrick3"),
    labels = c(
      "SSP1 RCP2.6",
      "SSP3 RCP7.0",
      "SSP5 RCP8.5")) +
  scale_x_continuous(breaks = c(2050, 2100)) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Year",
       y = "Index",
       colour = element_blank()) +
  facet_grid(ecol_sub~class,
             # scales = "free_y",
             labeller = labeller("class" = c("bird" = "Birds",
                                             "mammal" = "Mammals"),
                                 "ecol_sub" = c("All" = "All",
                                                "bm1" = "Small",
                                                "bm2" = "Medium",
                                                "bm3" = "Large"))) +
  theme_bw() +
  basic_thm +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        legend.position = c(0.175,0.7),
        legend.title = element_blank(),
        legend.margin = margin(t = 0, unit = "cm"),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))


# Combined figure 3
plot_grid(p_bm_pred_av,
          p_bm_proj_fut_incl_ranef_2100 +
            theme(legend.position = "bottom",
                  legend.direction = "vertical"),
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18
)

ggsave("../Results/Figs/bm_split_pred_av_proj_av_2100.pdf",
       width = 14, height = 11)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# And free scale for sm fig 18
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Mammals and birds in own columns for better scales...

plot_grid(ggplot(w_av_sc_fut_proj_df %>% filter(class == "bird")) +
            geom_hline(yintercept = 1, lty = 1, colour = "darkgrey", size = 0.75) +
            geom_rect(aes(xmin = 2010, xmax = 2010+min_lag,
                          ymin = 0, ymax = Inf),
                      alpha = 0.1,
                      data = lag_lim_df %>% filter(class == "bird")) +
            geom_rect(aes(xmin = 2010, xmax = 2010+max_lag,
                          ymin = 0, ymax = Inf),
                      alpha = 0.1,
                      data = lag_lim_df %>% filter(class == "bird")) +
            geom_vline(xintercept = 2100, lty = 2,
                       alpha = 0.5) +
            geom_line(aes(x = TS_end, y = pred_incl_ranef_Idx, colour = sc),
                      size = 1.2,
                      alpha = 0.9) +
            scale_color_manual(values = c("sc_126" = "dodgerblue3",
                                          "sc_370" = "orange",
                                          "sc_585" = "firebrick3"),
                               labels = c("SSP1 RCP2.6",
                                          "SSP3 RCP7.0",
                                          "SSP5 RCP8.5")) +
            scale_y_continuous(trans = "log10") +
            labs(x = "Year",
                 y = "Index",
                 colour = element_blank()) +
            facet_grid(ecol_sub~class,
                       scales = "free_y",
                       labeller = labeller("class" = c("bird" = "Birds",
                                                       "mammal" = "Mammals"),
                                           "ecol_sub" = c("All" = "All",
                                                          "bm1" = "Small",
                                                          "bm2" = "Medium",
                                                          "bm3" = "Large"))) +
            theme_bw() +
            # theme_classic() +
            basic_thm +
            theme(aspect.ratio = 1,
                  panel.grid = element_blank(),
                  legend.position = "none",
                  legend.title = element_blank(),
                  legend.margin = margin(t = 0, unit = "cm"),
                  panel.border = element_rect(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.text = element_text(colour = "black")),
          
          ggplot(w_av_sc_fut_proj_df %>% filter(class == "mammal")) +
            geom_hline(yintercept = 1, lty = 1, colour = "darkgrey", size = 0.75) +
            geom_rect(aes(xmin = 2010, xmax = 2010+min_lag,
                          ymin = 0, ymax = Inf),
                      alpha = 0.1,
                      data = lag_lim_df %>% filter(class == "mammal")) +
            geom_rect(aes(xmin = 2010, xmax = 2010+max_lag,
                          ymin = 0, ymax = Inf),
                      alpha = 0.1,
                      data = lag_lim_df %>% filter(class == "mammal")) +
            geom_vline(xintercept = 2100, lty = 2,
                       alpha = 0.5) +
            geom_line(aes(x = TS_end, y = pred_incl_ranef_Idx, colour = sc),
                      size = 1.2,
                      alpha = 0.9) +
            scale_color_manual(values = c("sc_126" = "dodgerblue3",
                                          "sc_370" = "orange",
                                          "sc_585" = "firebrick3"),
                               labels = c("SSP1 RCP2.6",
                                          "SSP3 RCP7.0",
                                          "SSP5 RCP8.5")) +
            scale_y_continuous(trans = "log10") +
            labs(x = "Year",
                 y = "",
                 colour = element_blank()) +
            facet_grid(ecol_sub~class,
                       scales = "free_y",
                       labeller = labeller("class" = c("bird" = "Birds",
                                                       "mammal" = "Mammals"),
                                           "ecol_sub" = c("All" = "All",
                                                          "bm1" = "Small",
                                                          "bm2" = "Medium",
                                                          "bm3" = "Large"))) +
            theme_bw() +
            # theme_classic() +
            basic_thm +
            theme(aspect.ratio = 1,
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  legend.direction = "vertical",
                  legend.title = element_blank(),
                  legend.margin = margin(t = 0, unit = "cm"),
                  panel.border = element_rect(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.text = element_text(colour = "black")),
          
          align = "hv",
          axis = "tblr"
)

ggsave("../Results/Figs/sm_fig_18.pdf",
       width = 7.5, height = 11)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SM Fig 6
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_bm_pred_av_interaction <- env_pred_df_av %>%
  filter(ecol_sub %in% c("All", "Small", "Medium", "Large")) %>%
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
                 filter(ecol_sub %in% c("All", "Small", "Medium", "Large")) %>%
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


plot_grid(p_bm_pred_av +
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
                           filter(ecol_sub %in% c("All", "Small", "Medium", "Large")) %>%
                           mutate(piggot_class = case_when(piggot_class == "Synergy" ~ 1,
                                                           piggot_class == "Antagonism" ~ 0))) ,
          p_bm_pred_av_interaction,
          ncol = 2,
          align = "hv",
          axis = "tblr",
          labels = c("a.", "b."),
          label_size = 18
)
ggsave("../Results/Figs/sm_fig_6.pdf",
       width = 14, height = 11)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##%%%%%
## SM Fig 16: Lag graphic 
##%%%%%

lpd <- read.csv("../Data/lpd.csv",
                stringsAsFactors = F, na.strings = c("NULL", ""))


pop_graphic_plttr(11886, lpd, dat$env_df, dat$mam_df) +
  theme_classic() +
  theme_bw() +
  basic_thm +
  theme(panel.grid = element_blank())

ggsave("../Results/Figs/sm_fig_16.pdf",
       device = "pdf", dpi = 300,
       width = 7.5, height = 8)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Model checks  
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mam_lmm_dat
bird_lmm_dat

# %%%%%%%
# General model checks
# %%%%%%%

mam_0_mod <- mam_lmm_dat$mod_dat[[which(mam_lmm_dat$mod_spec$delta == 0)]]$mod_obj
bird_0_mod <- bird_lmm_dat$mod_dat[[which(bird_lmm_dat$mod_spec$delta == 0)]]$mod_obj


ihs <- function(x){
  log(x + ( sqrt(x^2 + 1) ) )
}

mam_m_ihs <- lme4::lmer(ihs(lambda) ~ luc_s * cc_s + bm_s + pa + Managed + Utilised + 
                          (1 | spp_idx) + (1 | loc_idx) ,
                        data = mam_lmm_dat$mod_dat[[1]]$df)


bird_m_ihs <- lme4::lmer(ihs(lambda) ~ luc_s * cc_s + bm_s + pa + Managed + Utilised + 
                           (1 | spp_idx) + (1 | loc_idx) ,
                         data = bird_lmm_dat$mod_dat[[4]]$df)


mam_0_check <- check_model(mam_0_mod)
mam_ihs_check <- check_model(mam_m_ihs)

mam_0_check_plts <- plot(mam_0_check, return_list = T)
mam_ihs_check_plts <- plot(mam_ihs_check, return_list = T)

bird_0_check <- check_model(bird_0_mod)
bird_ihs_check <- check_model(bird_m_ihs)

bird_0_check_plts <- plot(bird_0_check, return_list = T)
bird_ihs_check_plts <- plot(bird_ihs_check, return_list = T)


# 1,6,2,3
plot_grid(#  Posterior Predictive Check
  mam_0_check_plts[[1]] + labs(x = expression(bar(lambda)),
                               title = element_blank(),
                               subtitle = "a.") +
    basic_thm +
    theme(legend.position = c(0.75,0.9)),
  NULL,
  mam_ihs_check_plts[[1]] + labs(x = expression(bar(lambda)),
                                 title = element_blank(),
                                 subtitle = "b.") +
    basic_thm +
    theme(legend.position = c(0.75,0.9)),
  #  Q-Q plot
  mam_0_check_plts[[6]] + labs(title = element_blank(),
                               subtitle = "   c.") +
    basic_thm,
  NULL,
  mam_ihs_check_plts[[6]] + labs(title = element_blank(),
                                 subtitle = "   d.") +
    basic_thm,
  #  Residuals v Fitted
  mam_0_check_plts[[2]] + labs(title = element_blank(),
                               subtitle = "   e.") +
    basic_thm,
  NULL,
  mam_ihs_check_plts[[2]] + labs(title = element_blank(),
                                 subtitle = "   f.") +
    basic_thm,
  
  align = "hv",
  axis = "tblr",
  labels = c("Untransformed response", "",
             "IHS transformed response",
             "","","",
             "","",""),
  label_size = 18,
  ncol = 3,
  rel_widths = c(0.475,0.05,0.475)
)

ggsave("../Results/Figs/mam_model_checks.pdf",
       width = 12, height = 12)

plot_grid(#  Posterior Predictive Check
  bird_0_check_plts[[1]] + labs(x = expression(bar(lambda)),
                                title = element_blank(),
                                subtitle = "a.") +
    basic_thm +
    theme(legend.position = c(0.275,0.9)),
  NULL,
  bird_ihs_check_plts[[1]] + labs(x = expression(bar(lambda)),
                                  title = element_blank(),
                                  subtitle = "b.") +
    basic_thm +
    theme(legend.position = c(0.275,0.9)),
  #  Q-Q plot
  bird_0_check_plts[[6]] + labs(title = element_blank(),
                                subtitle = "   c.") +
    basic_thm,
  NULL,
  bird_ihs_check_plts[[6]] + labs(title = element_blank(),
                                  subtitle = "   d.") +
    basic_thm,
  #  Residuals v Fitted
  bird_0_check_plts[[2]] + labs(title = element_blank(),
                                subtitle = "   e.") +
    basic_thm,
  NULL,
  bird_ihs_check_plts[[2]] + labs(title = element_blank(),
                                  subtitle = "   f.") +
    basic_thm,
  
  align = "hv",
  axis = "tblr",
  labels = c("Untransformed response", "",
             "IHS transformed response",
             "","","",
             "","",""),
  label_size = 18,
  ncol = 3,
  rel_widths = c(0.475,0.05,0.475)
)

ggsave("../Results/Figs/bird_model_checks.pdf",
       width = 12, height = 12)


# %%%%%%%
# Assess spatial autocorrelation
# %%%%%%%

# NOTE: won't run as anonymised data lacks coordinates.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mam_dists <- as.matrix(dist(mam_lmm_dat$mod_dat[[1]]$df[,c("Longitude", "Latitude")]))


calc_geogr_dist <- function(long_lat_df){
  n_row <- nrow(long_lat_df)
  matr_long <- data.frame("loc_id1" = rep(long_lat_df$loc_idx, each = n_row),
                          "long1" = rep(long_lat_df$Longitude, each = n_row),
                          "lat1" = rep(long_lat_df$Latitude, each = n_row),
                          "loc_id2" = rep(long_lat_df$loc_idx, n_row),
                          "long2" = rep(long_lat_df$Longitude, n_row),
                          "lat2" = rep(long_lat_df$Latitude, n_row))
  
  
  matr_long$dist_m <- geosphere::distHaversine(matr_long[,c("long1", "lat1")],
                                               matr_long[,c("long2", "lat2")])
  
  matr_ <- matrix(matr_long$dist_m, nrow = n_row, ncol = n_row)
  return(matr_)
}


# Mammals - geographic distance matrix
mam_dists_geog <- calc_geogr_dist(mam_lmm_dat$mod_dat[[1]]$df[,c("loc_idx", "Longitude", "Latitude")])
mam_w <- 1/mam_dists_geog
mam_w[mam_w>1] <- 1


# Spatial corr in residuals
Moran.I(resid(mam_0_mod), mam_w, alternative = "less")
Moran.I(resid(mam_0_mod), mam_w, alternative = "greater")
Moran.I(resid(mam_0_mod), mam_w)
diag(mam_w) <- 0


# Birds
bird_dists_geog <- calc_geogr_dist(bird_lmm_dat$mod_dat[[4]]$df[,c("loc_idx", "Longitude", "Latitude")])
bird_w <- 1/bird_dists_geog
bird_w[bird_w>1] <- 1
diag(bird_w) <- 0


Moran.I(resid(bird_0_mod), bird_w, alternative = "less")
Moran.I(resid(bird_0_mod), bird_w, alternative = "greater")
Moran.I(resid(bird_0_mod), bird_w)

# Negative spatial corr in residuals
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# %%%%%%%%%
# Assess phylogenetic autocorrelation
# %%%%%%%%%

# NOTE: code will not run properly with anonymised LPI data - Binomial is not present

get_spp_tree <- function(mod_df){
  spp_nms <- unique(mod_df$spp_idx)
  
  rotl_search <- tnrs_match_names(names = gsub("_", " ", spp_nms), 
                                  context_name = "All life")
  
  in_tree <- ott_id(rotl_search)[is_in_tree(ott_id(rotl_search))]
  spp_tree <- tol_induced_subtree(ott_ids = in_tree)
  
  spp_df <- data.frame(model_binomial = spp_nms)
  
  spp_df$ott_name <- rotl_search$unique_name
  spp_df$ott_id <- rotl_search$ott_id
  spp_df <- spp_df %>%
    mutate(tip_label = paste0(ott_name, "_ott", ott_id)) %>%
    mutate(tip_label = gsub(" ", "_", tip_label))
  
  return(list(spp_df = spp_df,
              spp_tree = spp_tree,
              missing_spp = setdiff(spp_df$tip_label, 
                                    spp_tree$tip.label))
  )
}

calc_phylosig <- function(mod_obj, spp_df, spp_tree){
  # Get error
  err_df <- data.frame(model_binomial = mod_obj@frame$spp_idx,
                       residual = residuals(mod_obj)) %>%
    left_join(spp_df)
  
  # Phylo signal
  tmp_vec = err_df$residual
  names(tmp_vec) = err_df$tip_label
  phylo_res <- phylosig(compute.brlen(spp_tree), tmp_vec, method = "lambda", test = T)
  
  return(list(error_df = err_df,
              phylo_res = phylo_res)
  )
  
}

bird_spp_tree <- get_spp_tree(bird_lmm_dat$mod_dat[[4]]$df)
mam_spp_tree <- get_spp_tree(mam_lmm_dat$mod_dat[[1]]$df)


compute.brlen(mam_spp_tree$spp_tree)

bird_phylo_1 <- calc_phylosig(bird_0_mod, 
                              bird_spp_tree$spp_df,
                              bird_spp_tree$spp_tree)

mam_phylo_1 <- calc_phylosig(mam_0_mod, 
                             mam_spp_tree$spp_df,
                             mam_spp_tree$spp_tree)
bird_phylo_1$phylo_res
# Birds, no phylo sig
mam_phylo_1$phylo_res
# Mammals, non-sig


# plot error on tree

mam_phylo_1$error_df %>% 
  filter(tip_label %in% mam_spp_tree$spp_tree$tip.label) %>%
  group_by(tip_label) %>% 
  dplyr::summarise(error = mean(residual)) %>%
  pull(error) %>% range()

bird_phylo_1$error_df %>% 
  filter(tip_label %in% bird_spp_tree$spp_tree$tip.label) %>%
  group_by(tip_label) %>% 
  dplyr::summarise(error = mean(residual)) %>%
  pull(error) %>% range()

bird_tree_plt <- ggtree(compute.brlen(bird_spp_tree$spp_tree), layout = "circular") %<+% 
  (bird_spp_tree$spp_df %>% 
     filter(tip_label %in% bird_spp_tree$spp_tree$tip.label)) +
  geom_fruit(geom=geom_tile,
             mapping = aes(y=tip_label, x=1, fill=error),
             width=0.1,offset=-0.14,
             data=bird_phylo_1$error_df %>% 
               filter(tip_label %in% bird_spp_tree$spp_tree$tip.label) %>%
               group_by(tip_label) %>% 
               dplyr::summarise(error = mean(residual))
  ) +
  scale_fill_gradient2("Error", midpoint = 0, 
                       limits = c(-0.3619195, 0.2336986),
                       high = "#D55E00", 
                       low = "#0072B2", 
                       mid = "gainsboro", na.value = "darkgrey", 
                       guide = "colourbar") 


bird_tree_plt

mam_tree_plt <- ggtree(compute.brlen(mam_spp_tree$spp_tree), layout = "circular") %<+% 
  (mam_spp_tree$spp_df %>% 
     filter(tip_label %in% mam_spp_tree$spp_tree$tip.label)) +
  geom_fruit(geom=geom_tile,
             mapping = aes(y=tip_label, x=1, fill=error),
             # pwidth=0.5,
             width=0.1,offset=-0.14,
             data=mam_phylo_1$error_df %>% 
               filter(tip_label %in% mam_spp_tree$spp_tree$tip.label) %>%
               group_by(tip_label) %>% 
               dplyr::summarise(error = mean(residual))) +
  scale_fill_gradient2("Error", midpoint = 0, 
                       limits = c(-0.3619195, 0.2336986),
                       high = "#D55E00", 
                       low = "#0072B2", 
                       mid = "gainsboro", na.value = "darkgrey", 
                       guide = "colourbar")
# offset=-0.46, 
mam_tree_plt


plot_grid(bird_tree_plt + theme(legend.position = "none"),
          mam_tree_plt,
          align = "hv",
          axis = "tblr",
          labels = c("a. Birds", "b. Mammals"),
          label_size = 18
)
ggsave("../Results/Figs/phylo_error.pdf",
       width = 15, height = 7)



# %%%%%%%%%%%
# Temporal autocorrelation checks
# %%%%%%%%%%%

mam_err_df <- mam_lmm_dat$mod_dat[[1]]$df
mam_err_df$error <- residuals(mam_0_mod)

bird_err_df <- bird_lmm_dat$mod_dat[[4]]$df
bird_err_df$error <- residuals(bird_0_mod)

# average residuals per year, lm test..

mam_av_temp_err_df <- mam_err_df[,c("TS_strt", "error")] %>%
  group_by(TS_strt) %>%
  summarise(av_error = mean(error))
bird_av_temp_err_df <- bird_err_df[,c("TS_strt", "error")] %>%
  group_by(TS_strt) %>%
  summarise(av_error = mean(error))

plot(head(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1),
     tail(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1))
summary(lm(tail(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1) ~ head(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1)))

plot(head(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1),
     tail(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1))
summary(lm(tail(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1) ~ head(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1)))


plot_grid(
  bind_rows(mam_err_df[,c("TS_strt", "error")] %>%
              mutate(class = "Mammals"),
            bird_err_df[,c("TS_strt", "error")] %>%
              mutate(class = "Birds")) %>%
    ggplot() +
    geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
    geom_point(aes(x = TS_strt, y = error),
               alpha = 0.4) +
    geom_blank(aes(x = x, y = y),
               data = data.frame(x = c(1950,2010, 1950,2010),
                                 y = c(-max(abs(range(mam_err_df$error))),
                                       max(abs(range(mam_err_df$error))),
                                       -max(abs(range(bird_err_df$error))),
                                       max(abs(range(bird_err_df$error)))),
                                 class = c("Mammals", "Mammals", "Birds", "Birds"))) +
    labs(x = "Time-series start",
         y = "Residual error") +
    facet_wrap(~class, scales = "free_y") +
    theme_bw() +
    basic_thm +
    theme(aspect.ratio = 1),
  
  bind_rows(data.frame(x = head(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1),
                       y = tail(mam_av_temp_err_df$av_error, nrow(mam_av_temp_err_df)-1),
                       class = "Mammals"),
            data.frame(x = head(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1),
                       y = tail(bird_av_temp_err_df$av_error, nrow(bird_av_temp_err_df)-1),
                       class = "Birds")) %>%
    ggplot() +
    # geom_abline(slope = 1) +
    geom_hline(yintercept = 0, colour = "darkgrey", lty = 2) +
    geom_point(aes(x = x, y = y)) +
    geom_smooth(aes(x = x, y = y), method = "lm") +
    geom_blank(aes(x = x, y = y),
               data = data.frame(x = c(-max(abs(range(mam_av_temp_err_df$av_error))),
                                       max(abs(range(mam_av_temp_err_df$av_error))),
                                       -max(abs(range(bird_av_temp_err_df$av_error))),
                                       max(abs(range(bird_av_temp_err_df$av_error)))),
                                 y = c(-max(abs(range(mam_av_temp_err_df$av_error))),
                                       max(abs(range(mam_av_temp_err_df$av_error))),
                                       -max(abs(range(bird_av_temp_err_df$av_error))),
                                       max(abs(range(bird_av_temp_err_df$av_error)))),
                                 class = c("Mammals", "Mammals", "Birds", "Birds"))) +
    facet_wrap(~class, scales = "free") +
    labs(x = expression("Average residual error"[(t)]),
         y = expression("Average residual error"[(t+1)])) +
    # coord_equal() +
    theme_bw() +
    basic_thm + 
    theme(aspect.ratio = 1),
  
  nrow = 2,
  labels = c("a.", "b."),
  align = "hv",
  axis = "tblr",
  label_size = 18
)

ggsave("../Results/Figs/temporal_error_.pdf",
       width = 10, height = 10)

