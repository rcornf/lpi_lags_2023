###%%%%%%%%%%
### Rscript to conduct analysis of initial models
### Assesses optimal lags based on AIC/Akaike weights, compares year- to 
### generation-based lags, and runs senitivity analyses of the top models
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
library(influence.ME)

# Specify some themes for ggplot
basic_thm <- theme(axis.text = element_text(size = 16),
                   axis.title = element_text(size = 20),
                   plot.title = element_text(size = 20),
                   plot.subtitle = element_text(size = 18),
                   strip.text.x = element_text(size = 18),
                   strip.text.y = element_text(size = 18),
                   strip.background = element_blank(),
                   legend.text = element_text(size = 16),
                   legend.title = element_text(size = 18),
                   legend.text.align = 0) 

clrs <- c("firebrick3", "dodgerblue3")
names(clrs) <- c("mammal", "bird")

clrs_ <- c("firebrick3", "dodgerblue3")
names(clrs_) <- c("Mammalia", "Aves")



###%%%%%%%%%%
# Functions
###%%%%%%%%%%

w_grd_plttr <- function(df, type_flag){
    
    if (type_flag == "year"){
        x_lab <- "CC lag (years)" 
        y_lab <- "LUC lag (years)" 
    }
    if (type_flag == "gen"){
        x_lab <- "CC lag (generations)" 
        y_lab <- "LUC lag (generations)" 
    }
    
    ggplot(df) + 
        geom_tile(aes(x = cc_lag, y = luc_lag, colour = w_sum, fill = w_sum)) +
        scale_colour_gradient(low = "gainsboro", high = "firebrick3") +
        scale_fill_gradient(low = "gainsboro", high = "firebrick3") +
        labs(x = x_lab,
             y = y_lab,
             fill = expression(Sigma*"(Akaike weights)"),
             colour = expression(Sigma*"(Akaike weights)")) +
        coord_equal(expand = F) +
        facet_grid(~class, 
                   labeller = labeller(class = c("bird" = "Birds", 
                                                 "mammal" = "Mammals"))) +
        theme(strip.background = element_blank()) +
        basic_thm
}


w_bar_plttr <- function(df, type_flag){
    
    if (type_flag == "year"){
        x_lab <- "Lag (years)" 
        # set diff between lags
        lag_diff <- 1
    }
    if (type_flag == "gen"){
        x_lab <- "Lag (generations)" 
        lag_diff <- 0.1
    }
    
    zero_line_df <- ddply(df,
                          .(class),
                          function(tmp){
                              data.frame(x = min(tmp$cc_lag)-lag_diff/2,
                                         xend = max(tmp$cc_lag)+lag_diff/2,
                                         y = 0, yend = 0)
                          })
    
    bar_df <- rbind.fill(ddply(df, 
                               .(class, cc_lag),
                               function(x){data.frame(w_sum = sum(x$w_sum),
                                                      env_var = "CC",
                                                      lag = x$cc_lag[1])}),
                         ddply(df, 
                               .(class, luc_lag),
                               function(x){data.frame(w_sum = -sum(x$w_sum),
                                                      env_var = "LUC",
                                                      lag = x$luc_lag[1])}))
    
    # bar_df$w_sum[bar_df$class == "bird"] <- -1*bar_df$w_sum[bar_df$class == "bird"]
    
    lim <- round(max(abs(bar_df$w_sum)), 1)
    if (lim == 0){
        lim <- round(max(abs(bar_df$w_sum)), 2)
    }
    
    if (lim>0.2){
        n_ <- 5
    }
    
    else{n_ <- 3}
    
    y_ <- max(abs(bar_df$w_sum))
    ggplot(bar_df) +
        geom_histogram(aes(x = lag, y = w_sum, fill = env_var), colour = NA,
                       stat = "identity", position = "identity", alpha = 0.9) +
        geom_segment(aes(x =x, xend = xend, 
                         y = y, yend = yend),
                     lty = 1, colour = "grey",
                     data = zero_line_df) +
        scale_fill_manual(values = c("firebrick3", "dodgerblue3"),
                          breaks = c("CC", "LUC"),
                          labels = c("CC", "LUC")) +
        labs(x = x_lab,
             y = expression(Sigma*"(Akaike weights)"),
             fill = "") +
        scale_y_continuous(breaks = seq(-1,1,length.out = 3),
                           labels = abs(seq(-1,1,length.out = 3)),
                           limits = c(-1, 1)) +
        facet_grid(class~., labeller = labeller(class = c("bird" = "Birds",
                                                          "mammal" = "Mammals")),
                   scales = "free_y") +
        basic_thm +
        theme(strip.background = element_blank())
}

akaike_w_calc_ <- function(metr_df, type_){
    w_df <- ddply(metr_df, 
                  .(class),
                  function(x){
                      x$delta <- x$AICc - min(x$AICc)
                      x$m_weights <- exp(-1/2*x$delta)/sum(exp(-1/2*x$delta))
                      x
                  })
    w_summ_df <- ddply(w_df,
                       .(class, type, cc_lag, luc_lag),
                       function(x){
                           data.frame(w_sum = sum(x$m_weights),
                                      w_mn  = mean(x$m_weights))
                       })
    
    w_grd_plt  <- w_grd_plttr(subset(w_summ_df, type == type_), type_)
    w_bar_plt  <- w_bar_plttr(subset(w_summ_df, type == type_), type_)
    
    return(list(w_df = w_df, w_summ_df = w_summ_df,
                w_bar_plt = w_bar_plt,
                w_grd_plt = w_grd_plt))
}


mod_rank_plttr <- function(df){
    # get ranks of models with delta<thr
    m_rank_delta_df <- ddply(df,
                             .(class),
                             function(x){
                                 d2 <- subset(x, delta <2)
                                 d2 <- subset(d2, Rank == max(Rank))
                                 d6 <- subset(x, delta <6)
                                 d6 <- subset(d6, Rank == max(Rank))
                                 d10 <- subset(x, delta <10)
                                 d10 <- subset(d10, Rank == max(Rank))
                                 
                                 # w95 <- subset(x, cum_w < 0.95)
                                 # w95 <- subset(w95, Rank == max(Rank))
                                 
                                 d2$delt_thr <- 2
                                 d6$delt_thr <- 6
                                 d10$delt_thr <- 10
                                 
                                 # w95$cum_w_thr <- 0.95
                                 
                                 return(rbind.fill(d2,d6,d10)) #,w95
                             })
    # identify max n models needed
    n_max <- max(m_rank_delta_df$Rank)
    # round to next 5 (or 10)
    x_max <- round((n_max+3)/5)*5
    
    if (x_max>=10){rank_breaks <- seq(10,x_max, 10)}
    if (x_max<10){rank_breaks <- c(5)}
    
    plt_df <- subset(df, Rank <=x_max)
    
    mod_rank_delta_plt <- ggplot(plt_df) +
        geom_line(aes(x = Rank, y = delta),
                  size = 1) +

        geom_segment(aes(x = Rank, xend = Rank,
                         y = 0, yend = delta),
                     lty = 2, colour = "grey",
                     size = 0.75,
                     data = m_rank_delta_df) +
        geom_segment(aes(x = 0.5, xend = Rank,
                         y = delta, yend = delta),
                     lty = 2, colour = "grey",
                     size = 0.75,
                     data = m_rank_delta_df) +

        labs(x = "Model rank",
             y = expression(Delta*" AICc")) +
        scale_y_continuous(breaks = seq(0,max(plt_df$delta), 5),
                           expand = expansion(add = c(0,1))) +
        scale_x_continuous(breaks = c(1,rank_breaks),
                           expand = expansion(add = c(0,0))) +
        facet_grid(class~., scales = "free_y",
                   labeller = labeller(class = c("bird" = "Birds",
                                                 "mammal" = "Mammals"))) +
        theme_classic() +
        basic_thm
    
    mod_rank_delta_plt
}

analyses_wrap_ <- function(df){
    
    # ignore loo/loobsr for now...
    coef_df <- subset(df, ecol_sub != "loo" & ecol_sub != "loobsr") 
    
    coef_df <- ddply(df, .(class, type, model, ecol_sub), # , dp_min
                     function(x){
                         x$delta <- x$AICc - min(x$AICc)
                         x
                     })
    
    # Plotting
    plt_ls <- dlply(coef_df, .(model, ecol_sub), #, dp_min
                    function(x){
                        # unique metric rows
                        metr_cn <- c("id","class","cc_lag","luc_lag","dur","off","type","model",
                                     "R2m","R2c","AICc","cc_col","luc_col","ecol_sub", "delta")
                        # "nolag_dAICc",
                        metr_df <- unique(x[,metr_cn])
                        
                        metr_df <- ddply(metr_df, .(class, type), # , dp_min
                                         function(x){
                                             x <- x[order(x$delta),]
                                             x$Rank <- 1:nrow(x)
                                             x
                                         })
                        # plot nolag_dAICc, per type, overall and per (cc_col, luc_col)
                        metr_df$type_class <- paste(metr_df$type, metr_df$class, sep = " ")
                        metr_df$env_comb <- paste(metr_df$cc_col, metr_df$luc_col, sep = " ")
                        
                        # Drop no-lag data before summarising
                        metr_sub_df <- subset(metr_df, !(cc_lag == 0 & luc_lag == 0))
                        
                        
                        aicc_plt <- ggplot(metr_sub_df) +
                            geom_violin(aes(x = type, y = AICc
                            )) +
                            geom_point(aes(x = type, y = AICc,
                                           group = paste(cc_col, luc_col)),
                                       position = position_dodge(0.9),
                                       data = subset(metr_df,
                                                     cc_lag == 0 & luc_lag == 0)) +
                            facet_grid(class~luc_col, scales = "free",
                                       labeller = labeller(class = c("bird" = "Birds",
                                                                     "mammal" = "Mammals"),
                                                           luc_col = c("luh2_norng" = "Excl. rangeland",
                                                                       "luh2_rng" = "Incl. rangeland"))) +
                            basic_thm + 
                            theme(strip.background = element_blank()) +
                            labs(x = "Lag type",
                                 y = "AICc",
                                 colour = "Env data",
                                 fill = "Env data")
                        
                        # calc overall model AIC - i.e. per class but across type and env comb
                        # plot w support for lags..
                        # split for yr and gen separately
                        w_yr <- akaike_w_calc_(subset(metr_df, type == "year"), "year")
                        w_gen <- akaike_w_calc_(subset(metr_df, type == "gen"), "gen")
                        
                        # also plot rank v delta - upto delta < 10
                        # yr
                        yr_rank_plt <- mod_rank_plttr(subset(metr_df, type == "year"))
                        # gen
                        gen_rank_plt <- mod_rank_plttr(subset(metr_df, type == "gen"))
                        
                        # return(list of plots)
                        return(list(aicc_plt = aicc_plt,
                                    w_yr  = w_yr,
                                    w_gen = w_gen,
                                    yr_rank_plt = yr_rank_plt,
                                    gen_rank_plt = gen_rank_plt,
                                    metr_df = metr_df))
                    })
    
    return(list(coef_df = coef_df,
                plt_ls = plt_ls))
}

# Code coming from @drob: https://gist.github.com/dgrtwo/eb7750e74997891d7c20#file-geom_flat_violin-r
"%||%" <- function(a, b) {
    if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomFlatViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            trim = trim,
            scale = scale,
            ...
        )
    )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
    ggproto("GeomFlatViolin", Geom,
            setup_data = function(data, params) {
                data$width <- data$width %||%
                    params$width %||% (resolution(data$x, FALSE) * 0.9)
                
                # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
                data %>%
                    group_by(group) %>%
                    mutate(ymin = min(y),
                           ymax = max(y),
                           xmin = x,
                           xmax = x + width / 2)
                
            },
            
            draw_group = function(data, panel_scales, coord) {
                # Find the points for the line to go all the way around
                data <- transform(data, xminv = x,
                                  xmaxv = x + violinwidth * (xmax - x))
                
                # Make sure it's sorted properly to draw the outline
                newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                                 plyr::arrange(transform(data, x = xmaxv), -y))
                
                # Close the polygon: set first and last point the same
                # Needed for coord_polar and such
                newdata <- rbind(newdata, newdata[1,])
                
                ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
            },
            
            draw_key = draw_key_polygon,
            
            default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                              alpha = NA, linetype = "solid"),
            
            required_aes = c("x", "y")
    )


geom_flat_violin_L <- function(mapping = NULL, data = NULL, stat = "ydensity",
                               position = "dodge", trim = TRUE, scale = "area",
                               show.legend = NA, inherit.aes = TRUE, ...) {
    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomFlatViolin_L,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            trim = trim,
            scale = scale,
            ...
        )
    )
}
GeomFlatViolin_L <-
    ggproto("GeomFlatViolin_L", Geom,
            setup_data = function(data, params) {
                data$width <- data$width %||%
                    params$width %||% (resolution(data$x, FALSE) * 0.9)
                
                # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
                data %>%
                    group_by(group) %>%
                    mutate(ymin = min(y),
                           ymax = max(y),
                           xmin = x,
                           xmax = x - width / 2)
                
            },
            
            draw_group = function(data, panel_scales, coord) {
                # Find the points for the line to go all the way around
                data <- transform(data, xminv = x,
                                  xmaxv = x - violinwidth * (x - xmax ))
                
                # Make sure it's sorted properly to draw the outline
                newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                                 plyr::arrange(transform(data, x = xmaxv), -y))
                
                # Close the polygon: set first and last point the same
                # Needed for coord_polar and such
                newdata <- rbind(newdata, newdata[1,])
                
                ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
            },
            
            draw_key = draw_key_polygon,
            
            default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                              alpha = NA, linetype = "solid"),
            
            required_aes = c("x", "y")
    )


## Env rate calc functions
cc_calc <- function(strt, end, loc, env_df){
    tmp_env <- subset(env_df, loc_idx == loc &
                          Year %in% (strt:end))
    strt_end <- cc1 <- cc2 <- cc3 <- NA
    if (nrow(tmp_env)>1){
        tmp_env <- tmp_env[order(tmp_env$Year),]
        strt_end <- paste(range(tmp_env$Year), collapse = "_")
        
        # calc rates
        if ("cru4.04" %in% colnames(tmp_env)){
            if (sum(!is.na(tmp_env$cru4.04))>=2){
                cc2 <- lm(cru4.04~Year, data = tmp_env)$coefficients[2]
            }
        }
        if (sum(!is.na(tmp_env$ipsl))>=2){
            cc3 <- lm(ipsl~Year, data = tmp_env)$coefficients[2]
        }
    }
    return(list(strt_end = strt_end, 
                cc1 = cc1, cc2 = cc2, cc3 = cc3))
}


luc_calc <- function(strt, end, loc, env_df){
    tmp_env <- subset(env_df, loc_idx == loc &
                          Year %in% (strt:end))
    strt_end <- luc1 <- luc2 <- luc3 <- luc4 <- luc5 <- NA
    if (nrow(tmp_env)>1){
        # order env dfs by year
        tmp_env <- tmp_env[order(tmp_env$Year),]
        strt_end <- paste(range(tmp_env$Year), collapse = "_")
        
        # calc rates
        if ("hyde3.2_crng" %in% colnames(tmp_env)){
            luc2 <- mean(diff(tmp_env$hyde3.2_crng))
            luc3 <- mean(diff(tmp_env$hyde3.2_grz))
        }
        luc4 <- mean(diff(tmp_env$luh2_norng))
        luc5 <- mean(diff(tmp_env$luh2_rng))
        
    }
    return(list(strt_end = strt_end, 
                luc1 = luc1, luc2 = luc2, luc3 = luc3, luc4 = luc4, luc5 = luc5))
}

env_rate_calc <- function(df, spec, env){
    df$cru4.04 <- NA
    df$ipsl <- NA
    
    df$hyde3.2_crng <- NA
    df$hyde3.2_grz <- NA
    df$luh2_norng <- NA
    df$luh2_rng <- NA
    
    # need to know lag type, dur, off, and lag 
    if (spec$type == "gen"){
        df$cc_gen_lag <- round(spec$cc_lag*df$GenLength)
        df$luc_gen_lag <- round(spec$luc_lag*df$GenLength)
    }
    
    # Calc cc/luc lags
    for (i in 1:nrow(df)){
        
        if (spec$dur == "ts"){
            if (spec$type == "gen"){
                cc_strt <- df$TS_strt[i] - df$cc_gen_lag[i]
                cc_end <- df$TS_end[i] - df$cc_gen_lag[i]
                luc_strt <- df$TS_strt[i] - df$luc_gen_lag[i]
                luc_end <- df$TS_end[i] - df$luc_gen_lag[i]
            }
            if (spec$type == "year"){
                cc_strt <- df$TS_strt[i] - spec$cc_lag 
                cc_end <- df$TS_end[i] - spec$cc_lag
                luc_strt <- df$TS_strt[i] - spec$luc_lag 
                luc_end <- df$TS_end[i] - spec$luc_lag
            }
        }
        
        # CC and luc rate calc
        cc_dat <- cc_calc(cc_strt, cc_end, df$loc_idx[i], env)
        luc_dat <- luc_calc(luc_strt, luc_end, df$loc_idx[i], env)
        
        df$cru3.23[i] <- cc_dat$cc1
        df$cru4.04[i] <- cc_dat$cc2
        df$ipsl[i] <- cc_dat$cc3
        
        df$hyde3.1[i] <- luc_dat$luc1
        df$hyde3.2_crng[i] <- luc_dat$luc2
        df$hyde3.2_grz[i] <- luc_dat$luc3
        df$luh2_norng[i] <- luc_dat$luc4
        df$luh2_rng[i] <- luc_dat$luc5
    }
    
    df
}

infl_assess_wrap <- function(pop_df, mod_spec, env_df){
    # prep model data
    m_df <- env_rate_calc(pop_df, mod_spec, env_df)
    
    m_df$cc_s <- scale(m_df[,mod_spec$cc_col])
    m_df$luc_s <- scale(m_df[,mod_spec$luc_col]) 
    m_df$bm_s <- scale(m_df$bm) 
    
    m_df$Managed <- as.factor(m_df$Managed)
    m_df$Utilised <- as.factor(m_df$Utilised)
    
    print("Fitting model to full data")
    ### 
    ## M1 models
    m_nm <- mod_spec$model
    if (m_nm == "m1a"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+(1|spp_idx)+(1|loc_idx),data=m_df,REML=F))
    }
    if (m_nm == "m1b"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+(1|spp_idx)+(1|loc_idx),data=m_df,REML=F))
    }
    
    # M2 models, add realm intercept
    if (m_nm == "m2a"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=m_df,REML=F))
    }
    if (m_nm == "m2b"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=m_df,REML=F))
    }
    
    
    # Edit model call to allow for influence calculation
    tmp_m@call[3] <- "data.update"
    
    # add col for row idx
    m_df$row_idx <- seq(1:nrow(m_df))
    # infl - for pop, spp and loc
    print("Calculating influence")
    pop_infl <- influence.ME::influence(tmp_m, obs = T)
    spp_infl <- influence.ME::influence(tmp_m, group = "spp_idx")
    loc_infl <- influence.ME::influence(tmp_m, group = "loc_idx")
    
    pop_cook <- influence.ME::cooks.distance.estex(pop_infl)
    spp_cook <- influence.ME::cooks.distance.estex(spp_infl)
    loc_cook <- influence.ME::cooks.distance.estex(loc_infl)
    
    pop_cook <- as.data.frame(pop_cook)
    spp_cook <- as.data.frame(spp_cook)
    loc_cook <- as.data.frame(loc_cook)
    
    pop_cook$drop_id <- rownames(pop_cook)
    spp_cook$drop_id <- rownames(spp_cook)
    loc_cook$drop_id <- rownames(loc_cook)
    
    pop_cook <- pop_cook[order(-pop_cook$V1),]
    spp_cook <- spp_cook[order(-spp_cook$V1),]
    loc_cook <- loc_cook[order(-loc_cook$V1),]
    
    pop_cook$idx <- 1:nrow(pop_cook)
    spp_cook$idx <- 1:nrow(spp_cook)
    loc_cook$idx <- 1:nrow(loc_cook)
    
    return(list(
        mod_spec = mod_spec,
        pop_cook = pop_cook,
        spp_cook = spp_cook,
        loc_cook = loc_cook
    ))
}

cooks_plttr <- function(res_ls, lag_type, drop_type){
    plt_df <- bind_rows(lapply(res_ls, function(x, lag_, drop_){
        if (!is.null(x)){
            if (x$mod_spec$type == lag_type){
                tmp <- x[[drop_type]]
                tmp$type <- x$mod_spec$type
                tmp$model <- x$mod_spec$model
                tmp$class <- x$mod_spec$class
                tmp$cc_lag <- x$mod_spec$cc_lag
                tmp$luc_lag <- x$mod_spec$cc_lag
                tmp
            }
        }
        
    }, lag_ = lag_type, drop_ = drop_type))
    ggplot(plt_df) +
        geom_hline(yintercept = 1,
                   lty = 1, colour = "red") +
        geom_hline(yintercept = 0.5,
                   lty = 2, colour = "red") +
        geom_point(aes(x = idx, y = V1)) +
        labs(x = "Rank",
             y = "Cook's distance") +
        facet_grid(class~model)
}



###%%%%%%%%%%
# Main Code
###%%%%%%%%%%

# Load data
# lapply output dir, read rds
# initial_models_test
lag_res <- lapply(list.files("../Results/initial_models",
                             full.names = T), 
                  FUN = readRDS)
# bind together
lag_res <- do.call(rbind.fill, lag_res)

# Run initial analysis
fig_ls <- analyses_wrap_(lag_res)

# Note: during model fitting, m3a and m3b were found to be overparameterised (rank-deficient) 
# for birds.
# We therefore focus on 1a-1b throughout the below.



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extended Data Fig 1
## Distribution of AIC scores - yr v gen
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aic_sub <- lag_res[,c("model", "AICc", "cc_col", "luc_col", "cc_lag", "luc_lag", "class", "type")] %>%
    unique() %>%
    mutate(type_ = case_when(type == "year" ~ "Year",
                             type == "gen" ~ "Generation"),
           model_ = gsub("m", "", model))

ggplot() +
    geom_segment(aes(x = x, xend = x,
                     y = low, yend = high),
                 colour = "grey30",
                 size = 0.9,
                 arrow = arrow(length = unit(0.07, "inch")),
                 data = aic_sub %>%
                     filter(model_ %in% c("1a", "1b", "2a", "2b")) %>%
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
                  filter(model_ %in% c("1a", "1b", "2a", "2b")) %>%
                  group_by(class) %>%
                  summarise(y = (min(AICc) + max(AICc))/2) %>%
                  mutate(x = 0.3,
                         label = "Better")
    ) +
   
    geom_boxplot(aes(x = as.factor(model_), y = AICc, colour = type_),
                 alpha = 0.25,
                 outlier.size = 1,
                 width = 0.1,
                 position = position_nudge(x = -0.075),
                 show.legend = F,
                 data = subset(aic_sub, model_ %in% c("1a", "1b", "2a", "2b") &
                                   type_ == "Year")) +
    
    geom_boxplot(aes(x = as.factor(model_), y = AICc, colour = type_),
                 alpha = 0.25,
                 outlier.size = 1,
                 width = 0.1,
                 position = position_nudge(x = 0.075),
                 show.legend = F,
                 data = subset(aic_sub, model_ %in% c("1a", "1b", "2a", "2b") &
                                   type_ == "Generation")) +
    
    geom_flat_violin_L(aes(x = as.factor(model_), y = AICc, 
                           colour = type_,
                           fill = type_,
                           group = interaction(class, model_)),
                       alpha = 0.25,
                       width = 0.4,
                       position = position_nudge(x = -0.175),
                       data = subset(aic_sub, model_ %in% c("1a", "1b", "2a", "2b") &
                                         type_ == "Year")) +
    geom_flat_violin(aes(x = as.factor(model_), y = AICc, 
                         colour = type_,
                         fill = type_,
                         group = interaction(class, model_)),
                     alpha = 0.25,
                     width = 0.4,
                     position = position_nudge(x = 0.177),
                     data = subset(aic_sub, model_ %in% c("1a", "1b", "2a", "2b") &
                                       type_ == "Generation")) +
    
    facet_grid(class~.,
               scales = "free_y",
               labeller = labeller("class" = c("bird" = "Birds",
                                               "mammal" = "Mammals"))
    ) + 
    
    scale_colour_manual(values = c("Year" = "firebrick3",
                                   "Generation" = "dodgerblue3")) +
    scale_fill_manual(values = c("Year" = "firebrick3",
                                 "Generation" = "dodgerblue3")) +
    
    
    labs(x = "Model structure",
         y = "AICc",
         colour = element_blank(),
         fill = element_blank()) +
    scale_x_discrete(expand = expansion(add = c(0.9,0.6)),
                     breaks = c("1a","1b","2a","2b"),
                     labels = c("Base", "+MU", "+R", "+MUR")
                     ) +
    
    theme_bw() +
    basic_thm +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1/2)

ggsave("../Results/Figs/sm_fig_1.pdf",
       device = "pdf", dpi = 300,
       width = 7.5, height = 5)
# Litte/no benefit of generation-based lags
# Including management and use seems beneficial, especially for mammals

# Some basic plots: Akaike weights per cc/luc lag across model structures
fig_ls$plt_ls$m1a.main$w_yr$w_bar_plt
fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt
fig_ls$plt_ls$m2a.main$w_yr$w_bar_plt
fig_ls$plt_ls$m2b.main$w_yr$w_bar_plt

plot_grid(fig_ls$plt_ls$m1a.main$w_yr$w_bar_plt + 
              theme(legend.position = "none",
                    axis.title = element_blank(),
                    axis.text = element_blank()) +
              facet_grid(~class,
                         labeller = labeller("class" = c("bird" = "Birds",
                                                         "mammal" = "Mammals"))),
          fig_ls$plt_ls$m1b.main$w_yr$w_bar_plt + 
              theme(legend.position = "none",
                    strip.text.x = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank()) +
              facet_grid(~class),
          fig_ls$plt_ls$m2a.main$w_yr$w_bar_plt + 
              theme(legend.position = "none",
                    strip.text.x = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank()) +
              facet_grid(~class),
          fig_ls$plt_ls$m2b.main$w_yr$w_bar_plt + 
              theme(legend.position = "none",
                    strip.text.x = element_blank()
                    # axis.title = element_blank(),
                    # axis.text = element_blank()
                    ) +
              facet_grid(~class),
          # fig_ls_123$res1$plt_ls$m3a.main$w_yr$w_bar_plt + 
          #     theme(legend.position = "none",
          #           strip.text.x = element_blank(),
          #           axis.title = element_blank(),
          #           axis.text = element_blank()) +
          #     facet_grid(~class),
          # fig_ls_123$res1$plt_ls$m3b.main$w_yr$w_bar_plt + 
          #     theme(legend.position = "none",
          #           strip.text.x = element_blank()) +
          #     facet_grid(~class),
          align = "hv",
          ncol = 1,
          labels = c("m1a", "m1b", "m2a", "m2b")
)


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sensitivity analysis - top models
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

infl_assess_specs <- unique(fig_ls$coef_df[fig_ls$coef_df$delta == 0,
                               c("model","class","cc_lag","luc_lag","type",
                                 "cc_col","luc_col","dp_min","R2m","R2c",
                                 "AICc", "delta","id","dur","off","ecol_sub")]) %>%
    filter(type == "year" & model %in% c("m1a", "m1b", "m2a", "m2b"))

infl_assess_specs$id <- 1:nrow(infl_assess_specs)

# load model fitting data...
dat <- readRDS("../Data/anon_dat.rds")


mam_infl <- dlply(subset(infl_assess_specs, class == "mammal"),
                  .(model),
                  function(x, y, z){
                      infl_assess_wrap(y, x, z)
                  },  
                  y = subset(dat$mam_df, lambda_rsq > 0 & DP_to_2014 >= 3),
                  z = dat$env_df)
      
# Cooks D > 0.5 liked to population of a confidential population/species
# (spp_68)

lapply(mam_infl,
       function(x){
           print(x$mod_spec$model)
           print(as.character(x$mod_spec$class))
           print(as.character(x$mod_spec$type))
           
           print(x$pop_cook[x$pop_cook$V1 > 0.5,])
           print(x$spp_cook[x$spp_cook$V1 > 0.5,])
           print(x$loc_cook[x$loc_cook$V1 > 0.5,])
       })
subset(dat$mam_df, loc_idx == "562")
subset(dat$mam_df, lambda_rsq > 0 & DP_to_2014 >= 3)[865,]

bird_infl <- dlply(subset(infl_assess_specs, class == "bird"),
                  .(model),
                  function(x, y, z){
                      infl_assess_wrap(y, x, z)
                  },  
                  y = subset(dat$bird_df, lambda_rsq > 0 & DP_to_2014 >= 3),
                  z = dat$env_df)

# Cooks D > 0.5 liked to Gyps_bengalensis and Podiceps_nigricollis
lapply(bird_infl,
       function(x){
           print(x$mod_spec$model)
           print(as.character(x$mod_spec$class))
           print(as.character(x$mod_spec$type))
           
           print(x$pop_cook[x$pop_cook$V1 > 0.5,])
           print(x$spp_cook[x$spp_cook$V1 > 0.5,])
           print(x$loc_cook[x$loc_cook$V1 > 0.5,])
       })

# For final analysis, remove these influential spp/pops



# Cooks D > 0.5 figs...
cooks_plttr(c(mam_infl, bird_infl), "year", "pop_cook") + theme_bw() + basic_thm
cooks_plttr(c(mam_infl, bird_infl), "year", "spp_cook") + theme_bw() + basic_thm
cooks_plttr(c(mam_infl, bird_infl), "year", "loc_cook") + theme_bw() + basic_thm



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### END
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
