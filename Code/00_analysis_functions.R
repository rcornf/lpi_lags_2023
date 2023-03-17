# Functions used in the main analysis


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
  
  # round((20+6)/10)*10
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

# Wrapper function to rank models and produce preliminary figures related to lags
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
                    
                    # get best models for each class within data considered
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


# Function to plot a summary of model fitting data
dat_summary_plttr <- function(df){
  # Location summary data
  loc_dat <- ddply(df, .(loc_idx, Class), 
                   function(x){
                     x$Npops <- nrow(x)
                     return(x[1,c("loc_idx", "pa", "Class", "Location", "Country", 
                                  "Latitude", "Longitude", "Npops")])})
  loc_dat$Longitude[loc_dat$Longitude < -170] <- 180+abs(diff(c(-180, loc_dat$Longitude[loc_dat$Longitude < -170])))
  
  loc_plt <- ggplot(loc_dat) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = "lightgrey",
                 colour = "lightgrey",
                 data = world_map_sub) +
    geom_point(aes(x = Longitude, y = Latitude, size = Npops, colour = Class),
               alpha = 0.65) +
    scale_color_manual(values = clrs_,
                       breaks = c("Aves", "Mammalia"),
                       labels = c("Birds", "Mammals")) +
    scale_size_continuous(breaks = c(1,10,50), range = c(1.5,10)) +
    labs(size = "No. of \npopulations",
         colour = "",
         x = element_blank(),
         y = element_blank()) +
    # ylim(c(-58,84)) +
    basic_thm +
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(0.05,0.25),
          legend.text = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes= list(size = 4,
                                                    shape = 15,
                                                    alpha = 0.65)))
  
  mam_pops <- subset(df, Class == "Mammalia")
  bird_pops <- subset(df, Class == "Aves")
  
  # Distr of lambda
  lmbd_plt <- ggplot() +
    geom_histogram(aes(x = lambda3, fill = Class), 
                   alpha = 0.65,
                   bins = 50,
                   data = mam_pops) +
    geom_histogram(aes(x = lambda3, fill = Class), 
                   alpha = 0.65,
                   bins = 50,
                   data = bird_pops) +
    labs(x = expression(bar(lambda)),
         y = "N") +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3")) +
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  # Distr of TS
  ts_plt <- ggplot() +
    geom_histogram(aes(x = TS_length_to_2014, fill = Class), 
                   alpha = 0.65,
                   binwidth = 2,
                   data = mam_pops) +
    geom_histogram(aes(x = TS_length_to_2014, fill = Class), 
                   alpha = 0.65,
                   binwidth = 2,
                   data = bird_pops) +
    labs(x = "Time-series length (years)",
         y = "N") +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3")) +
    scale_x_continuous(limits = c(0,63),
                       breaks = c(0,20,40,60)) +
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  # DP
  dp_plt <- ggplot() +
    geom_histogram(aes(x = DP_to_2014, fill = Class), 
                   alpha = 0.65,
                   binwidth = 2,
                   data = mam_pops) +
    geom_histogram(aes(x = DP_to_2014, fill = Class), 
                   alpha = 0.65,
                   binwidth = 2,
                   data = bird_pops) +
    labs(x = "N data points",
         y = "N") +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3")) +
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  # Generataion length
  gl_plt <- ggplot() +
    geom_histogram(aes(x = GenLength, fill = Class), 
                   alpha = 0.65,
                   binwidth = 1,
                   data = mam_pops) +
    geom_histogram(aes(x = GenLength, fill = Class), 
                   alpha = 0.65,
                   binwidth = 1,
                   data = bird_pops) +
    labs(x = "Generation length (years)",
         y = "N") +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3")) +
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  # Body mass
  bm_plt <- ggplot() +
    geom_histogram(aes(x = bm, fill = Class), 
                   alpha = 0.65,
                   bins = 40,
                   data = mam_pops) +
    geom_histogram(aes(x = bm, fill = Class), 
                   alpha = 0.65,
                   bins = 40,
                   data = bird_pops) +
    labs(x = expression("log"[10]~"Body mass (g)"),
         y = "N") +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3")) +
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  pa_man_use_plt <- bind_rows(mam_pops,
                              bird_pops) %>%
    group_by(Class) %>%
    summarise(y = c(sum(pa == "Yes")/n(),  sum(Managed == 1)/n(),  sum(Utilised == 1)/n()),
              x = c("PA", "Man", "Use")) %>%
    mutate(x = factor(x, levels = c("PA", "Man", "Use"))) %>%
    ggplot() +
    geom_bar(aes(x = x, y = y, group = Class, fill = Class),
             width = 0.3,
             position = 'dodge',
             stat = "identity",
             alpha = 0.65) +
    scale_y_continuous(breaks = c(0,0.5,1), 
                       limits = c(0,1)) +
    scale_fill_manual(values = c("Aves" = "dodgerblue3",
                                 "Mammalia" = "firebrick3"), 
                      guide = "none") +
    labs(x = element_blank(),
         y = "Proportion of pops.") + 
    # theme_bw() +
    theme_classic() +
    basic_thm +
    theme(aspect.ratio = 1,
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
  
  summ_plt <- plot_grid(loc_plt,
                        plot_grid(lmbd_plt, ts_plt,  gl_plt, bm_plt, pa_man_use_plt,
                                  nrow = 1, align = "hv",
                                  labels = c("b.", "c.", "d.", "e.", "f."),
                                  label_size = 18),
                        nrow = 2, rel_heights = c(0.7,0.3),
                        labels = c("a.", ""),
                        label_size = 18)
  
  return(list(loc_plt = loc_plt,
              lmbd_plt = lmbd_plt,
              tsdp_plt = tsdp_plt,
              gen_plt = gen_plt,
              bm_plt = bm_plt,
              pa_man_use_plt = pa_man_use_plt,
              summ_plt = summ_plt))
}


# Function to plot variation in model coefficients over lags
plot_coef_row <- function(df){
  coef_lim <- range(df$coef_val.1)
  
  val1 <- round(min(abs(coef_lim)-1.5))
  breaks_ <- c(-val1, 0, val1)
  p <- ggplot(df) +
    geom_tile(aes(x = lag, y = as.numeric(coef_nm), 
                  colour = coef_val.1, fill = coef_val.1)) +
    scale_color_gradient2(low = "firebrick3",
                          mid = "gainsboro",
                          high = "dodgerblue3",
                          midpoint = 0,
                          limits = coef_lim,
                          breaks = breaks_,
                          guide = guide_colorbar(direction = "horizontal",
                                                 title.position = "left",
                                                 title.vjust = 1,
                                                 label.hjust = 0.5,
                                                 label.vjust = 0.5)) +
    scale_fill_gradient2(low = "firebrick3",
                         mid = "gainsboro",
                         high = "dodgerblue3",
                         midpoint = 0,
                         limits = coef_lim,
                         breaks = breaks_,
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = "left",
                                                title.vjust = 1,
                                                label.hjust = 0.5,
                                                label.vjust = 0.5)) +
    facet_grid(ecol_sub~class,
               labeller = labeller("class" = c("bird" = "Birds",
                                               "mammal" = "Mammals"),
                                   "ecol_sub" = c("main" = "All",
                                                  "bm1" = "Small",
                                                  "bm2" = "Medium",
                                                  "bm3" = "Large"))) + 
    labs(x = "Lag (years)",
         y = element_blank(),
         colour = "Annual pop.\nchange (%)",
         fill = "Annual pop.\nchange (%)") +
    scale_x_continuous(breaks = c(0,20,40)) +
    scale_y_continuous(breaks = seq(1,8),
                       labels = c("Int.", "CC", "LUC", "CC:LUC",
                                  "BM", "PA", "Man", "Use")) +
    coord_equal(expand = F) +
    theme_bw() +
    basic_thm +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          panel.border = element_rect(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black"),
          aspect.ratio = 1)
  
  return(p)
}


# Fit a null model.
m_fit_null <- function(df){
  
  mnull<- try(lmer(lambda ~ (1|spp_idx)+(1|loc_idx),data=df,REML=F))
  if(inherits(mnull, "try-error")){
    mnull <- NA
  }
  
  
  m_ls <- list(mnull = mnull
  )
  
  m_ls <- m_ls[!is.na(m_ls)]
  
  rsq <- as.data.frame(do.call(rbind,lapply(m_ls, r.squaredGLMM)))
  rsq$model <- names(m_ls)
  
  aicc <- as.data.frame(do.call(rbind,lapply(m_ls, AICc)))
  colnames(aicc) <- "AICc"
  aicc$model <- rownames(aicc)
  
  metr_df <- merge(rsq, aicc)
  return(metr_df)
}

## Env rate of cange calc functions
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
    
    df$cru4.04[i] <- cc_dat$cc2
    df$ipsl[i] <- cc_dat$cc3
    
    df$hyde3.2_crng[i] <- luc_dat$luc2
    df$hyde3.2_grz[i] <- luc_dat$luc3
    df$luh2_norng[i] <- luc_dat$luc4
    df$luh2_rng[i] <- luc_dat$luc5
  }
  
  df
}


# Assess inlfuential populations/species/locations
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

# Plot cook's distance
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



luh2_extract <- function(loc_df){
  luh2 <- raster::stack(brick("../Data/LUH2/states.nc", var = "c3ann"),
                        brick("../Data/LUH2/states.nc", var = "c4ann"),
                        brick("../Data/LUH2/states.nc", var = "c3per"),
                        brick("../Data/LUH2/states.nc", var = "c4per"),
                        brick("../Data/LUH2/states.nc", var = "c3nfx"),
                        brick("../Data/LUH2/states.nc", var = "pastr"),
                        brick("../Data/LUH2/states.nc", var = "range")
  )
  
  # Index,  name,   year
  # 1,      0,      850
  # 1052,   1051,   1901
  # 1166,   1165,   2015
  names(luh2) <- paste0(rep(seq(850,2015), 7), 
                        "_", 
                        rep(c("c3ann", "c4ann", "c3per", "c4per", "c3nfx",
                              "pastr", "range"), each = 1166))
  
  luh2_1900_2015 <- luh2[[as.vector(sapply(seq(1051, by = 1166, length.out = 7), 
                                           FUN = function(x){seq(x,x+115)}))]]
  
  luh2_df <- data.frame(extract(luh2_1900_2015, 
                                loc_df[,c("x", "y")],
                                ncol=2))
  rownames(luh2_df) <- loc_df$loc_idx
  
  luh2_df <- t(luh2_df)
  luh2_yr_var <- data.frame(do.call(rbind, strsplit(rownames(luh2_df), "_")))
  colnames(luh2_yr_var) <- c("Year", "LU")
  luh2_yr_var$Year <- as.numeric(gsub("X", "", luh2_yr_var$Year))
  
  luh2_long_df <- reshape2::melt(cbind(luh2_df, luh2_yr_var), 
                                 id.vars = c("Year", "LU"),
                                 variable.name = "loc_idx",
                                 value.name = "prop_cover")
  
  # make luh2 wide - col per lu type
  luh2_wide_df <- tidyr::spread(luh2_long_df, LU, prop_cover)
  
  luh2_wide_df_ <- na.omit(luh2_wide_df)
  luh2_wide_df_
}

# Extract temperature data
ipsl_extract <- function(loc_df){
  ipsl_f_ls <- list.files(path = "../Data/IPSL",
                          full.names = T)
  
  ipsl_ls <- purrr::map(ipsl_f_ls,
                        function(x, loc){
                          tmp_brk <- brick(x)
                          
                          tmp_terra <- terra::rast(tmp_brk)
                          names(tmp_terra) <- names(tmp_brk)
                          
                          # Option 1
                          # extract, mean
                          tmp_t <- data.frame(terra::extract(tmp_terra, as.matrix(loc[,c("x", "y")]))) 
                          
                          rownames(tmp_t) <- loc$loc_idx
                          # need to know unique years...
                          tmp_yrs <- sort(as.numeric(gsub("X", "", 
                                                          unique(sapply(strsplit(names(tmp_terra), "\\."), getElement, 1)))))
                          
                          tmp_t <- data.frame(sapply(tmp_yrs, 
                                                     FUN = function(x){rowMeans(tmp_t[,grep(x, names(tmp_t))])}))
                          
                          colnames(tmp_t) <- tmp_yrs
                          tmp_t$loc_idx <- rownames(tmp_t)
                          tmp_t <- reshape2::melt(tmp_t, id.vars = "loc_idx", 
                                                  variable.name = "Year",
                                                  value.name = "av_t_K")
                          tmp_t$Year <- as.numeric(as.character(tmp_t$Year))
                          tmp_t$ipsl <- tmp_t$av_t_K - 273.15
                          tmp_t
                        },
                        loc_df)
  
  ipsl_df <- do.call(rbind, ipsl_ls)
  ipsl_df
}

# Function to extraact land-use data given specified cell centres/locations
luh2_fut_extract <- function(loc_df, f_nm){
  fp <- paste0("../Data/LUH2/", f_nm)
  
  luh2 <- raster::stack(brick(fp, var = "c3ann"),
                        brick(fp, var = "c4ann"),
                        brick(fp, var = "c3per"),
                        brick(fp, var = "c4per"),
                        brick(fp, var = "c3nfx"),
                        brick(fp, var = "pastr"),
                        brick(fp, var = "range")
  )
  
  names(luh2) <- paste0(rep(seq(2015,2100), 7), 
                        "_", 
                        rep(c("c3ann", "c4ann", "c3per", "c4per", "c3nfx",
                              "pastr", "range"), each = 86))
  
  
  luh2_df <- data.frame(extract(luh2, 
                                loc_df[,c("x", "y")],
                                ncol=2))
  rownames(luh2_df) <- loc_df$loc_idx
  
  luh2_df <- t(luh2_df)
  luh2_yr_var <- data.frame(do.call(rbind, strsplit(rownames(luh2_df), "_")))
  colnames(luh2_yr_var) <- c("Year", "LU")
  luh2_yr_var$Year <- as.numeric(gsub("X", "", luh2_yr_var$Year))
  
  luh2_long_df <- reshape2::melt(cbind(luh2_df, luh2_yr_var), 
                                 id.vars = c("Year", "LU"),
                                 variable.name = "loc_idx",
                                 value.name = "prop_cover")
  
  # make luh2 wide - col per lu type
  luh2_wide_df <- tidyr::spread(luh2_long_df, LU, prop_cover)
  
  luh2_wide_df_ <- na.omit(luh2_wide_df)
  
  luh2_wide_df_$luh2_norng <- base::rowSums(luh2_wide_df_[,c("c3ann", "c3nfx", "c3per",
                                                             "c4ann", "c4per", "pastr")])
  
  luh2_wide_df_$luh2_rng <- base::rowSums(luh2_wide_df_[,c("c3ann", "c3nfx", "c3per",
                                                           "c4ann", "c4per", "pastr",
                                                           "range")])
  luh2_wide_df_
}

# Extract temperature data
ipsl_fut_extract <- function(loc_df, f_dir){
  ipsl_f_ls <- list.files(path = paste0("../Data/", f_dir),
                          full.names = T)
  
  ipsl_ls <- purrr::map(ipsl_f_ls,
                        function(x, loc){
                          tmp_brk <- brick(x)
                          
                          tmp_terra <- terra::rast(tmp_brk)
                          names(tmp_terra) <- names(tmp_brk)
                          
                          # extract, mean
                          tmp_t <- data.frame(terra::extract(tmp_terra, as.matrix(loc[,c("x", "y")]))) 
                          
                          rownames(tmp_t) <- loc$loc_idx
                          # need to know unique years...
                          tmp_yrs <- sort(as.numeric(gsub("X", "", 
                                                          unique(sapply(strsplit(names(tmp_terra), "\\."), getElement, 1)))))
                          
                          tmp_t <- data.frame(sapply(tmp_yrs, 
                                                     FUN = function(x){rowMeans(tmp_t[,grep(x, names(tmp_t))])}))
                          
                          colnames(tmp_t) <- tmp_yrs
                          tmp_t$loc_idx <- rownames(tmp_t)
                          tmp_t <- reshape2::melt(tmp_t, id.vars = "loc_idx", 
                                                  variable.name = "Year",
                                                  value.name = "av_t_K")
                          tmp_t$Year <- as.numeric(as.character(tmp_t$Year))
                          tmp_t$ipsl <- tmp_t$av_t_K - 273.15
                          tmp_t
                        },
                        loc_df)
  
  ipsl_df <- do.call(rbind, ipsl_ls)
  ipsl_df
}



# function to extract env rates based on pop df, fut_env_ls, and models (lags...)
fut_env_rate_calc_yrs <- function(lags, fut_env_ls, env_){
  # uniqe locs
  locs <- unique(fut_env_ls[[1]]$loc_idx)
  sc_ <- names(fut_env_ls)
  decs <- c("2010_2020", "2020_2030", "2030_2040", "2040_2050",
            "2050_2060", "2060_2070", "2070_2080", "2080_2090", "2090_2100",
            "2100_2110", "2110_2120", "2120_2130", "2130_2140")
  
  out_df <- expand.grid(loc_id = locs,
                        sc = sc_,
                        Dec = decs,
                        lag = lags,
                        type = "year",
                        stringsAsFactors = F)
  
  
  out_df$TS_strt <- as.numeric(sapply(strsplit(out_df$Dec, "_"), getElement, 1))
  out_df$TS_end <- as.numeric(sapply(strsplit(out_df$Dec, "_"), getElement, 2))
  
  if (env_ == "CC"){
    out_df$cc_strt <- out_df$TS_strt - out_df$lag
    out_df$cc_end <- out_df$TS_end - out_df$lag
    
    out_df$ipsl <- NA
    
    for (i in 1:nrow(out_df)){
      if (out_df$cc_end[i]<=2100){
        # cc calc
        cc_dat <- cc_calc(out_df$cc_strt[i],
                          out_df$cc_end[i],
                          out_df$loc_idx[i],
                          fut_env_ls[[out_df$sc[i]]])
        out_df$ipsl[i] <- cc_dat$cc3
      }
    }
  }
  if (env_ == "LUC"){
    
    out_df$luc_strt <- out_df$TS_strt - out_df$lag
    out_df$luc_end <- out_df$TS_end - out_df$lag
    
    out_df$luh2_norng <- NA
    out_df$luh2_rng <- NA
    for (i in 1:nrow(out_df)){
      if (out_df$luc_end[i]<=2100){
        # luc calc
        luc_dat <- luc_calc(out_df$luc_strt[i],
                            out_df$luc_end[i],
                            out_df$loc_idx[i],
                            fut_env_ls[[out_df$sc[i]]])
        out_df$luh2_norng[i] <- luc_dat$luc4
        out_df$luh2_rng[i] <- luc_dat$luc5
      }
    }
  }
  
  return(out_df)
}



### function to plot pop specific graphic given pop id...
pop_graphic_plttr <- function(pop_id){
  pop_row <- subset(lpd, ID == pop_id)
  tmp_pop_df <- data.frame("Year" = 1950:2014,
                           "abund" = unname(t(pop_row[1,c(which(colnames(pop_row) == "X1950"):
                                                            which(colnames(pop_row) == "X2014"))])))
  
  tmp_pop_df <- na.omit(tmp_pop_df)
  tmp_env_df <- subset(dat$env_df, loc_idx == mam_df2$loc_idx[mam_df2$ID == pop_id])
  
  bind_rows(tmp_pop_df %>%
              mutate(var = "Abundance",
                     val = abund),
            tmp_env_df %>%
              filter(Year %in% c((min(tmp_pop_df$Year)-50):max(tmp_pop_df$Year))) %>%
              mutate(var = "Temperature",
                     val = ipsl),
            tmp_env_df %>%
              filter(Year %in% c((min(tmp_pop_df$Year)-50):max(tmp_pop_df$Year))) %>%
              mutate(var = "Human land",
                     val = luh2_norng)
            
  ) %>%
    ggplot() +
    geom_point(aes(x = Year, y = val, colour = var),
               alpha = 0.5) +
    # pop vline and gam
    geom_smooth(aes(x = Year, y = val, colour = var),
                size = 2,
                method = "gam",
                se = F,
                data = tmp_pop_df %>%
                  mutate(val = abund,
                         var = "Abundance")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year),
                                xmax = max(tmp_pop_df$Year),
                                var = "Abundance"),
              show.legend = F
    ) +
    # temp vline and lm
    geom_smooth(aes(x = Year, y = val, colour = var),
                size = 2,
                method = "lm",
                se = F,
                data = subset(tmp_env_df, Year %in% c(tmp_pop_df$Year-45)) %>%
                  mutate(val = ipsl,
                         var = "Temperature")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year)-45,
                                xmax = max(tmp_pop_df$Year)-45,
                                var = "Temperature"),
              show.legend = F
    ) +
    # luc vline and line
    geom_line(aes(x = Year, y = val, colour = var),
              size = 2,
              data = subset(tmp_env_df, Year %in% c(tmp_pop_df$Year-9)) %>%
                mutate(val = luh2_norng,
                       var = "Human land")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year)-9,
                                xmax = max(tmp_pop_df$Year)-9,
                                var = "Human land"),
              show.legend = F
    ) +
    
    # temp and land use lag error bar
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, y = 15),
                   height = 3,
                   colour = "firebrick3",
                   size = 1.5,
                   data = data.frame(xmin = min(tmp_pop_df$Year)-45,
                                     xmax = min(tmp_pop_df$Year),
                                     var = c("Abundance"))
    ) +
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, y = 45),
                   height = 3,
                   colour = "dodgerblue3",
                   size = 1.5,
                   data = data.frame(xmin = min(tmp_pop_df$Year)-9,
                                     xmax = min(tmp_pop_df$Year),
                                     var = c("Abundance"))
    ) +
    
    geom_text(aes(x = x, y = y, label = label),
              colour = "firebrick3",
              size = 5.5,
              data = data.frame(x = c(1927.5),
                                y = c(20),
                                label = "45 years",
                                var = c("Abundance"))
    ) +
    # 
    geom_text(aes(x = x, y = y, label = label),
              colour = "dodgerblue3",
              size = 5.5,
              data = data.frame(x = c(1941),
                                y = c(40),
                                label = "9 years",
                                var = c("Abundance"))
    ) +
    
    scale_color_manual(values = c("Abundance" = "black",
                                  "Temperature" = "firebrick3",
                                  "Human land" = "dodgerblue3"),
                       guide = F) +
    scale_fill_manual(values = c("Abundance" = "black",
                                 "Temperature" = "firebrick3",
                                 "Human land" = "dodgerblue3"),
                      guide = F) +
    
    labs(x = element_blank(),
         y = element_blank()) +
    
    scale_x_continuous(breaks = c(1900,1950,2000)) +
    facet_grid(var~.,
               scales = "free_y",
               as.table = F,
               labeller = labeller(
                 "var" = c("Abundance" = "Abundance\n",
                           "Human land" = "Anthropogenic land\n",
                           "Temperature" = "Temperature (\u00b0C)\n"))
    ) +
    basic_thm 
}


# Predict pop trends in response to env change across a range of conditions
predict_env_response <- function(model_coef_df){
  pred_df <- expand.grid(cc_s = seq(-2.5,2.5, length.out = 100),
                         luc_s = seq(-2.5,2.5, length.out = 100)) %>%
    as.data.frame() %>%
    mutate("cc_s:luc_s" = cc_s*luc_s,
           "int" = 1)
  
  pred_df$pred <- pred_df[,c("cc_s", "luc_s", "cc_s:luc_s")] %>%
    as.matrix() %*% 
    c(subset(model_coef_df, coef_nm == "cc_s")[,"coef_val"],
      subset(model_coef_df, coef_nm == "luc_s")[,"coef_val"],
      subset(model_coef_df, coef_nm == "luc_s:cc_s")[,"coef_val"])
  
  pred_df$pred_int <- pred_df$pred + subset(model_coef_df, coef_nm == "(Intercept)")[,"coef_val"]
  
  pred_df$pred_nointer <- pred_df[,c("cc_s", "luc_s")] %>%
    as.matrix() %*% 
    c(subset(model_coef_df, coef_nm == "cc_s")[,"coef_val"],
      subset(model_coef_df, coef_nm == "luc_s")[,"coef_val"])
  
  pred_df$pred_cc <- pred_df$cc_s * subset(model_coef_df, coef_nm == "cc_s")[,"coef_val"]
  pred_df$pred_luc <- pred_df$luc_s * subset(model_coef_df, coef_nm == "luc_s")[,"coef_val"]
  
  return(pred_df)
}

# Get model-averaged population responses to environmental change combinations
env_response_pred_av_wrap <- function(coef_df, metr_df){
  
  # for each m, get sub df, 
  env_pred_ls <- vector("list", nrow(metr_df))
  
  for (i in 1:nrow(metr_df)){
    # predict 
    tmp_coef_df <- subset(coef_df, delta == metr_df$delta[i])
    tmp_pred_df <- predict_env_response(tmp_coef_df)
    # add delta, 
    tmp_pred_df$delta <- metr_df$delta[i]
    env_pred_ls[[i]] <- tmp_pred_df
  }
  
  env_pred_df <- bind_rows(env_pred_ls)
  
  # calc w
  w_sum <- sum( exp( -0.5 * metr_df$delta  ) )
  env_pred_df$w <- exp(-0.5 * env_pred_df$delta)/w_sum
  
  env_pred_df_av <- env_pred_df %>%
    dplyr::group_by(cc_s, luc_s) %>%
    dplyr::summarise(av_pred = stats::weighted.mean(pred, w),
                     av_pred_int = stats::weighted.mean(pred_int, w),
                     av_pred_nointer = stats::weighted.mean(pred_nointer, w),
                     av_pred_cc = stats::weighted.mean(pred_cc, w),
                     av_pred_luc = stats::weighted.mean(pred_luc, w)) %>%
    ungroup()
  
  return(env_pred_df_av)
}


# Need to fit models - all, bm1, bm2, bm3
lmm_bm_fit_store <- function(mod_df, pop_df, env_df){
  
  pop_df$Managed <- as.factor(pop_df$Managed)
  pop_df$Utilised <- as.factor(pop_df$Utilised)
  
  print(sprintf("N models: %d", nrow(mod_df)))
  
  out <- vector("list", nrow(mod_df))
  
  bm_quants <- quantile(pop_df$bm,
                        c(1/3, 2/3))
  
  # for each mod
  for (i in 1:nrow(mod_df)){
    print(sprintf("Fitting model %d", i))
    
    # subset to bm_split
    # Note: if ecol_sub == main, model for all pops is fitted
    if (mod_df$ecol_sub[i] == "bm1"){
      df <- subset(pop_df, bm <= bm_quants[1])
    }
    if (mod_df$ecol_sub[i] == "bm2"){
      df <- subset(pop_df, bm <= bm_quants[2] & 
                     bm > bm_quants[1])
    }
    if (mod_df$ecol_sub[i] == "bm3"){
      df <- subset(pop_df, bm > bm_quants[2])
    }
    # prep model data
    df <- env_rate_calc(df, mod_df[i,], env_df)
    
    # Full data and model...
    tmp_full <- cbind(df, mod_df[i,c("model", "class", "cc_lag", "luc_lag",
                                     "cc_col", "luc_col", "delta")])
    # scale
    tmp_full$cc_s <- scale(tmp_full[,mod_df$cc_col[i]])
    tmp_full$luc_s <- scale(tmp_full[,mod_df$luc_col[i]])
    tmp_full$bm_s <- scale(tmp_full$bm)
    
    m_full <- lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+(1|spp_idx)+(1|loc_idx),
                   data=tmp_full,REML=F)
    
    rsq <- as.data.frame(r.squaredGLMM(m_full))
    # colnames(rsq) <- paste0(colnames(rsq), "_tr")
    
    tmp_coef_full <- data.frame(coef_val = fixef(m_full))
    tmp_coef_full$coef_nm <- rownames(tmp_coef_full)
    
    # Should get confints...
    confints <- try(confint(m_full))
    if(inherits(confints, "try-error")){
      confints <- as.data.frame(matrix(NA, nrow = nrow(tmp_coef_full), ncol = 2))
      confints$coef_nm <- rownames(tmp_coef_full)
    }
    else {
      confints <- data.frame(confints)
      confints$coef_nm <- rownames(confints)
    }
    tmp_coef_full <- merge(tmp_coef_full, confints)
    colnames(tmp_coef_full) <- c("coef_nm", "coef_val", "loCI", "hiCI")
    # tmp_coef_full[,3] <- ((10^tmp_coef_full[,c("coef_val")]) - 1)*100
    tmp_coef_full[,5:7] <- ((10^tmp_coef_full[c("coef_val", "loCI", "hiCI")]) - 1)*100
    tmp_coef_full <- cbind(tmp_coef_full,
                           mod_df[i, c("model", "class", "cc_lag", "luc_lag",
                                       "cc_col", "luc_col", "delta")])
    
    tmp_coef_full$AICc <- AICc(m_full)
    tmp_coef_full <- cbind(tmp_coef_full, rsq)
    
    
    out[[i]]$df <- tmp_full
    out[[i]]$coef <- tmp_coef_full
    out[[i]]$mod_obj <- m_full
    out[[i]]$mod_spec <- mod_df[i,]
  }
  
  return(list(mod_spec = mod_df,
              mod_dat = out))
}


mod_proj <- function(lmm_dat, env_fut){
  
  sc_proj_ls <- vector("list", length(lmm_dat))
  
  print("Calculating future cc")
  fut_cc <- fut_env_rate_calc_yrs(unique(lmm_dat$mod_spec$cc_lag), 
                                  env_fut,
                                  "CC") 
  print("Calculating future luc")
  fut_luc <- fut_env_rate_calc_yrs(unique(lmm_dat$mod_spec$luc_lag), 
                                   env_fut,
                                   "LUC") 
  
  # for each m,
  # %%%
  print(sprintf("N models: %d", nrow(lmm_dat$mod_spec)))
  
  # for each mod
  for (i in 1:nrow(lmm_dat$mod_spec)){
    print(sprintf("Projecting model %d", i))
    mod_dat_ <- lmm_dat$mod_dat[[i]]
    
    # get cc/luc from fut/hist dfs
    # combine rates...
    # rbind hist and fut cc
    tmp_cc <- subset(fut_cc,
                     lag == mod_dat_$mod_spec$cc_lag)
    # rbind hist and fut luc
    tmp_luc <- subset(fut_luc,
                      lag == mod_dat_$mod_spec$luc_lag)
    # then join 
    tmp_env <- left_join(tmp_cc[,c("loc_idx", "Dec", "type", "TS_strt", "TS_end", 
                                   "cc_strt", "cc_end", "ipsl", "sc")],
                         tmp_luc[,c("loc_idx", "Dec", "type", "TS_strt", "TS_end", 
                                    "luc_strt", "luc_end", "luh2_norng", "luh2_rng", "sc")],
                         by = c("loc_idx", "Dec", "type", "TS_strt", "TS_end", "sc"))
    
    # get scaling factors for cc/luc and scale
    tmp_env$cc_s <- scale(tmp_env[,mod_dat_$mod_spec$cc_col],
                          center = attr(mod_dat_$df$cc_s, "scaled:center"),
                          scale = attr(mod_dat_$df$cc_s, "scaled:scale"))  
    
    tmp_env$luc_s <- scale(tmp_env[,mod_dat_$mod_spec$luc_col],
                           center = attr(mod_dat_$df$luc_s, "scaled:center"),
                           scale = attr(mod_dat_$df$luc_s, "scaled:scale"))
    
    
    # pop_df is: mod_dat_$df
    # model is: mod_dat_$mod_obj
    
    # then, use model to predict for sc data with pop info
    # combine pop_pred_df (reduced) with loc_id_key
    pop_proj_df <- mod_dat_$df[,c("ID", "loc_idx", "luc_lag", "cc_lag",
                                  "LPI_Realm", "spp_idx", "bm_s", "pa", "Managed",
                                  "Utilised", "class", "cc_col", "luc_col", "delta")]
    
    # then inner join change rates to pop_pred_df_
    pop_proj_df <- inner_join(pop_proj_df,
                              tmp_env)
    
    # pred with/without ranef
    pop_proj_df$pred_incl_ranef <- lme4:::predict.merMod(mod_dat_$mod_obj,
                                                         pop_proj_df) 
    pop_proj_df$pred_excl_ranef <- lme4:::predict.merMod(mod_dat_$mod_obj,
                                                         pop_proj_df,
                                                         re.form = NA) 
    # add dfs to output list
    sc_proj_ls[[i]] <- pop_proj_df
    
  }
  
  # Aggregate /weight outputs to give LPI-style indices
  print("Aggregating projections")
  # Assign akaike weights to model predictions
  w_sc_proj_df <- bind_rows(sc_proj_ls)
  
  w_sum <- sum( exp( -0.5 * lmm_dat$mod_spec$delta ) )
  w_sc_proj_df$w <- exp(-0.5 * w_sc_proj_df$delta)/w_sum
  
  
  # %%%%%%%%%%%%%%%%%
  
  # LPI-style index per model
  spp_sc_proj_df <- w_sc_proj_df %>%
    group_by(spp_idx, LPI_Realm, class, Dec, TS_strt, TS_end, sc, delta, w) %>%
    dplyr::summarise(pred_incl_ranef = mean(pred_incl_ranef),
                     pred_excl_ranef = mean(pred_excl_ranef)) %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(spp_idx, LPI_Realm, class, sc, delta, w) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  realm_sc_proj_df <- spp_sc_proj_df %>%
    group_by(LPI_Realm, class, Dec, TS_strt, TS_end, sc, delta, w) %>%
    dplyr::summarise(pred_incl_ranef = mean(pred_incl_ranef),
                     pred_excl_ranef = mean(pred_excl_ranef)) %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(LPI_Realm, class, sc, delta, w) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  sc_proj_df <- realm_sc_proj_df %>%
    group_by(class, Dec, TS_strt, TS_end, sc, delta, w) %>%
    dplyr::summarise(pred_incl_ranef = mean(pred_incl_ranef),
                     pred_excl_ranef = mean(pred_excl_ranef)) %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(class, sc, delta, w) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  # LPI-style index - av over models
  w_av_spp_sc_proj_df <- spp_sc_proj_df %>%
    group_by(spp_idx, LPI_Realm, class, Dec, TS_strt, TS_end, sc) %>%
    dplyr::summarise(pred_incl_ranef = stats::weighted.mean(pred_incl_ranef, w, na.rm = T),
                     pred_excl_ranef = stats::weighted.mean(pred_excl_ranef, w, na.rm = T)) %>%
    ungroup() %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(spp_idx, LPI_Realm, class, sc) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  w_av_realm_sc_proj_df <- realm_sc_proj_df %>%
    group_by(LPI_Realm, class, Dec, TS_strt, TS_end, sc) %>%
    dplyr::summarise(pred_incl_ranef = stats::weighted.mean(pred_incl_ranef, w, na.rm = T),
                     pred_excl_ranef = stats::weighted.mean(pred_excl_ranef, w, na.rm = T)) %>%
    ungroup() %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(LPI_Realm, class, sc) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  
  w_av_sc_proj_df <- sc_proj_df %>%
    group_by(class, Dec, TS_strt, TS_end, sc) %>%
    dplyr::summarise(pred_incl_ranef = stats::weighted.mean(pred_incl_ranef, w, na.rm = T),
                     pred_excl_ranef = stats::weighted.mean(pred_excl_ranef, w, na.rm = T)) %>%
    ungroup() %>%
    mutate(pred_incl_ranef_mult = 10^(pred_incl_ranef*10),
           pred_excl_ranef_mult = 10^(pred_excl_ranef*10)) %>%
    group_by(class, sc) %>%
    arrange(TS_end) %>%
    mutate(pred_incl_ranef_Idx = cumprod(pred_incl_ranef_mult),
           pred_excl_ranef_Idx = cumprod(pred_excl_ranef_mult))
  
  
  return(list(w_sc_proj_df = w_sc_proj_df,
              
              spp_sc_proj_df = spp_sc_proj_df,
              realm_sc_proj_df = realm_sc_proj_df,
              sc_proj_df = sc_proj_df,
              
              w_av_spp_sc_proj_df = w_av_spp_sc_proj_df,
              w_av_realm_sc_proj_df = w_av_realm_sc_proj_df,
              w_av_sc_proj_df = w_av_sc_proj_df))
}


# Plotter function for population data graphic
pop_graphic_plttr <- function(pop_id, lpd, env_df, mod_df){
  pop_row <- subset(lpd, ID == pop_id)
  tmp_pop_df <- data.frame("Year" = 1950:2014,
                           "abund" = unname(t(pop_row[1,c(which(colnames(pop_row) == "X1950"):
                                                            which(colnames(pop_row) == "X2014"))])))
  
  tmp_pop_df <- na.omit(tmp_pop_df)
  tmp_env_df <- subset(env_df, loc_id == mod_df$loc_id[mod_df$ID == pop_id])
  
  tmp_pop_pred = data.frame("Year" = seq(min(tmp_pop_df$Year), max(tmp_pop_df$Year)))
  smth_par <- max(round(nrow(tmp_pop_df)/2),3)
  s <- mgcv::`s`
  tmp_gam <- mgcv::gam(log10(abund)~s(Year, k = smth_par), fx = T, data = tmp_pop_df)
  # Predict from pop start to end_yr
  
  tmp_pop_pred$abund_int <- predict(tmp_gam, tmp_pop_pred)
  
  
  bind_rows(tmp_pop_df %>%
              mutate(var = "Abundance",
                     val = abund),
            tmp_env_df %>%
              filter(Year %in% c((min(tmp_pop_df$Year)-50):max(tmp_pop_df$Year))) %>%
              mutate(var = "Temperature",
                     val = ipsl),
            tmp_env_df %>%
              filter(Year %in% c((min(tmp_pop_df$Year)-50):max(tmp_pop_df$Year))) %>%
              mutate(var = "Human land",
                     val = luh2_norng)
            
  ) %>%
    ggplot() +
    geom_point(aes(x = Year, y = val, colour = var),
               alpha = 0.5) +
    # pop vline and gam
    geom_line(aes(x = Year, y = 10^val, colour = var),
              size = 2,
              method = "gam",
              se = F,
              data = tmp_pop_pred %>%
                mutate(val = abund_int,
                       var = "Abundance")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year),
                                xmax = max(tmp_pop_df$Year),
                                var = "Abundance"),
              show.legend = F
    ) +
    # temp vline and lm
    geom_smooth(aes(x = Year, y = val, colour = var),
                size = 2,
                method = "lm",
                se = F,
                data = subset(tmp_env_df, Year %in% c(tmp_pop_df$Year-45)) %>%
                  mutate(val = ipsl,
                         var = "Temperature")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year)-45,
                                xmax = max(tmp_pop_df$Year)-45,
                                var = "Temperature"),
              show.legend = F
    ) +
    # luc vline and line
    geom_line(aes(x = Year, y = val, colour = var),
              size = 2,
              data = subset(tmp_env_df, Year %in% c(tmp_pop_df$Year-9)) %>%
                mutate(val = luh2_norng,
                       var = "Human land")) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = var),
              alpha = 0.125,
              
              data = data.frame(xmin = min(tmp_pop_df$Year)-9,
                                xmax = max(tmp_pop_df$Year)-9,
                                var = "Human land"),
              show.legend = F
    ) +
    
    # temp and land use lag error bar
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, y = 15),
                   height = 6,
                   colour = "firebrick3",
                   size = 1.5,
                   data = data.frame(xmin = min(tmp_pop_df$Year)-45,
                                     xmax = min(tmp_pop_df$Year),
                                     var = c("Abundance"))
    ) +
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, y = 45),
                   height = 6,
                   colour = "dodgerblue3",
                   size = 1.5,
                   data = data.frame(xmin = min(tmp_pop_df$Year)-9,
                                     xmax = min(tmp_pop_df$Year),
                                     var = c("Abundance"))
    ) +
    
    geom_text(aes(x = x, y = y, label = label),
              colour = "firebrick3",
              size = 5.5,
              data = data.frame(x = c(1927.5),
                                y = c(23),
                                label = "45 years",
                                var = c("Abundance"))
    ) +
    # 
    geom_text(aes(x = x, y = y, label = label),
              colour = "dodgerblue3",
              size = 5.5,
              data = data.frame(x = c(1941),
                                y = c(55),
                                label = "9 years",
                                var = c("Abundance"))
    ) +
    
    scale_color_manual(values = c("Abundance" = "black",
                                  "Temperature" = "firebrick3",
                                  "Human land" = "dodgerblue3"),
                       guide = F) +
    scale_fill_manual(values = c("Abundance" = "black",
                                 "Temperature" = "firebrick3",
                                 "Human land" = "dodgerblue3"),
                      guide = F) +
    
    labs(x = element_blank(),
         y = element_blank()) +
    
    scale_x_continuous(breaks = c(1900,1950,2000)) +
    facet_grid(var~.,
               scales = "free_y",
               as.table = F,
               labeller = labeller(
                 "var" = c("Abundance" = "Abundance\n",
                           "Human land" = "Anthropogenic land\n",
                           "Temperature" = "Temperature (\u00b0C)\n"))
    ) +
    basic_thm 
}
