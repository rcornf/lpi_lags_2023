##%%%%%%%%%%
## R script to run analyses of lagged effects of climate and land use 
## Supplementary Analyses 
## Leaave-one-out cross-validation - investigate stability of top models
##%%%%%%%%%%

rm(list = ls())
graphics.off()

# Packages
library(lme4)           #~ LMMs
library(MuMIn)          #~ Model selection
library(car)            #~ Model selection

##%%%%%
## Functions
##%%%%%

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
    
    # df$hyde3.1 <- NA
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

m_fit <- function(df, CI, cc_col, luc_col, m_nm){
    
    # scale cc, luc, bm
    df$cc_s <- scale(df[,cc_col])
    df$luc_s <- scale(df[,luc_col]) 
    df$bm_s <- scale(df$bm) 
    
    df$Managed <- as.factor(df$Managed)
    df$Utilised <- as.factor(df$Utilised)
    
    ## M1 models
    if (m_nm == "m1a"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
        if(inherits(tmp_m, "try-error")){
            return(list(comb_df = NULL,
                        mod = NULL))
        }
    }
    if (m_nm == "m1b"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
        if(inherits(tmp_m, "try-error")){
            return(list(comb_df = NULL,
                        mod = NULL))
        }
    }
    
    # M2 models, add realm intercept
    if (m_nm == "m2a"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
        if(inherits(tmp_m, "try-error")){
            return(list(comb_df = NULL,
                        mod = NULL))
        }
    }
    if (m_nm == "m2b"){
        tmp_m <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
        if(inherits(tmp_m, "try-error")){
            return(list(comb_df = NULL,
                        mod = NULL))
        }
    }
    
    rsq <- as.data.frame(r.squaredGLMM(tmp_m))
    
    aicc <- as.data.frame(AICc(tmp_m))
    colnames(aicc) <- "AICc"
    
    metr_df <- cbind(rsq, aicc)
    
    coef_df <- data.frame(coef_val = fixef(tmp_m))
    coef_df$coef_nm <- rownames(coef_df)
    
    if (CI == T){
        confints <- try(confint(tmp_m))
        if(inherits(confints, "try-error")){
            confints <- as.data.frame(matrix(NA, nrow = nrow(coef_df), ncol = 2))
            confints$coef_nm <- rownames(coef_df)
        }
        else {
            confints <- data.frame(confints)
            confints$coef_nm <- rownames(confints)
        }
        coef_df <- merge(coef_df, confints)
        colnames(coef_df) <- c("coef_nm", "coef_val", "lowCI", "highCI")
    }
    else if (CI == F){
        colnames(coef_df) <- c("coef_val", "coef_nm")
    }
    # coef_df
    if (CI == T){
        coef_df[,5:7] <- ((10^coef_df[c("coef_val", "lowCI", "highCI")]) - 1)*100
    }
    else if (CI == F){
        coef_df[,3] <- ((10^coef_df[c("coef_val")]) - 1)*100
    }
    
    comb_df <- cbind(coef_df, metr_df)
    
    return(comb_df)
}


loo_fit <- function(df, mod_spec){
    out_ls <- vector("list", nrow(df))
    
    for (i in 1:nrow(df)){
        tmp <- df[-i,]
        # df, CI, cc_col, luc_col, m_nm
        tmp_out <- m_fit(tmp, F, mod_spec$cc_col, mod_spec$luc_col, mod_spec$model)
        tmp_out$ID <- df$ID[i]
        # Store m_dat in list
        out_ls[[i]] <- tmp_out
    }
    return(do.call(rbind, out_ls))
}


##

loobsr_fit <- function(df, mod_spec){ #cc_col, luc_col, dp_min){
    df$binom_sys_realm <- paste(df$spp_idx, 
                                df$Sys_realm,
                                sep = "_")
    
    out_ls <- vector("list", length(unique(df$binom_sys_realm)))
    idx_counter <- 1
    for (bsr in unique(df$binom_sys_realm)){
        tmp_df <- subset(df, binom_sys_realm != bsr)
        
        tmp_out <- m_fit(tmp_df, F, mod_spec$cc_col, mod_spec$luc_col, mod_spec$model)
        tmp_out$binom_sys_realm <- bsr
        
        out_ls[[idx_counter]] <- tmp_out
        idx_counter <- idx_counter+1
    }
    
    # join all coefs into single df...
    return(do.call(rbind, out_ls))
}

##

##%%%%%
## Main Code
##%%%%%

# Load data
dat <- readRDS("../Data/anon_dat.rds")
print("Data loaded.")
# Get cv specs
m_spec <- readRDS("../Data/cv_specs.rds")

# Testing
# j <- 2

for (j in 1:nrow(m_spec)){
    f = paste0("lpi_lag_res_", as.character(j), ".rds")
    fp = paste0("../Results/cv_models/", f)
    
    spec <- m_spec[j,]
    
    if (spec$class == "bird"){
        df <- dat$bird_df
        df <- subset(df, lambda_rsq > 0 & DP_to_2014 >=3)
        df <- subset(df, !(spp_idx %in% c("Gyps_bengalensis", "Podiceps_nigricollis")))
    } 
    if (spec$class == "mammal"){
        df <- dat$mam_df
        df <- subset(df, lambda_rsq > 0 & DP_to_2014 >=3)
        df <- subset(df, !(ID %in% c(23514)))
    } 
    
    m_df <- env_rate_calc(df, spec, dat$env_df)
    
    loo <- loo_fit(m_df, spec)
    loobsr <- loobsr_fit(m_df, spec)
    
    loo <- cbind(loo, spec[,c("model", "class", "cc_lag", "luc_lag", "type", 
                              "cc_col", "luc_col", "dp_min", "dur", "off", "delta")])
    loobsr <- cbind(loobsr, spec[,c("model", "class", "cc_lag", "luc_lag", "type", 
                                    "cc_col", "luc_col", "dp_min", "dur", "off", "delta")])
    
    out <- list(loo = loo, loobsr = loobsr)
    
    # Save
    saveRDS(out,
            fp)
}