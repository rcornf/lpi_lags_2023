##%%%%%%%%%%
## R script to run analyses of lagged effects of climate and land use
## Runs initial models considering all popuations with at least 3 monitoring years
## and a trend fit with Rsq > 0
## Uses LUH2 and IPSL environmenta data
##%%%%%%%%%%

#  NOTE: for actual analysis, an HPC was used. Code here is re-written to run on
#  a single, local core.


rm(list = ls())
graphics.off()

# Packages
library(lme4)           #~ LMMs
library(MuMIn)          #~ Model selection
library(car)            #~ Model selection


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cc_calc <- function(strt, end, loc){
    tmp_env <- subset(env_df, loc_idx == loc &
                          Year %in% (strt:end))
    strt_end <- cc1 <- cc2 <- cc3 <- NA
    if (nrow(tmp_env)>1){
        tmp_env <- tmp_env[order(tmp_env$Year),]
        strt_end <- paste(range(tmp_env$Year), collapse = "_")
        
        # calc rates
        if (sum(!is.na(tmp_env$cru4.04))>=2){
            cc2 <- lm(cru4.04~Year, data = tmp_env)$coefficients[2]
        }
        if (sum(!is.na(tmp_env$ipsl))>=2){
            cc3 <- lm(ipsl~Year, data = tmp_env)$coefficients[2]
        }
    }
    return(list(strt_end = strt_end, 
                cc1 = cc1, cc2 = cc2, cc3 = cc3))
}


luc_calc <- function(strt, end, loc){
    tmp_env <- subset(env_df, loc_idx == loc &
                          Year %in% (strt:end))
    strt_end <- luc1 <- luc2 <- luc3 <- luc4 <- luc5 <- NA
    if (nrow(tmp_env)>1){
        # order env dfs by year
        tmp_env <- tmp_env[order(tmp_env$Year),]
        strt_end <- paste(range(tmp_env$Year), collapse = "_")
        
        # calc rates
        luc2 <- mean(diff(tmp_env$hyde3.2_crng))
        luc3 <- mean(diff(tmp_env$hyde3.2_grz))
        luc4 <- mean(diff(tmp_env$luh2_norng))
        luc5 <- mean(diff(tmp_env$luh2_rng))
        
    }
    return(list(strt_end = strt_end, 
                luc1 = luc1, luc2 = luc2, luc3 = luc3, luc4 = luc4, luc5 = luc5))
}




m_fit <- function(df, CI, cc_col, luc_col, dp_min){
    
    # scale cc, luc, bm
    df$cc_s <- scale(df[,cc_col])
    df$luc_s <- scale(df[,luc_col]) 
    df$bm_s <- scale(df$bm) 
    
    df$Managed <- as.factor(df$Managed)
    df$Utilised <- as.factor(df$Utilised)
    
    
    ## M1 models
    m1a <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m1a, "try-error")){
        m1a <- NA
    }
    
    m1b <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m1b, "try-error")){
        m1b <- NA
    }
    ##
    
    # M2 models, add realm intercept
    ##
    m2a <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m2a, "try-error")){
        m2a <- NA
    }
    
    m2b <- try(lmer(lambda ~ luc_s*cc_s+bm_s+pa+Managed+Utilised+LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m2b, "try-error")){
        m2b <- NA
    }
    ##
    # M3 models, realm interactions
    m3a <- try(lmer(lambda ~ luc_s*cc_s*LPI_Realm+bm_s+pa*LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m3a, "try-error")){
        m3a <- NA
    }
    
    m3b <- try(lmer(lambda ~ luc_s*cc_s*LPI_Realm+bm_s+pa*LPI_Realm+Managed*LPI_Realm+Utilised*LPI_Realm+(1|spp_idx)+(1|loc_idx),data=df,REML=F))
    if(inherits(m3b, "try-error")){
        m3b <- NA
    }
    
    m_ls <- list(m1a=m1a, m1b=m1b, 
                 m2a=m2a, m2b=m2b, 
                 m3a=m3a, m3b=m3b
                 )
    
    m_ls <- m_ls[!is.na(m_ls)]
    
    rsq <- as.data.frame(do.call(rbind,lapply(m_ls, r.squaredGLMM)))
    rsq$model <- names(m_ls)
    
    aicc <- as.data.frame(do.call(rbind,lapply(m_ls, AICc)))
    colnames(aicc) <- "AICc"
    aicc$model <- rownames(aicc)
    
    metr_df <- merge(rsq, aicc)
    
    coef_df <- do.call(rbind, lapply(m_ls, function(m){

        coef_df <- as.data.frame(summary(m)$coefficients[,1:2])
        colnames(coef_df) <- c("coef_val", "coef_se")
        coef_df$coef_nm <- rownames(coef_df)
        
        if (CI == T){
            confints <- try(confint(m))
            if(inherits(confints, "try-error")){
                confints <- as.data.frame(matrix(NA, nrow = nrow(coef_df), ncol = 2))
                colnames(confints) <- c("lowCI", "highCI")
                confints$coef_nm <- rownames(coef_df)
            }
            else {
                confints <- data.frame(confints)
                colnames(confints) <- c("lowCI", "highCI")
                confints$coef_nm <- rownames(confints)
            }
            coef_df <- merge(coef_df, confints)
        }
        coef_df
    }))
    
    coef_df$model <- sapply(rownames(coef_df), FUN = function(x){strsplit(x, "\\.")[[1]][1]})
    
    if (CI == T){
        coef_df[,7:9] <- ((10^coef_df[c("coef_val", "lowCI", "highCI")]) - 1)*100
    }
    else if (CI == F){
        coef_df[,5] <- ((10^coef_df[c("coef_val")]) - 1)*100
    }
    
    
    # merge 
    comb_df <- merge(coef_df, metr_df, all.x = T)
    
    comb_df$cc_col <- cc_col
    comb_df$luc_col <- luc_col
    
    comb_df$dp_min <- dp_min
    
    return(comb_df)
}

m_fit_wrp <- function(df, cc_cl, luc_col, dp_min, spec){
    print("Main")
    main <- m_fit(df, T, cc_cl, luc_col, dp_min)
    # 20220426: add spec info to results dataframe
    main <- cbind(main, spec)
    main$ecol_sub <- "main"

    return(list(main = main))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main code
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data
dat <- readRDS("../Data/anon_dat.rds")
print("Data loaded.")
# Get rate specs
m_spec <- readRDS("../Data/lag_specs.rds")

# 6282 combinations
# Testing
# j <- 2

for (j in 1:nrow(m_spec)){
    f = paste0("lpi_lag_res_", as.character(j), ".rds")
    fp = paste0("../Results/initial_models/", f)
    
    spec <- m_spec[j,]
    
    if (spec$class == "bird"){
        df <- dat$bird_df
    } 
    if (spec$class == "mammal"){
        df <- dat$mam_df
    } 
    
    # df$cru3.23 <- NA
    df$cru4.04 <- NA
    df$ipsl <- NA
    
    # df$hyde3.1 <- NA
    df$hyde3.2_crng <- NA
    df$hyde3.2_grz <- NA
    df$luh2_norng <- NA
    df$luh2_rng <- NA
    
    env_df <- dat$env_df
    
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
        cc_dat <- cc_calc(cc_strt, cc_end, df$loc_idx[i])
        luc_dat <- luc_calc(luc_strt, luc_end, df$loc_idx[i])
        
        df$cru4.04[i] <- cc_dat$cc2
        df$ipsl[i] <- cc_dat$cc3
        
        df$hyde3.2_crng[i] <- luc_dat$luc2
        df$hyde3.2_grz[i] <- luc_dat$luc3
        df$luh2_norng[i] <- luc_dat$luc4
        df$luh2_rng[i] <- luc_dat$luc5
    }
    
    
    # Fit models
    # Combinations of env data and dpmin
    # Initially, just ipsl, luh2 and 3
    ccluc_comb <- data.frame(expand.grid(cc_col = c("ipsl"),
                                         luc_col = c("luh2_norng", "luh2_rng"),
                                         dp_min = c(3),
                                         stringsAsFactors = F))
    out_ls <- vector("list", nrow(ccluc_comb))
    
    for(i in 1:nrow(ccluc_comb)){
        df1 <- subset(df, DP_to_2014 >= ccluc_comb$dp_min[i] & lambda_rsq > 0)
        
        res <- m_fit_wrp(df1, ccluc_comb$cc_col[i], ccluc_comb$luc_col[i], ccluc_comb$dp_min[i],
                         spec)
        # Bind dfs in res tog
        res <- do.call(rbind.fill, res)
        out_ls[[i]] <- res
    }
    # Bind out 1 and 2
    out_df <- do.call(rbind.fill, out_ls)
    
    # Save
    saveRDS(out_df,
            fp)
}



