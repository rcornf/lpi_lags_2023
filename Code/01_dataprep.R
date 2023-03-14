##%%%%%%%%%%
## R script to prepare population data, land-use data and temperature date for 
## analyses of lagged effects of climate and land use on population abundance trends
##%%%%%%%%%%


rm(list = ls())
graphics.off()

# Packages
# Load packages
library(reshape2)       #~ Data wrangling
library(plyr)           #~ Data wrangling
library(mgcv)           #~ GAMs
library(MuMIn)          #~ Model selection
library(car)            #~ Model selection
library(zoo)            #~ Linear interpolation
library(raster)         #~ Spatial data
library(ncdf4)

## NOTE: This script requires a number of data files that are not provided (due to 
## large file sizes).
## We do however include a redacted version of the LPD ubset required for all 
## subsequent analysis.

## Required data specifications:

## LPD
## LUH2 - https://luh.umd.edu/data.shtml, v2h Release, states.nc
## IPSL - https://data.isimip.org/datasets/ffade380-1650-4179-b522-b185712ed9fd/
## HYDE 3.2 - https://geo.public.data.uu.nl/vault-hyde-data/HYDE%203.2[1648738557]
## CRU TS 4.04 - https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/tmp/
## Amniote database - https://figshare.com/collections/An_amniote_life-history_database_to_perform_comparative_analyses_with_birds_mammals_and_reptiles/3308127
## Elton traits - https://figshare.com/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
## Generation length data - 
##    - Birds: https://conbio.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fcobi.13486&file=cobi13486-sup-0004-TableS4.xlsx
##    - Mammals: https://datadryad.org/stash/dataset/doi:10.5061/dryad.gd0m3

## A user supplied token is also required for access to the IUCN API.


##~~~~~~~~~~~~~~~~~
## Prepare/filter LPD data
##~~~~~~~~~~~~~~~~~


# Load lpd
lpd <- read.csv("../Data/lpd.csv",
                stringsAsFactors = F, na.strings = c("NULL", ""))

# Only consider 1950-2014 

# Subset to ter/fw, mammal/bird, specific location,
# PA!="Both" or "Unknown"
# Not marine, not excluded
lpd_sub <- subset(lpd, Class %in% c("Mammalia", "Aves") & Excluded != 1 &
                      Specific_location == 1 & Protected_status != "Both" &
                      Protected_status != "Unknown" & System != "Marine")


## Get DP, TS_strt, TS_end, TS_length, non-zero abund
for (i in 1:nrow(lpd_sub)){
    tmp <- lpd_sub[i,c(which(colnames(lpd_sub) == "X1950"):
                         which(colnames(lpd_sub) == "X2014"))]
    tmp_lpi_df <- data.frame("year" = seq(1950,2014),
                             "abund" = unname(t(tmp)))
    
    # dp
    lpd_sub$DP_to_2014[i] <- sum(!is.na(tmp_lpi_df$abund))
    
    tmp_lpi_df <- na.omit(tmp_lpi_df)
    # ts strt/end
    lpd_sub$TS_strt[i] <- min(tmp_lpi_df$year)
    lpd_sub$TS_end[i] <- max(tmp_lpi_df$year)
    # non-zero
    lpd_sub$non_zero_abund[i] <- ifelse(sum(tmp_lpi_df$abund == 0) == nrow(tmp_lpi_df), 0, 1)
}


lpd_sub$TS_length_to_2014 <- lpd_sub$TS_end-lpd_sub$TS_strt+1


lpd_sub <- subset(lpd_sub, DP_to_2014 >= 2 & TS_length_to_2014 >= 5 & non_zero_abund == 1)
# 3392 pops
rm(i, tmp, tmp_lpi_df)

lpd_sub$loc_id <- paste(lpd_sub$Latitude, lpd_sub$Longitude, sep = "_")
sum(lpd_sub$Replicate)
# 924


# Combine sytem and realm into a single col
# Unite terrestrial and frewshwater realms
lpd_sub$Realm <- apply(cbind(lpd_sub$T_realm, lpd_sub$FW_realm), 1,
                   function(x) paste(x[!is.na(x)], collapse = ""))

# Combine realms to LPI-style values
lpd_sub$LPI_Realm <- lpd_sub$Realm
# combine Austrlasia and Indo-Malayan into Indo-Pacific
lpd_sub$LPI_Realm[lpd_sub$Realm=="Australasia"] <- "Indo-Pacific"
lpd_sub$LPI_Realm[lpd_sub$Realm=="Indo-Malayan"] <- "Indo-Pacific"
lpd_sub$LPI_Realm[lpd_sub$Realm=="Oceania"] <- "Indo-Pacific"


lpd_sub$Sys_realm <- paste(lpd_sub$System, lpd_sub$Realm, sep = "_")

lpd_repl <- subset(lpd_sub, Replicate == 1)
lpd_non_repl <- subset(lpd_sub, Replicate == 0)
# subset(lpd_sub, Binomial == "Ceratotherium_simum")[,c("ID", "Replicate", "System", "T_realm")]


# Drop replicates
# But check if another pop with same binom and sys_realm in data or not...
lpd_repl$realm_repl <- 0
lpd_repl$ts_repl <- 0
lpd_repl$loc_repl <- 0
lpd_repl$loc_ts_repl <- 0

for (i in 1:nrow(lpd_repl)){
    # Get binom, sys_realm
    tmp <- subset(lpd_non_repl, Binomial == lpd_repl$Binomial[i] &
                      Sys_realm == lpd_repl$Sys_realm[i])

    if (nrow(tmp)>0){
        lpd_repl$realm_repl[i] <- 1
        # look at time-series overlap too..
        # if end non-repl >= strt repl+1, or strt non_repl <= end repl-1: overlap >= 2 data yrs
        if (sum(tmp$TS_end >= lpd_repl$TS_strt[i] + 1 | tmp$TS_strt <= lpd_repl$TS_end[i] - 1)>0){
            lpd_repl$ts_repl[i] <- 1
        }
    }
    tmp1 <- subset(lpd_non_repl, Binomial == lpd_repl$Binomial[i] &
                       Sys_realm == lpd_repl$Sys_realm[i] & 
                       loc_id == lpd_repl$loc_id[i])
    # print(tmp1)
    if (nrow(tmp1)>0){
        lpd_repl$loc_repl[i] <- 1
        if (sum(tmp1$TS_end >= lpd_repl$TS_strt[i] + 1 | tmp1$TS_strt <= lpd_repl$TS_end[i] - 1)>0){
            lpd_repl$loc_ts_repl[i] <- 1
        }
    }  
}

rm(i, tmp, tmp1)

# even if acounting for overlap, still same number of replicates to drop...
sum(lpd_repl$realm_repl) # 605
sum(lpd_repl$ts_repl) # 605
sum(lpd_repl$loc_repl) # 35
sum(lpd_repl$loc_ts_repl) # 35

# only 35 pops are replicates at the location level....

# rm replicates at loc level
lpd_repl_sub <- subset(lpd_repl, loc_ts_repl == 0) 


lpd_non_repl$realm_repl <- 0
lpd_non_repl$ts_repl <- 0
lpd_non_repl$loc_repl <- 0
lpd_non_repl$loc_ts_repl <- 0


lpd_sub1 <- rbind(lpd_non_repl, 
                  lpd_repl_sub)

# Note: depending on replicate rm approcah, some more pops may need to be dropped...

nrow(subset(lpd_sub1, Class == "Aves")) # 1979
nrow(subset(lpd_sub1, Class == "Mammalia")) # 1378


rm(lpd_repl, lpd_non_repl, lpd_sub, lpd_repl_sub)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract body mass, generation length and trophic info for each species
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bird_lpd <- subset(lpd_sub1, Class == "Aves")
mam_lpd <- subset(lpd_sub1, Class == "Mammalia")

# Body mass dataset
# amniote_df <- read.csv("../Data/ECOL_96_269/Data_Files/Amniote_Database_Aug_2015.csv")
amniote_df <- read.csv("../Data/Amniote_Database_Aug_2015.csv")
amniote_df$Binom <- paste(as.character(amniote_df$genus), 
                          as.character(amniote_df$species),
                          sep = "_")

# Gen length data
bird_gen_df <- readxl::read_excel("../Data/bird_gen_length_Bird2020.xlsx")
mam_gen_df <- readxl::read_excel("../Data/mam_gen_length_Pacifici2013.xls")

# Elton traits
bird_elton <- read.delim("../Data/BirdFuncDat.txt", sep = "\t", 
                         stringsAsFactors = F,
                         header = T)
bird_elton$Herb_pcnt <- bird_elton$Diet.Fruit+bird_elton$Diet.Nect+bird_elton$Diet.Seed+bird_elton$Diet.PlantO 
bird_elton$Carn_pcnt <- bird_elton$Diet.Inv+bird_elton$Diet.Vend+bird_elton$Diet.Vect+bird_elton$Diet.Vfish+bird_elton$Diet.Vunk+bird_elton$Diet.Scav

mam_elton <- read.delim("../Data/MamFuncDat.txt", sep = "\t", 
                        stringsAsFactors = F,
                        header = T)
mam_elton$Herb_pcnt <- mam_elton$Diet.Fruit+mam_elton$Diet.Nect+mam_elton$Diet.Seed+mam_elton$Diet.PlantO 
mam_elton$Carn_pcnt <- mam_elton$Diet.Inv+mam_elton$Diet.Vend+mam_elton$Diet.Vect+mam_elton$Diet.Vfish+mam_elton$Diet.Vunk+mam_elton$Diet.Scav

##~~~~~
## Names in LPD do not always match those in the other databases used for trait data
## Use IUCN to find binomials/synonyms and match automatically where possible
## Otherwise match manually 
##~~~~~

# mam_bm_binoms 
setdiff(mam_lpd$Binomial, amniote_df$Binom) # 3 diffs
# mam_gl_binoms 
setdiff(gsub("_", " ", mam_lpd$Binomial), mam_gen_df$Scientific_name) # 22 diffs
# mam_tr_binoms 
setdiff(gsub("_", " ", mam_lpd$Binomial), mam_elton$Scientific) # 4 diffs
# bird_bm_binoms 
setdiff(bird_lpd$Binomial, amniote_df$Binom) # 49 diffs
# bird_gl_binoms 
setdiff(gsub("_", " ", bird_lpd$Binomial), bird_gen_df$`Scientific name`) # no diffs
# bird_tr_binoms 
setdiff(gsub("_", " ", bird_lpd$Binomial), bird_elton$Scientific) # 87 diffs

# Initially format for merging
mam_lpd$Binom_bm <- mam_lpd$Binomial
mam_lpd$Binom_gl <- gsub("_", " ", mam_lpd$Binomial)
mam_lpd$Binom_tr <- gsub("_", " ", mam_lpd$Binomial)

bird_lpd$Binom_bm <- bird_lpd$Binomial
bird_lpd$Binom_gl <- gsub("_", " ", bird_lpd$Binomial)
bird_lpd$Binom_tr <- gsub("_", " ", bird_lpd$Binomial)


library(httr)
tok = "<usr_defined>"

# Set up function to call IUCN API with missing names, and add useful synonyms to df
iucn_synonyms_check <- function(lpd_binoms, trait_binoms){
    match_binoms <- rep(NA, length(lpd_binoms))
    for (i in 1:length(lpd_binoms)){
        
        r <- GET(paste0("http://apiv3.iucnredlist.org/api/v3/species/synonym/",
                        lpd_binoms[i],
                        "?token=",
                        tok))
        
        j <- jsonlite::fromJSON(content(r, "text"))
        
        # strip white space
        # unique synonyms
        syns <- unique(trimws(j$result$synonym))
        print(syns)
        # if in data binom, keep
        # should only be 1 match!
        if (length(syns[syns %in% trait_binoms])==1){
            match_binoms[i] <- syns[syns %in% trait_binoms]
        }
        if (length(syns[syns %in% trait_binoms])>1){
            print(syns[syns %in% trait_binoms])
        }
    }
    return(data.frame(lpd_binom = lpd_binoms,
                      match_binom = match_binoms,
                      stringsAsFactors = F))
}

# mam_bm_binoms 
mam_bm_iucn_match <- iucn_synonyms_check(gsub("_", " ", setdiff(mam_lpd$Binom_bm, 
                                                                amniote_df$Binom)),
                                         gsub("_", " ", amniote_df$Binom))
# No NAs, add names in, replace " ", with "_"
for(i in 1:nrow(mam_bm_iucn_match)){
    mam_lpd$Binom_bm[mam_lpd$Binom_bm == gsub(" ", "_", mam_bm_iucn_match$lpd_binom[i])] <- gsub(" ", "_", mam_bm_iucn_match$match_binom[i])
}
setdiff(mam_lpd$Binom_bm, amniote_df$Binom)
# Checked

# mam_gl_binoms 
mam_gl_iucn_match <- iucn_synonyms_check(setdiff(mam_lpd$Binom_gl, 
                                                 mam_gen_df$Scientific_name),
                                         mam_gen_df$Scientific_name)
subset(mam_gl_iucn_match, is.na(match_binom))
# first, go through iucn matches
for(i in 1:nrow(mam_gl_iucn_match)){
    if (!is.na(mam_gl_iucn_match$match_binom[i])){
        mam_lpd$Binom_gl[mam_lpd$Binom_gl == mam_gl_iucn_match$lpd_binom[i]] <- mam_gl_iucn_match$match_binom[i]
    }
}

mam_gl_man_match <- list(
    "Damaliscus korrigum"       = "Damaliscus lunatus",
    "Equus burchellii"          = "Equus quagga",
    "Piliocolobus tephrosceles" = "Procolobus badius",
    "Taurotragus derbianus"     = "Tragelaphus derbianus",
    "Taurotragus oryx"          = "Tragelaphus oryx",
    "Bos frontalis"             = "Bos gaurus",
    "Bubalus bubalis"           = "Bubalus arnee",
    "Hemitragus hylocrius"      = "Nilgiritragus hylocrius",
    # "Ovis aries"               
    "Capra hircus"              = "Capra aegagrus",
    "Saiga borealis"            = "Saiga tatarica",
    "Galerella sanguinea"       = "Herpestes sanguineus",
    "Loxodonta cyclotis"        = "Loxodonta africana",
    "Bos grunniens"             = "Bos mutus",
    "Alcelaphus caama"          = "Alcelaphus buselaphus",
    "Hydrochoeris hydrochaeris" = "Hydrochoerus hydrochaeris",
    "Platanista minor"          = "Platanista gangetica",
    "Uncia uncia"               = "Panthera uncia",
    "Piliocolobus gordonorum"   = "Procolobus gordonorum")
                   
# then go through manual matches
for(i in 1:length(mam_gl_man_match)){
    mam_lpd$Binom_gl[mam_lpd$Binom_gl == names(mam_gl_man_match)[i]] <- mam_gl_man_match[[i]]
}
setdiff(mam_lpd$Binom_gl, mam_gen_df$Scientific_name)
# ovis aries

# mam_tr_binoms 
mam_tr_iucn_match <- iucn_synonyms_check(setdiff(mam_lpd$Binom_tr, 
                                                 mam_elton$Scientific),
                                         mam_elton$Scientific)
# Still an NA
mam_tr_iucn_match
mam_lpd$Binom_tr[mam_lpd$Binom_tr == "Hydrochoeris hydrochaeris"] <- "Hydrochoerus hydrochaeris"

for(i in 1:nrow(mam_tr_iucn_match)){
    if (!is.na(mam_tr_iucn_match$match_binom[i])){
        mam_lpd$Binom_tr[mam_lpd$Binom_tr == mam_tr_iucn_match$lpd_binom[i]] <- mam_tr_iucn_match$match_binom[i]
    }
}
setdiff(mam_lpd$Binom_tr, mam_elton$Scientific)



# bird_bm_binoms 
bird_bm_iucn_match <- iucn_synonyms_check(gsub("_", " ", setdiff(bird_lpd$Binomial, 
                                                                 amniote_df$Binom)),
                                          gsub("_", " ", amniote_df$Binom))
subset(bird_bm_iucn_match, is.na(match_binom))
# first, go through iucn matches
for(i in 1:nrow(bird_bm_iucn_match)){
    if (!is.na(bird_bm_iucn_match$match_binom[i])){
        bird_lpd$Binom_bm[bird_lpd$Binom_bm == gsub(" ", "_", bird_bm_iucn_match$lpd_binom[i])] <- gsub(" ", "_", bird_bm_iucn_match$match_binom[i])
    }
}
# Manual fixes where possible
bird_bm_man_match <- list(#"Amazona_arausiaca" ,
    "Antigone_canadensis"      = "Grus_canadensis",
    "Antigone_vipio"           = "Grus_vipio",
    "Ardea_intermedia"         = "Egretta_intermedia",
    "Ardea_plumifera"          = "Egretta_intermedia",
    "Basileuterus_hypoleucus"  = "Basileuterus_culicivorus",
    "Bonasa_bonasia"           = "Tetrastes_bonasia",
    "Ceratopipra_rubrocapilla" = "Dixiphia_rubrocapilla",
    "Corvus_monedula"          = "Coloeus_monedula",
    "Cyanoloxia_cyanoides"     = "Cyanocompsa_cyanoides",
    "Hylatomus_lineatus"       = "Dryocopus_lineatus",
    "Larus_audouinii"          = "Ichthyaetus_audouinii",
    "Larus_cirrocephalus"      = "Chroicocephalus_cirrocephalus",
    "Larus_genei"              = "Chroicocephalus_genei",
    "Larus_ichthyaetus"        = "Ichthyaetus_ichthyaetus",
    "Larus_maculipennis"       = "Chroicocephalus_maculipennis",
    "Larus_melanocephalus"     = "Ichthyaetus_melanocephalus",
    "Microcarbo_pygmaeus"      = "Microcarbo_pygmeus",
    # "Petroica_longipes"        = "Petroica_australis", #supsp?
    "Pygochelidon_cyanoleuca"  = "Notiochelidon_cyanoleuca",
    # "Saltator_fuliginosus"
    # "Todirostrum_poliocephalum"
    "Xenops_rutilus"           = "Xenops_rutilans")
    # "Xolmis_velatus" 

# then go through manual matches
for(i in 1:length(bird_bm_man_match)){
    bird_lpd$Binom_bm[bird_lpd$Binom_bm == names(bird_bm_man_match)[i]] <- bird_bm_man_match[[i]]
}
setdiff(bird_lpd$Binom_bm, amniote_df$Binom)

# bird_tr_binoms 
bird_tr_iucn_match <- iucn_synonyms_check(setdiff(bird_lpd$Binom_tr, 
                                                 bird_elton$Scientific),
                                         bird_elton$Scientific)
for(i in 1:nrow(bird_tr_iucn_match)){
    if (!is.na(bird_tr_iucn_match$match_binom[i])){
        bird_lpd$Binom_tr[bird_lpd$Binom_tr == bird_tr_iucn_match$lpd_binom[i]] <- bird_tr_iucn_match$match_binom[i]
    }
}
subset(bird_tr_iucn_match, is.na(match_binom))
# Manual checking - Using birdlife and BirdFuncDat
bird_tr_man_match <- list(
    "Aramides cajaneus"      = "Aramides cajanea",
    "Cecropis daurica"       = "Hirundo daurica",
    "Chlamydotis macqueenii" = "Chlamydotis undulata",
    "Circus hudsonius"       = "Circus cyaneus",
    "Hylatomus lineatus"     = "Dryocopus lineatus",
    "Ardea intermedia"       = "Mesophoyx intermedia",
    "Ardea plumifera"        = "Mesophoyx intermedia",
    "Gallinula galeata"      = "Gallinula chloropus",
    "Gelochelidon nilotica"  = "Sterna nilotica",
    "Antigone canadensis"    = "Grus canadensis",
    "Antigone vipio"         = "Grus vipio",
    "Microcarbo pygmaeus"    = "Phalacrocorax pygmeus",
    "Poecile montanus"       = "Parus montanus",
    "Turdus atrogularis"     = "Turdus ruficollis",
    "Xenops rutilus"         = "Xenops rutilans",
    "Cyanoloxia cyanoides"   = "Cyanocompsa cyanoides"
    # "Petroica longipes"      = "Petroica_australis" diff species
)
for(i in 1:length(bird_tr_man_match)){
    bird_lpd$Binom_tr[bird_lpd$Binom_tr == names(bird_tr_man_match)[i]] <- bird_tr_man_match[[i]]
}
setdiff(bird_lpd$Binom_tr, bird_elton$Scientific)

###
# Merge population data with bm, gen length and trophic info
###
mam_lpd1 <- merge(mam_lpd, amniote_df[,c("Binom", "adult_body_mass_g")], 
                 by.x = "Binom_bm", by.y = "Binom", all.x = T)

mam_lpd1 <- merge(mam_lpd1, mam_gen_df[,c("Scientific_name", "GenerationLength_d")], 
                 by.x = "Binom_gl", by.y = "Scientific_name", all.x = T)
mam_lpd1$GenLength <- mam_lpd1$GenerationLength_d/365

mam_lpd1 <- merge(mam_lpd1, mam_elton[,c("Scientific", "Herb_pcnt", "Carn_pcnt")], 
                 by.x = "Binom_tr", by.y = "Scientific", all.x = T)

bird_lpd1 <- merge(bird_lpd, amniote_df[,c("Binom", "adult_body_mass_g")], 
                 by.x = "Binom_bm", by.y = "Binom", all.x = T)

bird_lpd1 <- merge(bird_lpd1, bird_gen_df[,c("Scientific name", "GenLength")], 
                 by.x = "Binom_gl", by.y = "Scientific name", all.x = T)

bird_lpd1 <- merge(bird_lpd1, bird_elton[,c("Scientific", "Herb_pcnt", "Carn_pcnt")], 
                 by.x = "Binom_tr", by.y = "Scientific", all.x = T)


# Check for NAs
range(mam_lpd1$adult_body_mass_g)
range(bird_lpd1$adult_body_mass_g, na.rm = T)
# Convert <0 to NA

mam_lpd1$adult_body_mass_g[mam_lpd1$adult_body_mass_g < 0] <- NA
bird_lpd1$adult_body_mass_g[bird_lpd1$adult_body_mass_g < 0] <- NA

range(mam_lpd1$GenLength, na.rm = T)
range(bird_lpd1$GenLength, na.rm = T)
sum(is.na(mam_lpd1$GenLength))
sum(is.na(bird_lpd1$GenLength))


# Drop pops with na for bm or gen length
mam_lpd1 <- subset(mam_lpd1, !is.na(adult_body_mass_g) & !is.na(GenLength))
bird_lpd1 <- subset(bird_lpd1, !is.na(adult_body_mass_g) & !is.na(GenLength))


# Tidy
rm(amniote_df, bird_bm_iucn_match, bird_bm_man_match, bird_elton, bird_gen_df,
   bird_tr_iucn_match, bird_tr_man_match, 
   mam_bm_iucn_match, mam_elton, mam_gen_df, mam_gl_iucn_match, mam_gl_man_match,
   mam_tr_iucn_match)

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Hyde/LUH2 and CRU/IPSL data from 1901-2014 for each location
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loc_df <- unique(rbind(mam_lpd1[,c("loc_id", "Longitude", "Latitude")],
                       bird_lpd1[,c("loc_id", "Longitude", "Latitude")]))

points <- SpatialPoints(loc_df[, c("Longitude", "Latitude")])

###
# Hyde3.2
###

## Preprocess hyde data to proportional area values

# Function from Tim Newbold's github to calculate area of cells
DegreeCellAreaKM <- function(lat, height, width) {
    # Returns the area in km squared of a grid cell in degrees of arc
    # lat - the latitudinal centre of the cell
    # height, width - the size of the grid cell in degrees
    
    radians <- function(theta) theta*pi/180.0
    
    # Convert the latitude into radians
    lat.rad <- radians(lat)
    
    # The equatorial and polar radii of the Earth in km
    eq.radius <-  6378137
    pol.radius <- 6356752.3142
    
    # Calculate cell area
    angular.eccentricity <- acos(radians(pol.radius/eq.radius))
    ecc.sq <- sin(radians(angular.eccentricity))^2
    flattening <- 1-cos(radians(angular.eccentricity))
    temp.val <- (eq.radius*cos(lat.rad))^2+(pol.radius*sin(lat.rad))^2
    m.phi <- ((eq.radius*pol.radius)^2)/(temp.val^1.5)
    n.phi <- (eq.radius^2)/sqrt(temp.val)
    lat.length <- pi/180*m.phi/1000
    long.length <- pi/180*cos(lat.rad)*n.phi/1000
    return (lat.length*height*long.length*width)
}

hyde_f <- list.files(path = "../Data/HYDE3.2",
                     pattern = c("grazing|pasture|cropland|rangeland"),
                     recursive = T,
                     full.names = T)

hyde1 <- raster(hyde_f[1])
cellareas <- DegreeCellAreaKM(lat=coordinates(hyde1)[,2],height=res(hyde1)[2],width=res(hyde1)[1])

# Move files into single directory...
for (f in hyde_f){
    rast <- raster(f)
    prop<-signif(values(rast),digits=4)/cellareas
    matr_prop<-matrix(prop, nrow=nrow(rast), ncol=ncol(rast), byrow=TRUE)
    rast_prop<-raster(matr_prop, xmn=rast@extent[1], xmx=rast@extent[2], 
                      ymn=rast@extent[3], ymx=rast@extent[4])
    writeRaster(rast_prop, paste0("../Data/HYDE3.2/prop_area/", names(rast),"_raster.tif"),overwrite=TRUE)
}

rm(rast, rast_prop, prop, matr_prop, hyde1, hyde_f, f, cellareas)

crop3.2_s <- stack(list.files(path = "../Data/HYDE3.2/prop_area",
                              pattern = c("cropland"),
                              recursive = F,
                              full.names = T))

crng3.2_s <- stack(list.files(path = "../Data/HYDE3.2/prop_area",
                              pattern = c("conv_rangeland"),
                              recursive = F,
                              full.names = T))

pstr3.2_s <- stack(list.files(path = "../Data/HYDE3.2/prop_area",
                              pattern = c("pasture"),
                              recursive = F,
                              full.names = T))

rng3.2_f <- list.files(path = "../Data/HYDE3.2/prop_area",
                       pattern = c("rangeland"),
                       recursive = F,
                       full.names = T)
rng3.2_s <- stack(rng3.2_f[!grepl("conv", rng3.2_f)])

graz3.2_s <- stack(list.files(path = "../Data/HYDE3.2/prop_area",
                              pattern = c("grazing"),
                              recursive = F,
                              full.names = T))

hyde3.2_grid_ext <- data.frame("loc_id" = loc_df$loc_id,
                               "xmin" = colFromX(crop3.2_s[[1]], points@coords[, 1]) - 1,
                               "xmax" = colFromX(crop3.2_s[[1]], points@coords[, 1]) + 1,
                               "ymin" = rowFromY(crop3.2_s[[1]], points@coords[, 2]) - 1,
                               "ymax" = rowFromY(crop3.2_s[[1]], points@coords[, 2]) + 1)

library(doParallel)
cores<- 6
cl <- makeCluster(cores, output="") #output should make it spit errors
registerDoParallel(cl)

hyde3.2_df <- foreach(idx=1:nrow(loc_df), .packages = c("raster", "zoo"), 
                      .combine = rbind) %dopar% {
 
                          crop3.2_cells <- crop(crop3.2_s, extent(crop3.2_s,
                                                                  hyde3.2_grid_ext[idx,"ymin"], hyde3.2_grid_ext[idx,"ymax"],
                                                                  hyde3.2_grid_ext[idx,"xmin"], hyde3.2_grid_ext[idx,"xmax"]))
                          crng3.2_cells <- crop(crng3.2_s, extent(crng3.2_s,
                                                                  hyde3.2_grid_ext[idx,"ymin"], hyde3.2_grid_ext[idx,"ymax"],
                                                                  hyde3.2_grid_ext[idx,"xmin"], hyde3.2_grid_ext[idx,"xmax"]))
                          pstr3.2_cells <- crop(pstr3.2_s, extent(pstr3.2_s,
                                                                  hyde3.2_grid_ext[idx,"ymin"], hyde3.2_grid_ext[idx,"ymax"],
                                                                  hyde3.2_grid_ext[idx,"xmin"], hyde3.2_grid_ext[idx,"xmax"]))
                          rng3.2_cells <- crop(rng3.2_s, extent(rng3.2_s,
                                                                hyde3.2_grid_ext[idx,"ymin"], hyde3.2_grid_ext[idx,"ymax"],
                                                                hyde3.2_grid_ext[idx,"xmin"], hyde3.2_grid_ext[idx,"xmax"]))
                          graz3.2_cells <- crop(graz3.2_s, extent(graz3.2_s,
                                                                  hyde3.2_grid_ext[idx,"ymin"], hyde3.2_grid_ext[idx,"ymax"],
                                                                  hyde3.2_grid_ext[idx,"xmin"], hyde3.2_grid_ext[idx,"xmax"]))
                          
                          crop3.2_matr <- as.matrix(crop3.2_cells)
                          crng3.2_matr <- as.matrix(crng3.2_cells)
                          pstr3.2_matr <- as.matrix(pstr3.2_cells)
                          rng3.2_matr  <- as.matrix(rng3.2_cells)
                          graz3.2_matr <- as.matrix(graz3.2_cells)
                          
                          # Caclulate average crop/grass coverage per decade
                          hyde_av_df <- data.frame("Year" = c(1900,1910,1920,1930,1940,1950,1960,1970,
                                                              1980,1990,2000:2017),
                                                   "mean_crop3.2_fs" = colMeans(crop3.2_matr,na.rm = F),
                                                   "mean_crng3.2_fs" = colMeans(crng3.2_matr,na.rm = F),
                                                   "mean_pstr3.2_fs" = colMeans(pstr3.2_matr,na.rm = F),
                                                   "mean_rng3.2_fs" = colMeans(rng3.2_matr ,na.rm = F),
                                                   "mean_graz3.2_fs" = colMeans(graz3.2_matr,na.rm = F)
                          )
                         
                          # Set any hyde crop/grass values <0 >1 to 0, 1 respectively
                          hyde_av_df[,c(2:6)][hyde_av_df[,c(2:6)]>1] <- 1
                          hyde_av_df[,c(2:6)][hyde_av_df[,c(2:6)]<0] <- 0
                          
                          # Interpolate to get years of interest
                          hyde_interp_df <- data.frame("Year" = 1900:2017)
                          hyde_interp_df <- merge(hyde_interp_df, hyde_av_df[,c(1:6)], by = "Year", all.x = T)
                          
                          # for each column, if  >=2 !is.na, interpolate between first and end yr..
                          hyde_interp_df[,c(2:6)] <- apply(hyde_interp_df[,c(2:6)],
                                                           2,
                                                           function(x){
                                                               dat_idx <- !is.na(x)
                                                               n_decs <- sum(dat_idx)
                                                               strt_end <- range(which(dat_idx))
                                                               
                                                               if (n_decs>2){
                                                                   x[strt_end[1]:strt_end[2]] <- na.approx(x[strt_end[1]:strt_end[2]])
                                                               }
                                                               else {
                                                                   x <- rep(NA, length(x))
                                                               }
                                                               return(x)
                                                           })
                          
                          hyde_interp_df$loc_id <- loc_df$loc_id[idx]
                          hyde_interp_df
                      }

stopImplicitCluster()
saveRDS(hyde3.2_df,
        "../Data/hyde3.2_1900_2017.rds")

# hyde3.2_df <- readRDS("../Data/hyde3.2_1900_2017_20210217.rds")
rm(hyde3.2_grid_ext, cl, crop3.2_s, crng3.2_s, graz3.2_s, pstr3.2_s, rng3.2_s)
rm(rng3.2_f)

###
# LUH2
###
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

names(luh2_1900_2015)

# Extract

luh2_1900_2015_df <- data.frame(extract(luh2_1900_2015, loc_df[,2:3], ncol=2))
head(luh2_1900_2015_df)

rownames(luh2_1900_2015_df) <- loc_df$loc_id
colnames(luh2_1900_2015_df)

luh2_1900_2015_df <- t(luh2_1900_2015_df)
luh2_1900_2015_yr_var <- data.frame(do.call(rbind, strsplit(rownames(luh2_1900_2015_df), "_")))
colnames(luh2_1900_2015_yr_var) <- c("Year", "LU")
luh2_1900_2015_yr_var$Year <- as.numeric(gsub("X", "", luh2_1900_2015_yr_var$Year))

luh2_1900_2015_long_df <- reshape2::melt(cbind(luh2_1900_2015_df, luh2_1900_2015_yr_var), 
                                         id.vars = c("Year", "LU"),
                                         variable.name = "loc_id",
                                         value.name = "prop_cover")

# make luh2 wide? - col per lu type
luh2_1900_2015_wide_df <- tidyr::spread(luh2_1900_2015_long_df, LU, prop_cover)

head(luh2_1900_2015_wide_df)

rm(luh2, luh2_1900_2015, luh2_1900_2015_df, luh2_1900_2015_yr_var, luh2_1900_2015_long_df)
# Save luh2
saveRDS(luh2_1900_2015_wide_df,
        "../Data/luh2_1900_2015.rds")


###
# CRU4.04
###
cru_t <- brick("../Data/cru_ts4.04.1901.2019.tmp.dat.nc", var = "tmp")

cru_t <- cru_t[[min(grep(1901, names(cru_t))):max(grep(2017, names(cru_t)))]]
cru_t <- data.frame(extract(cru_t, loc_df[,2:3], ncol=2)) 
rownames(cru_t) <- loc_df$loc_id

cru_t <- data.frame(sapply(1901:2017, 
                           FUN = function(x){rowMeans(cru_t[,grep(x, names(cru_t))])}))
rownames(cru_t)
head(cru_t)
colnames(cru_t) <- 1901:2017
cru_t$loc_id <- rownames(cru_t)
cru_t <- reshape2::melt(cru_t, id.vars = "loc_id", 
                        variable.name = "Year",
                        value.name = "cru4.04")
cru_t$Year <- as.numeric(as.character(cru_t$Year))
head(cru_t)

saveRDS(cru_t[,c("loc_id", "Year", "cru4.04")], 
        "../Data/cru4.04_1901_2017.rds")

###
# IPSL
###
ipsl_f_ls <- list.files(path = "../Data/IPSL",
                        full.names = T)

ipsl_ls <- purrr::map(ipsl_f_ls,
                      function(x, loc){
                          tmp_brk <- brick(x)
                          
                          # Option 1
                          # extract, mean
                          tmp_t <- data.frame(extract(tmp_brk, loc[,2:3], ncol=2)) 
                          
                          rownames(tmp_t) <- loc$loc_id
                          # need to know unique years...
                          tmp_yrs <- sort(as.numeric(gsub("X", "", 
                                                          unique(sapply(strsplit(names(tmp_brk), "\\."), getElement, 1)))))
                          
                          tmp_t <- data.frame(sapply(tmp_yrs, 
                                                     FUN = function(x){rowMeans(tmp_t[,grep(x, names(tmp_t))])}))
                          
                          colnames(tmp_t) <- tmp_yrs
                          tmp_t$loc_id <- rownames(tmp_t)
                          tmp_t <- reshape2::melt(tmp_t, id.vars = "loc_id", 
                                                  variable.name = "Year",
                                                  value.name = "av_t_K")
                          tmp_t$Year <- as.numeric(as.character(tmp_t$Year))
                          tmp_t$ipsl <- tmp_t$av_t_K - 273.15
                          tmp_t
                      },
                      loc_df)

ipsl_df <- do.call(rbind, ipsl_ls)

rm(ipsl_f_ls, ipsl_ls)
saveRDS(ipsl_df, 
        "../Data/ipsl_1901_2014.rds")
# ipsl_df <- readRDS("../Data/ipsl_1901_2014.rds")

# Sum hyde and luh2 cols into the alternative types
hyde3.2_df$hyde3.2_crng <- base::rowSums(hyde3.2_df[,c("mean_crop3.2_fs",
                                                             "mean_crng3.2_fs",
                                                             "mean_pstr3.2_fs")])
range(hyde3.2_df$hyde3.2_crng, na.rm = T)

hyde3.2_df$hyde3.2_grz <- base::rowSums(hyde3.2_df[,c("mean_crop3.2_fs",
                                                            "mean_graz3.2_fs")])
range(hyde3.2_df$hyde3.2_grz, na.rm = T) # some >1
hyde3.2_df$hyde3.2_grz[hyde3.2_df$hyde3.2_grz > 1] <- 1
 
luh2_1900_2015_wide_df$luh2_norng <- base::rowSums(luh2_1900_2015_wide_df[,c("c3ann", "c3nfx", "c3per",
                                                     "c4ann", "c4per", "pastr")])
range(luh2_1900_2015_wide_df$luh2_norng, na.rm = T)
 
luh2_1900_2015_wide_df$luh2_rng <- base::rowSums(luh2_1900_2015_wide_df[,c("c3ann", "c3nfx", "c3per",
                                                   "c4ann", "c4per", "pastr",
                                                   "range")])
range(luh2_1900_2015_wide_df$luh2_rng, na.rm = T)

env_df <- plyr::join_all(list(
                              hyde3.2_df[,c("Year", "loc_id", "hyde3.2_crng", "hyde3.2_grz")],
                              luh2_1900_2015_wide_df[,c("Year", "loc_id", "luh2_norng", "luh2_rng")],
                              cru_t[,c("Year", "loc_id", "cru4.04")],
                              ipsl_df[,c("Year", "loc_id", "ipsl")])
)

range(env_df$Year)

# subset to 1901-2014?
env_df_sub <- subset(env_df, Year %in% c(1901:2014))

rm(cru_t, hyde3.2_df, ipsl_df, luh2_1900_2015_wide_df, points)


######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Final collation/organisation
######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset to usable pops
# calc lambdas
# save rds


# General checks
min(mam_lpd1$DP_to_2014)
min(bird_lpd1$DP_to_2014)

min(mam_lpd1$TS_length_to_2014)
min(bird_lpd1$TS_length_to_2014)


nrow(mam_lpd1) # 1369
length(unique(mam_lpd1$Binomial)) # 358
length(unique(mam_lpd1$loc_id)) # 624

nrow(bird_lpd1) # 1960
length(unique(bird_lpd1$Binomial)) # 712
length(unique(bird_lpd1$loc_id)) # 487


# Check for NAs in env time series
env_cov_check <- function(pop_df, env_df){
    pop_df$hyde_crng_nas <- NA
    pop_df$hyde_grz_nas <- NA
    pop_df$luh2_norng_nas <- NA
    pop_df$luh2_rng_nas <- NA
    
    pop_df$cru_nas <- NA
    pop_df$ipsl_nas <- NA
    
    for (i in 1:nrow(pop_df)){
        tmp_loc <- pop_df$loc_id[i]
        tmp_env <- subset(env_df, loc_id == tmp_loc & Year>=1901 & Year <= 2014)
        
        pop_df$hyde_crng_nas[i] <- sum(is.na(tmp_env$hyde3.2_crng))
        pop_df$hyde_grz_nas[i] <- sum(is.na(tmp_env$hyde3.2_grz)) 
        pop_df$luh2_norng_nas[i] <- sum(is.na(tmp_env$luh2_norng))
        pop_df$luh2_rng_nas[i] <- sum(is.na(tmp_env$luh2_rng))
        
        pop_df$cru_nas[i] <- sum(is.na(tmp_env$cru4.04))
        pop_df$ipsl_nas[i] <- sum(is.na(tmp_env$ipsl))
    }
    # row sums... of nas
    pop_df$env_nas <- base::rowSums(pop_df[,c("hyde_crng_nas", "hyde_grz_nas",
                                              "luh2_norng_nas", "luh2_rng_nas",
                                              "cru_nas", "ipsl_nas")])
    
    return(pop_df)
} 


mam_lpd1 <- env_cov_check(mam_lpd1, env_df_sub)
bird_lpd1 <- env_cov_check(bird_lpd1, env_df_sub)

sum(mam_lpd1$env_nas>0) # 64
unique(mam_lpd1$hyde_crng_nas) # 0 114
unique(mam_lpd1$hyde_grz_nas) # 0 114
unique(mam_lpd1$luh2_norng_nas) # 0 114
unique(mam_lpd1$luh2_rng_nas) # 0 114

unique(mam_lpd1$cru_nas) # 0 114
unique(mam_lpd1$ipsl_nas) # 0

sum(bird_lpd1$env_nas>0) # 730
unique(bird_lpd1$hyde_crng_nas) # 0 114
unique(bird_lpd1$hyde_grz_nas) # 0 114
unique(bird_lpd1$luh2_norng_nas) # 0 114
unique(bird_lpd1$luh2_rng_nas) # 0 114

unique(bird_lpd1$cru_nas)  # 0 114
unique(bird_lpd1$ipsl_nas) # 0
# Each pop either has all the data or no data for a given env var


# Subset to pops with all data 
mam_lpd1 <- subset(mam_lpd1, env_nas == 0) 
bird_lpd1 <- subset(bird_lpd1, env_nas == 0)


min(round(bird_lpd1$TS_strt - bird_lpd1$GenLength*3.1)) # 1902
min(round(mam_lpd1$TS_strt - mam_lpd1$GenLength*2.3)) # 1902


nrow(mam_lpd1) # 1305
length(unique(mam_lpd1$Binomial)) # 347
length(unique(mam_lpd1$loc_id)) # 574

nrow(bird_lpd1) # 1230
length(unique(bird_lpd1$Binomial)) # 596
length(unique(bird_lpd1$loc_id)) # 319


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate rates of population change (lambda)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Rate calculations
lambda_calc <- function(lpi_dat){
    # fix pop vals if 0s
    lpi_dat$abund[lpi_dat$abund==0] <- 0.01*mean(lpi_dat$abund[lpi_dat$abund!=0])
    # +1 if pop vals < 1
    if (sum(lpi_dat$abund<1)>0){
        lpi_dat$abund <- lpi_dat$abund + 1
    }
    
    lpi_dat$log_abund <- log10(lpi_dat$abund)
    
    # if n_dat >=6
    gam_rsq <- NA
    gam_lambda <- NA
    gam_pred <- NA
    
    if (nrow(lpi_dat)>=6){
        # GAM
        smth_par <- max(round(nrow(lpi_dat)/2),3)
        s <- mgcv::`s`
        tmp_gam <- mgcv::gam(log_abund~s(Year, k = smth_par), fx = T, data = lpi_dat)
        gam_rsq <- summary(tmp_gam)$r.sq
        # Predict from pop start to end_yr
        gam_pred <- data.frame("Year"=min(lpi_dat$Year):max(lpi_dat$Year))
        gam_pred$abund <- predict(tmp_gam, gam_pred)
        gam_lambda <- mean(diff(gam_pred$abund), na.rm = T)
    }
    # LM
    tmp_lm <- lm(log_abund~Year, data = lpi_dat)
    lm_rsq <- summary(tmp_lm)$adj.r.squared
    # lm lambda is just slope...
    lm_lambda <- as.numeric(tmp_lm$coefficients[2])
    # Also record se, useful in post_df
    lm_se <- summary(tmp_lm)$coefficients[2,2]
    # Set rsq of lm to 1 if only 2 data points
    if (nrow(lpi_dat)==2 & is.nan(lm_rsq)){
        lm_rsq <- 1
    }
    # Also set to 1 if lm_se == 0 and rsq is nan
    if (lm_se == 0 & is.nan(lm_rsq)){
        lm_rsq <- 1
    }
    
    # Interpolation based lambda...
    tmp_abund <- data.frame("Year" = min(lpi_dat$Year):max(lpi_dat$Year))
    tmp_abund <- merge(tmp_abund, lpi_dat,
                       by = "Year", all.x = T)
    tmp_abund$abund_int <- na.approx(tmp_abund$log_abund)
    interp_lambda <- mean(diff(tmp_abund$abund_int), na.rm = T)
    
    # return info
    return(list(gam_lambda = gam_lambda, 
                gam_pred = gam_pred,
                gam_rsq = gam_rsq, 
                lm_lambda = lm_lambda, 
                lm_rsq = lm_rsq,
                interp_lambda = interp_lambda,
                lm_pred = tmp_abund))
}

rate_calc_wrap_ <- function(df){
    
    for (i in 1:nrow(df)){
        
        if(i %% 100 == 0){
            print(paste0("Processing population: ", i))
        }
        
        tmp_loc <- df$loc_id_char[i]
        
        tmp_lpi_df <- data.frame("Year" = 1950:2014,
                                 "abund" = unname(t(df[i,c(which(colnames(df) == "X1950"):
                                                     which(colnames(df) == "X2014"))])))
        
        tmp_lpi_df <- na.omit(tmp_lpi_df)
        
        n_dat <- nrow(tmp_lpi_df)
        strt <- min(tmp_lpi_df$Year)
        end <- max(tmp_lpi_df$Year)
        
        # Calc new lambda as a check
        lambda_dat <- lambda_calc(tmp_lpi_df)
        
        df$gam_lambda[i] <- lambda_dat$gam_lambda
        df$gam_rsq[i] <- lambda_dat$gam_rsq
        df$lm_lambda[i] <- lambda_dat$lm_lambda
        df$lm_rsq[i] <- lambda_dat$lm_rsq
        df$interp_lambda[i] <- lambda_dat$interp_lambda
        
    }
    
    # tidy lambda - gam unless NA, then interp...
    df$lambda <- df$gam_lambda
    df$lambda[is.na(df$lambda)] <- df$interp_lambda[is.na(df$lambda)]
    
    # tidy rsq - gam unless NA, then lm...
    df$lambda_rsq <- df$gam_rsq
    df$lambda_rsq[is.na(df$lambda_rsq)] <- df$lm_rsq[is.na(df$lambda_rsq)]
    
    return(df)
}


bird_rate_df <- rate_calc_wrap_(bird_lpd1)
mam_rate_df <- rate_calc_wrap_(mam_lpd1)

sum(is.na(bird_rate_df$lambda))
sum(is.na(mam_rate_df$lambda))

bird_rate_df <- subset(bird_rate_df, !is.na(lambda))

min(round(bird_rate_df$TS_strt - bird_rate_df$GenLength*3.1)) # 1902


# Some quick overview plots
hist(mam_rate_df$lambda, 150)
hist(bird_rate_df$lambda, 150)

hist(mam_rate_df$DP_to_2014, 50)
hist(bird_rate_df$DP_to_2014, 50)

hist(mam_rate_df$TS_length_to_2014, 50)
hist(bird_rate_df$TS_length_to_2014, 50)

sum(is.na(mam_rate_df$lambda))
sum(is.na(bird_rate_df$lambda)) 

# Check bodymass etc NAs
sum(is.na(mam_rate_df$adult_body_mass_g))
sum(is.na(mam_rate_df$Protected_status))
sum(is.na(mam_rate_df$loc_id))
sum(is.na(mam_rate_df$Binomial))
sum(is.na(mam_rate_df$GenLength))
# All 0: no nas
unique(mam_rate_df$Protected_status)

sum(is.na(bird_rate_df$adult_body_mass_g))
sum(is.na(bird_rate_df$Protected_status))
sum(is.na(bird_rate_df$loc_id))
sum(is.na(bird_rate_df$Binomial))
sum(is.na(bird_rate_df$GenLength))
# All 0: non nas
unique(bird_rate_df$Protected_status)


# saveRDS(mam_rate_df, "../Data/mam_lpd_lambda_1950_2014.rds")
# saveRDS(bird_rate_df, "../Data/bird_lpd_lambda_1950_2014.rds")


# Subset env data to pops/locs to be used in analysis
env_df_sub1 <- subset(env_df_sub, loc_id %in% unique(c(bird_rate_df$loc_id,
                                                       mam_rate_df$loc_id)))

# process mam and bird dfs ready for modelling
# log10 bodymass
# pa: grep no <- no
data_prep <- function(df){
    # df <- subset(df, !is.na(lambda) & adult_body_mass_g != -999)
    df$bm <- log10(df$adult_body_mass_g)
    df$pa <- "Yes"
    df$pa[grepl("No", df$Protected_status)] <- "No"
    
    return(df)
}

bird_sub_df <- data_prep(bird_rate_df) # 1229
mam_sub_df <- data_prep(mam_rate_df) # 1305

length(unique(mam_sub_df$loc_id))
length(unique(bird_sub_df$loc_id))

# restrict mam/bird data to cols of interest
colnames(mam_sub_df)
keep_cols <- c("ID", "Binomial", "lambda", "lambda_rsq", "loc_id", "pa", "bm", 
               "Class", "Order", "Family", "Genus", "Species", "Location", 
               "Country", "Latitude", "Longitude", "Realm", "LPI_Realm",
               "GenLength", "DP_to_2014", "TS_length_to_2014", 
               "gam_lambda","gam_rsq", "lm_lambda", "lm_rsq", 
               "interp_lambda", "TS_strt", "TS_end", "Sys_realm", "realm_repl",
               "Herb_pcnt", "Carn_pcnt")

# saveRDS(list(mam_df = mam_sub_df[,keep_cols],
#              bird_df = bird_sub_df[,keep_cols],
#              env_df = env_df_sub1),
#         "../Data/lpi_lag_data_1901_2014.rds")


##%%%%%
## Generate anonymised data
##%%%%%

# drop country, Location, Management Management_type IPBES_region IPBES_subregion LPI_Region
keep_cols1 <- c("ID", "Binomial", "Class", "Order", "Family", "Genus", "Species", 
                "DP_to_2014", "TS_length_to_2014", "TS_strt", "TS_end", 
                "GenLength", "bm", "Herb_pcnt", "Carn_pcnt",  
                "Managed", "Utilised", "pa", 
                "lambda_rsq", "gam_rsq", "lm_rsq",
                "lambda", "gam_lambda", "interp_lambda", "lm_lambda",
                "loc_id", "Latitude", "Realm", "LPI_Realm", "Confidential")

mam_df_anon <- mam_sub_df[,keep_cols1]
bird_df_anon <- bird_sub_df[,keep_cols1]

# sort loc_id, Latitude Longitude for use in modelling..
unique_locs <- unique(c(mam_df_anon$loc_id, bird_df_anon$loc_id))
unique_locs_df <- data.frame(loc_id = unique_locs,
                             loc_idx = as.character(1:length(unique_locs)),
                             stringsAsFactors = F)

mam_df_anon <- left_join(mam_df_anon, unique_locs_df, all.x = T,
                         by = "loc_id")
bird_df_anon <- left_join(bird_df_anon, unique_locs_df, all.x = T,
                          by = "loc_id")

# saveRDS(unique_locs_df, "../Data/unique_loc_lookup.rds")


# and do for conf spp/pops
conf_mam <- unique(mam_df_anon$Binomial[mam_df_anon$Confidential==1]) 
nonconf_mam <- setdiff(unique(mam_df_anon$Binomial), 
                       unique(mam_df_anon$Binomial[mam_df_anon$Confidential==1]))

unique_mam_spp_df <- data.frame(Binomial = c(conf_mam,
                                             nonconf_mam),
                                spp_idx = c(paste0("spp_",as.character(1:length(conf_mam))),
                                            nonconf_mam),
                                stringsAsFactors = F)

conf_bird <- unique(bird_df_anon$Binomial[bird_df_anon$Confidential==1]) 
nonconf_bird <- setdiff(unique(bird_df_anon$Binomial), 
                        unique(bird_df_anon$Binomial[bird_df_anon$Confidential==1]))

unique_bird_spp_df <- data.frame(Binomial = c(conf_bird,
                                              nonconf_bird),
                                 spp_idx = c(paste0("spp_",as.character(1:length(conf_bird))),
                                             nonconf_bird),
                                 stringsAsFactors = F)

mam_df_anon <- left_join(mam_df_anon, unique_mam_spp_df, all.x = T,
                           by = "Binomial")
bird_df_anon <- left_join(bird_df_anon, unique_bird_spp_df, all.x = T,
                            by = "Binomial")

# saveRDS(unique_mam_spp_df, "../Data/unique_mam_lookup.rds")
# saveRDS(unique_bird_spp_df, "../Data/unique_bird_lookup.rds")

mam_df_anon$trop_temp <- NA
mam_df_anon$trop_temp[(mam_df_anon$Latitude < 23.5) & 
                            (mam_df_anon$Latitude > -23.5)] <- "trop"
mam_df_anon$trop_temp[(mam_df_anon$Latitude > 23.5) | 
                            (mam_df_anon1Latitude < -23.5)] <- "tmpr"


bird_df_anon$trop_temp <- NA
bird_df_anon$trop_temp[(bird_df_anon$Latitude < 23.5) & 
                             (bird_df_anon$Latitude > -23.5)] <- "trop"
bird_df_anon$trop_temp[(bird_df_anon$Latitude > 23.5) | 
                             (bird_df_anon$Latitude < -23.5)] <- "tmpr"

keep_cols2 <- c("ID", "spp_idx", "Class",
                "DP_to_2014", "TS_length_to_2014", "TS_strt", "TS_end", 
                "GenLength", "bm", "Herb_pcnt", "Carn_pcnt",  
                "Managed", "Utilised", "pa", 
                "lambda_rsq", "gam_rsq", "lm_rsq",
                "lambda", "gam_lambda", "interp_lambda", "lm_lambda",
                "loc_idx", "Realm", "LPI_Realm", "trop_temp")

mam_df_anon <- mam_df_anon[,keep_cols2]
bird_df_anon <- bird_df_anon[,keep_cols2]



env_df_anon <- left_join(env_df_sub1, unique_locs_df,
                         by = "loc_id")

saveRDS(list(mam_df = mam_df_anon,
             bird_df = bird_df_anon,
             env_df = env_df_anon[,c("Year", "loc_idx", "hyde3.2_crng", "hyde3.2_grz", 
                                     "luh2_norng", "luh2_rng", "cru4.04", "ipsl")]),
        "../Data/anon_dat.rds")



###%%%%%
# Set up specs for lag models
###%%%%%

class_vec <- c("bird", "mammal")
lag_vec <- 0:49
bird_gen_vec <- seq(0.3,3.1,0.1)
mam_gen_vec <- gen_vec <- seq(0.3,2.3,0.1)

spec1 <- expand.grid("class" = class_vec,
            "cc_lag" = lag_vec,
            "luc_lag" = lag_vec,
            "dur" = "ts",
            "off" = NA,
            "type" = "year")
spec2 <- expand.grid("class" = "mammal",
            "cc_lag" = mam_gen_vec,
            "luc_lag" = mam_gen_vec,
            "dur" = "ts",
            "off" = NA,
            "type" = "gen")
spec3 <- expand.grid("class" = "bird",
            "cc_lag" = bird_gen_vec,
            "luc_lag" = bird_gen_vec,
            "dur" = "ts",
            "off" = NA,
            "type" = "gen")
spec_df <- rbind(spec1, spec2, spec3)
spec_df$id <- 1:nrow(spec_df)
nrow(spec_df)

saveRDS(spec_df, "../Data/lag_specs.rds")

