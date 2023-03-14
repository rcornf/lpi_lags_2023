README

This repository contains the code needed to reproduce the analyses found within **Cornford et al. 2023, Ongoing over-exploitation and delayed responses to environmental change highlight the urgency for action to promote vertebrate recoveries by 2030**.


---

All analysis was conducted using R 4.1.3  
A copy of the session information and packages used can be found at the bottom of this page.

Figures relevant to the manuscript can be found in Results/Figs.  
All other results and intermediate data files are excluded due to their cumulative size.

---

Before running any code files, please download the relevant environmental data.  
The links below provide access to the data, required local file paths are provided in brackets.

LUH2:  
[Historical](https://luh.umd.edu/LUH2/LUH2_v2h/states.nc)
(Data/LUH2/states.nc)

[link](https://luh.umd.edu/LUH2/LUH2_v2f/IMAGE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc)
(Data/LUH2/ssp1_rcp2.6.nc)

[link](https://luh.umd.edu/LUH2/LUH2_v2f/AIM/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-AIM-ssp370-2-1-f_gn_2015-2100.nc)
(Data/LUH2/ssp3_rcp7.0.nc)

[link](https://luh.umd.edu/LUH2/LUH2_v2f/MAGPIE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc)
(Data/LUH2/ssp5_rcp8.5.nc)


IPSL: 
[link](https://data.isimip.org/datasets/b6ed4b89-43aa-4d79-9206-5b185a356468/)
all files from 1901 onwards are needed
(Data/IPSL/<files>)

[link](https://data.isimip.org/datasets/a7569f3c-1543-46e9-b1ef-cf791cf83859/)
(Data/IPSL_126/< files>)

[link](https://data.isimip.org/datasets/206170bc-41e4-46df-bde5-76a86b815f8d/)
(Data/IPSL_370/< files>)

[link](https://data.isimip.org/datasets/f0fb9bcf-36ec-420d-8590-bd8dd697c7d6/)
(Data/IPSL_585/< files>)


HYDE 3.2: 
[link](https://dataportaal.pbl.nl/downloads/HYDE/HYDE3.2/baseline.zip)
all files should be extracted
(Data/HYDE3.2/< files>)

CRU TS 4.04: 
[link](https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.04/cruts.2004151855.v4.04/tmp/cru_ts4.04.1901.2019.tmp.dat.nc.gz)
file needs to be extracted
(Data/cru_ts4.04.1901.2019.tmp.dat.nc)

An IUCN API token is also required for 01_dataprep.R, one can be requested here: [link](https://apiv3.iucnredlist.org/api/v3/token.)


---

**Overview of code files:**

* 01_analysis_functions.R - Assorted functions used in subsequent analysis.
* 01_dataprep.R – Prepare data for subsequent analyses. Note, the raw LPD file is not provided due to confidential population data. We include an anonymised version of the processed data (populations and environmental) that can be used in the following analysis.  
* 02_initial_models.R – Run initial, lagged models.  
* 03_initial_analysis.R – Analyse the initial models for optimal lags and conduct sensitivity tests.  
* 04_main_models.R – Run final models, having dropped influential populations/species.  
* 05_bm_models.R - Run models separately for species with different body-mass categories (small, medium, large).
* 06_final_analysis.R – Analyse the main and bm-split models to obtain optimal lags and produce main results/figures.  
* 07_cv_top_models.R, 08_dp_models.R, 09_ecol_models.R, 10_env_models.R – Run additional models to check stability of results under cross-validation, when using alternative data selection thresholds, for different ecological subsets, and across different environmental data, respectively.  
* 11_SI_analysis.R – Conduct supplementary analysis and produce SI results/figures.

**Note, model fitting can take a substantial amount of time on a local machine. Using multiple remote CPUs is advised but not accounted for in the provided code.**


---

**Overview of data files:**

* Amniote_Database_Aug_2015.csv - Source of body mass darta.
* anon_dat.rds - Anonymised version of data used to fit models.
* bird_gen_length_Bird2020.xlsx - Source of bird generation lengths.
* BirdFuncDat.txt - Source of bird diet types.
* cv_specs.rds - Specification for cross-validation of best lag-based models.
* env_126.rds - Environmental data for SSP1 RCP 2.6, pre extracted/procesed.
* env_370.rds - Environmental data for SSP3 RCP 7.0, pre extracted/procesed.
* env_585.rds - Environmental data for SSP5 RCP 8.5, pre extracted/procesed.
* env_hist.rds - Historical environemntal data.
* esa_data.rds - Agricultural land-cover estimates from ESA.
* lag_specs.rds - Specification for lag-based models.
* mam_gen_length_Pacifici2013.xls - Source of mammal generation lengths.
* MamFuncDat.txt - Source of mammal diet types.


---

**License**

cc...


---

**sessionInfo()**

R version 4.1.3 (2022-03-10)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] doParallel_1.0.17  iterators_1.0.14   foreach_1.5.2      ncdf4_1.19         zoo_1.8-10        
 [6] car_3.0-13         carData_3.0-5      mgcv_1.8-39        nlme_3.1-155       ggtreeExtra_1.4.2 
[11] ggtree_3.2.1       rotl_3.0.12        phytools_1.0-3     maps_3.4.0         ape_5.6-2         
[16] geosphere_1.5-14   sf_1.0-7           fasterize_1.0.3    terra_1.5-21       raster_3.5-15     
[21] sp_1.4-7           influence.ME_0.9-9 dplyr_1.1.0        ggpubr_0.4.0       cowplot_1.1.1     
[26] gtable_0.3.0       gridExtra_2.3      ggplot2_3.3.6      MuMIn_1.46.0       lme4_1.1-29       
[31] Matrix_1.4-0       plyr_1.8.7         reshape2_1.4.4    

loaded via a namespace (and not attached):
 [1] ggnewscale_0.4.7        minqa_1.2.4             colorspace_2.0-3        ggsignif_0.6.3         
 [5] ellipsis_0.3.2          class_7.3-20            rgdal_1.5-32            aplot_0.1.6            
 [9] rstudioapi_0.13         proxy_0.4-26            farver_2.1.0            fansi_1.0.3            
[13] codetools_0.2-18        splines_4.1.3           mnormt_2.1.0            jsonlite_1.8.0         
[17] nloptr_2.0.0            broom_0.8.0             rentrez_1.2.3           compiler_4.1.3         
[21] httr_1.4.3              backports_1.4.1         assertthat_0.2.1        lazyeval_0.2.2         
[25] cli_3.6.0               prettyunits_1.1.1       tools_4.1.3             igraph_1.3.4           
[29] coda_0.19-4             glue_1.6.2              clusterGeneration_1.3.7 fastmatch_1.1-3        
[33] Rcpp_1.0.8.3            cellranger_1.1.0        vctrs_0.5.2             stringr_1.4.0          
[37] lifecycle_1.0.3         phangorn_2.9.0          rstatix_0.7.0           XML_3.99-0.10          
[41] MASS_7.3-55             scales_1.2.0            hms_1.1.1               expm_0.999-6           
[45] ggfun_0.0.6             yulab.utils_0.0.5       stringi_1.7.6           plotrix_3.8-2          
[49] e1071_1.7-9             tidytree_0.4.0          boot_1.3-28             rlang_1.0.6            
[53] pkgconfig_2.0.3         rncl_0.8.6              lattice_0.20-45         purrr_0.3.4            
[57] treeio_1.18.1           patchwork_1.1.1         labeling_0.4.2          tidyselect_1.2.0       
[61] magrittr_2.0.3          R6_2.5.1                generics_0.1.2          combinat_0.0-8         
[65] DBI_1.1.2               pillar_1.7.0            withr_2.5.0             units_0.8-0            
[69] scatterplot3d_0.3-41    abind_1.4-5             tibble_3.1.7            crayon_1.5.1           
[73] KernSmooth_2.23-20      utf8_1.2.2              progress_1.2.2          readxl_1.4.0           
[77] digest_0.6.29           classInt_0.4-3          tidyr_1.2.0             numDeriv_2016.8-1.1    
[81] gridGraphics_0.5-1      stats4_4.1.3            munsell_0.5.0           ggplotify_0.1.0        
[85] quadprog_1.5-8   