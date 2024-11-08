Supplementary information / reproducible research files for the manuscript 
Title: "Landmarking for left-truncated competing risk data"

Authors: Theresa Unseld, Tobias Bluhmki, Jan Beyersmann, Evelin Beck, Stephanie Padberg, and Regina Stegherr
Code was written by Unseld, T.
In case of questions or comments please contact theresa.unseld@uni-ulm.de
An online version of the R code can be found at https://github.com/unseldth/LandmarkingCode

The code was written/evaluated in R with the following software versions:
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8    LC_MONETARY=German_Germany.utf8
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3       dplyr_1.1.4         compareGroups_4.8.0 mvna_2.0.1          latex2exp_0.9.6    
 [6] doParallel_1.0.17   iterators_1.0.14    foreach_1.5.2       RColorBrewer_1.1-3  scales_1.3.0       
[11] ggplot2_3.4.4       xtable_1.8-4        data.table_1.14.10  mstate_0.3.2        survival_3.5-7     

loaded via a namespace (and not attached):
 [1] writexl_1.5.0           rlang_1.1.3             magrittr_2.0.3          compiler_4.3.2         
 [5] systemfonts_1.0.5       vctrs_0.6.5             rvest_1.0.3             stringr_1.5.1          
 [9] httpcode_0.3.0          pkgconfig_2.0.3         shape_1.4.6             crayon_1.5.2           
[13] fastmap_1.1.1           backports_1.4.1         ellipsis_0.3.2          labeling_0.4.3         
[17] utf8_1.2.4              promises_1.2.1          rmarkdown_2.25          nloptr_2.0.3           
[21] ragg_1.2.7              purrr_1.0.2             xfun_0.41               glmnet_4.1-8           
[25] jomo_2.7-6              jsonlite_1.8.8          later_1.3.2             uuid_1.2-0             
[29] pan_1.9                 broom_1.0.5             R6_2.5.1                stringi_1.8.3          
[33] Rsolnp_1.16             boot_1.3-28.1           rpart_4.1.23            Rcpp_1.0.12            
[37] knitr_1.45              httpuv_1.6.13           Matrix_1.6-5            splines_4.3.2          
[41] nnet_7.3-19             tidyselect_1.2.0        rstudioapi_0.15.0       codetools_0.2-19       
[45] curl_5.2.0              lattice_0.22-5          tibble_3.2.1            shiny_1.8.0            
[49] withr_3.0.0             flextable_0.9.5         askpass_1.2.0           evaluate_0.23          
[53] zip_2.3.1               xml2_1.3.6              pillar_1.9.0            mice_3.16.0            
[57] generics_0.1.3          truncnorm_1.0-9         munsell_0.5.0           minqa_1.2.6            
[61] chron_2.3-61            glue_1.7.0              gdtools_0.3.6           tools_4.3.2            
[65] gfonts_0.2.0            lme4_1.1-35.1           webshot_0.5.5           grid_4.3.2             
[69] tidyr_1.3.0             colorspace_2.1-0        nlme_3.1-164            cli_3.6.2              
[73] kableExtra_1.3.4        textshaping_0.3.7       officer_0.6.5           fontBitstreamVera_0.1.1
[77] fansi_1.0.6             viridisLite_0.4.2       svglite_2.1.3           gtable_0.3.4           
[81] digest_0.6.34           fontquiver_0.2.1        crul_1.4.0              farver_2.1.1           
[85] htmltools_0.5.7         lifecycle_1.0.4         httr_1.4.7              mitml_0.4-5            
[89] mime_0.12               fontLiberation_0.1.0    openssl_2.1.1           MASS_7.3-60.0.1        
[93] HardyWeinberg_1.7.5    

All simulations were originally run on Ulm University's cluster "Pacioli" with R Version 3.2.3 (December, 2015) and the corresponding package versions 
data.table_1.9.6 mstate_0.2.7 doParallel_1.0.10

The random numbers of this earlier R version are obtained by adding "RNGkind(sample.kind = "Rounding")" to the simulation scripts.
Spot checks are implemented by setting the pacioli variable in the simulation scripts (SimuRun) to FALSE.

-------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------

This folder contains the following data and files that can be used to reproduce all analysis and figures of the manuscript.
It contains the following files:

-------------------------------------------------------------------------------------------------------------------------------
Data sets as versions of the anonymised Berlin fluoroquinolone data in the "./Data/" folder that mimic the main features of the real data analysed in the article and produce similar results: 

Input data set:
Data/imp_anonymized.Rds	Anonymised version of the Berlin fluoroquinolone data from the Berlin Institute for Clinical Teratology and Drug Risk Assessment in Pregnancy.
The original data can not be published online due to data privacy reasons.

Different versions of this dataset in different mult-state and landmark formats:
Data/dMSM_anonymized.Rds (separate currently and previously exposed state)
Data/dMSM2.Rds (combined currently and previously exposed state)
Data/LMlong.Rds (format for landmark anlaysis)
Data/dMSMexpand.Rds (format for nonparametric estimation)

-------------------------------------------------------------------------------------------------------------------------------
R scripts for the analysis of the anonymised Berlin fluoroquinolone data: 
01ConvertoMstate.R	Conversion from wide format into a multi-state model in long format
02MstateAnalysisGyrase.R Non-parametric estimation of transition probabilities for Figure 3
03TDCAnalysisGyrase.R Hazard Ratios from a Cox model with exposure as time-dependent covariate (Table 1), transition probabilites in this model for Figure 3
04LandmarkDatasetsDescriptionGyrase.R	Descriptive analysis of Berlin fluoroquinolone landmark data (Figure 2)
05LMCrude.R Hazard Ratios from a crude landmark analysis (Table 2), transition probabilites in this model for Figure 3
06SupermodelGyrase.R	Hazard Ratios from a landmark supermodel (Table 3), transition probabilites in this model for Figure 3
07CumulativeIncidencesComparison.R Comparison of estimated transition probabilities from the non-parametric model, the TDC model, and the landmark models (Figure 3, Figure S1)
GlobalVariables.R Global variables, like the paths to the directories, variable labels, and transition matrices

These scripts make use of additional functions which can be found in:
Programs/ConvertMSM.R	Functions for conversion from wide format into multi-state model in long format
Programs/HelpFunctionsAnaGyrase.R	Functions for definition and evaluation of the landmark data


-------------------------------------------------------------------------------------------------------------------------------
R scripts for the simulation in the subfolder "./Simulation/":
01SimulationSetUp.R	Definition and description of the parametric models for the simulation (Figure S2), definition of the parameters for the extended simulation with additional scenarios
02Probtrans.R	Calculation of the true transition probabilities for the simulation (with matrix-products adapted from the function probtrans from mstate where the matrices are the true parametric cumulative hazards)
03SimuRun.R	Executes the simulations by calling LMSimulation.R and saves the simulation results
04SimuRunResults.R Descriptive analysis of the landmark datasets and the estimated hazard ratios, and transition probabilities in the simulation(Figure 4, Figure S3, Figure S4, Table 5)
GlobalVariables.R Global variables, like the paths to the directories, variable labels, and transition matrices

R scripts for the simulation with multiple scenarios in the subfolder "./Simulation/MultipleScenarios/:
01SimuRunScs.R	Executes the simulations by calling LMSimulationScs.R and saves simulation results
02SimuRunResultsScs.R Descriptive analysis of the landmark datasets and the estimated hazard ratios, and transition probabilities in the additonal simulation scenarios (Table S3, Figures S5-S24)



These scripts use additional functions which can be found in the subfolders "./RunFiles/":
LMSimulation.R	Generates completely observed multi-state data with the function mssample from the mstate package. Creates left-truncated multi-state data from completely observed data. Creates landmark data 
	from the multi-state data. Estimates coefficients and transition probabilities in the landmark models, coefficient in the TDC model, and transition probabilities in the non-parametric model
AnaLandmark.R	Similar to Programs/HelpFunctionsAnaGyrase.R, but with adaptation for the distincition between the current and previously exposed state

LMSimulationScs.R	Generates left-truncation multi-state data with the function mssample from the mstate package for all combinations of the additional scenarios set= 1,2,...,9 and left-truncation distributions 
	LT= 1,2,3,4. Creates landmark data from multi-state data. Estimates coefficients and transition probabilities in the landmark models, coefficient in the TDC model, and transition probabilities in the non-parametric model
AnaLandmarkScs.R	Similar to AnaLandmark.R but with adaptation for the analysis of the additional scenarios


-------------------------------------------------------------------------------------------------------------------------------

The estimation results of all scripts are saved in the respective "./Results/" folders.
For the simulations, the results folder contains two subfolders:
./Cluster/  original simulation results where the code was run on a cluster
./Test/ test results from spot checks with the first 2 iterations

The estimates for the anonymised data vary slightly from the estimates for the original Berlin fluoroquinolone data, since noise was added to the event times and exposure times of the anonymised data for privacy reasons.

-------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------