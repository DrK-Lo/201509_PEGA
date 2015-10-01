---
title: "PEGA for Pine Analysis"
author: "Katie Lotterhos"
date: "August 24, 2015"
output: html_document
---

Does the degree to which SNPs evolve to be correlated with the environment depend on the phenotype-environment association of the phenotypes they are associated with?

Basic algorithm:

1. choose environment
2. choose phenotype
3. get PE association
4. get relevant column data for each SNP (bayenv rho, bayenv bf, gwas slope, gwas p-value) using system calls
5. correct rho for different coding from gwas using flip.cor using merge
6. correct rho for PE/GWAS direction (polarize)
    + (PE-correlation sign) * (GWAS-slope sign) = correction sign
7. (remove large values that are not significant)
8. plot GWAS slope (x axis) vs GEA (y axis)
    + color points in 1st and 3rd quadrants in green (signs match expected direction)
    + color points in 2nd and 4th quadrants in purple (signs don't match expected direction)
    
#### Note
Note that the code runs fastest if you work through R in a terminal window and upload the "results_pine" object first.  Once it is uploaded, when the markdown is created it won't try to load the "results_pine" dataframe, which is very large.  To make this markdown file:
```
library(knitr)
library(markdown)
setwd("~/Desktop/CurrResearch/1-AdaptreeData/201509_PEGA/analysis")
knit("AnalysisPEGA.Rmd")  # produces the md file
markdownToHTML("AnalysisPEGA.md", "AnalysisPEGA.html")  # converts an md file to html
```
  
Install packages  

```r
  if (!("hexbin" %in% installed.packages())){install.packages("hexbin")}
  library(hexbin)
  if (!("gridExtra" %in% installed.packages())){install.packages("ash")}
  library(gridExtra)
  if (!("ash" %in% installed.packages())){install.packages("ash")}
  library(ash)
  if (!("fields" %in% installed.packages())){install.packages("fields")}
  library(fields)
```

Load data, if not already loaded.  The phenotype and environment labels were not consistent among studies.  The `matchXXX` objects are used to match inconsistent labeling.

Also load "flip" data.  A summary of the flip (i.e., allele coding) problem can be found in the `Log.Rmd`.  `XXX_flip_gwas` is an indicator (-1,1) to multiply the gwas slopes by to make them match the ways alleles were coded for bayenv GEAs.  `XXX_flip_raw` is an indicator (-1,1) to multiply Sam's raw correlations by to match the ways alleles were coded for the bayenv GEAs.


```r
### load data
PE <- read.csv("../data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2.csv", header=TRUE)
#reorder to make largest effects on top
PE <- PE[rev(order(abs(PE$pine_correlation))),]
head(PE)
```

```
##     Environment                   Phenotype pine_correlation
## 190    Latitude      ColdInjuryFallS02_Mean            -0.42
## 36          CMD                  BudSet_day             0.40
## 275        MCMT      ColdInjuryFallS02_Mean             0.38
## 342         SHM                  BudSet_day             0.37
## 191    Latitude ColdInjuryMidWinterS02_Mean            -0.37
## 54          DD0      ColdInjuryFallS02_Mean            -0.37
##     spruce_correlation
## 190              -0.74
## 36                0.33
## 275               0.76
## 342               0.25
## 191              -0.62
## 54               -0.74
```

```r
plot(PE$pine_correlation, PE$spruce_correlation, pch=20, col="blue")
abline(0,1)
abline(lm(PE$spruce_correlation~PE$pine_correlation), col="blue")
```

![plot of chunk load data](figure/load data-1.png) 

```r
PE_matchGWAS <- read.csv("../data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2_matchGWAS.csv", header=TRUE)
PE_matchGEA <- read.csv("../data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2_matchGEA.csv", header=TRUE)

#### Load PINE Flip and RESULTS files
flip_pine_gwas <- read.table("../data/large_files/flip_pine.cor", header=TRUE)
  nrow(flip_pine_gwas)
```

```
## [1] 1251074
```

```r
  names(flip_pine_gwas)[2] <- "flip_pine_gwas"
  head(flip_pine_gwas)
```

```
##   gcontig__gcontig_pos flip_pine_gwas
## 1       C10572747__100             -1
## 2       C10572747__101             -1
## 3       C10572747__124             -1
## 4        C10572747__20              1
## 5        C10572747__25              1
## 6         C10572747__9             -1
```

```r
flip_pine_raw <- read.table("../data/large_files/var_out_GATK3_allhet_pine688_ALL.table_filt10_p95_het7_passFILT2.2_LDformat_TOFLIP_alphabet_lowest_eq_1")
  nrow(flip_pine_raw)
```

```
## [1] 1089468
```

```r
  names(flip_pine_raw)[1] <- names(flip_pine_gwas)[1]
  names(flip_pine_raw)[2] <- "flip_pine_raw"
  head(flip_pine_raw)
```

```
##   gcontig__gcontig_pos flip_pine_raw
## 1       C10572747__101             1
## 2        C10572747__20            -1
## 3        C10572747__25            -1
## 4         C10572747__9             1
## 5       C11705109__103            -1
## 6        C11705109__89            -1
```

```r
flip_pine2 <- merge(flip_pine_gwas, flip_pine_raw, all.x=TRUE)
  nrow(flip_pine2)
```

```
## [1] 1251134
```

```r
  head(flip_pine2)
```

```
##   gcontig__gcontig_pos flip_pine_gwas flip_pine_raw
## 1       C10572747__100             -1            NA
## 2       C10572747__101             -1             1
## 3       C10572747__124             -1            NA
## 4        C10572747__20              1            -1
## 5        C10572747__25              1            -1
## 6         C10572747__9             -1             1
```

```r
rm(flip_pine_gwas, flip_pine_raw)

#### read in results_pine if it doesn't exist
  if (!("results_pine" %in% ls())){
  results_pine <- read.table("../data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS", header=TRUE, comment.char=" ")
    }
  dim(results_pine)
```

```
## [1] 1256174     359
```

```r
  ### merge flip with results
  results_pine2 <- merge(flip_pine2, results_pine)
  dim(results_pine2)
```

```
## [1] 1332974     359
```

```r
  head(results_pine2[,1:10])
```

```
##   gcontig__gcontig_pos flip_pine_gwas flip_pine_raw   gcontig pos_gcontig
## 1       C10572747__100             -1            NA C10572747         100
## 2       C10572747__101             -1             1 C10572747         101
## 3       C10572747__124             -1            NA C10572747         124
## 4        C10572747__20              1            -1 C10572747          20
## 5        C10572747__25              1            -1 C10572747          25
## 6         C10572747__9             -1             1 C10572747           9
##   X_annotation paralogy_rank          tcontig pos_tcontig gmap_percentID
## 1  non_tcontig             1 ctg7180057699683         163           98.2
## 2  non_tcontig             1 ctg7180057699683         162           98.2
## 3  non_tcontig             1 ctg7180057699683         139           98.2
## 4  non_tcontig             1 ctg7180057699683         243           98.2
## 5  non_tcontig             1 ctg7180057699683         238           98.2
## 6  non_tcontig             1 ctg7180057699683         254           98.2
```

```r
  results_pine <- results_pine2
  rm(results_pine2, flip_pine2)

  results_pine_head <- names(results_pine)
  results_pine_head
```

```
##   [1] "gcontig__gcontig_pos"               
##   [2] "flip_pine_gwas"                     
##   [3] "flip_pine_raw"                      
##   [4] "gcontig"                            
##   [5] "pos_gcontig"                        
##   [6] "X_annotation"                       
##   [7] "paralogy_rank"                      
##   [8] "tcontig"                            
##   [9] "pos_tcontig"                        
##  [10] "gmap_percentID"                     
##  [11] "offtarget_dist"                     
##  [12] "num_alleles"                        
##  [13] "num_ind"                            
##  [14] "num_comm_gtype"                     
##  [15] "major_allele"                       
##  [16] "major_num"                          
##  [17] "major_freq"                         
##  [18] "minor_allele"                       
##  [19] "minor_num"                          
##  [20] "minor_freq"                         
##  [21] "freq_hetzyg"                        
##  [22] "MAP_dir1"                           
##  [23] "MAP_dir2"                           
##  [24] "MAP_MgGm"                           
##  [25] "noncode_in_covar"                   
##  [26] "MQ"                                 
##  [27] "haplotype_score"                    
##  [28] "FS"                                 
##  [29] "MQRankSum"                          
##  [30] "ReadPosRankSum"                     
##  [31] "allele_balance"                     
##  [32] "window_stat"                        
##  [33] "scaffold"                           
##  [34] "scaf_pos"                           
##  [35] "xtx"                                
##  [36] "LAT"                                
##  [37] "LAT_rho1"                           
##  [38] "LAT_rho2"                           
##  [39] "LONG"                               
##  [40] "LONG_rho1"                          
##  [41] "LONG_rho2"                          
##  [42] "ELEVATION"                          
##  [43] "ELEVATION_rho1"                     
##  [44] "ELEVATION_rho2"                     
##  [45] "MAT"                                
##  [46] "MAT_rho1"                           
##  [47] "MAT_rho2"                           
##  [48] "MWMT"                               
##  [49] "MWMT_rho1"                          
##  [50] "MWMT_rho2"                          
##  [51] "MCMT"                               
##  [52] "MCMT_rho1"                          
##  [53] "MCMT_rho2"                          
##  [54] "TD"                                 
##  [55] "TD_rho1"                            
##  [56] "TD_rho2"                            
##  [57] "MAP"                                
##  [58] "MAP_rho1"                           
##  [59] "MAP_rho2"                           
##  [60] "MSP"                                
##  [61] "MSP_rho1"                           
##  [62] "MSP_rho2"                           
##  [63] "AHM"                                
##  [64] "AHM_rho1"                           
##  [65] "AHM_rho2"                           
##  [66] "SHM"                                
##  [67] "SHM_rho1"                           
##  [68] "SHM_rho2"                           
##  [69] "DD_0"                               
##  [70] "DD_0_rho1"                          
##  [71] "DD_0_rho2"                          
##  [72] "DD5"                                
##  [73] "DD5_rho1"                           
##  [74] "DD5_rho2"                           
##  [75] "NFFD"                               
##  [76] "NFFD_rho1"                          
##  [77] "NFFD_rho2"                          
##  [78] "bFFP"                               
##  [79] "bFFP_rho1"                          
##  [80] "bFFP_rho2"                          
##  [81] "eFFP"                               
##  [82] "eFFP_rho1"                          
##  [83] "eFFP_rho2"                          
##  [84] "FFP"                                
##  [85] "FFP_rho1"                           
##  [86] "FFP_rho2"                           
##  [87] "PAS"                                
##  [88] "PAS_rho1"                           
##  [89] "PAS_rho2"                           
##  [90] "EMT"                                
##  [91] "EMT_rho1"                           
##  [92] "EMT_rho2"                           
##  [93] "EXT"                                
##  [94] "EXT_rho1"                           
##  [95] "EXT_rho2"                           
##  [96] "Eref"                               
##  [97] "Eref_rho1"                          
##  [98] "Eref_rho2"                          
##  [99] "CMD"                                
## [100] "CMD_rho1"                           
## [101] "CMD_rho2"                           
## [102] "LAT_z_lfmm_beagle"                  
## [103] "LAT_q_lfmm_beagle"                  
## [104] "LONG_z_lfmm_beagle"                 
## [105] "LONG_q_lfmm_beagle"                 
## [106] "ELEVATION_z_lfmm_beagle"            
## [107] "ELEVATION_q_lfmm_beagle"            
## [108] "MAT_z_lfmm_beagle"                  
## [109] "MAT_q_lfmm_beagle"                  
## [110] "MWMT_z_lfmm_beagle"                 
## [111] "MWMT_q_lfmm_beagle"                 
## [112] "MCMT_z_lfmm_beagle"                 
## [113] "MCMT_q_lfmm_beagle"                 
## [114] "TD_z_lfmm_beagle"                   
## [115] "TD_q_lfmm_beagle"                   
## [116] "MAP_z_lfmm_beagle"                  
## [117] "MAP_q_lfmm_beagle"                  
## [118] "MSP_z_lfmm_beagle"                  
## [119] "MSP_q_lfmm_beagle"                  
## [120] "AHM_z_lfmm_beagle"                  
## [121] "AHM_q_lfmm_beagle"                  
## [122] "SHM_z_lfmm_beagle"                  
## [123] "SHM_q_lfmm_beagle"                  
## [124] "DD_0_z_lfmm_beagle"                 
## [125] "DD_0_q_lfmm_beagle"                 
## [126] "DD5_z_lfmm_beagle"                  
## [127] "DD5_q_lfmm_beagle"                  
## [128] "NFFD_z_lfmm_beagle"                 
## [129] "NFFD_q_lfmm_beagle"                 
## [130] "bFFP_z_lfmm_beagle"                 
## [131] "bFFP_q_lfmm_beagle"                 
## [132] "eFFP_z_lfmm_beagle"                 
## [133] "eFFP_q_lfmm_beagle"                 
## [134] "FFP_z_lfmm_beagle"                  
## [135] "FFP_q_lfmm_beagle"                  
## [136] "PAS_z_lfmm_beagle"                  
## [137] "PAS_q_lfmm_beagle"                  
## [138] "EMT_z_lfmm_beagle"                  
## [139] "EMT_q_lfmm_beagle"                  
## [140] "EXT_z_lfmm_beagle"                  
## [141] "EXT_q_lfmm_beagle"                  
## [142] "Eref_z_lfmm_beagle"                 
## [143] "Eref_q_lfmm_beagle"                 
## [144] "CMD_z_lfmm_beagle"                  
## [145] "CMD_q_lfmm_beagle"                  
## [146] "bayenv_PCA1"                        
## [147] "bayenv_PCA1_rho1"                   
## [148] "bayenv_PCA1_rho2"                   
## [149] "bayenv_PCA2"                        
## [150] "bayenv_PCA2_rho1"                   
## [151] "bayenv_PCA2_rho2"                   
## [152] "bayenv_PCA3"                        
## [153] "bayenv_PCA3_rho1"                   
## [154] "bayenv_PCA3_rho2"                   
## [155] "bayenv_PCA4"                        
## [156] "bayenv_PCA4_rho1"                   
## [157] "bayenv_PCA4_rho2"                   
## [158] "bayenv_PCA5"                        
## [159] "bayenv_PCA5_rho1"                   
## [160] "bayenv_PCA5_rho2"                   
## [161] "bayenv_PCA6"                        
## [162] "bayenv_PCA6_rho1"                   
## [163] "bayenv_PCA6_rho2"                   
## [164] "bayenv_PCA7"                        
## [165] "bayenv_PCA7_rho1"                   
## [166] "bayenv_PCA7_rho2"                   
## [167] "bayenv_PCA8"                        
## [168] "bayenv_PCA8_rho1"                   
## [169] "bayenv_PCA8_rho2"                   
## [170] "bayenv_PCA9"                        
## [171] "bayenv_PCA9_rho1"                   
## [172] "bayenv_PCA9_rho2"                   
## [173] "bayenv_PCA10"                       
## [174] "bayenv_PCA10_rho1"                  
## [175] "bayenv_PCA10_rho2"                  
## [176] "bayenv_PCA11"                       
## [177] "bayenv_PCA11_rho1"                  
## [178] "bayenv_PCA11_rho2"                  
## [179] "bayenv_PCA12"                       
## [180] "bayenv_PCA12_rho1"                  
## [181] "bayenv_PCA12_rho2"                  
## [182] "bayenv_PCA13"                       
## [183] "bayenv_PCA13_rho1"                  
## [184] "bayenv_PCA13_rho2"                  
## [185] "bayenv_PCA14"                       
## [186] "bayenv_PCA14_rho1"                  
## [187] "bayenv_PCA14_rho2"                  
## [188] "bayenv_PCA15"                       
## [189] "bayenv_PCA15_rho1"                  
## [190] "bayenv_PCA15_rho2"                  
## [191] "bayenv_PCA16"                       
## [192] "bayenv_PCA16_rho1"                  
## [193] "bayenv_PCA16_rho2"                  
## [194] "bayenv_PCA17"                       
## [195] "bayenv_PCA17_rho1"                  
## [196] "bayenv_PCA17_rho2"                  
## [197] "bayenv_PCA18"                       
## [198] "bayenv_PCA18_rho1"                  
## [199] "bayenv_PCA18_rho2"                  
## [200] "bayenv_PCA19"                       
## [201] "bayenv_PCA19_rho1"                  
## [202] "bayenv_PCA19_rho2"                  
## [203] "bayenv_PCA20"                       
## [204] "bayenv_PCA20_rho1"                  
## [205] "bayenv_PCA20_rho2"                  
## [206] "bayenv_PCA21"                       
## [207] "bayenv_PCA21_rho1"                  
## [208] "bayenv_PCA21_rho2"                  
## [209] "bayenv_PCA22"                       
## [210] "bayenv_PCA22_rho1"                  
## [211] "bayenv_PCA22_rho2"                  
## [212] "Budset_p"                           
## [213] "Budset_snp_effect"                  
## [214] "Budbreak_p"                         
## [215] "Budbreak_snp_effect"                
## [216] "Height_season_1_p"                  
## [217] "Height_season_1_snp_effect"         
## [218] "Height_season_2_p"                  
## [219] "Height_season_2_snp_effect"         
## [220] "Diameter_p"                         
## [221] "Diameter_snp_effect"                
## [222] "Shoot_weight_p"                     
## [223] "Shoot_weight_snp_effect"            
## [224] "Root_weight_p"                      
## [225] "Root_weight_snp_effect"             
## [226] "Max_growth_rate_p"                  
## [227] "Max_growth_rate_snp_effect"         
## [228] "Linear_growth_days_p"               
## [229] "Linear_growth_days_snp_effect"      
## [230] "X5_growth_complete_days_p"          
## [231] "X5_growth_complete_days_snp_effect" 
## [232] "X95_growth_complete_days_p"         
## [233] "X95_growth_complete_days_snp_effect"
## [234] "X5_95_growth_days_p"                
## [235] "X5_95_growth_days_snp_effect"       
## [236] "Fall_cold_injury_p"                 
## [237] "Fall_cold_injury_snp_effect"        
## [238] "Winter_cold_injury_p"               
## [239] "Winter_cold_injury_snp_effect"      
## [240] "Spring_cold_injury_p"               
## [241] "Spring_cold_injury_snp_effect"      
## [242] "root_wt__shoot_wt_p"                
## [243] "root_wt__shoot_wt_snp_effect"       
## [244] "root_wt__shoot_wt_p_1"              
## [245] "root_wt__shoot_wt_snp_effect_1"     
## [246] "gcta_PC1_p"                         
## [247] "gcta_PC1_snp_effect"                
## [248] "gcta_PC2_p"                         
## [249] "gcta_PC2_snp_effect"                
## [250] "gcta_PC3_p"                         
## [251] "gcta_PC3_snp_effect"                
## [252] "gcta_PC4_p"                         
## [253] "gcta_PC4_snp_effect"                
## [254] "gcta_PC5_p"                         
## [255] "gcta_PC5_snp_effect"                
## [256] "gcta_PC6_p"                         
## [257] "gcta_PC6_snp_effect"                
## [258] "gcta_PC7_p"                         
## [259] "gcta_PC7_snp_effect"                
## [260] "gcta_PC8_p"                         
## [261] "gcta_PC8_snp_effect"                
## [262] "gcta_PC9_p"                         
## [263] "gcta_PC9_snp_effect"                
## [264] "gcta_PC10_p"                        
## [265] "gcta_PC10_snp_effect"               
## [266] "gcta_PC11_p"                        
## [267] "gcta_PC11_snp_effect"               
## [268] "gcta_PC12_p"                        
## [269] "gcta_PC12_snp_effect"               
## [270] "gcta_PC13_p"                        
## [271] "gcta_PC13_snp_effect"               
## [272] "gcta_PC14_p"                        
## [273] "gcta_PC14_snp_effect"               
## [274] "gcta_PC15_p"                        
## [275] "gcta_PC15_snp_effect"               
## [276] "gcta_PC16_p"                        
## [277] "gcta_PC16_snp_effect"               
## [278] "category_snp_array"                 
## [279] "sequenced_contig_length"            
## [280] "affx_name"                          
## [281] "ax_name"                            
## [282] "LAT_raw_rho"                        
## [283] "LONG_raw_rho"                       
## [284] "ELEVATION_raw_rho"                  
## [285] "MAT_raw_rho"                        
## [286] "MWMT_raw_rho"                       
## [287] "MCMT_raw_rho"                       
## [288] "TD_raw_rho"                         
## [289] "MAP_raw_rho"                        
## [290] "MSP_raw_rho"                        
## [291] "AHM_raw_rho"                        
## [292] "SHM_raw_rho"                        
## [293] "DD_0_raw_rho"                       
## [294] "DD5_raw_rho"                        
## [295] "NFFD_raw_rho"                       
## [296] "bFFP_raw_rho"                       
## [297] "eFFP_raw_rho"                       
## [298] "FFP_raw_rho"                        
## [299] "PAS_raw_rho"                        
## [300] "EMT_raw_rho"                        
## [301] "EXT_raw_rho"                        
## [302] "Eref_raw_rho"                       
## [303] "CMD_raw_rho"                        
## [304] "LAT_raw_p"                          
## [305] "LONG_raw_p"                         
## [306] "ELEVATION_raw_p"                    
## [307] "MAT_raw_p"                          
## [308] "MWMT_raw_p"                         
## [309] "MCMT_raw_p"                         
## [310] "TD_raw_p"                           
## [311] "MAP_raw_p"                          
## [312] "MSP_raw_p"                          
## [313] "AHM_raw_p"                          
## [314] "SHM_raw_p"                          
## [315] "DD_0_raw_p"                         
## [316] "DD5_raw_p"                          
## [317] "NFFD_raw_p"                         
## [318] "bFFP_raw_p"                         
## [319] "eFFP_raw_p"                         
## [320] "FFP_raw_p"                          
## [321] "PAS_raw_p"                          
## [322] "EMT_raw_p"                          
## [323] "EXT_raw_p"                          
## [324] "Eref_raw_p"                         
## [325] "CMD_raw_p"                          
## [326] "Budset_p_raw_rho"                   
## [327] "Budbreak_p_raw_rho"                 
## [328] "Height_season_1_p_raw_rho"          
## [329] "Height_season_2_p_raw_rho"          
## [330] "Diameter_p_raw_rho"                 
## [331] "Shoot_weight_p_raw_rho"             
## [332] "Root_weight_p_raw_rho"              
## [333] "Max_growth_rate_p_raw_rho"          
## [334] "Linear_growth_days_p_raw_rho"       
## [335] "X5_growth_complete_days_p_raw_rho"  
## [336] "X95_growth_complete_days_p_raw_rho" 
## [337] "X5_95_growth_days_p_raw_rho"        
## [338] "Fall_cold_injury_p_raw_rho"         
## [339] "Winter_cold_injury_p_raw_rho"       
## [340] "Spring_cold_injury_p_raw_rho"       
## [341] "root_wt_shoot_wt_p_raw_rho"         
## [342] "root_wt_shoot_wt_p_1_raw_rho"       
## [343] "Budset_p_raw_p"                     
## [344] "Budbreak_p_raw_p"                   
## [345] "Height_season_1_p_raw_p"            
## [346] "Height_season_2_p_raw_p"            
## [347] "Diameter_p_raw_p"                   
## [348] "Shoot_weight_p_raw_p"               
## [349] "Root_weight_p_raw_p"                
## [350] "Max_growth_rate_p_raw_p"            
## [351] "Linear_growth_days_p_raw_p"         
## [352] "X5_growth_complete_days_p_raw_p"    
## [353] "X95_growth_complete_days_p_raw_p"   
## [354] "X5_95_growth_days_p_raw_p"          
## [355] "Fall_cold_injury_p_raw_p"           
## [356] "Winter_cold_injury_p_raw_p"         
## [357] "Spring_cold_injury_p_raw_p"         
## [358] "root_wt_shoot_wt_p_raw_p"           
## [359] "root_wt_shoot_wt_p_1_raw_p"
```

```r
#### Load SPRUCE Flip and RESULTS files
flip_spruce_gwas <- read.table("../data/large_files/flip_spruce.cor", header=TRUE)
  nrow(flip_spruce_gwas)
```

```
## [1] 1064437
```

```r
  names(flip_spruce_gwas)[2] <- "flip_spruce_gwas"
  head(flip_spruce_gwas)
```

```
##               snp_id flip_spruce_gwas
## 1 143857680_251__102               -1
## 2 143857680_251__108               -1
## 3 143857680_251__127                1
## 4 143857680_251__143               -1
## 5 143857680_251__163               -1
## 6 143857680_251__181               -1
```

```r
flip_spruce_raw <- read.table("../data/large_files/var_out_GATK3_spruce_ALL.table_filt10_p95_het7_passFILT2.3_LDformat_TOFLIP_alphabet_lowest_eq_1")
  nrow(flip_spruce_raw)
```

```
## [1] 887774
```

```r
  names(flip_spruce_raw)[1] <- names(flip_spruce_gwas)[1]
  names(flip_spruce_raw)[2] <- "flip_spruce_raw"
  head(flip_spruce_raw)
```

```
##               snp_id flip_spruce_raw
## 1 143857680_251__102               1
## 2 143857680_251__108               1
## 3 143857680_251__143               1
## 4 143857680_251__163               1
## 5 143857680_251__181              -1
## 6  143857680_251__79               1
```

```r
flip_spruce2 <- merge(flip_spruce_gwas, flip_spruce_raw, all.x=TRUE)
  nrow(flip_spruce2)
```

```
## [1] 1064437
```

```r
rm(flip_spruce_gwas, flip_spruce_raw)
names(flip_spruce2)[1] <- "gcontig__gcontig_pos"
  head(flip_spruce2)
```

```
##   gcontig__gcontig_pos flip_spruce_gwas flip_spruce_raw
## 1   143857680_251__102               -1               1
## 2   143857680_251__108               -1               1
## 3   143857680_251__127                1              NA
## 4   143857680_251__143               -1               1
## 5   143857680_251__163               -1               1
## 6   143857680_251__181               -1              -1
```

```r
if (!("results_spruce" %in% ls())){
results_spruce <- read.table("../data/large_files/var_out_GATK3_spruce_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS", header=TRUE, comment.char=" ")
  }
  dim(results_spruce)
```

```
## [1] 1064437     316
```

```r
  results_spruce$gcontig__gcontig_pos <- paste(results_spruce$gcontig, results_spruce$pos_gcontig, sep="__")
  ### merge flip with results
  results_spruce2 <- merge(flip_spruce2, results_spruce)
  dim(results_spruce2)
```

```
## [1] 1064437     316
```

```r
  head(results_spruce2[,1:10])
```

```
##   gcontig__gcontig_pos flip_spruce_gwas flip_spruce_raw X.annotation
## 1   143857680_251__102               -1               1         <NA>
## 2   143857680_251__108               -1               1         <NA>
## 3   143857680_251__127                1              NA         <NA>
## 4   143857680_251__143               -1               1         <NA>
## 5   143857680_251__163               -1               1         <NA>
## 6   143857680_251__181               -1              -1         <NA>
##   paralogy_rank tcontig pos_tcontig gmap_percentID offtarget_dist
## 1            NA    <NA>        <NA>             NA             NA
## 2            NA    <NA>        <NA>             NA             NA
## 3            NA    <NA>        <NA>             NA             NA
## 4            NA    <NA>        <NA>             NA             NA
## 5            NA    <NA>        <NA>             NA             NA
## 6            NA    <NA>        <NA>             NA             NA
##         gcontig
## 1 143857680_251
## 2 143857680_251
## 3 143857680_251
## 4 143857680_251
## 5 143857680_251
## 6 143857680_251
```

```r
  results_spruce <- results_spruce2
  rm(results_spruce2, flip_spruce2)
  results_spruce_head <- names(results_spruce)
  results_spruce_head
```

```
##   [1] "gcontig__gcontig_pos"               
##   [2] "flip_spruce_gwas"                   
##   [3] "flip_spruce_raw"                    
##   [4] "X.annotation"                       
##   [5] "paralogy_rank"                      
##   [6] "tcontig"                            
##   [7] "pos_tcontig"                        
##   [8] "gmap_percentID"                     
##   [9] "offtarget_dist"                     
##  [10] "gcontig"                            
##  [11] "pos_gcontig"                        
##  [12] "num_alleles"                        
##  [13] "num_ind"                            
##  [14] "num_comm_gtype"                     
##  [15] "major_allele"                       
##  [16] "major_num"                          
##  [17] "major_freq"                         
##  [18] "minor_allele"                       
##  [19] "minor_num"                          
##  [20] "minor_freq"                         
##  [21] "freq_hetzyg"                        
##  [22] "MAP_dir1"                           
##  [23] "MAP_dir2"                           
##  [24] "MAP_MgGm"                           
##  [25] "noncode_in_covar"                   
##  [26] "MQ"                                 
##  [27] "haplotype_score"                    
##  [28] "FS"                                 
##  [29] "MQRankSum"                          
##  [30] "ReadPosRankSum"                     
##  [31] "allele_balance"                     
##  [32] "window_stat"                        
##  [33] "scaf"                               
##  [34] "pos_scaf"                           
##  [35] "xtx"                                
##  [36] "LAT"                                
##  [37] "LAT_rho1"                           
##  [38] "LAT_rho2"                           
##  [39] "LONG"                               
##  [40] "LONG_rho1"                          
##  [41] "LONG_rho2"                          
##  [42] "ELEVATION"                          
##  [43] "ELEVATION_rho1"                     
##  [44] "ELEVATION_rho2"                     
##  [45] "MAT"                                
##  [46] "MAT_rho1"                           
##  [47] "MAT_rho2"                           
##  [48] "MWMT"                               
##  [49] "MWMT_rho1"                          
##  [50] "MWMT_rho2"                          
##  [51] "MCMT"                               
##  [52] "MCMT_rho1"                          
##  [53] "MCMT_rho2"                          
##  [54] "TD"                                 
##  [55] "TD_rho1"                            
##  [56] "TD_rho2"                            
##  [57] "MAP"                                
##  [58] "MAP_rho1"                           
##  [59] "MAP_rho2"                           
##  [60] "MSP"                                
##  [61] "MSP_rho1"                           
##  [62] "MSP_rho2"                           
##  [63] "AHM"                                
##  [64] "AHM_rho1"                           
##  [65] "AHM_rho2"                           
##  [66] "SHM"                                
##  [67] "SHM_rho1"                           
##  [68] "SHM_rho2"                           
##  [69] "DD_0"                               
##  [70] "DD_0_rho1"                          
##  [71] "DD_0_rho2"                          
##  [72] "DD5"                                
##  [73] "DD5_rho1"                           
##  [74] "DD5_rho2"                           
##  [75] "NFFD"                               
##  [76] "NFFD_rho1"                          
##  [77] "NFFD_rho2"                          
##  [78] "bFFP"                               
##  [79] "bFFP_rho1"                          
##  [80] "bFFP_rho2"                          
##  [81] "eFFP"                               
##  [82] "eFFP_rho1"                          
##  [83] "eFFP_rho2"                          
##  [84] "FFP"                                
##  [85] "FFP_rho1"                           
##  [86] "FFP_rho2"                           
##  [87] "PAS"                                
##  [88] "PAS_rho1"                           
##  [89] "PAS_rho2"                           
##  [90] "EMT"                                
##  [91] "EMT_rho1"                           
##  [92] "EMT_rho2"                           
##  [93] "EXT"                                
##  [94] "EXT_rho1"                           
##  [95] "EXT_rho2"                           
##  [96] "Eref"                               
##  [97] "Eref_rho1"                          
##  [98] "Eref_rho2"                          
##  [99] "CMD"                                
## [100] "CMD_rho1"                           
## [101] "CMD_rho2"                           
## [102] "bayenv_PCA1"                        
## [103] "bayenv_PCA1_rho1"                   
## [104] "bayenv_PCA1_rho2"                   
## [105] "bayenv_PCA2"                        
## [106] "bayenv_PCA2_rho1"                   
## [107] "bayenv_PCA2_rho2"                   
## [108] "bayenv_PCA3"                        
## [109] "bayenv_PCA3_rho1"                   
## [110] "bayenv_PCA3_rho2"                   
## [111] "bayenv_PCA4"                        
## [112] "bayenv_PCA4_rho1"                   
## [113] "bayenv_PCA4_rho2"                   
## [114] "bayenv_PCA5"                        
## [115] "bayenv_PCA5_rho1"                   
## [116] "bayenv_PCA5_rho2"                   
## [117] "bayenv_PCA6"                        
## [118] "bayenv_PCA6_rho1"                   
## [119] "bayenv_PCA6_rho2"                   
## [120] "bayenv_PCA7"                        
## [121] "bayenv_PCA7_rho1"                   
## [122] "bayenv_PCA7_rho2"                   
## [123] "bayenv_PCA8"                        
## [124] "bayenv_PCA8_rho1"                   
## [125] "bayenv_PCA8_rho2"                   
## [126] "bayenv_PCA9"                        
## [127] "bayenv_PCA9_rho1"                   
## [128] "bayenv_PCA9_rho2"                   
## [129] "bayenv_PCA10"                       
## [130] "bayenv_PCA10_rho1"                  
## [131] "bayenv_PCA10_rho2"                  
## [132] "bayenv_PCA11"                       
## [133] "bayenv_PCA11_rho1"                  
## [134] "bayenv_PCA11_rho2"                  
## [135] "bayenv_PCA12"                       
## [136] "bayenv_PCA12_rho1"                  
## [137] "bayenv_PCA12_rho2"                  
## [138] "bayenv_PCA13"                       
## [139] "bayenv_PCA13_rho1"                  
## [140] "bayenv_PCA13_rho2"                  
## [141] "bayenv_PCA14"                       
## [142] "bayenv_PCA14_rho1"                  
## [143] "bayenv_PCA14_rho2"                  
## [144] "bayenv_PCA15"                       
## [145] "bayenv_PCA15_rho1"                  
## [146] "bayenv_PCA15_rho2"                  
## [147] "bayenv_PCA16"                       
## [148] "bayenv_PCA16_rho1"                  
## [149] "bayenv_PCA16_rho2"                  
## [150] "bayenv_PCA17"                       
## [151] "bayenv_PCA17_rho1"                  
## [152] "bayenv_PCA17_rho2"                  
## [153] "bayenv_PCA18"                       
## [154] "bayenv_PCA18_rho1"                  
## [155] "bayenv_PCA18_rho2"                  
## [156] "bayenv_PCA19"                       
## [157] "bayenv_PCA19_rho1"                  
## [158] "bayenv_PCA19_rho2"                  
## [159] "bayenv_PCA20"                       
## [160] "bayenv_PCA20_rho1"                  
## [161] "bayenv_PCA20_rho2"                  
## [162] "bayenv_PCA21"                       
## [163] "bayenv_PCA21_rho1"                  
## [164] "bayenv_PCA21_rho2"                  
## [165] "bayenv_PCA22"                       
## [166] "bayenv_PCA22_rho1"                  
## [167] "bayenv_PCA22_rho2"                  
## [168] "Budset_p"                           
## [169] "Budset_snp_effect"                  
## [170] "Budbreak_p"                         
## [171] "Budbreak_snp_effect"                
## [172] "Height_season_1_p"                  
## [173] "Height_season_1_snp_effect"         
## [174] "Height_season_2_p"                  
## [175] "Height_season_2_snp_effect"         
## [176] "Diameter_p"                         
## [177] "Diameter_snp_effect"                
## [178] "Shoot_weight_p"                     
## [179] "Shoot_weight_snp_effect"            
## [180] "Root_weight_p"                      
## [181] "Root_weight_snp_effect"             
## [182] "Max_growth_rate_p"                  
## [183] "Max_growth_rate_snp_effect"         
## [184] "Linear_growth_days_p"               
## [185] "Linear_growth_days_snp_effect"      
## [186] "X5_growth_complete_days_p"          
## [187] "X5_growth_complete_days_snp_effect" 
## [188] "X95_growth_complete_days_p"         
## [189] "X95_growth_complete_days_snp_effect"
## [190] "X5_95_growth_days_p"                
## [191] "X5_95_growth_days_snp_effect"       
## [192] "Fall_cold_injury_p"                 
## [193] "Fall_cold_injury_snp_effect"        
## [194] "Winter_cold_injury_p"               
## [195] "Winter_cold_injury_snp_effect"      
## [196] "Spring_cold_injury_p"               
## [197] "Spring_cold_injury_snp_effect"      
## [198] "root_wt__shoot_wt_p"                
## [199] "root_wt__shoot_wt_snp_effect"       
## [200] "root_wt__shoot_wt_p_1"              
## [201] "root_wt__shoot_wt_snp_effect_1"     
## [202] "gcta_PC1_p"                         
## [203] "gcta_PC1_snp_effect"                
## [204] "gcta_PC2_p"                         
## [205] "gcta_PC2_snp_effect"                
## [206] "gcta_PC3_p"                         
## [207] "gcta_PC3_snp_effect"                
## [208] "gcta_PC4_p"                         
## [209] "gcta_PC4_snp_effect"                
## [210] "gcta_PC5_p"                         
## [211] "gcta_PC5_snp_effect"                
## [212] "gcta_PC6_p"                         
## [213] "gcta_PC6_snp_effect"                
## [214] "gcta_PC7_p"                         
## [215] "gcta_PC7_snp_effect"                
## [216] "gcta_PC8_p"                         
## [217] "gcta_PC8_snp_effect"                
## [218] "gcta_PC9_p"                         
## [219] "gcta_PC9_snp_effect"                
## [220] "gcta_PC10_p"                        
## [221] "gcta_PC10_snp_effect"               
## [222] "gcta_PC11_p"                        
## [223] "gcta_PC11_snp_effect"               
## [224] "gcta_PC12_p"                        
## [225] "gcta_PC12_snp_effect"               
## [226] "gcta_PC13_p"                        
## [227] "gcta_PC13_snp_effect"               
## [228] "gcta_PC14_p"                        
## [229] "gcta_PC14_snp_effect"               
## [230] "gcta_PC15_p"                        
## [231] "gcta_PC15_snp_effect"               
## [232] "gcta_PC16_p"                        
## [233] "gcta_PC16_snp_effect"               
## [234] "gcta_PC17_p"                        
## [235] "gcta_PC17_snp_effect"               
## [236] "sequenced_contig_length"            
## [237] "affx_name"                          
## [238] "ax_name"                            
## [239] "LAT_raw_rho"                        
## [240] "LONG_raw_rho"                       
## [241] "ELEVATION_raw_rho"                  
## [242] "MAT_raw_rho"                        
## [243] "MWMT_raw_rho"                       
## [244] "MCMT_raw_rho"                       
## [245] "TD_raw_rho"                         
## [246] "MAP_raw_rho"                        
## [247] "MSP_raw_rho"                        
## [248] "AHM_raw_rho"                        
## [249] "SHM_raw_rho"                        
## [250] "DD_0_raw_rho"                       
## [251] "DD5_raw_rho"                        
## [252] "NFFD_raw_rho"                       
## [253] "bFFP_raw_rho"                       
## [254] "eFFP_raw_rho"                       
## [255] "FFP_raw_rho"                        
## [256] "PAS_raw_rho"                        
## [257] "EMT_raw_rho"                        
## [258] "EXT_raw_rho"                        
## [259] "Eref_raw_rho"                       
## [260] "CMD_raw_rho"                        
## [261] "LAT_raw_p"                          
## [262] "LONG_raw_p"                         
## [263] "ELEVATION_raw_p"                    
## [264] "MAT_raw_p"                          
## [265] "MWMT_raw_p"                         
## [266] "MCMT_raw_p"                         
## [267] "TD_raw_p"                           
## [268] "MAP_raw_p"                          
## [269] "MSP_raw_p"                          
## [270] "AHM_raw_p"                          
## [271] "SHM_raw_p"                          
## [272] "DD_0_raw_p"                         
## [273] "DD5_raw_p"                          
## [274] "NFFD_raw_p"                         
## [275] "bFFP_raw_p"                         
## [276] "eFFP_raw_p"                         
## [277] "FFP_raw_p"                          
## [278] "PAS_raw_p"                          
## [279] "EMT_raw_p"                          
## [280] "EXT_raw_p"                          
## [281] "Eref_raw_p"                         
## [282] "CMD_raw_p"                          
## [283] "Budset_p_raw_rho"                   
## [284] "Budbreak_p_raw_rho"                 
## [285] "Height_season_1_p_raw_rho"          
## [286] "Height_season_2_p_raw_rho"          
## [287] "Diameter_p_raw_rho"                 
## [288] "Shoot_weight_p_raw_rho"             
## [289] "Root_weight_p_raw_rho"              
## [290] "Max_growth_rate_p_raw_rho"          
## [291] "Linear_growth_days_p_raw_rho"       
## [292] "X5_growth_complete_days_p_raw_rho"  
## [293] "X95_growth_complete_days_p_raw_rho" 
## [294] "X5_95_growth_days_p_raw_rho"        
## [295] "Fall_cold_injury_p_raw_rho"         
## [296] "Winter_cold_injury_p_raw_rho"       
## [297] "Spring_cold_injury_p_raw_rho"       
## [298] "root_wt__shoot_wt_p_raw_rho"        
## [299] "root_wt__shoot_wt_p_1_raw_rho"      
## [300] "Budset_p_raw_p"                     
## [301] "Budbreak_p_raw_p"                   
## [302] "Height_season_1_p_raw_p"            
## [303] "Height_season_2_p_raw_p"            
## [304] "Diameter_p_raw_p"                   
## [305] "Shoot_weight_p_raw_p"               
## [306] "Root_weight_p_raw_p"                
## [307] "Max_growth_rate_p_raw_p"            
## [308] "Linear_growth_days_p_raw_p"         
## [309] "X5_growth_complete_days_p_raw_p"    
## [310] "X95_growth_complete_days_p_raw_p"   
## [311] "X5_95_growth_days_p_raw_p"          
## [312] "Fall_cold_injury_p_raw_p"           
## [313] "Winter_cold_injury_p_raw_p"         
## [314] "Spring_cold_injury_p_raw_p"         
## [315] "root_wt__shoot_wt_p_raw_p"          
## [316] "root_wt__shoot_wt_p_1_raw_p"
```



(Currently not implemented) Get the list of superoutliers.  There are 4 types for each species: based on transcript or gcontig, and based on raw or corrected for population structure.



Function used for plotting

```r
  plot_compare<- function(x,y, xlab, ylab, xlim, ylim){
  data1 <- cbind(x, y)
  data1b <- data1[complete.cases(data1),]
        binned <- bin2(data1b, 
                 matrix(c(-xlim,xlim,-ylim,ylim), 2,2, byrow=TRUE), 
                 nbin=c(100,100))
    binned$nc[binned$nc==0]=NA
    image(seq(-xlim,xlim,length.out = 100), seq(-ylim,ylim, length.out=100),binned$nc,
             xlab=xlab, ylab=ylab, add=FALSE, col=tim.colors(75))
  }

### Plotter function
plotter <- function(x, y, main, xlim, ylim, xlab, ylab){
  data1 <- cbind(x, y)
  data1b <- data1[complete.cases(data1),]
  score1 <- sum(data1b[,1]*data1b[,2]) #score 1 based on sum of products
  score1b <- round(score1/nrow(data1b)*10000,1)
  score2 <- round(sum((sign(data1b[,1])*sign(data1b[,2]))==1)/nrow(data1b),4) #score 2 based on number of points
    binned <- bin2(data1b, 
                 matrix(c(-xlim,xlim,-ylim,ylim), 2,2, byrow=TRUE), 
                 nbin=c(100,100))
    if(binned$nskip!=0){print("Error in binning, some values not counted. This affects score."); break}
    binned$nc[binned$nc==0]=NA
  
  slope <- lm(y~x)
  
    image.plot(seq(-xlim,xlim,length.out = 100), seq(-ylim,ylim, length.out=100),binned$nc,
             xlab=xlab, ylab=ylab, 
             main=paste(main, "\n, Score1 =",score1b, 
                        "Score2 =", score2,
                        "Slope =",round(slope$coef[2],5)))

  abline(slope)
  return(c(score1b, score2, slope$coef[2]))
  }
```

Next, loop through each phenotype-environment combo.  For each combo, flip GWAS to match the GEA and polarize the GEA scores so they are in matching direction to PE. (The 'flip' is going to be different for each SNP depending on how it was coded, while the 'polarization' multiplies all SNPs by -1 or by 1 depending on the slope between the phenotype and the environment.)

Then, standardize the GWAS so that each phenotype has a variance of 1 and all slopes will be on the same scale.  Plot the architecture and calculate a score:  

* score1 is `sum(GWAS_effect_pine_standardized * GEA_rho_pine_polarized)/n_snps`.  The null hypothesis for score1 is 0 (equal product of snp effects in all quadrants).  Positive scores indicate architectures that have more and/or larger effect SNPs for that phenotype in that environment.

* score2 is the proportion of points in the 1st and 3rd quadrants. Score2 varies between 0 and 1, with the null hypothesis of 0.5 (equal number of SNPs in all quadrants.)  Score2 > 0.5 indicates that the architecture is polarized in the same direction as the environment, while score2 < 0.5 indicates that the architecture is in the opposite direction as expected given the relationship between the phenotype and the environment.



```r
### Loop through PE file
for (i in 1:nrow(PE)){
### Match names
  PE_envi_name <- as.character(PE$Environment[i])
  GEA_envi_name <- as.character(PE_matchGEA$Environment_GEA[which(PE_matchGEA$Environment_PE==PE_envi_name)])
  
  PE_pheno_name <- as.character(PE$Phenotype[i])
  GWAS_pheno_name <- as.character(PE_matchGWAS$Phenotype_GWAS[which(PE_matchGWAS$Phenotype_PE==PE_pheno_name)])
  
  print(c(i, GEA_envi_name, GWAS_pheno_name))

###################
### Get GWAS relevant columns from large file (pine) 
###################
GWAS_effect_col_pine_gcta <- which(results_pine_head==gsub("_p","_snp_effect", GWAS_pheno_name))
  if (length(GWAS_effect_col_pine_gcta)!=1){print("ERROR: GWAS_effect_col_pine_corrected does not equal 1")}
  GWAS_effect_pine_gcta <- results_pine[,GWAS_effect_col_pine_gcta]

###################
### Get Bayenv relevant columns from large file (pine) 
###################
 GEA_effect_cols_pine_bayenv <- grep(paste("^",GEA_envi_name,"_rho", sep=""), results_pine_head)
  if (length(GEA_effect_cols_pine_bayenv)!=2){print("ERROR: GEA_effect_cols_pine_corrected does not equal 2")}
  GEA_rho_pine_bayenv <- rowMeans(results_pine[,GEA_effect_cols_pine_bayenv], na.rm=TRUE)

###################
### Get LFMM relevant columns from large file (pine)
###################
GEA_effect_col_pine_lfmm <- grep(paste("^",GEA_envi_name,"_z_lfmm_beagle", sep=""), results_pine_head)
  if (length(GEA_effect_col_pine_lfmm)!=1){print("ERROR: GEA_effect_col_pine_lfmm does not equal 1")}
  GEA_z_pine_lfmm <- results_pine[,GEA_effect_col_pine_lfmm]

###################
### Get RAW gwas relevant columns from large file (pine) 
###################
GWAS_pheno_name_raw <- GWAS_pheno_name  
# next lines fix a problem specific to pine df
  if(GWAS_pheno_name=="root_wt__shoot_wt_p"){
    GWAS_pheno_name_raw <- "root_wt_shoot_wt_p"}
  if(GWAS_pheno_name=="root_wt__shoot_wt_p_1"){
    GWAS_pheno_name_raw <- "root_wt_shoot_wt_p_1"} 
  GWAS_effect_col_pine_raw <- which(results_pine_head==paste(
    GWAS_pheno_name_raw, "_raw_rho", sep=""))
  if (length(GWAS_effect_col_pine_raw)!=1){
    print("ERROR: GWAS_effect_col_pine_raw does not equal 1")}
  GWAS_effect_pine_raw <- results_pine[,GWAS_effect_col_pine_raw]

  GEA_effect_cols_pine_raw <- which(results_pine_head==paste(
    GEA_envi_name,"_raw_rho", sep=""))
  if (length(GEA_effect_cols_pine_raw)!=1){
    print("ERROR: GEA_effect_cols_pine_raw does not equal 1")}
  GEA_rho_pine_raw <- results_pine[,GEA_effect_cols_pine_raw]

###################
### Plot PINE GEA raw vs corrected
###################
### This plot has the allele coding problem
  #plot_compare(GEA_rho_pine_bayenv, GEA_rho_pine_raw, 
  #             xlab = "GEA_rho_pine_bayenv", ylab = "GEA_rho_pine_raw",
  #             xlim = 1, ylim = 1)
  
  ### These two should match up
  plot_compare(GEA_rho_pine_bayenv, GEA_z_pine_lfmm, 
               xlab = "GEA_rho_pine_bayenv", ylab = "GEA_z_pine_lfmm",
               xlim = 1, ylim = 20)
  abline(lm(GEA_z_pine_lfmm~GEA_rho_pine_bayenv))
  mtext(paste(i, GEA_envi_name, GWAS_pheno_name))
  #plot_compare(abs(GEA_rho_pine_bayenv), abs(GEA_z_pine_lfmm), 
  #             xlab = "abs(GEA_rho_pine_bayenv)", 
  #             ylab = "abs(GEA_z_pine_lfmm)",
  #             xlim = 1, ylim = 3)
  #abline(lm(abs(GEA_z_pine_lfmm)~abs(GEA_rho_pine_bayenv)))
  
  ### This plot has the abs values, so you can see the positive relationship
  #plot_compare(abs(GEA_rho_pine_bayenv), abs(GEA_rho_pine_raw), 
  #             xlab = "abs(GEA_rho_pine_bayenv)", 
  #             ylab = "abs(GEA_rho_pine_raw)",
  #             xlim = 1, ylim = 1)
  #  abline(lm(abs(GEA_rho_pine_raw)~abs(GEA_rho_pine_bayenv)))
 
  #plot_compare(abs(GEA_z_pine_lfmm), abs(GEA_rho_pine_raw), 
  #             xlab = "abs(GEA_rho_pine_bayenv)", 
  #             ylab = "abs(GEA_rho_pine_raw)",
  #             xlim = 5, ylim = 1)
  #  abline(lm(abs(GEA_rho_pine_raw)~abs(GEA_z_pine_lfmm)))
  
    ### This plot should reflect the relationship between corrected_rho and raw_rho if we correctly flipped these variables
  GEA_rho_pine_raw_GOOD <- GEA_rho_pine_raw*results_pine$flip_pine_raw
  plot_compare(GEA_rho_pine_raw_GOOD, GEA_rho_pine_bayenv, 
               ylab = "GEA_rho_pine_bayenv", 
               xlab = "flip(GEA_rho_pine_raw)",
               xlim = 1, ylim = 1)
  abline(lm(GEA_rho_pine_bayenv~GEA_rho_pine_raw_GOOD))
  mtext(paste(i, GEA_envi_name, GWAS_pheno_name))

  GEA_z_pine_lfmm_STD <- GEA_z_pine_lfmm/sd(GEA_z_pine_lfmm, na.rm=TRUE)
   plot_compare(GEA_rho_pine_raw_GOOD, GEA_z_pine_lfmm,
               ylab = "GEA_z_pine_lfmm", 
               xlab = "flip(GEA_rho_pine_raw)",
               xlim = 1, ylim = 20)
  abline(lm(GEA_z_pine_lfmm~GEA_rho_pine_raw_GOOD))
    mtext(paste(i, GEA_envi_name, GWAS_pheno_name))

### Plot PINE GWAS raw vs. corrected
  ### This plot has the allele coding problem
  #plot_compare(GWAS_effect_pine_gcta, GWAS_effect_pine_raw, 
  #            xlab = "GWAS_effect_pine_gcta",
  #             ylab = "GWAS_effect_pine_raw",
  #             xlim = 5, ylim = 1)
  ### This plot has the abs values, so you can see the positive relationship
  #plot_compare(abs(GWAS_effect_pine_gcta), abs(GWAS_effect_pine_raw), 
  #             xlab = "abs(GWAS_effect_pine_corrected)",
  #             ylab = "abs(GWAS_effect_pine_raw)",
  #             xlim = 10, ylim = 1)
  #abline(lm(abs(GWAS_effect_pine_raw)~abs(GWAS_effect_pine_gcta)))
  ### This plot should reflect the relationship between corrected_gwas and raw_gwas if we correctly flipped these variables
  GWAS_gcta_pine_GOOD <- GWAS_effect_pine_gcta*results_pine$flip_pine_gwas*(-1)
  GWAS_gcta_pine_GOOD_STD <- GWAS_gcta_pine_GOOD/sd(GWAS_gcta_pine_GOOD, na.rm=TRUE)
    ### note multiplication by -1 above, think this is correct
  GWAS_raw_pine_GOOD <- GWAS_effect_pine_raw*results_pine$flip_pine_raw
  plot_compare(GWAS_raw_pine_GOOD,GWAS_gcta_pine_GOOD ,
               ylab = "flip(GWAS_effect_pine_gcta)",
               xlab = "flip(GWAS_effect_pine_raw)",
               xlim = 1, ylim = 20)
  abline(lm(GWAS_gcta_pine_GOOD~GWAS_raw_pine_GOOD))


###################
### Get GWAS relevant columns from large file (spruce) 
###################
GWAS_effect_col_spruce_gcta <- which(results_spruce_head==gsub("_p","_snp_effect", GWAS_pheno_name))
  if (length(GWAS_effect_col_spruce_gcta)!=1){print("ERROR: GWAS_effect_col_spruce_corrected does not equal 1")}
  GWAS_effect_spruce_gcta <- results_spruce[,GWAS_effect_col_spruce_gcta]

###################
### Get Bayenv relevant columns from large file (spruce) 
###################
 GEA_effect_cols_spruce_bayenv <- grep(paste("^",GEA_envi_name,"_rho", sep=""), results_spruce_head)
  if (length(GEA_effect_cols_spruce_bayenv)!=2){print("ERROR: GEA_effect_cols_spruce_corrected does not equal 2")}
  GEA_rho_spruce_bayenv <- rowMeans(results_spruce[,GEA_effect_cols_spruce_bayenv], na.rm=TRUE)

###################
### Get LFMM relevant columns from large file (spruce)
###################
#GEA_effect_col_spruce_lfmm <- grep(paste("^",GEA_envi_name,"_z_lfmm_beagle", sep=""), results_spruce_head)
#  if (length(GEA_effect_col_spruce_lfmm)!=1){print("ERROR: GEA_effect_col_spruce_lfmm does not equal 1")}
#  GEA_z_spruce_lfmm <- results_spruce[,GEA_effect_col_spruce_lfmm]

###################
### Get RAW GWAS and GEA relevant columns from large file (spruce)
###################

  GWAS_effect_col_spruce_raw <- which(results_spruce_head==paste(GWAS_pheno_name, "_raw_rho", sep=""))
  if (length(GWAS_effect_col_spruce_raw)!=1){print("ERROR: GWAS_effect_col_spruce_raw does not equal 1")}
  GWAS_effect_spruce_raw <- results_spruce[,GWAS_effect_col_spruce_raw]

  GEA_effect_cols_spruce_raw <- which(results_spruce_head==paste(GEA_envi_name,"_raw_rho", sep=""))
  if (length(GEA_effect_cols_spruce_raw)!=1){print("ERROR: GEA_effect_cols_spruce_raw does not equal 1")}
  GEA_rho_spruce_raw <- results_spruce[,GEA_effect_cols_spruce_raw]

###################
### Plot SPRUCE GEA raw vs corrected
###################
  ### Bayenv vs LFMM
#  plot_compare(GEA_rho_spruce_bayenv, GEA_z_spruce_lfmm, 
#               xlab = "GEA_rho_spruce_bayenv", ylab = "GEA_z_spruce_lfmm",
#               xlim = 1, ylim = 10)
#  abline(lm(GEA_z_spruce_lfmm~GEA_rho_spruce_bayenv))
#  mtext(paste(i, GEA_envi_name, GWAS_pheno_name))

  ### raw_GEA vs Bayenv
  GEA_rho_spruce_raw_GOOD <- GEA_rho_spruce_raw*results_spruce$flip_spruce_raw
  plot_compare(GEA_rho_spruce_raw_GOOD, GEA_rho_spruce_bayenv, 
               ylab = "GEA_rho_spruce_bayenv", 
               xlab = "flip(GEA_rho_spruce_raw)",
               xlim = 1, ylim = 1)
  abline(lm(GEA_rho_spruce_bayenv~GEA_rho_spruce_raw_GOOD))
  mtext(paste(i, GEA_envi_name, GWAS_pheno_name))

  ### raw_GEA vs lfmm
#  GEA_z_spruce_lfmm_STD <- GEA_z_spruce_lfmm/sd(GEA_z_spruce_lfmm, na.rm=TRUE)
#   plot_compare(GEA_rho_spruce_raw_GOOD, GEA_z_spruce_lfmm,
#               ylab = "GEA_z_spruce_lfmm", 
#               xlab = "flip(GEA_rho_spruce_raw)",
#               xlim = 1, ylim = 15)
#  abline(lm(GEA_z_spruce_lfmm~GEA_rho_spruce_raw_GOOD))
#    mtext(paste(i, GEA_envi_name, GWAS_pheno_name))

  ### raw_GWAS vs gcta
  GWAS_gcta_spruce_GOOD <- GWAS_effect_spruce_gcta*results_spruce$flip_spruce_gwas*(-1)
  GWAS_gcta_spruce_GOOD_STD <- GWAS_gcta_spruce_GOOD/sd(GWAS_gcta_spruce_GOOD, na.rm=TRUE)
    ### note multiplication by -1 above, think this is correct
  GWAS_raw_spruce_GOOD <- GWAS_effect_spruce_raw*results_spruce$flip_spruce_raw
  plot_compare(GWAS_raw_spruce_GOOD,GWAS_gcta_spruce_GOOD ,
               ylab = "flip(GWAS_effect_spruce_gcta)",
               xlab = "flip(GWAS_effect_spruce_raw)",
               xlim = 1, ylim = 20)
  abline(lm(GWAS_gcta_spruce_GOOD~GWAS_raw_spruce_GOOD))

###################
### Polarize GEA so sign aligns with PE and GWAS sign
###################
  PE_cor_pine <- PE$pine_correlation[i]
  PE_cor_spruce <- PE$spruce_correlation[i]

  polarize_pine <- sign(PE_cor_pine)
  if(polarize_pine==0){polarize_pine=1}
  polarize_spruce <- sign(PE_cor_spruce)
  if(polarize_spruce==0){polarize_spruce=1}

  GEA_rho_pine_bayenv_PLRZD <- GEA_rho_pine_bayenv * polarize_pine
  GEA_z_pine_lfmm_STD_PLRZD <- GEA_z_pine_lfmm_STD * polarize_pine
  GEA_rho_pine_raw_GOOD_PLRZD <- GEA_rho_pine_raw_GOOD * polarize_pine

  GEA_rho_spruce_bayenv_PLRZD <- GEA_rho_spruce_bayenv * polarize_spruce
 # GEA_z_spruce_lfmm_STD_PLRZD <- GEA_z_spruce_lfmm_STD * polarize_spruce
  GEA_rho_spruce_raw_GOOD_PLRZD <- GEA_rho_spruce_raw_GOOD * polarize_spruce

###############
### Plotting PEGA
###############
  png(paste("analysis/figure_byPheno/PEGA_" ,PE_pheno_name,"_",i, "_", PE_envi_name,".png",sep=""), width=10, height=15, units="in", res=300)
  par(mfrow=c(3,2), mar=c(4,4,4,1), oma=c(0,0,3,3), cex.main=0.8)

make_all <- function(){
  ### GCTA vs Bayenv PINE
  a.1p <- plotter(GWAS_gcta_pine_GOOD_STD, GEA_rho_pine_bayenv_PLRZD, 
                 main="All Corrected SNPs", 
                 xlim=40, ylim=1, 
                 xlab = "GWAS flipped & standardized PINE", 
                 ylab="Bayenv rho polarized PINE")
      PE$pine_score1_ALL_corrected_bayenv[i] <- a.1p[1]
      PE$pine_score2_ALL_corrected_bayenv[i] <- a.1p[2]
      PE$pine_slope_ALL_corrected_bayenv[i] <- a.1p[3]

  ### GCTA vs Bayenv SPRUCE 
  a.1s <-  plotter(GWAS_gcta_spruce_GOOD_STD, GEA_rho_spruce_bayenv_PLRZD, 
                 main="All Corrected SNPs", 
                 xlim=40, ylim=1, 
                 xlab = "GWAS flipped & standardized spruce", 
                 ylab="Bayenv rho polarized spruce")
      PE$spruce_score1_ALL_corrected_bayenv[i] <- a.1s[1]
      PE$spruce_score2_ALL_corrected_bayenv[i] <- a.1s[2]
      PE$spruce_slope_ALL_corrected_bayenv[i] <- a.1s[3]

  ### GCTA vs LFMM PINE
  a.2p <- plotter(GWAS_gcta_pine_GOOD_STD, GEA_z_pine_lfmm_STD_PLRZD, 
                 main="All Corrected SNPs", 
                 xlim=20, ylim=20, 
                 xlab = "GWAS flipped & standardized PINE", 
                 ylab="LFMM z standardized & polarized PINE")
      PE$pine_score1_ALL_corrected_lfmm[i] <- a.2p[1]
      PE$pine_score2_ALL_corrected_lfmm[i] <- a.2p[2]
      PE$pine_slope_ALL_corrected_lfmm[i] <- a.2p[3]

  ### GCTA vs LFMM SPRUCE
  #a.2s <- plotter(GWAS_gcta_spruce_GOOD_STD, GEA_z_spruce_lfmm_STD_PLRZD, 
  #               main="All Corrected SNPs", 
  #               xlim=20, ylim=20, 
  #               xlab = "GWAS flipped & standardized spruce", 
  #               ylab="LFMM z standardized & polarized spruce")
  a.2s <- c(NA, NA, NA)
  plot(0,0)
      PE$spruce_score1_ALL_corrected_lfmm[i] <- a.2s[1]
      PE$spruce_score2_ALL_corrected_lfmm[i] <- a.2s[2]
      PE$spruce_slope_ALL_corrected_lfmm[i] <- a.2s[3]

    ### GWAS RAW VS GEA RAW PINE
    a.3p <- plotter(GWAS_raw_pine_GOOD, GEA_rho_pine_raw_GOOD_PLRZD, 
          main="All SNPs RAW", xlim=1, ylim=1,
          xlab = "GWAS RAW PINE", 
          ylab="GEA RAW polarized PINE")
        PE$pine_score1_ALL_raw[i] <- a.3p[1]
        PE$pine_score2_ALL_raw[i] <- a.3p[2] 
        PE$pine_slope_ALL_raw[i] <- a.3p[3] 

    ### GWAS RAW VS GEA RAW spruce
    a.3s <- plotter(GWAS_raw_spruce_GOOD, GEA_rho_spruce_raw_GOOD_PLRZD, 
          main="All SNPs RAW", xlim=1, ylim=1,
          xlab = "GWAS RAW spruce", 
          ylab="GEA RAW polarized spruce")
        PE$spruce_score1_ALL_raw[i] <- a.3s[1]
        PE$spruce_score2_ALL_raw[i] <- a.3s[2] 
        PE$spruce_slope_ALL_raw[i] <- a.3s[3] 
    }#end make_all

  make_all()
  title("") #not sure why but needed to make mtext work
   mtext(paste(PE_pheno_name, PE_envi_name, "\nPine", PE$pine_correlation[i], 
               "\nSpruce", PE$spruce_correlation[i], sep=" "),
         side=3, outer=TRUE, line=-1)
  dev.off()

  make_all()
} #end loop
```

```
## [1] "1"                  "LAT"                "Fall_cold_injury_p"
```

![plot of chunk Loop through PE file](figure/Loop through PE file-1.png) ![plot of chunk Loop through PE file](figure/Loop through PE file-2.png) ![plot of chunk Loop through PE file](figure/Loop through PE file-3.png) ![plot of chunk Loop through PE file](figure/Loop through PE file-4.png) ![plot of chunk Loop through PE file](figure/Loop through PE file-5.png) 

```
## Error in dev.off(): QuartzBitmap_Output - unable to open file 'analysis/figure_byPheno/PEGA_ColdInjuryFallS02_Mean_1_Latitude.png'
```

![plot of chunk Loop through PE file](figure/Loop through PE file-6.png) 



Now that all the scores have been calculated, we can look at the relationship between the PE association and the underlying architecture.  I expect that higher/larger PEs will have higher scores.  The code plots and summarizes the results of a linear model 


```r
head(PE)
```

```
##     Environment                   Phenotype pine_correlation
## 190    Latitude      ColdInjuryFallS02_Mean            -0.42
## 36          CMD                  BudSet_day             0.40
## 275        MCMT      ColdInjuryFallS02_Mean             0.38
## 342         SHM                  BudSet_day             0.37
## 191    Latitude ColdInjuryMidWinterS02_Mean            -0.37
## 54          DD0      ColdInjuryFallS02_Mean            -0.37
##     spruce_correlation
## 190              -0.74
## 36                0.33
## 275               0.76
## 342               0.25
## 191              -0.62
## 54               -0.74
```

```r
  ylim=c(-70,70)

plotter_score <- function(x,y, xlab,ylab, ymin, ymax){
  plot(x, y,  ylim=c(ymin,ymax), 
       col=as.numeric(PE$Phenotype), pch=as.numeric(PE$Phenotype),
       xlab=xlab, ylab=ylab)
  m1 <- lm(y~x + PE$Environment + PE$Phenotype)
  print(summary(m1))
  abline(lm(y~x))
  abline(0,0, col="grey", lwd=2)
  abline(0.5,0, col="grey", lwd=2)
}

png(paste("analysis/Score1_Pine_Spruce_Envi.png"), width=10, height=10, units="in", res=300)
 par(mfrow=c(2,2), mar=c(4,4,4,1), oma=c(0,0,2,0))
  ### All results
  plotter_score(abs(PE$pine_correlation), PE$pine_score1_ALL_corrected_bayenv, 
                "PINE Abs(Pheno-Envi) Correlation", "PINE (corrected) Score 1 ALL SNPs", -40,40)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$spruce_correlation), PE$spruce_score1_ALL_corrected_bayenv, 
                "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE (corrected) Score 1 ALL SNPs", -40,40)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$pine_correlation), PE$pine_score1_ALL_raw, 
                "PINE Abs(Pheno-Envi) Correlation", "PINE (raw) Score 1 ALL SNPs", -30,30)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$spruce_correlation), PE$spruce_score1_ALL_raw, 
                "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE (raw) Score 1 ALL SNPs", -350, 350)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
dev.off()
```

```
## Error in dev.off(): QuartzBitmap_Output - unable to open file 'analysis/Score1_Pine_Spruce_Envi.png'
```

```r
png(paste("analysis/Score2_Pine_Spruce_Envi.png"), width=10, height=10, units="in", res=300)
 par(mfrow=c(2,2), mar=c(4,4,4,1), oma=c(0,0,2,0))
  plotter_score(abs(PE$pine_correlation), PE$pine_score2_ALL_corrected_bayenv, 
                "PINE Abs(Pheno-Envi) Correlation", "PINE (corrected) Proportion All SNPs", .4,.6)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$spruce_correlation), PE$spruce_score2_ALL_corrected_bayenv, 
                "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE (corrected) Proportion All SNPs", .4,.6)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$pine_correlation), PE$pine_score2_ALL_raw, 
                "PINE Abs(Pheno-Envi) Correlation", "PINE (raw) Proportion All SNPs", 0.1, 1)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
  plotter_score(abs(PE$spruce_correlation), PE$spruce_score2_ALL_raw, 
                "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE (raw) Proportion All SNPs", 0.1, 1)
```

```
## Error in model.frame.default(formula = y ~ x + PE$Environment + PE$Phenotype, : invalid type (NULL) for variable 'y'
```

```r
dev.off()
```

```
## Error in dev.off(): QuartzBitmap_Output - unable to open file 'analysis/Score2_Pine_Spruce_Envi.png'
```

```r
#   ### tcontig results
#   plotter_score(abs(PE$pine_correlation), PE$pine_score_super_t, 
#                 "PINE Abs(Pheno-Envi) Correlation", "PINE Score super_t")
#   plotter_score(abs(PE$spruce_correlation), PE$spruce_score_super_t, 
#                 "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE Score super_t")
#   ### tcontig results RAW
#   plotter_score(abs(PE$pine_correlation), PE$pine_score_super_t_raw, 
#                 "PINE Abs(Pheno-Envi) Correlation", "PINE Score super_t_raw")
#   plotter_score(abs(PE$spruce_correlation), PE$spruce_score_super_t_raw, 
#                 "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE Score super_t_raw")
# 
#   ### gcontig results
#   plotter_score(abs(PE$pine_correlation), PE$pine_score_super_g, 
#                 "PINE Abs(Pheno-Envi) Correlation", "PINE Score super_g")
#   plotter_score(abs(PE$spruce_correlation), PE$spruce_score_super_g, 
#                 "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE Score super_g")
#   ### gcontig results RAW
#   plotter_score(abs(PE$pine_correlation), PE$pine_score_super_g_raw, 
#                 "PINE Abs(Pheno-Envi) Correlation", "PINE Score super_g_raw")
#   plotter_score(abs(PE$spruce_correlation), PE$spruce_score_super_g_raw, 
#                 "SPRUCE Abs(Pheno-Envi) Correlation", "SPRUCE Score super_g_raw")
dev.off()
```

```
## quartz_off_screen 
##                 3
```

```
#plot allsnps vs super, expect super to be above 1:1 line

# color by environment is different than by phenotype

plot(PE$pine_score_super_t, PE$spruce_score_super_t,col=as.numeric(PE$Environment), pch=as.numeric(PE$Environment))
```
