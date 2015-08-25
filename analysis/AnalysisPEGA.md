---
title: "PEGA for Pine Analysis"
author: "Katie Lotterhos"
date: "August 24, 2015"
output: html_document
---

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
Note that the code runs fastest if you work through R in a terminal window and upload the "results" object first.  Once it is uploaded, when the markdown is created it won't try to load the "results" dataframe, which is very large.  To make this markdown file:
```
library(knitr)
library(markdown)
setwd("~/Desktop/CurrResearch/1-AdaptreeData/201509_PEGA/analysis")
knit("AnalysisPEGA.Rmd")  # produces the md file
markdownToHTML("AnalysisPEGA.md", "AnalysisPEGA.html")  # converts an md file to html
```
    

```r
setwd("~/Desktop/CurrResearch/1-AdaptreeData/201509_PEGA")
  if (!("hexbin" %in% installed.packages())){install.packages("hexbin")}
  library(hexbin)
library(gridExtra)
### load data
PE <- read.csv("data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2.csv", header=TRUE)
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
PE_matchGWAS <- read.csv("data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2_matchGWAS.csv", header=TRUE)
PE_matchGEA <- read.csv("data/EnvironmentPhenotypeCorrelations-PineSpruce-spearman_v2_matchGEA.csv", header=TRUE)

flip <- read.table("data/large_files/flip.cor", header=TRUE)
head(flip)
```

```
##   gcontig__gcontig_pos X1
## 1       C10572747__100 -1
## 2       C10572747__101 -1
## 3       C10572747__124 -1
## 4        C10572747__20  1
## 5        C10572747__25  1
## 6         C10572747__9 -1
```

```r
results_head <- scan("data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS", what = "character", nlines = 1)
results_head
```

```
##   [1] "gcontig__gcontig_pos"               
##   [2] "gcontig"                            
##   [3] "pos_gcontig"                        
##   [4] "X_annotation"                       
##   [5] "paralogy_rank"                      
##   [6] "tcontig"                            
##   [7] "pos_tcontig"                        
##   [8] "gmap_percentID"                     
##   [9] "offtarget_dist"                     
##  [10] "num_alleles"                        
##  [11] "num_ind"                            
##  [12] "num_comm_gtype"                     
##  [13] "major_allele"                       
##  [14] "major_num"                          
##  [15] "major_freq"                         
##  [16] "minor_allele"                       
##  [17] "minor_num"                          
##  [18] "minor_freq"                         
##  [19] "freq_hetzyg"                        
##  [20] "MAP_dir1"                           
##  [21] "MAP_dir2"                           
##  [22] "MAP_MgGm"                           
##  [23] "noncode_in_covar"                   
##  [24] "MQ"                                 
##  [25] "haplotype_score"                    
##  [26] "FS"                                 
##  [27] "MQRankSum"                          
##  [28] "ReadPosRankSum"                     
##  [29] "allele_balance"                     
##  [30] "window_stat"                        
##  [31] "scaffold"                           
##  [32] "scaf_pos"                           
##  [33] "xtx"                                
##  [34] "LAT"                                
##  [35] "LAT_rho1"                           
##  [36] "LAT_rho2"                           
##  [37] "LONG"                               
##  [38] "LONG_rho1"                          
##  [39] "LONG_rho2"                          
##  [40] "ELEVATION"                          
##  [41] "ELEVATION_rho1"                     
##  [42] "ELEVATION_rho2"                     
##  [43] "MAT"                                
##  [44] "MAT_rho1"                           
##  [45] "MAT_rho2"                           
##  [46] "MWMT"                               
##  [47] "MWMT_rho1"                          
##  [48] "MWMT_rho2"                          
##  [49] "MCMT"                               
##  [50] "MCMT_rho1"                          
##  [51] "MCMT_rho2"                          
##  [52] "TD"                                 
##  [53] "TD_rho1"                            
##  [54] "TD_rho2"                            
##  [55] "MAP"                                
##  [56] "MAP_rho1"                           
##  [57] "MAP_rho2"                           
##  [58] "MSP"                                
##  [59] "MSP_rho1"                           
##  [60] "MSP_rho2"                           
##  [61] "AHM"                                
##  [62] "AHM_rho1"                           
##  [63] "AHM_rho2"                           
##  [64] "SHM"                                
##  [65] "SHM_rho1"                           
##  [66] "SHM_rho2"                           
##  [67] "DD_0"                               
##  [68] "DD_0_rho1"                          
##  [69] "DD_0_rho2"                          
##  [70] "DD5"                                
##  [71] "DD5_rho1"                           
##  [72] "DD5_rho2"                           
##  [73] "NFFD"                               
##  [74] "NFFD_rho1"                          
##  [75] "NFFD_rho2"                          
##  [76] "bFFP"                               
##  [77] "bFFP_rho1"                          
##  [78] "bFFP_rho2"                          
##  [79] "eFFP"                               
##  [80] "eFFP_rho1"                          
##  [81] "eFFP_rho2"                          
##  [82] "FFP"                                
##  [83] "FFP_rho1"                           
##  [84] "FFP_rho2"                           
##  [85] "PAS"                                
##  [86] "PAS_rho1"                           
##  [87] "PAS_rho2"                           
##  [88] "EMT"                                
##  [89] "EMT_rho1"                           
##  [90] "EMT_rho2"                           
##  [91] "EXT"                                
##  [92] "EXT_rho1"                           
##  [93] "EXT_rho2"                           
##  [94] "Eref"                               
##  [95] "Eref_rho1"                          
##  [96] "Eref_rho2"                          
##  [97] "CMD"                                
##  [98] "CMD_rho1"                           
##  [99] "CMD_rho2"                           
## [100] "LAT_z_lfmm_beagle"                  
## [101] "LAT_q_lfmm_beagle"                  
## [102] "LONG_z_lfmm_beagle"                 
## [103] "LONG_q_lfmm_beagle"                 
## [104] "ELEVATION_z_lfmm_beagle"            
## [105] "ELEVATION_q_lfmm_beagle"            
## [106] "MAT_z_lfmm_beagle"                  
## [107] "MAT_q_lfmm_beagle"                  
## [108] "MWMT_z_lfmm_beagle"                 
## [109] "MWMT_q_lfmm_beagle"                 
## [110] "MCMT_z_lfmm_beagle"                 
## [111] "MCMT_q_lfmm_beagle"                 
## [112] "TD_z_lfmm_beagle"                   
## [113] "TD_q_lfmm_beagle"                   
## [114] "MAP_z_lfmm_beagle"                  
## [115] "MAP_q_lfmm_beagle"                  
## [116] "MSP_z_lfmm_beagle"                  
## [117] "MSP_q_lfmm_beagle"                  
## [118] "AHM_z_lfmm_beagle"                  
## [119] "AHM_q_lfmm_beagle"                  
## [120] "SHM_z_lfmm_beagle"                  
## [121] "SHM_q_lfmm_beagle"                  
## [122] "DD_0_z_lfmm_beagle"                 
## [123] "DD_0_q_lfmm_beagle"                 
## [124] "DD5_z_lfmm_beagle"                  
## [125] "DD5_q_lfmm_beagle"                  
## [126] "NFFD_z_lfmm_beagle"                 
## [127] "NFFD_q_lfmm_beagle"                 
## [128] "bFFP_z_lfmm_beagle"                 
## [129] "bFFP_q_lfmm_beagle"                 
## [130] "eFFP_z_lfmm_beagle"                 
## [131] "eFFP_q_lfmm_beagle"                 
## [132] "FFP_z_lfmm_beagle"                  
## [133] "FFP_q_lfmm_beagle"                  
## [134] "PAS_z_lfmm_beagle"                  
## [135] "PAS_q_lfmm_beagle"                  
## [136] "EMT_z_lfmm_beagle"                  
## [137] "EMT_q_lfmm_beagle"                  
## [138] "EXT_z_lfmm_beagle"                  
## [139] "EXT_q_lfmm_beagle"                  
## [140] "Eref_z_lfmm_beagle"                 
## [141] "Eref_q_lfmm_beagle"                 
## [142] "CMD_z_lfmm_beagle"                  
## [143] "CMD_q_lfmm_beagle"                  
## [144] "bayenv_PCA1"                        
## [145] "bayenv_PCA1_rho1"                   
## [146] "bayenv_PCA1_rho2"                   
## [147] "bayenv_PCA2"                        
## [148] "bayenv_PCA2_rho1"                   
## [149] "bayenv_PCA2_rho2"                   
## [150] "bayenv_PCA3"                        
## [151] "bayenv_PCA3_rho1"                   
## [152] "bayenv_PCA3_rho2"                   
## [153] "bayenv_PCA4"                        
## [154] "bayenv_PCA4_rho1"                   
## [155] "bayenv_PCA4_rho2"                   
## [156] "bayenv_PCA5"                        
## [157] "bayenv_PCA5_rho1"                   
## [158] "bayenv_PCA5_rho2"                   
## [159] "bayenv_PCA6"                        
## [160] "bayenv_PCA6_rho1"                   
## [161] "bayenv_PCA6_rho2"                   
## [162] "bayenv_PCA7"                        
## [163] "bayenv_PCA7_rho1"                   
## [164] "bayenv_PCA7_rho2"                   
## [165] "bayenv_PCA8"                        
## [166] "bayenv_PCA8_rho1"                   
## [167] "bayenv_PCA8_rho2"                   
## [168] "bayenv_PCA9"                        
## [169] "bayenv_PCA9_rho1"                   
## [170] "bayenv_PCA9_rho2"                   
## [171] "bayenv_PCA10"                       
## [172] "bayenv_PCA10_rho1"                  
## [173] "bayenv_PCA10_rho2"                  
## [174] "bayenv_PCA11"                       
## [175] "bayenv_PCA11_rho1"                  
## [176] "bayenv_PCA11_rho2"                  
## [177] "bayenv_PCA12"                       
## [178] "bayenv_PCA12_rho1"                  
## [179] "bayenv_PCA12_rho2"                  
## [180] "bayenv_PCA13"                       
## [181] "bayenv_PCA13_rho1"                  
## [182] "bayenv_PCA13_rho2"                  
## [183] "bayenv_PCA14"                       
## [184] "bayenv_PCA14_rho1"                  
## [185] "bayenv_PCA14_rho2"                  
## [186] "bayenv_PCA15"                       
## [187] "bayenv_PCA15_rho1"                  
## [188] "bayenv_PCA15_rho2"                  
## [189] "bayenv_PCA16"                       
## [190] "bayenv_PCA16_rho1"                  
## [191] "bayenv_PCA16_rho2"                  
## [192] "bayenv_PCA17"                       
## [193] "bayenv_PCA17_rho1"                  
## [194] "bayenv_PCA17_rho2"                  
## [195] "bayenv_PCA18"                       
## [196] "bayenv_PCA18_rho1"                  
## [197] "bayenv_PCA18_rho2"                  
## [198] "bayenv_PCA19"                       
## [199] "bayenv_PCA19_rho1"                  
## [200] "bayenv_PCA19_rho2"                  
## [201] "bayenv_PCA20"                       
## [202] "bayenv_PCA20_rho1"                  
## [203] "bayenv_PCA20_rho2"                  
## [204] "bayenv_PCA21"                       
## [205] "bayenv_PCA21_rho1"                  
## [206] "bayenv_PCA21_rho2"                  
## [207] "bayenv_PCA22"                       
## [208] "bayenv_PCA22_rho1"                  
## [209] "bayenv_PCA22_rho2"                  
## [210] "Budset_p"                           
## [211] "Budset_snp_effect"                  
## [212] "Budbreak_p"                         
## [213] "Budbreak_snp_effect"                
## [214] "Height_season_1_p"                  
## [215] "Height_season_1_snp_effect"         
## [216] "Height_season_2_p"                  
## [217] "Height_season_2_snp_effect"         
## [218] "Diameter_p"                         
## [219] "Diameter_snp_effect"                
## [220] "Shoot_weight_p"                     
## [221] "Shoot_weight_snp_effect"            
## [222] "Root_weight_p"                      
## [223] "Root_weight_snp_effect"             
## [224] "Max_growth_rate_p"                  
## [225] "Max_growth_rate_snp_effect"         
## [226] "Linear_growth_days_p"               
## [227] "Linear_growth_days_snp_effect"      
## [228] "X5_growth_complete_days_p"          
## [229] "X5_growth_complete_days_snp_effect" 
## [230] "X95_growth_complete_days_p"         
## [231] "X95_growth_complete_days_snp_effect"
## [232] "X5_95_growth_days_p"                
## [233] "X5_95_growth_days_snp_effect"       
## [234] "Fall_cold_injury_p"                 
## [235] "Fall_cold_injury_snp_effect"        
## [236] "Winter_cold_injury_p"               
## [237] "Winter_cold_injury_snp_effect"      
## [238] "Spring_cold_injury_p"               
## [239] "Spring_cold_injury_snp_effect"      
## [240] "root_wt__shoot_wt_p"                
## [241] "root_wt__shoot_wt_snp_effect"       
## [242] "root_wt__shoot_wt_p_1"              
## [243] "root_wt__shoot_wt_snp_effect_1"     
## [244] "gcta_PC1_p"                         
## [245] "gcta_PC1_snp_effect"                
## [246] "gcta_PC2_p"                         
## [247] "gcta_PC2_snp_effect"                
## [248] "gcta_PC3_p"                         
## [249] "gcta_PC3_snp_effect"                
## [250] "gcta_PC4_p"                         
## [251] "gcta_PC4_snp_effect"                
## [252] "gcta_PC5_p"                         
## [253] "gcta_PC5_snp_effect"                
## [254] "gcta_PC6_p"                         
## [255] "gcta_PC6_snp_effect"                
## [256] "gcta_PC7_p"                         
## [257] "gcta_PC7_snp_effect"                
## [258] "gcta_PC8_p"                         
## [259] "gcta_PC8_snp_effect"                
## [260] "gcta_PC9_p"                         
## [261] "gcta_PC9_snp_effect"                
## [262] "gcta_PC10_p"                        
## [263] "gcta_PC10_snp_effect"               
## [264] "gcta_PC11_p"                        
## [265] "gcta_PC11_snp_effect"               
## [266] "gcta_PC12_p"                        
## [267] "gcta_PC12_snp_effect"               
## [268] "gcta_PC13_p"                        
## [269] "gcta_PC13_snp_effect"               
## [270] "gcta_PC14_p"                        
## [271] "gcta_PC14_snp_effect"               
## [272] "gcta_PC15_p"                        
## [273] "gcta_PC15_snp_effect"               
## [274] "gcta_PC16_p"                        
## [275] "gcta_PC16_snp_effect"               
## [276] "category_snp_array"                 
## [277] "sequenced_contig_length"            
## [278] "affx_name"                          
## [279] "ax_name"
```

```r
#### read in results if it doesn't exist
if (!("results" %in% ls())){
results <- read.table("data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS", header=TRUE, comment.char=" ")
  }
dim(results)
```

```
## [1] 1251074     279
```

```r
### Loop through PE file

for (i in 1:5){#nrow(PE)){
  print(i)
### Match names
  PE_envi_name <- as.character(PE$Environment[i])
  GEA_envi_name <- as.character(PE_matchGEA$Environment_GEA[which(PE_matchGEA$Environment_PE==PE_envi_name)])
  
  PE_pheno_name <- as.character(PE$Phenotype[i])
  GWAS_pheno_name <- as.character(PE_matchGWAS$Phenotype_GWAS[which(PE_matchGWAS$Phenotype_PE==PE_pheno_name)])

### Get relevant columns from large file
  
  GWAS_effect_col <- which(results_head==gsub("p","snp_effect", GWAS_pheno_name))
  GWAS_effect <- results[,GWAS_effect_col]

  GEA_effect_cols <- grep(paste(GEA_envi_name,"_rho", sep=""), results_head)
  GEA_rho <- rowMeans(results[,GEA_effect_cols], na.rm=TRUE)

### check that "flip"" order matches "results" order
if (!(identical(flip$gcontig__gcontig_pos, results$gcontig__gcontig_pos))){
  print("ERROR!  SNP order in flip file does not match results file")
}

### Flip GEA_rho so sign matches GWAS_effect
  GEA_rho_flip <- GEA_rho * flip$X1

### Polarize GEA so sign aligns with PE and GWAS sign
  PE_cor <- PE$pine_correlation[i]

  GEA_rho_polarized <- GEA_rho_flip * sign(PE_cor) * sign(GWAS_effect)

#  plot(0,0, xlim=c(-3,3), ylim=c(-0.5,0.5), col=0, main=paste(PE_envi_name, PE_pheno_name, PE$pine_correlation[i]))
  #plot(hexbin(GWAS_effect, GEA_rho_polarized,xbnds=c(-5,5), ybnds=c(-1,1), xbins=100), add=TRUE)
       
hxb<-hexbin(GWAS_effect, GEA_rho_polarized,xbnds=c(-5,5), ybnds=c(-1,1), xbins=100)
#hxb@xcm    #gives the x co-ordinates of each hex tile
#hxb@ycm    #gives the y co-ordinates of each hex tile
#hxb@count  #gives the cell size for each hex tile
#points(x=hxb@xcm, y=hxb@ycm, col=grey(hxb@count*0.01))

  png(paste("analysis/PEGA_",PE_envi_name,"_",PE_pheno_name,"_",i, ".png",sep=""), width=6, height=5, units="in", res=300)
  par(mar=c(1,1,1,1))
  a<- hexbinplot(GEA_rho_polarized~GWAS_effect, xbnds=c(-5,5), ybnds=c(-1,1), xbins=100, main=paste(PE_envi_name, PE_pheno_name, PE$pine_correlation[i]) , xlim=c(-3,3), ylim=c(-0.5,0.5), aspect=1)
  plot(a)
  #do.call(grid.arrange, c(a, ncol=1))
  dev.off()

  plot(a)
} #end loop
```

```
## [1] 1
```

```
## [1] 2
```

```
## [1] 3
```

```
## [1] 4
```

```
## [1] 5
```
