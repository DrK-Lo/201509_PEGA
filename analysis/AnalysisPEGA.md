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
Note that the code runs fastest if you work through R in a terminal window and upload the "results_pine" object first.  Once it is uploaded, when the markdown is created it won't try to load the "results_pine" dataframe, which is very large.  To make this markdown file:
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
  if (!("gridExtra" %in% installed.packages())){install.packages("ash")}
  library(gridExtra)
  if (!("ash" %in% installed.packages())){install.packages("ash")}
  library(ash)

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

flip_pine <- read.table("data/large_files/flip_pine.cor", header=TRUE)
head(flip_pine)
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
#### read in results_pine if it doesn't exist
if (!("results_pine" %in% ls())){
results_pine <- read.table("data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS", header=TRUE, comment.char=" ")
  }
dim(results_pine)
```

```
## [1] 1251074     279
```

```r
### Upload superoutliers file
  pine_raw_tcontig <- read.table("data/large_files/super_outliers_per_variable_pine_both_methods_100x_raw.txt", header=TRUE)
  pine_tcontig <- read.table("data/large_files/super_outliers_per_variable_pine_both_methods_100x.txt", header=TRUE)
  pine_raw_gcontig <- read.table("data/large_files/super_outliers_per_variable_pine_both_methods_gcontigs_100x_raw.txt", header=TRUE)
  pine_gcontig <- read.table("data/large_files/super_outliers_per_variable_pine_both_methods_gcontigs_100x.txt", header=TRUE)

  head(pine_tcontig)
```

```
##                 pine snp_count outlier_count      ratio
## 1      comp0_c0_seq1        14             1 0.07142857
## 2 comp100159_c0_seq1         3             2 0.66666667
## 3  comp10016_c0_seq1        56             2 0.03571429
## 4  comp10018_c0_seq1        53             2 0.03773585
## 5   comp1001_c0_seq1        42             1 0.02380952
## 6  comp10020_c0_seq1        20             1 0.05000000
##   super_outlier_score            status test_name
## 1         -0.09765642 non-super-outlier       LAT
## 2          0.42145652 non-super-outlier       LAT
## 3         -0.18458752 non-super-outlier       LAT
## 4         -0.16906323 non-super-outlier       LAT
## 5         -0.23498452 non-super-outlier       LAT
## 6         -0.12775533 non-super-outlier       LAT
```

```r
  head(pine_gcontig)
```

```
##                pine snp_count outlier_count      ratio super_outlier_score
## 1         C24432054        10             2 0.20000000          0.13929921
## 2         C26044696         1             1 1.00000000         -0.02491853
## 3         C26870990        14             1 0.07142857         -0.09595732
## 4         C27752974         1             1 1.00000000         -0.02491853
## 5  C28162626_93_546         3             1 0.33333333         -0.03685996
## 6 C28932838_119_624        14             2 0.14285714          0.09011353
##              status test_name
## 1 non-super-outlier       LAT
## 2     super_outlier       LAT
## 3 non-super-outlier       LAT
## 4     super_outlier       LAT
## 5 non-super-outlier       LAT
## 6 non-super-outlier       LAT
```

```r
### Get list of superoutliers
  pine_super_t <- as.character(pine_tcontig$pine[pine_tcontig$super_outlier_score > 1])
  pine_indexes_super_t <- which(as.character(results_pine$tcontig) %in% pine_super_t)

  pine_super_t_raw <- as.character(pine_raw_tcontig$pine[pine_raw_tcontig$super_outlier_score > 1])
  pine_indexes_super_t_raw <- which(as.character(results_pine$tcontig) %in% pine_super_t_raw)

  pine_super_g <- as.character(pine_gcontig$pine[pine_gcontig$super_outlier_score > 1])
  pine_indexes_super_g <- which(as.character(results_pine$gcontig) %in% pine_super_g)

  pine_super_g_raw <- as.character(pine_raw_gcontig$pine[pine_raw_gcontig$super_outlier_score > 1])
  pine_indexes_super_g_raw <- which(as.character(results_pine$gcontig) %in% pine_super_g_raw)

  PE$pine_score_ALL <- PE$pine_score_super_t <- PE$pine_score_super_t_raw <- PE$pine_score_super_g <- PE$pine_score_super_g_raw <- NA

### Loop through PE file

for (i in 1:nrow(PE)){
### Match names
  PE_envi_name <- as.character(PE$Environment[i])
  GEA_envi_name <- as.character(PE_matchGEA$Environment_GEA[which(PE_matchGEA$Environment_PE==PE_envi_name)])
  
  PE_pheno_name <- as.character(PE$Phenotype[i])
  GWAS_pheno_name <- as.character(PE_matchGWAS$Phenotype_GWAS[which(PE_matchGWAS$Phenotype_PE==PE_pheno_name)])
  
  print(c(i, GEA_envi_name, GWAS_pheno_name))

### Get relevant columns from large file
  results_pine_head <- names(results_pine)
  GWAS_effect_col <- which(results_pine_head==gsub("_p","_snp_effect", GWAS_pheno_name))
  if (length(GWAS_effect_col)!=1){print("ERROR: GWAS_effect_col does not equal 1")}
  GWAS_effect <- results_pine[,GWAS_effect_col]

  GEA_effect_cols <- grep(paste("^",GEA_envi_name,"_rho", sep=""), results_pine_head)
  if (length(GEA_effect_cols)!=2){print("ERROR: GEA_effect_cols does not equal 2")}
  GEA_rho <- rowMeans(results_pine[,GEA_effect_cols], na.rm=TRUE)

### check that "flip_pine"" order matches "results_pine" order
if (!(identical(flip_pine$gcontig__gcontig_pos, results_pine$gcontig__gcontig_pos))){
  print("ERROR!  SNP order in flip_pine file does not match results_pine file")
}

### flip_pine GEA_rho so sign matches GWAS_effect
  GEA_rho_flip_pine <- GEA_rho * flip_pine$X1

### Polarize GEA so sign aligns with PE and GWAS sign
  PE_cor <- PE$pine_correlation[i]

  GWAS_sign <- sign(GWAS_effect)
  GWAS_sign[GWAS_sign==0]=1

  if (sum(GWAS_sign==0, na.rm=TRUE)>0){print("ERROR: Some GWAS have 0 sign for polarization, will skew results_pine")} #fixed this with code below

  if(PE_cor != 0){
  GEA_rho_polarized <- GEA_rho_flip_pine * sign(PE_cor) * GWAS_sign
  }else{
  GEA_rho_polarized <- GEA_rho_flip_pine * GWAS_sign  
  }

#  plot(0,0, xlim=c(-3,3), ylim=c(-0.5,0.5), col=0, main=paste(PE_envi_name, PE_pheno_name, PE$pine_correlation[i]))
  #plot(hexbin(GWAS_effect, GEA_rho_polarized,xbnds=c(-5,5), ybnds=c(-1,1), xbins=100), add=TRUE)


### Plotting 


GWAS_effect_standardized <- GWAS_effect/sd(GWAS_effect, na.rm=TRUE)

plotter <- function(GWAS_effect, GEA_rho_polarized, main, is.standard){
  data1 <- cbind(GWAS_effect, GEA_rho_polarized)
  data1b <- data1[complete.cases(data1),]
  score1 <- sum(data1b[,1]*data1b[,2])
  score <- score1/nrow(data1b)*10000
  if(is.standard==FALSE){
    binned <- bin2(data1b, 
                 matrix(c(-3,3,-0.5,0.5), 2,2, byrow=TRUE), 
                 nbin=c(100,100))
    binned$nc[binned$nc==0]=NA
    image.plot(seq(-3,3,length.out = 100), seq(-0.5, 0.5, length.out=100),binned$nc,
             xlab="GWAS effect", ylab="GEA rho (polarized)", main=main)
    #text(2.9,0.45, round(score,1))
  }
  if(is.standard==TRUE){
    binned <- bin2(data1b, 
                 matrix(c(-10,10,-0.5,0.5), 2,2, byrow=TRUE), 
                 nbin=c(100,100))
    binned$nc[binned$nc==0]=NA
    image.plot(seq(-10,10,length.out = 100), seq(-0.5, 0.5, length.out=100),binned$nc,
             xlab="GWAS effect standardized", ylab="GEA rho (polarized)", main=paste(main, ", Score = ",round(score,1)))
  }
  return(score)
  }

  png(paste("analysis/figure_byPheno/PEGA_" ,PE_pheno_name,"_",i, "_", PE_envi_name,".png",sep=""), width=5, height=18, units="in", res=300)
  par(mfrow=c(6,1), mar=c(4,4,4,1), oma=c(0,0,2,0))
  plotter(GWAS_effect, GEA_rho_polarized, main="All SNPs", is.standard=FALSE)
  PE$pine_score_ALL[i]<-plotter(GWAS_effect_standardized, GEA_rho_polarized, main="All SNPs Standardized", is.standard=TRUE)
  PE$pine_score_super_t[i]<-plotter(GWAS_effect_standardized[pine_indexes_super_t], GEA_rho_polarized[pine_indexes_super_t], main="Super tcontigs", is.standard=TRUE)
  PE$pine_score_super_t_raw[i]<-plotter(GWAS_effect_standardized[pine_indexes_super_t_raw], GEA_rho_polarized[pine_indexes_super_t_raw], main="Super tcontigs RAW", is.standard=TRUE)
  PE$pine_score_super_g[i]<-plotter(GWAS_effect_standardized[pine_indexes_super_g], GEA_rho_polarized[pine_indexes_super_g], main="Super gcontigs", is.standard=TRUE)
  PE$pine_score_super_g_raw[i]<-plotter(GWAS_effect_standardized[pine_indexes_super_g_raw], GEA_rho_polarized[pine_indexes_super_g_raw], main="Super gcontigs RAW", is.standard=TRUE)
  title("yo")
  mtext(paste(PE_pheno_name, PE_envi_name, "\nPine", PE$pine_correlation[i], sep=" "),
        side=3, outer=TRUE, line=-1)
  dev.off()

} #end loop
```

```
## [1] "1"                  "LAT"                "Fall_cold_injury_p"
```

```
## [1] "2"        "CMD"      "Budset_p"
```

```
## [1] "3"                  "MCMT"               "Fall_cold_injury_p"
```

```
## [1] "4"        "SHM"      "Budset_p"
```

```
## [1] "5"                    "LAT"                  "Winter_cold_injury_p"
```

```
## [1] "6"                  "DD_0"               "Fall_cold_injury_p"
```

```
## [1] "7"        "EMT"      "Budset_p"
```

```
## [1] "8"                  "EMT"                "Fall_cold_injury_p"
```

```
## [1] "9"                    "MCMT"                 "Winter_cold_injury_p"
```

```
## [1] "10"       "MSP"      "Budset_p"
```

```
## [1] "11"                "ELEVATION"         "Height_season_2_p"
```

```
## [1] "12"       "EXT"      "Budset_p"
```

```
## [1] "13"                "DD5"               "Height_season_2_p"
```

```
## [1] "14"       "DD_0"     "Budset_p"
```

```
## [1] "15"                 "TD"                 "Fall_cold_injury_p"
```

```
## [1] "16"                 "MAT"                "Fall_cold_injury_p"
```

```
## [1] "17"                "EXT"               "Height_season_2_p"
```

```
## [1] "18"                "Eref"              "Height_season_1_p"
```

```
## [1] "19"                   "EMT"                  "Winter_cold_injury_p"
```

```
## [1] "20"                "MWMT"              "Height_season_2_p"
```

```
## [1] "21"                "DD5"               "Height_season_1_p"
```

```
## [1] "22"                "NFFD"              "Height_season_2_p"
```

```
## [1] "23"       "MCMT"     "Budset_p"
```

```
## [1] "24"                   "DD_0"                 "Winter_cold_injury_p"
```

```
## [1] "25"                "MWMT"              "Height_season_1_p"
```

```
## [1] "26"       "MAT"      "Budset_p"
```

```
## [1] "27"                "EXT"               "Height_season_1_p"
```

```
## [1] "28"                   "TD"                   "Winter_cold_injury_p"
```

```
## [1] "29"                "CMD"               "Height_season_2_p"
```

```
## [1] "30"       "NFFD"     "Budset_p"
```

```
## [1] "31"                "MAT"               "Height_season_1_p"
```

```
## [1] "32"                "Eref"              "Height_season_2_p"
```

```
## [1] "33"                "bFFP"              "Height_season_2_p"
```

```
## [1] "34"                "SHM"               "Height_season_2_p"
```

```
## [1] "35"                "NFFD"              "Height_season_1_p"
```

```
## [1] "36"                "MAT"               "Height_season_2_p"
```

```
## [1] "37"                   "MAP"                  "Winter_cold_injury_p"
```

```
## [1] "38"                "FFP"               "Height_season_2_p"
```

```
## [1] "39"                 "Eref"               "Fall_cold_injury_p"
```

```
## [1] "40"       "Eref"     "Budset_p"
```

```
## [1] "41"                "bFFP"              "Height_season_1_p"
```

```
## [1] "42"                 "MAP"                "Fall_cold_injury_p"
```

```
## [1] "43"        "ELEVATION" "Budset_p"
```

```
## [1] "44"                "FFP"               "Height_season_1_p"
```

```
## [1] "45"         "EXT"        "Diameter_p"
```

```
## [1] "46"             "ELEVATION"      "Shoot_weight_p"
```

```
## [1] "47"                   "eFFP"                 "Winter_cold_injury_p"
```

```
## [1] "48"       "DD5"      "Budset_p"
```

```
## [1] "49"                   "PAS"                  "Spring_cold_injury_p"
```

```
## [1] "50"         "PAS"        "Budbreak_p"
```

```
## [1] "51"                         "NFFD"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "52"         "MWMT"       "Diameter_p"
```

```
## [1] "53"       "MWMT"     "Budset_p"
```

```
## [1] "54"                   "MAT"                  "Winter_cold_injury_p"
```

```
## [1] "55"            "EXT"           "Root_weight_p"
```

```
## [1] "56"                         "eFFP"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "57"                 "eFFP"               "Fall_cold_injury_p"
```

```
## [1] "58"                         "FFP"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "59"                         "EMT"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "60"         "ELEVATION"  "Diameter_p"
```

```
## [1] "61"         "DD5"        "Diameter_p"
```

```
## [1] "62"            "DD5"           "Root_weight_p"
```

```
## [1] "63"         "bFFP"       "Budbreak_p"
```

```
## [1] "64"                  "AHM"                 "root_wt__shoot_wt_p"
```

```
## [1] "65"                "AHM"               "Height_season_1_p"
```

```
## [1] "66"                   "PAS"                  "Winter_cold_injury_p"
```

```
## [1] "67"            "MWMT"          "Root_weight_p"
```

```
## [1] "68"       "LONG"     "Budset_p"
```

```
## [1] "69"                  "EXT"                 "root_wt__shoot_wt_p"
```

```
## [1] "70"            "Eref"          "Root_weight_p"
```

```
## [1] "71"       "eFFP"     "Budset_p"
```

```
## [1] "72"                "DD_0"              "Max_growth_rate_p"
```

```
## [1] "73"                "CMD"               "Height_season_1_p"
```

```
## [1] "74"                    "AHM"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "75"                "AHM"               "Height_season_2_p"
```

```
## [1] "76"            "AHM"           "Root_weight_p"
```

```
## [1] "77"                "MCMT"              "Max_growth_rate_p"
```

```
## [1] "78"                "MAT"               "Max_growth_rate_p"
```

```
## [1] "79"                "EMT"               "Max_growth_rate_p"
```

```
## [1] "80"                        "ELEVATION"                
## [3] "X5_growth_complete_days_p"
```

```
## [1] "81"            "ELEVATION"     "Root_weight_p"
```

```
## [1] "82"                   "ELEVATION"            "Spring_cold_injury_p"
```

```
## [1] "83"                "PAS"               "Height_season_1_p"
```

```
## [1] "84"                 "PAS"                "Fall_cold_injury_p"
```

```
## [1] "85"            "NFFD"          "Root_weight_p"
```

```
## [1] "86"                    "MAP"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "87"                   "LAT"                  "Spring_cold_injury_p"
```

```
## [1] "88"         "ELEVATION"  "Budbreak_p"
```

```
## [1] "89"                "eFFP"              "Height_season_2_p"
```

```
## [1] "90"                  "DD5"                 "root_wt__shoot_wt_p"
```

```
## [1] "91"                         "DD5"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "92"                "DD_0"              "Height_season_2_p"
```

```
## [1] "93"            "CMD"           "Root_weight_p"
```

```
## [1] "94"                 "NFFD"               "Fall_cold_injury_p"
```

```
## [1] "95"                  "MWMT"                "root_wt__shoot_wt_p"
```

```
## [1] "96"                         "MWMT"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "97"                         "MAT"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "98"                "MSP"               "Height_season_2_p"
```

```
## [1] "99"                  "Eref"                "root_wt__shoot_wt_p"
```

```
## [1] "100"        "Eref"       "Diameter_p"
```

```
## [1] "101"                 "ELEVATION"           "root_wt__shoot_wt_p"
```

```
## [1] "102"               "DD_0"              "Height_season_1_p"
```

```
## [1] "103"                 "CMD"                 "root_wt__shoot_wt_p"
```

```
## [1] "104"                        "bFFP"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "105"        "bFFP"       "Diameter_p"
```

```
## [1] "106"           "bFFP"          "Root_weight_p"
```

```
## [1] "107"               "SHM"               "Height_season_1_p"
```

```
## [1] "108"        "NFFD"       "Diameter_p"
```

```
## [1] "109"           "FFP"           "Root_weight_p"
```

```
## [1] "110"        "FFP"        "Budbreak_p"
```

```
## [1] "111"               "EMT"               "Height_season_2_p"
```

```
## [1] "112"        "DD5"        "Budbreak_p"
```

```
## [1] "113"                        "DD_0"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "114"        "CMD"        "Diameter_p"
```

```
## [1] "115"                       "bFFP"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "116"        "AHM"        "Diameter_p"
```

```
## [1] "117"                  "TD"                   "Spring_cold_injury_p"
```

```
## [1] "118"           "SHM"           "Root_weight_p"
```

```
## [1] "119"      "PAS"      "Budset_p"
```

```
## [1] "120"               "NFFD"              "Max_growth_rate_p"
```

```
## [1] "121"                       "NFFD"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "122"                  "NFFD"                 "Winter_cold_injury_p"
```

```
## [1] "123"                        "MCMT"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "124"           "MAT"           "Root_weight_p"
```

```
## [1] "125"                  "MAP"                  "Spring_cold_injury_p"
```

```
## [1] "126"        "FFP"        "Diameter_p"
```

```
## [1] "127"      "FFP"      "Budset_p"
```

```
## [1] "128"                  "Eref"                 "Winter_cold_injury_p"
```

```
## [1] "129"                        "ELEVATION"                 
## [3] "X95_growth_complete_days_p"
```

```
## [1] "130"               "eFFP"              "Height_season_1_p"
```

```
## [1] "131"            "DD5"            "Shoot_weight_p"
```

```
## [1] "132"                   "TD"                    "root_wt__shoot_wt_p_1"
```

```
## [1] "133"        "TD"         "Diameter_p"
```

```
## [1] "134"                 "SHM"                 "root_wt__shoot_wt_p"
```

```
## [1] "135"                 "NFFD"                "root_wt__shoot_wt_p"
```

```
## [1] "136"            "NFFD"           "Shoot_weight_p"
```

```
## [1] "137"        "MWMT"       "Budbreak_p"
```

```
## [1] "138"                 "MAP"                 "root_wt__shoot_wt_p"
```

```
## [1] "139"               "MAP"               "Max_growth_rate_p"
```

```
## [1] "140"            "LAT"            "Shoot_weight_p"
```

```
## [1] "141"                       "FFP"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "142"                   "EXT"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "143"                        "EXT"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "144"                "EXT"                "Fall_cold_injury_p"
```

```
## [1] "145"               "ELEVATION"         "Height_season_1_p"
```

```
## [1] "146"                       "DD5"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "147"                 "bFFP"                "root_wt__shoot_wt_p"
```

```
## [1] "148"      "TD"       "Budset_p"
```

```
## [1] "149"        "SHM"        "Diameter_p"
```

```
## [1] "150"                   "PAS"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "151"               "PAS"               "Height_season_2_p"
```

```
## [1] "152"               "MCMT"              "Height_season_2_p"
```

```
## [1] "153"               "MCMT"              "Height_season_1_p"
```

```
## [1] "154"                  "LAT"                  "Linear_growth_days_p"
```

```
## [1] "155"        "LAT"        "Budbreak_p"
```

```
## [1] "156"               "Eref"              "Max_growth_rate_p"
```

```
## [1] "157"                        "Eref"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "158"                   "ELEVATION"             "root_wt__shoot_wt_p_1"
```

```
## [1] "159"                "ELEVATION"          "Fall_cold_injury_p"
```

```
## [1] "160"               "eFFP"              "Max_growth_rate_p"
```

```
## [1] "161"               "DD5"               "Max_growth_rate_p"
```

```
## [1] "162"                   "CMD"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "163"            "bFFP"           "Shoot_weight_p"
```

```
## [1] "164"            "SHM"            "Shoot_weight_p"
```

```
## [1] "165"            "MWMT"           "Shoot_weight_p"
```

```
## [1] "166"        "MAT"        "Diameter_p"
```

```
## [1] "167"        "MSP"        "Budbreak_p"
```

```
## [1] "168"                 "FFP"                 "root_wt__shoot_wt_p"
```

```
## [1] "169"                  "FFP"                  "Winter_cold_injury_p"
```

```
## [1] "170"               "EMT"               "Height_season_1_p"
```

```
## [1] "171"                        "CMD"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "172"                "CMD"                "Fall_cold_injury_p"
```

```
## [1] "173"            "AHM"            "Shoot_weight_p"
```

```
## [1] "174"      "AHM"      "Budset_p"
```

```
## [1] "175"                 "TD"                  "root_wt__shoot_wt_p"
```

```
## [1] "176"               "TD"                "Max_growth_rate_p"
```

```
## [1] "177"                        "SHM"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "178"                 "PAS"                 "root_wt__shoot_wt_p"
```

```
## [1] "179"                  "PAS"                  "Linear_growth_days_p"
```

```
## [1] "180"        "PAS"        "Diameter_p"
```

```
## [1] "181"           "PAS"           "Root_weight_p"
```

```
## [1] "182"        "NFFD"       "Budbreak_p"
```

```
## [1] "183"                   "MWMT"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "184"                       "MWMT"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "185"                  "MCMT"                 "Spring_cold_injury_p"
```

```
## [1] "186"                 "MAT"                 "root_wt__shoot_wt_p"
```

```
## [1] "187"            "LONG"           "Shoot_weight_p"
```

```
## [1] "188"           "MAP"           "Root_weight_p"
```

```
## [1] "189"                       "LAT"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "190"                 "LAT"                 "X5_95_growth_days_p"
```

```
## [1] "191"               "FFP"               "Max_growth_rate_p"
```

```
## [1] "192"            "FFP"            "Shoot_weight_p"
```

```
## [1] "193"                 "EMT"                 "X5_95_growth_days_p"
```

```
## [1] "194"                  "EMT"                  "Spring_cold_injury_p"
```

```
## [1] "195"                       "eFFP"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "196"                "DD5"                "Fall_cold_injury_p"
```

```
## [1] "197"            "CMD"            "Shoot_weight_p"
```

```
## [1] "198"                  "AHM"                  "Winter_cold_injury_p"
```

```
## [1] "199"                   "SHM"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "200"                  "SHM"                  "Winter_cold_injury_p"
```

```
## [1] "201"                       "PAS"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "202"                       "MAT"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "203"            "MAT"            "Shoot_weight_p"
```

```
## [1] "204"               "MAP"               "Height_season_2_p"
```

```
## [1] "205"        "MAP"        "Diameter_p"
```

```
## [1] "206"      "LAT"      "Budset_p"
```

```
## [1] "207"                "FFP"                "Fall_cold_injury_p"
```

```
## [1] "208"            "EXT"            "Shoot_weight_p"
```

```
## [1] "209"                   "Eref"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "210"                  "ELEVATION"            "Linear_growth_days_p"
```

```
## [1] "211"                  "ELEVATION"            "Winter_cold_injury_p"
```

```
## [1] "212"                   "DD5"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "213"                  "DD5"                  "Spring_cold_injury_p"
```

```
## [1] "214"                  "CMD"                  "Winter_cold_injury_p"
```

```
## [1] "215"      "bFFP"     "Budset_p"
```

```
## [1] "216"                  "AHM"                  "Spring_cold_injury_p"
```

```
## [1] "217"           "TD"            "Root_weight_p"
```

```
## [1] "218"                 "PAS"                 "X5_95_growth_days_p"
```

```
## [1] "219"            "PAS"            "Shoot_weight_p"
```

```
## [1] "220"               "MWMT"              "Max_growth_rate_p"
```

```
## [1] "221"                 "MCMT"                "X5_95_growth_days_p"
```

```
## [1] "222"               "LONG"              "Height_season_1_p"
```

```
## [1] "223"                 "MSP"                 "root_wt__shoot_wt_p"
```

```
## [1] "224"           "MSP"           "Root_weight_p"
```

```
## [1] "225"                  "MAP"                  "Linear_growth_days_p"
```

```
## [1] "226"                        "MAP"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "227"               "LAT"               "Max_growth_rate_p"
```

```
## [1] "228"               "EXT"               "Max_growth_rate_p"
```

```
## [1] "229"                  "EXT"                  "Winter_cold_injury_p"
```

```
## [1] "230"                  "DD5"                  "Winter_cold_injury_p"
```

```
## [1] "231"                   "bFFP"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "232"                  "TD"                   "Linear_growth_days_p"
```

```
## [1] "233"        "TD"         "Budbreak_p"
```

```
## [1] "234"                  "MWMT"                 "Spring_cold_injury_p"
```

```
## [1] "235"               "LONG"              "Height_season_2_p"
```

```
## [1] "236"               "MSP"               "Height_season_1_p"
```

```
## [1] "237"            "MSP"            "Shoot_weight_p"
```

```
## [1] "238"                  "MSP"                  "Winter_cold_injury_p"
```

```
## [1] "239"               "MAP"               "Height_season_1_p"
```

```
## [1] "240"            "MAP"            "Shoot_weight_p"
```

```
## [1] "241"                   "LAT"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "242"               "LAT"               "Height_season_1_p"
```

```
## [1] "243"                 "EXT"                 "X5_95_growth_days_p"
```

```
## [1] "244"                 "Eref"                "X5_95_growth_days_p"
```

```
## [1] "245"                 "eFFP"                "X5_95_growth_days_p"
```

```
## [1] "246"           "eFFP"          "Root_weight_p"
```

```
## [1] "247"                 "DD_0"                "X5_95_growth_days_p"
```

```
## [1] "248"           "DD_0"          "Root_weight_p"
```

```
## [1] "249"               "bFFP"              "Max_growth_rate_p"
```

```
## [1] "250"                  "bFFP"                 "Spring_cold_injury_p"
```

```
## [1] "251"                  "AHM"                  "Linear_growth_days_p"
```

```
## [1] "252"                  "MWMT"                 "Winter_cold_injury_p"
```

```
## [1] "253"                   "MSP"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "254"                        "MSP"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "255"        "MSP"        "Diameter_p"
```

```
## [1] "256"                 "MAP"                 "X5_95_growth_days_p"
```

```
## [1] "257"           "EMT"           "Root_weight_p"
```

```
## [1] "258"        "eFFP"       "Budbreak_p"
```

```
## [1] "259"               "CMD"               "Max_growth_rate_p"
```

```
## [1] "260"        "CMD"        "Budbreak_p"
```

```
## [1] "261"                        "TD"                        
## [3] "X95_growth_complete_days_p"
```

```
## [1] "262"                 "TD"                  "X5_95_growth_days_p"
```

```
## [1] "263"            "TD"             "Shoot_weight_p"
```

```
## [1] "264"               "SHM"               "Max_growth_rate_p"
```

```
## [1] "265"                "SHM"                "Fall_cold_injury_p"
```

```
## [1] "266"                   "NFFD"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "267"                "MWMT"               "Fall_cold_injury_p"
```

```
## [1] "268"                        "LONG"                      
## [3] "X95_growth_complete_days_p"
```

```
## [1] "269"                        "LAT"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "270"        "LAT"        "Diameter_p"
```

```
## [1] "271"                   "FFP"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "272"            "Eref"           "Shoot_weight_p"
```

```
## [1] "273"                  "EMT"                  "Linear_growth_days_p"
```

```
## [1] "274"        "eFFP"       "Diameter_p"
```

```
## [1] "275"            "eFFP"           "Shoot_weight_p"
```

```
## [1] "276"                  "DD_0"                 "Spring_cold_injury_p"
```

```
## [1] "277"                 "CMD"                 "X5_95_growth_days_p"
```

```
## [1] "278"                  "bFFP"                 "Winter_cold_injury_p"
```

```
## [1] "279"                "AHM"                "Fall_cold_injury_p"
```

```
## [1] "280"               "TD"                "Height_season_1_p"
```

```
## [1] "281"                       "TD"                       
## [3] "X5_growth_complete_days_p"
```

```
## [1] "282"                       "SHM"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "283"        "SHM"        "Budbreak_p"
```

```
## [1] "284"               "PAS"               "Max_growth_rate_p"
```

```
## [1] "285"                 "NFFD"                "X5_95_growth_days_p"
```

```
## [1] "286"                 "MWMT"                "X5_95_growth_days_p"
```

```
## [1] "287"                   "MCMT"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "288"                  "MCMT"                 "Linear_growth_days_p"
```

```
## [1] "289"                 "MAT"                 "X5_95_growth_days_p"
```

```
## [1] "290"                       "LONG"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "291"                 "LAT"                 "root_wt__shoot_wt_p"
```

```
## [1] "292"                 "FFP"                 "X5_95_growth_days_p"
```

```
## [1] "293"                  "EXT"                  "Spring_cold_injury_p"
```

```
## [1] "294"                 "eFFP"                "root_wt__shoot_wt_p"
```

```
## [1] "295"                  "eFFP"                 "Spring_cold_injury_p"
```

```
## [1] "296"                 "DD_0"                "root_wt__shoot_wt_p"
```

```
## [1] "297"                       "AHM"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "298"               "TD"                "Height_season_2_p"
```

```
## [1] "299"                 "SHM"                 "X5_95_growth_days_p"
```

```
## [1] "300"                        "PAS"                       
## [3] "X95_growth_complete_days_p"
```

```
## [1] "301"           "MCMT"          "Root_weight_p"
```

```
## [1] "302"               "LONG"              "Max_growth_rate_p"
```

```
## [1] "303"                  "LONG"                 "Spring_cold_injury_p"
```

```
## [1] "304"                "MSP"                "Fall_cold_injury_p"
```

```
## [1] "305"                  "FFP"                  "Spring_cold_injury_p"
```

```
## [1] "306"        "EMT"        "Budbreak_p"
```

```
## [1] "307"               "ELEVATION"         "Max_growth_rate_p"
```

```
## [1] "308"                  "DD5"                  "Linear_growth_days_p"
```

```
## [1] "309"        "DD_0"       "Diameter_p"
```

```
## [1] "310"            "DD_0"           "Shoot_weight_p"
```

```
## [1] "311"                  "bFFP"                 "Linear_growth_days_p"
```

```
## [1] "312"        "MCMT"       "Budbreak_p"
```

```
## [1] "313"        "MAT"        "Budbreak_p"
```

```
## [1] "314"                  "LONG"                 "Winter_cold_injury_p"
```

```
## [1] "315"                "LONG"               "Fall_cold_injury_p"
```

```
## [1] "316"                 "MSP"                 "X5_95_growth_days_p"
```

```
## [1] "317"               "LAT"               "Height_season_2_p"
```

```
## [1] "318"           "LAT"           "Root_weight_p"
```

```
## [1] "319"                   "EMT"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "320"                 "EMT"                 "root_wt__shoot_wt_p"
```

```
## [1] "321"            "EMT"            "Shoot_weight_p"
```

```
## [1] "322"                 "DD5"                 "X5_95_growth_days_p"
```

```
## [1] "323"                  "DD_0"                 "Linear_growth_days_p"
```

```
## [1] "324"                       "DD_0"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "325"                       "CMD"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "326"                 "AHM"                 "X5_95_growth_days_p"
```

```
## [1] "327"                  "SHM"                  "Linear_growth_days_p"
```

```
## [1] "328"                  "NFFD"                 "Linear_growth_days_p"
```

```
## [1] "329"                  "NFFD"                 "Spring_cold_injury_p"
```

```
## [1] "330"                  "MWMT"                 "Linear_growth_days_p"
```

```
## [1] "331"                   "MAT"                   "root_wt__shoot_wt_p_1"
```

```
## [1] "332"                  "LONG"                 "Linear_growth_days_p"
```

```
## [1] "333"               "MSP"               "Max_growth_rate_p"
```

```
## [1] "334"                       "EXT"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "335"                       "Eref"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "336"        "Eref"       "Budbreak_p"
```

```
## [1] "337"        "EMT"        "Diameter_p"
```

```
## [1] "338"                  "eFFP"                 "Linear_growth_days_p"
```

```
## [1] "339"        "DD_0"       "Budbreak_p"
```

```
## [1] "340"                 "bFFP"                "X5_95_growth_days_p"
```

```
## [1] "341"                "bFFP"               "Fall_cold_injury_p"
```

```
## [1] "342"               "AHM"               "Max_growth_rate_p"
```

```
## [1] "343"                  "SHM"                  "Spring_cold_injury_p"
```

```
## [1] "344"                 "MCMT"                "root_wt__shoot_wt_p"
```

```
## [1] "345"            "MCMT"           "Shoot_weight_p"
```

```
## [1] "346"        "LONG"       "Diameter_p"
```

```
## [1] "347"        "LONG"       "Budbreak_p"
```

```
## [1] "348"                       "MSP"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "349"                       "MAP"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "350"        "MAP"        "Budbreak_p"
```

```
## [1] "351"                  "FFP"                  "Linear_growth_days_p"
```

```
## [1] "352"        "EXT"        "Budbreak_p"
```

```
## [1] "353"                  "Eref"                 "Linear_growth_days_p"
```

```
## [1] "354"                  "Eref"                 "Spring_cold_injury_p"
```

```
## [1] "355"                       "EMT"                      
## [3] "X5_growth_complete_days_p"
```

```
## [1] "356"                   "DD_0"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "357"        "AHM"        "Budbreak_p"
```

```
## [1] "358"                       "MCMT"                     
## [3] "X5_growth_complete_days_p"
```

```
## [1] "359"                  "MAT"                  "Linear_growth_days_p"
```

```
## [1] "360"                   "LONG"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "361"                 "LONG"                "root_wt__shoot_wt_p"
```

```
## [1] "362"           "LONG"          "Root_weight_p"
```

```
## [1] "363"                  "MSP"                  "Linear_growth_days_p"
```

```
## [1] "364"                  "MSP"                  "Spring_cold_injury_p"
```

```
## [1] "365"      "MAP"      "Budset_p"
```

```
## [1] "366"                  "EXT"                  "Linear_growth_days_p"
```

```
## [1] "367"                   "eFFP"                  "root_wt__shoot_wt_p_1"
```

```
## [1] "368"                  "CMD"                  "Linear_growth_days_p"
```

```
## [1] "369"        "MCMT"       "Diameter_p"
```

```
## [1] "370"                  "MAT"                  "Spring_cold_injury_p"
```

```
## [1] "371"                 "LONG"                "X5_95_growth_days_p"
```

```
## [1] "372"                 "ELEVATION"           "X5_95_growth_days_p"
```

```
## [1] "373"                  "CMD"                  "Spring_cold_injury_p"
```

```
## [1] "374"                        "AHM"                       
## [3] "X95_growth_complete_days_p"
```

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
##     spruce_correlation pine_score_super_g_raw pine_score_super_g
## 190              -0.74               24.15438          19.454100
## 36                0.33               18.97972          15.099421
## 275               0.76               22.03100          21.998084
## 342               0.25               12.25069           7.654826
## 191              -0.62               22.56486          14.376012
## 54               -0.74               21.27823          20.898606
##     pine_score_super_t_raw pine_score_super_t pine_score_ALL
## 190               29.81696           25.53098       9.436613
## 36                21.32336           16.67241      21.309085
## 275               28.07914           24.64387       8.481924
## 342               13.85500            9.37160      16.011148
## 191               27.64752           19.38109       7.823803
## 54                26.74571           24.30750      11.937639
```

```r
png(paste("analysis/Scores.png"), width=5, height=18, units="in", res=300)
 par(mfrow=c(5,1), mar=c(4,4,4,1), oma=c(0,0,2,0))
  plot(abs(PE$pine_correlation), PE$pine_score_ALL,  ylim=c(-50,50))
  m1 <- lm(PE$pine_score_ALL~abs(PE$pine_correlation) + PE$Environment + PE$Phenotype)
  summary(m1)
```

```
## 
## Call:
## lm(formula = PE$pine_score_ALL ~ abs(PE$pine_correlation) + PE$Environment + 
##     PE$Phenotype)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -31.717  -3.859   0.884   4.328  35.085 
## 
## Coefficients:
##                                         Estimate Std. Error t value
## (Intercept)                              -4.1871     3.2403  -1.292
## abs(PE$pine_correlation)                 26.8304     7.6852   3.491
## PE$EnvironmentbFFP                       -0.6296     3.4222  -0.184
## PE$EnvironmentCMD                         9.4724     3.4226   2.768
## PE$EnvironmentDD0                         3.1949     3.4218   0.934
## PE$EnvironmentDD5                         4.0659     3.4395   1.182
## PE$EnvironmenteFFP                        1.5178     3.4200   0.444
## PE$EnvironmentElevation                  -2.5606     3.4409  -0.744
## PE$EnvironmentEMT                         6.7718     3.4230   1.978
## PE$EnvironmentEref                       12.1159     3.4260   3.537
## PE$EnvironmentEXT                        12.1576     3.4289   3.546
## PE$EnvironmentFFP                         1.9058     3.4226   0.557
## PE$EnvironmentLatitude                   -4.0235     3.4222  -1.176
## PE$Environmentlog(MAP)                   -5.6080     3.4199  -1.640
## PE$Environmentlog(MSP)                    2.9951     3.4280   0.874
## PE$EnvironmentLongitude                  -2.8317     3.4565  -0.819
## PE$EnvironmentMAT                         9.4464     3.4268   2.757
## PE$EnvironmentMCMT                        1.2307     3.4202   0.360
## PE$EnvironmentMWMT                        3.4532     3.4286   1.007
## PE$EnvironmentNFFD                        2.7120     3.4325   0.790
## PE$EnvironmentPAS                        -8.2039     3.4236  -2.396
## PE$EnvironmentSHM                         4.5251     3.4200   1.323
## PE$EnvironmentTD                         -3.5039     3.4200  -1.025
## PE$PhenotypeBudSet_day                    8.5470     3.1965   2.674
## PE$PhenotypeColdInjuryFallS02_Mean        7.3358     3.0987   2.367
## PE$PhenotypeColdInjuryMidWinterS02_Mean   7.1861     3.0735   2.338
## PE$PhenotypeColdInjurySpringS02_Mean      0.4124     3.0067   0.137
## PE$PhenotypeFinalRootWeightS02_g          6.0699     3.0277   2.005
## PE$PhenotypeFinalShootWeightS02_g         6.8783     3.0082   2.286
## PE$PhenotypeFinalStemDiamS02_mm           5.2766     3.0196   1.747
## PE$PhenotypeGth5.95pctS02_Days            9.3599     3.0146   3.105
## PE$PhenotypeGth5pctS02_Day                7.3547     3.0074   2.446
## PE$PhenotypeGth95pctS02_Day               9.5682     3.0265   3.161
## PE$PhenotypeHtFinalS01_mm                 4.1294     3.0961   1.334
## PE$PhenotypeHtFinalS02_mm                 4.3212     3.1415   1.376
## PE$PhenotypeLinearGthS02_Days             0.1106     3.0245   0.037
## PE$PhenotypeMaxGthRate_mm_Day           -12.8789     3.0102  -4.278
## PE$PhenotypeRootShoot_Ratio               6.5459     3.0156   2.171
## PE$PhenotypeTotalWeight                   2.1784     3.0063   0.725
##                                         Pr(>|t|)    
## (Intercept)                             0.197178    
## abs(PE$pine_correlation)                0.000545 ***
## PE$EnvironmentbFFP                      0.854141    
## PE$EnvironmentCMD                       0.005961 ** 
## PE$EnvironmentDD0                       0.351129    
## PE$EnvironmentDD5                       0.237987    
## PE$EnvironmenteFFP                      0.657458    
## PE$EnvironmentElevation                 0.457299    
## PE$EnvironmentEMT                       0.048709 *  
## PE$EnvironmentEref                      0.000463 ***
## PE$EnvironmentEXT                       0.000447 ***
## PE$EnvironmentFFP                       0.578014    
## PE$EnvironmentLatitude                  0.240560    
## PE$Environmentlog(MAP)                  0.101983    
## PE$Environmentlog(MSP)                  0.382896    
## PE$EnvironmentLongitude                 0.413233    
## PE$EnvironmentMAT                       0.006160 ** 
## PE$EnvironmentMCMT                      0.719201    
## PE$EnvironmentMWMT                      0.314585    
## PE$EnvironmentNFFD                      0.430031    
## PE$EnvironmentPAS                       0.017109 *  
## PE$EnvironmentSHM                       0.186686    
## PE$EnvironmentTD                        0.306330    
## PE$PhenotypeBudSet_day                  0.007866 ** 
## PE$PhenotypeColdInjuryFallS02_Mean      0.018480 *  
## PE$PhenotypeColdInjuryMidWinterS02_Mean 0.019971 *  
## PE$PhenotypeColdInjurySpringS02_Mean    0.890996    
## PE$PhenotypeFinalRootWeightS02_g        0.045791 *  
## PE$PhenotypeFinalShootWeightS02_g       0.022849 *  
## PE$PhenotypeFinalStemDiamS02_mm         0.081473 .  
## PE$PhenotypeGth5.95pctS02_Days          0.002066 ** 
## PE$PhenotypeGth5pctS02_Day              0.014979 *  
## PE$PhenotypeGth95pctS02_Day             0.001713 ** 
## PE$PhenotypeHtFinalS01_mm               0.183202    
## PE$PhenotypeHtFinalS02_mm               0.169884    
## PE$PhenotypeLinearGthS02_Days           0.970858    
## PE$PhenotypeMaxGthRate_mm_Day           2.46e-05 ***
## PE$PhenotypeRootShoot_Ratio             0.030658 *  
## PE$PhenotypeTotalWeight                 0.469205    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9.971 on 335 degrees of freedom
## Multiple R-squared:  0.4336,	Adjusted R-squared:  0.3693 
## F-statistic: 6.748 on 38 and 335 DF,  p-value: < 2.2e-16
```

```r
  abline(lm(PE$pine_score_ALL~abs(PE$pine_correlation)))

  plot(abs(PE$pine_correlation), PE$pine_score_super_t,  ylim=c(-50,50))
  m2 <- lm(PE$pine_score_super_t~abs(PE$pine_correlation) + PE$Environment + PE$Phenotype)
  summary(m2)
```

```
## 
## Call:
## lm(formula = PE$pine_score_super_t ~ abs(PE$pine_correlation) + 
##     PE$Environment + PE$Phenotype)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -33.507  -7.514   0.197   6.465  57.339 
## 
## Coefficients:
##                                         Estimate Std. Error t value
## (Intercept)                              -9.5782     4.1449  -2.311
## abs(PE$pine_correlation)                 29.9422     9.8307   3.046
## PE$EnvironmentbFFP                       -2.5861     4.3776  -0.591
## PE$EnvironmentCMD                         5.0377     4.3781   1.151
## PE$EnvironmentDD0                         7.7027     4.3770   1.760
## PE$EnvironmentDD5                         0.9489     4.3997   0.216
## PE$EnvironmenteFFP                        1.6648     4.3747   0.381
## PE$EnvironmentElevation                  -2.9867     4.4015  -0.679
## PE$EnvironmentEMT                         8.6839     4.3786   1.983
## PE$EnvironmentEref                        7.8218     4.3824   1.785
## PE$EnvironmentEXT                         4.4277     4.3862   1.009
## PE$EnvironmentFFP                         1.5151     4.3781   0.346
## PE$EnvironmentLatitude                   -3.6078     4.3776  -0.824
## PE$Environmentlog(MAP)                   -1.5652     4.3747  -0.358
## PE$Environmentlog(MSP)                    0.0486     4.3850   0.011
## PE$EnvironmentLongitude                  -0.5977     4.4215  -0.135
## PE$EnvironmentMAT                        12.0030     4.3834   2.738
## PE$EnvironmentMCMT                        7.1064     4.3750   1.624
## PE$EnvironmentMWMT                        0.2073     4.3858   0.047
## PE$EnvironmentNFFD                        2.9310     4.3908   0.668
## PE$EnvironmentPAS                        -2.8626     4.3793  -0.654
## PE$EnvironmentSHM                         0.9007     4.3747   0.206
## PE$EnvironmentTD                         -3.1039     4.3748  -0.710
## PE$PhenotypeBudSet_day                   13.7412     4.0889   3.361
## PE$PhenotypeColdInjuryFallS02_Mean       15.2764     3.9637   3.854
## PE$PhenotypeColdInjuryMidWinterS02_Mean  14.6651     3.9315   3.730
## PE$PhenotypeColdInjurySpringS02_Mean     14.4442     3.8461   3.756
## PE$PhenotypeFinalRootWeightS02_g          7.4025     3.8730   1.911
## PE$PhenotypeFinalShootWeightS02_g         7.3695     3.8480   1.915
## PE$PhenotypeFinalStemDiamS02_mm           7.7904     3.8625   2.017
## PE$PhenotypeGth5.95pctS02_Days            3.9848     3.8561   1.033
## PE$PhenotypeGth5pctS02_Day                2.5415     3.8470   0.661
## PE$PhenotypeGth95pctS02_Day              -3.1831     3.8714  -0.822
## PE$PhenotypeHtFinalS01_mm                 6.1315     3.9605   1.548
## PE$PhenotypeHtFinalS02_mm                 3.7935     4.0185   0.944
## PE$PhenotypeLinearGthS02_Days            11.2609     3.8689   2.911
## PE$PhenotypeMaxGthRate_mm_Day            -0.4456     3.8506  -0.116
## PE$PhenotypeRootShoot_Ratio               2.8733     3.8575   0.745
## PE$PhenotypeTotalWeight                   1.3059     3.8456   0.340
##                                         Pr(>|t|)    
## (Intercept)                             0.021448 *  
## abs(PE$pine_correlation)                0.002505 ** 
## PE$EnvironmentbFFP                      0.555084    
## PE$EnvironmentCMD                       0.250690    
## PE$EnvironmentDD0                       0.079356 .  
## PE$EnvironmentDD5                       0.829375    
## PE$EnvironmenteFFP                      0.703782    
## PE$EnvironmentElevation                 0.497893    
## PE$EnvironmentEMT                       0.048151 *  
## PE$EnvironmentEref                      0.075193 .  
## PE$EnvironmentEXT                       0.313477    
## PE$EnvironmentFFP                       0.729519    
## PE$EnvironmentLatitude                  0.410447    
## PE$Environmentlog(MAP)                  0.720733    
## PE$Environmentlog(MSP)                  0.991164    
## PE$EnvironmentLongitude                 0.892559    
## PE$EnvironmentMAT                       0.006507 ** 
## PE$EnvironmentMCMT                      0.105248    
## PE$EnvironmentMWMT                      0.962334    
## PE$EnvironmentNFFD                      0.504884    
## PE$EnvironmentPAS                       0.513772    
## PE$EnvironmentSHM                       0.837007    
## PE$EnvironmentTD                        0.478507    
## PE$PhenotypeBudSet_day                  0.000868 ***
## PE$PhenotypeColdInjuryFallS02_Mean      0.000139 ***
## PE$PhenotypeColdInjuryMidWinterS02_Mean 0.000225 ***
## PE$PhenotypeColdInjurySpringS02_Mean    0.000204 ***
## PE$PhenotypeFinalRootWeightS02_g        0.056816 .  
## PE$PhenotypeFinalShootWeightS02_g       0.056326 .  
## PE$PhenotypeFinalStemDiamS02_mm         0.044500 *  
## PE$PhenotypeGth5.95pctS02_Days          0.302180    
## PE$PhenotypeGth5pctS02_Day              0.509296    
## PE$PhenotypeGth95pctS02_Day             0.411539    
## PE$PhenotypeHtFinalS01_mm               0.122530    
## PE$PhenotypeHtFinalS02_mm               0.345842    
## PE$PhenotypeLinearGthS02_Days           0.003848 ** 
## PE$PhenotypeMaxGthRate_mm_Day           0.907940    
## PE$PhenotypeRootShoot_Ratio             0.456881    
## PE$PhenotypeTotalWeight                 0.734383    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.75 on 335 degrees of freedom
## Multiple R-squared:  0.3051,	Adjusted R-squared:  0.2263 
## F-statistic: 3.871 on 38 and 335 DF,  p-value: 1.025e-11
```

```r
  abline(lm(PE$pine_score_super_t~abs(PE$pine_correlation)))

  plot(abs(PE$pine_correlation), PE$pine_score_super_t_raw,  ylim=c(-50,50))
  m3 <- lm(PE$pine_score_super_t_raw~abs(PE$pine_correlation) + PE$Environment + PE$Phenotype)
  summary(m3)
```

```
## 
## Call:
## lm(formula = PE$pine_score_super_t_raw ~ abs(PE$pine_correlation) + 
##     PE$Environment + PE$Phenotype)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -35.566  -8.276  -0.612   7.843  56.322 
## 
## Coefficients:
##                                         Estimate Std. Error t value
## (Intercept)                             -14.4876     4.9466  -2.929
## abs(PE$pine_correlation)                 32.7588    11.7321   2.792
## PE$EnvironmentbFFP                       -3.8030     5.2244  -0.728
## PE$EnvironmentCMD                         9.5557     5.2249   1.829
## PE$EnvironmentDD0                        13.1935     5.2236   2.526
## PE$EnvironmentDD5                         0.5029     5.2506   0.096
## PE$EnvironmenteFFP                        4.4544     5.2209   0.853
## PE$EnvironmentElevation                  -2.8135     5.2529  -0.536
## PE$EnvironmentEMT                        14.6652     5.2255   2.806
## PE$EnvironmentEref                        9.7430     5.2300   1.863
## PE$EnvironmentEXT                         5.6614     5.2346   1.082
## PE$EnvironmentFFP                         1.7146     5.2249   0.328
## PE$EnvironmentLatitude                   -2.3968     5.2244  -0.459
## PE$Environmentlog(MAP)                   -1.7864     5.2208  -0.342
## PE$Environmentlog(MSP)                    5.7854     5.2331   1.106
## PE$EnvironmentLongitude                   0.7187     5.2767   0.136
## PE$EnvironmentMAT                        15.4344     5.2313   2.950
## PE$EnvironmentMCMT                       13.5381     5.2212   2.593
## PE$EnvironmentMWMT                       -1.4774     5.2341  -0.282
## PE$EnvironmentNFFD                        4.4892     5.2400   0.857
## PE$EnvironmentPAS                        -3.8629     5.2264  -0.739
## PE$EnvironmentSHM                         4.9764     5.2209   0.953
## PE$EnvironmentTD                         -4.1381     5.2209  -0.793
## PE$PhenotypeBudSet_day                   17.6303     4.8798   3.613
## PE$PhenotypeColdInjuryFallS02_Mean       18.3102     4.7304   3.871
## PE$PhenotypeColdInjuryMidWinterS02_Mean  18.6953     4.6919   3.985
## PE$PhenotypeColdInjurySpringS02_Mean     20.7253     4.5900   4.515
## PE$PhenotypeFinalRootWeightS02_g          9.5772     4.6221   2.072
## PE$PhenotypeFinalShootWeightS02_g         9.5082     4.5923   2.070
## PE$PhenotypeFinalStemDiamS02_mm           8.6753     4.6096   1.882
## PE$PhenotypeGth5.95pctS02_Days           10.4988     4.6020   2.281
## PE$PhenotypeGth5pctS02_Day                4.7405     4.5911   1.033
## PE$PhenotypeGth95pctS02_Day               6.2532     4.6202   1.353
## PE$PhenotypeHtFinalS01_mm                 6.0232     4.7265   1.274
## PE$PhenotypeHtFinalS02_mm                 5.7435     4.7957   1.198
## PE$PhenotypeLinearGthS02_Days            15.6371     4.6172   3.387
## PE$PhenotypeMaxGthRate_mm_Day            -1.4162     4.5953  -0.308
## PE$PhenotypeRootShoot_Ratio               7.8893     4.6036   1.714
## PE$PhenotypeTotalWeight                  -0.2136     4.5894  -0.047
##                                         Pr(>|t|)    
## (Intercept)                             0.003636 ** 
## abs(PE$pine_correlation)                0.005535 ** 
## PE$EnvironmentbFFP                      0.467162    
## PE$EnvironmentCMD                       0.068306 .  
## PE$EnvironmentDD0                       0.012006 *  
## PE$EnvironmentDD5                       0.923753    
## PE$EnvironmenteFFP                      0.394160    
## PE$EnvironmentElevation                 0.592585    
## PE$EnvironmentEMT                       0.005301 ** 
## PE$EnvironmentEref                      0.063351 .  
## PE$EnvironmentEXT                       0.280237    
## PE$EnvironmentFFP                       0.742999    
## PE$EnvironmentLatitude                  0.646688    
## PE$Environmentlog(MAP)                  0.732433    
## PE$Environmentlog(MSP)                  0.269721    
## PE$EnvironmentLongitude                 0.891743    
## PE$EnvironmentMAT                       0.003397 ** 
## PE$EnvironmentMCMT                      0.009935 ** 
## PE$EnvironmentMWMT                      0.777913    
## PE$EnvironmentNFFD                      0.392216    
## PE$EnvironmentPAS                       0.460355    
## PE$EnvironmentSHM                       0.341185    
## PE$EnvironmentTD                        0.428575    
## PE$PhenotypeBudSet_day                  0.000349 ***
## PE$PhenotypeColdInjuryFallS02_Mean      0.000130 ***
## PE$PhenotypeColdInjuryMidWinterS02_Mean 8.30e-05 ***
## PE$PhenotypeColdInjurySpringS02_Mean    8.77e-06 ***
## PE$PhenotypeFinalRootWeightS02_g        0.039024 *  
## PE$PhenotypeFinalShootWeightS02_g       0.039175 *  
## PE$PhenotypeFinalStemDiamS02_mm         0.060702 .  
## PE$PhenotypeGth5.95pctS02_Days          0.023155 *  
## PE$PhenotypeGth5pctS02_Day              0.302567    
## PE$PhenotypeGth95pctS02_Day             0.176828    
## PE$PhenotypeHtFinalS01_mm               0.203428    
## PE$PhenotypeHtFinalS02_mm               0.231906    
## PE$PhenotypeLinearGthS02_Days           0.000792 ***
## PE$PhenotypeMaxGthRate_mm_Day           0.758127    
## PE$PhenotypeRootShoot_Ratio             0.087505 .  
## PE$PhenotypeTotalWeight                 0.962915    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 15.22 on 335 degrees of freedom
## Multiple R-squared:  0.3213,	Adjusted R-squared:  0.2443 
## F-statistic: 4.173 on 38 and 335 DF,  p-value: 4.964e-13
```

```r
  abline(lm(PE$pine_score_super_t_raw~abs(PE$pine_correlation)))

  plot(abs(PE$pine_correlation), PE$pine_score_super_g,  ylim=c(-50,50))
  m4<- lm(PE$pine_score_super_g~abs(PE$pine_correlation) + PE$Environment + PE$Phenotype)
  summary(m4)
```

```
## 
## Call:
## lm(formula = PE$pine_score_super_g ~ abs(PE$pine_correlation) + 
##     PE$Environment + PE$Phenotype)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -32.016  -8.164   0.241   6.561  61.523 
## 
## Coefficients:
##                                         Estimate Std. Error t value
## (Intercept)                             -10.3813     4.2961  -2.416
## abs(PE$pine_correlation)                 25.0596    10.1893   2.459
## PE$EnvironmentbFFP                       -1.1360     4.5373  -0.250
## PE$EnvironmentCMD                         5.1528     4.5378   1.136
## PE$EnvironmentDD0                         7.7059     4.5367   1.699
## PE$EnvironmentDD5                         0.9152     4.5601   0.201
## PE$EnvironmenteFFP                        2.6160     4.5343   0.577
## PE$EnvironmentElevation                  -0.7566     4.5621  -0.166
## PE$EnvironmentEMT                         9.4942     4.5383   2.092
## PE$EnvironmentEref                        5.9686     4.5422   1.314
## PE$EnvironmentEXT                         2.6503     4.5462   0.583
## PE$EnvironmentFFP                         2.7125     4.5378   0.598
## PE$EnvironmentLatitude                   -2.4756     4.5373  -0.546
## PE$Environmentlog(MAP)                   -0.7335     4.5343  -0.162
## PE$Environmentlog(MSP)                    0.8490     4.5449   0.187
## PE$EnvironmentLongitude                   0.3379     4.5828   0.074
## PE$EnvironmentMAT                        11.7943     4.5433   2.596
## PE$EnvironmentMCMT                        7.5600     4.5346   1.667
## PE$EnvironmentMWMT                       -0.1633     4.5458  -0.036
## PE$EnvironmentNFFD                        3.7528     4.5509   0.825
## PE$EnvironmentPAS                        -1.8165     4.5391  -0.400
## PE$EnvironmentSHM                         1.5458     4.5343   0.341
## PE$EnvironmentTD                         -2.2481     4.5344  -0.496
## PE$PhenotypeBudSet_day                   15.3850     4.2381   3.630
## PE$PhenotypeColdInjuryFallS02_Mean       14.2855     4.1083   3.477
## PE$PhenotypeColdInjuryMidWinterS02_Mean  13.9810     4.0749   3.431
## PE$PhenotypeColdInjurySpringS02_Mean     15.0554     3.9864   3.777
## PE$PhenotypeFinalRootWeightS02_g          7.5416     4.0142   1.879
## PE$PhenotypeFinalShootWeightS02_g         7.1191     3.9884   1.785
## PE$PhenotypeFinalStemDiamS02_mm           7.4656     4.0034   1.865
## PE$PhenotypeGth5.95pctS02_Days            2.1898     3.9968   0.548
## PE$PhenotypeGth5pctS02_Day                2.6400     3.9874   0.662
## PE$PhenotypeGth95pctS02_Day              -5.2822     4.0126  -1.316
## PE$PhenotypeHtFinalS01_mm                 5.3814     4.1050   1.311
## PE$PhenotypeHtFinalS02_mm                 4.6261     4.1651   1.111
## PE$PhenotypeLinearGthS02_Days            11.7898     4.0100   2.940
## PE$PhenotypeMaxGthRate_mm_Day             0.4463     3.9910   0.112
## PE$PhenotypeRootShoot_Ratio               4.4847     3.9982   1.122
## PE$PhenotypeTotalWeight                   0.7512     3.9859   0.188
##                                         Pr(>|t|)    
## (Intercept)                             0.016205 *  
## abs(PE$pine_correlation)                0.014421 *  
## PE$EnvironmentbFFP                      0.802459    
## PE$EnvironmentCMD                       0.256963    
## PE$EnvironmentDD0                       0.090331 .  
## PE$EnvironmentDD5                       0.841058    
## PE$EnvironmenteFFP                      0.564370    
## PE$EnvironmentElevation                 0.868379    
## PE$EnvironmentEMT                       0.037188 *  
## PE$EnvironmentEref                      0.189739    
## PE$EnvironmentEXT                       0.560310    
## PE$EnvironmentFFP                       0.550399    
## PE$EnvironmentLatitude                  0.585701    
## PE$Environmentlog(MAP)                  0.871584    
## PE$Environmentlog(MSP)                  0.851923    
## PE$EnvironmentLongitude                 0.941259    
## PE$EnvironmentMAT                       0.009848 ** 
## PE$EnvironmentMCMT                      0.096415 .  
## PE$EnvironmentMWMT                      0.971366    
## PE$EnvironmentNFFD                      0.410171    
## PE$EnvironmentPAS                       0.689277    
## PE$EnvironmentSHM                       0.733381    
## PE$EnvironmentTD                        0.620366    
## PE$PhenotypeBudSet_day                  0.000327 ***
## PE$PhenotypeColdInjuryFallS02_Mean      0.000573 ***
## PE$PhenotypeColdInjuryMidWinterS02_Mean 0.000677 ***
## PE$PhenotypeColdInjurySpringS02_Mean    0.000188 ***
## PE$PhenotypeFinalRootWeightS02_g        0.061153 .  
## PE$PhenotypeFinalShootWeightS02_g       0.075172 .  
## PE$PhenotypeFinalStemDiamS02_mm         0.063082 .  
## PE$PhenotypeGth5.95pctS02_Days          0.584132    
## PE$PhenotypeGth5pctS02_Day              0.508361    
## PE$PhenotypeGth95pctS02_Day             0.188944    
## PE$PhenotypeHtFinalS01_mm               0.190774    
## PE$PhenotypeHtFinalS02_mm               0.267503    
## PE$PhenotypeLinearGthS02_Days           0.003509 ** 
## PE$PhenotypeMaxGthRate_mm_Day           0.911025    
## PE$PhenotypeRootShoot_Ratio             0.262806    
## PE$PhenotypeTotalWeight                 0.850635    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13.22 on 335 degrees of freedom
## Multiple R-squared:  0.2787,	Adjusted R-squared:  0.1968 
## F-statistic: 3.406 on 38 and 335 DF,  p-value: 1.096e-09
```

```r
  abline(lm(PE$pine_score_super_g~abs(PE$pine_correlation)))

  plot(abs(PE$pine_correlation), PE$pine_score_super_g_raw,  ylim=c(-50,50))
  m5 <- lm(PE$pine_score_super_g_raw~abs(PE$pine_correlation) + PE$Environment + PE$Phenotype)
  summary(m5)
```

```
## 
## Call:
## lm(formula = PE$pine_score_super_g_raw ~ abs(PE$pine_correlation) + 
##     PE$Environment + PE$Phenotype)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -34.545  -7.973  -0.208   7.304  60.221 
## 
## Coefficients:
##                                         Estimate Std. Error t value
## (Intercept)                             -13.0573     4.5731  -2.855
## abs(PE$pine_correlation)                 25.2499    10.8464   2.328
## PE$EnvironmentbFFP                       -2.5372     4.8300  -0.525
## PE$EnvironmentCMD                         7.6844     4.8305   1.591
## PE$EnvironmentDD0                        11.1080     4.8293   2.300
## PE$EnvironmentDD5                         0.6928     4.8543   0.143
## PE$EnvironmenteFFP                        3.5849     4.8267   0.743
## PE$EnvironmentElevation                  -2.0131     4.8563  -0.415
## PE$EnvironmentEMT                        12.1972     4.8310   2.525
## PE$EnvironmentEref                        8.2439     4.8352   1.705
## PE$EnvironmentEXT                         4.2889     4.8394   0.886
## PE$EnvironmentFFP                         2.1403     4.8305   0.443
## PE$EnvironmentLatitude                   -2.1661     4.8300  -0.448
## PE$Environmentlog(MAP)                   -1.5897     4.8267  -0.329
## PE$Environmentlog(MSP)                    3.2833     4.8381   0.679
## PE$EnvironmentLongitude                   1.4884     4.8783   0.305
## PE$EnvironmentMAT                        13.8948     4.8364   2.873
## PE$EnvironmentMCMT                       11.1569     4.8271   2.311
## PE$EnvironmentMWMT                       -0.6536     4.8389  -0.135
## PE$EnvironmentNFFD                        4.3017     4.8444   0.888
## PE$EnvironmentPAS                        -3.2756     4.8318  -0.678
## PE$EnvironmentSHM                         3.5103     4.8267   0.727
## PE$EnvironmentTD                         -3.5192     4.8268  -0.729
## PE$PhenotypeBudSet_day                   16.7661     4.5114   3.716
## PE$PhenotypeColdInjuryFallS02_Mean       16.8340     4.3733   3.849
## PE$PhenotypeColdInjuryMidWinterS02_Mean  16.6667     4.3377   3.842
## PE$PhenotypeColdInjurySpringS02_Mean     19.0304     4.2435   4.485
## PE$PhenotypeFinalRootWeightS02_g          8.2730     4.2732   1.936
## PE$PhenotypeFinalShootWeightS02_g         7.5968     4.2456   1.789
## PE$PhenotypeFinalStemDiamS02_mm           7.8419     4.2616   1.840
## PE$PhenotypeGth5.95pctS02_Days            7.8953     4.2546   1.856
## PE$PhenotypeGth5pctS02_Day                4.0229     4.2445   0.948
## PE$PhenotypeGth95pctS02_Day               2.6115     4.2714   0.611
## PE$PhenotypeHtFinalS01_mm                 6.4126     4.3697   1.468
## PE$PhenotypeHtFinalS02_mm                 5.7747     4.4337   1.302
## PE$PhenotypeLinearGthS02_Days            15.1068     4.2686   3.539
## PE$PhenotypeMaxGthRate_mm_Day             0.3431     4.2484   0.081
## PE$PhenotypeRootShoot_Ratio               7.2352     4.2561   1.700
## PE$PhenotypeTotalWeight                  -0.3901     4.2429  -0.092
##                                         Pr(>|t|)    
## (Intercept)                             0.004569 ** 
## abs(PE$pine_correlation)                0.020510 *  
## PE$EnvironmentbFFP                      0.599718    
## PE$EnvironmentCMD                       0.112591    
## PE$EnvironmentDD0                       0.022055 *  
## PE$EnvironmentDD5                       0.886595    
## PE$EnvironmenteFFP                      0.458168    
## PE$EnvironmentElevation                 0.678754    
## PE$EnvironmentEMT                       0.012038 *  
## PE$EnvironmentEref                      0.089125 .  
## PE$EnvironmentEXT                       0.376121    
## PE$EnvironmentFFP                       0.657994    
## PE$EnvironmentLatitude                  0.654106    
## PE$Environmentlog(MAP)                  0.742097    
## PE$Environmentlog(MSP)                  0.497832    
## PE$EnvironmentLongitude                 0.760471    
## PE$EnvironmentMAT                       0.004325 ** 
## PE$EnvironmentMCMT                      0.021422 *  
## PE$EnvironmentMWMT                      0.892635    
## PE$EnvironmentNFFD                      0.375193    
## PE$EnvironmentPAS                       0.498283    
## PE$EnvironmentSHM                       0.467571    
## PE$EnvironmentTD                        0.466460    
## PE$PhenotypeBudSet_day                  0.000237 ***
## PE$PhenotypeColdInjuryFallS02_Mean      0.000142 ***
## PE$PhenotypeColdInjuryMidWinterS02_Mean 0.000146 ***
## PE$PhenotypeColdInjurySpringS02_Mean    1.01e-05 ***
## PE$PhenotypeFinalRootWeightS02_g        0.053705 .  
## PE$PhenotypeFinalShootWeightS02_g       0.074464 .  
## PE$PhenotypeFinalStemDiamS02_mm         0.066636 .  
## PE$PhenotypeGth5.95pctS02_Days          0.064371 .  
## PE$PhenotypeGth5pctS02_Day              0.343925    
## PE$PhenotypeGth95pctS02_Day             0.541351    
## PE$PhenotypeHtFinalS01_mm               0.143178    
## PE$PhenotypeHtFinalS02_mm               0.193659    
## PE$PhenotypeLinearGthS02_Days           0.000458 ***
## PE$PhenotypeMaxGthRate_mm_Day           0.935682    
## PE$PhenotypeRootShoot_Ratio             0.090065 .  
## PE$PhenotypeTotalWeight                 0.926801    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 14.07 on 335 degrees of freedom
## Multiple R-squared:  0.2971,	Adjusted R-squared:  0.2173 
## F-statistic: 3.726 on 38 and 335 DF,  p-value: 4.429e-11
```

```r
  abline(lm(PE$pine_score_super_g_raw~abs(PE$pine_correlation)))
dev.off()
```

```
## quartz_off_screen 
##                 7
```
