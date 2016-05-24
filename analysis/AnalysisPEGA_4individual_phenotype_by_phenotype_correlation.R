

ee <- read.table ("~/Desktop/pine_association_phenotype_data.csv",sep = ",", header = T)
head(ee)

ee$fallcoldmean <- (ee$FallColdHardS02_T1It + ee$FallColdHardS02_T3It + ee$FallColdHardS02_T2It)/3

ff <- ee[ee$ExperimentalClimate == "MAT06",]

plot (ff$Height2S02_mm, ff$fallcoldmean)
abline(lm(ff$fallcoldmean~ff$Height2S02_mm))

cor.test (ff$Height2S02_mm, ff$fallcoldmean)
