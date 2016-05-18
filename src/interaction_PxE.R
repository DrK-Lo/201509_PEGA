

ee <- read.table ("~/Dropbox/desktop/adaptree/Laura_MAT6_analysis/MAT06-BLUEs&Climate-pine.csv", sep = ",", header = T)

#unexpected direction of cold vs. height: taller individuals are more cold tolerant!!
plot (ee$ColdInjuryMidWinterS02_Mean, ee$HeightFinalS02_mm)


#expected correlations for MAT, less cold injury in cold environments:
plot (ee$ColdInjuryMidWinterS02_Mean, ee$MAT)

#overall higher cold injury is found in wetter environments
plot (ee$ColdInjuryMidWinterS02_Mean, ee$log.MAP.)


lm1 <- lm (ee$ColdInjuryMidWinterS02_Mean ~ ee$MAT * ee$log.MAP.)

lm2 <- lm (ee$ColdInjuryMidWinterS02_Mean ~ ee$DD0 * ee$log.MAP.)

anova (lm1)

anova (lm2)

summary (lm1)
summary (lm2)