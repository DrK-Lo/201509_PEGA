
ls()
dt <- load("data/large_files/ResultsWithPEGAScoresByLocusAndPhenEnviComboByPhenotype.Rdata")
ls()
### get head of pine file
names(results_pine2)
### get head of spruce file
names(results_spruce2)

### PINE
### determine which
PinePlieotropyMat <- matrix(NA, 374, 374)
ncolpine <- ncol(results_pine2)
results_pineSub <- results_pine2[,(ncolpine-374+1):ncolpine]
names(results_pineSub)
neword <- order(names(results_pineSub))
results_pineSub <- results_pineSub[,neword]
names(results_pineSub)
nlocpine <- nrow(results_pine2)

### Plieotropy mat Pine
  #17*22=374
  for (i in 1:374){
    if(i%%10==0){print(i)}
     focal <- results_pineSub[,i]
     focalind <- which(rank(-focal,na.last = TRUE)<nlocpine*0.001)
     PinePlieotropyMat[i,] <- colSums(results_pineSub[focalind,], na.rm = TRUE)
  }
  library(fields)




### SPRUCE
### determine which
SprucePlieotropyMat <- matrix(NA, 374, 374)
ncolspruce <- ncol(results_spruce2)
results_spruceSub <- results_spruce2[,(ncolspruce-374+1):ncolspruce]
names(results_spruceSub)
neword <- order(names(results_spruceSub))
results_spruceSub <- results_spruceSub[,neword]
names(results_spruceSub)
nlocspruce <- nrow(results_spruce2)

### Plieotropy mat Spruce
  #17*22=374
  for (i in 1:374){
    if(i%%10==0){print(i)}
     focal <- results_spruceSub[,i]
     focalind <- which(rank(-focal,na.last = TRUE)<nlocspruce*0.001)
     SprucePlieotropyMat[i,] <- colSums(results_spruceSub[focalind,], na.rm = TRUE)
  }

MakePlot <- function(PMmat, matsize=374){
  par(mar=c(8,6,2,6))
  image(PMmat, breaks = seq(-8,8,1), 
             col=two.colors(n=16,start="blue", middle="white", end="red"),
             xaxt="n", yaxt="n",
             axes=FALSE)
  if(matsize==374){
    byenvi <- seq(1,374, by=17)
    a<-seq(0,1,length.out=374)
    axis(1, at=a[seq(1,374,by=17)], tick = FALSE,
       labels=names(results_spruceSub)[byenvi],cex.axis=0.5, las=2)
    axis(2, at=a[seq(17,374,by=17)], tick = FALSE,
       labels=names(results_spruceSub)[byenvi],cex.axis=0.5, las=2)
  }else{
    yo <- seq(0,1, length.out=(nrow(PMmat)))
    axis(1, at=yo, tick = FALSE,
       labels=rownames(PMmat),cex.axis=0.5, las=2)
    axis(2, at=yo, tick = FALSE,
       labels=rownames(PMmat),cex.axis=0.5, las=2)
  }
  image.plot(PMmat, breaks = seq(-8,8,1), 
             col=two.colors(n=16,start="blue", middle="white", end="red"),
             legend.only=TRUE)
}
rownames(PinePlieotropyMat) <- names(results_pineSub)
colnames(PinePlieotropyMat) <- names(results_pineSub)
MakePlot(PinePlieotropyMat, 374)

rownames(SprucePlieotropyMat) <- names(results_spruceSub)
colnames(SprucePlieotropyMat) <- names(results_spruceSub)
MakePlot(SprucePlieotropyMat, 374)

focalEnvi <- c("LAT", "MAP", "MAT", "^FFP")
focalPheno <- c("Max_growth_rate_p", "Fall_cold_injury_p", "Diameter_p")
names(results_spruceSub)
focalind1 <- as.integer(sapply(focalEnvi, grep, x= names(results_spruceSub)))
focalind2 <- as.integer(sapply(focalPheno, grep, x= names(results_spruceSub)))
focalind <- intersect(focalind1, focalind2)
focalPine <- PinePlieotropyMat[focalind,focalind]
focalSpruce <- SprucePlieotropyMat[focalind, focalind]
MakePlot(focalPine, matsize=0)
MakePlot(focalSpruce, matsize=0)
