#### Load PINE Flip and RESULTS files
flip_pine_gwas <- read.table("../data/large_files/flip_pine_v2.cor", header=TRUE)
  nrow(flip_pine_gwas)
  names(flip_pine_gwas)[2] <- "flip_pine_gwas"
  head(flip_pine_gwas)

flip_pine_raw <- read.table("../data/large_files/var_out_GATK3_allhet_pine688_ALL.table_filt10_p95_het7_passFILT2.2_LDformat_TOFLIP_alphabet_lowest_eq_1")
  nrow(flip_pine_raw)
  names(flip_pine_raw)[1] <- names(flip_pine_gwas)[1]
  names(flip_pine_raw)[2] <- "flip_pine_raw"
  head(flip_pine_raw)

flip_pine2 <- merge(flip_pine_gwas, flip_pine_raw, all.x=TRUE)
  nrow(flip_pine2)
  head(flip_pine2)
rm(flip_pine_gwas, flip_pine_raw)

#### read in results_pine if it doesn't exist
  if (!("results_pine" %in% ls())){
  results_pine <- read.table("../data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS_altGCTABAYENV", header=TRUE, comment.char=" ")
    }
  dim(results_pine)
  ### merge flip with results
  results_pine2 <- merge(flip_pine2, results_pine)
  dim(results_pine2)
  head(results_pine2[,1:10])
  results_pine <- results_pine2
  rm(results_pine2, flip_pine2)

  results_pine_head <- names(results_pine)
  results_pine_head

#### Load SPRUCE Flip and RESULTS files
flip_spruce_gwas <- read.table("../data/large_files/flip_spruce_v2.cor", header=TRUE)
  nrow(flip_spruce_gwas)
  names(flip_spruce_gwas)[2] <- "flip_spruce_gwas"
  head(flip_spruce_gwas)

flip_spruce_raw <- read.table("../data/large_files/var_out_GATK3_spruce_ALL.table_filt10_p95_het7_passFILT2.3_LDformat_TOFLIP_alphabet_lowest_eq_1")
  nrow(flip_spruce_raw)
  names(flip_spruce_raw)[1] <- names(flip_spruce_gwas)[1]
  names(flip_spruce_raw)[2] <- "flip_spruce_raw"
  head(flip_spruce_raw)

flip_spruce2 <- merge(flip_spruce_gwas, flip_spruce_raw, all.x=TRUE)
  nrow(flip_spruce2)
rm(flip_spruce_gwas, flip_spruce_raw)
names(flip_spruce2)[1] <- "gcontig__gcontig_pos"
  head(flip_spruce2)

if (!("results_spruce" %in% ls())){
results_spruce <- read.table("../data/large_files/var_out_GATK3_spruce_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS_altGCTABAYENV", header=TRUE, comment.char=" ")
  }
  dim(results_spruce)
  results_spruce$gcontig__gcontig_pos <- paste(results_spruce$gcontig, results_spruce$pos_gcontig, sep="__")
  ### merge flip with results
  results_spruce2 <- merge(flip_spruce2, results_spruce)
  dim(results_spruce2)
  head(results_spruce2[,1:10])
  results_spruce <- results_spruce2
  rm(results_spruce2, flip_spruce2)
  results_spruce_head <- names(results_spruce)
  head(results_spruce_head)

###  WRITE TO FILE
write.table(results_pine,file = "../data/large_files/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS_altGCTABAYENV_flipped", row.names=FALSE)
write.table(results_spruce,file = "../data/large_files/var_out_GATK3_spruce_ALL.summary.ALL.annots.sorted.GOOD.window_RESULTS_altGCTABAYENV_flipped", row.names=FALSE)