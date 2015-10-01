
### the *.bim file was used to code alleles for gcta (Hatkan)
### So if the first column is a “G” and the second column is an “A" 
### it coded a GG as a 11 and a AA as a 00.
### So if the first letter is also first alphabetically, the coding is correct
### else, the coding is incorrect

gwas_spruce <- read.table("data/large_files/spruce_gwas.bim")
gwas_pine <- read.table("data/large_files/pine_gwas.bim")


head(gwas_pine)

### Make sure all levels are ACGT
  levels(gwas_spruce$V5)
  levels(gwas_spruce$V6)
  levels(gwas_pine$V5)
  levels(gwas_pine$V6)

  head(gwas_spruce$V5)
  head(as.numeric(gwas_spruce$V5))

### Make spruce indicator (is 1st allele also 1st in alphabet)
  spruce_ind <- as.numeric(as.numeric(gwas_spruce$V5) < 
                             as.numeric(gwas_spruce$V6))
  head(gwas_spruce)
  spruce_ind[spruce_ind==0] <- -1
  gwas_spruce$spruce_ind <- spruce_ind
  dim(gwas_spruce)
  names(gwas_spruce)[2] <- "snp_id"
  head(gwas_spruce)
  ### compare to old flip file
  flip_spruce <- read.table("data/large_files/flip_spruce.cor", header=TRUE)
  dim(flip_spruce)
  head(flip_spruce)
  compare_flip <- merge(flip_spruce, gwas_spruce)
  head(compare_flip)
  dim(compare_flip)
  plot(compare_flip$X1, compare_flip$spruce_ind)
  abline(lm(compare_flip$spruce_ind~compare_flip$X1))
    # so most alleles were coded backwards, and some coded wrong
  write.table(data.frame(snp_id=gwas_spruce$snp_id,
                         X1 = gwas_spruce$spruce_ind),
              file="flip_spruce_v2.cor", row.names=FALSE)

### Pine indicator (is 1st allele also 1st in alphabet)
  pine_ind <- as.numeric(as.numeric(gwas_pine$V5) < 
                             as.numeric(gwas_pine$V6))
  head(gwas_pine)
  head(pine_ind)
  pine_ind[pine_ind==0] <- -1
  gwas_pine$pine_ind <- pine_ind
  dim(gwas_pine)
  names(gwas_pine)[2] <- "gcontig__gcontig_pos"
  head(gwas_pine)
  ### compare to old flip file
  flip_pine <- read.table("data/large_files/flip_pine.cor", header=TRUE)
  dim(flip_pine)
  head(flip_pine)
  compare_flip_pine <- merge(flip_pine, gwas_pine)
  head(compare_flip_pine)
  dim(compare_flip_pine)
  plot(compare_flip_pine$X1, compare_flip_pine$pine_ind)
  abline(lm(compare_flip_pine$pine_ind~compare_flip_pine$X1))
    # so most alleles were coded backwards, and some coded wrong
  write.table(data.frame(gcontig__gcontig_pos=gwas_pine$gcontig__gcontig_pos,
                         X1 = gwas_pine$pine_ind),
              file="flip_pine_v2.cor", row.names=FALSE)
getwd()
