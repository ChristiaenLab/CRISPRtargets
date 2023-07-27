library(GenomicRanges)
library(rtracklayer)

cint <- import('cint.bed')
files <- lapply(list.files(pattern='C_int_.\\.bed'),import)
conserved2 <- Reduce(function(x,y) reduce(pintersect(findOverlapPairs(x,y))),files,cint)

snps2 <- lapply(list.files(pattern='C_int_.\\.snps.bed'),import)
