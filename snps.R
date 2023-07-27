library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(parallel)

snps <- import('cint.snps.bed')
snps <- split(snps,paste0(seqnames(snps),':',start(snps)))

ncore <- detectCores()-2
cl <- makeCluster(ncore,"FORK")
res <- parSapply(cl,snps,
		 function(x) paste0(x$name,collapse=''))
stopCluster(cl)
