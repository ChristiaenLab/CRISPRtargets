library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)

ky <- makeTxDbFromGFF("HT.Gene.gff3")
tpt <- transcripts(ky)
promoter <- promoters(tpt,1107,107)
gene <- punion(tpt,promoter)
mcols(gene) <- mcols(tpt)

cint <- import('cint.bed')
tmp <- findOverlapPairs(gene,cint)
conserved <- pintersect(tmp)

export(conserved,'conserved.gff3')
