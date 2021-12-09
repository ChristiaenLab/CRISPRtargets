library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Crobusta.HT.KY)
library(optparse)

source('fns.R')

parser <- OptionParser()
parser <- add_option(parser, '--genes', action = 'store',default = 'sgrna.csv')
parser <- add_option(parser, '--cutoff', action = 'store',default = 0)
opts <- parse_args(parser)

prefix <- sub('\\..*','',opts$genes)
bed <- paste0(prefix,'.bed')
dat <- read.csv(opts$genes)
sel <- dat[,1]

exon <- getExons()

cint <- import('cint.bed')
# snps <- import('cint.snps.bed')

peaks <- import('accessomeKY.bed')

cint <- intersect(cint,peaks)
tmp <- findOverlapPairs(exon,cint)

conserved <- GRanges(
	seqnames(first(tmp)),
	IRanges(
		mapply(max,start(first(tmp)),start(second(tmp))),
		mapply(min,end(first(tmp)),end(second(tmp)))
	),
	strand(first(tmp)),
	mcols(first(tmp))
)
names(conserved) <- paste0(conserved$tx_name,'_exon',as.character(conserved$exon_rank))

export(conserved,'genometargets.bed')

targetSites <- sapply(paste0(sel,'\\.v'),grep,names(conserved))

export(conserved[unlist(targetSites)],bed)

system2('./flashfry.sh',bed)

out <- readOut(paste0(prefix,'.output'),dat,opts$cutoff)

contig <- getContigs(out)
targetGene <- getOnTarget(contig, out,exon,prefix)
