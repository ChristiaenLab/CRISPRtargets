library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Crobusta.HT.KY)
library(optparse)
library(purrr)

source("fns.R")

parser <- OptionParser()
parser <- add_option(parser, '--genes', action = 'store',default = 'genes.txt')
parser <- add_option(parser, '--cutoff', action = 'store',default = 0.5)
opts <- parse_args(parser)

kyid <- readLines(opts$genes)

kh.ky <- read.csv('kydb.csv',row.names=1)

prefix <- paste0(Sys.Date(),'/')
dir.create(prefix)

prefix <- paste0(prefix,'targets')
bed <- paste0(prefix,'.bed')

snps <- import('cint.snps.bed')
snps2 <- lapply(list.files(pattern='C_int_.\\.snps.bed'),import)
snps <- unlist(do.call(GRangesList,append(snps2,snps)))

ky <- makeTxDbFromGFF("HT.Gene.gff3")
tpt <- transcripts(ky)
tpt$KYID <- sub('\\.v.*','',tpt$tx_name)
promoter <- promoters(tpt,1107,107)
gene <- punion(tpt,promoter)
mcols(gene) <- mcols(tpt)

cint <- import('cint.bed')
tmp <- findOverlapPairs(gene,cint)
conserved <- pintersect(tmp)

conslist <- split(conserved,conserved$KYID)
		  

export(unlist(conslist[kyid]),bed)
system2('./flashfry.sh',bed)

out <- readOut(paste0(prefix,'.output'),0)

contig <- getContigs(out,conserved)
targetGene <- onTarget(contig)


ntsymbols <- read.table('ntsymbols.txt',T,row.names=1)
peaks <- import('accessomeKY.bed')

targetSNPs <- findOverlaps(targetGene,snps)
ix <- split(targetSNPs,from(targetSNPs))
# targetSNPs <- split(targetSNPs,paste0(seqnames(first(targetSNPs)),':',start(first(targetSNPs))))

sel <- sapply(ix, function(i) from(i)[1])
targets <- mapply(getSNPs, split(targetGene[sel],1:length(sel)),
		  lapply(ix,function(i) snps[to(i)]))

h <- sapply(targets,seqEntropy)
gap <- sapply(targets,letterFrequency,'-')

targetGene$entropy <- 0
targetGene$gaps <- 0
targetGene$entropy[sel] <- h
targetGene$gaps[sel] <- gap
targetGene$target[sel] <- sapply(targets,as.character)

fwdOligo <- DNAStringSet(sub('^.{1}(.*).{3}$','G\\1',
			     targetGene$target))
revOligo <- reverseComplement(fwdOligo)

fwdOligo <- paste0('agat',fwdOligo)
revOligo <- paste0('aaac',revOligo)
targetGene$fwd <- fwdOligo
targetGene$rev <- revOligo

targetGene$accessibility <- overlapWidth(targetGene,peaks)

targetGene$transcript_overlap <- mapply(function(i,kyid){
			countOverlaps(targetGene[i],
				      gene[gene$KYID==kyid],
				      ignore.strand=T)
		  },1:length(targetGene),targetGene$KYID)

cutoff <- 0.5
targetGene$over_cutoff <- targetGene$Doench2014OnTarget>cutoff

files <- lapply(list.files(pattern='C_int_.\\.bed'),import)

conservation <- sapply(files,overlapWidth,x=targetGene)
targetGene$conservation <- apply(conservation,1,sum)

sortedky <- lapply(split(targetGene,targetGene$KYID),
		 function(x){
			 x <- x[order(x$otCount)]
			 x <- x[order(x$Doench2014OnTarget,decreasing=T)]
			 x <- x[order(x$accessibility,decreasing=T)]
			 x <- x[order(x$transcript_overlap,
				      decreasing=T)]
			 x <- x[order(x$conservation,decreasing=T)]
			 x <- x[order(x$over_cutoff,decreasing=T)]
			 x <- x[order(x$entropy)]
			 x <- x[order(x$gaps)]
			 return(x)
		 })
sorted <- unlist(do.call(GRangesList,sortedky))
export(sorted,paste0(prefix,'.gff3'))
sorted$KH2013 <- kh.ky[sorted$KYID,'KH2013']
sorted$GeneName <- kh.ky[sorted$KYID,'Gene_Name']

cols <- c('KYID','KH2013','GeneName','target', 'fwd','rev', "Doench2014OnTarget",
	  "otCount",'entropy', 'accessibility',
	  'conservation')

write.csv(sorted[,cols],
	  paste0(prefix,'.csv'),row.names=F)
# 
# mapply(function(x,y) write.csv(x[x$gaps==0,c('target','fwd','rev',
#                    "Doench2014OnTarget","otCount",'entropy','accessibility')],
#                                y,row.names=F),
#        sorted,c('tyr_sgrna.csv','depdc_sgrna.csv'))
# 
best3 <- unlist(do.call(GRangesList,
			lapply(sortedky,function(x) x[1:3])))
best3$sequence_name <- paste0(best3$KYID,
			      '_sgRNA',as.character(1:3))

write.csv(best3[,c('sequence_name',cols)],
	  paste0(prefix,'_sgrna.csv'), row.names=F)

