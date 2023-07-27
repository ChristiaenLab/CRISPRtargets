library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Crobusta.HT.KY)

source("fns.R")

kyid <- c('KY.Chr12.791','KY.Chr2.2230')
prefix <- 'tyr_depdc'
bed <- paste0(prefix,'.bed')
dat <- read.csv("tvc514.csv")
snps <- import('cint.snps.bed')

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

getSNPs <- function(r,snps){
	reversed <- as.character(strand(r)[1])=='-'
	
	targetSNPs <- split(snps,
			    paste0(seqnames(snps),':',
				   start(snps)))
	target <- DNAString(r$target[1])
	nts <- sapply(targetSNPs,
		      function(x) {
			      paste0(unique(unlist(strsplit(paste0(x$name,
						    collapse=''),''))),
		      collapse='')
		      })
	gap <- grepl('\\.',nts)
	nts <- sub('\\.','',nts)
	subst <- DNAStringSet(mapply(function(gap,nts){
					     if(gap) '-' else{
						     mergeIUPACLetters(nts)
					     }},gap,nts))
	ix <- unique(start(snps))-start(r)[1]+1
	if(reversed) {
		#                 subst <- complement(subst)
		target <- reverseComplement(target)
		ix <- ix
	}
	target[ix] <- unlist(subst)
	if(reversed) target <- reverseComplement(target)
	return(target)
}

seqEntropy <- function(x){
	nn <- c('M','R','W','S','Y','K')
	nnn <- c('V','H','D','B')
	h <- sum(letterFrequency(x,nn))+
		sum(letterFrequency(x,nnn))*-log2(1/3)+
		sum(letterFrequency(x,c('N')))*2
	return(h)
}

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

targetGene$accessibility <- 0
accessible <- findOverlaps(targetGene,peaks)
sel <- unique(from(accessible))
accessibleNTs <- width(pintersect(targetGene[from(accessible)],
				  peaks[to(accessible)]))
accessibleNTs <- sapply(split(accessibleNTs,from(accessible)),sum)
targetGene$accessibility[sel] <- accessibleNTs
targetGene$transcript_overlap <- mapply(function(i,kyid){
			countOverlaps(targetGene[i],
				      gene[gene$KYID==kyid],
				      ignore.strand=T)
		  },1:length(targetGene),targetGene$KYID)

cutoff <- 0.5
targetGene$over_cutoff <- targetGene$Doench2014OnTarget>cutoff

sorted <- lapply(split(targetGene,targetGene$KYID),
		 function(x){
			 x <- x[order(x$otCount)]
			 x <- x[order(x$Doench2014OnTarget,decreasing=T)]
			 x <- x[order(x$accessibility,decreasing=T)]
			 x <- x[order(x$transcript_overlap,
				      decreasing=T)]
			 x <- x[order(x$over_cutoff,decreasing=T)]
			 x <- x[order(x$entropy)]
			 x <- x[order(x$gaps)]
			 return(x)
		 })

mapply(function(x,y) write.csv(x[x$gaps==0,c('target','fwd','rev',
		   "Doench2014OnTarget","otCount",'entropy','accessibility')],
			       y,row.names=F),
       sorted,c('tyr_sgrna.csv','depdc_sgrna.csv'))

best3 <- unlist(do.call(GRangesList,
			lapply(sorted,function(x) x[1:3])))
best3$sequence_name <- paste0(best3$Uniq_Name,
			      '_sgRNA',as.character(1:3))

write.csv(best3[,c('sequence_name','target','fwd','rev',
		   "Doench2014OnTarget","otCount")],
	  'tyr_depdc_sgrna.csv', row.names=F)
