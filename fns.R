getExons <- function(gff='HT.Gene.gff3') {
	require(GenomicRanges)
	require(rtracklayer)
	require(Biostrings)
	require(GenomicFeatures)
	require(BSgenome.Crobusta.HT.KY)

	ky <- makeTxDbFromGFF(gff)

	tpt <- transcripts(ky)
	promoter <- promoters(tpt,1107,107)
	promoter$exon_rank <- 0
	names(promoter) <- promoter$tx_name

	exon <- exonsBy(ky,'tx')
	names(exon) <- tpt$tx_name
	exon <- unlist(exon)
	exon$tx_name <- names(exon)

	exon <- c(promoter,exon)
	exon$KYID <- sub('\\.v.*','',exon$tx_name)
	exon
}

#split matches for each sequence
parseTarget <- function(x) {
	require(GenomicRanges)
	x <- sub('\\|.*','',x)
	seqnames <- sub('.*<(.*):.*','\\1',x)
	strand <- sapply(grepl('\\^R',x), function(x) if(x) '-' else '+')
	coord <- as.numeric(sub('.*:([0-9]+).*','\\1',x))
	ranges <- IRanges(coord+1,coord+23)
	mismatches <- sub('.*_[0-9]+_([0-9]+)<.*','\\1',x)
	sequence <- sub('_.*','',x)
	return(GRanges(
		       seqnames,ranges,strand,
		       mismatches=as.numeric(mismatches),
		       sequence=sequence))
}

readOut <- function(x,cutoff=0.50){
	out <- read.delim(x)
	outscore <- read.delim(paste0(x,'.scored'))
	outscore <- outscore[
	     is.finite(outscore$Doench2014OnTarget)&
	     outscore$Doench2014OnTarget>cutoff&
	     outscore$dangerous_polyT=="NONE"&
	     outscore$dangerous_GC=="NONE"&
	     outscore$dangerous_in_genome=="IN_GENOME=1",
	]

	out <- merge(out,outscore)
	#         out$KYID=sub('\\.v.*','',out$contig)
	out$KYID=sub('::.*','',out$contig)

	#         out <- merge(out,dat,by.x="KYID",by.y="Gene.ID",all.x=T)
	#         out[is.na(out$Uniq_Name),"Uniq_Name"] <- paste0(
	#                 "KY2019:",out$KYID[is.na(out$Uniq_Name)]
	#         )

	return(out)
}

parseOut <- function(out){
	require(GenomicRanges)
	offTargets <- strsplit(out$offTargets,',')
	offTargets <- do.call(
		GRangesList,
		lapply(offTargets, parseTarget)
	)
	return(offTargets)
}

#select 3 best targets for each gene
best3 <- function(x,cutoff=0.55) {
	tmp <- x[x$Doench2014OnTarget>cutoff]
	if(length(tmp)<3) {
		x <- x[order(x$Doench2014OnTarget,decreasing=T)[1:3]]
	} else x <- tmp
	if(sum(x$conserved)>2){
		x <- x[x$conserved,]
	}
	x <- x[order(x$exon)]
	x <- x[order(x$otCount)]
	x <- x[order(x$SNPct)]
	x <- x[order(x$exon_overlap_ct,decreasing=T)]
	return(x[1:min(length(x),3)])
}

countSNPs <- function(ranges, snps='cint.snps.bed'){
	require(GenomicRanges)
	require(rtracklayer)
	snps <- import('cint.snps.bed')
	ranges$SNPct <- countOverlaps(ranges,snps)
	ranges
}

countTx <- function(targetGene,exon) lapply(
	names(targetGene),
	function(x){
		tmp <- targetGene[[x]]
		tmp$exon_overlap_ct <- countOverlaps(
			tmp,exon[exon$KYID==x],ignore.strand=T
		)
		tmp
	}
)


getContigs <- function(out,conserved){
	require(GenomicRanges)
	contig <- GRanges(
		sub('.*::([A-Za-z0-9]+):.*','\\1',out$contig),
		IRanges(
			as.numeric(sub(
				'.*:([0-9]+)-.*','\\1',out$contig
			))+1,
			as.numeric(sub(
				'.*-([0-9]+)\\(.*','\\1',out$contig
			))
		),
		sub('.*\\(([+-])\\)$','\\1',out$contig),
		#KYID=sub('\\.v.*','',out$contig),
		contig=sub('(.*)::.*','\\1',out$contig)
	)
	mcols(contig) <- cbind(mcols(contig),out[,-1:-3])
	contig$conserved <- overlapsAny(contig,conserved)
	contig
}

onTarget <- function(contig){
	require(GenomicRanges)
	require(Biostrings)
	targetGene <- GRanges(seqnames(contig),
		do.call(IRanges,as.data.frame(t(mapply(
			function(x.start,x.stop,y.start,y.stop,is.fwd){
				if(is.fwd) c(start=x.start+y.start,
					     end=x.start+y.stop-1) 
				else c(start=x.stop-y.stop+1,
				       end=x.stop-y.start)
			},
			start(contig), end(contig), out$start,out$stop,
			as.character(strand(contig))=='+'
		)))),
		mapply(
			function(x,y) {
				if(y) x else {
					if(x=='+') {
						'-' 
					}else '+'
				}
			},
			as.character(strand(contig)),
			out$orientation=="FWD"
		),
		mcols(contig)
	)

	targetGene <- targetGene[!duplicated(targetGene)]
	return(targetGene)
}

getOnTarget <- function(contig,out,exon,file, snps='cint.snps.bed'){
	require(GenomicRanges)
	require(Biostrings)
	onTarget <- GRanges(
		seqnames(contig),
		do.call(IRanges,as.data.frame(t(mapply(
			function(x.start,x.stop,y.start,y.stop,is.fwd){
				if(is.fwd) c(start=x.start+y.start,
					     end=x.start+y.stop-1) 
				else c(start=x.stop-y.stop+1,
				       end=x.stop-y.start)
			},
			start(contig), end(contig), out$start,out$stop,
			as.character(strand(contig))=='+'
		)))),
		mapply(
			function(x,y) {
				if(y) x else if(x=='+') '-' else '+'
			},
			as.character(strand(contig)),
			out$orientation=="FWD"
		),
		mcols(contig)
	)

	onTarget$exon <- as.numeric(sub('.*exon([0-9A-Z]+)::.*','\\1',onTarget$contig))
	onTarget$exon[is.na(onTarget$exon)] <- 0

	onTarget <- onTarget[!duplicated(onTarget)]
	onTarget <- countSNPs(onTarget)

	#split seqs by gene
	targetGene <- split(onTarget,onTarget$KYID)
	targetGene <- targetGene[sapply(targetGene,length)>2]
	targetGene <- countTx(targetGene,exon)

	#select 3 best seqs for each gene
	targetGene <- lapply(targetGene,best3)

	targetGene <- unlist(do.call(GRangesList,targetGene))

	fwdOligo <- DNAStringSet(sub('^.{1}(.*).{3}$','G\\1',
				     targetGene$target))
	revOligo <- reverseComplement(fwdOligo)

	fwdOligo <- paste0('agat',fwdOligo)
	revOligo <- paste0('aaac',revOligo)
	targetGene$fwd <- fwdOligo
	targetGene$rev <- revOligo

	targetGene$sequence_name <- paste0(
		targetGene$Uniq_Name,'_sgRNA',as.character(1:3)
	)

	write.csv(mcols(targetGene),
		  paste0(file,'SiteScore.csv'),
		  row.names=F)

	targetGene
}

getPlates <- function(targetGene,out){
	targetGene$plate <- 1:3
	well_position <- c(sapply(
		1:12, function(x) paste0(
			sapply(LETTERS[1:8],rep,3), 
			as.character(x)
		)
	))
	well_position <- factor(
		well_position,
		unique(well_position)
	)
	targetGene$well_position <- well_position[1:length(targetGene)]

	plates <- targetGene[,c(
		'plate','well_position','sequence_name','fwd','rev'
	)]
	plates <- reshape2::melt(
		as.data.frame(mcols(plates)),
		1:3,4:5,
		value.name='sequence'
	)
	plates$sequence_name <- paste0(
		plates$sequence_name,'_',
		plates$variable
	)
	plates <- plates[order(plates$well_position),]
	plates <- split(plates[,c(2,3,5)],plates$plate)

	mapply(write.csv,
	       plates,
	       paste0(out,'Plate',as.character(1:3),'.csv'),
	       row.names=F,quote=F)

}


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


overlapWidth <- function(x,y){
	res <- rep(0,length(x))

	overlaps <- findOverlaps(x,y)
	sel <- unique(from(overlaps))
	nts <- width(pintersect(x[from(overlaps)],
					  y[to(overlaps)]))
	nts <- sapply(split(nts,from(overlaps)),sum)
	res[sel] <- accessibleNTs
	return(res)
}

