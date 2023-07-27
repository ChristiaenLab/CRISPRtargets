source('fns.R')

system2('Rscript','runFlashfry.R --genes sgrna.csv')

sgrna <- read.csv('sgrnaSiteScore.csv')

dat <- read.csv('tvc514.csv',row.names=1)
dat <- dat[dat$logfoldchanges>1&
	   dat$pvals_adj<0.01& 
	   #dat$SNR>1&
	   !is.na(dat$gene_network_leiden),]
dat <- dat[order(dat$logfoldchanges,decreasing=T),]

sel <- setdiff(dat[dat$gene_network_leiden==2,'Gene.ID'],sgrna$KYID)[1:192]

write.csv(dat[sel,],'genes.csv')
system2('Rscript','runFlashfry.R --genes genes.csv')

targets <- read.csv('genesSiteScore.csv')
targets <- targets[order(targets$logfoldchanges,decreasing=T),]

sel <- intersect(names(targets),names(sgrna))
targets <- targets[,sel]
sgrna <- sgrna[,sel]

targets <- targets[1:(96*3-nrow(sgrna)),]

targets <- rbind(sgrna,targets)

getPlates(targets,'tvc')
