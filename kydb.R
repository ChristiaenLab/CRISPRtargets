library(purrr)
kh.ky1 <- read.table("KHtoKY.tsv")
names(kh.ky1) <- c("KY_ID","KH2013","Transcripts",'KHtranscript')
kh.ky1[,1] <- sub('\\.v.*','',kh.ky1[,3])
kh.ky1[,2] <- sub('\\.v.*','',kh.ky1[,4])
kh.ky2 <- read.table("KH.KYname.tsv",header=T)
names(kh.ky2) <- c("KH2013","KY_ID","KHname")
kh.ky3 <- read.delim("KH.KY.dictionary.txt")

kh.ky <- Reduce(partial(merge,all.x=T,all.y=T),
		list(kh.ky3,kh.ky2,kh.ky1))
kh.ky[is.na(kh.ky)] <- ''

kh.ky <- split(kh.ky[,c("KY_ID","KH2013","KYname","KHname","Gene_Name")],kh.ky$KY_ID)
kh.ky <- do.call(rbind,
		 lapply(kh.ky,sapply,
			compose(partial(sub,'^;',''),partial(sub,';$',''),
				partial(paste,collapse=';'),
				unique)))
write.csv(kh.ky,'kydb.csv')
