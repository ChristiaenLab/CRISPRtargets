#!/usr/bin/bash

run_nucmer () {
	PREFIX=$1
	FILE=$2
	nucmer -p $PREFIX HT.Ref.fasta $FILE
	dnadiff -p $PREFIX -d $PREFIX.delta
	awk '{print $12"\t"$1"\t"$2"\t*"}' $PREFIX.mcoords > $PREFIX.bed
	awk '{print $11"\t"$1-1"\t"$1"\t"$2$3"\t0\t*"}' $PREFIX.snps > $PREFIX.snps.bed
}

merge_bed () {
	A=$1
	B=$2
	OUT=$3
	bedtools intersect -a $A.bed -b $B.bed | bedtools sort > tmp.bed
	bedtools merge -i tmp.bed > $OUT.bed
	rm tmp.bed
}

while read l; do run_nucmer $l; done < genomes.txt
