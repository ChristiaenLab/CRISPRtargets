#blastn -query sgrnaTargetExamples/targets.fasta -subject BNJZ01.1.fsa_nt -out tmp.txt -outfmt 7 -word_size 23

#awk '/>/ {sub(">","",$1); sub(",","",$5);print $1"\t"$5}' ncbi-genomes-2021-10-13/GCA_009617815.1_Ciona_intestinalis_HT_Hoya_T-line_assembly_2019_genomic.fna > htky.txt

nucmer -p roscoff -t 6 HT.Ref.fasta ncbi-genomes-2021-10-13/GCA_018327825.1_Cint_typeB-Roscoff_1.0_genomic.fna
dnadiff -p roscoff -d roscoff.delta

nucmer -p plymouth -t 6 HT.Ref.fasta ncbi-genomes-2021-10-13/GCA_018327805.1_Cint_typeB-Plymouth_1.0_genomic.fna
dnadiff -p plymouth -d plymouth.delta

awk '{print $12"\t"$1"\t"$2"\t*"}' roscoff.mcoords > roscoff.bed
awk '{print $12"\t"$1"\t"$2"\t*"}' plymouth.mcoords > plymouth.bed

bedtools intersect -a roscoff.bed -b plymouth.bed | bedtools sort > tmp.bed
bedtools merge -i tmp.bed > cint.bed
rm tmp.bed

bedtools intersect -a HT.Gene.gff3 -b cint.bed > conserved.bed

awk '{print $11"\t"$1-1"\t"$1"\t*"}' roscoff.snps > roscoff.snps.bed
bedtools merge -i roscoff.snps.bed > tmp.bed
mv tmp.bed roscoff.snps.bed

awk '{print $11"\t"$1-1"\t"$1"\t*"}' plymouth.snps > plymouth.snps.bed
bedtools merge -i plymouth.snps.bed > tmp.bed
mv tmp.bed plymouth.snps.bed

cat plymouth.snps.bed roscoff.snps.bed > cint.snps.bed
bedtools sort -i cint.snps.bed > tmp.bed
bedtools merge -i tmp.bed > cint.snps.bed
rm tmp.bed

