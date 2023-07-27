awk '{print $11"\t"$1-1"\t"$1"\t"$2$3"\t0\t*"}' roscoff.snps > roscoff.snps.bed
#bedtools merge -i roscoff.snps.bed > tmp.bed
#mv tmp.bed roscoff.snps.bed

awk '{print $11"\t"$1-1"\t"$1"\t"$2$3"\t0\t*"}' plymouth.snps > plymouth.snps.bed
#bedtools merge -i plymouth.snps.bed > tmp.bed
#mv tmp.bed plymouth.snps.bed

cat plymouth.snps.bed roscoff.snps.bed > cint.snps.bed
bedtools sort -i cint.snps.bed > tmp.bed
#bedtools merge -i tmp.bed > cint.snps.bed
#rm tmp.bed
mv cint.snps.bed



