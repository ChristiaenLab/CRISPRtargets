#!/bin/bash
IN=$1
FILE=$(echo $IN | sed 's/\.[A-Za-z]*$//')

#create db
if ! test -f htky_db; then
  mkdir -p build
  java -Xmx4g -jar FlashFry-assembly-1.12.jar \
   index \
   --tmpLocation ./build \
   --database htky_db \
   --reference HT.Ref.fasta \
   --enzyme spcas9ngg
fi

if [[ $IN =~ .*\.bed ]]; then
  bedtools getfasta -name -fi HT.Ref.fasta -fo $FILE.fasta -bed $IN -s
  IN=$FILE.fasta
fi

#get targets
java -Xmx16g -jar FlashFry-assembly-1.12.jar \
 discover \
 --database htky_db \
 --fasta $IN \
 --positionOutput \
 --output $FILE.output

#get scores
java -Xmx16g -jar FlashFry-assembly-1.12.jar \
 score \
 --input $FILE.output \
 --output $FILE.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database htky_db \

