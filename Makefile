FASTA = HT.Ref.fasta
URL = http://ghost.zool.kyoto-u.ac.jp/datas/
GFF = HT.Gene.gff
ZIP = HT.KYGene.gff3.zip
FF = FlashFry-assembly-1.12.jar
DB = htky_db

$(DB): $(FASTA) $(FF)
	mkdir -p build
	java -Xmx4g -jar $(FF) \
	 index \
	 --tmpLocation ./build \
	 --database $(DB) \
	 --reference $(FASTA) \
	 --enzyme spcas9ngg

$(FF):
	wget https://github.com/mckennalab/FlashFry/releases/download/1.12/$(FF)

$(FASTA):
	wget -U firefox $(URL)$(FASTA).zip
	unzip -o $(FASTA).zip

$(GFF):
	wget -U firefox $(URL)$(ZIP)
	unzip -o $(ZIP)

cint.snps.bed:
	tar -zxvf cint.snps.bed.tar.gz

clean:
	rm -f $(ZIP)
	rm $(FASTA).zip
	rm $(DB)*
