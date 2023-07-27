nucmer -p sgrna_HTKY -t 6 conserved_sgrna.fasta HT.Ref.fasta 
dnadiff -p sgrna_HTKY -d sgrna_HTKY.delta
nucmer -p sgrna_roscoff -t 6 conserved_sgrna.fasta ncbi-genomes-2021-10-13/GCA_018327825.1_Cint_typeB-Roscoff_1.0_genomic.fna
dnadiff -p sgrna_roscoff -d sgrna_roscoff.delta

nucmer -p sgrna_plymouth -t 6 conserved_sgrna.fasta ncbi-genomes-2021-10-13/GCA_018327805.1_Cint_typeB-Plymouth_1.0_genomic.fna
dnadiff -p sgrna_plymouth -d sgrna_plymouth.delta

blastn -query conserved_sgrna.fasta -subject HT.Ref.fasta -word_size 7 -penalty -1 -reward 1 -gapopen 3 -gapextend 2 > sgrna_aligned_HTKY.txt
blastn -query conserved_sgrna.fasta -subject ncbi-genomes-2021-10-13/GCA_018327825.1_Cint_typeB-Roscoff_1.0_genomic.fna -word_size 7 -penalty -1 -reward 1 -gapopen 3 -gapextend 2 > sgrna_aligned_roscoff.txt
blastn -query conserved_sgrna.fasta -subject ncbi-genomes-2021-10-13/GCA_018327805.1_Cint_typeB-Plymouth_1.0_genomic.fna -word_size 7 -penalty -1 -reward 1 -gapopen 3 -gapextend 2 > sgrna_aligned_plymouth.txt
