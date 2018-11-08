#!/bin/bash

WD=/scratch/rtraborn/tsrchitect-figures
At_genome=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_genome.fasta
At_genome_file=TAIR10.genome

#AtPEAT
AtPEATinr="AtPEAT_inr.fasta"
AtPEAT_TATA="AtPEAT_TATA.fasta"
#AtCAGE
AtCAGEinr="AtCAGE_inr.fasta"
AtCAGE_TATA="AtCAGE_TATA.fasta"
#AtOligo
AtOligoinr="AtOligo_inr.fasta"
AtOligo_TATA="AtOligo_TATA.fasta"

cd $WD

echo "Creating TATA background"
bedtools random -l 16 -n 1000000 -g $At_genome_file > At_TATA_bg.bed
bedtools getfasta -fi $At_genome -bed At_TATA_bg.bed -fo At_TATA_bg.fasta

echo "Creating Inr background"
bedtools random -l 10 -n 1000000 -g $At_genome_file > At_Inr_bg.bed
bedtools getfasta -fi $At_genome -bed At_Inr_bg.bed -fo At_Inr_bg.fasta

echo "Finding motifs AtPEAT"
findMotifs.pl $AtPEATinr fasta AtPEATinr_results -fasta At_Inr_bg.fasta -len 6 -p 4
findMotifs.pl $AtPEAT_TATA fasta AtPEAT_TATA_results -fasta At_TATA_bg.fasta -len 6 -p 4

echo "Finding motifs AtOligo"
findMotifs.pl $AtCAGEinr fasta AtCAGEinr_results -fasta At_Inr_bg.fasta -len 6 -p 4
findMotifs.pl $AtCAGE_TATA fasta AtCAGE_TATA_results -fasta At_TATA_bg.fasta -len 6 -p 4

echo "Finding motifs AtCAGE"
findMotifs.pl $AtOligoinr fasta AtOligoinr_results -fasta At_Inr_bg.fasta -len 6 -p 4
findMotifs.pl $AtOligo_TATA fasta AtOligo_TATA_results -fasta At_TATA_bg.fasta -len 6 -p 4

echo "Job Complete!"
