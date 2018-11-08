#!/bin/bash

WD=/scratch/rtraborn/tsrchitect-figures
At_genome=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_genome.fasta
AtPEAT="AtPEAT_TSSset-1.bed"
AtCAGE="AtCAGE_TSSset-1.bed"
AtOligo="AtOligo_TSSset-1.bed"

AtPEATfilt="AtPEAT_TSSset-1_filt.bed"
AtCAGEfilt="AtCAGE_TSSset-1_filt.bed"
AtOligofilt="AtOligo_TSSset-1_filt.bed"

nTAGs_threshold=20

cd $WD

echo "Creating the TAIR10 .genome file"

samtools faidx $At_genome
cut -f1,2 ${At_genome}.fai > TAIR10.genome

awk '{ if ($5 >= 10) { print } }' $AtPEAT > $AtPEATfilt
awk '{ if ($5 >= 10) { print } }' $AtCAGE > $AtCAGEfilt
awk '{ if ($5 >= 10) { print } }' $AtOligo > $AtOligofilt

echo "Creating intervals reflecting the desired genomic segments"

bedtools slop -i $AtPEATfilt -g TAIR10.genome -l 5 -r 5 -s > AtPEAT_inr.bed

bedtools flank -i $AtPEATfilt -g TAIR10.genome -l 40 -r 0 -s > AtPEAT_TATAi.bed
bedtools flank -i $AtPEATfilt -g TAIR10.genome -l 15 -r 0 -s > AtPEAT_TATAii.bed
bedtools subtract -s -a AtPEAT_TATAi.bed -b AtPEAT_TATAii.bed > AtPEAT_TATA.bed

bedtools slop -i $AtCAGEfilt -g TAIR10.genome -l 5 -r 5 -s > AtCAGE_inr.bed

bedtools flank -i $AtCAGEfilt -g TAIR10.genome -l 40 -r 0 -s > AtCAGE_TATAi.bed
bedtools flank -i $AtCAGEfilt -g TAIR10.genome -l 15 -r 0 -s > AtCAGE_TATAii.bed
bedtools subtract -s -a AtCAGE_TATAi.bed -b AtCAGE_TATAii.bed > AtCAGE_TATA.bed

bedtools slop -i $AtOligofilt -g TAIR10.genome -l 5 -r 5 -s > AtOligo_inr.bed

bedtools flank -i $AtOligofilt -g TAIR10.genome -l 40 -r 0 -s > AtOligo_TATAi.bed
bedtools flank -i $AtOligofilt -g TAIR10.genome -l 15 -r 0 -s > AtOligo_TATAii.bed
bedtools subtract -s -a AtOligo_TATAi.bed -b AtOligo_TATAii.bed > AtOligo_TATA.bed

echo "Extracting sequence data from each bed file"

echo "AtPEAT"
bedtools getfasta -fi $At_genome -bed AtPEAT_inr.bed -fo AtPEAT_inr.fasta
bedtools getfasta -fi $At_genome -bed AtPEAT_TATA.bed -fo AtPEAT_TATA.fasta

echo "AtOligo"
bedtools getfasta -fi $At_genome -bed AtCAGE_inr.bed -fo AtCAGE_inr.fasta
bedtools getfasta -fi $At_genome -bed AtCAGE_TATA.bed -fo AtCAGE_TATA.fasta

echo "AtCAGE"
bedtools getfasta -fi $At_genome -bed AtOligo_inr.bed -fo AtOligo_inr.fasta
bedtools getfasta -fi $At_genome -bed AtOligo_TATA.bed -fo AtOligo_TATA.fasta

echo "Job Complete!"
