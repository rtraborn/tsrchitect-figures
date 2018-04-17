#!/bin/bash

WD=/projects/TSRplantsV2.0/At/AtCoveragePlots
At_genes=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_GFF3_genes.gff
At_genome=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_genome.fasta
AtPEAT="AtPEAT_TSRset-1.bed"
AtCAGE="AtCAGE_TSRset-1.bed"
AtOligo="AtOligo_TSRset-1.bed"

cd $WD

echo "Creating the TAIR10 .genome file"

samtools faidx $At_genome
cut -f1,2 ${At_genome}.fai > TAIR10.genome

echo "Creating intervals reflecting the desired genomic segments"

#gt gff3 -tidy -addintrons /projects/TSRplantsV2.0/At/AtGENOME/TAIR10_GFF3_genes.gff > At_genes_new.gff3

#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="CDS" {print }' | sortBed > At_cds.gff3
#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="intron" {print }' | sortBed > At_introns.gff3
#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="gene" {print }' | sortBed > At_genes.gff3
#bedtools flank -i At_genes.gff3 -g TAIR10.genome -l 500 -r 0 -s > At_upstream.gff3
#bedtools flank -i At_genes.gff3 -g TAIR10.genome -l 0 -r 500 -s > At_downstream.gff3
#cat At_genes.gff3 At_upstream.gff3 At_downstream.gff3 | sortBed > At_combined.gff3
#bedtools complement -i At_combined.gff3 -g TAIR10.genome > At_intergenic.bed

echo "Performing coverage analysis of all TSRs"

echo "AtPEAT"
bedtools coverage -s -counts -a $AtPEAT -b At_upstream.gff3 > AtPEAT_upstream.tsr
bedtools coverage -s -counts -a $AtPEAT -b At_downstream.gff3 > AtPEAT_downstream.tsr
bedtools coverage -s -counts -a $AtPEAT -b At_cds.gff3 > AtPEAT_cds.tsr
bedtools coverage -s -counts -a $AtPEAT -b At_introns.gff3 > AtPEAT_introns.tsr
bedtools coverage -counts -a $AtPEAT -b At_intergenic.bed > AtPEAT_intergenic.tsr

echo "AtOligo"
bedtools coverage -s -counts -a $AtOligo -b At_upstream.gff3 > AtOligo_upstream.tsr
bedtools coverage -s -counts -a $AtOligo -b At_downstream.gff3 > AtOligo_downstream.tsr
bedtools coverage -s -counts -a $AtOligo -b At_cds.gff3 > AtOligo_cds.tsr
bedtools coverage -s -counts -a $AtOligo -b At_introns.gff3 > AtOligo_introns.tsr
bedtools coverage -counts -a $AtOligo -b At_intergenic.bed > AtOligo_intergenic.tsr

echo "AtCAGE"
bedtools coverage -s -counts -a $AtCAGE -b At_upstream.gff3 > AtCAGE_upstream.tsr
bedtools coverage -s -counts -a $AtCAGE -b At_downstream.gff3 > AtCAGE_downstream.tsr
bedtools coverage -s -counts -a $AtCAGE -b At_cds.gff3 > AtCAGE_cds.tsr
bedtools coverage -s -counts -a $AtCAGE -b At_introns.gff3 > AtCAGE_introns.tsr
bedtools coverage -counts -a $AtCAGE -b At_intergenic.bed > AtCAGE_intergenic.tsr

echo "Job Complete!"
