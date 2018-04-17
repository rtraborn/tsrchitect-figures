#!/bin/bash

WD=/projects/TSRplantsV2.0/At/AtCoveragePlots
At_genes=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_GFF3_genes.gff
At_genome=/projects/TSRplantsV2.0/At/AtGENOME/TAIR10_genome.fasta
AtPEAT_pos=/projects/TSRplantsV2.0/At/AtPEAT/tsr/AtPEAT_tss.pos 
AtCAGE_pos=/projects/TSRplantsV2.0/At/AtCAGE/tsr/AtCAGE_tss.pos
AtOligo_pos=/projects/TSRplantsV2.0/At/AtOligo/tsr/AtOligo_tss.pos
AtOligo_bed=/projects/TSRplantsV2.0/At/AtCoveragePlots/AtOligotss_corr.bed
AtPEAT_thresh=/projects/TSRplantsV2.0/At/AtDistPlots/pos_files/AtPEATtss_filtered.pos
AtCAGE_thresh=/projects/TSRplantsV2.0/At/AtDistPlots/pos_files/AtCAGEtss_filtered.pos
AtOligo_thresh=/projects/TSRplantsV2.0/At/AtDistPlots/pos_files/AtOligo_tss_filtered.pos 

cd $WD

echo "Creating the TAIR10 .genome file"

samtools faidx $At_genome
cut -f1,2 ${At_genome}.fai > TAIR10.genome

#echo "Converting the .pos formatted TSSs to BED format"

#pos2bed.pl $AtPEAT_pos > AtPEATtss.bed
#pos2bed.pl $AtCAGE_pos > AtCAGEtss.bed
#pos2bed.pl $AtOligo_pos > AtOligotss.bed
#pos2bed.pl $AtPEAT_thresh > AtPEATtss_t.bed
#pos2bed.pl $AtCAGE_thresh > AtCAGEtss_t.bed
pos2bed.pl $AtOligo_thresh > AtOligotss_t.bed

echo "Creating intervals reflecting the desired genomic segments"

#gt gff3 -tidy -addintrons /projects/TSRplantsV2.0/At/AtGENOME/TAIR10_GFF3_genes.gff > At_genes_new.gff3

#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="CDS" {print }' | sortBed > At_cds.gff3
#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="intron" {print }' | sortBed > At_introns.gff3
#cat At_genes_new.gff3 | awk 'BEGIN{OFS="\t";} $3=="gene" {print }' | sortBed > At_genes.gff3
#bedtools flank -i At_genes.gff3 -g TAIR10.genome -l 500 -r 0 -s > At_upstream.gff3
#bedtools flank -i At_genes.gff3 -g TAIR10.genome -l 0 -r 500 -s > At_downstream.gff3
#cat At_genes.gff3 At_upstream.gff3 At_downstream.gff3 | sortBed > At_combined.gff3
#bedtools complement -i At_combined.gff3 -g TAIR10.genome > At_intergenic.bed

echo "Performing coverage analysis of all TSSs"

echo "AtPEAT"
#bedtools coverage -s -counts -a AtPEATtss.bed -b At_upstream.gff3 > AtPEAT_upstream.coverage
#bedtools coverage -s -counts -a AtPEATtss.bed -b At_downstream.gff3 > AtPEAT_downstream.coverage
#bedtools coverage -s -counts -a AtPEATtss.bed -b At_cds.gff3 > AtPEAT_cds.coverage
#bedtools coverage -s -counts -a AtPEATtss.bed -b At_introns.gff3 > AtPEAT_introns.coverage
#bedtools coverage -counts -a AtPEATtss.bed -b At_intergenic.bed > AtPEAT_intergenic.coverage

echo "AtOligo"
bedtools coverage -s -counts -a $AtOligo_bed -b At_upstream.gff3 > AtOligo_upstream.coverage
bedtools coverage -s -counts -a $AtOligo_bed -b At_downstream.gff3 > AtOligo_downstream.coverage
bedtools coverage -s -counts -a $AtOligo_bed -b At_cds.gff3 > AtOligo_cds.coverage
bedtools coverage -s -counts -a $AtOligo_bed -b At_introns.gff3 > AtOligo_introns.coverage
bedtools coverage -counts -a $AtOligo_bed -b At_intergenic.bed > AtOligo_intergenic.coverage

echo "AtCAGE"
#bedtools coverage -s -counts -a AtCAGEtss.bed -b At_upstream.gff3 > AtCAGE_upstream.coverage
#bedtools coverage -s -counts -a AtCAGEtss.bed -b At_downstream.gff3 > AtCAGE_downstream.coverage
#bedtools coverage -s -counts -a AtCAGEtss.bed -b At_cds.gff3 > AtCAGE_cds.coverage
#bedtools coverage -s -counts -a AtCAGEtss.bed -b At_introns.gff3 > AtCAGE_introns.coverage
#bedtools coverage -counts -a AtCAGEtss.bed -b At_intergenic.bed > AtCAGE_intergenic.coverage

echo "Performing coverage analysis of all above-threshold TSSs"

echo "AtPEAT"
#bedtools coverage -s -counts -a AtPEATtss_t.bed -b At_upstream.gff3 > AtPEAT_t_upstream.coverage
#bedtools coverage -s -counts -a AtPEATtss_t.bed -b At_downstream.gff3 > AtPEAT_t_downstream.coverage
#bedtools coverage -s -counts -a AtPEATtss_t.bed -b At_cds.gff3 > AtPEAT_t_cds.coverage
#bedtools coverage -s -counts -a AtPEATtss_t.bed -b At_introns.gff3 > AtPEAT_t_introns.coverage
#bedtools coverage -counts -a AtPEATtss_t.bed -b At_intergenic.bed > AtPEAT_t_intergenic.coverage

echo "AtOligo"
bedtools coverage -s -counts -a AtOligotss_t.bed -b At_upstream.gff3 > AtOligo_t_upstream.coverage
bedtools coverage -s -counts -a AtOligotss_t.bed -b At_downstream.gff3 > AtOligo_t_downstream.coverage
bedtools coverage -s -counts -a AtOligotss_t.bed -b At_cds.gff3 > AtOligo_t_cds.coverage
bedtools coverage -s -counts -a AtOligotss_t.bed -b At_introns.gff3 > AtOligo_t_introns.coverage
bedtools coverage -counts -a AtOligotss_t.bed -b At_intergenic.bed > AtOligo_t_intergenic.coverage

echo "AtCAGE"
#bedtools coverage -s -counts -a AtCAGEtss_t.bed -b At_upstream.gff3 > AtCAGE_t_upstream.coverage
#bedtools coverage -s -counts -a AtCAGEtss_t.bed -b At_downstream.gff3 > AtCAGE_t_downstream.coverage
#bedtools coverage -s -counts -a AtCAGEtss_t.bed -b At_cds.gff3 > AtCAGE_t_cds.coverage
#bedtools coverage -s -counts -a AtCAGEtss_t.bed -b At_introns.gff3 > AtCAGE_t_introns.coverage
#bedtools coverage -counts -a AtCAGEtss_t.bed -b At_intergenic.bed > AtCAGE_t_intergenic.coverage

echo "Job Complete!"
