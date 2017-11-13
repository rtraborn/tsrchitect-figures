## Objective: calculate the coverage of ESTs and TSRs calculated from ESTs relative to gene annotations

CAGE_dir=/projects/TSRplants/ZmCAGE/tsrOut
CAGE_m=/projects/TSRplants/ZmCAGE/tsrOut/TSRsetCombined.bed

ZmGenbankEST=/projects/TSRplants/ESTcDNA/Z_mays/tsrOut/TSRset-1_reheader.bed

outputDir=/scratch/rtraborn/tsrchitect-figures/output_files
Zm_gff=/projects/TSRplants/ZmCAGE/annotation/Zea_mays.AGPv3.26.gff3
Zm_bed=/projects/TSRplants/ZmCAGE/annotation/Zea_mays.AGPv3.26.bed
Zm_genes=/projects/TSRplants/ZmCAGE/annotation/Zm_genes.bed
Zm_genes_nuc=/projects/TSRplants/ZmCAGE/annotation/Zm_genes_nuc.bed
Zm_output_dir=/projects/TSRplants/ESTcDNA/Z_mays/tsrOut/

TAIR10_genome=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10.genome
TAIR10_fasta=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genome.fasta

PEAT_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT
PEAT_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1_sum.bed
PEAT_tss=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1.txt
PEATbam=/projects/TSRplants/AtPEAT/alignment/AtPEAT.bam
PEAT_out=PEAT_TSSset-1.bed

TokCAGE_tss=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1.txt
TokCAGE_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1_sum.bed
TokCAGEbam=/projects/TSRplants/AtCAGE/alignment/AtCAGE.aligned.bam
TokCAGE_out=TokCAGE_TSSset-1.bed

TokVecbam=/projects/TSRplants/AtOligo/alignment/OligoCap_align.bam
TokVec_tss=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1.txt
TokVec_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1_sum.bed
TokVec_out=TokVec_TSSset-1.bed

AtGenbankEST=/projects/TSRplants/ESTcDNA/A_thaliana/tsrOut/AtEST_genbank_out/TSSset-1.txt
AtEST_out=AtEST_TSSset-1.bed

At_GROseqBED=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed
At_GROseq=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/GROseq_tagDirectory/tss.txt

At_TAIR_gff=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_GFF3_genes.gff
At_TAIR_genes=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genes.bed
At_TAIR_genes_nuc=/projects/TSRplants/ESTcDNA/A_thaliana/TAIR10_genes_nuc.bed

TSS_threshold=5

echo "Converting the TSS files to bed format."

#cat $PEAT_tss | awk -F'\t' -v OFS="\t" '{ if($3=="+") print $1, $2, $2+1, ".", $4, $3 }' | sed '1d' > outfile.p
#cat $PEAT_tss | awk -F'\t' -v OFS="\t" '{ if($3=="-") print $1, $2-1, $2, ".", $4, $3 }' | sed '1d' > outfile.m
#cat outfile.p outfile.m | sortBed -chrThenSizeA -i - > $PEAT_out

#awk '(NR>1) && ($1 != 'chloroplast') ' $TokCAGE_tss > TokCAGE_TSSset-1.txt

#cat TokCAGE_TSSset-1.txt | awk -F'\t' -v OFS="\t" '{ if($3=="+") print $1, $2, $2+1, ".", $4, $3 }' | sed '1d' > outfile.p
#cat TokCAGE_TSSset-1.txt | awk -F'\t' -v OFS="\t" '{ if($3=="-") print $1, $2-1, $2, ".", $4, $3 }' | sed '1d' > outfile.m
#cat outfile.p outfile.m | bedtools sort -chrThenSizeA -i - > $TokCAGE_out
#cat outfile.p outfile.m > $TokCAGE_out

#cat $TokVec_tss | awk -F'\t' -v OFS="\t" '{ if($3=="+") print $1, $2, $2+1, ".", $4, $3 }' | sed '1d' > outfile.p
#cat $TokVec_tss | awk -F'\t' -v OFS="\t" '{ if($3=="-") print $1, $2-1, $2, ".", $4, $3 }' | sed '1d' > outfile.m
#cat outfile.p outfile.m | sortBed -chrThenSizeA -i - > $TokVec_out

#cat $AtGenbankEST | awk -F'\t' -v OFS="\t" '{ if($3=="+") print $1, $2, $2+1, ".", $4, $3 }' | sed '1d' > outfile.p
#cat $AtGenbankEST | awk -F'\t' -v OFS="\t" '{ if($3=="-") print $1, $2-1, $2, ".", $4, $3 }' | sed '1d' > outfile.m
#cat outfile.p outfile.m | sortBed -chrThenSizeA -i - > $AtEST_out

#rm outfile.p outfile.m

echo "Selecting TSS datasets above the defined threshold."

awk '(NR>1) && ($5 >= '10') ' $TokCAGE_out > $(basename $TokCAGE_out .bed)_cutoff.bed
awk '(NR>1) && ($5 >= '4') ' $TokVec_out > $(basename $TokVec_out .bed)_cutoff.bed
awk '(NR>1) && ($5 >= '10') ' $AtEST_out > $(basename $AtEST_out .bed)_cutoff.bed
awk '(NR>1) && ($5 >= '10') ' $PEAT_out > $(basename $PEAT_out .bed)_cutoff.bed

cat $(basename $TokCAGE_out .bed)_cutoff.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

echo "Creating +/- 5bp intervals around TSS"
bedtools slop -l 4 -r 4 -s -i $(basename $TokCAGE_out .bed)_cutoff.bed -g $TAIR10_genome | awk -F'\t' '{ SUM=$3-$2 } {if(SUM>=10) print }' > $(basename $TokCAGE_out _cutoff.bed)_10bp.bed
bedtools slop -l 4 -r 4 -i $(basename $TokVec_out .bed)_cutoff.bed -g $TAIR10_genome | awk -F'\t' '{ SUM=$3-$2 } {if(SUM>=10) print }' > $(basename $TokVec_out _cutoff.bed)_10bp.bed
bedtools slop -l 4 -r 4 -s -i $(basename $AtEST_out .bed)_cutoff.bed -g $TAIR10_genome | awk -F'\t' '{ SUM=$3-$2 } {if(SUM>=10) print }' > $(basename $AtEST_out _cutoff.bed)_10bp.bed
bedtools slop -l 4 -r 4 -s -i $(basename $PEAT_out .bed)_cutoff.bed -g $TAIR10_genome | awk -F'\t' '{ SUM=$3-$2 } {if(SUM>=10) print }' > $(basename $PEAT_out _cutoff.bed)_10bp.bed

echo "Extracting the sequence data for each TSS"
bedtools getfasta -s -fi $TAIR10_fasta -bed $(basename $TokCAGE_out _cutoff.bed)_10bp.bed -fo $(basename $TokCAGE_out _cutoff.bed)_10bp.fa
bedtools getfasta -s -fi $TAIR10_fasta -bed $(basename $TokVec_out _cutoff.bed)_10bp.bed -fo $(basename $TokVec_out _cutoff.bed)_10bp.fa
bedtools getfasta -s -fi $TAIR10_fasta -bed $(basename $AtEST_out _cutoff.bed)_10bp.bed -fo $(basename $AtEST_out _cutoff.bed)_10bp.fa
bedtools getfasta -s -fi $TAIR10_fasta -bed $(basename $PEAT_out _cutoff.bed)_10bp.bed -fo $(basename $PEAT_out _cutoff.bed)_10bp.fa

echo "Job Complete!"
