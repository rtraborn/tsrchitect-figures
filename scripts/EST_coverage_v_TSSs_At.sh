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


PEAT_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT
PEAT_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1_sum.bed
PEATbam=/projects/TSRplants/AtPEAT/alignment/AtPEAT.bam

TokCAGE_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE
TokCAGE_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1_sum.bed
TokCAGEbam=/projects/TSRplants/AtCAGE/alignment/AtCAGE.aligned.bam

TokVecbam=/projects/TSRplants/AtOligo/alignment/OligoCap_align.bam
TokVec_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping
TokVec_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1_sum.bed

AtGenbankEST=/projects/TSRplants/ESTcDNA/A_thaliana/tsrOut/AtEST_genbank_out/TSRset-1.bed

At_GROseqBED=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed
At_GROseq=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/GROseq_tagDirectory/tss.txt

PEAT_out=PEAT_TSSset-1.bed
TokCAGE_out=TokCAGE_TSSset-1.bed
TokVec_out=TokVec_TSSset-1.bed
AtEST_out=AtEST_TSSset-1.bed

At_TAIR_gff=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_GFF3_genes.gff
At_TAIR_genes=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genes.bed
At_TAIR_genes_nuc=/projects/TSRplants/ESTcDNA/A_thaliana/TAIR10_genes_nuc.bed
At_TAIR_genome=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10.genome

AtEST_coverage=AtEST.coverage
ZmEST_coverage=ZmEST.coverage

#gff2bed < $At_TAIR_gff > $At_TAIR_genes
#cat $At_TAIR_genes | awk 'BEGIN{OFS="\t";} $1!="choloroplast" {print }' > $At_TAIR_genes_nuc

#gff2bed < $Zm_gff > $Zm_bed
#cat $Zm_bed | awk 'BEGIN{OFS="\t";} $8=="gene" {print }' > $Zm_genes
#cat $Zm_genes | awk 'BEGIN{OFS="\t";} $1!="Mt" {print }' > genes_int
#cat genes_int | awk 'BEGIN{OFS="\t";} $1!="Pt" {print }' > $Zm_genes_nuc

#cat $Zm_bed | awk 'BEGIN{OFS="\t";} $8!="repeat_region" {print }' > Zm_test.bed

echo "Creating 5' end datasets from the TAIR10 annotation"
#bedtools flank -l 1 -r 0 -s -i $At_TAIR_genes -g $At_TAIR_genome > TAIRgenes_5prime.bed

echo "Computing coverage of TSS profiling datasets vs TSSs"

#makeTagDirectory ESTtagDir $AtEST_out -format bed
#makeTagDirectory CAGEtagDir $TokCAGEbam -format sam
#makeTagDirectory TokVecDir $TokVecbam -format sam
#makeTagDirectory PEATtagDir $PEATbam -format sam

echo "Computing coverage histograms of ESTs and TSS profiling data relative to tair10 annotated TSSs"

annotatePeaks.pl TAIRgenes_5prime.bed tair10_new -size -500,100 -hist 10 -d ESTtagDir CAGEtagDir PEATtagDir TokVecDir > $outputDir/At_coverage_hist_all_v_tss_annot.txt

echo "Job Complete!"
