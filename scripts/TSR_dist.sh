## Objective: calculate the distribution of TSSs from various datsets to annotation

At_fa=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genome.fasta
At_g=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10.genome

PEAT_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT
PEAT_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1_sum.bed

TokCAGE_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE
TokCAGE_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1_sum.bed

TokVec_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping
TokVec_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1_sum.bed

At_GROseq=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed
At_TAIR_gff=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_GFF3_genes.gff
At_TAIR_gene=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_GFF3_genes.bed

#gff2bed < $At_TAIR_Gff > $At_TAIR_gene

bedtools closest -s -D a  -a $PEAT_m -b $At_GROseq > PEAT_closest_GRO.txt
bedtools closest -s -D a  -a $TokCAGE_m -b $At_GROseq > TokCAGE_closest_GRO.txt
bedtools closest -s -D a  -a $TokVec_m -b $At_GROseq > TokVec_closest_GRO.txt

bedtools closest -s -D a  -a $PEAT_m -b $At_TAIR_gene > PEAT_closest_gene.txt
bedtools closest -s -D a  -a $TokCAGE_m -b $At_TAIR_gene > TokCAGE_closest_gene.txt
bedtools closest -s -D a  -a $TokVec_m -b $At_TAIR_gene > TokVec_closest_gene.txt

echo "Job Complete!"
