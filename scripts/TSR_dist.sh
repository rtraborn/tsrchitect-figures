## Objective: calculate the distribution of TSSs from various datsets to annotation

At_fa=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genome.fasta
At_g=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10.genome

PEAT_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT
PEAT_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1_sum.bed

TokCAGE_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE
TokCAGE_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1_sum.bed

TokVec_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping
TokVec_m=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1_sum.bed

GenbankEST=/projects/TSRplants/ESTcDNA/A_thaliana/tsrOut/AtEST_genbank_out/TSRset-1.bed

At_GROseq=/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed
At_TAIR_gff=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_GFF3_genes.gff
At_TAIR_genes=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genes.bed
At_TAIR_genes_nuc=TAIR10_genes_nuc.bed

#gff2bed < $At_TAIR_gff > $At_TAIR_genes

#head $At_TAIR_genes

cat $At_TAIR_genes | awk 'BEGIN{OFS="\t";} $1!="choloroplast" {print }' > $At_TAIR_genes_nuc

head $At_TAIR_genes_nuc

#bedtools closest -s -Da  -a $PEAT_m -b $At_GROseq > PEAT_closest_GRO.txt
#bedtools closest -s -Da  -a $TokCAGE_m -b $At_GROseq > TokCAGE_closest_GRO.txt
#bedtools closest -s -Da  -a $TokVec_m -b $At_GROseq > TokVec_closest_GRO.txt

#bedtools closest -s -Da  -a $PEAT_m -b $At_TAIR_gene > PEAT_closest_gene.txt
#bedtools closest -s -Da  -a $TokCAGE_m -b $At_TAIR_gene > TokCAGE_closest_gene.txt
#bedtools closest -s -Da  -a $TokVec_m -b $At_TAIR_gene > TokVec_closest_gene.txt

bedtools closest -s -iu -D a -a $GenbankEST -b $At_TAIR_genes_nuc > AtEST_closest_gene.txt
cat AtEST_closest_gene.txt | awk 'BEGIN{OFS="\t";} $17!="-1" {print }' > AtEST_closest_gene_i.txt
mv AtEST_closest_gene_i.txt /scratch/rtraborn/tsrchitect-figures/output_files

echo "Job Complete!"
