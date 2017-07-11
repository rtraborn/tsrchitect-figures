## Objective: calculate the PWM for consensus promoters vs all Tokizawa promoters

At_fa=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10_genome.fasta
At_g=/projects/rtraborn/plant_promoters/A_thaliana/genome/TAIR10.genome

PEAT_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT
PEATtsrs=TSRset-1.bed
PEATtss=TSSset-1.bed
PEAT_m=TSSset-1_sum.bed
PEAT_out=PEAT_out.bed
PEAT_out2=PEAT_out2.bed
PEAT_fa=PEAT_tsrs.fa

TokCAGE_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE
TokCAGE_tsrs=TSRset-1_nuc.bed
TokCAGE_tss=TSSset-1.bed
TokCAGE_m=TSSset-1_sum.bed
TokCAGE_out=TokCAGE_out.bed
TokCAGE_out2=TokCAGE_out2.bed
TokCAGE_fa=TokizawaCAGE_tsrs.fa

TokVec_dir=/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping
TokVec_tsrs=TSRset-1_nuc.bed
TokVec_tss=TSSset-1.bed
TokVec_m=TSSset-1_sum.bed
TokVec_out=TokVec_out.bed
TokVec_out2=TokVec_out2.bed
TokVec_fa=TokizawaVec_tsrs.fa

uniq -c $PEAT_dir/$PEATtss | awk '{ print $2,$3,$4,$5,$1,$7}' | perl -p -e 's/ /\t/g' | awk 'BEGIN{OFS="\t";} $5>=3 {print }' > $PEAT_dir/$PEAT_m
bedtools intersect -s -u -a $PEAT_dir/$PEAT_m -b $PEAT_dir/$PEATtsrs > $PEAT_dir/$PEAT_out
bedtools slop -l 2 -r 2 -s -i $PEAT_dir/$PEAT_out -g $At_g  > $PEAT_dir/$PEAT_out2
bedtools getfasta -s -fi $At_fa -bed $PEAT_dir/$PEAT_out2 -fo $PEAT_fa

echo "Finished job 1"

uniq -c $TokCAGE_dir/$TokCAGE_tss | awk '{ print $2,$3,$4,$5,$1,$7}' | perl -p -e 's/ /\t/g' | awk 'BEGIN{OFS="\t";} $5>=3 {print }' > $TokCAGE_dir/$TokCAGE_m
bedtools intersect -s -u -a $TokCAGE_dir/$TokCAGE_m -b $TokCAGE_dir/$TokCAGE_tsrs > $TokCAGE_dir/$TokCAGE_out
bedtools slop -l 2 -r 2 -s -i $TokCAGE_dir/$TokCAGE_out -g $At_g > $TokCAGE_dir/$TokCAGE_out2
bedtools getfasta -s -fi $At_fa -bed $TokCAGE_dir/$TokCAGE_out2 -fo $TokCAGE_fa

echo "Finished job 2"

uniq -c $TokVec_dir/$TokVec_tss | awk '{ print $2,$3,$4,$5,$1,$7}' | perl -p -e 's/ /\t/g' | awk 'BEGIN{OFS="\t";} $5>=3 {print }' > $TokVec_dir/$TokVec_m
bedtools intersect -s -u -a $TokVec_dir/$TokVec_m -b $TokVec_dir/$TokVec_tsrs > $TokVec_dir/$TokVec_out
bedtools slop -l 2 -r 2 -s -i $TokVec_dir/$TokVec_out -g $At_g  > $TokVec_dir/$TokVec_out2
bedtools getfasta -s -fi $At_fa -bed $TokVec_dir/$TokVec_out2 -fo $TokVec_fa

echo "Finished job 3"

echo "Task complete!"
