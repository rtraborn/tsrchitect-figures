## Objective: calculate the coverage of ESTs and TSRs calculated from ESTs relative to gene annotations

CAGE_dir=/projects/TSRplants/ZmCAGE/tsrOut
CAGE_m=/projects/TSRplants/ZmCAGE/tsrOut/TSRsetCombined.bed
CAGEtss=/projects/TSRplants/ZmCAGE/tsrOut/TSSset-1.txt
CAGEtss_out=TSSset-1.bed
CAGEbam_shoot1=/projects/TSRplants/ZmCAGE/alignment/B73-shoot-1.bam
CAGEbam_shoot2=/projects/TSRplants/ZmCAGE/alignment/B73-shoot-2.bam

ZmGenbankEST=/projects/TSRplants/ESTcDNA/Z_mays/tsrOut/TSRset-1_reheader.bed

outputDir=/scratch/rtraborn/tsrchitect-figures/output_files
Zm_gff=/projects/TSRplants/ZmCAGE/annotation/Zea_mays.AGPv3.26.gff3
Zm_bed=/projects/TSRplants/ZmCAGE/annotation/Zea_mays.AGPv3.26.bed
Zm_genes=/projects/TSRplants/ZmCAGE/annotation/Zm_genes.bed
Zm_genes_nuc=/projects/TSRplants/ZmCAGE/annotation/Zm_genes_nuc.bed
Zm_output_dir=/projects/TSRplants/ESTcDNA/Z_mays/tsrOut/

#cat $CAGEtss | awk -F'\t' -v OFS="\t" '{ if($3=="+") print $1, $2, $2+1, ".", $4, $3 }' | sed '1d' > outfile.p
#cat $CAGEtss | awk -F'\t' -v OFS="\t" '{ if($3=="-") print $1, $2-1, $2, ".", $4, $3 }' | sed '1d' > outfile.m
#cat outfile.p outfile.m | sortBed -chrThenSizeA -i - > $CAGEtss_out


#gff2bed < $Zm_gff > $Zm_bed
#cat $Zm_bed | awk 'BEGIN{OFS="\t";} $8=="gene" {print }' > $Zm_genes
#cat $Zm_genes | awk 'BEGIN{OFS="\t";} $1!="Mt" {print }' > genes_int
#cat genes_int | awk 'BEGIN{OFS="\t";} $1!="Pt" {print }' > $Zm_genes_nuc

#cat $Zm_bed | awk 'BEGIN{OFS="\t";} $8!="repeat_region" {print }' > Zm_genes_new.bed

echo "Computing coverage of TSS profiling datasets vs TSSs"

makeTagDirectory ZmESTtagDir $ZmGenbankEST -format bed
findPeaks ZmESTtagDir -style tss -o auto
#makeTagDirectory ZmCAGEtagDir_shoot1 $CAGEbam_shoot1 -format sam
#findPeaks ZmCAGEtagDir_shoot1 -style tss -o auto
#makeTagDirectory ZmCAGEtagDir_shoot2 $CAGEbam_shoot2 -format sam
#findPeaks ZmCAGEtagDir_shoot2 -style tss -o auto

echo "Computing coverage histograms of ESTs and TSS profiling data relative to tair10 annotated TSSs"

annotatePeaks.pl Zm_genes_new.bed corn.AGPv3 -size -500,100 -hist 10 -d ZmESTtagDir ZmCAGEtagDir_shoot1 ZmCAGEtagDir_shoot2 > $outputDir/Zm_coverage_hist_all_v_tss_annot.txt

echo "Job Complete!"
