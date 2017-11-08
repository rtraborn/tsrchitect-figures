## Objective: calculate the coverage of ESTs and TSRs calculated from ESTs relative to gene annotations

CAGE_dir="/projects/TSRplants/ZmCAGE/tsrOut"
CAGE_m="/projects/TSRplants/ZmCAGE/tsrOut/TSRsetCombined.bed"

PEAT_dir="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT"
PEAT_m="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSSset-1_sum.bed"
PEATbam="/projects/TSRplants/AtPEAT/alignment/AtPEAT_reheader.bam"

TokCAGE_dir="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE"
TokCAGE_m="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSSset-1_sum.bed"
TokCAGEbam="/projects/TSRplants/AtCAGE/alignment/AtCAGE.aligned.bam"

TokVecbam="/projects/TSRplants/AtOligo/alignment/OligoCap_align_sorted.bam"
TokVec_dir="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping"
TokVec_m="/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSSset-1_sum.bed"

AtGenbankEST="/projects/TSRplants/ESTcDNA/A_thaliana/tsrOut/AtEST_genbank_out/TSRset-1.bed"
At_GROseqBED="/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed"
At_GROseq="/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/GROseq_tagDirectory/tss.txt"

AtTAIR_genes="/projects/rtraborn/plant_promoters/A_thaliana/genome/AtTAIR_genes_reheader.gff"

PEAT_out="PEAT_TSSset-1.bed"
TokCAGE_out="TokCAGE_TSSset-1.bed"
TokVec_out="TokVec_TSSset-1.bed"
AtEST_out="AtEST_TSSset-1.bed"

library(genomation)
library(rtracklayer)
library(GenomicRanges)

AtESTbed <- import.bed(AtEST_out)
At_annot <- import.gff(AtTAIR_genes)
At_genes <- At_annot[At_annot$type=="gene",]
At_genes <- At_genes[as.character(At_genes$Note)=="protein_coding_gene",]

#cat $AtTAIR_genes | awk '{gsub(/ChrM/,"mitochondria")}; 1' > AtTAIR_genes_i
#cat AtTAIR_genes_i | awk '{gsub(/ChrC/,"chloroplast")}; 1' > AtTAIR_genes_reheader


extend <- function(x, upstream=100, downstream=100)
{
    if (any(strand(x) == "*"))
       warning("'*' ranges were treated as '+'")
       on_plus <- strand(x) == "+" | strand(x) == "*"
       new_start <- start(x) - ifelse(on_plus, upstream, downstream)
       new_end <- end(x) + ifelse(on_plus, downstream, upstream)
       ranges(x) <- IRanges(new_start, new_end)
       trim(x)
}

At_windows <- promoters(At_genes, upstream=200, downstream=200)

scores1 <- ScoreMatrix(target=TokCAGEbam, windows=At_windows, type='bam', strand.aware=TRUE)
scores2 <- ScoreMatrix(target=TokVecbam, windows=At_windows, type='bam', strand.aware=TRUE)
#scores3 <- ScoreMatrix(target=PEATbam, windows=At_windows, type='bam', strand.aware=TRUE)

sm1.scaled = scaleScoreMatrix(scores1)
sm2.scaled = scaleScoreMatrix(scores2)

#sml=new("ScoreMatrixList",list(CAGE=scores1, OligoCap=scores2, PEAT=scores3))
#multiHeatMatrix(sml,kmeans=TRUE,k=2, matrix.main=c("cage","OligoCap","peat"), cex.axis=0.8)

#sm = ScoreMatrix(target = cage, windows = promoters, strand.aware = TRUE)
#cpg.ind = which(countOverlaps(promoters, cpgi) > 0)
#nocpg.ind = which(countOverlaps(promoters, cpgi) == 0)
#heatMatrix(sm, xcoords = c(-1000, 1000), group = list(CpGi = cpg.ind, noCpGi = nocpg.ind))

sml=new("ScoreMatrixList",list(cage=sm1.scaled, oligo=sm2.scaled))
#sml=new("ScoreMatrixList",list(cage=scores1, oligo=scores2))
multiHeatMatrix(sml,matrix.main=c("cage","oligo"), cex.axis=0.8, clustfun = function(x) kmeans(x, 
    centers = 2)$cluster)
dev.off()
