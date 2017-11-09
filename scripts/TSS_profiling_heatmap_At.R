## Objective: calculate the coverage of ESTs and TSRs calculated from ESTs relative to gene annotations

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
chr.names <- c("Chr1","Chr2","Chr3","Chr4", "Chr5")
At.names <- as.character(seqnames(At_genes))
match.ind <- na.omit(match(At.names, chr.names))
new.ind <- which(is.na(match.ind)==FALSE)
At_genes <- At_genes[new.ind,]
seqlevels(At_genes) <- chr.names

At_windows <- promoters(At_genes, upstream=200, downstream=200)

scores1 <- ScoreMatrix(target=TokCAGEbam, windows=At_windows, type='bam', strand.aware=TRUE)
scores2 <- ScoreMatrix(target=TokVecbam, windows=At_windows, type='bam', strand.aware=TRUE)
scores3 <- ScoreMatrix(target=PEATbam, windows=At_windows, type='bam', strand.aware=TRUE)

sm1.scaled = scaleScoreMatrix(scores1)
sm2.scaled = scaleScoreMatrix(scores2)
sm3.scaled = scaleScoreMatrix(scores3)

sml=new("ScoreMatrixList",list(cage=sm1.scaled, oligo=sm2.scaled, peat=sm3.scaled))
multiHeatMatrix(sml,matrix.main=c("cage","oligo", "peat"), cex.axis=0.8, clustfun = function(x) kmeans(x, 
    centers = 2)$cluster)
dev.off()
