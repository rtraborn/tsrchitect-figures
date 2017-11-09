## Objective: calculate the coverage of TSSs calculated from ESTs relative to gene annotations: Z. mays

CAGE_dir <- "/projects/TSRplants/ZmCAGE/tsrOut"
CAGEbam_shoot1 <- "/projects/TSRplants/ZmCAGE/alignment/B73-shoot-1_sorted_cut.bam"
CAGEbam_shoot2 <- "/projects/TSRplants/ZmCAGE/alignment/B73-shoot-2_sorted_cut.bam"
CAGEbam_root1 <- "/projects/TSRplants/ZmCAGE/alignment/B73-root-1_sorted_cut.bam"
CAGEbam_root2 <- "/projects/TSRplants/ZmCAGE/alignment/B73-root-2_sorted_cut.bam"

outputDir <- "/scratch/rtraborn/tsrchitect-figures/output_files/"

Zm_gff <- "/projects/TSRplants/ZmCAGE/annotation/Zea_mays.AGPv3.26.gff3"
setwd(outputDir)

library(genomation)
library(rtracklayer)
library(GenomicRanges)

chr.names <- c("1","2","3","4", "5", "6", "7", "8", "9", "10", "Pt", "Mt")

Zm_annot <- import.gff(Zm_gff)
Zm_genes <- Zm_annot[Zm_annot$type=="gene",]
Zm_genes <- Zm_genes[as.character(Zm_genes$biotype)=="protein_coding",]
Zm.names <- as.character(seqnames(Zm_genes))
match.ind <- na.omit(match(Zm.names, chr.names))
new.ind <- which(is.na(match.ind)==FALSE)
Zm_genes <- Zm_genes[new.ind,] #match string of chromosomes
seqlevels(Zm_genes) <- chr.names

Zm_windows <- promoters(Zm_genes, upstream=200, downstream=200)

scores1 <- ScoreMatrix(target=CAGEbam_shoot1, windows=Zm_windows, type='bam', strand.aware=TRUE)
scores2 <- ScoreMatrix(target=CAGEbam_shoot2, windows=Zm_windows, type='bam', strand.aware=TRUE)
scores3 <- ScoreMatrix(target=CAGEbam_root1, windows=Zm_windows, type='bam', strand.aware=TRUE)
scores4 <- ScoreMatrix(target=CAGEbam_root2, windows=Zm_windows, type='bam', strand.aware=TRUE)

sm1.scaled <- scaleScoreMatrix(scores1)
sm2.scaled <- scaleScoreMatrix(scores2)
sm3.scaled <- scaleScoreMatrix(scores3)
sm4.scaled <- scaleScoreMatrix(scores4)

sml=new("ScoreMatrixList",list(shoot1=sm1.scaled, shoot2=sm2.scaled, root1=sm3.scaled, root2=sm4.scaled))
multiHeatMatrix(sml,matrix.main=c("Shoot-1","Shoot-2","Root-1", "Root-2"), cex.axis=0.8, clustfun = function(x) kmeans(x, 
    centers = 2)$cluster)
dev.off()
