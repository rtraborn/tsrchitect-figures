## Measuring the total gene coverage for GROseq and various TSS profiling datasets in Arabidopsis
library(genomation)
library(GenomicRanges)

TSS.bams <- list.files(path="/scratch/rtraborn/TSRchitect_plant_results/alignment_data", pattern="\\.bam$", include.dirs=TRUE, full.names=TRUE)
PEAT.bam <- TSS.bams[2]
At_CAGE <- TSS.bams[3]
At_vec_capping <- TSS.bams[4]

GRO.bams <- list.files(path="/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/alignment/", pattern="\\.bam", include.dirs=TRUE)

TAIR10_annot <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/annotation_data/TAIR10_GFF3_genes.gff"
annot.import <- gffToGRanges(TAIR10_annot, filter="CDS")

#filtering the non-nuclear chromosomes
stuff <- annot.import[seqnames(annot.import) != "ChrM"]
stuff2 <- stuff[seqnames(stuff) != "ChrC"]
annot.final <- dropSeqlevels(stuff2, c("ChrM", "ChrC"))

sm <- ScoreMatrixBin(target=PEAT.bam, windows=annot.final, bin.num=10)
heatMatrix(sm, xcoords = c(-500, 500))
plotMeta(sm, xcoords = c(-1000, 200))




