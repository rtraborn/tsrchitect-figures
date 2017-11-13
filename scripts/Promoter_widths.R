
library(rtracklayer)

library(GenomicRanges)

library("ggplot2")

library("reshape2")

setwd("/scratch/rtraborn/tsrchitect-figures/output_files/")

PEAT_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSRset-1.bed"

TokCAGE_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSRset-1.bed"

TokVec_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSRset-1.bed"

CAGEzm <- "/projects/TSRplants/ZmCAGE/tsrOut/TSRsetCombined.bed"

PEATtsr <- import.bed(PEAT_m)

TokCAGEtsr <- import.bed(TokCAGE_m)

TokVectsr <- import.bed(TokVec_m)

zmCAGEtsr <- import.bed(CAGEzm)

dim(PEATtsr)

AtPEAT_width <- cbind(width(PEATtsr))

colnames(AtPEAT_width) <- c("A_thaliana_PEAT")

TokVec_width <- cbind(width(TokVectsr))

colnames(TokVec_width) <- c("A_thaliana_Oligo")

TokCAGE_width <- cbind(width(TokCAGEtsr))

colnames(TokCAGE_width) <- c("A_thaliana_CAGE")

zmCAGE_width <- cbind(width(zmCAGEtsr))

colnames(zmCAGE_width) <- c("Zea_mays_CAGE")

AtPEAT_m <- melt(AtPEAT_width)

TokCAGE_melt <- melt(TokCAGE_width)

TokVec_melt <- melt(TokVec_width)

zmCAGE_m <- melt(zmCAGE_width)

dim(TokVec_melt)

dim(zmCAGE_m)

width_df <- rbind(AtPEAT_m, TokVec_melt, TokCAGE_melt, zmCAGE_m)

width_df <- width_df[,-1]

colnames(width_df) <- c("sampleID", "TSR_width")

g <- ggplot(width_df, aes(sampleID, TSR_width))

g + geom_boxplot(fill="lightblue") + scale_y_continuous(limits=c(0,100))

ggsave(filename="TSR_widths_boxplot.png")

t.test(AtPEAT_width,zmCAGE_width, alternative='greater')

t.test(TokVec_width,zmCAGE_width, alternative='greater')

t.test(TokCAGE_width,zmCAGE_width, alternative='greater')
