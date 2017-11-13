
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

PEATtsr_df <- as.data.frame(PEATtsr)

TokVectsr_df <- as.data.frame(TokVectsr)

TokCAGEtsr_df <- as.data.frame(TokCAGEtsr)

zmCAGE_df <- as.data.frame(zmCAGEtsr)

PEATtsr_df$dataset <- "AtPEAT"

TokVectsr_df$dataset <- "AtOligo"

TokCAGEtsr_df$dataset <- "AtCAGE"

zmCAGE_df$dataset <- "zmCAGE"

score.df <- rbind(PEATtsr_df[,7:8],TokVectsr_df[,7:8], TokCAGEtsr_df[,7:8], zmCAGE_df[,7:8])

colnames(score.df) <- c("SI","Dataset")

head(score.df)

g <- ggplot(score.df, aes(Dataset),stat="count")

g + geom_bar(position="dodge", fill="orange3") + coord_flip()

ggsave(filename="TSR_count_barplot.png")

g <- ggplot(score.df, aes(Dataset, SI))

g + geom_violin(fill="purple2") + scale_y_reverse() + coord_flip()

ggsave(filename="TSR_shape_violinplot.png")
