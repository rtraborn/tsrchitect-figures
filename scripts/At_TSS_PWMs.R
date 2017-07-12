
library(Biostrings)

library(seqLogo)

setwd("/scratch/rtraborn/TSRchitect_plant_results/analysis_scripts")

fa.list <- list.files(pattern="*.fa")

fa.list

PEAT_fa <- readDNAStringSet(fa.list[1], format="fasta")

PEAT_fa

PEAT_fa_i <- PEAT_fa[width(PEAT_fa)==5,]

PEAT_fa_i

PEAT_matrix <- consensusMatrix(PEAT_fa_i)

PEAT_matrix2 <- PEAT_matrix[1:4,]

dim(PEAT_matrix2)

head(PEAT_matrix2)

PEAT_sum <- apply(PEAT_matrix2, 2, "sum")

PEAT_ma <- rbind(PEAT_matrix2, PEAT_sum)

#colnames(PEAT_df) <- c("1","2","3","4","5")

is(PEAT_ma)

PEAT_sw <- sweep(PEAT_ma, 2, PEAT_ma[5,],'/')

PEAT_sw2 <- PEAT_sw[-5,]

PEATpwm <- makePWM(PEAT_sw2)

PEATpwm

is(PEATpwm)

png(file="PEATpwm.png")

seqLogo(PEATpwm, ic.scale=TRUE)

list.files()



TokCAGE_fa <- readDNAStringSet(fa.list[2], format="fasta")

TokCAGE_fa_i <- TokCAGE_fa[width(TokCAGE_fa)==5,]

TokCAGE_fa_i

TokCAGE_matrix <- consensusMatrix(TokCAGE_fa_i)

TokCAGE_matrix2 <- TokCAGE_matrix[1:4,]

dim(TokCAGE_matrix2)

head(TokCAGE_matrix2)

TokCAGE_sum <- apply(TokCAGE_matrix2, 2, "sum")

TokCAGE_ma <- rbind(TokCAGE_matrix2, TokCAGE_sum)

TokCAGE_sw <- sweep(TokCAGE_ma, 2, TokCAGE_ma[5,],'/')

TokCAGE_sw2 <- TokCAGE_sw[-5,]

TokCAGEpwm <- makePWM(TokCAGE_sw2)

TokCAGEpwm

is(TokCAGEpwm)

png(file="TokCAGEpwm.png")

seqLogo(TokCAGEpwm, ic.scale=TRUE)

TokVec_fa <- readDNAStringSet(fa.list[3], format="fasta")

TokVec_fa

TokVec_fa_i <- TokVec_fa[width(TokVec_fa)==5,]

TokVec_fa_i

TokVec_matrix <- consensusMatrix(TokVec_fa_i)

TokVec_matrix2 <- TokVec_matrix[1:4,]

dim(TokVec_matrix2)

head(TokVec_matrix2)

TokVec_sum <- apply(TokVec_matrix2, 2, "sum")

TokVec_ma <- rbind(TokVec_matrix2, TokVec_sum)

TokVec_sw <- sweep(TokVec_ma, 2, TokVec_ma[5,],'/')

TokVec_sw2 <- TokVec_sw[-5,]

TokVec_sw2

TokVecpwm <- makePWM(TokVec_sw2)

TokVecpwm

png(file="TokVec_pwm.png")

seqLogo(TokVecpwm, ic.scale=TRUE)
