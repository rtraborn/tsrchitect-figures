myWD <- "/scratch/rtraborn/tsrchitect-figures/"
setwd(myWD)

PEATtsr <- "/projects/TSRplantsV3.0/At/AtPEAT/tsr/TSRset-1.txt"
CAGEtsr <- "/projects/TSRplantsV3.0/At/AtCAGE/tsr/TSRset-1.txt"
Oligotsr <- "/projects/TSRplantsV3.0/At/AtOligo/tsr/TSRset-1.txt"

AtTSR <- read.table(PEATtsr, header=TRUE)                                     
AtTSR.p <- AtTSR[AtTSR$strand=="+",]
AtTSR.m <- AtTSR[AtTSR$strand=="-",]
tsr.p.start <- AtTSR.p$nTSSs - 1
tsr.m.end <- AtTSR.m$nTSSs + 1
tsr.p.start <- AtTSR.p$TSS - 1
tsr.m.end <- AtTSR.m$TSS + 1

n.peaks <- 1:nrow(AtTSR)
peak.ids <- paste("tsr", n.peaks, sep="_")
AtTSR_df <- data.frame(cbind(as.character(AtTSR$seq), AtTSR$start, AtTSR$end, peak.ids, AtTSR$nTAGs, as.character(AtTSR$strand)))
colnames(AtTSR_df) <- c("chr","start", "end", "tsrID", "nTAGs", "strand")

write.table(AtTSR_df, file="AtPEAT_TSRset-1.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

AtTSR <- read.table(CAGEtsr, header=TRUE)                                     

AtTSR.p <- AtTSR[AtTSR$strand=="+",]
AtTSR.m <- AtTSR[AtTSR$strand=="-",]
tsr.p.start <- AtTSR.p$nTSSs - 1
tsr.m.end <- AtTSR.m$nTSSs + 1
tsr.p.start <- AtTSR.p$TSS - 1
tsr.m.end <- AtTSR.m$TSS + 1

n.peaks <- 1:nrow(AtTSR)
peak.ids <- paste("tsr", n.peaks, sep="_")
AtTSR_df <- data.frame(cbind(as.character(AtTSR$seq), AtTSR$start, AtTSR$end, peak.ids, AtTSR$nTAGs, as.character(AtTSR$strand)))
colnames(AtTSR_df) <- c("chr","start", "end", "tsrID", "nTAGs", "strand")

write.table(AtTSR_df, file="AtCAGE_TSRset-1.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

AtTSR <- read.table(Oligotsr, header=TRUE)                                     

AtTSR.p <- AtTSR[AtTSR$strand=="+",]
AtTSR.m <- AtTSR[AtTSR$strand=="-",]
tsr.p.start <- AtTSR.p$nTSSs - 1
tsr.m.end <- AtTSR.m$nTSSs + 1
tsr.p.start <- AtTSR.p$TSS - 1
tsr.m.end <- AtTSR.m$TSS + 1

n.peaks <- 1:nrow(AtTSR)
peak.ids <- paste("tsr", n.peaks, sep="_")
AtTSR_df <- data.frame(cbind(as.character(AtTSR$seq), AtTSR$start, AtTSR$end, peak.ids, AtTSR$nTAGs, as.character(AtTSR$strand)))
colnames(AtTSR_df) <- c("chr","start", "end", "tsrID", "nTAGs", "strand")

write.table(AtTSR_df, file="AtOligo_TSRset-1.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
