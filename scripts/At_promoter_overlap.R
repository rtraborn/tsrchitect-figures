
library(rtracklayer)

library(GenomicRanges)

library(UpSetR)

library(ChIPpeakAnno)

setwd("/scratch/rtraborn/tsrchitect-figures/output_files/")

PEAT_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Morton_PEAT/TSRset-1.bed"

TokCAGE_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_CAGE/TSRset-1.bed"

TokVec_m <- "/scratch/rtraborn/TSRchitect_plant_results/Arabidopsis/Tokizawa_Vec_capping/TSRset-1.bed"

At_GROseq <- "/scratch/rtraborn/TSRchitect_plant_results/plant_genomic_data/A_thaliana/GROseq/At_GROseq_out.bed"

PEATtsr <- import.bed(PEAT_m)

peat.vec <- 1:(nrow(as.data.frame(PEATtsr)))

peat.names <- paste("tsr",peat.vec, sep ="_" )

names(PEATtsr) <- peat.names

PEATtsr

TokCAGEtsr <- import.bed(TokCAGE_m)

CAGE.vec <- 1:(nrow(as.data.frame(TokCAGEtsr)))

CAGE.names <- paste("tsr", CAGE.vec, sep="_")

names(TokCAGEtsr) <- CAGE.names

TokCAGEtsr

TokVectsr <- import.bed(TokVec_m)

Oligo.vec <- 1:(nrow(as.data.frame(TokVectsr)))

Oligo.names <- paste("tsr", Oligo.vec, sep="_")

names(TokVectsr) <- Oligo.names

TokVectsr

AtGROpeaks <- import.bed(At_GROseq)

GRO.vec <- 1:(nrow(as.data.frame(AtGROpeaks))) 

GRO.names <- paste("tsr", GRO.vec, sep="_")

names(AtGROpeaks) <- GRO.names

test <- findOverlapsOfPeaks(TokCAGEtsr, TokVectsr, PEATtsr, AtGROpeaks, connectedPeaks = "keepAll", ignore.strand = FALSE)

getVennCounts(TokCAGEtsr, PEATtsr, TokVectsr, AtGROpeaks, ignore.strand = FALSE)

vennOut <- makeVennDiagram(test, NameOfPeaks = c("CAGE","Oligo","PEAT", "GRO"),totalTest = 33026)

vennOut

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
    header = T, sep = ";")

is(movies)

head(movies)

new.df <- as.data.frame(vennOut$vennCounts[-1,1:5])

nrow(new.df)

new.df

expressionInput <- c(CAGE = new.df[8,5], Oligo = new.df[4,5], PEAT = new.df[2,5], GRO = new.df[1,5], `PEAT&GRO` = new.df[3,5], `PEAT&Oligo` = new.df[6,5], 
    `CAGE&Oligo` = new.df[12,5], `CAGE&GRO`= new.df[9,5], `CAGE&PEAT`= new.df[10,5], `Oligo&GRO` = new.df[5,5], `Oligo&GRO&PEAT&CAGE`= new.df[15,5])

new.imput <- fromExpression(expressionInput)

dim(new.imput)

png(filename="TSRintersections.png")

upset(new.imput, nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, main.bar.color = "darkblue", sets.bar.color= "darkred",
    mainbar.y.label = "TSR Intersections", mainbar.y.max=9000, sets.x.label = "Test x", order.by = "freq", scale.intersections="identity", scale.sets="identity")

dev.off()

list.files()
