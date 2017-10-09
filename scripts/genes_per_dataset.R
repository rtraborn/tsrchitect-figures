## Measuring the number of promoters identified using TSRchtect, and the total gene coverage across the 4 plant experiments (2 species)

library(ggplot2)
library(GenomicRanges)

setwd("/scratch/rtraborn/tsrchitect-figures/figures")

AtEST.file <- "/scratch/rtraborn/tsrchitect-figures/output_files/AtEST_closest_gene_i.txt"
AtESTdist <- read.table(file=AtEST.file, header=FALSE)

head(AtESTdist)

colnames(AtESTdist) <- c("chr","start1","end1", "ID1", "score", "strand1", "chr2", "start2", "end2", "ID2", "score2", "strand2", "assemblyID", "type", "score3", "ID_string", "distance")

head(AtESTdist)

a <- ggplot(AtESTdist, aes(distance))
a + geom_histogram()

ggsave("AtEST_closest.png")

n.prom <- matrix(NA, nrow=4, ncol=3)

under100 <- length(which(AtESTdist[,17]<100))
under500 <- length(which(AtESTdist[,17]<500))
under1000 <- length(which(AtESTdist[,17]<1000))
above1k  <- length(which(AtESTdist[,17]>1000))

n.prom[,1] <- c("AtEST")                   
n.prom[,2] <- c(under100, under500, under1000, above1k)
n.prom[,3] <- c("under100", "under500", "under1000", "above1k")

n.prom <- as.data.frame(n.prom)
colnames(n.prom) <- c("dataset", "class", "distance")
                   
head(n.prom)

b <- ggplot(n.prom, aes(distance))
b + geom_bar()

ggsave("distance_classes.png")
