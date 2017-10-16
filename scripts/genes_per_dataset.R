
library(ggplot2)

library(GenomicRanges)

library(reshape2)

setwd("/scratch/rtraborn/tsrchitect-figures/figures")

AtEST.file <- "/scratch/rtraborn/tsrchitect-figures/output_files/AtEST_closest_gene_i.txt"

AtESTdist <- read.table(file=AtEST.file, header=FALSE)

ZmEST.file <- "/scratch/rtraborn/tsrchitect-figures/output_files/ZmEST_closest_gene.txt"

ZmESTdist <- read.table(file=ZmEST.file, sep="\t",header=FALSE)

head(AtESTdist)

head(ZmESTdist)

colnames(AtESTdist) <- c("chr","start1","end1", "ID1", "score", "strand1", "chr2", "start2", "end2", "ID2", "score2", "strand2", "assemblyI\
D", "type", "score3", "ID_string", "distance")

colnames(ZmESTdist) <- c("chr","start1","end1", "ID1", "score", "strand1", "chr2", "start2", "end2", "ID2", "score2", "strand2", "assemblyI\
D", "type", "score3", "ID_string", "distance")

head(AtESTdist)

dim(AtESTdist)

is(AtESTdist$distance)

a <- ggplot(AtESTdist, aes(distance))

a + geom_histogram(binwidth=500, colour="white",fill="navy") + scale_x_continuous(limit=c(0,50000)) + scale_y_continuous(limit=c(0,500))

ggsave("AtEST_closest_hist.png")

dim(ZmESTdist)

z <- ggplot(ZmESTdist, aes(distance))

z + geom_histogram(binwidth=500, colour="black",fill="gold") + scale_x_continuous(limit=c(0,50000)) + scale_y_continuous(limit=c(0,200))

ggsave("ZmEST_closest_hist.png")

n.prom <- matrix(NA, nrow=4, ncol=2)

under200 <- length(which(AtESTdist[,17]<200))

under500 <- length(which(AtESTdist[,17]<500))

under1000 <- length(which(AtESTdist[,17]<1000))

above1k  <- length(which(AtESTdist[,17]>=1000))

under200.ind <- which(AtESTdist[,17]<200)

under500.ind <- which(AtESTdist[,17]<500)

under1k.ind <- which(AtESTdist[,17]<1000)

above1k.ind <- which(AtESTdist[,17]>=1000)

dim(AtESTdist)

prom.cat <- matrix(NA, nrow=nrow(AtESTdist), ncol=2)

prom.cat[above1k.ind,1] <- "above1k"

prom.cat[under1k.ind,1] <- "under1k"

prom.cat[under500.ind,1] <- "under500"

prom.cat[under200.ind,1] <- "under200"

colnames(prom.cat) <- c("group","distance")

prom.cat <- as.data.frame(prom.cat)

above1k.c <- length(which(as.character(prom.cat$group)=="above1k"))

under1k.c <- length(which(as.character(prom.cat$group)=="under1k"))

under500.c <- length(which(as.character(prom.cat$group)=="under500"))

under200.c <- length(which(as.character(prom.cat$group)=="under200"))

above1k.c

under1k.c

under500.c

under200.c

prom.cat[,2] <- AtESTdist$distance

head(prom.cat)

prom.cat$distance <- as.numeric(unlist(prom.cat$distance))

prom.cat$group <- factor(prom.cat$group, levels = prom.cat$group[order(prom.cat$distance)])

prom.cat$org <- "Arabidopsis"

is(prom.cat$distance)

n.prom.z <- matrix(NA, nrow=4, ncol=2)

under200.z <- length(which(ZmESTdist[,17]<200))

under500.z <- length(which(ZmESTdist[,17]<500))

under1000.z <- length(which(ZmESTdist[,17]<1000))

above1k.z  <- length(which(ZmESTdist[,17]>1000))

under200.ind.z <- which(ZmESTdist[,17]<200)

under500.ind.z <- which(ZmESTdist[,17]<500)

under1k.ind.z <- which(ZmESTdist[,17]<1000)

above1k.ind.z <- which(ZmESTdist[,17]>=1000)

prom.cat.z <- matrix(NA, nrow=nrow(ZmESTdist), ncol=2)

prom.cat.z[above1k.ind.z,1] <- "above1k"

prom.cat.z[under1k.ind.z,1] <- "under1k"

prom.cat.z[under500.ind.z,1] <- "under500"

prom.cat.z[under100.ind.z,1] <- "under200"

prom.cat.z[,2] <- ZmESTdist$distance

colnames(prom.cat.z) <- c("group","distance")

dim(prom.cat.z)

prom.cat.z <- as.data.frame(prom.cat.z)

prom.cat.z$org <- "Maize"

head(prom.cat.z)

above1k.z <- length(which(as.character(prom.cat.z$group)=="above1k"))

under1k.z <- length(which(as.character(prom.cat.z$group)=="under1k"))

under500.z <- length(which(as.character(prom.cat.z$group)=="under500"))

under200.z <- length(which(as.character(prom.cat.z$group)=="under200"))

length(na.omit(prom.cat.z$group))

prom.count.at <- matrix(NA, ncol=3, nrow=4)

colnames(prom.count.at) <- c("group","number","org")

prom.count.at[,1] <- c("under200", "under500", "under1k", "above1k")

prom.count.at[,2] <- c(under200.c, under500.c, under1k.c, above1k.c)

prom.count.at <- as.data.frame(prom.count.at)

prom.count.at$group <- factor(prom.count.at$group, levels = prom.count.at$group[order(prom.count.at$number)])

prom.count.at$org <- "Arabidopsis"

prom.count.at$number <- as.numeric(unlist(as.character(prom.count.at$number)))

head(prom.count.at)

prom.count.z <- matrix(NA, ncol=3, nrow=4)

colnames(prom.count.z) <- c("group","number","org")

prom.count.z[,1] <- c("under200", "under500", "under1k", "above1k")

prom.count.z[,2] <- c(under200.z, under500.z, under1k.z, above1k.z)

prom.count.z <- as.data.frame(prom.count.z)

prom.count.z$group <- factor(prom.count.z$group, levels = prom.count.z$group[order(prom.count.z$number)])

prom.count.z$org <- "Maize"

prom.count.all <- rbind(prom.count.at, prom.count.z)

prom.count.z

prom.count.z$group

prom.count.all

prom.count.all$number <- as.numeric(unlist(as.character(prom.count.all$number)))

is(prom.count.all$number)

prom.count.all$org

prom.count.z$number <- as.numeric(unlist(as.character(prom.count.z$number)))

c <- ggplot(prom.count.z, aes(y=number, x=org, fill=group))

cc <- c + geom_bar(stat="identity")                                                               

cc

ggsave(file="ZmEST_distance_classes.png")

d <- ggplot(prom.count.all, aes(y=number, x=org, fill=group))

dd <- d + geom_bar(stat="identity")                                                               

dd

ggsave(file="EST_distance_classes_combined.png")
