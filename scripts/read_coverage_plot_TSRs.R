#Parses the output of bedtools coverage and creates stacked barplots reflecting the genomic distribution of TSSs

library(ggplot2)
library(reshape2)
setwd("/projects/TSRplantsV2.0/At/AtCoveragePlots")

#AtPEAT
AtPEAT.upstream <- "AtPEAT_upstream.tsr"
AtPEAT.downstream <- "AtPEAT_downstream.tsr"
AtPEAT.cds <- "AtPEAT_cds.tsr"
AtPEAT.introns <- "AtPEAT_introns.tsr"
AtPEAT.intergenic <- "AtPEAT_intergenic.tsr"

#AtCAGE
AtCAGE.upstream <- "AtCAGE_upstream.tsr"
AtCAGE.downstream <- "AtCAGE_downstream.tsr"
AtCAGE.cds <- "AtCAGE_cds.tsr"
AtCAGE.introns <- "AtCAGE_introns.tsr"
AtCAGE.intergenic <- "AtCAGE_intergenic.tsr"

#AtOligo
AtOligo.upstream <- "AtOligo_upstream.tsr"
AtOligo.downstream <- "AtOligo_downstream.tsr"
AtOligo.cds <- "AtOligo_cds.tsr"
AtOligo.introns <- "AtOligo_introns.tsr"
AtOligo.intergenic <- "AtOligo_intergenic.tsr"

#importing the coverage data

AtPEAT.upstream.cov <- read.table(file=AtPEAT.upstream, header=FALSE)
AtPEAT.downstream.cov <- read.table(file=AtPEAT.downstream, header=FALSE)
AtPEAT.cds.cov <- read.table(file=AtPEAT.cds, header=FALSE)
AtPEAT.introns.cov <- read.table(file=AtPEAT.introns, header=FALSE)
AtPEAT.intergenic.cov <- read.table(file=AtPEAT.intergenic, header=FALSE)

AtCAGE.upstream.cov <- read.table(file=AtCAGE.upstream, header=FALSE)
AtCAGE.downstream.cov <- read.table(file=AtCAGE.downstream, header=FALSE)
AtCAGE.cds.cov <- read.table(file=AtCAGE.cds, header=FALSE)
AtCAGE.introns.cov <- read.table(file=AtCAGE.introns, header=FALSE)
AtCAGE.intergenic.cov <- read.table(file=AtCAGE.intergenic, header=FALSE)

AtOligo.upstream.cov <- read.table(file=AtOligo.upstream, header=FALSE)
AtOligo.downstream.cov <- read.table(file=AtOligo.downstream, header=FALSE)
AtOligo.cds.cov <- read.table(file=AtOligo.cds, header=FALSE)
AtOligo.introns.cov <- read.table(file=AtOligo.introns, header=FALSE)
AtOligo.intergenic.cov <- read.table(file=AtOligo.intergenic, header=FALSE)


AtPEAT.cov <- cbind(AtPEAT.downstream.cov[,10], AtPEAT.upstream.cov[,10], AtPEAT.cds.cov[,10], AtPEAT.introns.cov[,10], AtPEAT.intergenic.cov[,4])
colnames(AtPEAT.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
print(head(AtPEAT.cov))
print(tail(AtPEAT.cov))
AtCAGE.cov <- cbind(AtCAGE.downstream.cov[,10], AtCAGE.upstream.cov[,10], AtCAGE.cds.cov[,10], AtCAGE.introns.cov[,10], AtCAGE.intergenic.cov[,4])
colnames(AtCAGE.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtOligo.cov <- cbind(AtOligo.downstream.cov[,10], AtOligo.upstream.cov[,10], AtOligo.cds.cov[,10], AtOligo.introns.cov[,10], AtOligo.intergenic.cov[,4])
colnames(AtOligo.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")

AtPEAT.m <- melt(AtPEAT.cov, value.name="nTags")
AtCAGE.m <- melt(AtCAGE.cov, value.name="nTags")
AtOligo.m <- melt(AtOligo.cov, value.name="nTags")

AtPEAT.m <- AtPEAT.m[,-1]
AtCAGE.m <- AtCAGE.m[,-1]
AtOligo.m <- AtOligo.m[,-1]

colnames(AtPEAT.m) <- c("Region","nTags")
colnames(AtCAGE.m) <- c("Region","nTags")
colnames(AtOligo.m) <- c("Region","nTags")

print(head(AtPEAT.m))
print(head(AtOligo.m))
print(head(AtCAGE.m))

n.total <- sum(as.numeric(AtPEAT.m$nTags))
print(n.total)
p.downstream <- sum(AtPEAT.m[AtPEAT.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtPEAT.m[AtPEAT.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtPEAT.m[AtPEAT.m$Region=="CDS",2])/n.total
p.introns <- sum(AtPEAT.m[AtPEAT.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtPEAT.m[AtPEAT.m$Region=="Intergenic",2])/n.total

AtTSR.df <- data.frame(matrix(NA, nrow=5, ncol=3))
colnames(AtTSR.df) <- c("AtPEAT-tsr","AtCAGE-tsr","AtOligo-tsr")

AtTSR.df[,1] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtCAGE.m$nTags))
p.downstream <- sum(AtCAGE.m[AtCAGE.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtCAGE.m[AtCAGE.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtCAGE.m[AtCAGE.m$Region=="CDS",2])/n.total
p.introns <- sum(AtCAGE.m[AtCAGE.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtCAGE.m[AtCAGE.m$Region=="Intergenic",2])/n.total

AtTSR.df[,2] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtOligo.m$nTags))
p.downstream <- sum(AtOligo.m[AtOligo.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtOligo.m[AtOligo.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtOligo.m[AtOligo.m$Region=="CDS",2])/n.total
p.introns <- sum(AtOligo.m[AtOligo.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtOligo.m[AtOligo.m$Region=="Intergenic",2])/n.total

AtTSR.df[,3] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

AtTSR.df.t <- t(AtTSR.df)

colnames(AtTSR.df.t) <- c("Downstream","Upstream","CDS","Introns","Intergenic")

AtTSR.df.t

AtTSR.df.t.m <- melt(AtTSR.df.t)

colnames(AtTSR.df.t.m) <- c("Experiment", "Region", "TSR_fraction")

#setting the colors (4)
mycols <- c('#FFFD00', '#97CB00', '#3168FF', '#FF0200')

#plotting the figure

p <- ggplot(AtTSR.df.t.m, aes(x=Experiment, y=TSR_fraction, fill=Region)) + geom_bar(stat = 'identity')
p + theme(axis.title.x =NULL)
ggsave("AtTSR_coverage_plot.png")
