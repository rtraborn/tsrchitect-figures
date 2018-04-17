#Parses the output of bedtools coverage and creates stacked barplots reflecting the genomic distribution of TSSs

library(ggplot2)
library(reshape2)
setwd("/projects/TSRplantsV2.0/At/AtCoveragePlots")

#AtPEAT
AtPEAT.upstream <- "AtPEAT_upstream.coverage"
AtPEAT.downstream <- "AtPEAT_downstream.coverage"
AtPEAT.cds <- "AtPEAT_cds.coverage"
AtPEAT.introns <- "AtPEAT_introns.coverage"
AtPEAT.intergenic <- "AtPEAT_intergenic.coverage"

AtPEAT_t.upstream <- "AtPEAT_t_upstream.coverage"
AtPEAT_t.downstream <- "AtPEAT_t_downstream.coverage"
AtPEAT_t.cds <- "AtPEAT_t_cds.coverage"
AtPEAT_t.introns <- "AtPEAT_t_introns.coverage"
AtPEAT_t.intergenic <- "AtPEAT_t_intergenic.coverage"

#AtCAGE
AtCAGE.upstream <- "AtCAGE_upstream.coverage"
AtCAGE.downstream <- "AtCAGE_downstream.coverage"
AtCAGE.cds <- "AtCAGE_cds.coverage"
AtCAGE.introns <- "AtCAGE_introns.coverage"
AtCAGE.intergenic <- "AtCAGE_intergenic.coverage"

AtCAGE_t.upstream <- "AtCAGE_t_upstream.coverage"
AtCAGE_t.downstream <- "AtCAGE_t_downstream.coverage"
AtCAGE_t.cds <- "AtCAGE_t_cds.coverage"
AtCAGE_t.introns <- "AtCAGE_t_introns.coverage"
AtCAGE_t.intergenic <- "AtCAGE_t_intergenic.coverage"

#AtOligo
AtOligo.upstream <- "AtOligo_upstream.coverage"
AtOligo.downstream <- "AtOligo_downstream.coverage"
AtOligo.cds <- "AtOligo_cds.coverage"
AtOligo.introns <- "AtOligo_introns.coverage"
AtOligo.intergenic <- "AtOligo_intergenic.coverage"

AtOligo_t.upstream <- "AtOligo_t_upstream.coverage"
AtOligo_t.downstream <- "AtOligo_t_downstream.coverage"
AtOligo_t.cds <- "AtOligo_t_cds.coverage"
AtOligo_t.introns <- "AtOligo_t_introns.coverage"
AtOligo_t.intergenic <- "AtOligo_t_intergenic.coverage"

#importing the coverage data

AtPEAT.upstream.cov <- read.table(file=AtPEAT.upstream, header=FALSE)
AtPEAT.downstream.cov <- read.table(file=AtPEAT.downstream, header=FALSE)
AtPEAT.cds.cov <- read.table(file=AtPEAT.cds, header=FALSE)
AtPEAT.introns.cov <- read.table(file=AtPEAT.introns, header=FALSE)
AtPEAT.intergenic.cov <- read.table(file=AtPEAT.intergenic, header=FALSE)

AtPEAT.upstream.t <- read.table(file=AtPEAT_t.upstream, header=FALSE)
AtPEAT.downstream.t <- read.table(file=AtPEAT_t.downstream, header=FALSE)
AtPEAT.cds.t <- read.table(file=AtPEAT_t.cds, header=FALSE)
AtPEAT.introns.t <- read.table(file=AtPEAT_t.introns, header=FALSE)
AtPEAT.intergenic.t <- read.table(file=AtPEAT_t.intergenic, header=FALSE)

AtCAGE.upstream.cov <- read.table(file=AtCAGE.upstream, header=FALSE)
AtCAGE.downstream.cov <- read.table(file=AtCAGE.downstream, header=FALSE)
AtCAGE.cds.cov <- read.table(file=AtCAGE.cds, header=FALSE)
AtCAGE.introns.cov <- read.table(file=AtCAGE.introns, header=FALSE)
AtCAGE.intergenic.cov <- read.table(file=AtCAGE.intergenic, header=FALSE)

AtCAGE.upstream.t <- read.table(file=AtCAGE_t.upstream, header=FALSE)
AtCAGE.downstream.t <- read.table(file=AtCAGE_t.downstream, header=FALSE)
AtCAGE.cds.t <- read.table(file=AtCAGE_t.cds, header=FALSE)
AtCAGE.introns.t <- read.table(file=AtCAGE_t.introns, header=FALSE)
AtCAGE.intergenic.t <- read.table(file=AtCAGE_t.intergenic, header=FALSE)

AtOligo.upstream.cov <- read.table(file=AtOligo.upstream, header=FALSE)
AtOligo.downstream.cov <- read.table(file=AtOligo.downstream, header=FALSE)
AtOligo.cds.cov <- read.table(file=AtOligo.cds, header=FALSE)
AtOligo.introns.cov <- read.table(file=AtOligo.introns, header=FALSE)
AtOligo.intergenic.cov <- read.table(file=AtOligo.intergenic, header=FALSE)

AtOligo.upstream.t <- read.table(file=AtOligo_t.upstream, header=FALSE)
AtOligo.downstream.t <- read.table(file=AtOligo_t.downstream, header=FALSE)
AtOligo.cds.t <- read.table(file=AtOligo_t.cds, header=FALSE)
AtOligo.introns.t <- read.table(file=AtOligo_t.introns, header=FALSE)
AtOligo.intergenic.t <- read.table(file=AtOligo_t.intergenic, header=FALSE)

AtPEAT.cov <- cbind(AtPEAT.downstream.cov[,10], AtPEAT.upstream.cov[,10], AtPEAT.cds.cov[,10], AtPEAT.introns.cov[,10], AtPEAT.intergenic.cov[,4])
colnames(AtPEAT.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtPEAT.t <- cbind(AtPEAT.downstream.t[,10], AtPEAT.upstream.t[,10], AtPEAT.cds.t[,10], AtPEAT.introns.t[,10], AtPEAT.intergenic.t[,4])
colnames(AtPEAT.t) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtCAGE.cov <- cbind(AtCAGE.downstream.cov[,10], AtCAGE.upstream.cov[,10], AtCAGE.cds.cov[,10], AtCAGE.introns.cov[,10], AtCAGE.intergenic.cov[,4])
colnames(AtCAGE.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtCAGE.t <- cbind(AtCAGE.downstream.t[,10], AtCAGE.upstream.t[,10], AtCAGE.cds.t[,10], AtCAGE.introns.t[,10], AtCAGE.intergenic.t[,4])
colnames(AtCAGE.t) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtOligo.cov <- cbind(AtOligo.downstream.cov[,10], AtOligo.upstream.cov[,10], AtOligo.cds.cov[,10], AtOligo.introns.cov[,10], AtOligo.intergenic.cov[,4])
colnames(AtOligo.cov) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")
AtOligo.t <- cbind(AtOligo.downstream.t[,10], AtOligo.upstream.t[,10], AtOligo.cds.t[,10], AtOligo.introns.t[,10], AtOligo.intergenic.t[,4])
colnames(AtOligo.t) <- c("Downstream", "Upstream", "CDS", "Introns", "Intergenic")

AtPEAT.m <- melt(AtPEAT.cov, value.name="nTags")
AtCAGE.m <- melt(AtCAGE.cov, value.name="nTags")
AtOligo.m <- melt(AtOligo.cov, value.name="nTags")
AtPEAT.t.m <- melt(AtPEAT.t, value.name="nTags")
AtCAGE.t.m <- melt(AtCAGE.t, value.name="nTags")
AtOligo.t.m <- melt(AtOligo.t, value.name="nTags")

AtPEAT.m <- AtPEAT.m[,-1]
AtCAGE.m <- AtCAGE.m[,-1]
AtOligo.m <- AtOligo.m[,-1]
AtPEAT.t.m <- AtPEAT.t.m[,-1]
AtCAGE.t.m <- AtCAGE.t.m[,-1]
AtOligo.t.m <- AtOligo.t.m[,-1]

colnames(AtPEAT.m) <- c("Region","nTags")
colnames(AtCAGE.m) <- c("Region","nTags")
colnames(AtOligo.m) <- c("Region","nTags")
colnames(AtPEAT.t.m) <- c("Region","nTags")
colnames(AtCAGE.t.m) <- c("Region","nTags")
colnames(AtOligo.t.m) <- c("Region","nTags")

print(head(AtPEAT.m))
print(head(AtPEAT.t.m))
print(head(AtOligo.m))
print(head(AtOligo.t.m))

n.total <- sum(as.numeric(AtPEAT.m$nTags))
print(n.total)
p.downstream <- sum(AtPEAT.m[AtPEAT.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtPEAT.m[AtPEAT.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtPEAT.m[AtPEAT.m$Region=="CDS",2])/n.total
p.introns <- sum(AtPEAT.m[AtPEAT.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtPEAT.m[AtPEAT.m$Region=="Intergenic",2])/n.total

AtPEAT.df <- data.frame(matrix(NA, nrow=5, ncol=2))
colnames(AtPEAT.df) <- c("AtPEAT-all","AtPEAT-filtered")

AtPEAT.df[,1] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtPEAT.t.m$nTags))
print(n.total)
p.downstream <- sum(AtPEAT.t.m[AtPEAT.t.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtPEAT.t.m[AtPEAT.t.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtPEAT.t.m[AtPEAT.t.m$Region=="CDS",2])/n.total
p.introns <- sum(AtPEAT.t.m[AtPEAT.t.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtPEAT.t.m[AtPEAT.t.m$Region=="Intergenic",2])/n.total

AtPEAT.df[,2] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtCAGE.m$nTags))
p.downstream <- sum(AtCAGE.m[AtCAGE.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtCAGE.m[AtCAGE.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtCAGE.m[AtCAGE.m$Region=="CDS",2])/n.total
p.introns <- sum(AtCAGE.m[AtCAGE.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtCAGE.m[AtCAGE.m$Region=="Intergenic",2])/n.total

AtCAGE.df <- data.frame(matrix(NA, nrow=5, ncol=2))
colnames(AtCAGE.df) <- c("AtCAGE-all", "AtCAGE-filtered")

AtCAGE.df[,1] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtCAGE.t.m$nTags))
p.downstream <- sum(AtCAGE.t.m[AtCAGE.t.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtCAGE.t.m[AtCAGE.t.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtCAGE.t.m[AtCAGE.t.m$Region=="CDS",2])/n.total
p.introns <- sum(AtCAGE.t.m[AtCAGE.t.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtCAGE.t.m[AtCAGE.t.m$Region=="Intergenic",2])/n.total

AtCAGE.df[,2] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtOligo.m$nTags))
p.downstream <- sum(AtOligo.m[AtOligo.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtOligo.m[AtOligo.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtOligo.m[AtOligo.m$Region=="CDS",2])/n.total
p.introns <- sum(AtOligo.m[AtOligo.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtOligo.m[AtOligo.m$Region=="Intergenic",2])/n.total

AtOligo.df <- data.frame(matrix(NA, nrow=5, ncol=2))
colnames(AtOligo.df) <- c("AtOligo-all", "AtOligo-filtered")

AtOligo.df[,1] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

n.total <- sum(as.numeric(AtOligo.t.m$nTags))
p.downstream <- sum(AtOligo.t.m[AtOligo.t.m$Region=="Downstream",2])/n.total
p.upstream <- sum(AtOligo.t.m[AtOligo.t.m$Region=="Upstream",2])/n.total
p.cds <- sum(AtOligo.t.m[AtOligo.t.m$Region=="CDS",2])/n.total
p.introns <- sum(AtOligo.t.m[AtOligo.t.m$Region=="Introns",2])/n.total
p.intergenic <- sum(AtOligo.t.m[AtOligo.t.m$Region=="Intergenic",2])/n.total

AtOligo.df[,2] <- c(p.upstream, p.downstream, p.cds, p.introns, p.intergenic)

AtPEAT.df.t <- t(AtPEAT.df)
AtCAGE.df.t <- t(AtCAGE.df)
AtOligo.df.t <- t(AtOligo.df)

colnames(AtPEAT.df.t) <- c("Downstream","Upstream","CDS","Introns","Intergenic")
colnames(AtCAGE.df.t) <- c("Downstream","Upstream","CDS","Introns","Intergenic")
colnames(AtOligo.df.t) <- c("Downstream","Upstream","CDS","Introns","Intergenic")

AtPEAT.df.t
AtCAGE.df.t
AtOligo.df.t

AtPEAT.df.t.m <- melt(AtPEAT.df.t)
AtCAGE.df.t.m <- melt(AtCAGE.df.t)
AtOligo.df.t.m <- melt(AtOligo.df.t)

colnames(AtPEAT.df.t.m) <- c("Experiment", "Region", "TSS_fraction")
colnames(AtCAGE.df.t.m) <- c("Experiment", "Region", "TSS_fraction")
colnames(AtOligo.df.t.m) <- c("Experiment", "Region", "TSS_fraction")

#setting the colors (4)
mycols <- c('#FFFD00', '#97CB00', '#3168FF', '#FF0200')

#plotting the figure

p <- ggplot(AtPEAT.df.t.m, aes(x=Experiment, y=TSS_fraction, fill=Region)) + geom_bar(stat = 'identity')
p + theme(axis.title.x =NULL)
ggsave("AtPEAT_coverage_plot.png")

p <- ggplot(AtCAGE.df.t.m, aes(x=Experiment, y=TSS_fraction, fill=Region)) + geom_bar(stat = 'identity')
p + theme(axis.title.x =NULL)
ggsave("AtCAGE_coverage_plot.png")

p <- ggplot(AtOligo.df.t.m, aes(x=Experiment, y=TSS_fraction, fill=Region)) + geom_bar(stat = 'identity')
p + theme(axis.title.x =NULL)
ggsave("AtOligo_coverage_plot.png")
