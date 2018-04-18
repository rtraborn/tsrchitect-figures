tsrBinary <- c("~/Desktop/TSRscripts/R_binaries/AtTSR_data.RData")
myWD <- c("~/Desktop/TSRscripts/R_binaries/")

setwd(myWD)

library(GenomicRanges)
library(UpSetR)
library(ChIPpeakAnno)
library(dplyr)

load(tsrBinary)

AtOligo_df <- AtTSR_data$AtOligo
AtCAGE_df <- AtTSR_data$AtCAGE
AtPEAT_df <- AtTSR_data$AtPEAT

head(AtPEAT_df)

AtPEAT_gr <- makeGRangesFromDataFrame(df=AtPEAT_df, seqnames.field="seq", keep.extra.columns=TRUE)
AtCAGE_gr <- makeGRangesFromDataFrame(df=AtCAGE_df, seqnames.field="seq", keep.extra.columns=TRUE)
AtOligo_gr <- makeGRangesFromDataFrame(df=AtOligo_df, seqnames.field="seq", keep.extra.columns=TRUE)

AtGRs <- GRangesList(PEAT=AtPEAT_gr, Oligo=AtOligo_gr, CAGE=AtCAGE_gr)
save(AtGRs, file="AtGR_list.RData")

#### What is the distribution of TSR widths within the samples?

AtPEAT_df <- as.data.frame(AtPEAT_gr)
AtCAGE_df <- as.data.frame(AtCAGE_gr)
AtOligo_df <- as.data.frame(AtOligo_gr)

AtPEAT_df$experiment <- "AtPEAT"
AtCAGE_df$experiment <- "AtCAGE"
AtOligo_df$experiment <- "AtOligo"

AtTSR_combined <- rbind(AtPEAT_df, AtCAGE_df, AtOligo_df)

AtTSR_new <- as.data.frame(cbind(as.numeric(AtTSR_combined$tsrWdth), as.numeric(AtTSR_combined$nTSSs), as.numeric(AtTSR_combined$nTAGs), as.numeric(AtTSR_combined$tsrMSI), as.character(AtTSR_combined$experiment)))

colnames(AtTSR_new) <- c("tsrWidth","nTSSs", "nTAGs", "M_SI", "experiment")
rownames(AtTSR_new) <- rownames(AtTSR_combined)

df.m <- melt(AtTSR_new, id.vars = c("experiment"))
df.m$value <- as.numeric(df.m$value)
colnames(df.m) <- c("Experiment","variable", "value")
head(df.m)

df.m$value <- as.numeric(df.m$value)
df.m.nTAGs <- df.m[df.m$variable=="nTAGs",]
df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]
df.m.width <- df.m[df.m$variable=="tsrWidth",]
df.m.SI <- df.m[df.m$variable=="M_SI",]

head(df.m)
dim(df.m)
head(df.m.nTSSs)

######### TSR width ########

p <- ggplot(data = df.m.width, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="TSR width (in base pairs)", limits = c(0,100)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                                                       axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen"))
ggsave(filename="TSR_widths_by_experiment_violin.png")

summary(df.m.width)

######### TSR expression ########

p <- ggplot(data = df.m.nTAGs, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="TSR expression (number of tags)", limits = c(0,500)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                                                       axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen"))
ggsave(filename="TSR_expression_by_experiment_violin.png")

summary(df.m.nTAGs)

######### TSR expression ########

p <- ggplot(data = df.m.SI, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="modified TSR Shape Index", limits = c(0,1)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                           axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen"))
ggsave(filename="TSR_mSI_by_experiment_violin.png")

summary(df.m.SI)

######### number of unique TSS sites per TSR ########

p <- ggplot(data = df.m.nTSSs, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="number of distinct TSS positions per TSR", limits = c(0,15)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                  axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen"))
ggsave(filename="TSR_nTSSs_by_experiment_violin.png")

summary(df.m.nTSSs)

#### What are the consensus promoters? #####

# subsetting all three experiments by overlaps

AtTSR_ol <- BiocGenerics::Reduce(subsetByOverlaps, list(AtPEAT_gr, AtOligo_gr, AtCAGE_gr))
AtTSR_ol #6803 overlaps

AtTSRol_df <- as.data.frame(cbind(as.numeric(AtTSR_ol$tsrWdth), as.numeric(AtTSR_ol$nTSSs), as.numeric(AtTSR_ol$nTAGs), as.numeric(AtTSR_ol$tsrMSI)))
AtTSRol_df$experiment <- "consensus"
colnames(AtTSRol_df) <- c("tsrWidth","nTSSs", "nTAGs", "M_SI", "experiment")
rownames(AtTSRol_df) <- rownames(as.data.frame(AtTSR_ol))

head(AtTSRol_df)
dim(AtTSRol_df)

colnames(AtTSRol_df) <- c("tsrWidth","nTSSs", "nTAGs", "M_SI", "experiment")
rownames(AtTSRol_df) <- rownames(as.data.frame(AtTSR_ol))

df.ol.m <- melt(AtTSRol_df, id.vars = c("experiment"))
df.ol.m$value <- as.numeric(df.ol.m$value)
colnames(df.ol.m) <- c("Experiment","variable", "value")
head(df.ol.m)

df.ol.m$value <- as.numeric(df.ol.m$value)
df.ol.m.nTAGs <- df.ol.m[df.ol.m$variable=="nTAGs",]
df.ol.m.nTSSs <- df.ol.m[df.ol.m$variable=="nTSSs",]
df.ol.m.width <- df.ol.m[df.ol.m$variable=="tsrWidth",]
df.ol.m.SI <- df.ol.m[df.ol.m$variable=="M_SI",]

summary(df.ol.m$value)
summary(df.ol.m.nTAGs)
summary(df.ol.m.nTSSs)
summary(df.ol.m.width)
summary(df.ol.m.SI)

df.m <- rbind(df.m,df.ol.m)

df.m$value <- as.numeric(df.m$value)
df.m.nTAGs <- df.m[df.m$variable=="nTAGs",]
df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]
df.m.width <- df.m[df.m$variable=="tsrWidth",]
df.m.SI <- df.m[df.m$variable=="M_SI",]

head(df.m)
dim(df.m)
head(df.m.nTSSs)

######### TSR width with consensus ########

p <- ggplot(data = df.m.width, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="TSR width (in base pairs)", limits = c(0,100)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                     axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen","maroon"))
ggsave(filename="TSR_widths_by_experiment_consensus_violin.png")

summary(df.m.width)

######### TSR expression with consensus ########

p <- ggplot(data = df.m.nTAGs, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="TSR expression (number of tags)", limits = c(0,500)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                           axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen","maroon"))
ggsave(filename="TSR_expression_by_experiment_consensus_violin.png")

summary(df.m.nTAGs)

######### TSR expression with consensus ########

p <- ggplot(data = df.m.SI, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="modified TSR Shape Index", limits = c(0,1)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                  axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen","maroon"))
ggsave(filename="TSR_mSI_by_experiment_consensus_violin.png")

summary(df.m.SI)

######### number of unique TSS sites per TSR with consensus ########

p <- ggplot(data =  vnTSSs, aes(x=Experiment, y=value, fill=Experiment)) + coord_flip()
p <- p + geom_violin(alpha=0.6) + scale_y_continuous(name="number of distinct TSS positions per TSR", limits = c(0,25)) + scale_x_discrete(name="Experiment") + theme_bw() + theme(axis.text=element_text(size=18),
                                                                                                                                                                                   axis.title=element_text(size=24,face="bold"))
p + scale_fill_manual(values=c("navy","orange","darkgreen","maroon"))
ggsave(filename="TSR_nTSSs_by_experiment_consensus_violin.png")

######### Selection of broad and peaked TSRs ########

broad.df <- subset(AtTSRol_df, M_SI <= 0.1)
broad.df$shapeClass <- "broad"
dim(broad.df)

peaked.df <- subset(AtTSRol_df, M_SI >= 0.9) 
peaked.df$shapeClass <- "peaked"
dim(peaked.df)

unclass.df <- subset(AtTSRol_df, M_SI > 0.1 & M_SI < 0.9)
unclass.df$shapeClass <- "unclassified"
dim(unclass.df)

consensus.df <- rbind(peaked.df, broad.df, unclass.df)
head(consensus.df)

cons.m <- melt(consensus.df, id.vars = c("shapeClass"))
head(cons.m)
cons.m$value <- as.numeric(cons.m$value)
colnames(cons.m) <- c("shapeClass","variable", "value")
head(cons.m)

df.m.nTAGs <- cons.m[cons.m$variable=="nTAGs",]
df.m.nTAGs$shapeClass = factor(df.m.nTAGs$shapeClass, c("broad","unclassified","peaked"))

p <- ggplot(data = df.m.nTAGs, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=shapeClass)) + scale_y_continuous(name="Tag count per TSR", limits = c(0,500)) + theme_bw() + scale_x_discrete(name="Consensus TSRs") + theme(axis.text=element_text(size=18),
                                                                                                                                                          axis.title=element_text(size=24,face="bold"))
q + facet_wrap( ~ shapeClass, scales="free")
ggsave(filename="nTSSs_by_shape_consensus_boxplot.png")

