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

