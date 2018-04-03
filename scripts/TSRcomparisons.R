tsrBinary <- c("~/Desktop/TSRscripts/R_binaries/AtGR_list.RData")
myWD <- c("~/Desktop/TSRscripts/R_binaries/")

setwd(myWD)

library(GenomicRanges)
library(UpSetR)
library(ChIPpeakAnno)
library(reshape2)
library(ggplot2)
library(MASS)
library(viridis)
theme_set(theme_bw(base_size = 16))

load(tsrBinary)

AtPEAT_gr <- AtGRs$PEAT
AtCAGE_gr <- AtGRs$CAGE
AtOligo_gr <- AtGRs$Oligo

PEAT_df <- data.frame(AtPEAT_gr$nTSSs, AtPEAT_gr$tsrWidth)
colnames(PEAT_df) <- c("nTSSs_PEAT", "tsrWidth_PEAT")
f <- ggplot(PEAT_df, aes(nTSSs_PEAT, tsrWidth_PEAT))
f + geom_jitter()

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

PEAT_df$density <- get_density(PEAT_df$x, PEAT_df$y)
ggplot(PEAT_df) + geom_point(aes(nTSSs_PEAT, tsrWidth_PEAT, color = density)) + scale_color_viridis() + scale_x_continuous(limits=c(0,250)) + scale_y_continuous(limits=c(0,250))

##### Comparing TSRs of widths of 1 with others ######

AtPEAT_single <- as.data.frame(AtPEAT_gr[AtPEAT_gr$tsrWidth==1,])
AtCAGE_single <- as.data.frame(AtCAGE_gr[AtCAGE_gr$tsrWidth==1,])
AtOligo_single <- as.data.frame(AtOligo_gr[AtOligo_gr$tsrWidth==1,])

AtOligo_multiple <- as.data.frame(AtOligo_gr[AtOligo_gr$tsrWidth>1,])
AtCAGE_multiple <- as.data.frame(AtCAGE_gr[AtCAGE_gr$tsrWidth>1,])
AtPEAT_multiple <- as.data.frame(AtPEAT_gr[AtPEAT_gr$tsrWidth>1,])

AtTSR_single <- rbind(AtPEAT_single, AtCAGE_single, AtOligo_single)
AtTSR_multiple <- rbind(AtPEAT_multiple, AtCAGE_multiple, AtOligo_multiple)

AtTSR_single$class <- "Single"
AtTSR_multiple$class <- "Multiple"
AtTSR_combined <- rbind(AtTSR_single, AtTSR_multiple)

AtTSR_new <- as.data.frame(cbind(as.numeric(AtTSR_combined$tsrWidth), as.numeric(AtTSR_combined$nTSSs), as.numeric(AtTSR_combined$shapeIndex), as.character(AtTSR_combined$experiment), as.character(AtTSR_combined$class)))
colnames(AtTSR_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment", "class")
rownames(AtTSR_new) <- rownames(AtTSR_combined)

df.m <- melt(AtTSR_new, id.vars = c("experiment", "class"))
df.m$value <- as.numeric(df.m$value)
df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]

p <- ggplot(data = df.m.nTSSs, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_continuous(name="Tag count per TSR", limits = c(0,500)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="nTSSs_by_experiment_boxplot.png")

#### Comparing Peaked and Broad promoters

AtPEAT_peaked <- as.data.frame(AtPEAT_gr[AtPEAT_gr$shapeIndex==2,])
AtCAGE_peaked <- as.data.frame(AtCAGE_gr[AtCAGE_gr$shapeIndex==2,])
AtOligo_peaked <- as.data.frame(AtOligo_gr[AtOligo_gr$shapeIndex==2,])

AtPEAT_broad <- as.data.frame(AtPEAT_gr[AtPEAT_gr$shapeIndex<=-1,])
AtCAGE_broad <- as.data.frame(AtCAGE_gr[AtCAGE_gr$shapeIndex<=-1,])
AtOligo_broad <- as.data.frame(AtOligo_gr[AtOligo_gr$shapeIndex<=-1,])

AtPEAT_unclass <- as.data.frame(AtPEAT_gr[AtPEAT_gr$shapeIndex>-1 & AtPEAT_gr$shapeIndex < 2,])
AtCAGE_unclass <- as.data.frame(AtCAGE_gr[AtCAGE_gr$shapeIndex>-1 & AtCAGE_gr$shapeIndex < 2,])
AtOligo_unclass <- as.data.frame(AtOligo_gr[AtOligo_gr$shapeIndex>-1 & AtOligo_gr$shapeIndex < 2,])

AtTSR_peaked <- rbind(AtPEAT_peaked, AtCAGE_peaked, AtOligo_peaked)
AtTSR_broad <- rbind(AtPEAT_broad, AtCAGE_broad, AtOligo_broad)
AtTSR_unclass <- rbind(AtPEAT_unclass, AtCAGE_unclass, AtOligo_unclass)

AtTSR_peaked$class <- "Peaked"
AtTSR_broad$class <- "Broad"
AtTSR_unclass$class <- "Unclassified"

AtTSR_combined <- rbind(AtTSR_peaked, AtTSR_broad, AtTSR_unclass)

AtTSR_new <- as.data.frame(cbind(as.numeric(AtTSR_combined$tsrWidth), as.numeric(AtTSR_combined$nTSSs), as.numeric(AtTSR_combined$shapeIndex), as.character(AtTSR_combined$experiment), as.character(AtTSR_combined$class)))
colnames(AtTSR_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment", "class")
rownames(AtTSR_new) <- rownames(AtTSR_combined)

df.m <- melt(AtTSR_new, id.vars = c("experiment", "class"))
df.m$value <- as.numeric(df.m$value)

df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]
df.m.nTSSs$class = factor(df.m.nTSSs$class, c("Broad","Unclassified","Peaked"))

p <- ggplot(data = df.m.nTSSs, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_continuous(name="Tag count per TSR", limits = c(0,500)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="nTSSs_by_shape_experiment_boxplot.png")

df.m.width <- df.m[df.m$variable=="tsrWidth",]
df.m.width$class = factor(df.m.width$class, c("Broad","Unclassified","Peaked"))

p <- ggplot(data = df.m.width, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_continuous(name="TSR width in base pairs", limits = c(0,100)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="tsrWidths_by_shape_experiment_boxplot.png")

#### Promoter shape by expression quantile

AtPEAT_df <- as.data.frame(AtPEAT_gr)
AtPEAT.quant <- quantile(AtPEAT_gr$nTSSs, c(0,0.25,0.50,0.75,1.0))
AtPEAT_df$expressionClass <- "Bottom_quantile"
second.ind <- which(AtPEAT_df$nTSSs>=AtPEAT.quant[2])
AtPEAT_df$expressionClass[second.ind] <- "Second_quantile"
third.ind <- which(AtPEAT_df$nTSSs>=AtPEAT.quant[3])
AtPEAT_df$expressionClass[third.ind] <- "Third_quantile"
top.ind <- which(AtPEAT_df$nTSSs>=AtPEAT.quant[4])
AtPEAT_df$expressionClass[top.ind] <- "Top_quantile"

AtCAGE_df <- as.data.frame(AtCAGE_gr)
AtCAGE.quant <- quantile(AtCAGE_gr$nTSSs, c(0,0.25,0.50,0.75,1.0))
AtCAGE_df$expressionClass <- "Bottom_quantile"
second.ind <- which(AtCAGE_df$nTSSs>=AtCAGE.quant[2])
AtCAGE_df$expressionClass[second.ind] <- "Second_quantile"
third.ind <- which(AtCAGE_df$nTSSs>=AtCAGE.quant[3])
AtCAGE_df$expressionClass[third.ind] <- "Third_quantile"
top.ind <- which(AtCAGE_df$nTSSs>=AtCAGE.quant[4])
AtCAGE_df$expressionClass[top.ind] <- "Top_quantile"

AtOligo_df <- as.data.frame(AtOligo_gr)
AtOligo.quant <- quantile(AtOligo_gr$nTSSs, c(0,0.25,0.50,0.75,1.0))
AtOligo_df$expressionClass <- "Bottom_quantile"
second.ind <- which(AtOligo_df$nTSSs>=AtOligo.quant[2])
AtOligo_df$expressionClass[second.ind] <- "Second_quantile"
third.ind <- which(AtOligo_df$nTSSs>=AtOligo.quant[3])
AtOligo_df$expressionClass[third.ind] <- "Third_quantile"
top.ind <- which(AtOligo_df$nTSSs>=AtOligo.quant[4])
AtOligo_df$expressionClass[top.ind] <- "Top_quantile"

AtTSR_combined <- rbind(AtPEAT_df, AtOligo_df, AtCAGE_df)
AtTSR_new <- as.data.frame(cbind(as.numeric(AtTSR_combined$tsrWidth), as.numeric(AtTSR_combined$nTSSs), as.numeric(AtTSR_combined$shapeIndex), as.character(AtTSR_combined$experiment), as.character(AtTSR_combined$expressionClass)))
colnames(AtTSR_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment", "expressionClass")
rownames(AtTSR_new) <- rownames(AtTSR_combined)

df.m <- melt(AtTSR_new, id.vars = c("experiment", "expressionClass"))
df.m$value <- as.numeric(df.m$value)

df.m.SI <- df.m[df.m$variable=="shapeIndex",]
df.m.SI$expressionClass = factor(df.m.SI$expressionClass, c("Bottom_quantile","Second_quantile","Third_quantile","Top_quantile"))

p <- ggplot(data = df.m.SI, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=expressionClass)) + scale_y_reverse(name="TSR Shape Index (SI)") + scale_x_discrete(name="TSR Expression Quantile")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="tsrSI_by_expressionClass_experiment_boxplot.png")
