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

ChrC.ind <- which(AtPEAT_df$seq=="ChrC")
AtPEAT_df_new <- AtPEAT_df[-ChrC.ind,]
ChrM.ind <- which(AtPEAT_df_new$seq=="ChrM")
AtPEAT_df_final <- AtPEAT_df_new[-ChrM.ind,]

AtPEAT_gr <- makeGRangesFromDataFrame(df=AtPEAT_df_final, seqnames.field="seq", keep.extra.columns=TRUE)
AtCAGE_gr <- makeGRangesFromDataFrame(df=AtCAGE_df, seqnames.field="seq", keep.extra.columns=TRUE)
AtOligo_gr <- makeGRangesFromDataFrame(df=AtOligo_df, seqnames.field="seq", keep.extra.columns=TRUE)

AtGRs <- GRangesList(PEAT=AtPEAT_gr, Oligo=AtOligo_gr, CAGE=AtCAGE_gr)
save(AtGRs, file="AtGR_list.RData")

AtOL <- findOverlapsOfPeaks(AtPEAT_gr, AtOligo_gr, AtCAGE_gr, connectedPeaks = "keepAll", ignore.strand = FALSE)
vennOut <- makeVennDiagram(AtOL, NameOfPeaks = c("PEAT", "Oligo", "CAGE"), totalTest=50000)
new.df <- as.data.frame(vennOut$vennCounts[-1,1:4])
expressionInput <- c(CAGE = new.df[1,4], Oligo = new.df[2,4], PEAT = new.df[4,4], `PEAT&CAGE` = new.df[5,4], `CAGE&Oligo` = new.df[3,4],
`PEAT&Oligo` = new.df[6,4], `CAGE&Oligo&PEAT` = new.df[7,4])
new.input <- fromExpression(expressionInput)

png(filename="TSRintersections.png")
upset(new.input, nsets = 3, nintersects = NA, scale.sets="identity", main.bar.color="navy", sets.bar.color = "navy",  number.angles=30, point.size=3.5, line.size=2, sets.x.label="# of TSRs Identified", mainbar.y.label="Number of Intersecting TSRs", matrix.color = "darkred", text.scale=c(2, 2, 1, 1, 2, 2))

#### What percentage of TSRs are supported by annotation?

AtPEAT_df <- as.data.frame(AtPEAT_gr)
AtCAGE_df <- as.data.frame(AtCAGE_gr)
AtOligo_df <- as.data.frame(AtOligo_gr)

AtPEAT_annot <- AtPEAT_df[is.na(AtPEAT_df$featureID)==FALSE,]
AtPEAT_annot$Annot <- "Gene"
AtPEAT_nonAnnot <- AtPEAT_df[is.na(AtPEAT_df$featureID)==TRUE,]
AtPEAT_nonAnnot$Annot <- "NonGene"

AtCAGE_annot <- AtCAGE_df[is.na(AtCAGE_df$featureID)==FALSE,]
AtCAGE_annot$Annot <- "Gene"
AtCAGE_nonAnnot <- AtCAGE_df[is.na(AtCAGE_df$featureID)==TRUE,]
AtCAGE_nonAnnot$Annot <- "NonGene"

AtOligo_annot <- AtOligo_df[is.na(AtOligo_df$featureID)==FALSE,]
AtOligo_annot$Annot <- "Gene"
AtOligo_nonAnnot <- AtOligo_df[is.na(AtOligo_df$featureID)==TRUE,]
AtOligo_nonAnnot$Annot <- "NonGene"

AtTSR_combined <- rbind(AtPEAT_annot, AtPEAT_nonAnnot, AtCAGE_annot, AtCAGE_nonAnnot, AtOligo_annot, AtOligo_nonAnnot)
AtTSR_new <- as.data.frame(cbind(as.numeric(AtTSR_combined$tsrWidth), as.numeric(AtTSR_combined$nTSSs), as.numeric(AtTSR_combined$shapeIndex), as.character(AtTSR_combined$experiment), as.character(AtTSR_combined$Annot)))
colnames(AtTSR_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment", "Annot")
rownames(AtTSR_new) <- rownames(AtTSR_combined)

df.m <- melt(AtTSR_new, id.vars = c("experiment", "Annot"))
df.m$Annot = factor(df.m$Annot, c("NonGene","Gene"))

df.m$value <- as.numeric(df.m$value)
df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]
df.m.width <- df.m[df.m$variable=="tsrWidth",]
df.m.SI <- df.m[df.m$variable=="shapeIndex",]

p <- ggplot(data = df.m.width, aes(x=Annot, y=value))
q <- p + geom_boxplot(aes(fill=Annot)) + scale_y_continuous(name="TSR width (in base pairs)", limits = c(0,200)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="TSR_width_by_annotation_boxplot.png")

p <- ggplot(df.m[order(df.m$Annot),], aes(x=experiment, fill=Annot)) 
p + geom_bar(position="fill") + scale_y_continuous() +
  labs(x="Experiment", y="Proportion")
ggsave(filename="TSR_width_by_annotation_boxplot.png")

df.m.nTSSs %>% count(annot = factor(Annot), exp = factor(experiment)) %>% 
  ungroup() %>%    # drop if you want percentages per cylinder
  mutate(pct = prop.table(n) * 100)

#### Stability of promoter definitions

At_unique_peaks <- AtOL$uniquePeaks
At_unique_df <- as.data.frame(At_unique_peaks)

threeOL_1 <- Reduce(subsetByOverlaps, list(AtPEAT_gr, AtCAGE_gr, AtOligo_gr))
OL_1_df <- as.data.frame(threeOL_1)
threeOL_2 <- Reduce(subsetByOverlaps, list(AtCAGE_gr, AtOligo_gr, AtPEAT_gr))
OL_2_df <- as.data.frame(threeOL_2)
threeOL_3 <- Reduce(subsetByOverlaps, list(AtOligo_gr, AtPEAT_gr, AtCAGE_gr))
OL_3_df <- as.data.frame(threeOL_3)

OL_combined <- rbind(OL_1_df, OL_2_df, OL_3_df)
OL_new <- as.data.frame(cbind(as.numeric(OL_combined$tsrWidth), as.numeric(OL_combined$nTSSs), as.numeric(OL_combined$shapeIndex), as.character(OL_combined$experiment), as.character(OL_combined$class)))
colnames(OL_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment")
rownames(OL_new) <- rownames(OL_combined)
OL_new$class <- "stable"

At_unique_new <- as.data.frame(cbind(as.numeric(At_unique_df$tsrWidth), as.numeric(At_unique_df$nTSSs), as.numeric(At_unique_df$shapeIndex), as.character(At_unique_df$experiment), as.character(At_unique_df$class)))
colnames(At_unique_new) <- c("tsrWidth","nTSSs", "shapeIndex", "experiment")
rownames(At_unique_new) <- rownames(At_unique_df)
At_unique_new$class <- "unique"

At_intersect <- rbind(OL_new, At_unique_new)
df.m <- melt(At_intersect, id.vars = c("experiment", "class"))
df.m$value <- as.numeric(df.m$value)
df.m.nTSSs <- df.m[df.m$variable=="nTSSs",]

p <- ggplot(data = df.m.nTSSs, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_continuous(name="Tag count per TSR", limits = c(0,500)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="nTSSs_by_stability_experiment_boxplot.png")

df.m.width <- df.m[df.m$variable=="tsrWidth",]

p <- ggplot(data = df.m.width, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_continuous(name="TSR width (in base pairs)", limits = c(0,250)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="TSR_width_by_stability_experiment_boxplot.png")

df.m.SI <- df.m[df.m$variable=="shapeIndex",]

p <- ggplot(data = df.m.SI, aes(x=variable, y=value))
q <- p + geom_boxplot(aes(fill=class)) + scale_y_reverse(name="TSR Shape Index", limits=c(2,-4)) + scale_x_discrete(name="Experiment")
q + facet_wrap( ~ experiment, scales="free")
ggsave(filename="TSR_shape_by_stability_experiment_boxplot.png")
