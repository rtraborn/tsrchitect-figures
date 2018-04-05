setwd("~/Desktop/tsrchitect-figures/dist_datasets")
cat("Loading the libraries\n")

library("reshape2")
library("ggplot2")

##### AtCAGE ########

AtCAGE_coverage <- read.table("AtCAGE_tss_all_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtCAGE_coverage)
colnames(AtCAGE_coverage) <- c("position","AtCAGE-all","AtCAGE-filtered","AtCAGE-below-threshold")
AtCAGE_cov <- AtCAGE_coverage[,-1]
AtCAGE_scaled <- as.data.frame(scale(AtCAGE_cov, center=FALSE, scale=colSums(AtCAGE_cov)))
head(AtCAGE_scaled)
AtCAGE_scaled$position  <- AtCAGE_coverage$position
AtCAGE_cov_m <- melt(AtCAGE_scaled, id.vars = "position")
#head(AtCAGE_cov_m)
colnames(AtCAGE_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtCAGE_cov_m, aes(x=position, y=coverage, by=TSSset, linetype=TSSset))
p + geom_line(color="navy", size=1, alpha=0.5) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16)
ggsave("AtCAGE_coverage_v_TSSs.png")

##### AtPEAT ########

AtPEAT_coverage <- read.table("AtPEATtss_all_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtPEAT_coverage)
colnames(AtPEAT_coverage) <- c("position","AtPEAT-all","AtPEAT-filtered","AtPEAT-below-threshold")
AtPEAT_cov <- AtPEAT_coverage[,-1]
AtPEAT_scaled <- as.data.frame(scale(AtPEAT_cov, center=FALSE, scale=colSums(AtPEAT_cov)))
#head(AtAtPEAT_scaled)
AtPEAT_scaled$position  <- AtPEAT_coverage$position
AtPEAT_cov_m <- melt(AtPEAT_scaled, id.vars = "position")
#head(AtPEAT_cov_m)
colnames(AtPEAT_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtPEAT_cov_m, aes(x=position, y=coverage, by=TSSset, linetype=TSSset))
p + geom_line(color="darkgreen", size=1, alpha=0.7) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16)
ggsave("AtPEAT_coverage_v_TSSs.png")

#### AtOligo ####

AtOligo_coverage <- read.table("AtOligo_tss_all_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtOligo_coverage)
colnames(AtOligo_coverage) <- c("position","AtOligo-all","AtOligo-filtered","AtOligo-below-threshold")
AtOligo_cov <- AtOligo_coverage[,-1]
AtOligo_scaled <- as.data.frame(scale(AtOligo_cov, center=FALSE, scale=colSums(AtOligo_cov)))
#head(AtAtOligo_scaled)
AtOligo_scaled$position  <- AtOligo_coverage$position
AtOligo_cov_m <- melt(AtOligo_scaled, id.vars = "position")
#head(AtOligo_cov_m)
colnames(AtOligo_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtOligo_cov_m, aes(x=position, y=coverage, by=TSSset, linetype=TSSset))
p + geom_line(color="#E69F00", size=1, alpha=0.8) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16) #color is a light orange
ggsave("AtOligo_coverage_v_TSSs.png")

#### All filtered datasets ####

AtTSS_coverage <- read.table("AtTSS_filtered_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtTSS_coverage)
colnames(AtTSS_coverage) <- c("position","AtCAGE-filtered","AtPEAT-filtered","AtOligo-filtered")
AtTSS_cov <- AtTSS_coverage[,-1]
AtTSS_scaled <- as.data.frame(scale(AtTSS_cov, center=FALSE, scale=colSums(AtTSS_cov)))
#head(AtAtTSS_scaled)
AtTSS_scaled$position  <- AtTSS_coverage$position
AtTSS_cov_m <- melt(AtTSS_scaled, id.vars = "position")
#head(AtTSS_cov_m)
colnames(AtTSS_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtTSS_cov_m, aes(x=position, y=coverage, by=TSSset, color=TSSset))
p + geom_line(size=1, alpha=0.4) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16) + scale_colour_manual(values=c("navy","darkgreen","#E69F00"))
ggsave("AtTSS_filtered_coverage_v_TSSs.png")

#### All below-threshold datasets ####

AtTSS_coverage <- read.table("AtTSS_belowThresh_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtTSS_coverage)
colnames(AtTSS_coverage) <- c("position","AtCAGE-belowThresh","AtPEAT-belowThresh","AtOligo-belowThresh")
AtTSS_cov <- AtTSS_coverage[,-1]
AtTSS_scaled <- as.data.frame(scale(AtTSS_cov, center=FALSE, scale=colSums(AtTSS_cov)))
#head(AtAtTSS_scaled)
AtTSS_scaled$position  <- AtTSS_coverage$position
AtTSS_cov_m <- melt(AtTSS_scaled, id.vars = "position")
#head(AtTSS_cov_m)
colnames(AtTSS_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtTSS_cov_m, aes(x=position, y=coverage, by=TSSset, color=TSSset))
p + geom_line(size=1, alpha=0.4) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16) + scale_colour_manual(values=c("navy","darkgreen","#E69F00"))
ggsave("AtTSS_belowThresh_coverage_v_TSSs.png")

#### All raw-threshold datasets ####

AtTSS_coverage <- read.table("AtTSS_tss_raw_dist_norm_1kb.txt", header=FALSE, skip=1)
#head(AtTSS_coverage)
colnames(AtTSS_coverage) <- c("position","AtCAGE-all","AtPEAT-all","AtOligo-all")
AtTSS_cov <- AtTSS_coverage[,-1]
AtTSS_scaled <- as.data.frame(scale(AtTSS_cov, center=FALSE, scale=colSums(AtTSS_cov)))
#head(AtAtTSS_scaled)
AtTSS_scaled$position  <- AtTSS_coverage$position
AtTSS_cov_m <- melt(AtTSS_scaled, id.vars = "position")
#head(AtTSS_cov_m)
colnames(AtTSS_cov_m) <- c("position","TSSset", "coverage")
p <- ggplot(AtTSS_cov_m, aes(x=position, y=coverage, by=TSSset, color=TSSset))
p + geom_line(size=1, alpha=0.4) + xlab("Position relative to TSS (TAIR10)") + ylab("Normalized coverage") + theme_bw(base_size = 16) + scale_colour_manual(values=c("navy","darkgreen","#E69F00"))
ggsave("AtTSS_all_raw_coverage_v_TSSs.png")
