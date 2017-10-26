
library(ggplot2)

library(reshape2)

setwd("/scratch/rtraborn/tsrchitect-figures/figures/")

coverage.out <- read.table(file="/scratch/rtraborn/tsrchitect-figures/output_files/Zm_coverage_hist_all_v_tss_annot.txt", skip=1, header=FALSE)

head(coverage.out)

colnames(coverage.out) <- c("Distance","ESTcoverage","EST_plusTags","EST_minusTags","CAGE1coverage","CAGE1_plusTags","CAGE1_minusTags","CAGE2coverage","CAGE2_plusTags","CAGE2_minusTags")

head(coverage.out)

coverage.table <- cbind(coverage.out[,1:2],coverage.out[,5],coverage.out[,8])

colnames(coverage.table) <- c("Distance","ESTcoverage","CAGEcoverage-rep1", "CAGEcoverage-rep2")

head(coverage.table)

coverage.rs <- melt(coverage.table,id.vars = "Distance")

colnames(coverage.rs) <- c("Distance","dataset", "density")

head(coverage.rs)

a <- ggplot(coverage.rs)

a + geom_line(aes(x=Distance, y=density, colour=dataset))

ggsave(file="Zm_peaks_v_annotated_genes.png")
