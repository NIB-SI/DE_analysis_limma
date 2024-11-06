#analysis of SUSPHIRE RNAseq dataset with 2 groups (virgin females, mated females), 4 biol. repl.
#raw total counts exported from CLC Genomics Workbench 11.0.1

library ("limma")
library("edgeR")
library("stringr")


## Read counts table into object x
x1 <- read.table("./input/Expression Browser.txt", header=TRUE, sep="\t", row.names="Name", stringsAsFactors=FALSE, dec = ".")
colnames(x1)
rownames(x1)

x <- subset(x1, select=c(4,8,12,16,20,24,28,32))
colnames(x) <- str_sub(colnames(x), 9, 16)

group <- factor(c(1,1,1,1,2,2,2,2))

## Create a DGEList object for limma statistical analysis and specify library size i.e. number of reads sequenced per sample
y <- DGEList(counts=x, group=group, lib.size=c(41049930, 35952768, 42218346, 42697256, 41771538, 50015709, 54231914, 39663913))
y <- DGEList(counts=x, group=group)
#check density plot
pdf("other/density_plot_before_lowexpr_filter.pdf", onefile=TRUE, family= "Helvetica")
nsamples <- ncol(x)
col = c(rep("green", 4), rep("red", 4))

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()

## Filter genes based on expression ? keep genes (rows) that have over 50 counts in at least 3 samples
# filtering is necessary for voom method to work properly
keep.exprs <- rowSums(y$counts>50)>=3


## Normalization of dataset for different library sizes
y1 <- y[keep.exprs, , keep.lib.sizes=TRUE]
y1 <- calcNormFactors(y1)

## Plot QC plots using different functions e.g.:
col = c(rep("green", 4), rep("red", 4))

pdf("other/log10rawcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

pdf("other/log10filteredcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y1$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

#density plots before and after removing low expressed genes
pdf("other/norm_counts_raw&filtered_densityplots.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(x)
col = c(rep("green", 4), rep("red", 4))

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")

lcpm <- log(as.matrix(y1),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
title(main="B. AFTER REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")
par(opar)
dev.off()


#MDS (PCA-like) graph
pdf("other/norm_counts_filtered_MDS.pdf", onefile=TRUE, family= "Helvetica")
plotMDS(y1, labels=colnames(y1), col = col, cex = 0.6)
dev.off()

pdf("other/norm_cf_lcpm_MDS.pdf", onefile=TRUE, family= "Helvetica")
lcpm <- log(as.matrix(y1),10)
plotMDS(lcpm, labels=colnames(lcpm), col = col, cex = 0.6)
dev.off()

#PCA plots and genes that contribute most to PC1 and PC2
pc<- prcomp(t(as.matrix(y1)))

pdf("other/PCA_biplot.pdf", onefile=TRUE, family= "Helvetica")
biplot(pc, expand=1, cex=0.6)
dev.off()
#too many genes (variables) shown!!! ASK PROF. BLEJEC HOW TO RESTRICT THESE TO TOP10!

#top10 PC1 (control, infection)
PC1_top10 <- as.matrix(round(sort(pc$rotation[,1], decreasing = TRUE)[1:10],2))
colnames(PC1_top10) <- "PC1"

#top10 PC2 (June, Sept)
PC2_top10 <- as.matrix(round(sort(pc$rotation[,2], decreasing = TRUE)[1:10],2))
colnames(PC2_top10) <- "PC2"

write.table(PC1_top10, file = "output/PCA_PC1_top10genes.txt", sep = "\t")
write.table(PC2_top10, file = "output/PCA_PC2_top10genes.txt", sep = "\t")

# install.packages(c("FactoMineR", "factoextra"))
# library("FactoMineR")
# library("factoextra")
# 
# fviz_pca_biplot(pc, col.ind = col, palette = "jco", addEllipses = TRUE, label = "var",
#                 col.var = "black", repel = TRUE, legend.title = "PCA") 


## limma-voom protocol
# Create design matrix
design <- model.matrix(~0+group)
colnames(design) <- c("virgin", "mated")


# limma voom fit for filtered RNA-seq dataset (y1)
pdf("other/voom_mean-variance_trend.pdf", onefile=TRUE, family= "Helvetica")
v <- voom(y1,design,plot=TRUE)
dev.off()
fit <- lmFit(v, design)

## Define contrasts i.e. comparisons between groups
contrastMatrix = makeContrasts("virgin-mated", levels=design)
fit2 = contrasts.fit(fit, contrastMatrix)
## Check DEG in contrasts (adj p-val cutoff 0.05, |logFC| > 0)
tfit <- treat(fit2)
tfit
logFCcut <- 1
dt <- decideTests(tfit, lfc=logFCcut) #defaults adjust.method = "BH", p.value = 0.05, lfc=0
summary(dt)
colnames(dt)

## eBayes statistics calculation
fit2 <- eBayes(fit2)
pdf("other/SIGMA_vs_A_plot.pdf", onefile=TRUE, family= "Helvetica")
plotSA(fit2)
dev.off()

## make results table
results <- topTable(fit2, coef=1, number=1000000, sort.by="none")

# add raw expression data
length(rownames(y1))
length(results[,1])
summary(rownames(y1) == results[,1]) # all FALSE, have to do merge, not cbind

results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
write.table(results.raw, file="output/P_citri_RNAseq_logFC_padj_counts.txt", sep="\t", quote=TRUE, row.names=FALSE)


sessionInfo()
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Slovenian_Slovenia.1250  LC_CTYPE=Slovenian_Slovenia.1250    LC_MONETARY=Slovenian_Slovenia.1250
# [4] LC_NUMERIC=C                        LC_TIME=Slovenian_Slovenia.1250    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.3.1 edgeR_3.20.9  limma_3.34.9 
# 
# loaded via a namespace (and not attached):
#   [1] compiler_3.4.2  magrittr_1.5    tools_3.4.2     yaml_2.1.19     Rcpp_0.12.17    stringi_1.1.7   grid_3.4.2     
# [8] locfit_1.5-9.1  lattice_0.20-35