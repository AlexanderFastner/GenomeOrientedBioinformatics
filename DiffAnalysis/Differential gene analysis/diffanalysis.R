#differential gene analysis
#-------------------------------------------------------------------------------
rm(list = ls())
install.packages("dslabs")
install.packages("data.table")
install.packages("magrittr")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReportingTools")
BiocManager::install("DESeq2")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("ReportingTools")
#-------------------------------------------------------------------------------
#read in 
hisat <- as.matrix(read.csv("C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\gene.counts.hisat", sep = "\t", row.names = "Geneid"))
star <- as.matrix(read.csv("C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\gene.counts.star", sep = "\t", row.names = "Geneid"))
listannot <- readr::read_delim("C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\sample.list", delim="\t")

#check dim and output
dim(hisat)
head(hisat)
dim(star)
head(star)
head(listannot)

#-------------------------------------------------------------------------------
#change listannot into S1,S2,mock
listannot$condition[grep("S2+", listannot$condition)] <- "S2"
listannot$condition[grep("S1+", listannot$condition)] <- "S1"
listannot$condition[grep("mock+", listannot$condition)] <- "mock"

listannot <- mutate(listannot, condition = factor(condition, levels = c("S1", "S2", "mock")))
#listannot <- mutate(listannot, condition = factor(condition, levels = c("S1", "S2", "mock")),
#                    sample = factor(sample))
#istannot <- mutate(listannot, condition = factor(condition),
                   #sample = factor(sample))
with(listannot,
     table(condition, sample))

#-------------------------------------------------------------------------------
mt <- match(colnames(hisat), sub(".bam", "", colnames(star)))
#check na
stopifnot(!any(is.na(mt)))


cov <- DESeqDataSetFromMatrix(
  countData <- hisat,
  colData <- listannot,
  design = ~ condition)
class(cov)
is(cov, "SummarizedExperiment")
#check na
stopifnot(!any(is.na(cov)))


cov2 <- DESeqDataSetFromMatrix(
  countData <- star,
  colData <- listannot,
  design = ~ condition)
class(cov2)
is(cov2, "SummarizedExperiment")

cov2 <- DESeq(cov2)
res2 <- results(cov2)


cov <- DESeq(cov)
res <- results(cov)
#remove NA from DESeq object?

res[order(res$padj), ] | head()

#-------------------------------------------------------------------------------
#size factors
sizefact <- plot(estimateSizeFactorsForMatrix(hisat), estimateSizeFactorsForMatrix(star), lty = 1, type = "b", pch = 1, col = "red")
sizefact <- lines((colSums(hisat)/10000000), (colSums(star)/10000000), lty = 1, type = "b", pch = 0, col = "blue")

#-------------------------------------------------------------------------------
#p-value histogram
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\phist.png",
    width=600, height=350)
ggplot(as(res, "data.frame"), aes(x = pvalue)) + geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\phist2.png",
    width=600, height=350)
ggplot(as(res2, "data.frame"), aes(x = pvalue)) + geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)
dev.off()

#-------------------------------------------------------------------------------
#MA plot
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\Ma.png",
    width=600, height=350)
Ma <- plotMA(cov, ylim = c( -2, 2))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\Ma2.png",
    width=600, height=350)
Ma <- plotMA(cov2, ylim = c( -2, 2))
dev.off()

#-------------------------------------------------------------------------------
#PCA
cov_rlog = rlogTransformation(cov)
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\pcasimple.png",
    width=600, height=250)
plotPCA(cov_rlog, intgroup = "condition") + coord_fixed()
dev.off()
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\pca.png",
    width=600, height=250)
plotPCA(cov_rlog, intgroup = c("condition", "sample")) + coord_fixed()
dev.off()


cov_rlog2 = rlogTransformation(cov2)
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\pcasimple2.png",
    width=600, height=250)
plotPCA(cov_rlog2, intgroup = "condition") + coord_fixed()
dev.off()
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\pca2.png",
    width=600, height=250)
plotPCA(cov_rlog2, intgroup = c("condition", "sample")) + coord_fixed()
dev.off()

#-------------------------------------------------------------------------------
#heatmap
library("pheatmap")
png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\heatmap.png",
    width=600, height=600)
select <- order(rowMeans(assay(cov_rlog)), decreasing = TRUE)[1:30]
pheatmap( assay(cov_rlog)[select, ],
          scale = "row",
          annotation_col = as.data.frame(
            colData(cov_rlog)[, c("condition", "sample")] ))
dev.off()


png(file="C:\\Users\\alex\\gobi2021\\DiffAnalysis\\Differential gene analysis\\heatmap2.png",
    width=600, height=600)
select <- order(rowMeans(assay(cov_rlog2)), decreasing = TRUE)[1:30]
pheatmap( assay(cov_rlog2)[select, ],
          scale = "row",
          annotation_col = as.data.frame(
            colData(cov_rlog2)[, c("condition", "sample")] ))
dev.off()

#-------------------------------------------------------------------------------
#output results
#write.csv(as.data.frame(res), file = "cov-cov2.csv")

htmlRep <- HTMLReport(shortName = "cov-cov2.csv", title = "testing",
                       reportDirectory = "./reports")
publish("size factor estimation", htmlRep)
publish(sizefact, htmlRep, name = "Size estimate")
publish("p-value histogram", htmlRep)
publish(phist, htmlRep, name = "p-value histogram")
publish("MA plot", htmlRep)
publish(Ma, htmlRep, name = "MA plot")
publish("PCA with samples", htmlRep)
publish(pca2, htmlRep, name = "PCA with samples")
publish("PCA simplified", htmlRep)
publish(pca, htmlRep, name = "PCA simplified")

finish(htmlRep)








