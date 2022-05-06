rm(list = ls())
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dslabs")
install.packages("data.table")
install.packages("magrittr")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("rlang")
install.packages("ggrepel")
install.packages("rmarkdown")
install.packages("dabestr")
BiocManager::install("vsn")
BiocManager::install("recount3", force = TRUE)
BiocManager::install("ReportingTools")
BiocManager::install("DESeq2")
BiocManager::install("GSEABase")
BiocManager::install("apeglm")
BiocManager::install("clusterProfiler", version = "3.14", force = TRUE)
BiocManager::install("pathview", force = TRUE)
BiocManager::install("enrichplot", force = TRUE)
BiocManager::install("EnhancedVolcano")
library("recount3")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("ReportingTools")
library("pheatmap")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("tidyr")
library("magrittr")
library("stringr")
library("fgsea")
library("EnhancedVolcano")
require(reshape2)
library(clusterProfiler)
library(enrichplot)
require(DOSE)
library(dabestr)

#----------------------------

readcounts <- read.csv("C:\\Users\\alex\\gobi2021\\Block\\Project3\\project3_raw.readcounts.csv", sep = ",")
readcounts <- as.data.frame(readcounts)
dim(readcounts)
head(readcounts)

sampleinfo <- read.csv("C:\\Users\\alex\\gobi2021\\Block\\Project3\\project3.sample_info.csv", sep = ",")
sampleinfo <- as.data.frame(sampleinfo)


umicounts <- read.csv("C:\\Users\\alex\\gobi2021\\Block\\Project3\\project3_raw.umicounts.csv", sep = ",")
umicounts <- as.data.frame(umicounts)

geneNames <- read.csv("C:\\Users\\alex\\gobi2021\\Block\\Project3\\project3.gene_names.txt", sep = "\t")
geneNames <- as.data.frame(geneNames)

#merged readcounts
merged <- merge(readcounts, geneNames, by.x = "ï..ENSEMBL", by.y = "gene_id")
merged <- merged[,-1]
merged <- merged %>% relocate("gene_name")

x <- make.names(merged$gene_name,unique=T)
merged$gene_name <- x
sum(duplicated(merged$gene_name)) 


#merged umicounts
mergedumi <- merge(umicounts, geneNames, by.x = "ï..ENSEMBL", by.y = "gene_id")
mergedumi <- mergedumi[,-1]
mergedumi <- mergedumi %>% relocate("gene_name")

y <- make.names(mergedumi$gene_name, unique=T)
mergedumi$gene_name <-y
sum(duplicated(mergedumi$gene_name))

#DESeq2 readcounts
#----------------------------
norm <- DESeqDataSetFromMatrix(
  countData <- merged,
  colData <- sampleinfo[,c(2,5)],
  design = ~ Condition, tidy = TRUE)
class(norm)
is(norm, "SummarizedExperiment")

reg <- DESeq(norm)
head(reg)

#add cols to reg
colData(reg)

reg$type1 <- grepl("M1",reg$Condition)
reg$type2 <- grepl("WT",reg$Condition)
colData(reg)

reg$type1 <- as.character(reg$type1)
reg$type2 <- as.character(reg$type2)

reg$type1 <- str_replace_all(reg$type1, "TRUE", "M1")
reg$type1 <- str_replace_all(reg$type1, "FALSE", "M2")
reg$type2 <- str_replace_all(reg$type2, "TRUE", "WT")
reg$type2 <- str_replace_all(reg$type2, "FALSE", "KO")
colData(reg)


res <- results(reg)
summary(res)

res <- res[order(res$padj),]
head(res)


#split into M1 and M2
m1in <- readcounts[,c(1,3,4,6,11,12,13,14,17,20,21,22,23)]
m2in <- readcounts[,c(1,2,5,7,8,9,10,15,16,18,19,24,25)]

m1sample <- sampleinfo[grepl("M1", sampleinfo$Sample_ID), ]
m2sample <- sampleinfo[grepl("M2", sampleinfo$Sample_ID), ]

m1des <- DESeqDataSetFromMatrix(
  countData <- m1in,
  colData <- m1sample[,c(2,5)],
  design = ~ Condition, tidy = TRUE)

m2des <- DESeqDataSetFromMatrix(
  countData <- m2in,
  colData <- m2sample[,c(2,5)],
  design = ~ Condition, tidy = TRUE)

m1reg <- DESeq(m1des)
m2reg <- DESeq(m2des)

m1res <- results(m1reg)
m2res <- results(m2reg)


#volcano Plot

# add a column of NAs
res$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"
#res$diffexpressed

row.names(res)

res$reslabel <- NA
res$reslabel[res$diffexpressed != "NO"] <- row.names(res)[res$diffexpressed != "NO"]
res$reslabel
head(res)

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Volcano.png",
    width=1000, height=1000)
ggplot(as(res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) + geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept= 10e-6, col="grey") + 
  #geom_text() + 
  theme(text = element_text(size = 30))
dev.off()



png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Volcano.png",
    width=1000, height=1000)
EnhancedVolcano(res, lab=rownames(res), x = "log2FoldChange", y = "pvalue" )
dev.off()





# add a column of NAs
m1res$diffexpressed <- "NO"
m2res$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
m1res$diffexpressed[m1res$log2FoldChange > 2 & m1res$pvalue < 0.05] <- "UP"
m2res$diffexpressed[m2res$log2FoldChange > 2 & m2res$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
m1res$diffexpressed[m1res$log2FoldChange < -2 & m1res$pvalue < 0.05] <- "DOWN"
m2res$diffexpressed[m2res$log2FoldChange < -2 & m2res$pvalue < 0.05] <- "DOWN"

m1res$reslabel <- NA
m2res$reslabel <- NA
m1res$reslabel[m1res$diffexpressed != "NO"] <- row.names(m1res)[m1res$diffexpressed != "NO"]
m2res$reslabel[m2res$diffexpressed != "NO"] <- row.names(m2res)[m2res$diffexpressed != "NO"]



png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\VolcanoM1.png",
    width=1000, height=1000)
ggplot(as(m1res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) +      geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=10e-6, col="grey") +
  theme(text = element_text(size = 30))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\VolcanoM2.png",
    width=1000, height=1000)
ggplot(as(m2res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) +      geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=10e-6, col="grey") +
  theme(text = element_text(size = 30))
dev.off()

#Deseq umi
#----------------------------

umiDeseq <- DESeqDataSetFromMatrix(
  umicountData <- mergedumi,
  umicolData <- sampleinfo[,c(2,5)],
  design = ~ Condition, tidy = TRUE)
class(umiDeseq)
is(umiDeseq, "SummarizedExperiment")

umi <- DESeq(umiDeseq)
head(umi)

resUMI <- results(umi)
summary(resUMI)

resUMI <- resUMI[order(resUMI$padj),]
head(resUMI)

#volcano plot
# add a column of NAs
resUMI$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
resUMI$diffexpressed[resUMI$log2FoldChange > 2 & resUMI$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
resUMI$diffexpressed[resUMI$log2FoldChange < -2 & resUMI$pvalue < 0.05] <- "DOWN"

#row.names(resUMI)

resUMI$reslabelUMI <- NA
resUMI$reslabelUMI[resUMI$diffexpressed != "NO"] <- row.names(resUMI)[resUMI$diffexpressed != "NO"]
resUMI$reslabelUMI

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\VolcanoUMI.png",
    width=1000, height=1000)
ggplot(as(resUMI, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabelUMI)) + geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="grey") +
  theme(text = element_text(size = 30))
dev.off()

head(res)
head(resUMI)

#----------------------------

is.data.frame(readcounts)
colnames(readcounts)
readcounts[readcounts == 0] <- NA
readcounts
is.na(readcounts)
typeof(readcounts)

readcountsMelted <- melt(readcounts, id.vars = "ï..ENSEMBL", na.rm = TRUE)

colnames(readcountsMelted)


#TODO make faster by pre calculating?
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\readcountDistributions.png",
    width=2000, height=1000)
ggplot(readcountsMelted, aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))
#ggplot(readcounts, aes(y = readcounts$let7b_KO_M2_2)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)
dev.off()


png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\readcountDistributionsHist.png",
    width=2000, height=1000)
ggplot(as(readcountsMelted, "data.frame"), aes(x=value)) + geom_histogram(binwidth = 1000) + scale_x_continuous() + scale_y_continuous()
dev.off()


#----------------------------
#plot how many counts per sample total

#need 
#sum each col in readcounts
typeof(readcounts)
names <- colnames(readcounts[,-1])
countSums <- colSums(na.rm = TRUE, readcounts[,-1])
countSums

#split into WT or KO and label
categories <- c("Knockout", "Knockout", "Wild Type","Knockout", "Knockout", "Knockout", "Wild Type", "Knockout", "Wild Type","Wild Type", "Wild Type","Knockout", "Knockout","Wild Type", "Knockout", "Wild Type", "Knockout","Wild Type","Wild Type", "Knockout", "Knockout","Wild Type","Wild Type","Wild Type")
str(readcounts)
counts <- data.frame(categories,names, countSums)

WTcounts <- counts %>% filter(categories == "Wild Type")
KOcounts <- counts %>% filter(categories == "Knockout")

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\WTcountsPerSample.png",
    width=2000, height=1000)
ggplot(data = WTcounts, aes(x = names, y = countSums, fill = categories)) + geom_bar(stat = "identity") + theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + theme(text = element_text(size = 40), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\KOcountsPerSample.png",
    width=2000, height=1000)
ggplot(data = KOcounts, aes(x = names, y = countSums, fill = categories)) + geom_bar(stat = "identity") + theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + theme(text = element_text(size = 40), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#----------------------------
#UMI readcount distribution

umicounts[umicounts == 0] <- NA
colnames(umicounts[,-1])
uminames <- colnames(umicounts[,-1])
umicountSums <- colSums(na.rm = TRUE, umicounts[,-1])
#split into WT or KO and label
categories <- c("Knockout", "Knockout", "Wild Type","Knockout", "Knockout", "Knockout", "Wild Type", "Knockout", "Wild Type","Wild Type", "Wild Type","Knockout", "Knockout","Wild Type", "Knockout", "Wild Type", "Knockout","Wild Type","Wild Type", "Knockout", "Knockout","Wild Type","Wild Type","Wild Type")
str(umicounts)
countsUMI <- data.frame(categories, uminames, umicountSums)

umiWTcounts <- countsUMI %>% filter(categories == "Wild Type")
umiKOcounts <- countsUMI %>% filter(categories == "Knockout")

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\umiWTcountsPerSample.png",
    width=2000, height=1000)
ggplot(data = umiWTcounts, aes(x = uminames, y = umicountSums, fill = categories)) + geom_bar(stat = "identity") + theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + theme(text = element_text(size = 40), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\umiKOcountsPerSample.png",
    width=2000, height=1000)
ggplot(data = umiKOcounts, aes(x = uminames, y = umicountSums, fill = categories)) + geom_bar(stat = "identity") + theme_minimal() + scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + theme(text = element_text(size = 40), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#----------------------------
colData(reg)
resultsNames(reg)

reg$Condition <- factor(reg$Condition, levels = c("KO_M1", "KO_M2", "WT_M1", "WT_M2"))
reg$Condition
colData(reg)

res <- results(reg, contrast = c("Condition", "KO_M1", "KO_M2"))
res

resultsNames(reg)

resLFC <- lfcShrink(reg, coef = "Condition_KO_M2_vs_KO_M1", type = "apeglm")
resLFC
resNorm <- lfcShrink(reg, coef = "Condition_KO_M2_vs_KO_M1", type = "normal")

#MA lfc plot
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Ma.png",
    width=1000, height=500)
Ma <- plotMA(res, ylim = c( -2, 2))
dev.off()

#MA apeglm plot
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\MaApeglm.png",
    width=1000, height=500)
Ma <- plotMA(resLFC, ylim = c( -2, 2))
dev.off()

#MA normal 
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\MaNormal.png",
    width=1000, height=500)
Ma <- plotMA(resNorm, ylim = c( -2, 2))
dev.off()


#----------------------------


vsd <- vst(reg, blind=FALSE)
rld <- rlog(reg, blind=FALSE)
head(assay(vsd), 3)


# this gives log2(n + 1)
ntd <- normTransform(norm)
library("vsn")
meanSdPlot(assay(ntd))

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\meanSdNormal.png",
    width=1000, height=1000)
meanSdPlot(assay(ntd))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\meanSdVSD.png",
    width=1000, height=1000)
meanSdPlot(assay(vsd))
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\meanSdRLD.png",
    width=1000, height=1000)
meanSdPlot(assay(rld))
dev.off()

#----------------------------
#basemean
colnames(res)

base <- data.frame(res[,c(1,7)])
head(base)

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Basemean.png",
    width=1000, height=1000)
ggplot(base, aes(x = baseMean, color = diffexpressed)) + geom_density(kernel = "gaussian")
dev.off()

#----------------------------
#PCA 

KOcounts <- count_rlog[,c(1,2,4,5,6,8,12,13,15,17,20,21)]
colnames(KOcounts)
WTcounts <- count_rlog[,c(3,7,9,10,11,14,16,18,19,22,23,24)]
colnames(WTcounts)

colnames(res)
head(res)

count_rlog = rlogTransformation(reg)
dim(count_rlog)
colnames(count_rlog)
row.names(count_rlog)
colData(count_rlog)

png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\PCAreadcounts.png",
    width=1000, height=500)
pcaplot <- plotPCA(count_rlog, intgroup = "Sample_ID", returnData=TRUE)
percentVar <- round(100 * attr(pcaplot, "percentVar"))
ggplot(pcaplot, aes(PC1, PC2, color=count_rlog$type1, shape=count_rlog$type2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

colData(count_rlog)


#----------------------------
#Heatmap

#TODO change colors of groups

#NORM
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\heatmapNorm.png",
    width=1000, height=1000)
select <- order(rowMeans(counts(reg,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(reg)[,c("type1","type2")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 20)
dev.off()

#VSD
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\heatmapVSD.png",
    width=1000, height=1000)
select <- order(rowMeans(counts(reg,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(reg)[,c("type1","type2")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 20)
dev.off()

#RLD
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\heatmapRLD.png",
    width=1000, height=1000)
select <- order(rowMeans(counts(reg,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(reg)[,c("type1","type2")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, fontsize = 20)
dev.off()


#----------------------------
# Sample distances

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type1, vsd$type2, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$type1, vsd$type2, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\sample2sampleheatmap.png",
    width=1000, height=1000)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, fontsize = 20)
dev.off()

#----------------------------
#----------------------------
#Cluster profiling Data prep

readcounts$ï..ENSEMBL <- sub("\\..*$", "",readcounts$ï..ENSEMBL)
readcounts$ï..ENSEMBL

readcountids <- DESeqDataSetFromMatrix(
  countData <- readcounts,
  colData <- sampleinfo[,c(2,5)],
  design = ~ Condition, tidy = TRUE)

readids <- DESeq(readcountids)

gseares <- results(readids)

# add a column of NAs
gseares$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gseares$diffexpressed[gseares$log2FoldChange > 2 & gseares$pvalue < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gseares$diffexpressed[gseares$log2FoldChange < -2 & gseares$pvalue < 0.05] <- "DOWN"

gseares$reslabel <- NA
gseares$reslabel[gseares$diffexpressed != "NO"] <- row.names(gseares)[gseares$diffexpressed != "NO"]
dim(gseares)

gseares$diffexpressed <- na_if(gseares$diffexpressed, "NO")
#gseares$diffexpressed
filtered <- subset(gseares, !is.na(diffexpressed))
#is.na(filtered$diffexpressed)
dim (filtered)

#----------------------------
#Cluster profiling

# Mouse 
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#need log2fold as svector
original_gene_list <-filtered$log2FoldChange

#name 
names(original_gene_list) <- row.names(filtered)

#omit NAs
gene_list <- na.omit(original_gene_list)

#sort in decreasing order !!!
gene_list = sort(gene_list, decreasing = TRUE)
gene_list

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 1000,
             minGSSize = 5, 
             maxGSSize = 100, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")

#TODO save these plots


dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
emapplot(gse, showCategory = 10)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
ridgeplot(gse) + labs(x = "enrichment distribution")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)



#----------------------------
#----------------------------
#----------------------------
#----------------------------
#----------------------------
#TODO dabest

plot_dabest <- function(input_df_dabest) {
  input_df_dabest <- na.omit(input_df_dabest) # remove NA
  
  # Performing unpaired (two independent groups) analysis.
  unpaired_mean_diff <- dabest(input_df_dabest, macrophage, log2FoldChange,
                               idx = c("M1", "M2"),
                               paired = FALSE)
  
  # Display the results in a user-friendly format.
  #unpaired_mean_diff
  # Compute the mean difference.
  #mean_diff(unpaired_mean_diff)
  
  # Plotting the mean differences.
  mean_diff(unpaired_mean_diff) %>% plot()
  
}


## 1) select all significant genes from the entire data set
m1res
m2res
#select only genes that were classified as significant earlier
m1res$diffexpressed <- na_if(m1res$diffexpressed, "NO")
m2res$diffexpressed <- na_if(m2res$diffexpressed, "NO")
filteredm1 <- subset(m1res, !is.na(diffexpressed))
filteredm2 <- subset(m2res, !is.na(diffexpressed))
dim (filteredm1)
dim (filteredm2)
head(filteredm1)

da_best_M1 <- data.frame(filteredm1$reslabel, filteredm1$log2FoldChange)
da_best_M1$macrophage <- "M1"
colnames(da_best_M1) <- c("gene", "log2FoldChange", "macrophage")
da_best_M2 <- data.frame(filteredm2$reslabel, filteredm2$log2FoldChange)
da_best_M2$macrophage <- "M2"
colnames(da_best_M2) <- c("gene", "log2FoldChange", "macrophage")

# save data from M1 and M2 in one dataframe
dabest_M1M2 <- merge(da_best_M1, da_best_M2, by=c("gene", "log2FoldChange", "macrophage"), all.x = T, all.y = T)

plot_dabest(dabest_M1M2)


##----------------------------------------------
## 2) select all significant genes from M1 
##	-> get logfoldchange of these genes in M1 and M2 separately and plot them
## get M1 lfc
norm <- DESeqDataSetFromMatrix(
  countData <- only_m1,
  colData <- sampleinfo_M1[,c(2,5)],
  design = ~ Condition, tidy = TRUE)


reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]

res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)

# fc_M1 contains all up- or down-regulated genes from M1
fc_M1 <- data.frame(res[res$diffexpressed == 'DOWN' | res$diffexpressed == 'UP' & res$reslabel != 'NA', ]$reslabel,
                    res[res$diffexpressed == 'DOWN' | res$diffexpressed == 'UP' & res$reslabel != 'NA', ]$log2FoldChange)

colnames(fc_M1) <- c("gene", "log2FoldChange")
fc_M1$macrophage <- "M1"


## get M2 lfc
norm <- DESeqDataSetFromMatrix(
  countData <- only_m2,
  colData <- sampleinfo_M2[,c(2,5)],
  design = ~ Condition, tidy = TRUE)


reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]

res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)


#save only genes that are significant in M1 ("focusM1")
bla <- res[rownames(res) %in% fc_M1$gene, ]
fc_M2 <- data.frame(bla$reslabel, bla$log2FoldChange)
fc_M2$macrophage <- "M2"
colnames(fc_M2) <- c("gene", "log2FoldChange", "macrophage")


# dabest_focusM1 contains all the log2foldchange for M1 and M2 separately
# (for all genes that were significant in M1, size = 248*2)
dabest_focusM1 <- merge(fc_M1, fc_M2, by=c("gene", "log2FoldChange", "macrophage"), all.x = T, all.y = T)


plot_dabest(dabest_focusM1)


##----------------------------------------------
## 3) select all significant genes from M2
##	-> get logfoldchange of these genes in M1 and M2 separately and plot them
## get M2 lfc
norm <- DESeqDataSetFromMatrix(
  countData <- only_m2,
  colData <- sampleinfo_M2[,c(2,5)],
  design = ~ Condition, tidy = TRUE)


reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]

res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)

#save all genes that are significant in M2 (for focusM2) -> 245
fc_M2 <- data.frame(res[res$diffexpressed == 'DOWN' | res$diffexpressed == 'UP' & res$reslabel != 'NA', ]$reslabel,
                    res[res$diffexpressed == 'DOWN' | res$diffexpressed == 'UP' & res$reslabel != 'NA', ]$log2FoldChange)
fc_M2$macrophage <- "M2"
colnames(fc_M2) <- c("gene", "log2FoldChange", "macrophage")


## get M1 lfc
norm <- DESeqDataSetFromMatrix(
  countData <- only_m1,
  colData <- sampleinfo_M1[,c(2,5)],
  design = ~ Condition, tidy = TRUE)


reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]

res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)


#for focusM2: get lfc for all genes that are significant in M2 (saved in fc_M2)
bla <- res[rownames(res) %in% fc_M2$gene, ]
fc_M1 <- data.frame(bla$reslabel, bla$log2FoldChange)
fc_M1$macrophage <- "M1"
colnames(fc_M1) <- c("gene", "log2FoldChange", "macrophage") 

# dabest_focusM2 contains all the log2foldchange for M1 and M2 separately
# (for all genes that were significant in M2, size = 245*2)
dabest_focusM2 <- merge(fc_M1, fc_M2, by=c("gene", "log2FoldChange", "macrophage"), all.x = T, all.y = T)


plot_dabest(dabest_focusM2)




#----------------------------




#----------------------------



#----------------------------
