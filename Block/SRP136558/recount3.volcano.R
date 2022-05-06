## recount3 SRP136558
# RNAseq to determine gene expression changes in adipose tissue macrophages (ATM1 and ATM2) and adipocyte progenitors (AP) following cold exposure or FGF21-mimetic antibody administration


readcounts_recount <- read.csv("C:\\Users\\kolle\\OneDrive - tum.de\\GOBI\\recount3_dataset\\51943c1e241e_sra.gene_sums.SRP136558_counts.M023", sep = "\t", skip = 2)
readcounts_recount <- as.data.frame(readcounts)



sampleinfo_recount <- read.csv("C:\\Users\\kolle\\OneDrive - tum.de\\GOBI\\recount3_dataset\\recount3.sampleInfo.txt", header = T, sep = "\t")
#remove some samples ("AP" and "BFK")
sampleinfo_recount <- sampleinfo_recount[sampleinfo_recount$condition != "AP" & sampleinfo_recount$condition != "BFK", ]
readcounts_recount <- readcounts_recount[, ! names(readcounts_recount) %in% sampleinfo_recount$external_id]


#----------------------------
# DESeq2
norm <- DESeqDataSetFromMatrix(
  countData <- readcounts_recount,
  colData <- sampleinfo_recount[,c(1,3)],
  design = ~ condition, tidy = TRUE)

reg <- DESeq(norm)
res <- results(reg, contrast = c("condition", "Cold_M1", "Cold_M2"))


res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)


ggplot(as(res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) + geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept= 10e-6, col="grey") + 
  #geom_text() + 
  theme(text = element_text(size = 30))

