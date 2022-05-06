
sampleinfo_M1 <- read.csv("C:\\Users\\kolle\\OneDrive - tum.de\\GOBI\\example_input\\project3.sample_info_M1.csv", sep = ",")
sampleinfo_M1 <- as.data.frame(sampleinfo_M1)

sampleinfo_M2 <- read.csv("C:\\Users\\kolle\\OneDrive - tum.de\\GOBI\\example_input\\project3.sample_info_M2.csv", sep = ",")
sampleinfo_M2 <- as.data.frame(sampleinfo_M2)




##todo: only M1 bzw M2 data frames
only_m1 <- data.frame(merged$gene_name, 
                      merged$let7b_KO_M1_1, merged$let7b_KO_M1_2, merged$let7b_KO_M1_3, merged$let7b_KO_M1_4, merged$let7b_KO_M1_5, merged$let7b_KO_M1_6,
                      merged$let7b_WT_M1_1, merged$let7b_WT_M1_2, merged$let7b_WT_M1_3, merged$let7b_WT_M1_4, merged$let7b_WT_M1_5, merged$let7b_WT_M1_6)
row.names(only_m1) <- only_m1$merged.gene_name

head(only_m1)

##todo: only M1 bzw M2
only_m2 <- data.frame(merged$gene_name, 
                      merged$let7b_KO_M2_1, merged$let7b_KO_M2_2, merged$let7b_KO_M2_3, merged$let7b_KO_M2_4, merged$let7b_KO_M2_5, merged$let7b_KO_M2_6,
                      merged$let7b_WT_M2_1, merged$let7b_WT_M2_2, merged$let7b_WT_M2_3, merged$let7b_WT_M2_4, merged$let7b_WT_M2_5, merged$let7b_WT_M2_6)
row.names(only_m2) <- only_m2$merged.gene_name



## ----------------------------------- DESeq M1
norm <- DESeqDataSetFromMatrix(
  countData <- only_m1,
  colData <- sampleinfo_M1[,c(2,5)],
  design = ~ Condition, tidy = TRUE)
reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]



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


#png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Volcano.png",
 #   width=1000, height=1000)
ggplot(as(res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) + geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept= 10e-6, col="grey") + 
  #geom_text() + 
  theme(text = element_text(size = 30))
dev.off()



## ----------------------------------- DESeq M2
norm <- DESeqDataSetFromMatrix(
  countData <- only_m2,
  colData <- sampleinfo_M2[,c(2,5)],
  design = ~ Condition, tidy = TRUE)


reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]



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


#png(file="C:\\Users\\alex\\gobi2021\\Block\\Project3\\Volcano.png",
#   width=1000, height=1000)
ggplot(as(res, "data.frame"), aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label = reslabel)) + geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept= 10e-6, col="grey") + 
  #geom_text() + 
  theme(text = element_text(size = 30))
dev.off()



