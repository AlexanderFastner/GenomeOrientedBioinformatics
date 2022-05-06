## dabest
install.packages("dabestr")
library(dabestr)

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
##	-> get logfoldchange of these genes in M1 and M2 separately and plot them


norm <- DESeqDataSetFromMatrix(
  countData <- merged,
  colData <- sampleinfo[,c(2,5)],
  design = ~ Condition, tidy = TRUE)
class(norm)
is(norm, "SummarizedExperiment")

reg <- DESeq(norm)
res <- results(reg)
res <- res[order(res$padj),]

# significant genes (up- or down-regulated)
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 2 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -2 & res$pvalue < 0.05] <- "DOWN"

res$reslabel <- row.names(res)

# save only significant genes
signif_genes <- data.frame(res[res$diffexpressed == 'DOWN' | res$diffexpressed == 'UP' & res$reslabel != 'NA', ]$reslabel)
colnames(signif_genes) <- c("gene")

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

#select only genes that were classified as significant earlier
signif_M1 <- res[rownames(res) %in% signif_genes$gene, ]

da_best_M1 <- data.frame(signif_M1$reslabel, signif_M1$log2FoldChange)
da_best_M1$macrophage <- "M1"
colnames(da_best_M1) <- c("gene", "log2FoldChange", "macrophage")


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

#select only genes that were classified as significant earlier
signif_M2 <- res[rownames(res) %in% signif_genes$gene, ]

da_best_M2 <- data.frame(signif_M2$reslabel, signif_M2$log2FoldChange)
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





