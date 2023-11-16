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
BiocManager::install("recount3", force = TRUE)
BiocManager::install("ReportingTools")
BiocManager::install("DESeq2")
library("recount3")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("ReportingTools")
#----------------------------


gene <- recount3::create_rse_manual(
  project = "SRP119232",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)

dim(gene)
head(gene)
rowData(gene)
colData(gene)


#----------------------------

exon <- recount3::create_rse_manual(
  project = "SRP119232",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "exon"
)

dim(exon)
head(exon)
rowData(exon)
colData(exon)

#----------------------------
# gene data plots
#head(gene)
#rowData(gene)
#colData(gene)
#assays(gene)
#colnames(bp_length)
#colData(bp_length)

genes <- data.frame(rowRanges(gene)) 


png(file="C:\\Users\\alex\\gobi2021\\Block\\geneWidth.png",
    width=600, height=350)
#ggplot(BPlength, aes(x = bp_length)) + geom_histogram(binwidth = 100, fill = "Royalblue", boundary = 0) + expand_limits(x = 0, count = 0) + scale_y_sqrt()
ggplot(genes, aes(x = width)) + geom_area(stat = "bin", fill = "Royalblue") + expand_limits(x = 0, count = 0) +# + scale_y_sqrt()
  labs(title="length of genes", x="length in Bp", y="Count") + 
dev.off()

#dotplot
png(file="C:\\Users\\alex\\gobi2021\\Block\\geneWidthdotplot.png",
    width=1500, height=500)
#ggplot(BPlength, aes(x = bp_length)) + geom_histogram(binwidth = 100, fill = "Royalblue", boundary = 0) + expand_limits(x = 0, count = 0) + scale_y_sqrt()
ggplot(genes, aes(x = seqnames, y = width)) + geom_dotplot(binaxis = "y" ,binwidth = 1/500, fill = "Royalblue") + # + expand_limits(x = 0, count = 0)# + scale_y_sqrt()
  labs(title="length of genes", x="length in Bp", y="Count") + 
dev.off()


#----------------------------
#exons data

exons <- data.frame(rowRanges(exon))
dim(exons)
head(exons)


png(file="C:\\Users\\alex\\gobi2021\\Block\\exonWidth.png",
    width=600, height=350)
ggplot(exons, aes(x = width)) + geom_histogram(binwidth = 100, fill = "Royalblue") + 
  labs(title="length of exons", x="length in Bp", y="Count") + 
  scale_x_continuous() + 
  scale_y_continuous()
dev.off()


#----------------------------
#gene type - width






#----------------------------




