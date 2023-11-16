

rm(list = ls())
install.packages("ggplot2")
install.packages("dplyr")
library("ggplot2")
library("dplyr")
library("stringr")
#-------------------------------------------------------------------------------

#1
#No of genes, no of gene sets, no of leafs, shortest and longest path to root, number of ''shortcuts" in the DAG

#provide an overall visualization and stats for the DAGnode tree
#this can be a table


#-------------------------------------------------------------------------------
#2 Distribution of number of genes in the gene sets 

#read in 
allnodes <- as.data.frame(read.csv("C:\\Users\\alex\\gobi2021\\goenrich\\nomaxmin", sep = "\t", row.names = "term"))
dim(allnodes)
head(allnodes)
adjusted  <- allnodes %>% filter(size > 0, size < 500)

#go through all dagnodes and count num genes. then make a distribution from that
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\allgenesDistribution.png",
    width=800, height=500)
#ggplot(as(dagnodes, "data.frame"), aes(x = size)) + geom_area(stat = "bin")
ggplot(as(adjusted, "data.frame"), aes(x = size)) + geom_histogram(binwidth = 10, color = "blue", fill=I("blue"))+
    labs(x = "Size of Geneset", y = "Count", title = "Distribution of Geneset sizes for all Genesets")
dev.off()
#-------------------------------------------------------------------------------
#3 Same, for gene sets between minsize and maxsize

#dagnodes %>% filter(size > 50, size < 500)
#filtered  <- dagnodes %>% filter(size > 50, size < 500)
filtered <- as.data.frame(read.csv("C:\\Users\\alex\\gobi2021\\goenrich\\testing", sep = "\t", row.names = "term"))

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\genesBetweenMinandMax.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size)) + geom_histogram(binwidth = 10, color = "blue", fill=I("blue"))+
    labs(x = "Size of Geneset", y = "Count", title = "Distribution of Geneset sizes with size > 50 && size < 500")
dev.off()
#-------------------------------------------------------------------------------
#4 of path lengths from all leafs to the root

#get length of each path and plot those

pathlength <- str_count(filtered$shortest_path_to_a_true, "|")
df <- data.frame(val = pathlength)

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\PathLengthDistributions.png",
    width=800, height=500)
#todo add pathlength calculation
ggplot(as(df, "data.frame"), aes(x = pathlength)) + geom_histogram(binwidth = 10, color = "black", fill=I("blue")) + 
    labs(x = "Pathlength", y = "Count", title = "Distribution of path lengths to Root by geneset size")
dev.off()

#-------------------------------------------------------------------------------
#5 size of set differences between child and parent sets
#?? no idea


#-------------------------------------------------------------------------------
#6 Distribution of significant gens (SGs) in the all gene sets

#how many significant genes for each dag node
#noverlap / total genes in Dagnode


png(file="C:\\Users\\alex\\gobi2021\\goenrich\\distofSigGenes.png",
    width=800, height=500)
ggplot(as(adjusted, "data.frame"), aes(x = size, y = noverlap)) + geom_point(color = "blue", fill=I("black")) +
    geom_abline(intercept=0, slope=1) +
    labs(x = "Geneset Size", y = "Number of Significant Genes", title = "Distribution of significant genes over Genesets for all genes") 
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\distofSigGenesMinMax.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = noverlap)) + geom_point(color = "blue", fill=I("black")) +
    geom_abline(intercept=0, slope=1) +
    labs(x = "Geneset Size", y = "Number of Significant Genes", title = "Distribution of significant genes over Genesets for genes 50<size<500") 
dev.off()

#-------------------------------------------------------------------------------
#7 Same, for gene sets between minsize and maxsize

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\distofSigGenesHistminmax.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = noverlap)) + geom_col(binwidth = 10, color = "blue", fill=I("blue")) +
    labs(x = "Geneset Size", y = "Number of Significant Genes", title = "Distribution of significant genes over Genesets for genes 50<size<500")
dev.off()

#-------------------------------------------------------------------------------
#8  the enrichment scores and p-values

overlaps <- as.data.frame(read.csv("C:\\Users\\alex\\gobi2021\\goenrich\\overlap_out_tsv", sep = "\t"))
dim(overlaps)
head(overlaps)

#todo add cutoff at p = 0.05

#hg pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\hypergeometricPvalues.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = hg_pval)) + geom_area(stat = "bin", color = "black", fill=I("blue")) + 
    labs(x = "hypergeometric P-values", y = "Count", title = "Hypergeometric P-value counts")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\hypergeometricPvaluescorrected.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = hg_fdr)) + geom_area(stat = "bin", color = "black", fill=I("blue")) +
    labs(x = "hypergeometric P-values fdr adjusted size", y = "Count", title = "Hypergeometric P-value counts 50<size<500")
dev.off()

#Fishers exact pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\FishersExactPval.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = fej_pval)) + geom_area(stat = "bin", color = "black", fill=I("blue")) +
    labs(x = "Fishers Exact P-values", y = "Count", title = "Fishers Exact P-value counts")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\FishersExactPvalcorrected.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = fej_fdr)) + geom_area(stat = "bin", color = "black", fill=I("blue")) +
    labs(x = "Fishers Exact P-values fdr adjusted size", y = "Count", title = "Fishers Exact P-values counts 50<size<500")
dev.off()

#ks pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\KolmogorovSmirnovPval.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = ks_pval)) + geom_area(stat = "bin", color = "black", fill=I("blue")) +
    labs(x = "KolmogorovSmirnov P-values", y = "Count", title = "KolmogorovSmirnov P-values counts")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\KolmogorovSmirnovPvalcorrected.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = ks_fdr)) + geom_area(stat = "bin", color = "black", fill=I("blue")) +
    labs(x = "KolmogorovSmirnov P-values fdr adjusted size", y = "Count", title = "KolmogorovSmirnov P-values counts 50<size<500")
dev.off()

#-------------------------------------------------------------------------------
#9 scatter score vs size




#-------------------------------------------------------------------------------
#10 Scatter of p-value against size

#hg pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\hypergeometricPvaluesScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = hg_pval)) + geom_point(color = "black")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\hypergeometricPvaluescorrectedScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = hg_fdr)) + geom_point( color = "black")
dev.off()

#Fishers exact pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\FishersExactPvalScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = fej_pval)) + geom_point(color = "black")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\FishersExactPvalcorrectedScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = fej_fdr)) + geom_point(color = "black")
dev.off()

#ks pval
png(file="C:\\Users\\alex\\gobi2021\\goenrich\\KolmogorovSmirnovPvalScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = ks_pval)) + geom_point(color = "black")
dev.off()

png(file="C:\\Users\\alex\\gobi2021\\goenrich\\KolmogorovSmirnovPvalcorrectedScatter.png",
    width=800, height=500)
ggplot(as(filtered, "data.frame"), aes(x = size, y = ks_fdr)) + geom_point(color = "black")
dev.off()


#-------------------------------------------------------------------------------
#make table with genesets, number of sig, pvalues

newtable <- filtered[, c("size", "noverlap", "hg_fdr", "fej_fdr", "ks_fdr")]
dim(newtable)
head(newtable)

write.table(newtable,"C:\\Users\\alex\\gobi2021\\goenrich\\newReduced.tsv", quote=FALSE, sep='\t', col.names = NA)



#-------------------------------------------------------------------------------









