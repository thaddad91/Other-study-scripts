#!/usr/bin/env Rscript
#Script to run DESEQ2 on Transcriptome data

#For this step we need a count matrix called cts and a sample information called coldata
#The design variable corresponds to 

#dds <- DESeqDataSetFromMatrix(countData = cts,
#                              colData = coldata,
#                              design= ~ batch + condition)
#dds <- DESeq(dds)
#resultsNames(dds) # lists the coefficients
#res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

#These commands process the dds object, the factors are organized as vectors
#condition = factor(c("ST","ST","ST","HT","HT","HT"),c("ST","HT"))
#col_data = data.frame(condition)
#dds = DESeqDataSetFromMatrix(expression_data, col_data, ~condition)

#Tximport can be used to import the data and make a txi object that can be
#converted into the matrix that's needed
#library("tximport")
#library("readr")
#library("tximportData")
#dir <- system.file("extdata", package="tximportData")
#samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
#samples$condition <- factor(rep(c("A","B"),each=3))
#rownames(samples) <- samples$run

unzipped=unzip("/SRP095740/SRR5133628_fastq.gz")
expression_data <- read.table(unzipped, row.names=1,header=TRUE, sep =",", stringsAsFactors=FALSE)
dim(expression_data)
class(expression_data)
colnames(expression_data)
rownames(expression_data)
summary(expression_data)
apply(expression_data,2,sum)
mx = apply(expression_data,1,max)
expression_data = expression_data[ mx>10, ]
dim(expression_data)
condition = factor(c("ST","ST","ST","HT","HT","HT"),c("ST","HT"))
col_data = data.frame(condition)
dds = DESeqDataSetFromMatrix(expression_data, col_data, ~condition)
dds = estimateSizeFactors(dds)
sizeFactors(dds)
norm_versus_non_norm( dds, 1, 2, left = 2, right = 8 )

#From here on it's cluster analysis
rld = rlog(dds)
plot(density(assay(dds)[,1]), main="counts")
plot(density(assay(rld)[,1]), main="log counts")
dists = dist(t(assay(rld)))
plot(hclust(dists))
dds = estimateDispersions(dds)
plotDispEsts(dds)
dds = nbinomWaldTest(dds)
res = results(dds)
head(res)
res$padj = ifelse(is.na(res$padj), 1, res$padj)
write.table(res, col.names=NA, row.names=T, file ="expressions.tsv", sep ="\t")
plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)
View(res)
View(res)
View(res)
print(res)
head(res,10)
