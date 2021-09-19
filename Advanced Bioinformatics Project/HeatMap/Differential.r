# get  the  table of read  counts
read.counts=read.table(
  "MIAgenes.txt", row.names=1, header=TRUE, sep ="\t", stringsAsFactors=FALSE) 

# the  gene  names are stored  as rownames
rownames(read.counts)

orig_names <- names(read.counts)

# make a data.frame where  rownames should  match  the sample  names

sample_info  <- data.frame(condition = gsub( "_[0 -9]+", "", names(read.counts)
),

row.names = names(read.counts) )

#install DESeq2
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
source("http://www.bioinformatics.nl/courses/RNAseq/DEseq2Exercise.R")

#generate  the  DESeqDataSet
mx = apply( read.counts, 1, max )
read.counts = read.counts[ mx > 5, ] 
condition = factor(c("CK","CK","CK","ET","ET","ET","MJ","MJ","MJ"),c("CK","ET","MJ"))
col_data = data.frame(condition)
DESeq.ds = DESeqDataSetFromMatrix(read.counts, col_data, ~condition)

# test  what  counts ()  returns
counts(DESeq.ds)

# remove  genes  without  any  counts
DESeq.ds <- DESeq.ds[ rowSums(counts(DESeq.ds)) > 0, ]

#different library  sizes
colSums(counts(DESeq.ds))

# calculate  the  size  factor  and add it to the  data  set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# counts gets the normalized counts
counts.sf_normalized  <- counts(DESeq.ds, normalized = TRUE)

# transform  size factor normalized  read  counts  to log2  scale
log.norm.counts  <- log2(counts.sf_normalized + 1)

#box plots  of non -transformed  read  counts (one  per  sample)
boxplot(counts.sf_normalized , notch = TRUE ,
        
        main = "untransformed  read  counts", ylab = "read  counts")

#box plots of log2 transformed read  counts
boxplot(log.norm.counts , notch = TRUE ,
        
        main = "log2 -transformed  read  counts",
        
        ylab = "log2(read  counts)")
plot(log.norm.counts [,1:2], cex=.1, main = "Normalized log2(read counts)")

#install  the vsn  package
source("http://bioconductor.org/biocLite.R")
biocLite("vsn")
library(vsn)
library("ggplot2")

#mean sd plot
msd_plot  <- meanSdPlot(log.norm.counts ,
                        
                        ranks=FALSE , # show  the  data on the  original  scale
                        
                        plot = FALSE)

msd_plot$gg +
  
  ggtitle("sequencing  depth  normalized  log2(read  counts)") +
  
  ylab("standard  deviation")

#obtain  regularized  log transformed values
DESeq.rlog  <- rlog(DESeq.ds , blind = TRUE)
rlog.norm.counts  <- assay(DESeq.rlog)

# mean sd plot for rlog transformed data
msd_plot  <- meanSdPlot(rlog.norm.counts ,
                        
                        ranks=FALSE , # show  the  data on the  original  scale
                        
                        plot = FALSE)
msd_plot$gg +
  
  ggtitle("rlog -transformed  read  counts") +
  
  ylab("standard  deviation")

#Plotting the counts and log counts
plot(density(assay(DESeq.ds)[,1]), main="counts")
plot(density(assay(DESeq.rlog)[,1]), main="log counts")

#Calculating the Euclidean distance and hierarchical clustering
dists = dist(t(assay(DESeq.rlog)))
plot(hclust(dists)) 

#Dispersion of gene expression
dds = estimateDispersions(DESeq.ds)
plotDispEsts(dds)

#Differential expression analysis
dds = nbinomWaldTest(dds)
res = results(dds)
res$padj = ifelse(is.na(res$padj), 1, res$padj)
write.table(res, col.names=NA, row.names=T, file ="expressions.tsv", sep ="\t")
plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)

#Volcano plot
tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$pvalue))
head(tab)
par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
lfc = 2
pval = 0.05
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

#Printing the list of differentially expressed genes
(resOrdered <- res[order(res$padj), ])
write.table(resOrdered, file = "genes.txt", sep ="\t")

#making the aheatmap()
#Install
install.packages("NMF")
#Load
library(NMF)

#sort the results according to the adjusted p-value
res.sorted <- res[order(res$padj), ]

#identify genes with desireed adjusted p-value cut-off
DGEgenes <- rownames(subset(res.sorted, padj < 0.05))

#extract the normalized read counts for DE genes into a matrix
hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

#plot the normalized read counts of DE genes sorted by the adjusted p-value 
aheatmap(hm.mat_DGEgenes,Rowv = NA, Colv = NA)

#Combine the heatmap with hierchical clustering
aheatmap(hm.mat_DGEgenes,Rowv = TRUE, Colv =TRUE, distfun = "euclidean", hclustfun = "average")

#scale the read counts per gene 
aheatmap(hm.mat_DGEgenes,Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row")

#Principal component analysis graph
z <- plotPCA(vsd)
#points with same sample name
z + geom_label(aes(label = name))
nudge <- position_nudge(y = 1)
z + geom_label(aes(label = name), position = nudge)
#different labels
z + geom_text(aes(label = name), position = nudge)

