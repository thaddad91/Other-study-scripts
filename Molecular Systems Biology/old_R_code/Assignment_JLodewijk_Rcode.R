###################################################
# Authors: Jeroen Lodewijk <jeroen.lodewijk@wur.nl> #
# Student numbers: 930203525030    #
# Subject: Molecular Systems Biology              #
#                               #
###################################################

#############
# Libraries #

# Library used for getting the colors in the heatmap
library(RColorBrewer)

# Is used to perform the network inference
library(minet)


####################
# Data processing  #

# Read the csv in as table, this is done  to get the proper formatting. Then transform it as a matrix.
#dataset <- as.matrix(read.table("set_17.csv", header = T))

# Print the results to the screen
print(paste( paste("Number of genes is",  dim(dataset)[1])," and the number of conditions is", dim(dataset)[2]) )
print(paste("Are there NA values in the dataset:", any(is.na(dataset))) )

# Check if there are an missing values in the data matrix
any(is.na(dataset))

# Get all rows(genes) , with NA values
new_ds <- dataset[rowSums(is.na(dataset))>0, ]

# 2 genes in total.
nrow(new_ds)

# Check which values are NA in both rows
which(is.na(new_ds[1,])==TRUE)
which(is.na(new_ds[2,])==TRUE)

# Delete all rows with NA values
data.cleared <- na.omit(dataset)

# Check dimensions of matrix without missing values
dim(data.cleared)

# Transpose the matrix.
t.data.cleared <- t(data.cleared)

#######################
# Check normalization #

# Make the all the text larger for the report, otherwise you can not read it.
par(cex = 1.50)

# Create a histogram to check if the  data is normalized
h <- hist(data.cleared,
          main = "Dataset number 17 Histogram\nWith a normal curve",
          xlab = "Fold Changes between studied and reference")

# Calculates the normal curve for the x-axis
xfit <- seq(min(data.cleared), 
            max(data.cleared),
            length = 40)
# Calculates the normal curve for the y-axis
yfit <- dnorm(xfit, 
              mean = mean(data.cleared), 
              sd = sd(data.cleared) )
yfit <- yfit*diff(h$mids[1:2])*length(data.cleared)

# Add the normal curve in the histogram
lines(xfit, 
      yfit, 
      col = "red", 
      lwd = 3) 

#########################
# Preliminary analysis  #

###########################
#   Clustering of genes   #
###########################

# Divide big cluster in two groups for more clear dendrogram without overlapping gene names
# Divide tree in two main groups
groups <- cutree(hclust(dist(data.cleared)), k = 2, h = NULL)

# Fet all genes in group 1
group1 <- which(groups == 1)

# Get all genes in group 2
group2 <- which(groups == 2)

# Get all rows in first cluster out of big dataset
group1 <- data.cleared[group1,]

# Cluster first group of data.
hc_g1 <- hclust(dist(group1))

# To change the size of the labels in the dendrogram.
par(cex = 1.35)

# Get the labels for the hc group 1 and remove some of them, otherwise the plot is unreadable in the report.
labels.group.1 <- hc_g1$labels
labels.group.1[seq(1, length(hc_g1$labels), 2)] <- ""
labels.group.1[seq(2, length(hc_g1$labels), 4)] <- ""
labels.group.1[match( c("Rv0280", "Rv1087", "Rv2685"), labels.group.1)] <- ""

# Plot cluster 1 of dataset 17
plot(hc_g1, 
     hang = -1, 
     main = "Dendrogram of dataset 17 cluster 1", 
     sub = NA,
     labels = labels.group.1,
     xlab = "genes")

hc_g1$labels[seq(1, length(hc_g1$labels), 2)]

# Get all rows in second cluster out of big dataset
group2 <- data.cleared[group2,]

# Cluster second group of data
hc_g2 <- hclust(dist(group2))

# Get the labels for the hc group 1 and remove some of them, otherwise the plot is unreadable in the report.
labels.group.2 <- hc_g2$labels
labels.group.2[seq(1, length(hc_g2$labels), 2)] <- ""

# Plot cluster 2 of dataset 17
plot(hc_g2, 
     hang = -1, 
     main = "Dendrogram of dataset 17 cluster 2", 
     sub = NA,
     labels = labels.group.2,
     xlab = "genes")

# Calculate the distance of correlations between the genes.
distcorr<-as.dist(1-cor(group1))

# Cluster first group of data.
hc_g1_cor <- hclust(as.dist(1-cor(t(group1))))

# Plot cluster 1 of dataset 17
plot(hc_g1_cor, 
     hang = -1, 
     main = "Dendrogram of dataset 17 cluster 1\nBased on correlation between genes", 
     xlab = "genes")

# Get all rows in second cluster out of big dataset
group2 <- data.cleared[group2,]

# Cluster second group of data
hc_g2_cor <- hclust(as.dist(1-cor(t(group2))))

# Plot cluster 2 of dataset 17
plot(hc_g2_cor, 
     hang = -1, 
     main = "Dendrogram of dataset 17 cluster 2\nBased on correlation between genes", 
     xlab = "genes")

# Calculate the distance of correlations between the genes.
# No seperate groups for this dendrogram
plot(hclust(as.dist(1-cor(t.data.cleared))),
     main = "Dendrogram based on correlation between genes",
     xlab = "",
     hang = -1)

# Prepare hierarchical cluster
hc <- hclust(dist(data.cleared))

# Plot a 'normal' dendrogram
plot(hc,
     xlab = "",
     hang = -1,
     main = "Dendrogram of dataset 17")

# Restore to default
par(cex = 1)

##################################
#  Principal Component Analysis  #
##################################

# Do the PCA analysis
prcomp <- prcomp(data.cleared)

# Calculate the explained variance of the first two principal components
var <- (prcomp$sdev)^2
sumvar <- sum(var)
expl1 <- 100*var[1]/sumvar
expl2 <- 100*var[2]/sumvar
explboth <- sum(expl1, expl2)

summary(prcomp)

# Set labels for the biplot
xlab<-paste('PC1 ',format(expl1, digits = 3),'%')
ylab<-paste('PC2 ',format(expl2, digits = 3),'%')
mainlab<-paste('PCA biplot dataset 17 (',format(explboth, digits = 3),'%)',sep = "")

# Plot the pca in a biplot, to shoq the variance of PC1 and PC2
biplot(prcomp,
       xlab = xlab,
       ylab = ylab, 
       main = mainlab)

#############
# Heatmaps  #
#############

# Choose the colors for the heatmap, range is from red to blue
heatmap.colors <- brewer.pal(10,"RdBu")

# Set the order for the colors from high (red) to low values (blue)
heatmap.colors<-heatmap.colors[10:1]

# Heatmap of the data without any NA's
heatmap(data.cleared, 
        col = heatmap.colors, 
        scale = "none")

# This heatmap is based on the correlations between the genes. Since the correlation is symmetrical, no scaling will be needed
heatmap(cor(t.data.cleared),
        symm = T,
        col = heatmap.colors)

#####################
# Network inference #

# Set the number of bins that are going to be for the discretization of the dataset
nbins <- sqrt(nrow(data.cleared))  

# Uses mutual information matrix as input, then tries to infer the network using maximum relevance/minimum redundancy
mrnet.network <- minet(t.data.cleared, 
                       method = "mrnet", 
                       estimator = "mi.empirical", 
                       disc = "equalfreq", 
                       nbins = nbins)

# Uses mutual information matrix as input, then tries to infer the network using Aracne algorithm
aracne.network <- minet(t.data.cleared, 
                        method = "aracne", 
                        estimator = "mi.empirical",
                        disc = "equalfreq", 
                        nbins = nbins)

# Uses mutual information matrix as input, then tries to infer the network using CLR algorithm
clr.network <- minet(t.data.cleared, 
                     method = "clr", 
                     estimator = "mi.empirical", 
                     disc = "equalfreq", 
                     nbins = nbins)

# Write the results for the clr.network so that Cytoscape can use it
clr.th <- clr.network

# Set all diagonal of a matrix to 0
diag(clr.th) <- 0

# Replace all values which entries return TRUE for the lowe triangle
clr.th[ lower.tri(clr.th) ] <- 0

# Replace all values with 1 if the entrie is higher than 0
clr.th[ which(clr.th > 0) ] <- 1

# Get all the 1 out of clr.th.
xx <- which(clr.th == 1, 
            arr.ind = T) # Return the array indices for the positons

# Write the results in a csv file for Cytoscape
write.table(
  cbind(
    rownames(clr.th)[xx[,1]],
    rownames(clr.th)[xx[,2]]),
  file = "clr_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)

# Write the results for the aracne.network so that Cytoscape can use it
aracne.th <- aracne.network

# Set all diagonal of a matrix to 0
diag(aracne.th) <- 0

# Replace all values which entries return TRUE for the lowe triangle
aracne.th[ lower.tri(aracne.th ) ] <- 0

# Replace all values with 1 if the entrie is higher than 0
aracne.th[ which(aracne.th > 0) ] <- 1

# Get all the 1 out of clr.th.
xx <- which(aracne.th == 1, arr.ind = T)

# Write the results in a csv file for Cytoscape
write.table(
  cbind(
    rownames(aracne.th)[xx[,1]],
    rownames(aracne.th)[xx[,2]]),
  file = "aracne_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)

# Write the results for the mrnet.network so that Cytoscape can use it
mrnet.th <- mrnet.network

# Set all diagonal of a matrix to 0
diag(mrnet.th) <- 0 

# Replace all values which entries return TRUE for the lowe triangle
mrnet.th[ lower.tri(mrnet.th ) ] <- 0

# Replace all values with 1 if the entrie is higher than 0
mrnet.th[ which(mrnet.th > 0) ] <- 1

# Get all the 1 out of clr.th.
xx <- which(mrnet.th == 1, arr.ind = T)

# Write the results in a csv file for Cytoscape
write.table(
  cbind(
    rownames(mrnet.th)[xx[,1]],
    rownames(mrnet.th)[xx[,2]]),
  file = "mrnet_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)

#######################
# Fixing a threshold  #

# Set the threshold
thrsh <- 0.25

# Replace the values of clr network results
clr.th2 <- clr.network
clr.th2[which(clr.network <= thrsh)] <- 0
clr.th2[which(clr.network > thrsh)] <- 1

# Set all diagonal of a matrix to 0
diag(clr.th2) <- 0

# Replace all values which entries return TRUE for the lowe triangle
clr.th2[lower.tri(clr.th2 )] <- 0

# Get all the 1 out of clr.th2
xx <- which(clr.th2 == 1, arr.ind = T)

# Write the results in csv file for analysis
write.table(
  cbind(
    rownames(clr.th2)[xx[,1]],
    rownames(clr.th2)[xx[,2]]),
  file = "clr_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)

# Replace the values of mrnet network results
mrnet.th2 <- mrnet.network
mrnet.th2[which(mrnet.network <= thrsh)] <- 0
mrnet.th2[which(mrnet.network > thrsh)] <- 1
diag(mrnet.th2) <- 0

# Replace all values which entries return TRUE for the lowe triangle
mrnet.th2[lower.tri(mrnet.th2 )] <- 0

# Get all the 1 out of mrnet.th2
xx <- which(mrnet.th2 == 1, arr.ind = T)

# Write the results in csv file for analysis
write.table(
  cbind(
    rownames(mrnet.th2)[xx[,1]],
    rownames(mrnet.th2)[xx[,2]]
  ),
  file="mrnet_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)

# Replace the values of aracne network results
aracne.th2 <- aracne.network
aracne.th2[which(aracne.network<= thrsh)] <- 0
aracne.th2[which(aracne.network> thrsh)] <- 1

# Set all diagonal of a matrix to 0
diag(aracne.th2) <- 0
aracne.th2[lower.tri(aracne.th2 )] <- 0

# Get all the 1 out of aracne.th2
xx <- which(aracne.th2 == 1, arr.ind = T)
write.table(
  cbind(
    rownames(aracne.th2)[xx[,1]],
    rownames(aracne.th2)[xx[,2]]
  ),
  file = "aracne_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = "\t", 
  quote = F)
