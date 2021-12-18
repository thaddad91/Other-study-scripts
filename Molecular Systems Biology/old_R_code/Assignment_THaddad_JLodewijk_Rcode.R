###################################################
# Authors: Jeroen Lodewijk <jeroen.lodewijk@wur.nl> #
#          Thierry Haddad <thierry.haddad@wur.nl>
# Student numbers: 930203525030    #
# Subject: Molecular Systems Biology              #
# Dataset:                               #
###################################################

# RECOMMENT threshold CODE!!!!

####################
#### Libraries #####
####################

source("http://bioconductor.org/biocLite.R")
biocLite("minet")
biocLite("impute")
biocLite("pheatmap") #to build nice heatmaps. 

library("minet")
library("impute")

####################
#### Functions #####
####################

#' create.dendogram.plot
#' Creates dendogram plots for a cluster
#' @param group: A cluster produced by cluster()
#' @param dataset: Dataset containing all the results
#' @param plot.title: String that is part of the main of plot()
#' @return None
create.dendogram.plot <- function(group, dataset, plot.title = ""){
  
  group <- dataset[group,]
  
  hc_g <- hclust(dist(group))
  
  labels.group <- hc_g$labels
  labels.group[seq(1, length(hc_g$labels), 2)] <- ""
  labels.group[seq(2, length(hc_g$labels), 4)] <- ""
  
  plot(hc_g, 
       hang = -1, 
       main = plot.title, 
       sub = NA,
       labels = labels.group,
       xlab = "Genes names")
  
}
#' create.dendogram.plot
#' Creates correlation dendogram plots
#' @param group: A cluster produced by cluster()
#' @param dataset: Dataset containing all the results
#' @param plot.title: String that is part of the main of plot()
#' @return None

create.dendogram.correlation.plot <- function(group, dataset, plot.title = ""){
  
  group <- dataset[group,]
  
  hc.g.cor <- hclust(as.dist(1-cor(t(group))))
  
  labels.group <- hc.g.cor$labels
  labels.group[seq(1, length(hc.g.cor$labels), 2)] <- ""
  labels.group[seq(2, length(hc.g.cor$labels), 4)] <- ""
  
  plot(hc.g.cor, 
       hang = -1, 
       main = plot.title, 
       sub = NA,
       labels = labels.group,
       xlab = "Genes names")
  
}


#########################
#### Data processing ####
#########################

# Read in the CSV as a data.frame and with a proper header.
raw.dataset.22 <- as.matrix(read.table("set_22.csv", header = T))


# Print the results to the screen
print(paste(
  paste("Total number genes is",  nrow(raw.dataset.22)),
  " and we have so many conditions",
  ncol(raw.dataset.22)
))

# checking for missing data points in the data set
print(paste("Are there NA values in the dataset:", any(is.na(raw.dataset.22))))
print(paste("How many NA values in the dataset:", sum(is.na(raw.dataset.22)) ))

# impute.knn(raw.dataset.22, k=10, rowmax=0.3, colmax=0.3, maxp=1500)
# knn.result <- impute.knn(as.matrix(raw.dataset.22))
# 
# subset(knn.result, rownames(knn.result) %in% c("Rv1660", "Rv1801"))
# 
# 
# test <- raw.dataset.22[rowSums(is.na(raw.dataset.22)) > 0,]


# Remove all NAs
knn.result <- impute.knn(raw.dataset.22, k=3, rowmax=0.3, colmax=0.3, maxp=1500)
dataset.22 <- knn.result$data

# Debug code to find and compare before and after values for KNN-changed values
#na.vals <- raw.dataset.22[rowSums(is.na(raw.dataset.22)) > 0,]
#na.names <- row.names(na.vals)

# dataset.22 <- na.omit(raw.dataset.22)


#######################
# Check normalization #

# Increase font size.
par(cex = 1.50)

# Transpose dataset, otherwise it can not be used for normalization.
dataset.22 <- t(dataset.22)

# Create a histogram to check if the  data is normalized
hist.22 <- hist(dataset.22,
                breaks=50,
          main = "Dataset number 22 Histogram\nWith a normal curve",
          xlab = "Fold Changes between studied and reference")
# Calculates the normal curve for the x-axis
xfit <- seq(min(dataset.22), 
            max(dataset.22))
# Calculates the normal curve for the y-axis
yfit <- dnorm(xfit, 
              mean = mean(dataset.22), 
              sd = sd(dataset.22) )
yfit <- yfit*diff(hist.22$mids[1:2])*length(dataset.22)
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
groups <- cutree(hclust(dist(dataset.22)), k = 2, h = NULL)

group.1 <- which(groups == 1)
group.2 <- which(groups == 2)

create.dendogram.plot(group.1, dataset.22, plot.title = "Group 1 dataset 22")
create.dendogram.plot(group.2, dataset.22, plot.title = "Group 2 dataset 22")


create.dendogram.correlation.plot(group.1, dataset.22, plot.title = "Group 1 dataset 22\nBased on correlation between genes")
create.dendogram.correlation.plot(group.2, dataset.22, plot.title = "Group 2 dataset 22\nBased on correlation between genes")


# Prepare hierarchical cluster
hc <- hclust(dist(dataset.22))

# Plot a 'normal' dendrogram
plot(hc,
     xlab = "",
     hang = -1,
     main = "Dendrogram of dataset 22")


###########################
# Redundant PCA           #
###########################
# Due to the high dimensionality of the data and the lack of a reference vs condition set,
# the PCA is deemed non-informational at this point.
# Perform Principal Component Analysis on the log-transformed dataset and calculate the proportions of explained
# variance by the first two principal components, including a summary.
data.pca <- prcomp(dataset.22)
pca.var <- (data.pca$sdev)^2
pca.sumvar <- sum(pca.var)
pcas <- c(100*pca.var[1]/pca.sumvar, 100*pca.var[2]/pca.sumvar)
var.of.pcas <- sum(pcas)
summary(data.pca)
# Eigenvalues of the principal components
data.pca$rotation

# Set labels for the biplot
xlab<-paste('PC1 ',format(expl1, digits = 3),'%')
ylab<-paste('PC2 ',format(expl2, digits = 3),'%')
mainlab<-paste('PCA biplot dataset 17 (',format(explboth, digits = 3),'%)',sep = "")

# Plot the pca in a biplot, to shoq the variance of PC1 and PC2
biplot(data.pca,
       xlab = xlab,
       ylab = ylab, 
       main = mainlab,
       pc.biplot = TRUE)

#############
# Heatmaps  #
#############

# Choose the colors for the heatmap, range is from red to blue
heatmap.colors <- brewer.pal(10,"RdBu")

# Set the order for the colors from high (red) to low values (blue)
heatmap.colors<-heatmap.colors[10:1]

# Heatmap of the data without any NA's
heatmap(dataset.22, 
        col = heatmap.colors, 
        scale = "none")

# This heatmap is based on the correlations between the genes. Since the correlation is symmetrical, no scaling will be needed
heatmap(cor(dataset.22),
        symm = T,
        col = heatmap.colors)

#####################
# Network inference #

# Set the number of bins that are going to be for the discretization of the dataset
nbins <- sqrt(nrow(dataset.22))

mrnet.network <- minet(dataset.22, 
                       method = "mrnet", 
                       nbins = nbins)

# Uses mutual information matrix as input, then tries to infer the network using Aracne algorithm
aracne.network <- minet(dataset.22, 
                        method = "aracne", 
                        nbins = nbins)

# Uses mutual information matrix as input, then tries to infer the network using CLR algorithm
clr.network <- minet(dataset.22, 
                     method = "clr", 
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
  sep = ",", 
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
  sep = ",", 
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
  sep = ",", 
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
  sep = ",", 
  quote = F)

# Replace the values of mrnet network results
mrnet.th2 <- mrnet.network
mrnet.th2[which(mrnet.network <= thrsh)]  <- 0
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
  sep = ",", 
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
  sep = ",", 
  quote = F)
