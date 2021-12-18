#####################################################
# Authors: Jeroen Lodewijk <jeroen.lodewijk@wur.nl> #
#          Thierry Haddad <thierry.haddad@wur.nl>   #
# Student numbers: 930203525030                     #
#                  911031296110                     #
# Subject: Molecular Systems Biology                #
# Dataset: 22                                       #
#####################################################

####################
#### Libraries #####
####################

source("http://bioconductor.org/biocLite.R")
biocLite("minet")    # Network creation
biocLite("impute")   # KNN imputation
biocLite("pheatmap") # Nice heatmaps. 
install.packages("RColorBrewer") # Heatmap colors

library("minet")
library("impute")
library("pheatmap")
library("RColorBrewer") 

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

# Raw fold-change data used for network creation/inference, loaded as matrix
raw.dataset.22 <- as.matrix(read.table("set_22.csv", header = T))

# Meta-information about dimensions
nr.genes <- nrow(raw.dataset.22)
nr.cons <- ncol(raw.dataset.22)
print(paste(
  paste("Total number genes is",  nr.genes),
  " and we have so many conditions", nr.cons
))
# Small-n,  large-p?
if (nr.genes < nr.cons) {
  print("This seems to be a small-n, large-p problem.")
}

################################
# KNN-imputation for NA values #

# Quantify presence of NA values in dataset
na.presence <- any(is.na(raw.dataset.22))
print(paste("Presence of NA values:", na.presence))
if (na.presence) {
  print(paste("Number of NA values:", sum(is.na(raw.dataset.22)) ))
}

# Impute the NA values based on KNN with neighbours = 3
knn.result <- impute.knn(raw.dataset.22, k=3, rowmax=0.3, colmax=0.3, maxp=1500)
dataset.22 <- knn.result$data  # Save as new dataset, leave raw set intact

# Debug code to find and compare before and after values for KNN-changed values
na.vals <- raw.dataset.22[rowSums(is.na(raw.dataset.22)) > 0,]  # List of rows with NA values
na.names <- row.names(na.vals)  # Row names
new.vals <- dataset.22[na.names,]  # The individual rows can be checked for change


###########################
# Check for normalization #

# Increase font size for readability
par(cex = 1.50)

# Transpose dataset for histogram and normal cruves
dataset.22 <- t(dataset.22)

# Create a histogram to check if the  data is normalized
hist.22 <- hist(dataset.22,
                breaks=50,
                main = "Dataset number 22 Histogram\nWith a normal curve",
                xlab = "Fold Changes between conditions and reference")

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
      type = "b",
      lwd = 3) 

# Re-transpose the matrix to its original state
dataset.22 <- t(dataset.22)


#########################
# Preliminary analysis  #
#########################

###########################
#   Clustering of genes   #

# Separate dendrogram in two groups for a more clear representation of genes

# Divide into two main groups
d.groups <- cutree(hclust(dist(dataset.22)), k = 2, h = NULL)

# Assign groups to respective var
group.1 <- which(d.groups == 1)

# But there is not really a second group, so this idea was dropped
#group.2 <- which(d.groups == 2)

# Group 1 dendrogram + correlation
create.dendogram.plot(group.1, 
                      dataset.22, 
                      plot.title = "Dataset 22"
)
create.dendogram.correlation.plot(group.1, 
                                  dataset.22, 
                                  plot.title = "Dataset 22\nBased on correlation between genes"
)


# Perform ierarchical clustering on the dataset
hierar.cl <- hclust(dist(dataset.22))

# Plot a 'normal' dendrogram
# Remains hard to read due to the excessive amount of gene names
plot(hierar.cl,
     xlab = "",
     hang = -1,
     main = "Dendrogram of dataset 22")


################################
# Principal Component Analysis #

# NOTE:
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
xlab<-paste('PC1 ',format(pcas[1], digits = 3),'%')
ylab<-paste('PC2 ',format(pcas[-1], digits = 3),'%')
mainlab<-paste('PCA on dataset 22 (',format(sum(pcas), digits = 2),'% var explained)',sep = "")

# Biplot plot of the PCA to show the variance explained by PC1 and PC2
biplot(data.pca,
       xlab = xlab,
       ylab = ylab, 
       main = mainlab,
       var.axes = TRUE,
       pc.biplot = TRUE)

#############
# Heatmaps  #

heatmap.colors <- brewer.pal(10,"RdBu")  # Red to blue colour range
heatmap.colors<-heatmap.colors[10:1]  # Colour order, red=high, blue=low

# Heatmap based on log2 fold-change
heatmap(t(dataset.22), 
        col = heatmap.colors, 
        scale = "none")

# Heatmap based on correlation between genes, scaling not needed due to symmetrical correlation
heatmap(cor(t(dataset.22)),
        symm = T,
        col = heatmap.colors)


#####################
# Network inference #
#####################

# NOTE: 3 network inference methods will be created and tested
# These are the mrnet, aracne and clr methods:
# - mrnet is based on maximum relevance/minimum redundancy (MRMR) gene selection
# - Aracne tests gene triplets for data inequality through thresholds
# - CLR performs pair-wise comparisons of mutual information (MI) values

nr.bins <- sqrt(nrow(dataset.22))  # nr of bins for discrete dataset
mrnet.net <- minet(dataset.22, method = "mrnet", nbins = nr.bins)
arac.net <- minet(dataset.22, method = "aracne", nbins = nr.bins)
clr.net <- minet(dataset.22, method = "clr", nbins = nr.bins)

###################################################################################
# Below the various networks will be written to csv files for usage in Ctyoscape  #

### CLR
clr.th <- clr.net
diag(clr.th) <- 0  # Set diagonal to 0
clr.th[ lower.tri(clr.th) ] <- 0  # Nullify entries according to lower triangle
clr.th[ which(clr.th > 0) ] <- 1  # Rest becomes 1
clr.ones <- which(clr.th == 1, arr.ind = T) # Return the array indices for the positons

# Write to csv file
write.table(
  cbind(
    rownames(clr.th)[clr.ones[,1]],
    rownames(clr.th)[clr.ones[,2]]),
  file = "clr_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

### ARACNE
arac.th <- arac.net
diag(arac.th) <- 0  # Set diagonal to 0
arac.th[ lower.tri(arac.th ) ] <- 0  # Nullify entries according to lower triangle
arac.th[ which(arac.th > 0) ] <- 1  # Rest becomes 1
arac.ones <- which(arac.th == 1, arr.ind = T)  # Return the array indices for the positons

# Write to csv file
write.table(
  cbind(
    rownames(arac.th)[arac.ones[,1]],
    rownames(arac.th)[arac.ones[,2]]),
  file = "aracne_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

### MRNET
mrnet.th <- mrnet.net
diag(mrnet.th) <- 0  # Set diagonal to 0
mrnet.th[ lower.tri(mrnet.th ) ] <- 0  # Nullify entries according to lower triangle
mrnet.th[ which(mrnet.th > 0) ] <- 1  # Rest becomes 1
mrnet.ones <- which(mrnet.th == 1, arr.ind = T)  # Return the array indices for the positons

# Write to csv file
write.table(
  cbind(
    rownames(mrnet.th)[mrnet.ones[,1]],
    rownames(mrnet.th)[mrnet.ones[,2]]),
  file = "mrnet_network.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

#######################
# Fixing a threshold  #

# NOTE: Here we try a different threshold to test the effect on the networks rigidness
#       The performed steps are identical to above, just with a different threshold

# Set the threshold
thrsh <- 0.25

# Replace the values of clr network results
clr.th2 <- clr.net
clr.th2[which(clr.net <= thrsh)] <- 0
clr.th2[which(clr.net > thrsh)] <- 1
diag(clr.th2) <- 0
clr.th2[lower.tri(clr.th2 )] <- 0
clr.ones2 <- which(clr.th2 == 1, arr.ind = T)

# Write to csv file
write.table(
  cbind(
    rownames(clr.th2)[clr.ones2[,1]],
    rownames(clr.th2)[clr.ones2[,2]]),
  file = "clr_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

# Replace the values of mrnet network results
mrnet.th2 <- mrnet.net
mrnet.th2[which(mrnet.net <= thrsh)]  <- 0
mrnet.th2[which(mrnet.net > thrsh)] <- 1
diag(mrnet.th2) <- 0
mrnet.th2[lower.tri(mrnet.th2 )] <- 0
mrnet.ones2 <- which(mrnet.th2 == 1, arr.ind = T)

# Write to csv file
write.table(
  cbind(
    rownames(mrnet.th2)[mrnet.ones2[,1]],
    rownames(mrnet.th2)[mrnet.ones2[,2]]
  ),
  file="mrnet_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

# Replace the values of aracne network results
arac.th2 <- arac.net
arac.th2[which(arac.net<= thrsh)] <- 0
arac.th2[which(arac.net> thrsh)] <- 1
diag(arac.th2) <- 0
arac.th2[lower.tri(arac.th2 )] <- 0
arac.ones2 <- which(arac.th2 == 1, arr.ind = T)

# Write to csv file
write.table(
  cbind(
    rownames(arac.th2)[arac.ones2[,1]],
    rownames(arac.th2)[arac.ones2[,2]]
  ),
  file = "aracne_th.csv", 
  col.names = F, 
  row.names = F, 
  sep = ",", 
  quote = F)

