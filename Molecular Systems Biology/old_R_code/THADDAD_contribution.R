###################################################
# Authors: Jeroen Lodewijk <jeroen.lodewijk@wur.nl> #
#          Thierry Haddad <thierry.haddad@wur.nl>
# Student numbers: 930203525030    #
# Subject: Molecular Systems Biology              #
# Dataset:                               #
###################################################

####################
#### Libraries #####
####################

source("http://bioconductor.org/biocLite.R")
biocLite("minet")
biocLite("impute")
library("minet")
library("impute")

####################
#### Functions #####
####################



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

# Debug code to find and compare before and after values for KNN-changed values
#na.vals <- raw.dataset.22[rowSums(is.na(raw.dataset.22)) > 0,]
#na.names <- row.names(na.vals)

# Replace NA values by imputing with KNN with k=3
knn.result <- impute.knn(raw.dataset.22, k=3, rowmax=0.3, colmax=0.3, maxp=1500)
dataset.22 <- knn.result$data


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
