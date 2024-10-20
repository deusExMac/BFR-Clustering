


###############################################################################################
#
#
#  Script executes BFR algorithm on some synthetic data files in directory
#  testdata/. All features of the synthetic data are used for clustering. File names
#  indicate the distribution the data were drawn from (N) number of rows (500, 4000) and
#  number of features (x2, x5).
#
#  The results of BFR clustering are evaluated using the Silhouette sore [1, 2, 3] and 
#  displayed with the Silhouette score of a clustering using the traditional K-means algorithm. 
#  Also, and only for two dimensional data only, the centers of the resulting clusters are 
#  shown on a graph  along with the centers of clusters when kmeans is executed on the same data.
#  
#  Hyper-parameters for BFR clustering can be found in file Clustering-BFRLibrary.R .
#
#  References:
#     1) Silhouette_(clustering), https://en.wikipedia.org/wiki/Silhouette_(clustering)
#     2) How to Evaluate the Performance of Clustering Algorithms Using Silhouette Coefficient, 
#        Available from: https://medium.com/@MrBam44/how-to-evaluate-the-performance-of-clustering-algorithms-3ba29cad8c03
#     3) How do you manually compute for silhouette, cohesion and separation of Cluster, 
#        https://stackoverflow.com/questions/23387275/how-do-you-manually-compute-for-silhouette-cohesion-and-separation-of-cluster   
# 
#  Required libraries: 
#     1) LaF (for chunk-based reading of files)
#     2) cluster (for silhouette calculation)
#     3) BFR_Implementation.R (implementation of BFR)
#
#
#  v0.2/mmt/Sep 2024
#
###############################################################################################



# Close any graphic device left open
graphics.off()

# Cleanup the environment
# all=TRUE in ls() is optional
rm(list = ls())

# Clear console
cat("\014")

options(scipen=8)

# Fix this
setwd("MyCode/BFR-Clustering/")




library('LaF')
library(cluster)

# Contains the BFR implementation used by this script 
source('BFR_Implementation.R')


# Number of clusters. The J
N_CLUSTERS <- 14

# Reading this number of lines during
# each read of file
CHUNK_SIZE <- 1000

# File containing data to cluster using BFR.
# Relative path.
DATA_FILE <- "testdata/syntheticDataN4000x2.csv"





#' Returns the next chunk of data, as a data frame, to be passed to BFT for clustering. 
#' 
#' @details the function reads the file on disk in chunks. The file, and specifically LAFCON,
#'          needs to be properly initialized for reading using the LaF library. 
#'          
#' @param chunkSize integer. The number of lines to read from the csv file.                      
#' 
#' @return a DATA FRAME of size at most CHUNK_SIZE of the source data representing the next chunk to
#'         be clustered using BFR. If no more chunks can be read, the function returns NULL.
#'
NextFileChunk <- function(chunkSize=CHUNK_SIZE){
  if (is.null(LAFCON)){
      message('No valid file connection. Cannot read file')
      stop('Terminating.')
  } 
  cdata <- next_block(LAFCON, nrows=chunkSize)
  if (nrow(cdata)==0){
    return(NULL)
  }
  
  return(cdata)
}


cat('####################################################################\n')
cat('Executing BFR with following hyperparameters:\n')
cat('\tNumber of clusters:', N_CLUSTERS, '\n')
cat('\tMinimum cluster points (CS generation):', MINIMUM_CLUSTER_POINTS, '\n')
cat('\tData file:', DATA_FILE, '\n')
model <- detect_dm_csv(DATA_FILE, sep=",", header=TRUE)
cat('\tNumber of features:', nrow(model$columns), '\n')
cat('\tChunk size:', CHUNK_SIZE, '\n')
cat('####################################################################\n')


cat('Will start in 5 seconds\n')
for (i in 1:5){
     cat('.')
     Sys.sleep(1)
}
Sys.sleep(1)


##########################################################################################
# Executing BFR to cluster data  
##########################################################################################



# Preparing for chunked reading of file residing on disk. The idea is
# that file is too large to fit into main memory.
# First detect a data model for your file
model <- detect_dm_csv(DATA_FILE, sep=",", header=TRUE)
LAFCON <- laf_open(model)
# Make sure we are at the start of file
goto(LAFCON, 1)

# Execute BFR.
# Returned value is a list containing item named DS that is a list containing
# the summaries of all DS clusters. Cluster summaries are indexed by sequential
# numbers starting from 1. Index name is a character data type, not a numeric 
# one (i.e. [['1']] nor [[1]])
# Print res[['DS']] to see structure.
res <- BFR(callback=NextFileChunk, K=N_CLUSTERS)

cat('\nBFR finished.\n')

close(LAFCON)


##########################################################################################
# Evaluation: comparing BFR centers with K-means centers. Calculating Silhouette scores 
# and charting BFR and K-means centroids only for 2D data.  
##########################################################################################


# Executing K-means on source data. Entire data must be read into main memory.
# Using same hyperparameters/settings with BFR. 
cat('Executing K-means over same data\n')

data <- read.csv(DATA_FILE, sep=",", header=T)

# K-means with data in memory.
classicKmeans <- kmeans(data, centers=as.matrix(data[sample(nrow(data), N_CLUSTERS), ]), iter.max=MAX_ITERATIONS)


cat('Calculating silhouette scores (K-means, BFR)\n')


# Calculate silhouette score for K-means clustering
kmeansSilhouetteScore <- silhouette(classicKmeans$cluster, dist(data))
cat('\tK-means silhouette score:', mean(kmeansSilhouetteScore[,3]), '\n')

# Calculate silhouette score for BFR clustering
sc <- c()
for (i in 1:N_CLUSTERS){
  cat('Silhouette score for cluster', i, '\n')
  # Use only 6 other clusters to calculate Silhouette score
  sc <- c( sc, ..DBGClusterSilhouetteScore(i, 6))
}

cat('\n\tBFR silhouette score:', mean(sc))




# Plot centers of BFR and result of plain simple Kmeans applied on the same data.
# This is done only if data dimension is 2.
if (res[['DS']]$dim == 2){
    plot(-4:4, -4:4, type = "n")
    # centers from BFR
    cent <- ClusterCenters(res[['DS']])
    points(cent[, 1], cent[, 2], col='orange', pch=15, cex=1.2)
    
    # kmeans centers
    points(classicKmeans$centers[, 1], classicKmeans$centers[, 2], col='grey', pch=16, cex=1.2)
    
    legend("topright", legend=c("BFR centers", "K-means centers"),
           col=c("orange", "grey"), fill = c("orange", "grey") , cex=0.57)
}




