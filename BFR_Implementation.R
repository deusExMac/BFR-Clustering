


################################################################################################
#
#  Contains functions implementing from scratch the Bradley, Fayyad and Reina (BFR) clustering   
#  algorithm [1, 2]. BFR is a variant of the K-means algorithm that works for multidimensional
#  data that cannot fit into main memory.
#  The script makes the assumptions that the data is read in chunks from an external file or 
#  resource. 
#
#  To call BFR in this implementation, call the BFR function. The function returns the summarized 
#  discard set (DS) of clusters. The implemented algorithm follows the steps
#  outlined in [2, 5]. Some implementation details were based on [3, 4].
#   
#  NOTE: Impmentation details of the BFR algorithm may differ. These details are related
#        to:
#         1) How to determine the criteria for merging CS clusters
#         2) The way to distribute data in the (old) RS to the CS and (new) RS
#
#
#  References: 
#     1) BFR algorithm, Available from https://en.wikipedia.org/wiki/BFR_algorithm
#     2) Rajaraman, Anand; Ullman, Jeffrey; Leskovec, Jure (2011). Mining of Massive Datasets. 
#        New York, NY, USA: Cambridge University Press. pp. 257–258. ISBN 1107015359. 
#        (section 7.3.4, pg 269). Available online from http://www.mmds.org/
#     3) Python BFR implementation, Available from https://github.com/artisan1218/BFR-Clustering/tree/main
#     4) Week 4: Clustering - Part 5: Implementation Details of BFR, Available from  
#        https://www.youtube.com/watch?v=UXUc7CfHsBU
#     5) BFR Algorithm: Extension of K-means to Large Data, Available from
#        https://chih-ling-hsu.github.io/2017/09/01/Clustering#bfr-algorithm-extension-of-k-means-to-large-data
#  
# 
#
#  v0.2/mmt/Sep 2024
#
################################################################################################



# The number of clusters to create
NUMBER_OF_CLUSTERS <- 8

# DEFAULT number of clusters for clustering outliers
# then generating the compression set (CS).
NUMBER_OF_OUTLIER_CLUSTERS <- 3*NUMBER_OF_CLUSTERS


# For kmeans
MAX_ITERATIONS <- 100

# Percentage of intial batch to 
# use to generate DS, CS. RS sets
# Must be in range (0, 1]
INIT_BATCH_PCT <- 1


# Minimum number of points a cluster must have
# in order to not consider the points as outliers
MINIMUM_CLUSTER_POINTS <- 8

# Minimum number of points a cluster must have to add the 
# cluster into the Compression set
MINIMUM_CS_CLUSTER_POINTS <- 80



# How many standard deviations away
# from a cluster a data point is considered close
# to cluster. For Mahalanobis distance.
STDEV_THRESHOLD <- 2.8

..DEBUGVAR <- list('DS'=NULL, 'CS'=NULL)

# Set to FALSE to disable keeping debug info for
# calculating scores at the end.
DEBUG_ENABLED <- TRUE






##########################################################################################
# DEBUG functions  
# The next functions are only used to store ALSO the data inside clusters in order to support 
# evaluation of the final result using various methods such as the Silhouette coefficient.
# 
# The rows in clusters are stored in a similar data structure (list) used to store
# summarized clusters. This has been done since this evaluation part of code is only
# for development purposes; the debug functions in this section will/should be removed 
# once development reaches a stable state/version.
##########################################################################################

..DBGAddDSCluster <- function(df){
  if (!DEBUG_ENABLED){
      return()
  }  
  
  if (is.null(..DEBUGVAR[['DS']])){
      idx <- '1'
  } else{
      idx <- as.character(length(..DEBUGVAR[['DS']])+1)
  }
  ..DEBUGVAR[['DS']][[idx]] <<- rbind(..DEBUGVAR[['DS']][[idx]], df)
  
}



..DBGAddCSCluster <- function(df){
  if (!DEBUG_ENABLED){
    return()
  }
  
  if (is.null(..DEBUGVAR[['CS']])){
    idx <- '1'
  } else{
    idx <- as.character(length(..DEBUGVAR[['CS']])+1)
  }
  ..DEBUGVAR[['CS']][[idx]] <<- rbind(..DEBUGVAR[['CS']][[idx]], df) 
}



..DBGRemoveCSCluster <- function(cid){
  if (!DEBUG_ENABLED){
    return()
  }
  
  ..DEBUGVAR[['CS']][[cid]] <<- NULL
}



..DBGReindexCS <- function(){
  if (!DEBUG_ENABLED){
    return()
  }
  
  names(..DEBUGVAR[['CS']]) <<- 1:length(..DEBUGVAR[['CS']])
}



..DBGAddDataToDS <- function(v, cid){
  if (!DEBUG_ENABLED){
    return()
  }
  ..DEBUGVAR[['DS']][[cid]][nrow(..DEBUGVAR[['DS']][[cid]]) + 1, ] <<- v
}


..DBGAddDataToCS <- function(v, cid){
  if (!DEBUG_ENABLED){
    return()
  }
  ..DEBUGVAR[['CS']][[cid]][nrow(..DEBUGVAR[['CS']][[cid]]) + 1, ] <<- v
}



..DBGAppendDFToDS <- function(srcId, destId){
  if (!DEBUG_ENABLED){
    return()
  }
  
  #..DEBUGVAR[['DS']][[cid]][nrow(..DEBUGVAR[['DS']][[cid]]) + 1, ] <<- v
  
  cat('\tBEFORE DEST', nrow( ..DEBUGVAR[['DS']][[destId]]), '\n')
  cat('\tBEFORE SRC', nrow( ..DEBUGVAR[['DS']][[srcId]]), '\n')
  ..DEBUGVAR[['DS']][[destId]] <<- rbind(..DEBUGVAR[['DS']][[destId]], ..DEBUGVAR[['CS']][[srcId]])
  cat('\tAFTER', nrow( ..DEBUGVAR[['DS']][[destId]]), '\n')
}



..DBGAddDataCS <- function(v, cid){
  if (!DEBUG_ENABLED){
    return()
  }
  ..DEBUGVAR[['CS']][[cid]][nrow(..DEBUGVAR[['CS']][[cid]]) + 1, ] <<- v
}


..DBGMergeCS <- function(srcId, destId){
  if (!DEBUG_ENABLED){
    return()
  }
  ..DEBUGVAR[['CS']][[destId]] <<- rbind(..DEBUGVAR[['CS']][[destId]], ..DEBUGVAR[['CS']][[srcId]])
}


...DBGDataCount <- function(){
   if (!DEBUG_ENABLED){
      return()
   }
  
   total <- 0
   for (i in 1:length(..DEBUGVAR[['DS']])){
        total <- total + nrow(..DEBUGVAR[['DS']][[as.character(i)]])
   }
   return(total)
}



..DBGEuclideanDistance <- function(v, cntrd){
  if (!DEBUG_ENABLED){
    return()
  }
  return( sqrt(sum((v - cntrd)^2)))
}


# TODO: Finish this...
..DBGWithinSS <- function(clusterId){
  if (!DEBUG_ENABLED){
    return()
  }
   centroid <- colSumns(..DEBUGVAR[['DS']][[as.character(clusterId)]]) / nrow(..DEBUGVAR[['DS']][[as.character(clusterId)]])
   return(apply(..DEBUGVAR[['DS']][[as.character(clusterId)]], 1, ..DBGEuclideanDistance, cntrd=centroid))
}






#' Calculates the Silhouette score for one specific data point. 
#' 
#'
#' @details  Takes a point and calculates its silhouette score from the
#'           clusters it is not part of.
#'           Calculation of silhouette score for point follows these steps:
#'           https://medium.com/@MrBam44/how-to-evaluate-the-performance-of-clustering-algorithms-3ba29cad8c03
#'           
#'
#' @param p numeric vector representing point 
#'                  
#' @param cT List of clusters from which to calculate average distances 
#'           from.
#' @return a numeric value (scalar) representing silhouette value of one
#'         data point.
#' 
..DBGPointSihlouetteScore <- function(p, cT, avgDist){
  if (!DEBUG_ENABLED){
    return()
  }
  
  #print(cT)
  cpD <- c()
  for (i in 1:length(cT)){
       #cat(i, ') point distance from', cT[i], '\n')
       cpD <- c(cpD, mean(apply(..DEBUGVAR[['DS']][[cT[i]]], 1, ..DBGEuclideanDistance, cntr = p)))
  }
  
  mDist <- min(cpD)
  return( (mDist - avgDist) / max(mDist, avgDist) )
}





#' Calculates the Silhouette score for a specific cluster. 
#' TODO: Calculation is slow. Need to use vectorized expressions to speed
#'       thnings up.
#' 
#'
#' @details  Takes all points of the cluster and calculates for each
#'           one the distance from the other clusters. Euclidean distance
#'           is used. Number of clusters to calculate distance of points
#'           is given as argument due to its potential large size.
#'           Calculation of silhouette score follows these steps:
#'           https://medium.com/@MrBam44/how-to-evaluate-the-performance-of-clustering-algorithms-3ba29cad8c03
#'           
#'
#' @param clusterId number representing the id of the cluster for which
#'                  to calculate the silhouette score.
#' @param nClust a numeric value indicating the number of other clusters
#'               from which to calculate the distance of cluster points.
#'               Negative value indicates all clusters. 
#' @return a numeric value (scalar) representing the mean value of the 
#'         distances of points to the other clusters.
#' 
..DBGClusterSilhouetteScore <- function(clusterId, nClust=-1){
       if (!DEBUG_ENABLED){
          message('No debug data available')
          return()
       }
    
       if (nrow(..DEBUGVAR[['DS']][[as.character(clusterId)]]) ==1){
           return(0)
       }
  
       #cat('Checking against', nClust, 'clusters. This may take a while....\n')
       avgInterClusterDistance <- mean( dist(..DEBUGVAR[['DS']][[as.character(clusterId)]]))
       if (nClust <= 0){
           targets <- names(..DEBUGVAR[['DS']])[!names(..DEBUGVAR[['DS']]) %in% c(as.character(clusterId))] 
       }else{
           targets <- sample(1:length(names(..DEBUGVAR[['DS']])), nClust + 1)
           if (!clusterId %in% targets){
               targets <- as.character(targets[1:nClust])  
           }else{
               targets <- as.character(targets[!targets %in% c(clusterId)])
           }
           
       }
       srcDF <- ..DEBUGVAR[['DS']][[as.character(clusterId)]]
       
       sScores <- c()
       for (i in 1:nrow(srcDF)){
            #cat(i, ' ')
            if (i%%300==0){
                cat('[', i,'/', nrow(srcDF), ']', sep='')
            }
            #print('')
            #print(names(..DEBUGVAR[['DS']])[!names(..DEBUGVAR[['DS']]) %in% c(as.character(clusterId))])
            pss <- ..DBGPointSihlouetteScore(as.numeric(srcDF[i,]), 
                                             targets, 
                                             avgInterClusterDistance)
            sScores <- c(sScores, pss)
       }
       
       return(mean(sScores))
      
}






##########################################################################################
# Distance function metrics  
##########################################################################################


#' Calculates the Mahalanobis distance between a cluster and a vector. 
#' 
#'
#' @details  Uses the formula in [2, section 7.3.5, page 273]
#'           MDistiance = sqrt()
#'           for each element in the vector.
#'                            
#'
#' @param clst a CLUSTER object. 
#' @param vec a numeric VECTOR . 
#' @return a numeric value (scalar) representing Mahalanobis distance
#'         between cluster and vector.
#'          
MDistance <- function(clst, vec){
  centroid <- clst[['SUM']] / clst$n
  normalizedVector <- (vec - centroid) / ClusterVariances(clst)
  return(sqrt(sum(normalizedVector^2)))
}




#' Calculates the average Mahalanobis distance between two clusters. 
#' 
#'
#' @details  Calculates Mahalanobis distance between cluster1 and
#'           cluster2 and Mahalabobis distance between cluster2 and cluster1.
#'           Returns their average.
#'                            
#'
#' @param cluster1 a nonempty CLUSTER object. 
#' @param cluster2 a nonempty CLUSTER object. 
#'              
#' @return a numeric value (scalar) with the average distance of the
#'         Mahalanobis distances.
#'  
MClusterDistance <- function(cluster1, cluster2){
  centroid1 <- cluster1[['SUM']] / cluster1$n
  centroid2 <- cluster2[['SUM']] / cluster2$n
  
  d1 <- MDistance(cluster2, centroid1)
  d2 <- MDistance(cluster1, centroid2)
  return((d1+d2)/2)
}


# TODO: Do we need this?
CSClustersShouldMerge <- function(c1, c2, threshold){
  
}

##########################################################################################
# Cluster set and cluster management functions  
##########################################################################################


#' Returns the closest cluster in a cluster set to a vector that meets
#' a minimum distance threshold. 
#' 
#'
#' @details  Calculates distance between all clusters in cluster set and vector  
#'           using the Mahalanobis distance. 
#'           Returns the closest cluster (smallest distance) among those whose 
#'           distance is smaller than a specific threshold. 
#'           Threshold value is STDEV_THRESHOLD*sqrt(dimension_of_data).
#'           
#'                            
#'
#' @param cst a nonempty CLUSTER SET object. 
#' @param vec a numeric VECTOR. 
#'              
#' @return the cluster (cluster id) in the cluster set with the smallest distance smaller
#'         than the threshold to the vector. If no such cluster is found (e.g. all distances
#'         larger than threshold) Inf is returned.
#' 
ClosestCluster <- function(cst, vec){
  cidx <- Inf
  cmd <- Inf
  
  if (is.null(cst)){
      return(cidx)
  }
  
  for (i in 1:cst$N){
    md <- MDistance(cst[[as.character(i)]], vec)
    if (is.na(md)){
        md <- Inf
    }
    #cat('Cluster ', i, ') MD=', md, 'dimensions=', cst[['dim']], '\n')
    if (md <= STDEV_THRESHOLD*sqrt(cst[['dim']])){
      #cds <- append(cds, i)
      if (md < cmd){
        #cat(md, 'smaller than', cmd, '\n')
        cidx <- i
        cmd <- md
      }
    }
  } # for
  # Return index
  return(cidx)
}





#' Returns the variances of all features in a cluster.  
#' 
#'
#' @details  Calculates variances using the cluster's summary statistics using
#'           the formula:
#'           (sumsquares/numberOfObservations) - (sums/numberOfObservations)^2
#'           which is actually the known: E(X^2) - (E(X))^2
#'           
#'                            
#'
#' @param clst a nonempty CLUSTER object. 
#'              
#' @return a VECTOR with the variance for each feature. 
#' 
ClusterVariances <- function(clst){
  return((clst$SUMSQ / clst$n) - (clst$SUM / clst$n)^2)
}






#' Returns a new, empty CLUSTER SET  
#' 
#'
#' @details  a cluster set contains a set of clusters along with some metadata.
#'           A cluster set is represented by list containing a number of named items
#'           that include:
#'           * N: number of clusters in the cluster set
#'           * dim: the dimension of the data (number of features)
#'           * an arbitrary number of summarized clusters, each cluster indexed by a unique
#'           increasing name 1,2,3 etc in the list, in character format (e.g. '1', '2', '3', etc). 
#'           Cluster sets are used to represent the Discard and Compression sets (DS, CS).
#'           The RS is not a cluster set; it's a data.frame
#'                            
#'              
#' @return a LIST containing only an initialized named field for the number of clusters (N), 
#'         and dimentionality of data (dim) representing the cluster set. Cluster set does
#'         not contain any cluster. 
#' 
ClusterSet <- function(){
  return(list('N'=0, 'dim'=0))
}






#' Returns a new summarized CLUSTER object that is initialized by the number of
#' points in the cluster, a vector representing the sum in each dimension and a
#' vector of the sum of squares in each dimension of cluster.  
#' 
#'
#' @details  Clusters objects are represented as lists with the following
#'           3 named items:
#'           * n: the number of items in the cluster 
#'           * SUM: the sum of values in each dimension as a vector
#'           * SUMSQ: the sum of squared values in each dimension as a vector  
#'
#' @param n An integer value specifying the number of observations/points in the
#'          cluster
#' @param s a VECTOR with integers each representing the sum of values of each dimension 
#'        in the cluster data
#' @param sSq a VECTOR with integers each representing the sum of squared values of 
#'        each dimension in the cluster data
#' 
#' @return a named LIST representing a fully summarized cluster (items: n, SUM and SUMSQ) 
#'         based on the centroid data.
#' 
CentroidCluster <- function(n, s=NULL, sSq=NULL){
  #cat('>>> New centroid cluster with', n, 'items:\n')
  return(list('n'=n,
              'SUM'=s,
              'SUMSQ'=sSq))
}


#' Adds a new summarized cluster to an existing cluster set. If the cluster set does not
#' exist, the cluster set is created and the summarized cluster is added to the 
#' newly created cluster.  
#' 
#'
#' @details The new cluster is added as a named item to the list. The new number of
#'          clusters in list (N) is used as the new cluster's name (index).
#'           
#' @param cst a list representing the cluster set where the new cluster should be added
#' @param nData an integer representing the number of items/points in the newly added cluster             
#' @param s a VECTOR with integers each representing the sum of values of each dimension 
#'        in the cluster data
#' @param sSq a VECTOR with integers each representing the sum of squared values of 
#'        each dimension in the cluster data           
#' @return a CLUSTER SET (list) with the new summarized cluster added.  
#' 
AddClusterByCentroid <- function(cst, nData, s=NULL, sSq=NULL){
  
  # Check if cluster set needs to be created
  if (is.null(cst)){
      ncs <- ClusterSet() 
      ncs[['dim']] <- length(s)
      ncs[['N']] <- ncs[['N']] + 1
      ncs[[as.character(ncs[['N']])]] <- CentroidCluster(nData, s, sSq)
      return(ncs)
  }
  
  cst[['N']] <- cst[['N']] + 1
  cst[[as.character(cst[['N']])]] <- CentroidCluster(nData, s, sSq)
  return(cst)
}



#' Adds new data (a point) to an existing cluster in a cluster set. The added data   
#' is summarized.
#'
#' @details Adding the new data will also update properly the named fields SUM
#'          and SUMSQ.
#'           
#' @param cst a list representing the cluster set where the new data should be added
#' @param clusterId an integer specifying the index of the cluster inside the cluster set
#'        to add (and summarize) the data
#' @param newData a vector specifying the data to be added to the cluster
#'               
#' @return a CLUSTER SET (list) with the data added to the specified cluster. Cluster is summarized
#'         with the new data.   
#' 
AddDataToCluster <- function(cst, clusterId, newData){
  
  if (clusterId > cst[['N']]){
    stop('Out of bounds.')
  }
  
  dsId <- as.character(clusterId)
  cst[[dsId]]$n <- cst[[dsId]]$n + 1
  cst[[dsId]][['SUM']] <- cst[[dsId]][['SUM']] + newData
  cst[[dsId]][['SUMSQ']] <- cst[[dsId]][['SUMSQ']] + newData^2
  
  if (any(cst[[dsId]][['SUM']]==0)){
      readline("Zero in SUM vector>")
  }
  
  return(cst)
}



### NEW! For testing!
AddToClosestCluster <- function(df, resultSet){
  
  ds <- resultSet[['DS']]
  cs <- resultSet[['CS']]
  rs <- resultSet[['RS']]
  
  cIdx <- ClosestCluster(ds, df)
  if (!is.infinite(cIdx)){
    ds <- AddDataToCluster(ds, cIdx, df)
    return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  }
  
  cIdx <- ClosestCluster(cs, df)
  if (!is.infinite(cIdx)){
    cs <- AddDataToCluster(cs, cIdx, df)
    return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  }
  
  if (!'F1' %in% colnames(rs) ){
    print(colnames(rs))
    stop('INCONSISTENT HEADER')
  }
  
  
  rs <- rbind(rs, df)
  
  return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  
}






#' Adds new data (a point) to the closest cluster in a cluster set. Closeness is determined using   
#' a distance metric. Searches closest cluster in the clusters sets given as arguments. 
#'
#' @details Cluster sets are scanned for closest cluster in specific order. First, closest
#'          cluster in the ds cluster set is located. If found, data is added and summarized
#'          in this cluster. If no such cluster is found, cluster set cs is scanned for closest
#'          cluster. If found, data is added and summarized in located cluster. If not found, data
#'          is added in the rs data frame. 
#'           
#' @param v a VECTOR representing the data to be added to the closest cluster
#' @param ds an LIST representing a CLUSTER SET. Specifies the Discard Set.
#' @param cs an LIST representing a CLUSTER SET. Specifies the Compression Set.
#' @param rs a DATA FRAME representing the Retained Set.
#'               
#' @return a LIST representing a NEW CLUSTER SET with the data added to the appropriate (closest) cluster.
#'         Cluster sets/data frames not updated with the data are returned as given.   
#' 
AddDataToClosestCluster <- function(v, ds, cs, rs){
   
  # Add to closest DS cluster
  
  cIdx <- ClosestCluster(ds, v)
  if (!is.infinite(cIdx)){
    ds <- AddDataToCluster(ds, cIdx, v)
    #write.table(as.data.frame(t(v)), file=paste0(cIdx,'.csv'), append=TRUE, sep = ",", row.names=FALSE, col.names=!file.exists(paste0(cIdx,'.csv')))
    ..DBGAddDataToDS(v, as.character(cIdx))
    return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  }
  
  
  cIdx <- ClosestCluster(cs, v)
  if (!is.infinite(cIdx)){
    cs <- AddDataToCluster(cs, cIdx, v)
    ..DBGAddDataToCS(v, as.character(cIdx))
    return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  }
  
  if (!'F1' %in% colnames(rs) ){
      print(colnames(rs))
      stop('INCONSISTENT HEADER')
  }
  
  #cat(' before ', nrow(rs), ' ', colnames(rs), '\n')
  #rs <- rbind(rs, v)
  #print('adding')
  #print(as.numeric(v))
  #print(rs)
  #print(str(rs))
  
  rs[nrow(rs)+1, ] <- v
  
  #cat(' AFTER ', nrow(rs), ' ', colnames(rs), '\n')
  #cat(' AFTER rbind',colnames(rs), '\n')
  
  
  
  return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
}




##########################################################################################
# Merging functions  
##########################################################################################


#' Merges two cluster sets: all clusters in the source cluster set will be merged into the 
#' closest cluster in the destination cluster set.    
#'
#' @details Uses the Euclidean distance between the cluster centroids to find the 
#'          closest cluster.
#'           
#' @param sourceC a list representing the source CLUSTER SET whose clusters will be merged into
#'                another cluster. 
#' @param ds a list representing the destination CLUSTER SET where the source clusters will be
#'           merged with.
#'               
#' @return a list which is the destination cluster set where the source cluster set has been
#'         merged into.    
#'
MergeClusterSets <- function(sourceC, destC){
     if (is.null(sourceC)){
         return(destC)
     }
  
     for (i in 1:sourceC$N){
          closestClusterDist <- Inf
          cIdx <- -Inf
          for (j in 1:destC$N){
               # Euclidean distance between centroids
               d <- sqrt(sum(((sourceC[[as.character(i)]]$SUM/sourceC[[as.character(i)]]$n) - (destC[[as.character(i)]]$SUM/destC[[as.character(i)]]$n) )^2))
               if ( d <= closestClusterDist){
                    closestClusterDist <- d
                    cIdx <- as.character(j)
               }
          }
          
          # Merge into closest destination cluster, updating the summaries
          destC[[cIdx]]$n <- destC[[cIdx]]$n + sourceC[[as.character(i)]]$n
          destC[[cIdx]]$SUM <- destC[[cIdx]]$SUM + sourceC[[as.character(i)]]$SUM
          destC[[cIdx]]$SUMSQ <- destC[[cIdx]]$SUMSQ + sourceC[[as.character(i)]]$SUMSQ
          
          ..DBGAppendDFToDS(as.character(i), cIdx)
     }
  
  
     return(destC)
}



#' Merges the points in data frame into a cluster set. Does this by finding the closest cluster in the
#' cluster set to a point/data item in the data frame. Adds the point to the closest cluster and updates
#' the summaries.     
#'
#' @details Calculates the Euclidean distance between the cluster centroid and the data item. Adds the  
#'          point/data item to the closest cluster in cluster set. Cluster summaries are updated. 
#'           
#' @param rs a data frame containing the rows/data/points to be merged
#' @param destC a list representing the cluster set containing the clusters
#'              where the points in rs should be merged with.
#'               
#' @return the list representing the cluster set where all items of the 
#'          data frame have been merged into the closest cluster.
#'
MergeRSToDS <- function(rs, destC){
              
              nMerged <- 0
              for (i in 1:nrow(rs)){
                   cD <- Inf
                   cIdx <- -Inf
                   if (i%%1000==0){
                       cat(' ', i, ' ')
                   }
                   for (j in 1:destC$N){
                     v1 <- as.numeric(rs[i,])
                     v2 <- destC[[as.character(j)]]$SUM/destC[[as.character(j)]]$n
                     d <- sqrt(sum(( v1-v2 )^2))
                     if (d <= cD){
                         cD <- d
                         cIdx <- j
                     } 
                   } # for j
                
                   if (cIdx == -Inf){
                       stop('-Inf index in MergeRSToDS. You Should never see this.')
                   }
                   
                   #write.table(data.frame(t(as.numeric(rs[i,]))), file = paste0(cIdx,'.csv'), append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(paste0(cIdx,'.csv')))
                   
                   # Merge point to choosen Cluster in cluster set and summarize 
                   destC[[as.character(cIdx)]]$n <- destC[[as.character(cIdx)]]$n + 1
                   destC[[as.character(cIdx)]]$SUM <- destC[[as.character(cIdx)]]$SUM + as.numeric(rs[i,])
                   destC[[as.character(cIdx)]]$SUMSQ <- destC[[as.character(cIdx)]]$SUMSQ + as.numeric(rs[i,])^2
                   nMerged <- nMerged + 1
                   #cat('\t\tAdding RS index', i, 'to cluster id:', cIdx, '\n')
                   ..DBGAddDataToDS(as.numeric(rs[i, ]), as.character(cIdx))
              } # for i 
              
              # TODO: This next for loop can be removed
              tot <- 0
              for (k in 1:destC$N){
                   tot <- tot + destC[[as.character(k)]]$n
              }
              
              cat('\n\t\tTotal number of observations in RS merged:', tot, '(ver:', nMerged, ')\n')
              
              return(destC)
}




#' Merges clusters in the cluster set given as argument. 
#' 
#' @details Calculates the Mahalanobis distance between two clusters in the same cluster set and merges
#'          the clusters that have the smallest Mahalanobis distance below a threshold. For all pairs of
#'          clusters, the average Mahalanobis distance is returned. If distances do not meet threshold
#'          criteria, clusters are not merged. Only pairs of clusters are merged: Clusters resulting
#'          from a merge, are ignored during merge operations. 
#'          When two clusters are merged, summaries are updated. Names of clusters in the cluster set 
#'          are updated so that they constitute increasing indexes starting from 1.
#'            
#' @param cs a LIST representing a cluster set
#'               
#' @return the list representing the cluster set (Discard Set) with the merged clusters.
#'
MergeCSCluster <- function(cs){
       
       if (is.null(cs)){
           return(NULL)
       }
  
       if (cs$N < 2){
           return(cs)
       }
  
  
       mergedPairs <- c()
       # Generate all pairs
       allPairs <- combn(1:cs$N, 2)
       for (i in 1:dim(allPairs)[2]){
         
         if ((allPairs[1,i] %in% mergedPairs) || (allPairs[2,i] %in% mergedPairs)){
           #print(mergedPairs)
           #cat('Clusters already merged\n')
           next
         }
         
         # Average Mahalanobis distance between clusters
         d <- MClusterDistance(cs[[as.character(allPairs[1,i])]], cs[[as.character(allPairs[2,i])]])
         if (d < STDEV_THRESHOLD*cs$dim){
             cat('\tMERGING', allPairs[1,i], 'and', allPairs[2,i], '\n')
             mergedPairs <- c(mergedPairs, allPairs[1,i])
             mergedPairs <- c(mergedPairs, allPairs[2,i])
             cs[[as.character(allPairs[1,i])]]$n <- cs[[as.character(allPairs[1,i])]]$n +
                                                    cs[[as.character(allPairs[2,i])]]$n
             cs[[as.character(allPairs[1,i])]]$SUM <- cs[[as.character(allPairs[1,i])]]$SUM +
                                                      cs[[as.character(allPairs[2,i])]]$SUM
             cs[[as.character(allPairs[1,i])]]$SUMSQ <- cs[[as.character(allPairs[1,i])]]$SUMSQ +
                                                        cs[[as.character(allPairs[2,i])]]$SUMSQ
             # Remove merged list item
             cs[[as.character(allPairs[2,i])]] <- NULL
             cs$N <- cs$N - 1
             
              #cat('') 
             ..DBGMergeCS(as.character(allPairs[2,i]), as.character(allPairs[1,i]))
             ..DBGRemoveCSCluster(as.character(allPairs[2,i]))
         } # if
       } # for
       
       
       # Check if clusters need reindexing i.e. new numbers starting from 1
       if (length(mergedPairs) > 0){
           newNames <- c('N', 'dim') 
           
           # TODO: Change Loop to newNames <- c(newNames, as.character(1:cs$N)) 
           #for (i in 1:cs$N) {
           #      newNames <- c(newNames, as.character(i))
           #}
           
           # TODO: Check if next is correct!
           newNames <- c(newNames, as.character(1:cs$N))
           
           oldNames <- names(cs)
           names(cs) <- newNames
           
           cat('\tTotal of [', cs$N, '] clusters\n')
           if (cs$N > 0){
               for (k in 3:cs$N){
                    if (oldNames[k] == newNames[k]){
                        cat('\t\t>>> pos ', k, ': Old index:', oldNames[k], 'New index:', newNames[k], ' UNCHANGED\n')
                    }else{
                       cat('\t\t>>> pos ', k, ': Old index:', oldNames[k], 'New index:', newNames[k], '\n')
                       # TODO: uncomment next and make it work
                       #file.rename(paste0('CS', oldNames[k], '.csv'), paste0('CS', newNames[k], '.csv'))
                    }
                    
               } #for
             ..DBGReindexCS()
           }#if
       }
       
       return(cs)
}


#' Determines if the cluster originating from K-means is considered a Compression Set or
#' a Retained Set. Does this by checking the number of data points in the cluster. if
#' number of points is greater or equal than MINIMUM_CS_CLUSTER_POINTS, all data points are considered
#' close together. A number of points smaller than  MINIMUM_CS_CLUSTER_POINTS indicate
#' outliers.
#' 
#' @details If number of points in cluster is greater or equal than MINIMUM_CS_CLUSTER_POINTS, cluster 
#' will be considered a Compression Set. Number of points smaller than  MINIMUM_CS_CLUSTER_POINTS will
#' add cluster points to the Retained Set.
#'            
#' @param clusterData a data frame with the points belonging to a cluster resulting from
#'                    k-means clustering 
#' @param kmeansClusters the kmeans object returned by kmeans that generating the cluster clusterData 
#'                       belongs to. Unused in this implementation. Kept for future use to make function
#'                       compliant with interface. In future versions, a dynamic way to specify the CS/RS
#'                       generation process is planned.
#'               
#' @return a list with two named fields: CS with the input data frame or NULL, RS with the input
#'         data frame or NULL (depending on cluster criteria). Only one of CS / RS can be NULL.  
#'
CountBasedCSRSGeneration <- function(clusterData, kmeansClusters){
   cat('\t[count-based] Number of rows:', nrow(clusterData), '\n')
   if (nrow(clusterData) >= MINIMUM_CS_CLUSTER_POINTS){
       # this cluster will belong to the CS 
       return(list('CS'=clusterData, 'RS'=NULL))
   }else{
       # add it to the RS
       return(list('CS'=NULL, 'RS'=clusterData))
   }  
      
}




#' Merges data in a data frame into clusters of a new cluster set (cs), if clusters meet criteria. New
#' clusters will be added Compression Set. Clusters that do not meet criteria, remain in the RS.
#' 
#' @details Executes K-means on the data on the data frame. If K-means clusters have more points than 
#'          a threshold, clusters are added to a new cluster set. If K-means clusters have less data
#'          than a threshold, data remains in the data frame. A large K value is used to force more 
#'          detailed/compact clusters.
#'            
#' @param rsdf a DATA FRAME containing the data to be added to clusters
#' @param K number of clusters to generate using K-means. A different and larger value than the number
#'          of clusters to group the source data.
#'               
#' @return a new list, representing a (partial) cluster set with the summarized clusters generated 
#'         and a data frame containing the data that could not be clustered (are outliers). 
#'
MergeDFIntoClusters <- function(rsdf, K=NUMBER_OF_OUTLIER_CLUSTERS){
  cs <- NULL
  rs <- rsdf[0, ]
  
  
  if (nrow(rsdf) == 1){
    rs <- rbind(rs, rsdf)
    return(list('CS'=NULL, 'RS'=rs))
  }
  
  # TODO: Is the next correct?
  # Determine optimal number of clusters
  k <- K
  if (nrow(rsdf) <= k){
    k <-  as.integer(nrow(rsdf)/2)
  }
  
  
  
  csClusters <- kmeans(rsdf, iter.max=MAX_ITERATIONS, centers=k)
  for (i in 1:length(unique(csClusters$cluster))){
    
    clusterD <- rsdf[which(csClusters$cluster==i), ]
    r <- CountBasedCSRSGeneration(clusterD, csClusters)
    # NOTE: CS and RS in r are not clusters; they are data frames.
    if (is.null(r[['RS']])){
        cat('\t[count-based] CS:', nrow(r[['CS']]), '. Adding CS cluster.\n')
        cs <- AddClusterByCentroid(cs, nrow(r[['CS']]), 
                                 colSums(r[['CS']]),
                                 colSums(r[['CS']])^2)
        
        ..DBGAddCSCluster(r[['CS']])
        #write.table(clusterD, file = paste0('CS', i,'.csv'), append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(paste0('CS', i, '.csv')))
    } # if
    else {
          cat('\t[count-based] RS:', nrow(r[['RS']]), '. Adding RS\n')
          rs <- rbind(rs, r[['RS']])
    }
    
  } #for
  
  return(list('CS'=cs, 'RS'=rs))          
}




##########################################################################################
# Support functions for displaying cluster sets and clusters 
##########################################################################################


#' Displays on console a single cluster set
#' 
#' @details Displays for the cluster set the named fields N, dim and for each cluster its index/id,
#'          n, SUM and SUMSQ vectors. Total number of data inside all clusters in cluster set is
#'          also displayed.
#'            
#' @param cs a CLUSTER SET
#'               
#' @return 
#'
show <- function(cs){
  if (is.null(cs)){
    cat('Number of clusters: 0\n')
    cat('Data dimension: 0\n')
    cat('**null cluster set**\n')
    return()
  } 
  
  totalPoints <- 0
  cat('Number of clusters:', cs[['N']], '\n')
  cat('Data dimension:', cs[['dim']], '\n')
  for (i in 1:cs[['N']]){
    cat('> Cluster ', i, '\n', sep='')
    cat('\tn=', cs[[as.character(i)]]$n, '\n', sep='')
    totalPoints <- totalPoints + cs[[as.character(i)]]$n
    cat('\tSUM=', cs[[as.character(i)]]$SUM, '\n', sep=' ')
    cat('\tSUMSQ=', cs[[as.character(i)]]$SUMSQ, '\n', sep=' ')
  }
  cat('Total points in all clusters in cluster set:', totalPoints, '\n')
}


#' Returns the total number of data in all clusters of a cluster set. 
#' 
#' @details Goes over all clusters in cluster set and adds up the number of data/items in each
#'          cluster.
#'            
#' @param cst a CLUSTER SET
#'               
#' @return integer representing the total number of data/items in cluster set (all clusters)
#'
ClusterSetTotalPoints <- function(cst){
   if (is.null(cst)){
       return(0)   
   }
  
   tot <- 0
   for (i in 1:cst[['N']]){
        tot <- tot + cst[[as.character(i)]]$n
   }   
   
   return(tot)
}





#' Displays on console for each cluster resulting from K-means, the number of data items 
#' belonging to that cluster. 
#' 
#' @details 
#'            
#' @param cResult an object of class kmeans that is the returned value of the K-means 
#'                algorithm
#' @param df the data frame clustered with K-means 
#'               
#' @return 
#'
showClusterDistribution <- function(cResult, df){
    for (i in 1:length(unique(cResult$cluster))){
         pnts <- df[ which(cResult$cluster==i), ]
         cat('\tCluster [', i, '] number of points ', nrow(pnts), '\n', sep='')
    }
}



#' Displays on console information about the cluster sets and data frame in the result set returned
#' by functions of the BFR implementation.  
#' 
#' @details Displays data of the cluster sets DS, CS and the data frame RS found in the
#'          list result set.
#'            
#' @param resultSet a LIST containing the named items DS, CS and RS standing for Discard Set,
#'                  Compression Set and Retained Set respectively. DS and CS are cluster sets, RS
#'                  is a data frame.
#'               
#' @return 
#'
display <- function(resultSet){
    cat('###### Discard set:\n')
    show(resultSet[['DS']])
    cat('\n###### Compression set:\n')
    show(resultSet[['CS']])
    cat('\n###### Retained set:\n')
    print(nrow(resultSet[['RS']]))
    
    if (is.null(resultSet[['RS']]))
        nrs <- 0
    else
        nrs <- nrow(resultSet[['RS']])
    
    
    cat('\n[Total number of points in ALL cluster sets:', ClusterSetTotalPoints(resultSet[['DS']]) +
                                           ClusterSetTotalPoints(resultSet[['CS']]) +
                                           nrs, ']\n') 
}




#' Returns the cluster centroids of all clusters in a cluster set.  
#' 
#' @details The cluster centroid is calculated using the SUM and n fields of the cluster.
#'            
#' @param cst a LIST representing a cluster set. 
#'               
#' @return a Nxdim MATRIX with each row representing the centroids of clusters. Row i of 
#'         resulting matrix corresponds to centroid of cluster with index i in cluster
#'         set.
#'
ClusterCenters <- function(cst){
  
        m <- matrix(nrow=cst$N, ncol=cst$dim)
        for (i in 1:cst$N){
             #print(cst[[as.character(i)]]$SUM/cst[[as.character(i)]]$n)
             m[i,] <- cst[[as.character(i)]]$SUM/cst[[as.character(i)]]$n
        }   
        
        return(m)
}







##########################################################################################
# Intialisation and batch processing functions  
##########################################################################################



#' Processes a batch of the source file containing the data. Returns an LIST containing
#' the updated Discard and Compression Sets (DS, CS cluster sets) and Retained Set (RS), a data frame.
#' 
#' @details Does the following in order: 
#'           I) adds data of the batch to the sufficiently close cluster, 
#'          II) clusters using K-means the remaining data generating CS and updated RS and 
#'         III) merges the resulting CS clusters. 
#'         This function is executed for every new batch read from the source file; except the 
#'         first batch.
#'          
#'            
#' @param batch a DATA FRAME representing a batch of the source data to be clustered 
#' @param ds a LIST representing the current Discard Set
#' @param cs a LIST representing the current Compression Set
#' @param rs a DATA FRAME representing the data that is not into any  Discard or Compression Set
#' @param K integer representing the number of clusters to generate
#' 
#' @return a named LIST (result set) with named items representing the updated DS, CS and RS.
#'
processBatch <- function(batch, ds, cs, rs, K=NUMBER_OF_CLUSTERS){
  res <- list('DS'=ds, 'CS'=cs, 'RS'=rs)
  
  # Step 3: Find those points that are “sufficiently close” to a cluster centroid; add those
  #         points to that cluster and the DS.
  
  cat('Step 3: Adding sufficiently close points to clusters ')
  
  #res <- lapply(batch, AddToClosestCluster, resultSet=res)
  
  for (r in 1:nrow(batch)){
    
    #cat('\tprocess Batch: FOR LOOP RS column names (before):', colnames(rs), '\n') 
    if (r%%100==0){
        cat(r, ' ')
    }
    res <- AddDataToClosestCluster(as.numeric(batch[r, ]), res[['DS']], res[['CS']], res[['RS']])
    #cat('\tprocess Batch: FOR LOOP RS column names (after):', colnames(res[['RS']]), '\n') 
  } # for

  cat('\n')  
  
  # Step 4: Use any main-memory clustering algorithm to cluster the remaining points and the old RS.
  #         Clusters go to the CS; outlying points to the RS.
  
  cat('Step 4: Using main-memory clustering algorithm to cluster remaining points and old RS (K=', 3*K, ')\n', sep='')
  
  # If CS is empty
  if (is.null(res[['CS']])){
    cat('\tEMPTY CS set. Calling MergeDFIntoClusters RS data size: ', nrow(res[['RS']]), '\n')
    if (nrow(res[['RS']]) > 0){
       # Merge RS into CS and RS
       
       if (!'F1' %in% colnames(res[['RS']])){
           print(colnames(res[['RS']]))
           stop('INCORECT HEADER')
       }
       
       
       rsRes <- MergeDFIntoClusters(res[['RS']], 3*K)
       cat('\tGenerated', rsRes[['CS']]$N, ' CS clusters and added ', nrow(rsRes[['RS']]), 'to RS\n')
       res <- list('DS'=res[['DS']], 'CS'=rsRes[['CS']], 'RS'=rsRes[['RS']]) 
    }
  } 
  
  
  # Step 5: Consider merging compression sets in the CS.
  cat('Step 5: Merging CS clusters\n')
  
  # Merge CS
  nCS <- MergeCSCluster(res[['CS']])
  res[['CS']] <- nCS
  
  return(list('DS'=res[['DS']], 'CS'=res[['CS']], 'RS'=res[['RS']]))
}




#' Initializes the BFT algorithm by processing the first batch read from the source file.
#' 
#' @details Does the following:
#'          I) does a first clustering with K-means to separate inliers from outliers. This
#'             separation is based on cluster size (in terms of data)
#'         II) clusters with K-means the inliers of step I) above. The resulting clusters will
#'             constitute the Discard set (DS)
#'        III) Outliers of step I) are clustered using K-means (wit a greater value of K) in order
#'             to generate the Compression and Retained sets. Clusters with more than one item will
#'             constitute a cluster in the Compression set. Data in clusters with only one item will
#'             be added to the Retained set.                      
#'            
#' @param initBatch a DATA FRAME representing first batch of the source data to be clustered 
#' @param K integer representing the number of clusters to generate
#' 
#' @return a named LIST (result set) with named items representing the initialized Discard (DS), Compression 
#'         (CS) and Retained (RS) Sets.
#'
initializeBFR <- function(initBatch, K=NUMBER_OF_CLUSTERS){
  
  # Get a fraction of batch.
  initRows <-  sample(1:nrow(initBatch), nrow(initBatch)*INIT_BATCH_PCT)
  batch <- initBatch[initRows, ]
  remainingBatch <- initBatch[-initRows, ]
  
  
  outliers <- batch[0, ]
  inliers <- batch[0, ]
  ds <- NULL # summarized
  cs <- NULL # summarized
  rs <- batch[0, ] # not summarized. Simply observations as a data frame
  
  if (!'F1' %in% colnames(rs)){
      print(colnames(rs))
      stop('initializeBFR: inconsistend header')
  }
  
  
  # Initialization: Make a first attempt to cluster the batch using K-means. This will
  # separate clusters into outliers and inliers.
  # This is to make sure that the initial set of K clusters have enough data. 
  # Clusters with more than MINIMUM_CLUSTER_POINTS points will be considered for
  # initial generation of Discard Set. Clusters with less than MINIMUM_CLUSTER_POINTS 
  # points will be added to to the Retained Set.
  
  cat(sprintf("Step 0: Initialization of initial batch. Batch size:%d", nrow(batch)), '\n')
  dataClusters <- kmeans(batch, iter.max=MAX_ITERATIONS, centers=as.matrix(batch[sample(nrow(batch), K), ]))
  for (i in 1:K){
    clusteredItems <- batch[which(dataClusters$cluster==i), ]
    #cat('Cluster with ', nrow(clusteredItems), 'items. Limit:', MINIMUM_CLUSTER_POINTS, '\n')
    # Step 2:
    if (nrow(clusteredItems) < MINIMUM_CLUSTER_POINTS){
      cat('\tAdding ', nrow(clusteredItems), ' rows to ourliers\n')
      outliers <- rbind(outliers, clusteredItems)
    }else{
      inliers <- rbind(inliers, clusteredItems)
    }
  } #for   
  
  
  # Cluster inliers.
  cat(sprintf("\tClustering INLIERS. Data size:%d", nrow(inliers)), '\n')
  if (K >= nrow(inliers)){
    message(sprintf('[BFR] Cannot generate initial DS sets. Inliers too few: expected > K=%d, found:%d', K, nrow(inliers)) )
    stop('Terminating.')
  }
  
  
  # Step 1: Take a small random sample and cluster optimally. Here the inliers resulting from above 
  #        are clustered. The resulting clusters will generate the Discard Set.
  
  cat('Step 1: Clustering optimally random sample (inliers)\n')
  inliersClusters <- kmeans(inliers, iter.max=MAX_ITERATIONS, centers=as.matrix(inliers[sample(nrow(inliers), K), ]))
  for (i in 1:K){
    iCluster <- inliers[which(inliersClusters$cluster==i), ]
    sumSq <- colSums(iCluster^2)
    ds <- AddClusterByCentroid(ds, nrow(iCluster), 
                               nrow(iCluster)*as.vector(inliersClusters$centers[i, ]), sumSq)
    
    #write.table(iCluster, file = paste0(i,'.csv'), append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(paste0(i,'.csv')))
    ..DBGAddDSCluster(iCluster)
  } # for
  
  
  
  # Step 2: Use any main-memory clustering algorithm to cluster the remaining points and the old RS.
  #         Clusters go to the CS; outlying points to the RS.
  #         Outliers, the first rows of the RS, are clustered using K-means. Based on cluster criteria,
  #         this clustering will generate the CS (if number of rows > 1) and the other clusters are added
  #         to the RS.
  
  cat('Step 2: Clustering remaining points (initial outliers)\n')
  
  if (K >= nrow(outliers) ){
    # if clustering is not possible, add rows to the RS set
    if (nrow(outliers) > 0){
        #cat('before rbind clusters>outliers', colnames(rs), '\n')
        rs <- rbind(rs, outliers)
        #cat('after', colnames(rs), '\n')
    }
  } else {
    outliersClusters <- kmeans(outliers, iter.max=MAX_ITERATIONS, centers=K)
    
    #showClusterDistribution(outliersClusters, outliers)
    # Outliers will generate the CS and RS
    for (i in 1:length(unique(outliersClusters$cluster))){
      
      oCluster <- outliers[which(outliersClusters$cluster==i), ]
      if (nrow(oCluster) > 1){
        sumSq <- colSums(oCluster^2)
        cs <- AddClusterByCentroid(cs, nrow(oCluster), 
                                   nrow(oCluster)*as.vector(outliersClusters$centers[i, ]),
                                   sumSq)
        ..DBGAddCSCluster(oCluster)
      } else{
        # Add to Retained Set (RS)
        cat('before rbind # data in cluster ==1', colnames(rs), '\n')
        rs <- rbind(rs, oCluster)
        cat('AFTER rbind # data in cluster ==1', colnames(rs), '\n')
      } # nrow(oCluster)
      
    } # for
    
  } # else
  
  
  # Add remaining rows to the closest cluster or RS.
  if (nrow(remainingBatch) > 0){
    cat('Doing remaining batch: size=', nrow(remainingBatch))
    for (i in 1:nrow(remainingBatch)){
      res <- AddDataToClosestCluster(as.numeric(remainingBatch[i,]), ds, cs, rs)
    } #for
    
  } # if
  
  return(list('DS'=ds, 'CS'=cs, 'RS'=rs))
  
} # initialize BFR





##########################################################################################
# Main function to call to execute BFR
##########################################################################################



#' Starts execution of the BFR algorithm. Implements the main loop for processing the 
#' chunks of data. Data is processed in chunks. Chunks of source data to be clustered 
#' are returned by a callback function.
#' 
#' @details BFR calls the function specified by the callback argument and processes
#'          the data accordingly based on the BFR algorithm. The aim of the callback function is 
#'          to model the chunk-based reading and processing of the data. The function specified
#'          in the callback argument may read a source file in chunks or cut a main memory
#'          
#'                      
#' @param callback a FUNCTION that returns a data frame representing the next chunk of source data
#'                 that should be clustered.The function may take as argument the chunk size. The function 
#'                 must return NULL when no more chunks are available.  
#' @param K integer representing the number of clusters to generate
#' 
#' @return a named LIST (result set) with named items representing the initialized Discard (DS), Compression 
#'         (CS) and Retained (RS) Sets. The returned value is the final clustering.
#'
BFR <- function(callback, K=NUMBER_OF_CLUSTERS){
  
    nBatches <- 0
    totalRows <- 0
    res <- list('DS'=NULL, 'CS'=NULL, 'RS'=NULL)
    while(TRUE){
          batch <- callback() # TODO: add arguments to callback to avoid globals?
          if (is.null(batch)){
              cat('........NULL batch (total:', totalRows, ')\n', sep='')
            
              # Step 6: If this is the last round, merge all compression sets in the 
              #         CS and all RS points into their nearest cluster
            
              cat('Step 6: If this is the last round, merge all CS and all RS into nearest cluster\n')
              cat('\tMerging CS into DS (', (if (!is.null(res[['CS']])) res[['CS']]$N else '0'), ' clusters)...\n', sep='')
              res[['DS']] <- MergeClusterSets(res[['CS']], res[['DS']] )
              # Not needed anymore
              res[['CS']] <- NULL 
              
              cat('\tMerging RS into DS (RS rows:', nrow(res[['RS']]), ')...\n', sep='')
              
              if (!is.null(res[['RS']]) && nrow(res[['RS']]) > 0 ){
                  res[['DS']] <- MergeRSToDS(res[['RS']], res[['DS']])
                  res[['RS']] <- NULL
              }
              break
          }
          cat('........New batch (batch size:', nrow(batch), ', total:', totalRows+nrow(batch), ')\n', sep='')
          
          nBatches <- nBatches + 1
          totalRows <- totalRows + nrow(batch)
          if (nBatches == 1){
              res <- initializeBFR(batch, K)
          }else{
             res <- processBatch(batch, res[['DS']], res[['CS']], res[['RS']], K)
          }
    } # while
    
    return(res)
} # BFR


