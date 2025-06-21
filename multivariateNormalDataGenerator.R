


###############################################################################################
#
#
#  Generates multivariate normally distributed data (synthetic). Subsets may have different
#  parameters.
#
#  v0.1/mmt/Jun 2024
#
###############################################################################################


library(MASS)

# Prefix of file name to save data frame.
DESTINATION_FILE_NAME_PREFIX <- 'testdata/syntheticData'



# random int between 10 and 51: floor(runif(1, min = 10, max = 51))
# TODO: Finish this. Incomplete...
generateRandomMultivariateNormalDistribution<-function(n=100, 
                                                       d=2, 
                                                       mu=c(5, 6), 
                                                       sigMin=0.3,
                                                       sigMax=3.4){
  
  # Standard deviation (a covariance matrix since we want a multivariate
  # normal distribution). Must be a symmetric matrix  
  s<-matrix(runif(d*d, min=sigMin, max=sigMax), d, d)
  sigma <- s %*% t(s)
  
  return( mvrnorm(n = n, mu = mu, Sigma = sigma) )
}


data <- generateRandomMultivariateNormalDistribution(n=200)
data<-rbind(data, generateRandomMultivariateNormalDistribution(n=200, mu=c(20, 18)) )
data<-rbind(data, generateRandomMultivariateNormalDistribution(n=200, mu=c(50, 30)))
data<-rbind(data, generateRandomMultivariateNormalDistribution(n=200, 
                                                              mu=c(35, 5), 
                                                              sigMin=0.04, 
                                                              sigMax=8.34))

colVector <- rep( c("red",  "blue", "yellow", "orange"), each=200)
plot(data, main = "2D Normal Distribution", xlab = "X", ylab = "Y", col=colVector)

# Save to csv file
dimnames(data)<-list(NULL, c("F1", "F2"))
fname <- paste0(DESTINATION_FILE_NAME_PREFIX, 'N200x2-MultiVariate-4-subsets', '.csv')
write.table(data, file=fname, row.names=FALSE, col.names=TRUE, sep=",")


