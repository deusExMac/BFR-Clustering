

PREFIX = 'F'
NUM_DIMENSIONS = 20
NUM_ROWS = 25400

# Prefix of file name to save data frame.
DESTINATION_FILE_NAME_PREFIX <- 'testdata/syntheticData'
  
options(scipen=5)
tryCatch(
  {
    HOME_DIRECTORY <- getSrcDirectory(function(){})[1]
    setwd(getSrcDirectory(function(){})[1])
  },
  error = function(err) {
    message(paste("\nERROR changing directory.\nHINT: Are you executing script via RStudio and not source()?\n>>> Original error follows:"))
    err$message <- paste(err, sep = "", end='')
    # Terminate
    stop(err)
  }
)

#setwd(HOME_DIRECTORY)




# Generate a data frame with random data in the range from 0 to 1.
# Can save the data frame to a csv file.
# NOTE: using write.table and not write.csv since write.csv has issues
generateRandomDataFrameUNIF <- function(nr=NUM_ROWS, nc=NUM_DIMENSIONS){

      nms <- c()
      for (i in 1:nc){
           nms <- c(nms, paste0(PREFIX, as.character(i)))
      }
      
      
      # Initialize
      randomDF <- data.frame(matrix(ncol=length(nms), nrow=0) )
      
      # TODO: change this for loop
      for (i in 1:nr){
        randomDF <- rbind(randomDF, runif(n=nc, min=0, max=1))
        if (i%%1000==0){
            cat(i, ' ')
        }
          
      }
      colnames(randomDF) <- nms
      
      return(randomDF)
}





generateRandomDataFrameNORM <- function(nr=NUM_ROWS, nc=NUM_DIMENSIONS, m=0, sdev=1){
  
  randomDF <- data.frame(matrix( rnorm(nr*nc, mean=m, sd=sdev), nr, nc))
  
  nms <- c()
  for (i in 1:nc){
    nms <- c(nms, paste0(PREFIX, as.character(i)))
  }
  
  colnames(randomDF) <- nms
  return(randomDF)
} 




generateDF<-function(nr, nc, distr='U', saveCSV=FALSE){
    if (distr=='U'){
        print('Generating UNIFORM DISTRIBUTED data')
        df <- generateRandomDataFrameUNIF(nr, nc)
    }else{
        print('Generating NORMAL DISTRIBUTED data')
        df <- generateRandomDataFrameNORM(nr, nc, 0, 1)
    }
    
    if (saveCSV){
        fname <- paste0(DESTINATION_FILE_NAME_PREFIX, distr, as.integer(nr), 'x', as.integer(nc), '.csv')
        write.table(df, file=fname, row.names=FALSE, col.names=TRUE, sep=",")
        cat('\n', nr, ' rows ', nc, ' columns saved to [', fname, ']\n', sep='')
    }
    
    return(df)
}

# Generate data frame
df<-generateDF(nr=20000, nc=25, 'N', TRUE)






