# This function can be used to write the following files:
# - Pearson_listwise.csv - Listwise Pearson covariances
# - Pearson_pairwise.csv - Pairwise Pearson covariances
# - Spearman_listwise.csv - Listwise Spearman covariances
# - descriptives.txt - Some general descriptive measures
# If there is a column that contains "ID" in its name, it is assumed this is leader ID, 
# and the above mentioned files will be written for samples taken such that
# each leader is rated only by one participant.

# simply run the function below to read it into your environment, and then apply the function to your dataset (see at the bottom)
# please make sure that the MLQ5x items are ordered and/or labeled in the same way as they are in the MLQ5x official distribution.

getCovariances <- function (data,nochecks=F,index=NULL) {
  # check if directory exists and if not, create it using dir.create(paste0(getwd(),"/export"),showWarnings = FALSE)
  if (!dir.exists(paste0(getwd(),"/export"))) {
    dir.create(paste0(getwd(),"/export"),showWarnings = FALSE)
  }
  # compute sample sizes
  nSample_listwise <- sum(apply(data, 1, function (x) all(!is.na(x))))
  nSample_full <- sum(apply(data, 1, function (x) any(!is.na(x))))
  
  # Compute average sample size for pairwise correlations:
  nomisdata <- !is.na(as.matrix(data))
  nMat <- t (nomisdata) %*% nomisdata
  nSample_pairwise <- mean(nMat[lower.tri(nMat, diag=FALSE)],na.omit=TRUE)
  if (!nochecks) {
  # Check if this is a matrix or data frame:
  if (!is.matrix(data) && !is.data.frame(data) ) {
    stop("Input is not a matrix or data frame.")
  }
  # If it is a matrix, make it a data frame:
  if (is.matrix(data) ) {
      data <- as.data.frame(data)
  }
    
  # check if it contains a leader ID column
  if (!any(grepl("ID",colnames(data)))) {
    cat("No leader ID column found. Continuing with the assumption that each participant rated a different leader")
  } else {
    # extract the leader ID column
    leader.grouped <- data
    data <- data[,!grepl("ID",colnames(data))]
  }
  # Descriptives:
  
  # Means:
  means <- colMeans(data,na.rm = TRUE)
  # SDs:
  SDs <- sapply(data, sd, na.rm = TRUE)
  # Number of levels:
  nLevels <- sapply(data,function(x) length(unique(x)[!is.na(unique(x))]))
  
# Set names:
if (is.null(colnames (data) )) {
  cat("No column labels found, assuming the MLQ5x items are ordered from 1 to 36, starting with \n
      'Provides me with assistance in exchange for my efforts...'\n
      if this is not the case, please label the columns accordingly, or contact us")
  colnames (data) <- paste0("MLQ", seq_len(ncol(data)))
}
  # Write descriptives to a file:
  descriptivesFile <- paste0(getwd(),"/export/descriptives.txt")
  write(paste0(
    "Sample size (full): ", nSample_full, "\n",
    "Sample size (listwise): ", nSample_listwise, "\n",
    "Sample size (pairwise average): ", nSample_pairwise,"\n",
    "Name: ", paste0(colnames(data), collapse = "; "), "\n",
    "Means:", paste0(means, collapse ="; "), "\n",
    "Standard deviations: ", paste0(SDs, collapse ="; "),"\n",
    "Number of levels: ", paste0(nLevels, collapse = "; ")),
    file = descriptivesFile)
# Covariances:
  }
  try({
      pearsonCovsFile_listwise <- paste0(getwd(),"/export/Pearson_listwise","_n",nSample_listwise,"i",index,".csv")
      write.csv(cov(data, use = "complete.obs"), file = pearsonCovsFile_listwise)
  })
  try({
    pearsonCovsFile_pairwise <- paste0(getwd(),"/export/Pearson_pairwise","_n",nSample_pairwise,"i",index,".csv")
    write.csv(cov(data, use = "pairwise.complete.obs"), file = pearsonCovsFile_pairwise)
  })
  try({
    spearmanCovsFile_listwise <- paste0(getwd(),"/export/Spearman_listwise","_n",nSample_listwise,"i",index,".csv")
    write.csv(cov(data, use = "complete.obs", method ="spearman"), file = spearmanCovsFile_listwise)
  })
  
  # check if object leader.grouped exists
  if (exists("leader.grouped")) {
    cat("\ngrouping ID detected, computing covariances for randomly subsetted unique leader raters\n")
    index <- 1 
    
    repeat {
      # shuffle the order of responses in dataframe randomly
      leader.grouped <- leader.grouped[sample(nrow(leader.grouped)),]
      # remove all cases that have the same leader ID 
      leader.grouped.unique <- leader.grouped[!duplicated(leader.grouped[,grepl("ID",colnames(leader.grouped))]),]
      # save the removed cases to sample from later
      leader.grouped <- leader.grouped[duplicated(leader.grouped[,grepl("ID",colnames(leader.grouped))]),]
      # remove the leader ID column
      leader.grouped.unique <- leader.grouped.unique[,!grepl("ID",colnames(leader.grouped.unique))]
      
      
      # repeat this as long as there are enough cases 
      sample.obs <- nrow(leader.grouped.unique)
      if(sample.obs <= 36){
        break
      }
      getCovariances(leader.grouped.unique,nochecks = T,index=index)
      index <- index + 1
    }
    cat("\n\nComputed successfully! Please mail us all the files in the folder:\n",
        paste0(getwd(),"/export/"),
        "\n\nThank you for your assistance!","\n In case of any problems, do not hesitate to contact us!")
  }
  if (!nochecks && !exists("leader.grouped")) {
  cat("\n\nComputed successfully! Please mail us the following files:\n\n1. ",
      pearsonCovsFile_listwise,
      "\n2. ",pearsonCovsFile_pairwise,
      "\n3. ", spearmanCovsFile_listwise,
      "\n4. ",descriptivesFile, 
      "\n\nThank you for your assistance!","\n In case of any problems, do not hesitate to contact us!")
  }
}


# how to use the function?
# first load your data using the appropriate function depending on the data type 
## eg. my.data <- read.csv(filename)
# then subset only the MLQ5x columns and the ID column of the leader (if applicable) 
## eg. if the mlq5x items start at column 20 and end at column 56, and the leader ID is the column 1, then use
## my.data <- my.data[,c(20:56,1)]
# then run the function on the subsetted data
## getCorrelations(my.data)

# acknowledgements: This script is based on one supplied by 
## Isvoranu, A.-M., Epskamp, S., & Cheung, M. W.-L. (2021). 
## Network models of posttraumatic stress disorder: A meta-analysis. 
## Journal of Abnormal Psychology, 130(8), 841–861. https://doi.org/10.1037/abn0000704



