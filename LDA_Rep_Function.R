###################################################################################
# "Detection of LDA representative with LDA Detection function"
#
# About: This R script describes the function for detection of LDA representative.  
# 
#
###################################################################################
  

# Defining function LDA Representative
#
# Function: LDA_Rep
#
# Arguments:
# LDAReplList            list of LDA replications
# gamma_threshold        per-document-per-topic probabilities (gamma) threshold
#

LDA_Rep <- function(LDAReplList, gamma_threshold = .25) {
  
  # Logical control of arguments. 

  if(gamma_threshold > 1 | gamma_threshold < 0) stop(   
    "Gamma threshold should be in interval [0,1]."
    ,call. = TRUE
  ) 
  
  if (!is.list(LDAReplList)) stop(
    "LDA replications are not in the form of list. Please, provide list of LDA replications in the form of list. See example in /data/SRL_lda.RData. " 
    ,call. = TRUE
    )
  
  if(!all(sapply(LDAReplList, FUN = class) == "LDA_Gibbs")) stop(
    "Elements of LDA replications are not of class LDA_Gibbs."
    ,call. =  TRUE
  )
    
  if(var(sapply(LDAReplList, FUN = function(x) slot(x,"Dim"))[1,]) > 0) stop(
    "Numbers of documents in replications are not equal."
    ,call. = TRUE
  ) 
  
  if(var(sapply(LDAReplList, FUN = function(x) slot(x,"k"))) > 0) stop(
    "Numbers of topics in replications are not equal."
    ,call. = TRUE
  ) 
  
  # Loading packages
  
  require("topicmodels")
  require("tidyverse")
  require("flextable")
  
  # Extraction of parameters 
  
  n.documents <- LDAReplList[[1]]@Dim[1]
  n.topics <- LDAReplList[[1]]@k
  n.replications <- length(LDAReplList)
  
  # Defining Jaccard function
 
  my.jaccard <- function(x,y) {
    n <- max(x)
    m <- max(y)
    d <- matrix(0, nrow = n, ncol = m)
    for(i in 1:n) {
      for(j in 1:m){
        a <- which(x == i)
        b <- which(y == j)
        d[i,j] <- length(intersect(a,b)) / length(union(a,b))
        
      }
    }
    return(d)
  }
  
  
  # Initialization of temporary matrix for Manhattan distance
  
  temp_man <- matrix(0, 
      nrow=n.documents, 
      ncol=n.topics*n.replications)
  
  for(m in 1:n.replications){
    for(d in 1:n.documents){
      temp_man[d, (m-1)*n.topics + 
                    topicmodels::topics(LDAReplList[[m]])[d]] <- 1
    }}
  
  
  # Manhattan distance between documents based on dominant topics from all replications
  
  d_man <- dist(temp_man, method = "manhattan")
  
  # Clustering
  
  viz_man <- hclust(d_man, method = "ward.D2") 
  
  
  # Initialization of temporary matrix for Euclidean distance
  
  temp_euc <- lapply(LDAReplList, function(x) x@gamma)
  temp_euc <- do.call(cbind, temp_euc)
  
  temp_euc[temp_euc < gamma_threshold] <- 0 
  
  
  # Euclidean distance between documents based on topics where gamma >= gamma_threshold
  
  d_euc <- dist(temp_euc, method = "euclidean")
  
  # Clustering
  
  viz_euc <- hclust(d_euc, method = "ward.D2")

  
  # Initialization of data frame for storing results
  
  result <- data.frame(matrix(0, 
                              nrow = n.replications, 
                              ncol = 3))
  
  colnames(result) <- c("rep", paste0("repl_man_",n.topics, "topics"), 
                        paste0("repl_euc_",n.topics, "topics"))
  
  # Results
  
  for(i in 1:n.replications) {
     
    # For each replication, calculating average Jaccard similarity between clusters of documents based on the main topic and clustering based on all replications
    
    replication_manhatan <- 
      mean(apply(my.jaccard(cutree(viz_man, k = n.topics),
                            topicmodels::topics(LDAReplList[[i]])),1,max)) + 
      mean(apply(my.jaccard(cutree(viz_man, k = n.topics), 
                            topicmodels::topics(LDAReplList[[i]])),2,max))
    
    replication_euc <- 
      mean(apply(my.jaccard(cutree(viz_euc, k = n.topics), 
                            topicmodels::topics(LDAReplList[[i]])),1,max)) + 
      mean(apply(my.jaccard(cutree(viz_euc, k = n.topics), 
                            topicmodels::topics(LDAReplList[[i]])),2,max))
    
    result[i, 1] <- i
    result[i, 2] <- replication_manhatan
    result[i, 3] <- replication_euc
  }
 
  # Selecting replication with the highest Jaccard similarity 
  
  result_man <- result |> dplyr::slice_max(result[,2]) |> dplyr::select(c(1,2))
  result_euc <- result |> dplyr::slice_max(result[,3]) |> dplyr::select(c(1,3))
  
  # Returning results
  
  return(list(n.documents = n.documents,
              n.topics = n.topics,
              n.replications = n.replications,
              viz_man = viz_man,
              result_man = result_man,
              viz_euc = viz_euc,
              result_euc = result_euc))
    
}



