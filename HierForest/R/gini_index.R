gini_index <- function(class_sort, splitMat, rMin = FALSE){
  
  len <- length(class_sort)
  class_sort <- as.numeric(class_sort)-1
  
  # Find gini index for each split
  ps1 <- (class_sort %*% splitMat)/(1:len)
  ps2 <- (class_sort %*% (1-splitMat))/((len:1)-1)
  p1 <- (1:len/len)
  giniVec <- (1 - (ps1^2 + (1-ps1)^2)) * (p1) + 
    (1 - (ps2^2 + (1-ps2)^2)) * (1-p1)
  
  # It's often only necessary to get the minimum
  if(rMin) return(min(giniVec, na.rm = TRUE)) else return(giniVec)
}