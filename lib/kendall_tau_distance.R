inversionNumber <- function(x){
  mergeSort <- function(x){
    if(length(x) == 1){
      inv <- 0
      #printind(' base case')
    } else {
      n <- length(x)
      n1 <- ceiling(n/2)
      n2 <- n-n1
      y1 <- mergeSort(x[1:n1])
      y2 <- mergeSort(x[n1+1:n2])
      inv <- y1$inversions + y2$inversions
      x1 <- y1$sortedVector
      x2 <- y2$sortedVector
      i1 <- 1
      i2 <- 1
      while(i1+i2 <= n1+n2+1){
        if(i2 > n2 || (i1 <= n1 && x1[i1] <= x2[i2])){ # ***
          x[i1+i2-1] <- x1[i1]
          i1 <- i1 + 1
        } else {
          inv <- inv + n1 + 1 - i1
          x[i1+i2-1] <- x2[i2]
          i2 <- i2 + 1
        }
      }
    }
    return (list(inversions=inv,sortedVector=x))
  }
  r <- mergeSort(x)
  return (r$inversions)
}

kendallTauDistance <- function(x,y){
  return(inversionNumber(order(x)[rank(y)]))
}