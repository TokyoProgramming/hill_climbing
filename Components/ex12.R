ex12 <- function(x){
  
  
  sum1 <- function(y){
    y <- y + 1
    
    return(y)
  }
  
  
  sum2 <- function(z){
    
    z <- z + 1
    return(z)
  }
  
  result_1 <- sum1(x)
  result <- sum2(result_1)
  
  return(result)
  
}