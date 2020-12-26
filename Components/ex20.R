ex20 <- function(){
  
  even_sum <- 0

  n <- 1
  repeat {
    if((n %% 2) == 0){
      even_sum <- even_sum + n

      
      
      n <- n + 1
    }else{
      n <- n + 1
    }
    if (even_sum >= 12) break

  }
  
  mat1 <- matrix(c(1:10),nrow=2,ncol=5)
  mat2 <- matrix(c(1:10),nrow=2,ncol=5)
  if(identical(mat1,mat2) == TRUE){
    return(0)
  }else{
    
  
  return(n)
  }
}



