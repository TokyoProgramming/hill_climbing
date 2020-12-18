ex4 <- function(x){
  
  q <- 0
  
  for (z in 1:x+1) {
    
    
    if(z %% 2 == 1){
      
    }else{
      
      y <- array(c(1), dim = c(x,x, z))
      q <- q + y[,,z]
      
    }
    
    p <- 0

    
    
  }
  
  
  return(q)

  
  
}