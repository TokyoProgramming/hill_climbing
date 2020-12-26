ex5 <- function(){
  
  q <- 0
  p <- array(c(1), dim = c(3, 3, 3))
  
  for (z in 1:3) {
    

    my_matrix <- p[,,z]
    
    for (row in 1:3) {
      for (col in 1:3) {
        
        if(row==col){
          my_matrix[row,col] = 0
        } else {
          my_matrix[row,col] = 1


      }
      
    }
    
  }
  
  
  
  p[,,z] <- my_matrix
  
  
  }
  
  return(p)
}