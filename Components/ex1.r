ex3 <- function(data){
  
  

  col_num <- ncol(data)
  
  path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)  

  
  for (row in 1:col_num) {
    for (col in 1:col_num) {
       
        path[row,col] = 1        
    
    }
    
  }
  
  return(path)
  
  
}