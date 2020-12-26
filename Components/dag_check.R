dag_check <- function(path){
  
  n_row <- nrow(path)
  n_col <- ncol(path) 
  q <- matrix(c(0), nrow=n_row, ncol=n_col) 
 
  check <- path
  
  N <- n_row - 1

  for (x in 1:N) {
    check <- check %*% path
  }

  for (row in 1:n_row) {
    for (col in 1:n_col) {
      if(check[row,col] == q[row,col]){
      }else{
        return(0)
        next
      }
    }
    
  }
  return(1)
}