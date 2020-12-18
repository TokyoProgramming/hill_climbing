dag_check <- function(path){
  
  n_row <- nrow(path)
  n_col <- ncol(path) 
#  q <- matrix(c(0), nrow=n_row, ncol=n_col) 
  
  check <- path
  
  N <- n_row - 1
  
  for (x in 1:N) {
    check <- check %*% path
  }
  
#  if(identical(all.equal(check,q),TRUE) == FALSE){
  if(identical(all.equal(sum(check),0),TRUE) == FALSE){
    return(0)
  
  }else{
    
    return(1)
  }
}