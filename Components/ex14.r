ex14 <- function(x){
  

  Path     <- matrix(c(1), nrow = x, ncol = x, byrow=T)
  sum_matrix <- rowSums(Path)
  sum <- sum(sum_matrix)
  

  
  return(sum)
  
}