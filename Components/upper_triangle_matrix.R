urm <- function(data){
  
  col_num <- ncol(data)
  
  urm_num <- ((col_num * col_num) - col_num) /2  
  c_vector <- sequence(urm_num)
  random_order <- sample(c_vector)

  mat1 <- matrix(c(0),nrow=col_num,ncol=col_num, byrow=T)
  i = 0
  row_num <- col_num - 1 
  
  for (row in 1:row_num) {
    k <- row + 1
    
    for(col in k:col_num){
      
      i = i + 1
      
      mat1[row,col] = random_order[i]
      
      
      
      
      
    }
      
  }
  
  return(mat1)


  
}
