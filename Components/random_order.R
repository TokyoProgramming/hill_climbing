random <- function(data){
  
  col_num <- ncol(data)
  
  
  urm_num <- ((col_num * col_num) - col_num) /2  
    
  check_order <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  
  c_vector <- sequence(urm_num)
  
  random_order <- sample(c_vector)
  
  return(random_order)

  
}
