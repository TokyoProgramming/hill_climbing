ex8<- function(data){
  
  
  col_num <- ncol(data)
  row_num <- col_num
  
  
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
  
  ########################################################
  #################################     ランダム作成 配列割り当て
  ########################################################
  check_order <- array(0, dim=c(1, 2, urm_num))
  
  
  i = 1
  for (order_num in 1:urm_num) {
    
    
    for (row in 1:row_num) {
      for(col in 1:col_num) {
        
        if(mat1[row,col] == i){
          
          check_order[1,1,i] = row
          check_order[1,2,i] = col
          
          i = i + 1
          
          if(i==urm_num) break
          
          
          
        }
        
        
        
        
      }
      
    }
  }
    
  return(check_order)
  
 
}

