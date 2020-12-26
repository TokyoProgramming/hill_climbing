hairetsu <- function(data){
  
  col_num <- ncol(data)
  row_num <- col_num
  
  urm_num <- ((col_num * col_num) - col_num) /2  

  check_order <- array(0, dim=c(1, 2, urm_num))
  

  i = 1
  for (order_num in 1:urm_num) {
    
  
    for (row in 1:row_num) {
      for(col in 1:col_num) {
        
        if(data[row,col] == i){
          
          check_order[1,1,i] = row
          check_order[1,2,i] = col
          
          i = i + 1
          
          if(i==urm_num) break
          
  
          
        }
        
        
        
        
      }
      
    }
   
  }
  return(list(check_order=check_order,data=data))
  
  
}