hairetsu_o <- function(data){


  urm <- function(Data){
    col_num <- ncol(Data)
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
    
    return(mat1)
    
    
    
  }

  
  
  hairetsu <- function(order){
    
    col_num <- ncol(order)
    row_num <- col_num

    urm_num_k <- ((col_num * col_num) - col_num) /2  
    
    
    check_order <- array(0, dim=c(1, 2, urm_num_k))
    
    
    
    i = 1
    for (order_num in 1:urm_num_k) {
      
      
      for (row in 1:row_num) {
        for(col in 1:col_num) {
          
          if(order[row,col] == i){
            
            check_order[1,1,i] = row
            check_order[1,2,i] = col
            
            i = i + 1
            
            if(i==urm_num_k) break
            
            
            
          }
          
          
          
          
        }
        
      }
      
    }
    return(check_order)
    
    
  }
  
  hairetsu_order <- hairetsu(urm(data))

  return(hairetsu_order)  
  
  
  
}