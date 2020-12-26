bz <- function(x){
  
  
  set_bic = 1000
  col_num <- x
  
  path <- matrix(c(0), nrow = x, ncol = x, byrow=T)
  
  

    
  for(row in 1:col_num) {
    for(col in 1:col_num) {
      if(row==col){
        path[row,col] = 0
      } else {
        path[row,col] = 1
        check_bic = rnorm(1,1000,10)
        if(check_bic < set_bic){
          set_bic = check_bic
          new_path <- matrix(c(0), nrow = x, ncol = x, byrow=T)
          new_path[row,col] = 1
          path[row,col] = 0
        } else {
          path[row,col] = 0
        }
          
          
      }
        
        
  
      }
      
      
        
    }
    


  path = new_path

  create_path <- matrix(c(0), nrow = 1, ncol = 1, byrow=T)
  
  for (row in 1:col_num) {
    for (col in 1:col_num) {
      if(row==col){
        path[row,col]= 0
          
      }else if(path[row,col] == 1){
        path[row,col] = 1
      
        }else{
        
          path[row,col] = 1
        
          check_bic = rnorm(1,1000,10)
        
          if(check_bic < set_bic){
          
            set_bic = check_bic
          
          
            create_path <- matrix(c(0), nrow = x, ncol = x, byrow=T)
          
          
            create_path[row,col] = 1
          
          
            path[row,col] = 0
        
          
            } else {
          
          
              path[row,col] = 0
        
      
                  }      
       
      
             }    
  
          }
    

  }
  
  
  if(nrow(create_path) > 2){
    p <- new_path + create_path
    
  }else{
    p <- new_path
  }
  
  

  print(p)
  
}