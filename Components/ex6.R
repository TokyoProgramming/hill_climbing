ex6 <- function(data){
  
  
  col_num <- ncol(data)

  
  set_bic = 1000
  



  new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)

    for (z in 1:3) {
      for (row in 1:col_num) {
        for (col in 1:col_num) {
          
          if(row==col){
            Path[row,col] = 0
            
           
          } else if(Path[row,col] == 0) {
            Path[row,col] = 1
            check_bic = rnorm(1,1000,10)
            if(check_bic < set_bic){
              set_bic <- check_bic
              new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
              new_path[row,col] = 1
              Path[row,col] = 0
              
              
              
              
              
            } else {
              
              Path[row,col] = 0
            }
            
          }
          
        }

        
      }

      Path <- Path + new_path
      
      
      
      new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)


    
  }
  
  return(Path)
}