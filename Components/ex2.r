ex2 <- function(){
  
  path <- matrix(c(0), nrow = 3, ncol = 3, byrow=T)

  for (k in 1:3) {
    for (i in 1:3) {
      for (j in 1:3) {


            
            
            if(i==j){
              path[i,j]=1
            
              }else{
              path[i,j] = 0
            
              }
          
            }
      print(paste("Row", i, "and column",j, "have values of", path[i,j]))   
    
          }
    
        
  
      

    }
  
  
}