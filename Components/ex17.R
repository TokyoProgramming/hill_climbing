ex17 <- function(){
  
  

  non_dag =  matrix(c(0,1,0,0,
                      0,0,0,1,
                      0,0,0,1,
                      1,0,0,0), ncol=4,nrow=4, byrow=T)
  
  
  dag = matrix(c(0,0,1,0,
                 0,0,0,1,
                 0,0,0,1,
                 0,0,0,0), nrow= 4, ncol=4, byrow = T)
  
  


  
  p <- dag
  
  p_check <- dag
  
  for (x in 1:4) {
    p_check <- p_check %*% p
  }
  
  
  q <- non_dag
  
  q_check <- non_dag
  
  for (x in 1:3) {
    q_check <- q_check %*% q
  }
  

  

  
  return(q_check)
}