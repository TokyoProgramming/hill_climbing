ex18 <- function(path){
  
  
  source("dag_check.r")
  result <- dag_check(path)
  
  if(result == 1){
    print('This path is dag')
  }else{
    print('This path is not dag')
  }

}