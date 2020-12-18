main <- function(A,D=100, kurikaeshi = 20, show=15, set_bic = 1000, whiteList=0, blackList=0){
  
  
  input.show <- show
  
  if(kurikaeshi < show){
    stop('error: "kurikaeshi" must larger than "show" ') 
  } 
  
  
  # create Data frame 
  source('create_df04.r')
  
  result <- create_df04(A, D=100, kurikaeshi, show,set_bic, whiteList, blackList)
  
  res.df <- result$result_df
  res.path_array <- result$res.path_array
  
  # NA check
  for (p in 1:show) {
    if(is.na(res.df$BIC[p]) == TRUE){
      
      show <- (p - 1) 
      break
    }
    
  }
  
  
  
  # 
  # path_id 
  # 
  
  path_id <- rep(0,show)
  for(row in 1:show){
    path_id[row] = res.df[row,1]
  }
  
  cat("--------------------------------------------\n")
  
  for (x in 1:show) {
    
    cat("No.	=",x," \n")
    cat("---------------\n")
    cat("AIC	=",res.df[x,2],"	BIC	=",res.df[x,3],"\n")
    
    cat("DUP	=", res.df[x,4]," \n")
    cat("\n")
    
    print(res.path_array[,,path_id[x]])
    
    cat("--------------------------------------------\n")
  }
  if(show < input.show){
    
    
    
    cat("Available Path matrices : ", show," \n")
  }
  
  
  
  
  
  
  
}