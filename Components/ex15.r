ex15 <- function(A,D=100, kurikaeshi = 5, show=3, set_bic = 1000, whiteList=0, blackList=0){
  
  
  if(kurikaeshi < show){
    stop('error: "kurikaeshi" must larger than "show" ') 
  } 
  
  source("b_hill_climbing_03.r")
  
  col_num <- ncol(A)
  row_num <- col_num
  


  res.BIC_array <- numeric(kurikaeshi)
  res.AIC_array <- numeric(kurikaeshi)

  res.path_array <- matrix(c(0), nrow = row_num, ncol = col_num, byrow=T)

  for (k in 1:kurikaeshi) {


  
    hill_climbing <- b_hill_climbing_03(A,D,set_bic, whiteList, blackList)
    
    res.BIC <- hill_climbing$pas_hc.res$BIC
    res.AIC <- hill_climbing$pas_hc.res$AIC
    res.path <- hill_climbing$Path_hc
    

    res.path_array <- cbind(res.path_array, res.path)
   
    
    res.BIC_array[k] <- res.BIC
    
    res.AIC_array[k] <- res.AIC

  }
  re <- res.path_array[,-(1:col_num)]
  
  res <- array(re, dim = c(row_num,col_num,kurikaeshi))
  
  ###################################################
  ### Crete DataFrame for BIC 
  ###################################################
  
  
  id <- c(1:kurikaeshi)
  df <- data.frame(id=id, BIC=res.BIC_array,AIC=res.AIC_array)
  res.df <-df[order(df$BIC),]


  

    
  ###################################################
  ### show Path and BIC
  ###################################################
  
  bic_num <- numeric(show)
  for (row in 1:show) {
    
    bic_num[row] <- res.df[row,1] 
    
  }

  
  
  
  
  cat("--------------------------------------------\n")
  
  for (x in 1:show) {
    
  
  cat("AIC	=",res.df[x,2],"	BIC	=",res.df[x,3],"\n")
  print(res[,,bic_num[x]])
  cat("--------------------------------------------\n")
  }
  

  
  
  

  ###################################################
  ### return(invisible(list(res.df=res.df,res=res)))
  ###################################################
  
  

  
  return(list(res.df=res.df, res=res))
}