ex21 <- function(num){
  
  df <- data.frame(id = numeric(num), AIC = numeric(num), BIC =numeric(num), dup=numeric(num))
  
  for (i in 1:num) {
    
    id <- 4-i
    AIC <- 4-i
    BIC <- 4-i
    
    if((BIC %in% df$BIC) == FALSE){ 
    
    df$id[i] <- id
    df$AIC[i] <- AIC
    df$BIC[i] <- BIC
    df1 <- df[order(df$BIC),]
    }else{
      path_id <- match(BIC,df$BIC)
      if(identical(res.path,res.path_array[,,path_id]) == True){
        
        
      }
    }
  }

  
  return(df1)
  
}