create_df02 <- function(A, D=100, kurikaeshi= 10, show=4,set_bic=1000, 
                        whiteList=0, blackList=0){
  
  
  col_num <- ncol(A)
  row_num <- col_num
  
  res.BIC_array <- rep(0, kurikaeshi)
  res.AIC_array <- rep(0, kurikaeshi)
  
  res.path_array <- array(0, dim=c(row_num,col_num,kurikaeshi))
  
  
  # Read b_hill_climbing03s
  source('b_hill_climbing_03s.r')
  
  df <- data.frame(id = numeric(), AIC = numeric(), BIC =numeric(), dup=numeric())
  k <- 1
  repeat{
    
    hill_climbing <- b_hill_climbing(A,D,set_bic, whiteList, blackList,silent=T)
    
    res.BIC <- hill_climbing$res$BIC
    res.AIC <- hill_climbing$res$AIC
    res.path <- hill_climbing$Path_hc
    
    res.path_array[,,k] <- res.path
    res.BIC_array[k] <- res.BIC
    res.AIC_array[k] <- res.AIC

    # Read dag_check02
    source('dag_check02.r')
    
    if(dag_check(res.path) == 1){
      
      large_BIC <- df$BIC[show]

      
      if(is.na(large_BIC) == TRUE){
        # Data frame's row number is smaller than show Number 
        
        if((res.BIC %in% df$BIC) == TRUE){
          # there is the identical BIC Number in the data frame
          
          # 
          # 1. check : if there are multiple identical BIC numbers
          # 
            df_BIC_v <- as.vector(df$BIC)
            df_id    <- which(res.BIC == df_BIC_v)
            df_id_len <- length(df_id)
            
            for (x in 1:df_id_len) {
              
              id_location <- df_id[x]
              if(identical(res.path,res.path_array[,,df_id]) == True){
                # If path_matrices are identical
                
                # add dup (duplicate Number)
                df$dup[df_id[x]] <- df$dup[df_id[x]] + 1
                
              }
            }
            k <- k + 1 

            }else{
              # Identical BIC BUT Non-identical Path
              
              # Add the data to the df
              last_row <- nrow(df) + 1 
              df[last_row, ] = c(k, res.AIC, res.BIC, 0)
              
              
              
              # sort the df  
              df1 <- df[order(df$BIC),]
              df <- df1
              
              
              # if nrow(df) is larger than show delete the last row
              row_Num <- nrow(df)
              if(row_Num < show){
                k <- k + 1
              }else{
                df1 <- df[c(-(show + 1)), ]
                df <- df1
                k <- k + 1 
                
                
              }

            }
        
        



        }else{
          # there is No-identical BIC Number in the data frame
          
          # Append the data to the df
          last_row <- nrow(df) + 1 
          df[last_row, ] = c(k, res.AIC, res.BIC, 0)
          
          # sort the df  
          df1 <- df[order(df$BIC),]
          df <- df1
          
          # if nrow(df) is larger than show delete the last row
          row_Num <- nrow(df)
          if(row_Num < show){
            k <- k + 1
          }else{
            df1 <- df[c(-(show + 1)), ]
            df <- df1
            k <- k + 1 
            
            
          }

  
      }else{
        # Data frame's row number is larger than show Number 
        
        if(res.BIC < large_BIC ){
          # res.BIC is smaller than large BIC 
          # replace res.BIC and large_BIC 
          # sort
          df[large_BIC, ] = c(k, res.AIC, res.BIC, 0)
          # sort the df  
          df1 <- df[order(df$BIC),]
          df <- df1
          # if nrow(df) is larger than show delete the last row
          row_Num <- nrow(df)
          if(row_Num < show){
            k <- k + 1
          }else{
            df1 <- df[c(-(show + 1)), ]
            df <- df1
            k <- k + 1 
            
            
          }
        }else{
          # res.BIC is larger than large_BIC
          k <- k + 1
          
        }

        
      }
      
}

    

    
    if(k > kurikaeshi){
      break
    } 

  }
  

  return(list(k=k,df=df, large_BIC=large_BIC, res.path_array=res.path_array))
  
  print('Until now No Stack')
  
}