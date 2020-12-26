create_df <- function(A, D, kurikaeshi, show,set_bic, whiteList, blackList, k){
  
  
  col_num <- ncol(A)
  row_num <- col_num
  
  res.BIC_array <- rep(0,kurikaeshi)
  res.AIC_array <- rep(0,kurikaeshi)

  res.path_array <- array(0, dim=c(row_num,col_num,kurikaeshi))
  

  # Read b_hill_climbing_03s
  source("b_hill_climbing_03s.r")
  
  hill_climbing <- b_hill_climbing(A,D,set_bic, whiteList, blackList,silent=T)
  
  res.BIC <- hill_climbing$res$BIC
  res.AIC <- hill_climbing$res$AIC
  res.path <- hill_climbing$Path_hc
  
  
  res.path_array[,,k] <- res.path
  res.BIC_array[k] <- res.BIC
  res.AIC_array[k] <- res.AIC
  
  # Read dag_check02
  source('dag_check02.r')
  
  if(dag_check(hill_climbing$Path_hc) == 1){
    # Create Empty data frame
    df <- data.frame(id = numeric(), AIC = numeric(), BIC =numeric(), dup=numeric())
    
    
    # : New_BIC
    # : Last_BIC : BIC in the last row
    large_BIC <- df$BIC[show]
    
    # If the df row number is less than show
    if(is.na(large_BIC) == TRUE){
      
      
      # if BIC is not the unique Number
      # Then check path matrix
      if((res.BIC %in% df$BIC) == FALSE){
  
        if(nrow(df) != 0){
          df_id <- match(res.BIC,df$BIC)
          
          if(identical(res.path,res.path_array[,,df_id]) == True){
            # If path_matrices are identical
            # Then add dup (duplicate Number)
            
            df$dup[df_id] <- df$dup[df_id] + 1
            k <- k + 1 
          
             
          }else{
            # Same BIC value But different Path 
            # If path_matrices are identical
            
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
              df1 <- df[-c(show + 1), ]
              df <- df1
              k <- k + 1 

              
            }
          }
          
        }else{
          # Append the data to the df
          last_row <- nrow(df) + 1 
          df[last_row, ] = c(k, res.AIC, res.BIC, 0)
          
          # sort the df  
          df1 <- df[order(df$BIC),]
          df <- df1
          
          k <- k + 1
        }
  
      # if BIC is the unique Number
      }else{
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
          df1 <- df[-c(show + 1), ]
          df <- df1
          k <- k + 1 

          
        }
  
      }
    
    # If the df row number is larger than show
    # Check the BIC is in the df or not 
    }else{
      # If there is the same BIC Number 
      # Add to dup(duplicate number)
      if((res.BIC %in% df$BIC) == TRUE){
        df_id <- match(res.BIC, df$BIC)
        df$dup[df_id] <- df$dup[df_id] + 1
        k <- k + 1 
      
      # Compare large_BIC with BIC
        
      }else{
        # Append the data to the df
        # Drop the last row data in the df 
        if(res.BIC < large_BIC){

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
            df1 <- df[-c(show + 1), ]
            df <- df1
            k <- k + 1 

            
          }
          
        }else{
          k <- k + 1

        }
        
        }
      
      }
    
  }
  
  return(list(df=df,k=k))
  
}






































