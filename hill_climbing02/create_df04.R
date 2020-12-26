create_df04 <- function(A, D=100, kurikaeshi= 20, show=5,set_bic=1000, whiteList=0, blackList=0){
  
  col_num <- ncol(A)
  row_num <- col_num
  
  res.BIC_array <- rep(0, kurikaeshi)
  res.AIC_array <- rep(0, kurikaeshi)
  
  res.path_array <- array(0, dim=c(row_num,col_num,kurikaeshi))
  
  
  # Read b_hill_climbing03s
  source('b_hill_climbing_03s.r')
  
  df <- data.frame(id = as.integer(), AIC = numeric(), BIC =numeric(), dup=as.integer())
  k <- 1


  repeat{
    hill_climbing <- b_hill_climbing(A,D,set_bic, whiteList, blackList,silent=T)
    
    res.BIC <- hill_climbing$res$BIC
    res.AIC <- hill_climbing$res$AIC
    res.path <- hill_climbing$Path_hc
    
    res.path_array[,,k] <- res.path
    res.BIC_array[k] <- res.BIC
    res.AIC_array[k] <- res.AIC
  
    large_BIC <- df$BIC[show]
    
    # Read dag_check02
    source('dag_check02.r')
    if(dag_check(res.path) == 1){
      # If res.path == DAG
      
      if((res.BIC %in% as.double(df$BIC)) == TRUE){
        # there is identical BIC Number in the data frame

        # 
        # 1. check : if there are multiple identical BIC numbers
        # 
        
        df_BIC_v <- as.double(df$BIC)
        # data frame row
        df_row    <- which(res.BIC == df_BIC_v)
        
        # get ID in the df_id
        df_id <- df[df_row, 1]
          
        df_id_len <- length(df_id)

        j <- 1
        for (x in 1:df_id_len){
          
          # Check path => identical => add duplicate number   
          id_location <- df_id[x]
          
          if(identical(res.path,res.path_array[,,id_location]) == TRUE){
            
            # Get row number of id_location in the data frame
            df.id_num <- which(grepl(id_location, df$id))
            
            # add dup (duplicate Number)
            new_dup <- as.double(df$dup[df.id_num]) + 1
            
            # Replace the dup element 
            dup_col <- ncol(df)
            df[df.id_num, dup_col] = new_dup

            j <- j + 1 
          }else{
            # path is not identical 
            
            if(x == df_id_len && j == 1){
              
              # if nrow(df) < show => add data to df  
              if(nrow(df) < show){
                add_row <- nrow(df) + 1
                df[add_row,] = c(k, res.AIC, res.BIC, 0)
                
                
                # sort the df  
                df <- df[order(df$BIC),]

                
              }else{
                
                # if nrow(df) > show =>
                # res.BIC < large_BIC =>
                # df[show, ] = c(,,,,)
                if(res.BIC < large_BIC){
                  df[show, ] = c(k,res.AIC, res.BIC, 0)
                  # sort the df  
                  df <- df[order(df$BIC),]


                }

              }

            }

          }

        }

      }else{
        # there is No-identical BIC Number in the data frame
      
        # if nrow(df) < show => add data to df  
        if(nrow(df) < show){
          add_row <- nrow(df) + 1
          df[add_row,] = c(k, res.AIC, res.BIC, 0)
          
          
          # sort the df  
          df <- df[order(df$BIC),]

          
        }else{
          
          # if nrow(df) > show =>
          # res.BIC < large_BIC =>
          # df[show, ] = c(,,,,)
          if(res.BIC < large_BIC){
            df[show, ] = c(k,res.AIC, res.BIC, 0)
            # sort the df  
            df <- df[order(df$BIC),]

            
          }
          
        }

      }


    }
      
    
    if(k == kurikaeshi){
      break
    } 
    
    k <- k + 1
    
    
    
  }
  
  
  result_df <- df
  invisible(list(result_df=result_df, res.path_array=res.path_array))
}




























