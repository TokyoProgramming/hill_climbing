b_hill_climbing_03<-function(A,D=100, set_bic = 1000, whiteList=0, blackList=0){
  
  col_num <- ncol(A)
  row_num <- col_num
  urm_num <- ((col_num * col_num) - col_num) /2  

  ########################################################
  #################################     Hill-climb　
  #################################     ランダム経路
  #################################     作成 & 配列格納
  #################################     upper_triangle_matrix_hairetsu.r
  ########################################################
  
  source('upper_triangle_matrix_hairetsu.r')

  check_order <- hairetsu_o(A)

  ########################################################
  #################################     Path & Path_taikaku 比較
  ########################################################
  
  source("pas_add_bic3.r")

  Path     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path_taikaku     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)

  w_row <- whiteList[1]
  w_col <- whiteList[2]
  Path[w_row,w_col] = 1
  Path_taikaku[w_row,w_col] = 1
  
  for (path_order in 1:urm_num) {
    row = check_order[1,1,path_order]
    col = check_order[1,2,path_order]
    
    ########################################################
    #################################     上三角行列がすべて1でないかの確認
    ########################################################
    
    if(path_order == urm_num-1){
      Path_check <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
      Path_check <- Path + Path_taikaku
      sum_matrix <- rowSums(Path_check)
      sum <- sum(sum_matrix)
      if(urm_num == sum){
        Path[row,col] = 0
      }
      
      ########################################################
      #################################    確認 END
      ########################################################

    
    }else if(row == blackList){
      Path[row,col] = 0
    }else if(row == w_row && col == w_col){
      Path[w_row,w_col] = 1
      Path_taikaku <- Path
    
    }else{
      
      if(col == blackList){
        Path[row,col] = 1
        Path_taikaku <- Path
        Path_taikaku[blackList,row] = 0

      }else{
        Path[row,col] = 1
        Path_taikaku[col,row] = 1

      }

      ########################################################
      ################################## Path
      ########################################################

      pas.res <- pas(A,Path,D,silent=T)

      ########################################################
      ################################## Path_taikaku 
      ########################################################

      pas_taikaku.res <- pas(A,Path_taikaku,D,silent=T)

      ########################################################
      #################################      hill-climbing check add, delete
      ########################################################

      BIC <- pas.res$BIC
      BIC_taikaku <- pas_taikaku.res$BIC
      
      if(is.na(BIC) && !is.na(BIC_taikaku)){

        BIC <- BIC_taikaku

        if(BIC <= set_bic){
          set_bic <- BIC
          Path <- Path_taikaku
        }else{
          Path[row,col] = 0
          Path_taikaku <- Path
        }
        
      }else if(is.na(BIC_taikaku) &&!is.na(BIC)){

        BIC_taikaku <- BIC
        if(BIC_taikaku < set_bic){
          set_bic <- BIC_taikaku
          Path_taikaku <- Path
          
        }else{
          Path[row,col] = 0
          Path_taikaku <- Path
        }
        
      }else if(is.na(BIC) && is.na(BIC_taikaku)){

        stop('error ') 
        
        
      }else{
      
        if(BIC <= BIC_taikaku){
          if(BIC <= set_bic){
            set_bic <- BIC
            Path_taikaku <- Path
          }else{
            Path[row,col] = 0
            Path_taikaku <- Path
          }
        }else{
          if(BIC_taikaku < set_bic){
            set_bic <- BIC_taikaku
            Path <- Path_taikaku
            
          }else{
            Path[row,col] = 0
            Path_taikaku <- Path
          }
        }
      }
    }
  }
  
  ########################################################
  #################################       作成したPathを用いたPath解析
  ########################################################
  
  Path_hc <- Path
  pas_hc.res <- pas(A,Path_hc,D,silent=T)
  
  return(list(Path_hc=Path_hc, pas_hc.res=pas_hc.res))
  





  #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
  #
  #  >>> input <<<
  #	pas(A,Path,D)
  #          A		: item.category_data matrix.
  #          Path		: This is a procession who expresses the pass diagram that consists of 0 and 1. 
  #	   D		: It is the greatest repetition number of times for method of maximum likelihood estimation.
  #					  default value is 100.
  #  >>> output <<<
  #	    $path	: It is a matrix of the pass coefficient. 
  #	    $BVar	: It is presumed covariance matrix. 
  #
  # Example
  # (pass diagram)	
  #       1 → 2 →4
  #            ↑
  #            3
  #
  # p=matrix(c(0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0),4,4)
  # pas(data,p)
  # pas(data,p,D=50)
  

  #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
  #                                                1/07/05  by Sakakibara
  #                                               12/22/06  by matsu
  #                                               10/23/20  by Yoshioka
}