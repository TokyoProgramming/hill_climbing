b_hill_climbing<-function(A,D=100, set_bic = 1000, whiteList=0, blackList=0, silent=F){
  
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

  Path_old     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)

  w_row <- whiteList[1]
  w_col <- whiteList[2]
  Path_old[w_row,w_col] = 1
  
  for (path_order in 1:urm_num) {
    row = check_order[1,1,path_order]
    col = check_order[1,2,path_order]
    Path <- Path_old
    Path_taikaku <- Path_old
    
    if(path_order == urm_num){
      sum_matrix <- rowSums(Path)
      sum <- sum(sum_matrix)
      if(urm_num-1 == sum){
        Path[row,col] = 0
      }
    
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
        if(BIC_taikaku < set_bic){
          set_bic <- BIC_taikaku
          Path_old <- Path_taikaku
        }else{
          Path[row,col] = 0
          Path_old <- Path
        }
        
      }else if(is.na(BIC_taikaku) &&!is.na(BIC)){
        if(BIC < set_bic){
          set_bic <- BIC
          Path_old <- Path
        }else{
          Path[row,col] = 0
          Path_old <- Path
        }
        
      }else if(is.na(BIC) && is.na(BIC_taikaku)){
        print(Path)
        print(Path_taikaku)
        stop('error ') 
        
      }else{
        if(BIC <= BIC_taikaku){
          if(BIC < set_bic){
            set_bic <- BIC
            Path_old <- Path
          }else{
            Path[row,col] = 0
            Path_old <- Path
          }
        }else{
          if(BIC_taikaku < set_bic){
            set_bic <- BIC_taikaku
            Path_old <- Path_taikaku
          }else{
            Path[row,col] = 0
            Path_old <- Path
          }
        }
      }
    
    }
  }

  Path_hc <- Path_old

    ########################################################
    #################################       作成したPathを用いたPath解析
    ########################################################
  
    pas_hc.res <- pas(A,Path_hc,D,silent)
  
    invisible(list(Path_hc=Path_hc,res=pas_hc.res))


  #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
  #
  #  >>> input <<<
  #	b_hill_climbing(A,D,set_bic,whiteList,blackList)
  #          A		: item.category_data matrix.
  #          Path	: This is a procession who expresses the pass diagram that consists of 0 and 1. 
  #	     D		: It is the greatest repetition number of times for method of maximum likelihood estimation.
  #					  default value is 100.
  #          set_bic	: 
  #					  default value is 1000.
  #          whiteList	: 
  #					  default value is 0.
  #          blackList	: 
  #					  default value is 0.
  #          silent	: 
  #					  default value is F.
  #  >>> output <<<
  #	    $Path_hc	: It is the best BIC pass matrix.. 
  #	    $res	: It is the result of pas function. 
  #
  #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
  #                                                1/07/05  by Sakakibara
  #                                               12/22/06  by matsu
  #                                               10/23/20  by Yoshioka
}