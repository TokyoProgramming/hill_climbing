b_hill_climbing_02<-function(A,D=100, set_bic = 1000, whiteList=0, blackList=0){
  
  col_num <- ncol(A)
  row_num <- col_num
  urm_num <- ((col_num * col_num) - col_num) /2  

  source('pas_calculation.r')
  
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
  
  
  Path     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path_taikaku     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)

  w_row <- whiteList[1]
  w_col <- whiteList[2]
  Path[w_row,w_col] = 1
  Path_taikaku[w_row,w_col] = 1
  
  
  for (path_order in 1:urm_num) {
    row = check_order[1,1,path_order]
    col = check_order[1,2,path_order]
    
    
    if(row == blackList){
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

      if ((ncol(A)!=ncol(Path))||(ncol(A)!=nrow(Path))) stop("Path matrix error!")
      colnames(Path)<-colnames(A)
      rownames(Path)<-colnames(A)
      saiteki<-saiyu(scale(A),Path,1,D)
      ki<-Kai(nrow(A),saiteki$BVar,cor(A))
      n<-ncol(A)
      nd<-nrow(A)
      k<-0
      
      for(i in 1:ncol(A)){
        if(sum(Path[,i]) == 0){
          for(j in i:ncol(A)){
            k<-k+(i != j && sum(Path[,j]) == 0)
          }
        }
      }

      ########################################################
      ################################## Path_taikaku 
      ########################################################


      if ((ncol(A)!=ncol(Path_taikaku))||(ncol(A)!=nrow(Path_taikaku))) stop("Path_taikaku matrix error!")
      colnames(Path_taikaku)<-colnames(A)
      rownames(Path_taikaku)<-colnames(A)
      saiteki_taikaku<-saiyu(scale(A),Path_taikaku,1,D)
      ki_taikaku<-Kai(nrow(A),saiteki_taikaku$BVar,cor(A))
      
      k_taikaku<-0
      
      for(i in 1:ncol(A)){
        if(sum(Path_taikaku[,i]) == 0){
          for(j in i:ncol(A)){
            k_taikaku<-k_taikaku+(i != j && sum(Path_taikaku[,j]) == 0)
          }
        }
      }

      ########################################################
      ################################## Path AIC & BIC
      ########################################################

      
      N <-n*(n+1)/2-(k+n+sum(Path))
      BIC <- ki-N*log(nd)
      AIC	<- ki-2*N

      ########################################################
      #################################   Path_taikaku AIC & BIC
      ########################################################
      
      N_taikaku <-n*(n+1)/2-(k_taikaku+n+sum(Path_taikaku))
      BIC_taikaku <- ki_taikaku-N_taikaku*log(nd)
      AIC	<- ki-2*N

      ########################################################
      #################################      hill-climbing check add, delete
      ########################################################
      

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
  
  
  ########################################################
  #################################       作成したPathを用いたPath解析
  ########################################################
  
  Path_hc <- Path

  ########################################################
  #################################     パス解析　メイン
  ########################################################

  if ((ncol(A)!=ncol(Path_hc))||(ncol(A)!=nrow(Path_hc))) stop("Path_hc matrix error!")
  colnames(Path_hc)<-colnames(A)
  rownames(Path_hc)<-colnames(A)
  saiteki_hc<-saiyu(scale(A),Path_hc,1,D)
  ki_hc<-Kai(nrow(A),saiteki_hc$BVar,cor(A))
  n<-ncol(A)
  nd<-nrow(A)
  k_hc<-0
  
  for(i in 1:ncol(A)){
    if(sum(Path_hc[,i]) == 0){
      for(j in i:ncol(A)){
        k_hc<-k_hc+(i != j && sum(Path_hc[,j]) == 0)
      }
    }
  }
  
  cat("--------------------------------------------\n")
  N_hc<-n*(n+1)/2-(k_hc+n+sum(Path_hc))

  cat("AIC	=",ki_hc-2*N_hc,"	BIC	=",ki_hc-N_hc*log(nd),"\n")

  cat("--------------------------------------------\n")

  return(Path_hc)


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