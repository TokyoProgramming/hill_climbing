pas<-function(A,D=100){
  
  source('others.r')
  

  col_num <- ncol(A)
  row_num <- col_num
  

  Path    <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path[1,2] = 1
  Path[1,4] = 1
  Path[2,4] = 1
  Path[3,4] = 1
  
  Path_taikaku     <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path_taikaku[1,2] = 1 
  

  
  
  ########################################################
  #################################     Path解析　
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
  #################################     Path_taikaku解析　
  ########################################################
  
  
  
  
 
  
  if ((ncol(A)!=ncol(Path_taikaku))||(ncol(A)!=nrow(Path_taikaku))) stop("Path_taikaku matrix error!")
  colnames_taikaku(Path_taikaku)<-colnames_taikaku(A)
  rownames_taikaku(Path_taikaku)<-colnames_taikaku(A)
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
  

  cat("--------------------------------------------\n")
  ########################################################
  #################################     BIC 追加　＆　配置変更
  ########################################################
  
  

  N_taikaku<-n*(n+1)/2-(k_taikaku+n+sum(Path_taikaku))
  
  N<-n*(n+1)/2-(k+n+sum(Path))
  cat("AIC	=",ki-2*N,"	BIC	=",ki-N*log(nd),"\n")
  cat("AIC_taikaku	=",ki_taikaku-2*N_taikaku,"	BIC_taikaku	=",ki_taikaku-N_taikaku*log(nd),"\n")

  cat("--------------------------------------------\n")
  
  return(Path)
  
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