b_zouka_pas<-function(A,D=100, set_bic = 1000){
  
  
  
  kaiki<-function(A){
    ########################################################
    ##################################  回帰係数を求める
    ########################################################
    
    Ax   <- ncol(A)
    Ay   <- nrow(A)
    
    a<-var(A)
    a<-a[-c(Ax),-c(Ax)]
    b<-var(A)
    b<-b[-c(Ax),Ax]
    return(solve(a)%*%b)
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used only in 'paskei'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  Rr<-function(n,Data,Pathkeis){
    ########################################################
    #############################　　構造上のnの分散を求める
    ########################################################
    z<-t(t(Data)*Pathkeis[,n])
    z<-apply(z,1,sum)
    return(var(z))
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'bvar', 'saiyu'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  col_num <- ncol(data)
  
  new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path<-matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  
  
  
  for (Zz in 1:7) {
    for(row in 1:col_num) {
      for(col in 1:col_num) {
        if(row==col){
          Path[row,col] = 0
        } else if(Path[row,col] == 0) {
          Path[row,col] = 1
          
          
          paskei<-function(A,Path){
            ########################################################
            #################################       パス係数を求める
            ########################################################
            
            Ax   <- ncol(Path)
            Ay   <- nrow(Path)
            #	Solve<- matrix(0,Ax,Ay)
            Solve<- Path
            #	flag <- matrix(0,Ay)
            
            for (x in 1:Ax){     	
              part<-NULL
              flag<-0;	
              for(y in 1:Ay)		
                if(Path[y,x]!=0){
                  part<-cbind(part,A[,y])
                  flag[y]<-1
                }else{
                  flag[y]<-0
                }
              
              if (is.null(part)==0) {
                part<-cbind(part,A[,x])
                kai<-kaiki(part)
                
                j<-1
                for(i in 1:Ay){
                  if(flag[i]==1) {
                    Solve[i,x]<-kai[j]
                    j<-j+1
                  }
                }
              }
            }
            return(Solve)
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used in 'bvar', 'saiyu'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara
            #                   12/22/06  by matsu
          }
          
          BVar<-function(data,pathkeis,zan){
            ########################################################
            #################################    母共分散行列を返す
            ########################################################
            a<-solve(diag(ncol(data))-t(pathkeis))
            return(a %*% zan %*% t(a) )
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used in 'saiyu', 'pas'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara	
          }
          
          saiyu<-function(data,pathkei,d=1,kurikaesi){
            ########################################################
            #################################    最尤推定法
            ########################################################
            D<-0.000000002
            s<-cor(data);ssp<-pathkei
            pathkei<-paskei(data,pathkei)
            zan<-diag(ncol(data))
            for(i in 1:ncol(data)){
              zan[i,i]<-1-Rr(i,data,pathkei)
            }
            DEzan=!(zan==1)
            for(i in 1:ncol(data))
              if(sum(pathkei[,i])==0)
                for(k in 1:ncol(data))
                  if(i != k && sum(pathkei[,k])==0){
                    #					if(VarPath[i,k] == 1){
                    DEzan[i,k]<-1
                    DEzan[k,i]<-1
                    zan[i,k]  <-cor(data[,k],data[,i])
                    zan[k,i]  <-cor(data[,k],data[,i])
                    #				 	}
                  }
            gi  <-matrix(0,ncol(data),ncol(data))
            giz <-matrix(0,ncol(data),ncol(data))
            
            for(j in 1:kurikaesi){
              path<-pathkei
              Dzan<-zan
              sig1<-BVar(data,pathkei,zan)
              
              for(i in 1:ncol(data)){
                for(k in 1:ncol(data)){
                  if(pathkei[i,k] != 0){
                    path[i,k]<-pathkei[i,k]+D
                    sig2<-BVar(data,path,zan)
                    gi[i,k]<-bibun(sig1,sig2,s)/D
                  }
                }
              }	
              
              for(i in 1:ncol(data)){
                Dzan[i,i]<-zan[i,i]+D
                sig2<-BVar(data,pathkei,Dzan)
                giz[i,i]<-bibun(sig1,sig2,s)/D
              }
              
              
              sigms<-solve(sig1)%*%s
              f1<-tr(sigms)-log(det(sigms)); 
              
              ggi<-gi*d
              pathkei<-pathkei-( ( ggi/ (  (abs(ggi)<1)+abs((abs(ggi)>1)*(ggi)))))
              giz<-giz*DEzan
              ggiz<-giz*d
              zan<-zan-( ( ggiz/ (  (abs(ggiz)<1)+abs((abs(ggiz)>1)*(ggiz)))))
              
              sig<-BVar(data,pathkei,zan)
              sigms<-solve(sig)%*%s
              f2<-tr(sigms)-log(det(sigms))
              
              if(f1 <= f2){
                d<-d/2
              }
              if(max(abs(giz),abs(gi))<0.0001 ){
                break
              }
              
            }
            
            return(list(path=pathkei,BVar=BVar(data,pathkei,zan)))
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used only in 'pas'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara	
          }
          
          tr<-function(a){
            ########################################################
            #################################       対角行列の和
            ########################################################
            a<-sum(diag(a))
            return(a)
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used in 'Kai', 'GFI', 'saiyu', 'bibun'.
            #
            #  >>> input <<<
            #    tr(a)
            #      a : Square matrix
            #  >>> output <<<
            #	    The sum of a diagonal matrix.
            #
            # Example
            #   a<-matrix(c(1:3,2:4,4:6),3,3)
            #   tr(a)
            #
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #                                                11/25/05  by Sakakibara
          }
          
          bibun <- function(sig1,sig2,s){
            ########################################################
            #################################     　　　　微分
            ########################################################
            sigs1<-solve(sig1)%*%s
            sigs2<-solve(sig2)%*%s
            return(tr(sigs2)-log(det(sigs2))-(tr(sigs1)-log(det(sigs1))) )
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used only in 'saiyu'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara	
          }
          
          Kai<-function(n,sig,s){
            ########################################################
            #################################     　　　　χ^2
            ########################################################
            return((n-1)*(tr(solve(sig)%*%s)-log(det(solve(sig)))-log(det(s))-ncol(s)))
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used only in 'pas'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara	
          }
          
          Pvalue <- function(a,n){
            ########################################################
            #################################     　　　　Ｐ値
            ########################################################
            if(n==0) return(1)
            return(1-pchisq(a,n))
            
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            # This function is used only in 'pas'.
            #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
            #   		     1/07/05  by Sakakibara	
          }
          
          
          
          
          
          ########################################################
          #################################     パス解析　メイン
          ########################################################
          
          Vark<-0
          
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
          
          
          
          
          
          
          
          N <-n*(n+1)/2-(k+n+sum(Path))
          BIC <- ki-N*log(nd)
          
          
          
          if(BIC < set_bic){
            set_bic <- BIC
            new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
            new_path[row,col] = 1
            Path[row,col] = 0
            
            
            
            
            
          } else {
            
            Path[row,col] = 0
          }
          
        }
        
      }
      
      
    }
    
    Path <- Path + new_path
    new_path <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
    
    
    
  }
  
  
  ########################################################
  #################################       追加
  ########################################################
  
  
  
  
  paskei<-function(A,Path){
    ########################################################
    #################################       パス係数を求める
    ########################################################
    
    Ax   <- ncol(Path)
    Ay   <- nrow(Path)
    #	Solve<- matrix(0,Ax,Ay)
    Solve<- Path
    #	flag <- matrix(0,Ay)
    
    for (x in 1:Ax){     	
      part<-NULL
      flag<-0;	
      for(y in 1:Ay)		
        if(Path[y,x]!=0){
          part<-cbind(part,A[,y])
          flag[y]<-1
        }else{
          flag[y]<-0
        }
      
      if (is.null(part)==0) {
        part<-cbind(part,A[,x])
        kai<-kaiki(part)
        
        j<-1
        for(i in 1:Ay){
          if(flag[i]==1) {
            Solve[i,x]<-kai[j]
            j<-j+1
          }
        }
      }
    }
    return(Solve)
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'bvar', 'saiyu'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara
    #                   12/22/06  by matsu
  }
  
  BVar<-function(data,pathkeis,zan){
    ########################################################
    #################################    母共分散行列を返す
    ########################################################
    a<-solve(diag(ncol(data))-t(pathkeis))
    return(a %*% zan %*% t(a) )
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'saiyu', 'pas'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  saiyu<-function(data,pathkei,d=1,kurikaesi){
    ########################################################
    #################################    最尤推定法
    ########################################################
    D<-0.000000002
    s<-cor(data);ssp<-pathkei
    pathkei<-paskei(data,pathkei)
    zan<-diag(ncol(data))
    for(i in 1:ncol(data)){
      zan[i,i]<-1-Rr(i,data,pathkei)
    }
    DEzan=!(zan==1)
    for(i in 1:ncol(data))
      if(sum(pathkei[,i])==0)
        for(k in 1:ncol(data))
          if(i != k && sum(pathkei[,k])==0){
            #					if(VarPath[i,k] == 1){
            DEzan[i,k]<-1
            DEzan[k,i]<-1
            zan[i,k]  <-cor(data[,k],data[,i])
            zan[k,i]  <-cor(data[,k],data[,i])
            #				 	}
          }
    gi  <-matrix(0,ncol(data),ncol(data))
    giz <-matrix(0,ncol(data),ncol(data))
    
    for(j in 1:kurikaesi){
      path<-pathkei
      Dzan<-zan
      sig1<-BVar(data,pathkei,zan)
      
      for(i in 1:ncol(data)){
        for(k in 1:ncol(data)){
          if(pathkei[i,k] != 0){
            path[i,k]<-pathkei[i,k]+D
            sig2<-BVar(data,path,zan)
            gi[i,k]<-bibun(sig1,sig2,s)/D
          }
        }
      }	
      
      for(i in 1:ncol(data)){
        Dzan[i,i]<-zan[i,i]+D
        sig2<-BVar(data,pathkei,Dzan)
        giz[i,i]<-bibun(sig1,sig2,s)/D
      }
      
      
      sigms<-solve(sig1)%*%s
      f1<-tr(sigms)-log(det(sigms)); 
      
      ggi<-gi*d
      pathkei<-pathkei-( ( ggi/ (  (abs(ggi)<1)+abs((abs(ggi)>1)*(ggi)))))
      giz<-giz*DEzan
      ggiz<-giz*d
      zan<-zan-( ( ggiz/ (  (abs(ggiz)<1)+abs((abs(ggiz)>1)*(ggiz)))))
      
      sig<-BVar(data,pathkei,zan)
      sigms<-solve(sig)%*%s
      f2<-tr(sigms)-log(det(sigms))
      
      if(f1 <= f2){
        d<-d/2
      }
      if(max(abs(giz),abs(gi))<0.0001 ){
        break
      }
      
    }
    
    
    
    return(list(path=pathkei,BVar=BVar(data,pathkei,zan)))
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used only in 'pas'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  tr<-function(a){
    ########################################################
    #################################       対角行列の和
    ########################################################
    a<-sum(diag(a))
    return(a)
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'Kai', 'GFI', 'saiyu', 'bibun'.
    #
    #  >>> input <<<
    #    tr(a)
    #      a : Square matrix
    #  >>> output <<<
    #	    The sum of a diagonal matrix.
    #
    # Example
    #   a<-matrix(c(1:3,2:4,4:6),3,3)
    #   tr(a)
    #
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #                                                11/25/05  by Sakakibara
  }
  
  bibun <- function(sig1,sig2,s){
    ########################################################
    #################################     　　　　微分
    ########################################################
    sigs1<-solve(sig1)%*%s
    sigs2<-solve(sig2)%*%s
    return(tr(sigs2)-log(det(sigs2))-(tr(sigs1)-log(det(sigs1))) )
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used only in 'saiyu'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  Kai<-function(n,sig,s){
    ########################################################
    #################################     　　　　χ^2
    ########################################################
    return((n-1)*(tr(solve(sig)%*%s)-log(det(solve(sig)))-log(det(s))-ncol(s)))
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used only in 'pas'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  
  
  
  
  ########################################################
  #################################     パス解析　メイン
  ########################################################
  
  Vark<-0
  
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
  
  cat("--------------------------------------------\n")
  N<-n*(n+1)/2-(k+n+sum(Path))
  
  
  ########################################################
  #################################     BIC 追加　＆　配置変更
  ########################################################
  
  cat("AIC	=",ki-2*N,"	BIC	=",ki-N*log(nd),"\n")
  
  
  
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
  
  
  
  
  
  
  cat("--------------------------------------------\n")
  N<-n*(n+1)/2-(k+n+sum(Path))
  
  
  ########################################################
  #################################     BIC 追加　＆　配置変更
  ########################################################
  
  cat("AIC	=",ki-2*N,"	BIC	=",ki-N*log(nd),"\n")
  
  
  
  cat("--------------------------------------------\n")
  
  
  
  
  
  #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
  #                                                1/07/05  by Sakakibara
  #                                               12/22/06  by matsu
  #                                               10/23/20  by Yoshioka
}