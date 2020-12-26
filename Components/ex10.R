ex10<-function(A,D=100, set_bic = 1000){
  
  
  
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
  row_num <- col_num
  
  Path    <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  Path_taikaku <- matrix(c(0), nrow = col_num, ncol = col_num, byrow=T)
  
  row = 5
  col = 6
  
  Path[row,col] = 1
  Path_taikaku[col,row] = 1 

  
  
  
  ########################################################
  ########################################################
  ################################## Path
  ########################################################
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
  
  ########################################################
  ########################################################
  ################################## Path END
  ########################################################
  ########################################################
  
  N <-n*(n+1)/2-(k+n+sum(Path))
  BIC <- ki-N*log(nd)
  AIC	<- ki-2*N
  
  
  ########################################################
  ########################################################
  ################################## Path_taikaku 
  ########################################################
  ########################################################
  

  

  
  
  ########################################################
  #################################   Path_taikaku_taikaku 
  ########################################################
  
  
  paskei_taikaku<-function(A,Path_taikaku){
    ########################################################
    #################################       パス係数を求める
    ########################################################
    
    Ax   <- ncol(Path_taikaku)
    Ay   <- nrow(Path_taikaku)
    #	Solve<- matrix(0,Ax,Ay)
    Solve_taikaku<- Path_taikaku
    #	flag <- matrix(0,Ay)
    
    for (x in 1:Ax){     	
      part<-NULL
      flag<-0;	
      for(y in 1:Ay)		
        if(Path_taikaku[y,x]!=0){
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
            Solve_taikaku[i,x]<-kai[j]
            j<-j+1
          }
        }
      }
    }
    return(Solve_taikaku)
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'bvar', 'saiyu'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara
    #                   12/22/06  by matsu
  }
  
  Bvar_taikaku<-function(data,pathkeis_taikaku,zan_taikaku){
    ########################################################
    #################################    母共分散行列を返す
    ########################################################
    a<-solve(diag(ncol(data))-t(pathkeis_taikaku))
    return(a %*% zan_taikaku %*% t(a) )
    
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    # This function is used in 'saiyu', 'pas'.
    #<<<<<>>>>><<<<<>>>>><<<<<>>>>><<<<<>>>>>
    #   		     1/07/05  by Sakakibara	
  }
  
  saiyu_taikaku<-function(data,pathkei_taikaku,d=1,kurikaesi){
    ########################################################
    #################################    最尤推定法
    ########################################################
    D<-0.000000002
    s<-cor(data);ssp<-pathkei_taikaku
    pathkei_taikaku<-paskei_taikaku(data,pathkei_taikaku)
    zan_taikaku<-diag(ncol(data))
    for(i in 1:ncol(data)){
      zan_taikaku[i,i]<-1-Rr(i,data,pathkei_taikaku)
    }
    DEzan_taikaku=!(zan_taikaku==1)
    for(i in 1:ncol(data))
      if(sum(pathkei_taikaku[,i])==0)
        for(k in 1:ncol(data))
          if(i != k && sum(pathkei_taikaku[,k])==0){
            #					if(VarPath_taikaku[i,k] == 1){
            DEzan_taikaku[i,k]<-1
            DEzan_taikaku[k,i]<-1
            zan_taikaku[i,k]  <-cor(data[,k],data[,i])
            zan_taikaku[k,i]  <-cor(data[,k],data[,i])
            #				 	}
          }
    gi  <-matrix(0,ncol(data),ncol(data))
    giz <-matrix(0,ncol(data),ncol(data))
    
    for(j in 1:kurikaesi){
      path<-pathkei_taikaku
      Dzan_taikaku<-zan_taikaku
      sig1<-Bvar_taikaku(data,pathkei_taikaku,zan_taikaku)
      
      for(i in 1:ncol(data)){
        for(k in 1:ncol(data)){
          if(pathkei_taikaku[i,k] != 0){
            path[i,k]<-pathkei_taikaku[i,k]+D
            sig2<-Bvar_taikaku(data,path,zan_taikaku)
            gi[i,k]<-bibun(sig1,sig2,s)/D
          }
        }
      }	
      
      for(i in 1:ncol(data)){
        Dzan_taikaku[i,i]<-zan_taikaku[i,i]+D
        sig2<-Bvar_taikaku(data,pathkei_taikaku,Dzan_taikaku)
        giz[i,i]<-bibun(sig1,sig2,s)/D
      }
      
      
      sigms<-solve(sig1)%*%s
      f1<-tr(sigms)-log(det(sigms)); 
      
      ggi<-gi*d
      pathkei_taikaku<-pathkei_taikaku-( ( ggi/ (  (abs(ggi)<1)+abs((abs(ggi)>1)*(ggi)))))
      giz<-giz*DEzan_taikaku
      ggiz<-giz*d
      zan_taikaku<-zan_taikaku-( ( ggiz/ (  (abs(ggiz)<1)+abs((abs(ggiz)>1)*(ggiz)))))
      
      sig<-Bvar_taikaku(data,pathkei_taikaku,zan_taikaku)
      sigms<-solve(sig)%*%s
      f2<-tr(sigms)-log(det(sigms))
      
      if(f1 <= f2){
        d<-d/2
      }
      if(max(abs(giz),abs(gi))<0.0001 ){
        break
      }
      
    }
    
    
    
    return(list(path=pathkei_taikaku,Bvar_taikaku=Bvar_taikaku(data,pathkei_taikaku,zan_taikaku)))
    
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
    # This function is used in 'Kai', 'GFI', 'saiyu_taikaku', 'bibun'.
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
    # This function is used only in 'saiyu_taikaku'.
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
  
  if ((ncol(A)!=ncol(Path_taikaku))||(ncol(A)!=nrow(Path_taikaku))) stop("Path_taikaku matrix error!")
  colnames(Path_taikaku)<-colnames(A)
  rownames(Path_taikaku)<-colnames(A)
  saiteki<-saiyu_taikaku(scale(A),Path_taikaku,1,D)
  ki<-Kai(nrow(A),saiteki$Bvar_taikaku,cor(A))
  n<-ncol(A)
  nd<-nrow(A)
  k<-0
  
  for(i in 1:ncol(A)){
    if(sum(Path_taikaku[,i]) == 0){
      for(j in i:ncol(A)){
        k<-k+(i != j && sum(Path_taikaku[,j]) == 0)
      }
    }
  }
  

  #######################################################
  ################################   Path_taikaku_taikaku end
  #######################################################
  
  
  
  
  ########################################################
  #################################   Path_taikaku_BIC
  ########################################################
  
  
  N_taikaku<-n*(n+1)/2-(k+n+sum(Path_taikaku))
  AIC_taikaku <- ki-2*(N_taikaku)
  BIC_taikaku <- ki-(N_taikaku)*log(nd)

  
  ########################################################
  #################################      hill-climbing check add, delete
  ########################################################
  
  
  if(BIC < BIC_taikaku){
    
    
  }else{

    Path <- Path_taikaku
  }
  
  
  
  
  
  
  
  
  cat("--------------------------------------------\n")

  
  
  ########################################################
  #################################     BIC 追加　＆　配置変更
  ########################################################
  cat("AIC	=",AIC,"	BIC=",BIC,"\n")
  cat("AIC_taikaku	=",AIC_taikaku,"	BIC_taikaku=",BIC_taikaku,"\n")

  
  
  
  cat("--------------------------------------------\n")
  

  
  
  return(Path)
  
  

}