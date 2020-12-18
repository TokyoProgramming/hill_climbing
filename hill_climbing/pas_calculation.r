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

paskei<-function(A,k){
  ########################################################
  #################################       パス係数を求める
  ########################################################
  
  Ax   <- ncol(k)
  Ay   <- nrow(k)
  #	Solve<- matrix(0,Ax,Ay)
  Solve<- k
  #	flag <- matrix(0,Ay)
  
  for (x in 1:Ax){     	
    part<-NULL
    flag<-0;	
    for(y in 1:Ay)		
      if(k[y,x]!=0){
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


