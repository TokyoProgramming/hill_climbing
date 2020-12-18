f2 <- function(){
  
  
 x <- matrix(1, nrow = 3, ncol=3, dimnames = list(c("X","Y","Z"), c("A","B","C")))


  y <- x[,1]
  a <- c(3, 2, -7, 3, 5, 2)
  y_1 <- as.vector(y)
  z <- which(y_1 == 1)
  z_1 <- as.vector(z)
  
  
  x <- 1
  repeat{
    
    

    for (i in 1:3) {

      
    }
    x <- x + 1
    if(x == 5){
      break
    }
  }
  
  
  sex    <- c("F","F","M","M","M")
  height <- c(158,162,177,173,166)
  weight <- c(51,55,72,57,64)
  x    <- data.frame(SEX=sex, HEIGHT=height, WEIGHT=weight)
  id <- tail(x,1)

  
  df <- rbind(x,c('M', 1, 1))
  x <- df
  
  return(x)


}