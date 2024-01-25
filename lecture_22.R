

rm(list = ls())
setwd("C:/Users/ASUS/OneDrive - Università degli Studi di Padova/uni/Magistrale/1° anno/2° - Machine learning for bioengineering/R_codes")


data <- read.table("BreastCancer_complete.dat")
head(data)
data$class <- as.factor(data$class)
detach(data)
attach(data)


library(FNN)
library(caret)


set.seed( 19 )
idx <- sample( nrow( data ) , round( nrow( data )*0.8 ) )
trn <- data[ idx , ]
tst <- data[ -idx , ]

res <- matrix( rep( 0 , 5*20 ) , ncol = 5 , nrow = 20 )

# colnames( trn )[ 10 ] <- c( "class" )
# colnames( tst )[ 10 ] <- c( "class" )

colnames( res ) <- c("k","trn_error","tst_error","trn_error_PCA","tst_error_PCA")
rownames( res ) <- rep("",20)


for ( i in 1:20 ) {
  
  # pred trn w\ trn

  knn.trn <- knn( scale( data[ idx , -c( 1 , 11 ) ] ) ,
                  scale( data[ idx , -c( 1 , 11 ) ] ) ,
                  trn$class ,
                  k = i )
  
  res[ i , 1 ] <- i
  res[ i , 2 ] <- sum( knn.trn != trn$class ) / dim( trn )[ 1 ]
  
  
  # pred tst w\ trn
  
  knn.tst <- knn( scale( data[ idx , -c( 1 , 11 ) ] ) ,
                  scale( data[ -idx , -c( 1 , 11 ) ] ) ,
                  trn$class ,
                  k = i )
  
  res[ i , 3 ] <- sum( knn.tst != tst$class ) / dim( tst )[ 1 ]

}

plot( res[ , 1 ] , res[ , 3 ] , ylim = c( 0 , 0.1 ) , col = "red" , type = "b" )
lines( res[ , 1 ] , res[ , 2 ] , ylim = c( 0 , 0.1 ) , col = "green" , type = "b" )



############################## predizione con pca ##############################


trn.pca <- prcomp( trn[ , 2:10 ] , scale = T , center = T )

trn.PC <- trn.pca$x

tst.PC <- predict( trn.pca , tst[ , 2:10 ] )


for ( i in 1:20 ) {
  
  # pred trn w\ trn
  
  knn.trn <- knn( trn.PC[ , 1:2 ] ,
                  trn.PC[ , 1:2 ] ,
                  trn$class ,
                  k = i )
  
  res[ i , 4 ] <- sum( knn.trn != trn$class ) / dim( trn )[ 1 ]
  
  
  # pred tst w\ trn
  
  knn.tst <- knn( trn.PC[ , 1:2 ] ,
                  tst.PC[ , 1:2 ] ,
                  trn$class ,
                  k = i )
  
  res[ i , 5 ] <- sum( knn.tst != tst$class ) / dim( tst )[ 1 ]
  
}

lines( res[ , 1 ] , res[ , 4 ] , ylim = c( 0 , 0.1 ) , col = "orange" , type = "b" )
lines( res[ , 1 ] , res[ , 5 ] , ylim = c( 0 , 0.1 ) , col = "blue" , type = "b" )


