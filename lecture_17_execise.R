


rm(list = ls())

setwd("C:/Users/ASUS/OneDrive - Università degli Studi di Padova/uni/Magistrale/1° anno/2° - Machine learning for bioengineering/R_codes")


# ###
hd <- read.table("HeartDisease.dat" , stringsAsFactors = T )
head(hd)
hd$num <- as.factor(hd$num)
hd$slope <- as.factor(hd$slope)
detach(hd)
attach(hd)
library(caret)
library(gbm)
# summary(hd)


#######

# interaction.depth
# n.trees
# shrinkage

set.seed(1001)
idx <- sample( nrow( hd ) , 200 )
trn <- hd[ idx ,]
tst <- hd[ -idx ,]


gbm1 <- gbm( num ~ . ,
             data = trn ,
             train.fraction = 0.5 ,
             n.tree = 500 , 
             cv.folds = 5 ,
             bag.fraction = 0.8 ,
             shrinkage = 0.1 )

best.iter <- gbm.perf( gbm1 , method = "test" )
print(best.iter)
# red curve : validation error
# black  error : training error

# points( 1:500 , gbm1$valid.error )

l <- 5 * length( seq( from = 0.001 , by = 0.01 , to = 0.11 ) )

res <- matrix( rep( 0 , 3*l ) , ncol = 7 , nrow = l )

c = 1

names(res) <- c( "shr" , "dpt" , "err" , "valerr" , "oob" , "test" , "cv" )

# aggiungere ciclo e colonna per maximal depth
for ( shr in seq( from = 0.001 , by = 0.01 , to = 0.11 ) ) {
  
  for ( dpt in seq( from = 1 , by = 1 , to = 5 ) ) {
    
    gbm1 <- gbm( num ~ . ,
               data = trn ,
               train.fraction = 0.5 ,
               n.tree = 500 , 
               cv.folds = 5 ,
               bag.fraction = 0.8 ,
               interaction.depth = dpt ,
               shrinkage = 0.1 )
    
    error <- min( gbm1$valid.error )
    
    # calcolare n.trees dell'albero con i tre metodi
    
    ntrees.valerr <- which( gbm1$valid.error == min( gbm1$valid.error ) )
    
    best.iter.oob <- gbm.perf( gbm1 , method = "OOB" )
    
    best.iter.test <- gbm.perf( gbm1 , method = "test" )
    
    best.iter.cv <- gbm.perf( gbm1 , method = "cv" )
    
    res[c,1] <- shr
    res[c,2] <- dpt
    res[c,3] <- error
    res[c,4] <- ntrees.valerr
    res[c,5] <- best.iter.oob
    res[c,6] <- best.iter.test
    res[c,7] <- best.iter.cv
    
    c = c+1
  
}

}


print(res)

best <- which( res[,4] == min( res[,4] ) )
best.shr <- res[best,2]
best.ntrees <- res[best,1]
best.dpt <- res[]

# costruisco il modello con il minor numero di alberi
best.gbm <- gbm( num ~ . ,
                 data = trn ,
                 train.fraction = 0.5 ,
                 cv.folds = 5 ,
                 bag.fraction = 0.8 ,
                 n.tree = best.ntrees ,
                 shrinkage = best.shr )

summary(best.gbm)


pred.gbm <- predict( best.gbm,
                     newdata = tst ,
                     n.trees = best.ntrees ,
                     type = "response" )
pred.gbm
pred.gbm.mat <- matrix( pred.gbm , ncol = 5 )
pred.class <- as.factor( max.col( pred.gbm.mat ) - 1 )
# sottraggo 1 perchè le colonne vanno da 1 a 5 mentre le classi da 0 a 4
pred.class


# library(caret)
confusionMatrix( pred.class , tst$num )

confusionMatrix( as.factor( pred.class == 0 ) , as.factor( tst$num == 0 ) )






# predizione con due sole classi


gbm2 <- gbm( num != 0 ~ . ,
             data = trn ,
             train.fraction = 0.5 ,
             n.tree = 500 ,
             cv.folds = 5 ,
             bag.fraction = 0.8 ,
             shrinkage = 0.1 )


gbm3 <- gbm( num == 0 ~ . ,
             data = trn ,
             n.trees = 500 ,
             train.fraction = 0.5 ,
             cv.folds = 5 ,
             bag.fraction = 0.8 ,
             shrinkage = 0.1 ,
             interaction.depth = 1 )

best.iter1 <- gbm.perf( gbm3 , method = "OOB" )
print(best.iter1)

best.iter2 <- gbm.perf( gbm3 , method = "cv" )
print(best.iter2)

best.iter3 <- gbm.perf( gbm3 , method = "test" )
print(best.iter3)



pred.gbm1.1 <- predict( gbm3 ,
                        newdata = tst ,
                        n.trees = best.iter1 ,
                        type = "response" )

confusionMatrix( as.factor( pred.gbm1.1 > 0.5 ) , as.factor( tst$num == 0 ) )

x <- table( pred.gbm1.1 > 0.5 ,
            tst$num == 0 )# Accuracy: 0.81

accuracy <- sum( diag( x ) ) / sum( rowSums( x ) )
TP = x[2,2]/sum(x[,2]) #0.833 Sensitivity
TN = x[1,1]/sum(x[,1]) #0.79 Specificity




