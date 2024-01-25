

rm(list = ls())
data("PimaIndiansDiabetes2",package="mlbench")
pima.data <- na.omit(PimaIndiansDiabetes2)

pima.data$diabetes <- as.factor(ifelse(pima.data$diabetes=='pos',1,0))
detach(pima.data)
attach(pima.data)
library(pROC)
library(klaR)
library(e1071)
library(caret)


##################################### ROC #####################################

# for reproducibility: every time I run the code I'll get the same
set.seed(42) # generates random number starting from a seed
idx <- sample(nrow(PimaIndiansDiabetes2),500)
# generates 500 values for the data to use for the training
trn <- PimaIndiansDiabetes2[idx,]
tst <- PimaIndiansDiabetes2[-idx,]

log.roc <- glm(diabetes~pregnant+glucose+mass+age,family=binomial(),data=trn)
test.prob <- predict(log.roc,newdata=tst,type='response')
table(test.prob>0.5,tst$diabetes)

tpr = 53/(53+29)
tnr = 153/(153+29)

test.roc <- roc(tst$diabetes~test.prob,plot=T,print.auc=T,
                print.thres="best",add=F)
# TOP alpha is 0: predict everyone  is diabetic: sens = 1, spec = 0
# BOTTOM alpha is 1: predict no one is diabetic: sens = 0, spec = 1
table(test.prob>0.403,tst$diabetes) # valore di alpha che vedo dal plot
# nella tabella la colonna T F è riferita al prob>0.403: T->diab F->nodiab

TPR = 57/(57+29)
TNR = 147/(147+29)

log.roc.1 <- glm(diabetes~mass,family=binomial(),data=trn)
test.prob.1 <- predict(log.roc.1,newdata=tst,type="response")
test.roc.1 <- roc(tst$diabetes~test.prob.1,plot=T,print.auc=T,
                print.thres="best",add=T,col="pink")
# TOP aplha is 0: predict everyone  is diabetic: sens = 1, spec = 0
# BOTTOM alpha is 1: predict no one is diabetic: snes = 0, spec = 1
table(test.prob>0.403,tst$diabetes)


nb.roc <- NaiveBayes(formula(log.roc),data=trn,na.action=na.omit)
test.prob.nb <- predict(nb.roc,newdata=tst,type="response")$posterior[,2]
test.roc.nb <- roc(tst$diabetes~test.prob.nb,plot=T,print.auc=T,add=T,col="red")


log.roc.tot <- glm(diabetes~.,family =binomial(),data=trn)
test.prob.log.tot <- predict(log.roc.tot,newdata=tst,type="response")
test.roc.log.tot <- roc(tst$diabetes~test.prob.log.tot,plot=T,print.auc=T,
                        add=T,col="green")


nb.roc.tot <- NaiveBayes(formula(log.roc.tot),data=trn,na.action=na.omit)
test.prob.nb.tot <- predict(nb.roc.tot,newdata=tst,type="response")$posterior[,2]
test.roc.nb.tot <- roc(tst$diabetes~test.prob.nb.tot,plot=T,print.auc=T,
                       add=T,col="blue")

# bayes does worse than the one including few ones because of the assumption that
# the covariates are independent, which is probably not the case



############################### REGRESSION TREES ###############################


# pred glu as fun of all other var except diab
library(rpart)
library(tree)
library(rattle)

regtree.1 <- rpart(glucose~.-diabetes,data=pima.data,
                   control=rpart.control(maxdepth=30,cp=0.01))
# control with such default parameters
# fine bc i can modify them
# maxdepth: Set the maximum depth of any node of the final tree, with the root
# node counted as depth 0. Values greater than 30 rpart will give nonsense
# results on 32-bit machines
# cp: complexity parameters: Any split that does not decrease the overall lack
# of fit by a factor of cp is not attempted. The lower, the more complex tree
par(mfrow=c(1,2))
fancyRpartPlot(regtree.1)


regtree.2 <- tree(glucose~.-diabetes,data=pima.data)
plot(regtree.2)
text(regtree.2)





### Boosting
N=20		# number of trees
set.seed(22)
idx = sample(nrow(pima.data),250)
trn = pima.data[idx,]
tst = pima.data[-idx,]

trn$eold=trn$glucose
# I add this new column to the data which is the residuals that are updated every cycle
# at 1st cycle corresponds to y (glucose)

dev=rep(NA,N)     # deviance / MSE
ll <- list(NULL)  # trees
res <- list(NULL) # residuals

for(j in 1:N){
  regT <-  rpart(eold~pregnant+pressure+triceps+insulin+mass+pedigree+age,
                 data =trn, control = rpart.control(maxdepth=1,cp=0.005)) # put in list
  # this is predicting the residuals
  ll[[j]] <- regT
  #	data$enew <- data$eold - predict(regT, newdata = trn)
  #	data$eold <- data$enew
  trn$eold <- res[[j]] <- residuals(regT)
  dev[j] <- norm(res[[j]], type="2")^2/length(res[[j]])
  # implementa la formula da minimizzare (slide 7)
}
plot(dev,type='b')



# verification
# prediction is the sum of all the results of the trees
dev_ver <-rep(NA,N)
pred=0
par(mfrow=c(3,7))
for(j in 1:N){
  pred <- pred + unlist(predict(ll[j], newdata = trn))
  # unlist così diventa vettore
  plot(pred~ trn$glucose, col= trn$diabetes, pch=20)
  abline(0,1)
  # straight line where the points should fall if the predictions
  # were perfect 
  dev_ver[j] <- norm(pred-trn$glucose, type="2")^2/length(pred)
}
plot(dev_ver,type='b') # =dev
# deviance is lowering so predictive power is improving

# predict
dev_pred <-rep(NA,N)
pred=0
par(mfrow=c(3,7))
for(j in 1:N){
  pred <- pred + unlist(predict(ll[j], newdata = tst))
  plot(pred~tst$glucose, col=tst$diabetes, pch=20)
  abline(0,1)
  dev_pred[j] <- norm(pred-tst$glucose, type="2")^2/length(pred)
}

plot(dev,type='b', ylim=c(000,1100))
lines(dev_pred, type='l',col=2)
# we can see the test deviance has a minimum and than does not improve
# anymore, above i am overfitting
# if i change the maxdepth the test prediction worsens bc of overfitting


# LM for comparison
fit <- lm(glucose~pregnant+pressure+triceps+insulin+mass+pedigree+age, data=trn)
summary(fit)
pred_lm <- predict(fit, newdata=tst)
dev_lm <- norm(pred_lm-tst$glucose, type="2")^2/length(pred_lm)
dev_lm
dev_pred


par(mfrow=c(1,3))
plot(pred~pred_lm, col=tst$diabetes)
abline(0,1)
plot(pred~tst$glucose, col=tst$diabetes, pch=20)
abline(0,1)
plot(pred_lm~tst$glucose, col=tst$diabetes, pch=20)
abline(0,1)



################### GENERALIZED GRADIENT BOOSTING REGRESSION ###################

# gbm to create the model
# needs to have 0 or 1 in the var to predict
# gbl.perf to evaluate performance of the model finding best num of 
#boosting iterations

library(gbm)

gbm1.a <- gbm( glucose ~ .-diabetes ,
             data = pima.data ,
             n.trees = 20 ,
             train.fraction = 0.7 ,
             cv.folds = 10 ,
             bag.fraction = 0.8 ,
             shrinkage = 1 )
dev.new()
summary(gbm1.a)
# Distribution not specified, assuming gaussian ...
# capisce che glucose è una var continua e assume una distribuzione gaussiana
# come viene fatto nello stimatore ML

gbm1.b <- gbm( glucose ~ pregnant+pressure+triceps+insulin+mass+pedigree+age ,
               data = pima.data ,
               n.trees = 20 ,
               train.fraction = 0.7 ,
               cv.folds = 10 ,
               bag.fraction = 1 ,
               shrinkage = 1 )
summary(gbm1.b)


# how good is the model changing the number of trees
# can be used with different methods
bestiter.aOOB <- gbm.perf( gbm1.a , method = "OOB" )
print( bestiter.aOOB )

bestiter.atest <- gbm.perf( gbm1.a , method = "test" )
print( bestiter.atest )

bestiter.acv <- gbm.perf( gbm1.a , method = "cv" )
print( bestiter.acv )


bestiter.bOOB <- gbm.perf( gbm1.b , method = "OOB" )
print( bestiter.bOOB )

bestiter.btest <- gbm.perf( gbm1.b , method = "test" )
print( bestiter.btest )

bestiter.bcv <- gbm.perf( gbm1.b , method = "cv" )
print( bestiter.bcv )



### want to predict diabetes
# needs to be converted to 0 1 first

gbm2 <- gbm(diabetes~. , data = pima.data , n.trees = 500 ,
            train.fraction = 0.75 , cv.folds = 5 , bag.fraction = 0.8 , shrinkage = 1 )

best.iter2 <- gbm.perf(gbm2,method="OOB")
print(best.iter2)
# measures how good is the model changing the number of trees


# can be used with different methods

best.iter2 <- gbm.perf(gbm2,method="test")
print(best.iter2)

best.iter2 <- gbm.perf(gbm2,method="cv")
print(best.iter2)


# which one is the best? == tuning boosting algorithm
# split data in to trn and tst, tune on trn and test on tst

set.seed(22)
idx <- sample( nrow( pima.data ) , 250 )
trn <- pima.data[ idx ,]
tst <- pima.data[ -idx ,]


gbm3 <- gbm( as.integer(diabetes) ~ . , data = trn , n.trees = 1000 ,
            train.fraction = 0.5 , cv.folds = 5 , bag.fraction = 0.5 , shrinkage = 0.01 )
plot( gbm3 , i.var = "insulin" )
# lowering shrinkage means we are taking more steps moving down the gradient

best.iter <- gbm.perf(gbm3, method = "cv")
print(best.iter)

pred.gbm <- predict( gbm3 , newdata = tst , n.trees = best.iter , type = "response" )
# if n.trees not specified it uses the one from test
# asking type response gives a probability
pred.gbm

 
x <- table( pred.gbm>0.5 , tst$diabetes )
x

TP = x[2,2]/sum(x[,2]) # specificity 66% non viene
TN = x[1,1]/sum(x[,1]) # sensitivity 89% non viene

library(pROC)
test_roc = roc( tst$diabetes ~ pred.gbm , plot = T , print.auc = T,
                print.thres = "best" , add = F )
# AUC = 0.86


