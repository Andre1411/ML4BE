

### EXERCISE
# Plot the data: Ccell vs V2h with B-cells and non-B-cells having different colors and symbols
# Add a legend to the figure
# Build logistic models with V2h, Ccell, kh as covariates
# Use different tools to evaluate the different models, to see which of V2h, Ccell, kh predicts the cell type (B/non-B) "best" 
# Build models that combine the different predictive covariates (V2h, Ccell, kh).
# Is there an optimal model?
#
# If you have time: Add a new column to the dataset indicating if the cell is an A-cell:
# ct$Acell <- ifelse(ct$CellType ... ) # fill in, write ?ifelse to see the help	
# Repeat the exercise above for A-cells


################################### INTRO ######################################

rm(list = ls())
setwd("C:/Users/ASUS/OneDrive - Università degli Studi di Padova/uni/Magistrale/1° anno/2° - Machine learning for bioengineering/R_codes")

ct <- read.table("lez_04_CellType.dat")
ct$Bcell <- as.factor(ct$Bcell)
ct$Acurrent <- as.factor(ct$Acurrent)
ct$TailCurrent <- as.factor(ct$TailCurrent)
ct$Acell <- as.factor(ifelse(ct$CellType=="A","yes","no"))
ct$Dcell <- as.factor(ifelse(ct$CellType=="D","yes","no"))
ct$CellType <- as.factor(ct$CellType)
head(ct)
detach(ct)
attach(ct)
library(caret)
library(nnet)
library(klaR)
library(e1071)
library(tree)
library(randomForest)

table(Bcell)
table(Acell)
table(Dcell)

train.control <- trainControl(method="LOOCV")

set.seed(10)
ct_idx <- sample(nrow(ct),100) # samples 100 data points
ct_trn <- ct[ct_idx,]
ct_tst <- ct[-ct_idx,]


################################# B CELLS ######################################

plot(Ccell~V2h,col=Bcell,pch=ifelse(Bcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(Ccell~kh,col=Bcell,pch=ifelse(Bcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(kh~V2h,col=Bcell,pch=ifelse(Bcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(V2h~Bcell) # dal grafo si nota che è la var migliore per separare le Bcell
plot(Ccell~Bcell)
plot(kh~Bcell)

plot(V2h~CellType)
plot(Ccell~CellType)
plot(kh~CellType)
plot(Imax~CellType)


###

logi0B <- glm( Bcell ~ kh ,
               binomial( link = "logit" ) ,
               data = ct )

summary(logi0B)

model0B <- train( formula( logi0B ) ,
                  data = ct ,
                  method = "glm" ,
                  family = binomial() ,
                  trControl = train.control )

print(model0B)


###
logi1B <- glm(Bcell~V2h,binomial(link="logit"),data=ct)
summary(logi1B)
model1B <- train(formula(logi1B),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model1B)

###
logi2B <- glm(Bcell~Ccell,binomial(link="logit"),data=ct)
summary(logi2B)
model2B <- train(formula(logi2B), data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model2B)

###
logi3B <- glm(Bcell~V2h+kh+Ccell,binomial(link="logit"),data=ct)
summary(logi3B)
model3B <- train(formula(logi3B),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model3B)

newdata3 <- data.frame( V2h = -60 , Ccell = 4 , kh = 10 )
ppB = predict( logi3B , newdata3  ) # se aggiungo type = "response" calcola direttamente
1/(1+exp(-ppB))

newdata4 <- data.frame(V2h=-80,Ccell=6,kh=13)
ppB = predict(logi3B,newdata4)
1/(1+exp(-ppB))

###

logiB <- glm( Bcell ~ Imax*Ccell ,
              binomial( link = "logit" ) ,
              data = ct )

summary( logiB )

modelB <- train( formula( logiB ) ,
                 data = ct ,
                 method = "glm" ,
                 family = binomial() ,
                 trControl = train.control )

print( modelB )
# probabilità che un nuovo dato sia una cellula B
newdata3 <- data.frame(Imax=-900,Ccell=5)
ppB = predict(logiB,newdata3)
1/(1+exp(-ppB))




################################# A CELLS ######################################

plot(Ccell~V2h,col=Acell,pch=ifelse(Acell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(Ccell~kh,col=Acell,pch=ifelse(Acell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(kh~V2h,col=Acell,pch=ifelse(Acell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(V2h~Acell) # V2h è la miglior var per separare anche Acell
plot(Ccell~Acell)
plot(kh~Acell)

###
logi0A <- glm(Acell~kh,binomial(link="logit"),data=ct)
summary(logi0A)
model0A <- train(formula(logi0A), data=ct, method="glm",family=binomial(),
                 trControl=train.control)
print(model0A)

###
logi1A <- glm(Acell~V2h, binomial(link="logit"),data=ct)
summary(logi1A)
model1A <- train(formula(logi1A),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model1A)

###
logi2A <- glm(Acell~Ccell,binomial(link="logit"),data=ct)
summary(logi2A)
model2A <- train(formula(logi2A),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model2A)

###
logi3A <- glm(Acell~V2h+kh+Ccell,binomial(link="logit"),data=ct)
summary(logi3A)
model3A <- train(formula(logi3A),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model3A)

###
logiA <- glm(Acell~Imax*Ccell,binomial(link="logit"),data=ct)
summary(logiA)
modelA <- train(formula(logiA),data=ct,method="glm",family=binomial(),
                trControl=train.control)
print(modelA)
# probabilità che un nuovo dato sia una cellula A
ppA = predict(logiA,newdata3)
1/(1+exp(-ppA))

logiA.1 <- glm(Acell~Imax+Ccell, binomial(link="logit"), data=ct)

# probabilità della griglia di valori dei Ccell e Imax
CC <- seq(2.5,15,by=0.5)
II <- seq(-1700,-100,by=50)
gr <- expand.grid(Imax=II,Ccell=CC) # crea tuttle le possibili combinazioni di II e CC
head(gr)

gr.predict.1 <- 1/(1+exp(-predict(logiA.1,gr)))
gr.classA.1 <- as.logical(gr.predict>0.5)
gr.classA.1 <- matrix(gr.classA.1,nrow=length(II)) # crea una matrice di nrow righe
# riempendole con i valori di gr.classA
image(II,CC,gr.classA.1,xlab="Imax",ylab="Ccell")
points(Ccell~Imax,col=Acell,pch=19)

gr.predict <- 1/(1+exp(-predict(logiA,gr)))
gr.classA <- as.logical(gr.predict>0.5)
gr.classA <- matrix(gr.classA,nrow=length(II))
image(II,CC,gr.classA,xlab="Imax",ylab="Ccell")
points(Ccell~Imax,col=Acell,pch=19)




################################# D CELLS ######################################

plot(Ccell~V2h,col=Dcell,pch=ifelse(Dcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(Ccell~kh,col=Dcell,pch=ifelse(Dcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(kh~V2h,col=Dcell,pch=ifelse(Dcell=="yes",19,2))
legend("topright",c("B","A or D"),col=c(2,1),pch=c(19,2))

plot(V2h~Dcell)
plot(Ccell~Dcell)
plot(kh~Dcell)

# per Dcells anche kh permette di dividere abbastaza bene, ma forse resta
# migliore la variabile V2h


###
logi0D <- glm(Dcell~kh,binomial(link="logit"),data=ct)
summary(logi0D)
model0D <- train(formula(logi0D),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model0D)

###
logi1D <- glm(Dcell~V2h, binomial(link="logit"), data=ct)
summary(logi1D)
model1D <- train(formula(logi1D),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model1D)

###
logi2D <- glm(Dcell~Ccell,binomial(link="logit"),data=ct)
summary(logi2D)
model2D <- train(formula(logi2D),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model2D)

###
logi3D <- glm(Dcell~V2h+kh+Ccell,binomial(link="logit"),data=ct)
summary(logi3D)
model3D <- train(formula(logi3B),data=ct,method="glm",family=binomial(),
                 trControl=train.control)
print(model3D)

###
logiD <- glm(Dcell~Imax*Ccell, binomial(link="logit"), data=ct)
summary(logiD)
modelD <- train(formula(logiD),data=ct,method="glm",family=binomial(),
                trControl=train.control)
print(modelD)
# probabilità che un nuovo dato sia una cellula D
ppD = predict(logiD,newdata3)
1/(1+exp(-ppD))




############################ MULTINOMIAL REGRESSION ############################


# richiede la fattorizzazione della var CellType

plot(Acurrent~CellType)
plot(TailCurrent~CellType)
plot(Leak~CellType)
plot(Ccell~CellType)
plot(Raccess~CellType)
plot(V2h~CellType)
plot(kh~CellType)
plot(Imax~CellType)


###

ct$CellType<-relevel(ct$CellType,"A")
# stabilisco il celltype B come classe di riferimento
# la odds che ottengo è riferita a tale classe

controlvar <- trainControl(method="LOOCV")

fit1 <- multinom(CellType~V2h,data=ct)
summary(fit1)
control1 <- train(formula(fit1), data=ct,method="multinom",trControl=controlvar,
                  trace=FALSE)
print(control1)

fit2 <- multinom(CellType~V2h+kh)
summary(fit2)
control2 <- train(formula(fit2), data=ct,method="multinom",trControl=controlvar,
                  trace=FALSE)
print(control2)

fit3 <- multinom(CellType~V2h*kh)
summary(fit3)
control3 <- train(formula(fit3), data=ct,method="multinom",trControl=controlvar,
                  trace=FALSE)
print(control3)

newdata <- data.frame(kh=6,V2h=-90,kh=12)
ppn <- predict(fit3,newdata)
ppn


########################### NAIVE BAYES CLASSIFIER #############################

# !!!!!!!!! CANNOT INCLUDE INTERACTIONS !!!!!!!!!

nb0B <- NaiveBayes( formula(logi0B) , data=ct )
plot(nb0B)
nb0B.test <- train( formula(logi0B) ,
                    data=ct,
                    method="nb",
                    trControl=controlvar)
print(nb0B.test)

nb1B <- NaiveBayes(formula(logi1B),data=ct)
plot(nb1B)
nb1B.test <- train(formula(logi1B),data=ct,method="nb",trControl=controlvar)
print(nb1B.test)

nb2B <- NaiveBayes(formula(logi2B),data=ct)
plot(nb2B)
nb2B.test <- train(formula(logi2B),data=ct,method="nb",trControl=controlvar)
print(nb2B.test)

nb3B <- NaiveBayes(formula(logi3B),data=ct)
plot(nb3B)
nb3B.test <- train(formula(logi3B),data=ct,method="nb",trControl=controlvar)
print(nb3B.test)
# errori a cazzo

nb <- NaiveBayes(Bcell~kh+V2h,data=ct)
nb.train <- train(Bcell~kh+V2h,data=ct,methd="nb",trControl=controlvar)
print(nb.train)
# errori a cazzo anche qua

nbb <- NaiveBayes(Acell~Imax+Ccell,data=ct)
names(nbb)
nbb$apriori # restituissce la probabilità no yes di essere Acell a priori (ratio)
nbb$tables # restituisce avg e std dei fit con var gaus delle variab  considerate
# (Imax e Ccell) sei per la Acell a priori che per le notAcell a priori
plot(nbb) # plotta num della formula di Bayes per ciascuna variabile

newdata <- data.frame(Imax=-400,Ccell=5)
pp = predict(nbb,newdata)

gr.predict.nbb <- predict(nbb,gr)
gr.classA.nbb <- as.logical(gr.predict.nbb$class=="yes")
gr.classA.nbb <- matrix(gr.classA.nbb,nrow=length(II))
image(II,CC,gr.classA.nbb,xlab="Imax",ylab="Ccell")
points(Ccell~Imax,col=Acell,pch=19)


# posso predire anche più classi

nbmulti1 <-NaiveBayes(CellType~V2h,data=ct)
#plot(nbmulti1)
nbmulti1.train <- train(CellType~V2h,data=ct,method="nb",trControl=controlvar)
nbmulti1.train

newdata <- data.frame(V2h=-90)
pp = predict(nbmulti1,newdata)

nbmulti2 <-NaiveBayes(CellType~V2h+kh,data=ct)
#plot(nbmulti2)
nbmulti2.train <- train(CellType~V2h+kh,data=ct,method="nb",trControl=controlvar)
nbmulti2.train



############################ SUPPORT VECTOR MACHINE ############################


svm.Acell <- svm(Acell~Ccell+Imax,data=ct,kernel="linear")
plot(svm.Acell,data=ct,Ccell~Imax,grid=100)

# grid tells how fine you want the grid
# red ones Acells, crosses are data that contribute to loss function: if they are
# removed the separation line changes
# prima predizione corretta che non ha un peso determina il limite del margine
# all'interno del quale si hanno pesi piccoli nella hinge loss

# HOW TO READ THE PLOT
# black --> not Acell
# red --> Acell
# X --> data for which correspond a loss inside the Hinge Loss
# so X represent MISCLASSIFICATION or OBSERVATION inside the MARGIN
# O --> correct prediction that have 0 loss

# confrontare le predizioni dei tre modelli per le Acells: logis, bayes, svm

cv.logiA <- train(formula(logiA),data=ct,method="glm",family="binomial",
                  trControl=train.control)
cv.logiA

cv.logiA.1 <- train(formula(logiA.1),data=ct,method="glm",family="binomial",
                    trControl=train.control)
cv.logiA.1

cv.nbb <- train(Acell~Imax+Ccell,data=ct,method="nb",trControl=train.control)
cv.nbb

cv.svm.Acell <- train(formula(svm.Acell),data=ct,method="svmLinear",
                      trControl=train.control,
                      preProcess = c("center","scale"),tuneLength=10)
cv.svm.Acell


# tune_out <- tune(svm, Acell~Ccell+Imax, data=ct, kernel="polynomial",
# ranges = list(cost = c(0,1,10,100,1000), degree = c(1,2,3,4)))
# Error in svm.default(x, y, scale = scale, ..., na.action = na.action) : 
# C <= 0!
  


##################################### ROC #####################################

set.seed(43)
idx <- sample(nrow(ct),150)
trn <- ct[idx,]
tst <- ct[-idx,]

fit.multi <- multinom(CellType~V2h,data=trn)
pred <- predict(fit.multi,tst)
table(pred,tst$CellType)

confusionMatrix(pred,tst$CellType)

# class specific sens and spec can be calculated from the table:
# sensitivity class A: 10/(10+1+11)
# specificity class A: (10+12)/(6+10+12)
# so on so forth for the other classes
# rivedere la definizione teorica


##################################### SVM #####################################

# cost is lambda, degree is for polynomial kernels

svm.lin <- svm(Acell~Ccell+Imax,data=ct,kernel="linear")
plot(svm.lin,data=ct,Ccell~Imax,grid=100)
tst.svm.lin <- train(formula(svm.lin),data=ct,method="svmLinear",
                     trControl=train.control,
                     preprocess=c('centre','scale'),tuneLength=10)

svm.poli.1 <- svm(Acell~Ccell+Imax,data=ct,kernel="polynomial",cost=1,degree=3)
plot(svm.poli.1,data=ct,Ccell~Imax,grid=100)
tst.svm.poli.1 <- train(formula(svm.poli.1),data=ct,method="svmPoly",
                     trControl=train.control,
                     preprocess=c('centre','scale'),tuneLength=10)

svm.poli.2 <- svm(Acell~Ccell+Imax,data=ct,kernel="polynomial",cost=10,degree=3)
plot(svm.poli.2,data=ct,Ccell~Imax,grid=100)

svm.poli.3 <- svm(Acell~Ccell+Imax,data=ct,kernel="polynomial",cost=1,degree=10)
plot(svm.poli.3,data=ct,Ccell~Imax,grid=100)

# aumentando il costo, una predizione sbagliata è "più grave"
# e stiamo forzando i confini ad includere più punti possibile al costo di aver
# una norma di w più grande con il proble di non predire benissimi nuovi dati
# se lo diminuisco, una predizione sbagliata o l'essere vicino al margine viene
# consoderato meno grave e penalizzato poco
# cambiando il grado, aumentandolo spingo per predire tutte le cellule, anche le
# più lontante facendo overfitting del training set e predicendo male nuovi dati


# per valutare la bontà di predizione dei valori di costo e grado uso tune
tune_out.poly <- tune( svm , Acell ~ Ccell + Imax , data = ct , kernel = "polynomial" ,
                       ranges = list( cost = c(0.1,1,10,100,100) , degree = c(0.5,1,2,3,4)))
# ranges: valori da provare per i parametri specificati, prova le combinazioni
# e restituisce quella con la crossvalidation minore

svm.rad.1 <- svm(Acell~Ccell+Imax,data=ct,kernel="radial",cost=0.5,gamma=10)
plot(svm.rad.1,data=ct,Ccell~Imax,grid=100)

svm.rad.2 <- svm(Acell~Ccell+Imax,data=ct,kernel="radial",cost=1,gamma=10)
plot(svm.rad.2,data=ct,Ccell~Imax,grid=100)

svm.rad.3 <- svm(Acell~Ccell+Imax,data=ct,kernel="radial",cost=1,gamma=1)
plot(svm.rad.3,data=ct,Ccell~Imax,grid=100)

svm.rad.4 <- svm(Acell~Ccell+Imax,data=ct,kernel="radial",cost=0.1,gamma=1)
plot(svm.rad.4,data=ct,Ccell~Imax,grid=100)

# gamma = 1/(2*sigma) sigma measures what I mean by "close"

tune_out.radial <- tune(svm,Acell~Ccell+Imax,data=ct,kernel="radial",
                    ranges=list(cost=c(0.1,1,10,100,100),gamma=c(0.5,1,2,3,4)))

train.control <- trainControl(method="LOOCV")

svm.lin.cv <- train(Acell~Ccell+Imax,data=ct,method="svmLinear",
                    trControl=train.control,preProcess=c("center","scale"))
svm.lin.cv

svm.rad.cv <- train(Acell~Ccell+Imax,data=ct,method="svmRadial",
                    trControl=train.control,preProcess=c("center","scale"),
                    tuneLength=10)
svm.rad.cv

# testa sui valori del costo calcolando il valore di sigma con metodo derivativo



################################## MULTI SVM ##################################


multisvm.0 <- svm(CellType~Ccell+Imax,data=ct,kernel="linear")
plot(multisvm.0,data=ct,Ccell~Imax,grid=100)
# uses the voting system, aka one-vs-one
# i want to try to build this on my own

ctAB <- ct[CellType!="D",]
multisvm.AB <- svm(CellType~Ccell+Imax,data=ctAB,kernel="linear")
dev.new()
plot(multisvm.AB,data=ctAB,Ccell~Imax,grid=100)


ctAD <- ct[CellType!="B",]
multisvm.AD <- svm(CellType~Ccell+Imax,data=ctAD,kernel="linear")
dev.new()
plot(multisvm.AD,data=ctAD,Ccell~Imax,grid=100)


ctBD <- ct[CellType!="A",]
multisvm.BD <- svm(CellType~Ccell+Imax,data=ctBD,kernel="linear")
dev.new()
plot(multisvm.BD,data=ctBD,Ccell~Imax,grid=100)



multisvm <- svm(CellType~.-Bcell-Acell-Dcell,data=ct,kernel="linear")
summary(multisvm)

set.seed(22)
#tune_out.multi <- tune(svm,CellType~.-Bcell-Acell-Dcell,data=ct,kernel="linear",
#                       ranges=list(cost=c(0.1,1,3,7,10,100)))
#summary(tune_out.multi)
#plot(tune_out.multi$best.model,data=ct,V2h~Ccell,grid=100) 
# per qualche ragione sconosciuta tale comando non funge
# per suddetto motivo siamo costretti a fare trick illegali del tipo:

tune_out.multi <- tune(svm,CellType~.,data=ct[,1:9],kernel="linear",
                       ranges=list(cost=c(0.1,1,2,3,4,5,7,10,100)),
                       tunecontrol=tune.control(cross=nrow(ct)))
# cross dice quanti folds voglio nella cross validazione
summary(tune_out.multi)
plot(tune_out.multi$best.model,data=ct[,1:9],V2h~Ccell,grid=100)


# can do the same using train
multisvm.cv <- train(CellType~.-Acell-Bcell-Dcell,data=ct,method="svmPoly",
                     trControl=train.control)
multisvm.cv



# usando caret e tuneGrid

svmtune.1 <- train(CellType~.-Acell-Bcell-Dcell,data=ct,method="svmPoly",
                   trControl=trainControl(method="cv",number=10),tunelength=4)
svmtune.1

svmtune.2 <- train(CellType~.-Acell-Bcell-Dcell,data=ct,method="svmPoly",
                   trControl=trainControl(method="cv",number=10),
                   tuneGrid=expand.grid(C=c(0.5,0.8,1.2),degree=c(1,2,3),scale=1))
svmtune.2


################################ DECISION TREES ################################


cttree <- tree(CellType~.-Acell-Bcell-Dcell,data=ct)
summary(cttree)
# residual mean deviance is similar to residuals sum of squares
plot(cttree)
text(cttree,pretty=0)
cttree
# if condition is true the tree goes to the left branch
# V2h good at separating beta cells from the rest

# CV simulation with for loop
tot1 <- 0
for (i in seq(from=10,to=nrow(ct),by=10)) {
  txt_idx <- (i-9):i
  train <- ct[-txt_idx,]
  test <- ct[txt_idx,]
  cv.tree <- tree(CellType~.-Acell-Bcell-Dcell,data=train)
  prediction <- predict(cv.tree,test,type="class")
  misclassi <- prediction != test$CellType
  tot1 <- tot1+sum(misclassi==T)
}

tot2 <- 0
for (i in 1:nrow(ct)) {
  train <- ct[-i,]
  test <- ct[i,]
  cv.tree <- tree(CellType~.-Acell-Bcell-Dcell,data=train)
  prediction <- predict(cv.tree,test,type="class")
  if(prediction!=test$CellType) {tot2 <- tot2+1}
}


# this is a way to find the best pruning tree based on crossvalidation
# find the best pruning tree based on 9 of the folds
cv.cttree <- cv.tree(cttree,FUN=prune.misclass)
plot(cv.cttree,type="b")
cv.cttree
# misclass: if I cut away a branch, will the misclassification (error I make)
# suffer? If not the pruning is good
cv.cttree2 <- cv.tree(cttree,FUN=prune.misclass,K=nrow(ct))
# in questo modo faccio una loocv con k num di folds pari al num di dati
cv.cttree2
plot(cv.cttree2,type="b")
# nel grafico cambia perchè si basa su cv
# size tells the num of nodes in the three
# dev tells how good it is at predicting
# looking at the misclassification graph there is a number of nodes at which
# the value drops

# this way i prune without doing cross validation
# get sth similar
misclass <- prune.misclass(cttree)
plot(misclass)
misclass
lines(misclass$size,misclass$dev,col="red",type="b")
# nel grafico resta sempre lo stesso
# this one instead tries to  prune the tree based on the complete training data
# not via cross-validation
# we can see the fit is better on the complete training data (final part)
# for the misclass object

# PRUNING: cutting away some branches ending elements of the tree -> leaf
prune.cttree <- prune.misclass(cttree,best=4)
summary(prune.cttree)
# the misclassification seen here (last row) is the one on the training data!
# not very informative for future predictions
dev.new()
plot(prune.cttree)
text(prune.cttree,pretty=0)

# to solve this since there is no cross validation in this package,
# we split the dataset into train and test sets
set.seed(431)


# train on the training data
cttree.trn <- tree(CellType~.-Acell-Bcell-Dcell,data=ct_trn)
plot(cttree.trn)
text(cttree.trn,pretty=0)

# try to predict the test data
cttree.pred <- predict(cttree.trn,ct_tst,type="class")
confusionMatrix(cttree.pred,ct_tst$CellType)

# cerca il punto migliore di pruning valutando la misclassification sul test data
misclass.trn <- prune.misclass(cttree.trn,newdata=ct_tst)
plot(misclass.trn)


cv.cttree <- cv.tree(cttree.trn,FUN=prune.misclass)
# cv.cttree <- cv.tree(cttree.trn,FUN=prune.misclass,K=nrow(ct))
# implemetns loocv
cv.cttree
plot(cv.cttree)

prune.cttree <- prune.misclass(cttree.trn,best=4)
plot(prune.cttree)
text(prune.cttree,pretty=0)
prune.cttree.pred <- predict(prune.cttree,ct_tst,type="class")
confusionMatrix(prune.cttree.pred,ct_tst$CellType)
# any improvement but tree easier to interpret

# con caret package

caret.cttree <- train(CellType~.-Acell-Bcell-Dcell,data=ct,method="rpart",
                      trControl=train.control,tuneLength=10)
caret.cttree
par(mfrow=c(1,2),xpd=NA)
plot(caret.cttree$finalModel)
text(caret.cttree$finalModel,all=T,cex=.8)
library(rattle)
fancyRpartPlot(caret.cttree$finalModel)


# impurity
ctrl <- tree.control(nobs=nrow(ct),mincut=10,minsize=20,mindev=0.01)
# mincut: min num of obs to include in either child node
# minsize: smallest allowed node size. default=10
cttree.gini <- tree(CellType~.-Acell-Bcell-Dcell,data=ct,split="gini",control=ctrl)
plot(cttree.gini)
text(cttree.gini,pretty=0)
summary(cttree.gini)

cttree.dev <- tree(CellType~.-Acell-Bcell-Dcell,data=ct,split="deviance",control=ctrl)
plot(cttree.dev)
text(cttree.dev,pretty=0)
summary(cttree.dev)



################################ RANDOM FOREST ################################

# choosing too many random variable to do each split we can fall into the usual
# problem of overfitting

set.seed(111)


ctrf <- randomForest(CellType~.-Acell-Bcell-Dcell,data=ct_trn,importance=T)
ctrf

ctrf.predict <- predict(ctrf,ct_tst)

table( pred = ctrf.predict , truth = ct_tst$CellType )


### tuning on all dataset

rf.1 <- train(CellType~.-Acell-Bcell-Dcell,data=ct,tuneLength=7,method="rf",
              trControl=train.control)
# this up builds a rf using train with caret pckg
rf.1
# the parameter changing in the cv is mtry, num of tree doesn't change much

rf.2 <- train(CellType~.-Acell-Bcell-Dcell,data=ct,tuneGrid=expand.grid(.mtry=(1:8)),
              method="rf",trControl=train.control)
rf.2


# to evaluate the best value for mtry
best.mtry <- tuneRF(ct[,2:8],ct[,1],stepFactor=1.5,improve=1e-5,ntree=500)
# it needs the column of predictors and responses, not the rf formula :(
# stepFactor -> how much I increase mtry at each attempt
# improve ->



############################## NEAREST NEIGHBOUR ##############################

library(FNN)


set.seed(101)
ct_idx <- sample(nrow(ct),100) # samples 100 data points
ct_trn <- ct[ct_idx,]
ct_tst <- ct[-ct_idx,]


prednn.1 <- knn( ct_trn[ , 4:9] ,
                 ct_tst[ , 4:9] ,
                 ct_trn$CellType )
# fun knn needs numerical covariates to work with so I have not to take the
# others into consideration
# retituisce direttamente la classe scelta in base al vicino, l'indice del vicino
# e la distanza a cui si trova

confusionMatrix( prednn.1 , ct_tst$CellType ) # 0.6 acc

# poblems: I had to ignore tail current and A current which were two of the best
# predictors and I din't scale the data

ctnum <- cbind( as.numeric( ct$Acurrent == "yes" ) ,
                as.numeric( ct$TailCurrent == "yes" ) ,
                ct[ , 4:9 ] )

names( ctnum ) <- names( ct[ 2:9 ] )

ctsc <- scale(ctnum)


prednn.2 <- knn( ctsc[ ct_idx , ] ,
                 ctsc[ -ct_idx , ] ,
                 ct_trn$CellType )

confusionMatrix( prednn.2 , ct_tst$CellType ) # 0.82 acc





knn.1 <- train( CellType ~ . - Acell - Bcell- Dcell ,
                data = ct_trn ,
                tuneLength = 20 ,
                method = "knn" ,
                preProcess = c( "center" , "scale" ) ,
                trControl = train.control )
knn.1
plot( knn.1 )


knn.2 <- train( CellType ~ . - Acell - Bcell - Dcell ,
                data = ct_trn ,
                trControl = train.control ,
                method = "knn" ,
                preProcess = c( "center" , "scale" ) ,
                tuneGrid = expand.grid( k = 1:20 ) )
knn.2
plot( knn.2 )


prednn.2 <- knn( ctsc[ ct_idx , ] ,
                 ctsc[ -ct_idx , ] ,
                 ct_trn$CellType, k = 4 )

confusionMatrix(prednn.2,ct_tst$CellType)



# if i don't scale the covariates with bigger values would dominate even if they are not good predicors


################################ CLUSTERING ####################################



ct <- read.table("lez_04_CellType.dat")
ct$Acell <- as.factor(ifelse(ct$CellType=="A","yes","no"))
detach(ct)
attach(ct)


library(NbClust)
library(factoextra)
library(cluster)

# package factoextra has fviz for clustering visualization
# fviz_nbclust plotta il valore di un parametro (tipo wrss) al variare del
# numero di clusters (#max come input)



head(ct)

ctnum <- cbind( as.numeric( ct$Acurrent == "yes" ),
                as.numeric( ct$TailCurrent == "yes" ),
                ct[,4:8] )

summary(ctnum)

ctsc <- as.data.frame( scale( ctnum ) )
# all mean = 0

summary( ctsc )


normalize <- function(x){ ( x-min(x) ) / ( max(x)-min(x) ) }

ctnorm <- as.data.frame( sapply( ctnum,normalize ) )
# all min = 0, all max = 1
# sapply to avoid for loops to apply the normalize function to the numeric 
# dataframe column by column. la fun non lavora sulle singole colonne ma prendebbe
# min e max su tutto il data frame, mentre noi vogliamo scalare le colnne
# ?apply

names(ctsc)[1:2] <- names(ctnorm)[1:2] <- c("Acurrent","TailCurrent")

summary(ctnorm)


CellTypeNum <-  ifelse( ct$CellType == "A" , 1 , ifelse( ct$CellType == "B", 2 , 4 ) )
# per assegnare un valore numerico al tipo di cellula per poterlo vedere come cluster

plot( ctsc, col = CellTypeNum , pch = CellTypeNum )
plot( ctnorm, col = CellTypeNum , pch = CellTypeNum )


### K means

# 1 cluster
# A vector of integers (1:k) indicating the cluster to which each point is allocated.
# 2 centers
# A matrix of cluster centres.
# 3 totss
# The total sum of squares. ( no clusters )
# 4 withinss
# Vector of within-cluster sum of squares, one component per cluster.
# 5 tot.withinss
# Total within-cluster sum of squares, i.e. sum(withinss).
# 6 betweenss
# The between-cluster sum of squares, i.e. =totss-tot.withinss.
# 7 size
# The number of points in each cluster.
# 8 iter
# The number of (outer) iterations.
# 9 ifault
# integer: indicator of a possible algorithm problem - for experts.

k_sc <- kmeans( ctsc , centers = 3 , nstart = 1 )

k_norm <- kmeans( ctnorm , centers = 3 , nstart = 1 )

table(k_sc$cluster,k_norm$cluster)

x <- table(k_sc$cluster,ct$CellType)

x

# accuracy
sum( apply( x , 1 , max ) ) / sum( x )
x <- table( k_norm$cluster , ct$CellType )
x
sum( apply( x , 1 , max ) ) / sum( x )

k_sc$betweenss / k_sc$totss
k_norm$betweenss / k_norm$totss


###

plot( ct$Ccell ~ ct$V2h, col = CellTypeNum , pch = k_norm$cluster )

fviz_cluster( k_norm , ctnorm[,c(6,4)] , ellipse.type = "norm" )

sil <- silhouette( k_norm$cluster , dist( ctnorm ) )
# dist calculates all the differnent distances between all the data points
fviz_silhouette( sil )


# how many clusters?

# within sum of squares
fviz_nbclust( ctnorm , kmeans , method = "wss", k.max = 30 )
fviz_nbclust( ctsc , kmeans , method = "wss" , k.max = 30 )

fviz_nbclust( ctnorm , kmeans , method = "silhouette" , k.max = 30 )
fviz_nbclust( ctsc , kmeans , method = "silhouette" , k.max = 30 )

fviz_nbclust( ctnorm , kmeans , method = "gap_stat" , k.max = 30 )
fviz_nbclust( ctsc , kmeans , method = "gap_stat" , k.max = 30 )



k_norm4 <- kmeans( ctnorm , centers = 4 , nstart = 10 )
table( k_norm4$cluster , ct$CellType )
k_norm4$betweenss / k_norm4$totss # 0.75

k_norm5<-kmeans( ctnorm , centers = 5 , nstart = 5 )
table( k_norm5$cluster , ct$CellType )
k_norm5$betweenss / k_norm5$totss # 0.79

fviz_cluster( k_norm5 , ctnorm[,c(6,4)] , ellipse.type = "norm" )
sil <- silhouette( k_norm5$cluster , dist( ctnorm ) )
fviz_silhouette( sil )

##
k_Ccell_V2h <- kmeans( ctnorm[,c(4,6)] , centers = 3 )
table( k_Ccell_V2h$cluster , ct$CellType )

plot( ctnorm$V2h ~ ctnorm$Ccell , col = CellTypeNum , pch = k_Ccell_V2h$cluster )
points( k_Ccell_V2h$center , cex = 3, pch = 20, col = 3 )

fviz_cluster( k_Ccell_V2h , ctnorm[,c(4,6)] , ellipse.type = "norm" )


k_Tail_V2h <- kmeans( ctnorm[,c(2,6)] , centers = 3 )
table( k_Tail_V2h$cluster , ct$CellType )
plot( ctnorm$TailCurrent , ctnorm$V2h , col = CellTypeNum , pch = k_Tail_V2h$cluster )
points( k_Tail_V2h$center , cex = 3 , pch = 20 , col = 3 )



#### Agglomerate (hierarchical clustering)

str( ctnorm )
dist_mat <- dist( ctnorm )

hclust_avg <- hclust( dist_mat , method = 'average' )
plot( hclust_avg )
# quando ci sono distanze verticai molto elevate, significa che l'accoppiamento
# dei cluster sta creando molto cambiamento ( sto unendo dati che in realtà
# sarebbero molto distanti )

rect.hclust( hclust_avg , k = 4 , border = 2:6 )
# disegna k  rettangoli che corrispondono ai cluster
cut_avg <- cutree( hclust_avg , k = 4 )
# assegna i diversi elementi ai rispettivi cluster in base a dove taglio (k)
cut_avg
table( cut_avg , ct$CellType )

hclust_sgl <- hclust( dist_mat , method = 'single' )
plot( hclust_sgl )
# cambiando la distanza a min tra due cluster si creano dei cluster che hanno
# forma allungata
rect.hclust( hclust_sgl , k = 5 , border = 2:4 )
cut_avg <- cutree( hclust_sgl , k=5 )
table( cut_avg , ct$CellType )

hclust_cplt <- hclust( dist_mat , method = 'complete' )
plot( hclust_cplt )
rect.hclust( hclust_cplt , k=5, border=2:4 )
cut_avg <- cutree( hclust_cplt , k=5 )
table( cut_avg , ct$CellType )

### k-Medoids

# in k means we risk our center does not correspond to a data point because we
# are performing averaging (example we get a point with sex btw female and male)

# partition around midoids
kmed <- pam( ctnorm , k = 3 )
# here the centers correspond to data items which minimizes the distance from
# other data points of the cluster as the centers in k means

table( kmed$cluster , ct$CellType )
# fviz_cluster(kmed)





