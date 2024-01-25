
#########################################################################

setwd("C:/Users/ASUS/OneDrive - Università degli Studi di Padova/uni/Magistrale/1° anno/2° - Machine learning for bioengineering/R_codes")

# per capire le covariate significative per distinguere la classe si può fare
# modello logistico che prende in considerazione tutte le covariate e selezionare
# quelle che hanno pvalue significativo

data <- read.table("BreastCancer_complete.dat")
head(data)
data$class <- as.factor(data$class)

detach(data)
attach(data)
library(caret)
library(nnet)
library(klaR)
library(e1071)


################################################################################

svm.fit1 <- svm(class~adhesion+thickness, data=data, kernel="linear")
summary(svm.fit1)
plot(svm.fit1,data=data,adhesion~thickness,grid=100)


svm.fit2 <- svm(class~adhesion+thickness,data=data,
                kernel="radial",cost=0.5,gamma=10)
plot(svm.fit2,data=data,adhesion~thickness,grid=100)


svm.fit3 <- svm(class~adhesion+thickness,data=data,
                kernel="radial",cost=1,gamma=0.5)
plot(svm.fit3,data=data,adhesion~thickness,grid=100)


svm.fit3.2 <- svm(class~adhesion+thickness,data=data,
                  kernel="radial",cost=1,gamma=4.023792)
plot(svm.fit3.2,data=data,adhesion~thickness,grid=100)


svm.fit4 <- svm(class~adhesion+thickness,data=data,
                kernel="radial",cost=.1,gamma=1)
plot(svm.fit4,data=data,adhesion~thickness,grid=100)


svm.fit5 <- svm(class~adhesion+thickness,data=data,
                kernel="polynomial",cost=1,degree=3)
plot(svm.fit5,data=data,adhesion~thickness,grid=100)


svm.fit6 <- svm(class~adhesion+thickness,data=data,
                kernel="polynomial",cost=.1,degree=5)
plot(svm.fit6,data=data,adhesion~thickness,grid=100)




tune_out = tune(svm, class~adhesion+thickness, data = data, kernel = "radial",
                ranges = list(cost = c(0.1,1,10,100,1000), gamma = c(0.5,1,2,3,4)))
bestmod = tune_out$best.model
tune_out
# scegliamo il modello svm.fit3 (cost=1, gamma=0.5)
summary(bestmod)

################################ CV ################################

train.control <- trainControl(method="CV",number=10)

cv.svm.fit3 <- train(class~adhesion+thickness, data = data, method = "svmRadial",
                     trControl = train.control)



################################################################################


svm.fit7 <- svm(class~adhesion+thickness+unif.shape+chromatin, data = data,
                kernel="radial",cost=1, gamma=1.017803)

cv.svm.fit7 <- train(class~adhesion+thickness+unif.shape+chromatin, data = data,
                     method = "svmRadial", trControl = train.control)
cv.svm.fit7




