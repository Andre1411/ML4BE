
rm(list = ls())
setwd("C:/Users/ASUS/OneDrive - Università degli Studi di Padova/uni/Magistrale/1° anno/2° - Machine learning for bioengineering/R_codes")

data <- read.table("lez_02_bbd_data.dat")

data$diab <- as.factor(data$diab)
data$gender <- as.factor(data$gender)
library(caret)
detach(data)
attach(data)
summary(data)
plot(data)
nrow(data)

table(diab,gender) # relatively more males are diabetic

# par(mfcol=c(1,2))   # serve per fare subplots
plot(age ~ waist, col=diab, pch=ifelse(gender==1,2,20)) # circles female (gender=2)
plot(height~weight,col=diab, pch=ifelse(gender==1,2,20)) # circles female (gender=2)

plot(waist~age)
abline(a=0,b=1,col=2)

# fit <- lm(height~gender/weight-1)
fit <- lm(height~gender+weight)
summary(fit)

plot(BMI~diab)
plot(age~diab)

### outlier example
plot(weight~height,data=data[100:150,]) 
points(weight~height,data=data[c(147,120),],pch=20,col="red")

### fitting only data btw 100 and 150
fit <- lm(weight~height,data=data[100:150,])
summary(fit)
abline(fit)
anova(fit)
plot(fit)

### fitting data excluding red points (outliers) I improve R2
fit2 <- lm(weight~height,data=data[c(100:119,121:146,148:150),])
summary(fit2)
abline(fit2,col="blue")
anova(fit2)
plot(fit2)

### Blood pressure
plot(sysBP~weight,col=gender,pch=20,xlim=c(30,210), ylab="systolic blood pressure")

BPfit <- lm(sysBP~weight+weight:gender,data=data)
summary(BPfit)
abline(BPfit)
abline(a=coef(BPfit)[1], b=coef(BPfit)[2]+coef(BPfit)[3], col=2)
plot(BPfit)


###

bpfit1 <- lm( sysBP ~ gender + weight + height ,
              data = data )
summary(bpfit1)

controlvar <- trainControl( method = "LOOCV" )
model1 <- train( formula( bpfit1 ) ,
                 data = data[ complete.cases( data ) , ] ,
                 method = "lm",
                 trControl = controlvar )
print(model1)
summary(model1)


###

bpfit2 <- lm(sysBP~weight*gender+height)
summary(bpfit2)

model2 <- train(formula(bpfit2),
                data=data[complete.cases(data),],
                method="lm",
                trControl=controlvar)
print(model2)
summary(model2)


###

bpfit1.1 <- lm(sysBP~weight+gender+height-1)
summary(bpfit1.1)

model1.1 <- train(formula(bpfit1.1),
                  data=data[complete.cases(data),],
                  method="lm",
                  trControl=controlvar)
print(model1.1)
summary(model1.1)


### modelli logistici

logistic <- glm(diab~.,
                binomial(link="logit"),
                data=data[complete.cases(data),])
summary(logistic)
anova(logistic,test="Chisq")

logistic0 <- glm(diab~gender+age+waist+height,
                 binomial(link="logit"),
                 data=data[complete.cases(data),])
summary(logistic0)
anova(logistic0,logistic,test="Chisq")


