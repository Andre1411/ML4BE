
### Multiple imputation with mice
rm(list = ls())
library(mice)
df<- nhanes
summary(df)
df
md.pattern(df)

fit <- lm(age~bmi, data=df)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  3.76718    1.31945   2.855   0.0127 *
# bmi         -0.07359    0.04910  -1.499   0.1561 

imp <- mice(df, method="mean", m=1, maxit=1)
complete(imp)

fit <- lm(age~bmi, data=complete(imp))
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  3.71468    1.32901   2.795   0.0103 *
# bmi         -0.07359    0.04966  -1.482   0.1520  
fit <- with(imp, lm(age~bmi))
summary(fit)
# term          estimate std.error statistic p.value  nobs
# 1 (Intercept)   3.71      1.33        2.80  0.0103    25
# 2 bmi          -0.0736    0.0497     -1.48  0.152     25

imp <- mice(df, method = "norm.predict", m = 1, maxit = 1)  # use predicted values
fit <- with(imp, lm(age~bmi))
summary(fit) # BMI has significant effect!
# term          estimate std.error statistic   p.value  nobs
# 1 (Intercept)    4.84     0.899       5.39 0.0000180    25
# 2 bmi           -0.115    0.0331     -3.47 0.00208      25

plot(chl~bmi,data=df)
points(chl~bmi, data = complete(imp), pch=4)


imp <- mice(df, method = "norm.nob", m = 1, maxit = 1)		# use predicted values + noise
fit <- with(imp, lm(age~bmi))
summary(fit)
# term          estimate std.error statistic  p.value  nobs
# 1 (Intercept)    4.53     1.02        4.46 0.000179    25
# 2 bmi           -0.104    0.0377     -2.76 0.0112      25

imp <- mice(df, method = "norm", m = 1, maxit = 1)		# use predicted values + noise + uncertainty
fit <- with(imp, lm(age~bmi))
summary(fit)

imp <- mice(df, method = "norm.boot", m = 1, maxit = 1) # use predicted values + noise + uncertainty
fit <- with(imp, lm(age~bmi))
summary(fit)

set.seed(123)
imp <- mice(df, method = "pmm", m=10, maxit=10, print=F)		# use "similar" observed values
imp$data # original data
imp$imp # imputated data
plot(imp)
c3 <- complete(imp, 3)
md.pattern(c3)
c.long <- complete(imp, "long")
c.long
c.broad <- complete(imp, "broad")
c.broad

fit <- with(imp, lm(age~bmi))
fit
summary(fit$analyses[[2]])

pool.fit <- pool(fit)
summary(pool.fit)
#          term    estimate  std.error statistic       df     p.value
# 1 (Intercept)  3.87035955 1.19968231  3.226154 15.95480 0.005294927
# 2         bmi -0.07928106 0.04459533 -1.777788 16.07811 0.094355143

## 
df <- nhanes2 # just age as a factor now and ht into factor
md.pattern(df) # same as before
library(rpart)
tree1 <- rpart(hyp~., df, cp=0.01, minsplit=1)
plot(tree1)
text(tree1)

imp <- mice(df, m=1, maxit=1)

tree2 <- rpart(hyp~., complete(imp), cp=0.01, minsplit=1)
plot(tree2)
text(tree2)


trn <- nhanes2[is.na(nhanes2$hyp)==F,]
tst <- nhanes2[is.na(nhanes2$hyp)==T,]
trn
tst
tree3 <- rpart(hyp~., trn, cp=0.01, minsplit=1)
plot(tree3)
text(tree3)
predict(tree3,tst, type="class")


df <- nhanes2
imp <- mice(df[,c(1,2,4)], m=1, maxit=1) # non voglio imputare ht
imp$pred
df[,c(1,2,4)] <- complete(imp)
df
train <- df[is.na(df$hyp)==F,]
test <- df[is.na(df$hyp)==T,]
train
test
tree4 <- rpart(hyp~., train, cp=0.01, minsplit=1)
plot(tree4)
text(tree4)
predict(tree4,test, type="class")
# facendo diverse run imputo valori diversi e le predizioni cambiano per questo


df <- nhanes2
imp <- mice(df[,c(1,2,4)], m=1, maxit=1, method="cart") # predice usando tree
df[,c(1,2,4)] <- complete(imp)
df
train <- df[is.na(df$hyp)==F,]
test <- df[is.na(df$hyp)==T,]
tree4 <- rpart(hyp~., train, cp=0.01, minsplit=1)
plot(tree4)
text(tree4)
predict(tree4,test, type="class")

##
fit <- with(train, rpart(hyp~., train, cp=0.01, minsplit=1))
fit

##
df <- nhanes2
m <- 10 		# no. imputations
imp <- mice(df[,c(1,2,4)], m=m, maxit=3, method="cart")
# imp <- mice(df[,c(1,2,4)], m=m, maxit=5, print=F)
df[,c(1,2,4)] <- complete(imp)
#df
train <- df[is.na(df$hyp)==F,]
test <- df[is.na(df$hyp)==T,]
tree4 <- rpart(hyp~., train, cp=0.01, minsplit=1)
plot(tree4)
text(tree4)
predict(tree4,test, type="class")

analyses <- as.list(rep(0,m))
for(i in 1:m){
  data <- as.data.frame(cbind(complete(imp,i), df$hyp))
  names(data) <- c("age","bmi","chl","hyp")
  train <- data[is.na(df$hyp)==F,]
  test <- data[is.na(df$hyp)==T,]
  tree <- rpart(hyp~., train, cp=0.01, minsplit=1)
  analyses[[i]] <- predict(tree4,test, type="class")
}

res <- apply(matrix(unlist(analyses), ncol=m), 1, function(x){table(factor(x,levels=c("no","yes")))})

res<- as.data.frame(t(res))
pred <- ifelse(res$no > res$yes, "no", "yes")
res

## using "with"
df <- nhanes2
imp1 <- mice(df, m=m, maxit=5, print=F, method=c("","pmm","","pmm")) # do not impute for age and hyp
p.matrix <- imp1$pred
p.matrix[,"hyp"] <- 0
imp1 <- mice(df, m=m, maxit=5, print=F, method=c("","pmm","","pmm"), predictorMatrix=p.matrix) # do not impute for age and hyp
complete(imp1,1)

expr <- expression(mean(bmi, na.rm=T))
fit1 <- with(imp1, eval(expr)) # applies function as written in "expr" on each imputed dataset

# fitting a tree and predict outcome where hyp is not known, for each imputed dataset
expr <- expression(
  data<- as.data.frame(cbind(age,bmi,hyp,chl)), 
  names(data) <- c("age","bmi","hyp","chl"),
  data$age <- as.factor(data$age),
  data$hyp <- as.factor(data$hyp),
  model <- rpart(hyp~. , data=data[is.na(data$hyp)==F,], cp=0.01, minsplit=1),
  predict(model,newdata=data[is.na(data$hyp)==T,], type="class")
)

fit <- with(imp1, eval(expr))
fit$analyses

# fit
analyses <- fit$analyses
res<- apply(matrix(unlist(analyses), ncol=m), 1, function(x){table(factor(x,levels=c("1","2")))})
res<- as.data.frame(t(res))
names(res) <- c("no", "yes")
pred <- ifelse(res$no > res$yes, "no", "yes")
res