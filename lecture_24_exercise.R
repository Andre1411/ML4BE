
rm(list = ls())
data("PimaIndiansDiabetes2", package = "mlbench")
df <- PimaIndiansDiabetes2[,2:9]
detach()
attach(df)
df
summary(df)
md.pattern(df)

set.seed(123)
idx <- sample(nrow(df),round(nrow(df)*90/100))
trn <- df[idx,]
tst <- df[-idx,]

fit1 <- glm(diabetes~.,data=df,binomial(link=logit))
summary(fit1)


imp <- mice(df,m=10,maxit=10)
summary(complete(imp))

fit2 <- with(complete(imp),glm(diabetes~.,data=df,binomial(link=logit)))
summary(fit2)



