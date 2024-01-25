# Read the data from the file "lez_02_bbd_data.dat"
data <- read.table("lez_02_bbd_data.dat")

# Convert 'diab' and 'gender' columns to factors
data$diab <- as.factor(data$diab)
data$gender <- as.factor(data$gender)

# Load required libraries
library(caret)

# Detach the data frame to avoid conflicts
detach(data)

# Attach the data frame for easy access to variables
attach(data)

# Display a summary of the data
summary(data)

# Create a plot of the data
plot(data)

# Display the number of rows in the data
nrow(data)

# Create a contingency table for 'diab' and 'gender'
table(diab, gender)  # relatively more males are diabetic

# Set up for subplot creation (not executed)
# par(mfcol=c(1,2))   # used to create subplots

# Create a scatter plot of age vs. waist, colored by 'diab' and shaped by 'gender'
plot(age ~ waist, col = diab, pch = ifelse(gender == 1, 2, 20))  # circles female (gender=2)

# Create a scatter plot of height vs. weight, colored by 'diab' and shaped by 'gender'
plot(height ~ weight, col = diab, pch = ifelse(gender == 1, 2, 20))  # circles female (gender=2)

# Create a scatter plot of waist vs. age and add a diagonal line
plot(waist ~ age)
abline(a = 0, b = 1, col = 2)

# Fit a linear model of height on gender and weight
fit <- lm(height ~ gender + weight)
summary(fit)

# Create a scatter plot of BMI vs. diab
plot(BMI ~ diab)

# Create a scatter plot of age vs. diab
plot(age ~ diab)

# Outlier example: Scatter plot of weight vs. height for a subset of data
plot(weight ~ height, data = data[100:150, ])
# Highlight specific points in red
points(weight ~ height, data = data[c(147, 120), ], pch = 20, col = "red")

# Fit a linear model for weight on height using data between 100 and 150
fit <- lm(weight ~ height, data = data[100:150, ])
summary(fit)
abline(fit)
anova(fit)
plot(fit)

# Fit a linear model for weight on height excluding red points (outliers)
fit2 <- lm(weight ~ height, data = data[c(100:119, 121:146, 148:150), ])
summary(fit2)
abline(fit2, col = "blue")
anova(fit2)
plot(fit2)

# Blood pressure: Scatter plot of systolic blood pressure vs. weight, colored by gender
plot(sysBP ~ weight, col = gender, pch = 20, xlim = c(30, 210), ylab = "systolic blood pressure")

# Fit a linear model for systolic blood pressure on weight and the interaction term between weight and gender
BPfit <- lm(sysBP ~ weight + weight:gender, data = data)
summary(BPfit)
abline(BPfit)
abline(a = coef(BPfit)[1], b = coef(BPfit)[2] + coef(BPfit)[3], col = 2)
plot(BPfit)

# Linear model for systolic blood pressure on gender, weight, and height
bpfit1 <- lm(sysBP ~ gender + weight + height, data = data)
summary(bpfit1)

# Perform LOOCV using linear model 'bpfit1'
controlvar <- trainControl(method = "LOOCV")
model1 <- train(formula(bpfit1),
                data = data[complete.cases(data), ],
                method = "lm",
                trControl = controlvar)
print(model1)
summary(model1)

# Linear model for systolic blood pressure on the interaction terms between weight and gender, and height
bpfit2 <- lm(sysBP ~ weight * gender + height)
summary(bpfit2)

# Perform LOOCV using linear model 'bpfit2'
model2 <- train(formula(bpfit2),
                data = data[complete.cases(data), ],
                method = "lm",
                trControl = controlvar)
print(model2)
summary(model2)

# Linear model for systolic blood pressure on weight, gender, and height excluding intercept
bpfit1.1 <- lm(sysBP ~ weight + gender + height - 1)
summary(bpfit1.1)

# Perform LOOCV using linear model 'bpfit1.1'
model1.1 <- train(formula(bpfit1.1),
                  data = data[complete.cases(data), ],
                  method = "lm",
                  trControl = controlvar)
print(model1.1)
summary(model1.1)

# Logistic models for diabetes prediction
logistic <- glm(diab ~ .,
                binomial(link = "logit"),
                data = data[complete.cases(data), ])
summary(logistic)
anova(logistic, test = "Chisq")

logistic0 <- glm(diab ~ gender + age + waist + height,
                 binomial(link = "logit"),
                 data = data[complete.cases(data), ])
summary(logistic0)
anova(logistic0, logistic, test = "Chisq")
