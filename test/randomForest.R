library(dplyr)
data_train <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/train.csv")
data_train <- data_train[complete.cases(data_train), ]
data_train$Survived <- factor(data_train$Survived)
glimpse(data_train)
data_test <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/test.csv") 
data_test <- data_test[complete.cases(data_test), ]
glimpse(data_test)

library(randomForest)
library(caret)
library(e1071)

# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")

# Run the model
set.seed(716)
rf_default <- train(Survived~.,
                    data = data_train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)

### simpler iris test
data(iris)
TrainData <- iris[,1:4]
TrainClasses <- iris[,5]

knnFit1 <- train(TrainData, TrainClasses,
                 method = "knn",
                 preProcess = c("center", "scale"),
                 tuneLength = 10,
                 trControl = trainControl(method = "cv"))









