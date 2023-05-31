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
library(randomForest)
library(tree)
library(ggplot2)
library(GGally)
library(dplyr)

index_row <- sample(2,
                    nrow(iris), 
                    replace = T, 
                    prob = c(0.7, 0.3)
                    )

train_data <- iris[index_row == 1, ]
test_data <- iris[index_row == 2, ]

iris_classifier <- randomForest(Species ~., 
                                data = train_data, #train data set 
                                importance = T)

plot(iris_classifier)

importance(iris_classifier)   #Petal features are more important

qplot(Petal.Width, Petal.Length, data = iris, color = Species)

predicted_table <- predict(iris_classifier, test_data[,-5])
table(observed = test_data[,5], predicted = predicted_table)

