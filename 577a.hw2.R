

library(DNAshapeR)
library(caret)
library(ROCR)
library(Biostrings)

#Q4
fn1="/Users/luna/Documents/577a/Mad.txt.fa"
fn2="/Users/luna/Documents/577a/Max.txt.fa"
fn3="/Users/luna/Documents/577a/Myc.txt.fa"

featureType1 <- c("1-mer") 
featureType2 <- c("1-mer", "1-shape") 

featureVector1 <- encodeSeqShape(fn1, getShape(fn1) , featureType1)
featureVector2 <- encodeSeqShape(fn1, getShape(fn1), featureType2)
experimentalData_mad <- read.table("/Users/luna/Documents/577a/Mad.txt") 
df_mad1 <- data.frame(affinity=experimentalData_mad$V2, featureVector1)
df_mad2 <- data.frame(affinity=experimentalData_mad$V2, featureVector2)


trainControl<- trainControl(method = "cv", number = 10, savePredictions = TRUE) 

model_mad1 <- train (affinity~ ., data = df_mad1, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_mad2 <- train (affinity~ ., data = df_mad2, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

model_mad1$result$Rsquared[1]
model_mad2$result$Rsquared[1]


featureVector1 <- encodeSeqShape(fn2, getShape(fn2) , featureType1)
featureVector2 <- encodeSeqShape(fn2, getShape(fn2), featureType2)
experimentalData_max<- read.table("/Users/luna/Documents/577a/Max.txt") 
df_max1 <- data.frame(affinity=experimentalData_max$V2, featureVector1)
df_max2 <- data.frame(affinity=experimentalData_max$V2, featureVector2)


model_max1 <- train (affinity~ ., data = df_max1, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_max2 <- train (affinity~ ., data = df_max2, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

model_max1$result$Rsquared[1]
model_max2$result$Rsquared[1]


featureVector1 <- encodeSeqShape(fn3, getShape(fn3) , featureType1)
featureVector2 <- encodeSeqShape(fn3, getShape(fn3), featureType2)
experimentalData_myc<- read.table("/Users/luna/Documents/577a/Myc.txt") 
df_myc1 <- data.frame(affinity=experimentalData_myc$V2, featureVector1)
df_myc2 <- data.frame(affinity=experimentalData_myc$V2, featureVector2)

model_myc1 <- train (affinity~ ., data = df_myc1, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model_myc2 <- train (affinity~ ., data = df_myc2, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

model_myc1$result$Rsquared[1]
model_myc2$result$Rsquared[1]

#Q5
x<-c(model_mad1$result$Rsquared[1], model_max1$result$Rsquared[1], model_myc1$result$Rsquared[1])
y<-c(model_mad2$result$Rsquared[1], model_max2$result$Rsquared[1], model_myc2$result$Rsquared[1])
att<-c(1,2,3)
plot(x,y, pch=as.integer(att),ylim=c(0.75,1), xlim=c(0.75,1),xlab= expression('1-mer R' ^ 2), ylab= expression('1-mer+shape R' ^ 2))
abline(a = 0, b=1)


#Q7
bound="/Users/luna/Documents/577a/bound_500.fa"
unbound="/Users/luna/Documents/577a/unbound_500.fa"
pred1 <- getShape(bound)
pred2 <- getShape(unbound)

plotShape(pred1$MGW, main="Bound500 MGW")
plotShape(pred2$MGW, main="Unbound500 MGW")

plotShape(pred1$ProT, main="Bound500 ProT")
plotShape(pred2$ProT, main="Unbound500 ProT")

plotShape(pred1$Roll, main="Bound500 Roll")
plotShape(pred2$Roll, main="Unbound500 Roll")

plotShape(pred1$HelT, main="Bound500 HelT")
plotShape(pred2$HelT, main="Unbound500 HelT")


#Q8
workingPath<- "/Users/luna/Documents/577a/"

#bound
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape(paste0(workingPath, "ctcf.fa"))


## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)

featureType2 <- c("1-mer")
featureVector2 <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType2)
df2 <- data.frame(isBound = exp_data$isBound, featureVector2)

## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~., data = df, trControl = trainControl, method = "glm", family = binomial, metric ="ROC")
summary(model)

model2 <- train(isBound~., data = df2, trControl = trainControl, method = "glm", family = binomial, metric ="ROC")
summary(model2)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance, main="ROC for feature+shape")


prediction2 <- prediction( model2$pred$Y, model2$pred$obs )
performance2 <- performance( prediction2, "tpr", "fpr" )
plot(performance2, main="ROC for feature only")

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
text(x=0.5, y = 0.5, labels = paste('AUC=', round(auc,3)) )

auc2 <- performance(prediction2, "auc")
auc2 <- unlist(slot(auc2, "y.values"))
text(x=0.5, y = 0.5, labels = paste('AUC=', round(auc2,3)) )




