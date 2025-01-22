library(Seurat)
library(Boruta)
library(randomForest)
library(dplyr)
library(caret)
library(ggRandomForests)

setwd("~/GitHub/MACanalyzeR_MaterialsMethods/FoamSpotteR/Z_RandomForest/")

###################################
## SETTING DATA
###################################
mac <- readRDS("../A_DatasetSC_ATHEROSPECTRUM/mac.rds")
DimPlot(mac)

## count matrix
mac_expr <- as.data.frame(t(as.data.frame(mac@assays$RNA@data)))
mac_expr$Sample <- factor(mac$Cells, levels = c("fMAC-", "fMAC+"))
table(mac_expr$Sample)

##  training and test set 
training_set <- mac_expr %>% slice_sample(prop=0.75, replace = F)
test_set <- anti_join(mac_expr, training_set)

nrow(training_set) + nrow(test_set) == nrow(mac_expr)

#############################
## boruta
#############################
set.seed(20)
boruta <- Boruta(Sample ~ ., data = training_set, doTrace = 2, maxRuns = 500)
print(boruta)

boruta$finalDecision

att <- attStats(boruta)
att <- att[which(att$decision=="Confirmed"),]
att <- att[order(att$meanImp, decreasing = T),]

mac_gene <- training_set[,c(rownames(att[1:20,]), "Sample")]

##############################
##  K fold cross validation
##############################
# 5 fold - means that 20% of the data is used for testing (pretty accurate)
repeat_cv <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
tunegrid <- expand.grid(.mtry = as.integer(sqrt(20))) 

cv <- train(Sample ~. , 
           data=mac_gene,
           method="rf",
           trControl=repeat_cv,
           tuneGrid=tunegrid,
           metric="Accuracy"
)

cv$finalModel$confusion

#############################
## random forest
#############################
rf <- randomForest(Sample~. -1,  # response
                   data=mac_gene, #predictors
                   ntree=500, #total num of trees in forest
                   mtry= 4,   #mtry
                   nodesize= 5, #nodesize
                   importance=T, #computing the importance of predictive variables? Y/N
                   keep.forest=T,  #retaining the forest in the output object? Y/N
)

rf
summary(rf)
str(rf)

#################
##  validation
#################
# PREDICT THE TEST SET DATA
pred_test <- as.data.frame(predict(rf, newdata=test_set[,rownames(att[1:20,])], , type = "prob"))

pred_test$Sample <- "fMAC-"
pred_test[which(pred_test$`fMAC+`>0.5), "Sample"] <- "fMAC+"


confmatrix <- table(Actual_value=test_set$Sample, Predicted_Value = pred_test$Sample)
confmatrix

confusionMatrix(as.factor(pred_test$Sample),test_set$Sample)


gg_dta <- gg_roc(rf, which_outcome = 2)
roc <- plot(gg_dta) +
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        legend.key=element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black"))


gg_dta<- gg_error(rf)
colnames(gg_dta) <- c("OOB", "fMAC-", "fMAC+", "ntree")
err <- plot(gg_dta) +
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black"))

roc | err
