# Clear environment
##############
rm(list=ls())
gc()
########

# Import libraries
##############
library(dplyr)
library(ROCR)
library(ggplot2)
library(gbm)
########

# Import data
##############
setwd("C:/Users/bhavi/Desktop/Research")

# Import data
setwd("C:/Users/bhavi/Desktop/Research")

# Import data
cells <- read.csv("GSE71585_RefSeq_TPM.csv")
gene <- cells$gene
# Transpose and change to data frame
cells <- data.frame(t(data.matrix(cells[,-1])))
colnames(cells) <- gene

class1 <- cells[703:781,]
class2 <- cells[782:838,]
########

# Data Pre-processing
##############
# Shuffle rows
set.seed(100)
class1<- class1[sample(nrow(class1)),]
set.seed(100)
class2<- class2[sample(nrow(class2)),]

class1['y'] <- 1
class2['y'] <- 0
l1 <- nrow(class1)
l2 <- nrow(class2)

class_all <- rbind(class1,class2)

# Remove genes not expressed or equally expressed in all cells
removeZeroVar <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}
data <- removeZeroVar(class_all)

data['y'] <- class_all['y']

gc()
########

# Pre-allocate variables
##############
df_genes <- data.frame(genes=colnames(data[,-ncol(data)]), beta = numeric(ncol(data)-1))
df_genes$beta <- 0
row.names(df_genes) <- df_genes$genes

k=1
j=l1+1
nfolds = 3
l1fold = floor(l1/nfolds)
l2fold = floor(l2/nfolds)


# create parameter list
params <- list(
  eta = .1,
  max_depth = 5,
  min_child_weight = 2,
  subsample = .8,
  colsample_bytree = .9
)
#########



# K-Fold Cross validation starts here
##############
for (i in 1:nfolds){
  
  if(i==nfolds){
    test_rows <- c(k:l1,j:(l1+l2))
  } else {
    test_rows <- c(k:(k+l1fold-1),j:(j+l2fold-1))
  }
  
  test_data <- data[test_rows,]
  train_data <- data[-test_rows,]
  test_data$y <- as.factor(test_data$y)
  train_data$y <- as.factor(train_data$y)
  
  k=k+l1fold
  j=j+l2fold
  
  
# gbm
##############
  start.time <- Sys.time()
  set.seed(100)
  gbm.fit <- gbm( formula = y~.
                  , data = train_data
                  , distribution="bernoulli"
                  , n.trees =50
                  , interaction.depth = 4
                  , shrinkage = 0.1
                  , bag.fraction = .65 
                  , train.fraction = 1
                  , cv.folds = 5
                  , verbose = FALSE
  )
  end.time <- Sys.time()
  comp.time <- end.time - start.time
##################  
  
# AUC
#######
  gbm.probabilities <- predict(gbm.fit, newdata = test_data)
  gbm.y_pred <- ifelse(gbm.probabilities > 0.5, 1, 0)
  gbm.pred <- prediction(gbm.y_pred, y_test)
  gbm.auc.perf <- performance(gbm.pred, measure = "auc")
  gbm.val[i] <- gbm.auc.perf@y.values
#######

}
#########


# get number of trees that minimize error
#######
gbm.fit$evaluation_log %>%
  dplyr::summarise(
    ntrees.train = which(train_logloss_mean == min(train_logloss_mean))[1],
    rmse.train   = min(train_logloss_mean),
    ntrees.test  = which(test_logloss_mean == min(test_logloss_mean))[1],
    rmse.test   = min(test_logloss_mean),
  )
#######

# create hyperparameter grid
###########
hyper_grid <- expand.grid(
  eta = c(.01, .05, .1, .3),
  max_depth = c(1, 3, 5, 7),
  min_child_weight = c(1, 3, 5, 7),
  subsample = c(.65, .8, 1), 
  colsample_bytree = c(.8, .9, 1),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)
#############

# grid search
###########
for(i in 1:nrow(hyper_grid)) {
  
  # create parameter list
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i],
    min_child_weight = hyper_grid$min_child_weight[i],
    subsample = hyper_grid$subsample[i],
    colsample_bytree = hyper_grid$colsample_bytree[i]
  )
  
  # reproducibility
  set.seed(100)
  
  # train model
  xgb.tune <- xgb.cv(
    params = params,
    data = x_train,
    label = y_train,
    nrounds = 1000,
    nfold = 10,
    objective = "binary:logistic",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_logloss_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_logloss_mean)
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)
##########