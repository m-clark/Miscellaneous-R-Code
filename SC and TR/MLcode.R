#-------------------------------------------------------------------------------#
# See http://www.nd.edu/~mclark19/learn/ML.pdf for the handout on machine       #  
# learning this code regards.                                                   # 
#-------------------------------------------------------------------------------#

sqerrloss = function(beta, X, y){
  mu = X%*%beta
  sum((y-mu)^2)
}

set.seed(123)
X = cbind(1, rnorm(100), rnorm(100))
y = rowSums(X[,-1] + rnorm(100))
out1 = optim(par=c(0,0,0), fn=sqerrloss, X=X, y=y)
out2 = lm(y ~ X[,2] + X[,3])  # check with lm 
rbind(c(out1$par, out1$value), c(coef(out2),sum(resid(out2)^2)) )




sqerrloss_reg = function(beta, X, y, lambda=.1){
  mu = X%*%beta
  sum((y-mu)^2) + lambda*sum(abs(beta[-1]))
}

out3 = optim(par=c(0,0,0), fn=sqerrloss_reg, X=X, y=y)
rbind(c(out1$par, out1$value), 
      c(coef(out2),sum(resid(out2)^2)), 
      c(out3$par, out3$value) )






wine = read.csv('http://www.nd.edu/~mclark19/learn/data/goodwine.csv')
summary(wine)



library(doSNOW)
registerDoSNOW(makeCluster(11, type = "SOCK"))



library(corrplot)
corrplot(cor(wine[,-c(13,15)]), method="number", tl.cex=.5)



library(caret)
set.seed(1234) #so that the indices will be the same when re-run
trainIndices = createDataPartition(wine$good, p=.8, list=F)
wanted = !colnames(wine) %in%  c("free.sulfur.dioxide", "density", "quality", "color", "white" )
wine_train = wine[trainIndices, wanted]  #remove quality and color, as well as density and others

wine_test = wine[-trainIndices, wanted]
# prep_test = preProcess(wine_test[,-10], method="range")
# wine_test = data.frame(predict(prep_test, wine_test[,-10]), good=wine_test[ ,10])



wine_trainplot = predict(preProcess(wine_train[,-10], method="range"), 
                         wine_train[,-10])
featurePlot(wine_trainplot, wine_train$good, "box")



set.seed(1234)
cv_opts = trainControl(method="cv", number=10)
knn_opts = data.frame(.k=c(seq(3, 11, 2), 25, 51, 101)) #odd to avoid ties
results_knn = train(good~., data=wine_train, method="knn", 
                    preProcess="range", trControl=cv_opts, 
                    tuneGrid = knn_opts)

results_knn



preds_knn = predict(results_knn, wine_test[,-10])
confusionMatrix(preds_knn, wine_test[,10], positive='Good')
conf_knn = confusionMatrix(preds_knn, wine_test[,10], positive='Good') #create an object to use in Sexpr



dotPlot(varImp(results_knn))



set.seed(1234)
results_nnet = train(good~., data=wine_train, method="avNNet", 
                     trControl=cv_opts, preProcess="range", 
                     tuneLength=5, trace=F, maxit=1000)
results_nnet



preds_nnet = predict(results_nnet, wine_test[,-10])
confusionMatrix(preds_nnet, wine_test[,10], positive='Good')
conf_nnet = confusionMatrix(preds_nnet, wine_test[,10], positive='Good') #create an object to use in Sexpr



set.seed(1234)
rf_opts = data.frame(.mtry=c(2:6))
results_rf = train(good~., data=wine_train, method="rf", 
                   preProcess='range',trControl=cv_opts, tuneGrid=rf_opts, 
                   n.tree=1000)
results_rf



preds_rf = predict(results_rf, wine_test[,-10])
confusionMatrix(preds_rf, wine_test[,10], positive='Good')
conf_rf = confusionMatrix(preds_rf, wine_test[,10], positive='Good') #create an object to use in Sexpr



set.seed(1234)
results_svm = train(good~., data=wine_train, method="svmLinear", 
                    preProcess="range", trControl=cv_opts, tuneLength=5)
results_svm



preds_svm = predict(results_svm, wine_test[,-10])
confusionMatrix(preds_svm, wine_test[,10], positive='Good')
conf_svm = confusionMatrix(preds_svm, wine_test[,10], positive='Good') #create an object to use in Sexpr


stopCluster()