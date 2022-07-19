
{
  .libPaths( c( "/home/rodrig/docker/" , .libPaths() ) )
  packlist <- c( "parallel", "foreach", "doParallel",
                 "dplyr", "tidyr", "purrr", "caTools", "doSNOW", "data.table", "caret")

  package.check <- lapply(
    packlist,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )

  rm(package.check, packlist)

}

pathIn = "/home/rodrig/repository/"

maxDistYear=readRDS(paste0(pathIn, "processing_results/full_LTS_ICP_Data.RDS"))

cNames= colnames(maxDistYear)

selNames = do.call(c,lapply(c("time","plot_id", "Defoliation", "delta", "chn", "var", "acf",
                              "tp", "t2m", "w", "ssr", "dem", "asp","slp", "seg", "len", "Dist"), function(i)
                                Filter(function(x) grepl(i, x), cNames)))

maxDistYear = as(maxDistYear, "data.frame")

#maxDistYear = maxDistYear[, selNames]

varSelection = lapply(seq(5,30,5), function(x){

  trainingDF = maxDistYear[maxDistYear$TreeCount >=x,selNames[c(6:40,42,44)]]#; str(trainingDF)

  trainingDF=trainingDF[!trainingDF$genDist %in% c("NoData"),]

  stable = maxDistYear[maxDistYear$TreeCount <5,selNames[c(6:40,42,44)]]

  stable=stable[!stable$genDist %in% c("NoData"),]

  stable$genDist = "LowDisturbance"

  trainingDF= rbind(trainingDF,stable)


  trainingDF$dummy[trainingDF$genDist=="Abiotic"] <- 1
  trainingDF$dummy[trainingDF$genDist=="Biotic"] <- 2
  trainingDF$dummy[trainingDF$genDist=="Management"] <- 3
  trainingDF$dummy[trainingDF$genDist=="Others"] <- 4
  trainingDF$dummy[trainingDF$genDist=="LowDisturbance"] <- 5

  model= lm(dummy~., data = trainingDF[,-c(length(trainingDF)-1)])#;  summary(model)
  vif = olsrr::ols_vif_tol(model)
  # corr = olsrr::ols_correlations(model)
  var=vif[vif$VIF <=4,1]

  var= data.frame(variables = var, model = as.factor(x))


})


varSelection = do.call(rbind, varSelection)

varSelection = unique(varSelection$variables)


newData = lapply(seq(5,30, 5), function(k){

  trainingDF = maxDistYear[maxDistYear$TreeCount >=k,c(varSelection, "genDist")]

  trainingDF=trainingDF[!trainingDF$genDist %in% c("NoData"),]

  stable = maxDistYear[maxDistYear$TreeCount <5,c(varSelection, "genDist")]

  stable=stable[!stable$genDist %in% c("NoData"),]

  stable$genDist = "LowDisturbance"

  trainingDF= rbind(trainingDF,stable)

  trainingDF0 = trainingDF



  bio = trainingDF0[trainingDF0$genDist =="Biotic",]

  bio$id = 1:nrow(bio)

  n = length(trainingDF0[trainingDF0$genDist =="Abiotic",1])

  #set.seed(123456)
  bio0= lapply(1:100, function(x) bio[sample(bio$id, size = n),-length(bio)])

  newData= bio0 %>% map(~rbind(trainingDF0[trainingDF0$genDist !="Biotic",], .x))

  inTrain = newData %>% map(~createDataPartition(.x$genDist, p = 0.60)[[1]])
  training = map2(newData, inTrain, ~.x[.y,])
  rest = map2(newData, inTrain, ~.x[-.y,])
  inTest = rest %>% map(~createDataPartition(.x$genDist, p = 0.50)[[1]])
  testing=map2(rest, inTest, ~.x[.y,])
  validation=map2(rest, inTest, ~.x[-.y,])

  list(training, testing, validation)

})

saveRDS(newData, paste0(pathIn,"processing_results/LTS_ICP_training_resample_100_seq_5_30_trees_all_lowDisturbance.RDS"))



modelTraining_treeThreshold=lapply(1:6, function(k){

  data= newData[[k]][[1]] %>% map(~na.omit(.x))

  validation = newData[[k]][[2]] %>% map(~na.omit(.x))

  iterations <- 100
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cl=parallel::makeCluster(detectCores()-2)
  doSNOW::registerDoSNOW(cl)

  trControl <- trainControl(method = "cv",
                            number = 10,
                            search = "grid")





  clusterEvalQ(cl,  .libPaths() )



  modelTraining =foreach(i=1:iterations, .packages = c("caret", "xgboost", "kernlab",  "randomForest"), .options.snow = opts,.verbose = T)%dopar%{

    #set.seed(1234)
    rf_default <- train(genDist ~.,
                        data = data[[i]],
                        method = "rf",
                        metric = "Accuracy",
                        trControl = trControl)


    svm_default <- train(genDist ~.,
                         data = data[[i]],
                         method = "svmRadial",
                         metric = "Accuracy",
                         trControl = trControl)


    knn_default <- train(genDist ~.,
                         data = data[[i]],
                         method = "kknn",
                         metric = "Accuracy",
                         trControl = trControl)


    gps_default <- train(genDist ~.,
                         data =  data[[i]],
                         method = "gaussprRadial",
                         metric = "Accuracy",
                         trControl = trControl)

    xb_default <- train(genDist ~.,
                        data =  data[[i]],
                        method = "xgbTree",
                        metric = "Accuracy",
                        trControl = trControl)

    predRF <- predict(rf_default, newdata = validation[[i]])


    predSVM <- predict(svm_default, newdata = validation[[i]])



    predKNN <- predict(svm_default, newdata = validation[[i]])

    predGPS <- predict(gps_default, newdata = validation[[i]])

    predXGB <- predict(xb_default, newdata = validation[[i]])

    validation[[i]]$genDist = as.factor(validation[[i]]$genDist)

    cmRF =(confusionMatrix(predRF, validation[[i]]$genDist))

    cmSVM = (confusionMatrix(predSVM, validation[[i]]$genDist))

    cmKKNN = confusionMatrix(predKNN, validation[[i]]$genDist)

    cmGPS = confusionMatrix(predGPS, validation[[i]]$genDist)

    cmXGB = confusionMatrix(predXGB, validation[[i]]$genDist)

    list(print(rf_default), print(svm_default),
         print(knn_default), print(gps_default), print(xb_default),
         rf_default, svm_default, knn_default, gps_default,xb_default,cbind(cmRF$overall[1], cmRF$overall[2]),
         cbind(cmSVM$overall[1], cmSVM$overall[2]), cbind(cmKKNN$overall[1], cmKKNN$overall[2]), cbind(cmGPS$overall[1], cmGPS$overall[2]),
         cbind(cmXGBS$overall[1], cmXGB$overall[2])
    )

  }


  stopCluster(cl)
  rm(cl)

  {
    aRf = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[1]])))#9
    aSvm = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[2]])))#11
    aKnn = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[3]])))#13
    aGPS =  data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[4]])))
    aXGB =  data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[5]])))

    aRf$model = "RF"
    aSvm$model= "SVM"
    aKnn$model="KNN"
    aGPS$model ="GPS"
    aXGB$model ="XGB"


    #


    mAcu = bind_rows(aRf, aSvm,  aKnn,aGPS, aXGB)

    mAcu$model = as.factor(mAcu$model)


    mAcu$Accuracy = as.numeric(mAcu$Accuracy)

    mAcu$Kappa = as.numeric(mAcu$Kappa)


    write.table(mAcu, paste0(pathIn,"processing_results/modelAccu_training_CV_", seq(5,30, 5)[k], ".txt"), quote = F, row.names = F, sep = "\t")
  }


  {
    aRf = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[11]])))#9
    aSvm = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[12]])))#11
    aKnn = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[13]])))#13
    aGPS =  data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[14]])))
    aXGB =  data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[15]])))

    aRf$model = "RF"
    aSvm$model= "SVM"
    aKnn$model="KNN"
    aGPS$model ="GPS"
    aXGB$model ="XGB"


    #


    mAcu = bind_rows(aRf, aSvm,  aKnn,aGPS, aXGB)

    mAcu$model = as.factor(mAcu$model)


    mAcu$Accuracy = as.numeric(mAcu$Accuracy)

    mAcu$Kappa = as.numeric(mAcu$Kappa)


    write.table(mAcu, paste0(pathIn,"processing_results/modelAccu_validation_CV_0", seq(5,30, 5)[k], ".txt"), quote = F, row.names = F, sep = "\t")



  }

  #boxplot(Accuracy ~ model, mAcu, main=paste0("ICP-Forest Data Classification Accuracy with ", length(var), " variables"))

  saveRDS(modelTraining, file =paste0(pathIn,"processing_results/trained_model_0", seq(5,30, 5)[k], "_trees.RDS"))

})

unregister_dopar()
accu =  list.files("Z:/trained_LTS_to_ICP_4_classes/", pattern = "validation_CV_0", full.names = T) %>% map(~read.csv(.x, sep="\t")) %>%
  map(~setnames(.x, c("Accuracy", "Kappa", "model")))

accu =accu %>% map(~mutate(.x, model = as.factor(.x$model)))

accu =accu %>% map(~mutate(.x, Accuracy = as.numeric(.x$Accuracy)))

accu %>% map(~boxplot(Accuracy ~ model, .x, main=paste0("ICP-Forest Data Classification Accuracy with ", length(var), " variables")))

accu = lapply(1:6, function(x) mutate( accu[[x]], sampled_trees=seq(5,30, 5)[x]))

accu= do.call(rbind, accu)

library(reshape2)
library(ggplot2)

# Assume your data frame is named dat
dat.m = melt(accu, id.var=c("model","sampled_trees"))
dat.m$model = as.factor(dat.m$model)
dat.m$sampled_trees = as.factor(dat.m$sampled_trees)
# If you want the two levels of event plotted side by side
ggplot(dat.m, aes(sampled_trees, value, colour=model)) +
  facet_grid(. ~ variable) +
  geom_boxplot(width=0.7)

modelList = gsub(".RDS", "", list.files("D:/trained_LTS_to_ICP_4_classes/", full.names = F, pattern = "trained_model_0")[c(20,15:19)])



h=6

modelTraining = readRDS(paste0("D:/trained_LTS_to_ICP_4_classes/", modelList[[h]], ".RDS"))
iterations=100


{

  aRf = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[7]])))#9
  aSvm = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[8]])))#11
  aKnn = data.frame(do.call(rbind, lapply(1:iterations, function (x) modelTraining[[x]][[9]])))#13

  aRf$model = "RF"
  aSvm$model= "SVM"
  aKnn$model="KNN"

  #


  mAcu = bind_rows(aRf, aSvm,  aKnn)

  colnames(mAcu) = c("Accuracy", "Kappa", "model")

  mAcu$model = as.factor(mAcu$model)


  mAcu$Accuracy = as.numeric(mAcu$Accuracy)

  mAcu$Kappa = as.numeric(mAcu$Kappa)

  #b=boxplot(Accuracy ~ model, mAcu, main=paste0("ICP-Forest Data Classification Accuracy with ", length(var), " variables"))


}
b= b$stats

colnames(b) =c("KNN", "RF", "SVM")

print(b)


{

  names(newData) = modelList
  bestM = c(ceiling(which.max(aRf[,1])/3), #ceiling(which.max(aPls[,3])/3),
            ceiling(which.max(aSvm[,1])/3), #ceiling(which.max(aNb[,3])/2),
            ceiling(which.max(aKnn[,1])/3))


  models = modelTraining#_treeThreshold[[6]][[4]]#[[18]]


  bestRF= models[[bestM[1]]][[4]]

  #bestPLS = models[[ceiling(which.max(aPls[,3])/3)]][[7]]

  bestSVM = models[[bestM[2]]][[5]]
  #bestNB = models[[ceiling(ceiling(which.max(aNb[,3])/2))]][[9]]
  bestKNN = models[[bestM[3]]][[6]]


  newData0= bestM %>% map(~newData[[modelList[[h]]]][[1]][[.x]])

  newData0 = newData0 %>% map(~mutate(.x, id = as.numeric(rownames(.x))))

  newData0 = do.call(rbind, newData0)

  training = na.omit(newData0[!duplicated(newData0$id),-length(newData0)])

  newData1= bestM %>% map(~newData[[modelList[[h]]]][[2]][[.x]])

  newData1 = newData1 %>% map(~mutate(.x, id = as.numeric(rownames(.x))))

  newData1 = do.call(rbind, newData1)

  testing = na.omit(newData1[!duplicated(newData1$id),-length(newData1)])

  newData2= bestM %>% map(~newData[[modelList[[h]]]][[3]][[.x]])

  newData2 = newData2 %>% map(~mutate(.x, id = as.numeric(rownames(.x))))

  newData2 = do.call(rbind, newData2)

  validation = na.omit(newData2[!duplicated(newData2$id),-length(newData2)])

}
{
  trControl <- trainControl(method = "cv",
                            number = 10,
                            search = "grid")
  #,allowParallel=T)

  # rf_default <- train(genDist ~.,
  #                     data = training,
  #                     method = "rf",
  #                     metric = "Accuracy",
  #                     trControl = trControl)
  #
  # # pls_default <- train(genDist ~.,
  # #                      data = training,
  # #                      method = "kernelpls",
  # #                      metric = "Accuracy",
  # #                      trControl = trControl)
  #
  # svm_default <- train(genDist ~.,
  #                      data = training,
  #                      method = "svmRadial",
  #                      metric = "Accuracy",
  #                      trControl = trControl)
  #
  # # nb_default <- train(genDist ~.,
  # #                     data = training,
  # #                     method = "naive_bayes",
  # #                     metric = "Accuracy",
  # #                     trControl = trControl)
  #
  # knn_default <- train(genDist ~.,
  #                      data = training,
  #                      method = "kknn",
  #                      metric = "Accuracy",
  #                      trControl = trControl)

  predRF <- predict(bestRF, newdata = validation)

  #predPLS <- predict(pls_default, newdata = validation)

  predSVM <- predict(bestSVM, newdata = validation)

  #predNB <- predict(nb_default, newdata = validation)

  predKNN <- predict(bestKNN, newdata = validation)


  validation$genDist = as.factor(validation$genDist)

  print(confusionMatrix(predRF, validation$genDist))
  #confusionMatrix(predPLS, validation$genDist)
  print(confusionMatrix(predSVM, validation$genDist))
  #confusionMatrix(predNB, validation$genDist)
  print(confusionMatrix(predKNN, validation$genDist))


  predDF <- data.frame(predRF,  predSVM,  predKNN, genDist = validation$genDist, stringsAsFactors = F)
}
{# Train the ensemble
  modelStack <- train(genDist ~ ., data = predDF, method = "gbm")

  tibble::as_tibble(summary(modelStack))

  # Generate predictions on the test set
  predRF <- predict(bestRF, newdata = testing)

  #predPLS <- predict(pls_default, newdata = testing)

  predSVM <- predict(bestSVM, newdata = testing)

  #predNB <- predict(nb_default, newdata = testing)

  predKNN <- predict(bestKNN, newdata = testing)


  # Using the base learner test set predictions,
  # create the level-one dataset to feed to the ensemble
  testPredLevelOne <- data.frame(predRF,  predSVM,  predKNN, genDist = testing$genDist, stringsAsFactors = F)
  combPred <- predict(modelStack, newdata=testPredLevelOne)

  testing$genDist = as.factor(testing$genDist)

  # Evaluate ensemble test performance
  confusionMatrix(combPred, testing$genDist)#$overall[1]
}
# Evaluate base learner test performance
confusionMatrix(testPredRF, testing$diagnosis)$overall[1]
confusionMatrix(testPredGBM, testing$diagnosis)$overall[1]
confusionMatrix(testPredLDA, testing$diagnosis)$overall[1]

















trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid",
                          savePredictions="final",
                          classProbs=F,
                          index= createResample(newData[[31]]$genDist, 10)
                          ,allowParallel=T)



model_list = caretList(genDist ~.,
                       data=newData[[31]],
                       trControl = trControl,
                       methodList = c("rf", "svmRadial", "naive_bayes", "kknn", "gaussprRadial"))

pre=as.data.frame(predict(model_list, newData[[1]]))

gn = as.factor(newData[[1]]$genDist)

prRF = as.factor(pre$gaussprRadial)

confusionMatrix(gn, prRF)
table(newData[[1]]$genDist, pre$svmRadial)
table(newData[[1]]$genDist, pre$naive_bayes)


xyplot(resamples(model_list))

modelCor(resamples(model_list))

greedy_ensemble <- caretEnsemble(
  model_list,
  metric="Accuracy",
  trControl=trainControl(
    number=10,

    classProbs=TRUE
  ))
summary(greedy_ensemble)




