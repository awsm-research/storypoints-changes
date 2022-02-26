# This file measures the performance of our classifiers and save to /data/rq4_performance 
#For each bootstrap dataset, we first load the data and the pre-built model (from RQ4.r).

#We first determine the threshold (cutoff; getCutOff()) for threshold dependent metrics (i.e., Balanced Accuracy and F1-Score).
#Then, we use the classifiers to predict whether a work item in testing dataset will have an SP change.
#After that, we calculate AUC, BACC, F1-score, and save prediction performance.
#Prediction performance of for each bootstrap is saved to data/RQ4_performance/[PROJECT]_LIST.csv
#Average prediction performance is saved to data/RQ4_performance/[PROJECT].csv

########################################################################################
# SETUP
########################################################################################
this_file = gsub("--file=", "", commandArgs()[grepl("--file", commandArgs())])
if (length(this_file) > 0){
  wd <- paste(head(strsplit(this_file, '[/|\\]')[[1]], -1), collapse = .Platform$file.sep)
}else{
  wd <- dirname(rstudioapi::getSourceEditorContext()$path)
}
setwd(wd)
source("util.R")
source("RQ4_util.R")
########################################################################################


library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fpc)
library(Rnalytica)
library(rms)
library(randomForest)
library(e1071)
library(car)
library(ROCR)
library(caret)
library(OneR)
library(data.table)
library(OptimalCutpoints)

# determine cutoff based on work items in training dataset.
getCutOff = function(data) {
  spChanged = data[data$y == TRUE,]
  spNotChanged = data[data$y == FALSE,]
  
  spChangedFirstQ = quantile(spChanged$y_predicted, 0.25)
  spNotChangedFirstQ = quantile(spNotChanged$y_predicted, 0.25)
  
  spChangedThirdQ = quantile(spChanged$y_predicted, 0.75)
  spNotChangedThirdQ = quantile(spNotChanged$y_predicted, 0.75)
  
  if (spChangedFirstQ > spNotChangedFirstQ & spChangedThirdQ > spNotChangedThirdQ) {
    cutoff = abs(spChangedFirstQ + spNotChangedThirdQ)/2
    return(cutoff)
  } else {
    cutoff = abs(median(spChanged$y_predicted) + median(spNotChanged$y_predicted))/2
    return(cutoff)
  }
}


print(paste0("Analyze model for RQ4"))

getCutOff = function(data) {
  changedSP = data[data$y == TRUE,]
  unchangedSP = data[data$y == FALSE,]
  
  changedSpFirstQ = quantile(changedSP$y_predicted, 0.25)
  if (is.na(changedSpFirstQ)) {
    changedSpFirstQ = 0
  }
  unchangedSpFirstQ = quantile(unchangedSP$y_predicted, 0.25)
  
  changedSpThirdQ = quantile(changedSP$y_predicted, 0.75)
  if (is.na(changedSpThirdQ)) {
    changedSpThirdQ = 0
  }
  unchangedSpThirdQ = quantile(unchangedSP$y_predicted, 0.75)
  
  if (changedSpFirstQ > unchangedSpFirstQ & changedSpThirdQ > unchangedSpThirdQ) {
    cutoff = abs(changedSpFirstQ + unchangedSpThirdQ)/2
    return(cutoff)
  } else {
    cutoff = abs(median(changedSP$y_predicted) + median(unchangedSP$y_predicted))/2
    return(cutoff)
  }
}

run <- function(projectKey, i, testingData, trainingData, ldaTrainingData, ldaTestingData) {
  load(paste0("../data/rq4_model/", projectKey, "_lrm_", i, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_rf_", i, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_svm_", i, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_cart_", i, ".rda"))
  
  # -Calculate cutoff for each model based on the prediction of trainingData (for randomGuessing, we use 0.5)
  # -predict work items in testing data with the cutoff
  # -Calculate performance measure.
  
  #LRM AUC2
  trainingData$y_predicted = predict(lrm, trainingData, type="fitted")
  cutoff = getCutOff(trainingData)
  lrmPredicted = predict(lrm, testingData, type="fitted")
  lrmResult = calculatePRF(lrmPredicted, testingData$y, cutoff, "lrm")
  names(lrmResult) = paste0("lrm_", names(lrmResult))
  #RF AUC
  trainingData$y_predicted = data.frame(predict(rf, trainingData, type = "prob"))$`TRUE.`
  cutoff = getCutOff(trainingData)
  rfPredicted=data.frame(predict(rf, testingData, type = "prob"))$`TRUE.`
  rfResult = calculatePRF(rfPredicted, testingData$y, cutoff, "rf")
  names(rfResult) = paste0("rf_", names(rfResult))
  #SVM AUC
  tempSelfPredicted = predict(svm, newdata = trainingData, probability = TRUE, decision.values = TRUE)
  trainingData$y_predicted = data.frame(attr(tempSelfPredicted, "probabilities"))$'TRUE.'
  cutoff = getCutOff(trainingData)
  svmPredicted = predict(svm, newdata = testingData, probability = TRUE, decision.values = TRUE)
  svmPredicted = data.frame(attr(svmPredicted, "probabilities"))$'TRUE.'
  svmResult = calculatePRF(svmPredicted, testingData$y, cutoff, "svm")
  names(svmResult) = paste0("svm_", names(svmResult))
  #CART AUC
  trainingData$y_predicted = data.frame(predict(cart, trainingData, type = 'prob'))$`TRUE.`
  cutoff = getCutOff(trainingData)
  cartPredicted=data.frame(predict(cart, testingData, type = 'prob'))$`TRUE.`
  cartResult = calculatePRF(cartPredicted, testingData$y, cutoff, "rf")
  names(cartResult) = paste0("cart_", names(cartResult))
  
  # LOAD BASELINE
  load(paste0("../data/rq4_model/", projectKey, "_oner_", i, ".rda"))
  ruleDf = read.csv(paste0("../data/rq4_model/", projectKey, "_majority_", i, ".csv"), stringsAsFactors = TRUE)
  trainingData = trainingData %>% relocate(y, .after = last_col())
  
  #random guessing
  randomPredicted = randomGuess(length(testingData$y))
  randomResult = calculatePRF(randomPredicted, testingData$y, 0.5, "random")
  names(randomResult) = paste0("random_", names(randomResult))
  
  #oneR
  data.bin <- optbin(testingData)
  predictedOneR = as.integer(predict(model.OneR, testingData) == TRUE)
  onerResult = calculatePRF(predictedOneR, testingData$y, 0.5, "oner")
  names(onerResult) = paste0("oner_", names(onerResult))
  
  #MAJORITY VOTE AUC
  trainingData$y_predicted = predictMajorityVote(ruleDf, trainingData, ind_vars)
  cutoff = getCutOff(trainingData)
  majorityVotePredicted=predictMajorityVote(ruleDf, testingData, ind_vars)
  majorityVoteResult = calculatePRF(majorityVotePredicted, testingData$y, cutoff, "majority")
  names(majorityVoteResult) = paste0("majorityVote_", names(majorityVoteResult))
  
  ############ LDA #############
  load(paste0("../data/rq4_model/", projectKey, "_ldaRf_", i, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_ldaCart_", i, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_ldaNb_", i, ".rda"))
  #LDA RF
  topicCount = ldaRf$numTopic
  data = loadLdaDataFrame(i, projectKey, topicCount)
  ldaTrainingData = data$trainingData
  ldaTestingData = data$testingData
  ldaTrainingData$y_predicted = data.frame(predict(ldaRf, ldaTrainingData, type = "prob"))$`TRUE.`
  cutoff = getCutOff(ldaTrainingData)
  ldaRfPredicted=data.frame(predict(ldaRf, ldaTestingData, type = "prob"))$`TRUE.`
  ldaRfResult = calculatePRF(ldaRfPredicted, ldaTestingData$y, cutoff, "rf")
  names(ldaRfResult) = paste0("ldaRf_", names(ldaRfResult))
  #LDA Cart
  topicCount = ldaCart$numTopic
  data = loadLdaDataFrame(i, projectKey, topicCount)
  ldaTrainingData = data$trainingData
  ldaTestingData = data$testingData
  ldaTrainingData$y_predicted = data.frame(predict(ldaCart, ldaTrainingData, type = "prob"))$`TRUE.`
  cutoff = getCutOff(ldaTrainingData)
  ldaCartPredicted=data.frame(predict(ldaCart, ldaTestingData, type = "prob"))$`TRUE.`
  ldaCartResult = calculatePRF(ldaCartPredicted, ldaTestingData$y, cutoff, "rf")
  names(ldaCartResult) = paste0("ldaCart_", names(ldaCartResult))
  # LDA NB
  topicCount = ldaNb$numTopic
  data = loadLdaDataFrame(i, projectKey, topicCount)
  ldaTrainingData = data$trainingData
  ldaTestingData = data$testingData
  ldaTrainingData$y_predicted = data.frame(predict(ldaNb, ldaTrainingData, type = "prob"))$`TRUE.`
  cutoff = getCutOff(ldaTrainingData)
  ldaNbPredicted=data.frame(predict(ldaNb, ldaTestingData, type = "prob"))$`TRUE.`
  ldaNbResult = calculatePRF(ldaNbPredicted, ldaTestingData$y, cutoff, "rf")
  names(ldaNbResult) = paste0("ldaNb_", names(ldaNbResult))
  #Text NB
  nbDf = read.csv(paste0("../data/rq4_feature/textdimension/", projectKey, "_", i, ".csv"), stringsAsFactors = TRUE)
  nbDf$y = as.factor(nbDf$sp_change == "True" | nbDf$sp_change == TRUE)
  trainNbDf = nbDf[nbDf$GROUP != "S",]
  testNbDf = nbDf[nbDf$GROUP == "S",]
  
  trainNbDf$y_predicted = trainNbDf$text_bernoulli
  cutoff = getCutOff(trainNbDf)
  bernoulliResult = calculatePRF(testNbDf$text_bernoulli, testNbDf$y, cutoff, "bernoulliNb")
  names(bernoulliResult) = paste0("bernoulliNb_", names(bernoulliResult))
  
  trainNbDf$y_predicted = trainNbDf$text_gaussian
  cutoff = getCutOff(trainNbDf)
  gaussianResult = calculatePRF(testNbDf$text_gaussian, testNbDf$y, cutoff, "gaussianNb")
  names(gaussianResult) = paste0("gaussianNb_", names(gaussianResult))
  
  trainNbDf$y_predicted = trainNbDf$text_complement
  cutoff = getCutOff(trainNbDf)
  complementResult = calculatePRF(testNbDf$text_complement, testNbDf$y, cutoff, "complementNb")
  names(complementResult) = paste0("complementNb_", names(complementResult))
  
  trainNbDf$y_predicted = trainNbDf$text_multinomial
  cutoff = getCutOff(trainNbDf)
  multinomialResult = calculatePRF(testNbDf$text_multinomial, testNbDf$y, cutoff, "multinomialNb")
  names(multinomialResult) = paste0("multinomialNb_", names(multinomialResult))
  
  
  return(c(lrmResult, rfResult, svmResult, cartResult,
           randomResult, onerResult, majorityVoteResult,
           ldaRfResult, ldaCartResult, ldaNbResult,
           bernoulliResult, gaussianResult, complementResult, multinomialResult))
}


parentPath = "RQ4"
projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
for (projectKey in projectList) {
  resultDf = data.frame()
  # load datasets, use similar metrics as when we train the model (RQ4.r)
  load(paste0("../data/rq4_model/", projectKey, "_rf_", 0, ".rda"))
  ind_vars = rownames(data.frame(rf$importance))
  
  for (i in 0:99) {
    print(paste0(projectKey, "-", i))
    df = loadDataFrame(i, projectKey)
    
    finalColnames = c(ind_vars,'y')
    testingData = dplyr::select(df[df$type == 'test',], finalColnames)
    trainingData = dplyr::select(df[df$type == 'train',], finalColnames)
    
    tryCatch({ 
      results = run(projectKey, i, testingData, trainingData, ldaTrainingData, ldaTestingData)
      resultDf = rbind(resultDf, results)
    }, error=function(cond){
      print(cond)
      cat(paste0(projectKey, "-", i, "\n", cond), file="error.txt",sep="\n", append=TRUE)
    })
  }
  
  writeCsv(resultDf, paste0('../data/rq4_performance/',projectKey, '_LIST.csv'))
  outputDf <- data.frame()
  outputList = list(projectKey=projectKey)
  for (col in colnames(resultDf)) {
    outputList[col] = mean(na.omit(resultDf[,col]))
  }
  writeCsv(outputList, paste0('../data/rq4_performance/', projectKey, '.csv'))
}


