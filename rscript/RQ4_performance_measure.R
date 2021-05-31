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

run <- function(projectKey, fileNum, testingData, trainingData) {
  # load models
  load(paste0("../data/rq4_model/", projectKey, "_lrm_", fileNum, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_rf_", fileNum, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_svm_", fileNum, ".rda"))
  load(paste0("../data/rq4_model/", projectKey, "_oner_", fileNum, ".rda"))
  trainingData = trainingData %>% relocate(y, .after = last_col())
  
  
  # -Calculate cutoff for each model based on the prediction of trainingData (for randomGuessing, we use 0.5)
  # -predict work items in testing data with the cutoff
  # -Calculate performance measure.
  
  #random guessing
  randomPredicted = randomGuess(length(testingData$y))
  randomResult = calculatePRF(randomPredicted, testingData$y, 0.5, "random")
  names(randomResult) = paste0("random_", names(randomResult))
  #oneR
  trainingData$y_predicted = data.frame(na.zero(predict(model.OneR, trainingData, type="prob")))$'TRUE.'
  cutoff = getCutOff(trainingData)
  predictedOneR = data.frame(na.zero(predict(model.OneR, testingData, type = "prob")))$'TRUE.'
  onerResult = calculatePRF(predictedOneR, testingData$y, cutoff, "oner")
  names(onerResult) = paste0("oner_", names(onerResult))
  #LRM
  trainingData$y_predicted = predict(lrm, trainingData, type="fitted")
  cutoff = getCutOff(trainingData)
  lrmPredicted = predict(lrm, testingData, type="fitted")
  lrmResult = calculatePRF(lrmPredicted, testingData$y, cutoff, "lrm")
  names(lrmResult) = paste0("lrm_", names(lrmResult))
  #RF
  trainingData$y_predicted = data.frame(predict(rf, trainingData, type = "prob"))$`TRUE.`
  cutoff = getCutOff(trainingData)
  rfPredicted=data.frame(predict(rf, testingData, type = "prob"))$`TRUE.`
  rfResult = calculatePRF(rfPredicted, testingData$y, cutoff, "rf")
  names(rfResult) = paste0("rf_", names(rfResult))
  #SVM
  tempSelfPredicted = predict(svm, newdata = trainingData, probability = TRUE, decision.values = TRUE)
  trainingData$y_predicted = data.frame(attr(tempSelfPredicted, "probabilities"))$'TRUE.'
  cutoff = getCutOff(trainingData)
  svmPredicted = predict(svm, newdata = testingData, probability = TRUE, decision.values = TRUE)
  svmPredicted = data.frame(attr(svmPredicted, "probabilities"))$'TRUE.'
  svmResult = calculatePRF(svmPredicted, testingData$y, cutoff, "svm")
  names(svmResult) = paste0("svm_", names(svmResult))
  return(c(lrmResult, rfResult, svmResult, randomResult, onerResult))
}


resultDf = data.frame()
projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "XD")
for (projectKey in projectList) {
  # load datasets, use similar metrics as when we train the model (RQ4.r)
  df = loadDataFrame(1, projectKey)
  ind_vars = removeVars(df, projectKey)
  
  for (i in 0:99) {
    print(paste0(projectKey, "-", i))
    # load each bootstrap datasets
    df = loadDataFrame(i, projectKey)
    
    finalColnames = c(ind_vars,'y')
    testingData = dplyr::select(df[df$GROUP == 'S',], finalColnames) # 'S' means the work items that will be used for testing
    trainingData = dplyr::select(df[df$GROUP != 'S',], finalColnames) # Hence, this line is for training datasets
    
    dd <<- datadist(testingData)
    options(datadist="dd") 
    
    tryCatch({ 
      results = run(projectKey, i, testingData, trainingData)
      resultDf = rbind(resultDf, results)
    }, error=function(cond){
      print(cond)
      cat(paste0(projectKey, "-", i, "\n", cond), file="error.txt",sep="\n", append=TRUE)
    })
  }
  
  writeCsv(resultDf, paste0('../data/RQ4_performance/', projectKey, '_LIST.csv'))
  outputDf <- data.frame()
  outputList = list(projectKey=projectKey)
  for (col in colnames(resultDf)) {
    outputList[col] = mean(na.omit(resultDf[,col]))
  }
  writeCsv(outputList, paste0('../data/RQ4_performance/', projectKey, '.csv'))
}

