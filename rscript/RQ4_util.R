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
library(ggplot2)
########################################################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(Rnalytica)
library(rms)
library(randomForest)
library(e1071)
library(car)
library(ROCR)
library(caret)
library(OneR)
library(DEoptim)

source("util.R")
source("MajorityVote.R")

buildLrmCommand <- function(ind_vars, enable_df, spearmanResult, buildLrmData) {
  formula = "y ~ " 
  i = 1
  for (var in ind_vars) {
    rho2 = spearmanResult[i]
    
    if (nrow(unique(buildLrmData[var])) <= 2) {
      enable_df = FALSE
    }
    if (i > 1) {
      formula = paste0(formula, " + ")
    }
    if (enable_df & rho2 > 0.3) {
      # five degrees of freedom if ρ2 > 0.3
      formula = paste0(formula, "rcs(", var, ", ", 5, ")")
    } else if (enable_df & rho2 > 0.15) {
      # three degrees of freedom if ρ2 > 0.3
      formula = paste0(formula, "rcs(", var, ", ", 3, ")")
    } else {
      # one degree of freedom for linear relationship
      formula = paste0(formula, var)
    }
    i = i + 1
  }
  return(formula)
}


getLrm <- function(projectKey, ind_vars, trainingData) {
  # Calculate Spearman multiple ρ2
  spearmanResult <- spearman2(formula(paste("y ~ ",paste0(ind_vars, collapse=" + "))), data=trainingData, p=2)
  tryCatch({ 
    formula = buildLrmCommand(ind_vars, TRUE, spearmanResult, trainingData)
    lrm = lrm(as.formula(formula), data=trainingData, x=TRUE, y=TRUE, maxit = 1000)
  }, error=function(cond){
    print("LRM ERROR")
    formula = buildLrmCommand(ind_vars, FALSE, spearmanResult, trainingData)
    lrm = lrm(as.formula(formula), data=trainingData, x=TRUE, y=TRUE, maxit = 1000)
  })
  return(lrm)
}


calculatePRF <- function(predicted, result, cutoff, name) {
  predictedLabels = factor(predicted > cutoff, levels=c(FALSE,TRUE))
  if (sum(predictedLabels==FALSE) == length(predicted)) {
    acc = 0
  } else {
    xtab = table("a" = factor(predicted > cutoff, levels=c(FALSE,TRUE)), "b" = result)
    matrix = caret::confusionMatrix(xtab, positive = "TRUE")
    acc = matrix$byClass['Balanced Accuracy']
  }
  
  if (name != "oner") {
    auc = ModelMetrics::auc(actual = result, predicted = predicted)
  } else {
    auc = 0
  }
  
  library(caret)
  p <- caret::posPredValue(predictedLabels, result, positive=TRUE)
  r <- caret::sensitivity(predictedLabels, result, positive=TRUE)
  f1 <- 2 * (p * r) / (p + r)
  
  matrix = caret::confusionMatrix(data=predictedLabels, reference=result, dnn = c("Prediction", "Reference"), positive = 'TRUE')
  fp = matrix$table[2,1]
  tn = matrix$table[1,1]
  falseAlarmRate = fp/(fp+tn)
  distanceToHeaven = sqrt(
    (
      (1-r)^2 + (0-falseAlarmRate)^2
    ) / 2
  )
  
  return(list(auc = auc, 
              acc = acc,
              f = f1,
              precision = p, 
              recall = r,
              far = falseAlarmRate,
              d2h = distanceToHeaven))
}




randomGuess = function(size) {
  return(sample(0:1, size, replace=TRUE))
}

removeVars_new <- function(df) {
  binary_vars = c('has_TESTCASE', 'has_ATTACHMENT', 'has_CODE', 'has_STACK', 'STEP', 'OB', 'EB', 'TASK', 'LINK')
  numeric_vars = c('firststorypoints', 'flesch', 'fog', 
                   'lix', 'kinkaid', 'ari', 'colemanlieu', 'smog', 'text_multinomial', 
                   'text_gaussian', 'text_bernoulli', 'text_complement', 
                   'reporter_bugnum', 'reporter_recent_bugnum', 'reporter_stable_sp_rate', 
                   'assignee_bugnum', 'assignee_recent_bugnum', 'assignee_stable_sp_rate', 
                   #'estimator_bugnum', 'estimator_recent_bugnum', 'estimator_stable_sp_rate', 
                   'count_comment_before', 'count_summ_before', 'wordsDiff', 'sp_diff', 'sp_sprint',
                   'collab_betweenness', 'collab_closeness', 'collab_clustering_coeff', 'collab_eigen_vector', 
                   'collab_in_degree', 'collab_kcore', 'collab_out_degree', 'collab_total_degree')
  ind_vars = AutoSpearman(dataset = df, metrics = numeric_vars)
  ind_vars = c(ind_vars, binary_vars)
  ind_vars_copy = ind_vars
  for (var in ind_vars_copy) {
    count = nrow(unique(df[var]))
    if (count <= 1) {
      ind_vars = ind_vars[ind_vars != var]
    }
  }
  red <- redun(~., data=df[,ind_vars], nk=0) 
  reject_vars <- red$Out
  ind_vars <- ind_vars[!(ind_vars %in% reject_vars)]
  budgetted_DF = floor(min(nrow(df[df$y == TRUE,]), nrow(df[df$y == FALSE,]) )/15)
  return(ind_vars)
}


stringToBoolean <- function(df) {
  df$y = as.factor(df$sp_change == "True" | df$sp_change == TRUE)
  df$has_ATTACHMENT = as.factor(df$has_ATTACHMENT == "True" | df$has_ATTACHMENT == TRUE)
  df$has_CODE = as.factor(df$has_CODE == "True" | df$has_CODE == TRUE)
  df$has_STACK = as.factor(df$has_STACK == "True" | df$has_STACK == TRUE)
  df$has_TESTCASE = as.factor(df$has_TESTCASE == "True" | df$has_TESTCASE == TRUE)
  
  
  df$EB = as.factor(df$EB > 0)
  df$LINK = as.factor(df$LINK > 0)
  df$OB = as.factor(df$OB > 0)
  df$STEP = as.factor(df$STEP > 0)
  df$TASK = as.factor(df$TASK > 0)
  
  return(df)
}

loadDataFrame <- function(i, projectKey) {
  # load dataframe from data/rq4_feature/merge and data/reverted 
  revertedDf <- read.csv(paste0("../data/reverted/", projectKey, ".csv"), stringsAsFactors = TRUE)
  featureDf <- read.csv(paste0("../data/rq4_feature/merge/", projectKey, "_", i, ".csv"), stringsAsFactors = TRUE)
  df = merge(featureDf, dplyr::select(revertedDf, issue_key, isValid_sprint), by = "issue_key")
  df = stringToBoolean(df)
  df <- filterIssueValidSprint(df)
  # check no data leak, 
  trainAndValidateInTest = df[df[df$type != 'test',]$issue_key %in% df[df$type == 'test',]$issue_key,]
  assertthat::assert_that(nrow(trainAndValidateInTest) == 0)
  return(df)
}



loadLdaDataFrame <- function(i, projectKey, topicCount) {
  ldaDf = read.csv(paste0("../data/LdaResult/", topicCount, "/", projectKey, "_", i, ".csv"), stringsAsFactors = TRUE)
  ldaDf$y = as.factor(ldaDf$sp_change == "True")
  ind_vars = c()
  for (t in 0:(topicCount-1)) {
    ind_vars = c(ind_vars, paste0("t", t))
  }
  ind_vars <- ind_vars
  finalColNames <- c(ind_vars,'y')
  
  trainingData <- dplyr::select(ldaDf[ldaDf$type == 'train',], finalColNames)
  validatingData <- dplyr::select(ldaDf[ldaDf$type == 'valid',], finalColNames)
  testingData <- dplyr::select(ldaDf[ldaDf$type == 'test',], finalColNames)
  ldaData = c()
  ldaData$ind_vars = ind_vars
  ldaData$finalColNames = finalColNames
  ldaData$trainingData = trainingData
  ldaData$validatingData = validatingData
  ldaData$testingData = testingData
  return(ldaData)
}


setDataDist <- function(dataset) {
  dd <<- datadist(dataset)
  options(datadist="dd") 
}

filterCompleteCaseAndFixGaussianDataType <- function(dataset) {
  if("text_gaussian" %in% colnames(dataset)){
    dataset$text_gaussian = as.integer(round(dataset$text_gaussian))
  }
  return(dataset[complete.cases(dataset), ])
}

roundOptim = function(x){ x <- round(x)}

###############################
kLower = 1 # K of smote
kUpper = 20
pLower = 100 # percent of smote (e.g., 100 -> creating N samples), 0.1 for performanceEstimation package
pUpper = 300 # percent of smote (e.g., 300 -> creating N * 2 samples), 3 for performanceEstimation package
parVars = c("trainingData", "validatingData", "finalColNames", "ind_vars", 
            "projectKey", "getLrm", "buildLrmCommand", "optimDataPart", 
            "optimAucPart", "setDataDist", "i", "projectKey", "datadist", 
            "options", "predictMajorityVote", "str_split",
            "createMajorityVoteRule", "findMajorityFactorRule", 
            "findMajorityNumericalRule", "maximize_metric",
            "sens_constrain", "optbin", "OneR", "na.zero", "getFormula")
mNum = 30
mIter = 30
pType = 1
###############################

getFormula <- function(ind_vars) {
  return(formula(paste("y ~ ",paste0(ind_vars, collapse=" + "))))
}

normalDeOptim <- function(optimFunction) {
  outDEoptim = NA
  tryCatch({
    control = DEoptim.control(NP = mNum, itermax = mIter, parallelType = pType, parVar = parVars, packages = c("rms", "DMwR", "e1071", "cutpointr"))
    outDEoptim <- DEoptim(optimFunction, lower=c(kLower, pLower), upper=c(kUpper, pUpper), control = control, fnMap = roundOptim)
  }, error = function(e) {
    print(e)
  }
  )
  if (is.na(outDEoptim)) {
    # this only happened for few bootstrap iterations of TISTUD (AS). 
    # Error occurs here as there are few positive cases. We go for K = 5 and percent = 200
    output = c(5, 200)
    outDEoptim = c()
    outDEoptim$optim$bestmem = output
  }
  return(outDEoptim)
}

optimDataPart <- function(trainingData, kOptim, percentOptim, ind_vars) {
  newTrainingData = SMOTE(getFormula(ind_vars), trainingData, perc.over = percentOptim, k = kOptim, perc.under = 0)
  newTrainingData = rbind(trainingData, newTrainingData)
  newTrainingData = newTrainingData[sample(nrow(newTrainingData)),]
  
  numberCompleteCases = nrow(newTrainingData[complete.cases(newTrainingData), ])
  numberAll = nrow(newTrainingData)
  if (numberCompleteCases != numberAll) {
    cat(paste0(projectKey, "-", i, " na cases = ", numberAll - numberCompleteCases), file="error.txt",sep="\n", append=TRUE)
  }
  newTrainingData = newTrainingData[complete.cases(newTrainingData), ]
  return (newTrainingData)
}

optimAucPart <- function(validatingOptim, validatingData) {
  auc = ModelMetrics::auc(actual = validatingData$y, predicted = validatingOptim)
  return(auc*(-1))
}

optimRf <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  ntrees = params[3]
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  setDataDist(newTrainingData)
  rf <- randomForest(getFormula(ind_vars), data = newTrainingData, importance = FALSE, type = "classification", ntree=ntrees)
  return(optimAucPart(data.frame(predict(rf, validatingData, type = "prob"))$`TRUE.`, validatingData))
}

optimLrm <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  setDataDist(newTrainingData)
  lr = lrm(getFormula(ind_vars), data=newTrainingData, x=TRUE, y=TRUE, maxit = 1000)
  return(optimAucPart(predict(lr, validatingData, type="fitted"), validatingData))
}


optimSvm <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  costOptim = params[3] #1-6 convert to 0.001 - 100
  gammaOptim = params[4] #1-6 convert to 0.001 - 100
  mappingList = c(0.001, 0.01, 0.1, 1, 10, 100)
  
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  
  setDataDist(newTrainingData)
  svm <- svm(getFormula(ind_vars), newTrainingData, type = "C-classification", probability=TRUE, cost=mappingList[costOptim], gamma=mappingList[gammaOptim])
  tempPredicted = predict(svm, newdata = validatingData, probability = TRUE, decision.values = TRUE)
  return(optimAucPart(data.frame(attr(tempPredicted, "probabilities"))$'TRUE.', validatingData))
}



optimCart <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  complexity = params[3] # 1-11
  mappingList = c(0.001, 0.05, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  
  setDataDist(newTrainingData)
  cart <- rpart(getFormula(ind_vars), data = newTrainingData, control=rpart.control(cp=mappingList[complexity], minsplit=0))
  return(optimAucPart(data.frame(predict(cart, validatingData, type = 'prob'))$`TRUE.`, validatingData))
}




optimMajorityVote <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  
  ruleDf = createMajorityVoteRule(df = newTrainingData, ind_vars = ind_vars)
  predictedResult = predictMajorityVote(ruleDf = ruleDf, predictingDf = validatingData, ind_vars = ind_vars)
  return(optimAucPart(predictedResult, validatingData))
}


optimOneR <- function(params) {
  kOptim = params[1]
  percentOptim = params[2]
  newTrainingData = optimDataPart(trainingData, kOptim, percentOptim, ind_vars)
  
  data.bin <- optbin(newTrainingData)
  model.OneR <- OneR(formula = getFormula(ind_vars), data = newTrainingData, ties.method = c("first", "chisq")) #, verbose = TRUE)
  
  predictedResult = data.frame(na.zero(predict(model.OneR, validatingData, type="prob")))$'TRUE.'
  return(optimAucPart(predictedResult, validatingData))
}


