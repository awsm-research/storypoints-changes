# This file builds OneR (baseline), Logistic Regression, Support Vector Machine, and Random Forest classifiers

# It first load data from data/rq4_feature combined with data/reverted 
# For each project, we determine the metrics that will be used using "removeVars()" function from RQ4_util.
# In removeVars(), we use AutoSpearman() and redun() to eliminate highly correlated and redundant metrics.
# Then, we separate testing and training dataset for each of 100 bootstrap datasets (feature_text_bootstrap).
# We build the classifiers using createModel() function and save them to data/rq4_model

# To run this, call "R RQ4.r <project key> <from iteration> <to iteration>
# For example, "R RQ4.r DM 0 99" to build and save models for 100 bootstrap datasets

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
source("MajorityVote.R")
########################################################################################

set.seed(1234)
home_data = "../data/"
modelFolderPath = paste0(home_data, "rq4_model/")
#install.packages('devtools')
#install.packages("../DMwR_0.4.1.tar.gz", repos = NULL, type="source")
#devtools::install_github('software-analytics/Rnalytica')


# build classifiers
createModel <- function(projectKey, i) {
  dir.create(file.path(modelFolderPath), showWarnings = FALSE)
  
  # RF ################## 
  pathRf = paste0(modelFolderPath, projectKey, "_rf_", i, ".rda")
  if (!file.exists(pathRf)) {
    print(paste0(projectKey, " - ", i, " RF"))
    rfControl = DEoptim.control(NP = mNum, itermax = mIter, parallelType = pType, parVar = parVars, packages = c("randomForest", "DMwR"))
    rfOutDEoptim <- DEoptim(optimRf, lower=c(kLower, pLower, 50), upper=c(kUpper, pUpper, 300), control = rfControl, fnMap = roundOptim)
    
    smoteK = rfOutDEoptim$optim$bestmem[1]
    smotePercent = rfOutDEoptim$optim$bestmem[2]
    nTree = round(rfOutDEoptim$optim$bestmem[3], digits = 0)
    rfSmotedData = SMOTE(getFormula(ind_vars), rbind(trainingData, validatingData), perc.over = smotePercent, k = smoteK, perc.under = 0)
    rfTrainingData = shuffleDf(rbind(trainingData, validatingData, rfSmotedData))
    rfTrainingData = filterCompleteCaseAndFixGaussianDataType(rfTrainingData)
    setDataDist(rfTrainingData) 
    
    rf <- randomForest(getFormula(ind_vars), data = rfTrainingData, importance = TRUE, type = "classification", ntree=nTree)
    save(rf, file = pathRf)
  }
  
  #LRM ################## 
  pathLrm = paste0(modelFolderPath, projectKey, "_lrm_", i, ".rda")
  if (!file.exists(pathLrm)) {
    print(paste0(projectKey, " - ", i, " LRM"))
    lrOutDEoptim = normalDeOptim(optimLrm)
    smoteK = lrOutDEoptim$optim$bestmem[1]
    smotePercent = lrOutDEoptim$optim$bestmem[2]
    
    lrSmotedData = SMOTE(y ~ ., rbind(trainingData, validatingData), perc.over = smotePercent, k = smoteK, perc.under = 0)
    lrTrainingData = shuffleDf(rbind(trainingData, validatingData, lrSmotedData))
    lrTrainingData = filterCompleteCaseAndFixGaussianDataType(lrTrainingData)
    setDataDist(lrTrainingData)
    
    lrm = getLrm(projectKey = projectKey, ind_vars = ind_vars, trainingData = lrTrainingData)
    save(lrm, file = pathLrm)
  }
  
  #SVM ################## 
  pathSvm = paste0(modelFolderPath, projectKey, "_svm_", i, ".rda")
  if (!file.exists(pathSvm)) {
    print(paste0(projectKey, " - ", i, " SVM"))
    svmControl = DEoptim.control(NP = 40, itermax = mIter, parallelType = pType, parVar = parVars, packages = c("e1071", "DMwR"))
    svmOutDEoptim <- DEoptim(optimSvm, lower=c(kLower, pLower, 1, 1), upper=c(kUpper, pUpper, 6, 6), control = svmControl, fnMap = roundOptim)
    
    smoteK = svmOutDEoptim$optim$bestmem[1]
    smotePercent = svmOutDEoptim$optim$bestmem[2]
    costOptim = svmOutDEoptim$optim$bestmem[3]
    gammaOptim = svmOutDEoptim$optim$bestmem[4]
    mappingList = c(0.001, 0.01, 0.1, 1, 10, 100)
    
    svmSmotedData = SMOTE(y ~ ., rbind(trainingData, validatingData), 
                          perc.over = smotePercent, k = smoteK, perc.under = 0, 
                          cost=mappingList[costOptim], gamma=mappingList[gammaOptim])
    svmTrainingData = shuffleDf(rbind(trainingData, validatingData, svmSmotedData))
    svmTrainingData = filterCompleteCaseAndFixGaussianDataType(svmTrainingData)
    setDataDist(svmTrainingData)
    
    svm <- svm(y ~ ., svmTrainingData, type = "C-classification", probability=TRUE)
    save(svm, file = pathSvm)
  }
  
  #CART ######### http://www.sthda.com/english/articles/35-statistical-machine-learning-essentials/141-cart-model-decision-tree-essentials/
  pathCart = paste0(modelFolderPath, projectKey, "_cart_", i, ".rda")
  if (!file.exists(pathCart)) {
    print(paste0(projectKey, " - ", i, " CART"))
    cartControl = DEoptim.control(NP = mNum, itermax = mIter, parallelType = pType, parVar = parVars, packages = c("rpart", "DMwR"))
    cartOutDEoptim <- DEoptim(optimCart, lower=c(kLower, pLower, 1), upper=c(kUpper, pUpper, 11), control = cartControl, fnMap = roundOptim)
    
    smoteK = cartOutDEoptim$optim$bestmem[1]
    smotePercent = cartOutDEoptim$optim$bestmem[2]
    cartComplexity = cartOutDEoptim$optim$bestmem[3]
    mappingList = c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    cartSmotedData = SMOTE(getFormula(ind_vars), rbind(trainingData, validatingData), perc.over = smotePercent, k = smoteK, perc.under = 0)
    cartTrainingData = shuffleDf(rbind(trainingData, validatingData, cartSmotedData))
    cartTrainingData = filterCompleteCaseAndFixGaussianDataType(cartTrainingData)
    setDataDist(cartTrainingData) 
    
    cart <- rpart(getFormula(ind_vars), data = cartTrainingData, control=rpart.control(cp=mappingList[cartComplexity], minsplit=0))
    #cart <- prune.rpart(cart, mappingList[cartComplexity])
    # pruning is not needed as we optimize the CP already -> https://stackoverflow.com/questions/37721047/selecting-cp-value-for-decision-tree-pruning-using-rpart
    par(xpd = NA) # otherwise on some devices the text is clipped
    save(cart, file = pathCart)
  }
  
  # One R #################################### 
  pathOneR = paste0(modelFolderPath, projectKey, "_oner_", i, ".rda")
  if (!file.exists(pathOneR)) {
    print(paste0(projectKey, " - ", i, " OneR"))
    oneROutDEoptim = normalDeOptim(optimOneR)
    
    smoteK = oneROutDEoptim$optim$bestmem[1]
    smotePercent = oneROutDEoptim$optim$bestmem[2]
    
    oneRSmotedData = SMOTE(y ~ ., rbind(trainingData, validatingData), perc.over = smotePercent, k = smoteK, perc.under = 0)
    oneRTrainingData = shuffleDf(rbind(trainingData, validatingData, oneRSmotedData))
    oneRTrainingData = filterCompleteCaseAndFixGaussianDataType(oneRTrainingData)
    data.bin <- optbin(trainingData)
    
    model.OneR <- OneR(formula = y~., data = oneRTrainingData, ties.method = c("first", "chisq")) #, verbose = TRUE)
    save(model.OneR, file = pathOneR)
  }
  
  # majorityVote ###################
  pathMajority = paste0(modelFolderPath, projectKey, "_majority_", i, ".csv")
  if (!file.exists(pathMajority)) {
    print(paste0(projectKey, " - ", i, " Majority"))
    majorityVoteOutDEoptim = normalDeOptim(optimMajorityVote)
    
    smoteK = majorityVoteOutDEoptim$optim$bestmem[1]
    smotePercent = majorityVoteOutDEoptim$optim$bestmem[2]
    
    majorityVoteSmotedData = SMOTE(y ~ ., rbind(trainingData, validatingData), perc.over = smotePercent, k = smoteK, perc.under = 0)
    majorityVoteTrainingData = shuffleDf(rbind(trainingData, validatingData, majorityVoteSmotedData))
    majorityVoteTrainingData = filterCompleteCaseAndFixGaussianDataType(majorityVoteTrainingData)
    
    majorityVoteRuleDf = createMajorityVoteRule(majorityVoteTrainingData, ind_vars)
    writeCsv(majorityVoteRuleDf, pathMajority)
  }
}



run <- function(projectKey, from, to) {
  # determine what metrics will be used for this project
  survivedColumns = c()
  for (i in 0:99) {
    trainTempDf = loadDataFrame(i, projectKey)
    trainTempDf = trainTempDf[trainTempDf$type == 'train',]
    
    cols = removeVars_new(trainTempDf)
    survivedColumns = c(survivedColumns, cols)
  }
  tempTable = table(survivedColumns)
  highestCount = max(tempTable)
  ind_vars <<- names(tempTable[tempTable == highestCount])

  # For each bootstrap dataset
  for (i in from:to) {
    projectKey <<- projectKey
    i <<- i
    print(paste0(projectKey, "-", i))
    # load dataset of that bootstrap
    df = loadDataFrame(i, projectKey)
    
    # R has a bug somehow these variables cannot passed through parVars
    # as this script is running sequentially, we use global vars
    finalColNames <<- c(ind_vars,'y')
    trainingData <<- dplyr::select(df[df$type == 'train',], finalColNames)
    validatingData <<- dplyr::select(df[df$type == 'valid',], finalColNames)
    testingData <<- dplyr::select(df[df$type == 'test',], finalColNames)
    ###############################################################
    # https://stats.stackexchange.com/questions/423230/problems-with-smote-optimizing-function
    ###############################################################
    tryCatch({
      createModel(projectKey, i)
    }, error=function(cond){
      print(cond)
      cat(paste0(projectKey, "-", i, "\n", cond), file="error.txt",sep="\n", append=TRUE)
    })
  }
}


########################################################################################
### To run this, call "R RQ4.r <project key> <from iteration> <to iteration>
### For example, "R RQ4.r DM 0 99" to build and save models for 100 bootstrap datasets
########################################################################################

#projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
#args<-commandArgs(TRUE)
#projectKey = args[1]
#from = args[2]
#to = args[3]
#print(projectKey)
#print(from)
#print(to)
#run(projectKey, from, to)



########################################################################################
### run straight for 0-99 bootstraps dataset
########################################################################################
projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
for (projectKey in projectList) {
  from = 0
  to = 99
  run(projectKey, from, to)
}

########################################################################################
### for printing avilable variable in models
########################################################################################

#for (projectKey in projectList) {
#  df = loadDataFrame(1, projectKey)
#  ind_vars = removeVars(df, projectKey)
#  print(projectKey)
#  print(ind_vars)
#}

