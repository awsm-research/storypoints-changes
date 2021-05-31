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
########################################################################################

set.seed(1234)

#install.packages('devtools')
#install.packages("../DMwR_0.4.1.tar.gz", repos = NULL, type="source")
#devtools::install_github('software-analytics/Rnalytica')


# build classifiers
createModel <- function(projectKey, i, testingData, trainingData, ind_vars) {
  dir.create("../data/rq4_model/", showWarnings = FALSE)
  # One R
  data.bin <- optbin(trainingData)
  model.OneR <- OneR(formula = as.factor(y)~., data = trainingData, ties.method = c("first", "chisq")) #, verbose = TRUE)
  save(model.OneR, file = paste0("../data/rq4_model/", projectKey, "_oner_", i, ".rda"))
  #Linear Regression Model, for this one, we allocate additional degress of freedom to the explanatory variables ((MC-3) Constructing Classiﬁers)
  lrm = getLrm(projectKey = projectKey, ind_vars = ind_vars, trainingData = trainingData)
  save(lrm, file = paste0("../data/rq4_model/", projectKey, "_lrm_", i, ".rda"))
  #Random Forest
  rf <- randomForest(as.factor(y) ~ ., data = trainingData, importance = TRUE, type = "classification", ntree=10)
  save(rf, file = paste0("../data/rq4_model/", projectKey, "_rf_", i, ".rda"))
  #Support Vector Machine
  svm <- svm(as.factor(y) ~ ., trainingData, type = "C-classification", probability=TRUE)
  save(svm, file = paste0("../data/rq4_model/", projectKey, "_svm_", i, ".rda"))
}



run <- function(projectKey, from, to) {
  # determine what metrics will be used for this project
  df = loadDataFrame(1, projectKey)
  ind_vars = removeVars(df, projectKey)
  
  # For each bootstrap dataset
  for (i in from:to) {
    print(paste0(projectKey, "-", i))
    # load dataset of that bootstrap
    df = loadDataFrame(i, projectKey)
    finalColnames = c(ind_vars,'y')
    # 'S' means that, for this bootstrap dataset, the work item is for testing, otherwise it is for training
    testingData = dplyr::select(df[df$GROUP == 'S',], finalColnames)
    trainingData = dplyr::select(df[df$GROUP != 'S',], finalColnames)
    
    dd <<- datadist(trainingData)
    options(datadist="dd") 
    
    tryCatch({
      results = createModel(projectKey, i, testingData, trainingData, ind_vars)
    }, error=function(cond){
      cat(paste0(projectKey, "-", i, "\n", cond), file="error.txt",sep="\n", append=TRUE)
    })
  }
}


########################################################################################
### To run this, call "R RQ4.r <project key> <from iteration> <to iteration>
### For example, "R RQ4.r DM 0 99" to build and save models for 100 bootstrap datasets
########################################################################################

#projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "XD")
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
projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "XD")
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

