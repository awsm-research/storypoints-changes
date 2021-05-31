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
library(fpc)
library(Rnalytica)
library(rms)
library(randomForest)
library(e1071)
library(car)
library(ROCR)
library(caret)
library(OneR)


buildLrmCommand <- function(ind_vars, enable_df, spearmanResult, buildLrmData) {
  formula = "as.factor(y) ~ " 
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
    print(paste0("LRM ERROR", str(cond)))
    formula = buildLrmCommand(ind_vars, FALSE, spearmanResult, trainingData)
    lrm = lrm(as.formula(formula), data=trainingData, x=TRUE, y=TRUE, maxit = 1000)
  })
  return(lrm)
}

calculatePRF <- function(predicted, result, cutoff, name) {
  temp_predicted = predicted > cutoff
  if (sum(temp_predicted==FALSE) == length(predicted)) {
    acc = 0
  } else {
    xtab = table("a" = (predicted > cutoff), "b" = result)
    matrix = confusionMatrix(xtab, positive = "TRUE")
    acc = matrix$byClass['Balanced Accuracy']
  }
  
  auc = ModelMetrics::auc(actual = result, predicted = predicted)
  p = ModelMetrics::precision(actual = result, predicted = predicted, cutoff = cutoff)
  r = ModelMetrics::recall(actual = result, predicted = predicted, cutoff = cutoff)
  f1 = ModelMetrics::f1Score(actual = result, predicted = predicted, cutoff = cutoff)
  return(list(auc = auc, 
              acc = acc,
              f = f1,
              precision = p, 
              recall = r))
}


randomGuess = function(size) {
  return(sample(0:1, size, replace=TRUE))
}

removeVars <- function(df, projectKey) {
  
  # COMPLETENESS DIMENSION
  binary_vars = c('has_TESTCASE', # has-testcases
                  'has_ATTACHMENT', # has-attachment
                  'has_CODE', # has-code
                  'has_STACK', # has-stack
                  'STEP', # has-steps
                  'OB', # has-observed-behaviors
                  'EB', # has-expected-behaviors
                  'TASK', # has-bullets-task
                  'LINK') # has-link
  
  numeric_vars = c( # ACTIVITY DIMENSION
                   'count_summ_before', # activity-infochange-count
                   'wordsDiff', # activity-word-count
                   'count_comment_before', # activity-comment-count
                   'sp_sprint', # activity-sprint-sp
                   'firststorypoints', # activity-initial-sp
                   'sp_diff', # activity-sp-size-change
                   
                   # READABILITY DIMENSION
                   'flesch', # read-flesch
                   'fog', # read-fog
                   'lix', # read-lix
                   'kinkaid', # read-kincaid
                   'ari', # read-ari
                   'colemanlieu', #read-coleman-liau
                   'smog', # read-smog
                   
                   # TEXT DIMENSION
                   'text_multinomial', 
                   'text_gaussian', 
                   'text_bernoulli', 
                   'text_complement', 
                   
                   # EXPERIENCE DIMENSION
                   'reporter_bugnum', # reporter-workitem-num
                   'reporter_recent_bugnum', # reporter-recent-workitem-num
                   'reporter_stable_sp_rate', # reporter-stable-sp-rate
                   'assignee_bugnum', # assignee-workitem-num
                   'assignee_recent_bugnum', # assignee-recent-workitem-num
                   'assignee_stable_sp_rate', # assignee-stable-sp-rate
                   
                   # COLLABORATION DIMENSION
                   'collab_betweenness', 
                   'collab_closeness', 
                   'collab_clustering_coeff', 
                   'collab_eigen_vector', 
                   'collab_in_degree', 
                   'collab_kcore', 
                   'collab_out_degree', 
                   'collab_total_degree'
                   )
  # has-patch is not here because we did not detect any patch in any work items
  
                  
  survive_vars = AutoSpearman(dataset = df, metrics = numeric_vars)
  ind_vars = append(binary_vars, survive_vars)
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
  
  constantVars = colnames(Filter(function(x) length(unique(x))==1, df))
  ind_vars = ind_vars[!(ind_vars %in% constantVars)]
  
  # These metrics are excluded due to errors produced by LR models (lrm())
  if (projectKey == "XD") {
    ind_vars <- ind_vars[!(ind_vars %in% c("EB"))]
  } else if (projectKey == "DM") {
    ind_vars <- ind_vars[!(ind_vars %in% c("OB", "EB"))]
  } else if (projectKey == "TISTUD") {
    ind_vars <- ind_vars[!(ind_vars %in% c("sp_diff", "text_gaussian", "collab_kcore", "OB", "collab_closeness", "collab_betweenness", "collab_in_degree", "wordsDiff", "assignee_stable_sp_rate"))]
  } else if (projectKey == "MESOS") {
    ind_vars <- ind_vars[!(ind_vars %in% c("OB"))]
  } else if (projectKey == "MULE") {
    ind_vars <- ind_vars[!(ind_vars %in% c("EB", "OB"))]
  }
  return(ind_vars)
}

stringToBoolean <- function(df) {
  df$y = df$sp_change == "True"
  df$has_ATTACHMENT = df$has_ATTACHMENT == "True"
  df$has_CODE = df$has_CODE == "True"
  df$has_STACK = df$has_STACK == "True"
  df$has_TESTCASE = df$has_TESTCASE == "True"
  
  
  df$EB = df$EB > 0
  df$LINK = df$LINK > 0
  df$OB = df$OB > 0
  df$STEP = df$STEP > 0
  df$TASK = df$TASK > 0
  
  return(df)
}

loadDataFrame <- function(i, projectKey) {
  # load dataframe from data/rq4_feature/merge and data/reverted 
  revertedDf <- read.csv(paste0("../data/reverted/", projectKey, ".csv"))
  featureDf <- read.csv(paste0("../data/rq4_feature/merge/", projectKey, "_", i, ".csv"))
  df = merge(featureDf, dplyr::select(revertedDf, issue_key, isValid_sprint), by = "issue_key")
  df = stringToBoolean(df)
  df <- filterIssueValidSprint(df)
  return(df)
}