# We evaluate the prediction performance of our classifiers with baseline approaches.
# To do so, we use Wilcoxon signed-rank test to compare the performance of our classifiers with OneR.
# We also measure the eﬀect size. A larger effect size represents a better performance of our models over OneR.


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


library(coin)

getEffSize = function(z, n) {
  return(z/sqrt(n))
}

getEffSizeText = function(size) {
  # According to effect size's rule of thumb, effect size 0.8 = Large, 0.5-0.8 = Medium, 0.2-0.5 = Small, else Negligible.
  if (size >= 0.8) {
    return("L") 
  } else if (size >= 0.5) {
    return("M")
  } else if (size >= 0.2) {
    return("S")
  } else{
    return("N")
  }
}

getValue = function(projectName, testName, random, oner, lrm, svm, rf) {
  # compare Lotistic Regression with oneR
  oner_lrmWilcox = wilcoxsign_test(lrm ~ oner, distribution="exact", alternative="greater")
  oner_lrmEff = getEffSize(attr(oner_lrmWilcox@statistic, "teststatistic"), length(lrm))
  oner_lrmP = coin::pvalue(oner_lrmWilcox)[1]
  
  # compare Support Vector Machine with oneR
  oner_svmWilcox = wilcoxsign_test(svm  ~ oner, distribution="exact", alternative="greater")
  oner_svmEff = getEffSize(attr(oner_svmWilcox@statistic, "teststatistic"), length(svm))
  oner_svmP = coin::pvalue(oner_svmWilcox)[1]
  
  # compare Random Forest with oneR
  oner_rfWilcox = wilcoxsign_test(rf ~ oner, distribution="exact", alternative="greater")
  oner_rfEff = getEffSize(attr(oner_rfWilcox@statistic, "teststatistic"), length(rf))
  oner_rfP = coin::pvalue(oner_rfWilcox)[1]
  
  # compare Lotistic Regression with random guessing
  random_lrmWilcox = wilcoxsign_test(lrm ~ random, distribution="exact", alternative="greater")
  random_lrmEff = getEffSize(attr(random_lrmWilcox@statistic, "teststatistic"), length(lrm))
  random_lrmP = coin::pvalue(random_lrmWilcox)[1]
  
  # compare Support Vector Machine with random guessing
  random_svmWilcox = wilcoxsign_test(svm  ~ random, distribution="exact", alternative="greater")
  random_svmEff = getEffSize(attr(random_svmWilcox@statistic, "teststatistic"), length(svm))
  random_svmP = coin::pvalue(random_svmWilcox)[1]
  
  # compare Random Forest with random guessing
  random_rfWilcox = wilcoxsign_test(rf ~ random, distribution="exact", alternative="greater")
  random_rfEff = getEffSize(attr(random_rfWilcox@statistic, "teststatistic"), length(rf))
  random_rfP = coin::pvalue(random_rfWilcox)[1]
  
  # Convert effect size into 3 classes according to Rule of Thumb (getEffSizeText())
  # Convert p-value as well, i.e., *** p < 0.001, ** p < 0.01, * p < 0.05, ◦ p ≥ 0.05
  # Noted that, this will convert into (my)Latex command, i.e., \three = ***, \two = **, \one = *, \zero = ◦
  lrmOneR = paste0(getEffSizeText(oner_lrmEff), getPStar(oner_lrmP))
  rfOneR = paste0(getEffSizeText(oner_rfEff), getPStar(oner_rfP))
  svmOneR = paste0(getEffSizeText(oner_svmEff), getPStar(oner_svmP))
  
  lrmRandom = paste0(getEffSizeText(random_lrmEff), getPStar(random_lrmP))
  rfRandom = paste0(getEffSizeText(random_rfEff), getPStar(random_rfP))
  svmRandom = paste0(getEffSizeText(random_svmEff), getPStar(random_svmP))
  
  df = data.frame(
    Proj = projectName,
    Measure = testName,
    lrmOneR = lrmOneR,
    lrmRandom = lrmRandom,
    rfOneR = rfOneR,
    rfRandom = rfRandom,
    svmOneR = svmOneR,
    svmRandom = svmRandom
    )
  
  return(df)
}

projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "XD")
resultDf = data.frame()
for (projectKey in projectList) {
  print(paste0("calculating performance for ", projectKey)) 
  performance <- read.csv(paste0("../data/RQ4_performance/", projectKey, "_LIST.csv"))
  performance[is.na(performance)] <- 0
  
  # compare prediction performance for each performance metrics
  auc = getValue(projectKey, "AUC", performance$random_auc, performance$oner_auc, performance$lrm_auc, performance$svm_auc, performance$rf_auc)
  resultDf = rbind(resultDf, auc)
  
  acc = getValue(projectKey, "BACC", performance$random_acc, performance$oner_acc, performance$lrm_acc, performance$svm_acc, performance$rf_acc)
  resultDf = rbind(resultDf, acc)
  
  f = getValue(projectKey, "F1-score", performance$random_f, performance$oner_f, performance$lrm_f, performance$svm_f, performance$rf_f)
  resultDf = rbind(resultDf, f)
  
  precision= getValue(projectKey, "precision", performance$random_precision, performance$oner_precision, performance$lrm_precision, performance$svm_precision, performance$rf_precision)
  resultDf = rbind(resultDf, precision)
  
  recall = getValue(projectKey, "recall", performance$random_recall, performance$oner_recall, performance$lrm_recall, performance$svm_recall, performance$rf_recall)
  resultDf = rbind(resultDf, recall)
}

writeCsv(dataFrame = resultDf, "../data/RQ4_performance/stattest_paper.csv")

