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
  absSize = abs(size)
  if (absSize >= 0.8) {
    return("L") 
  } else if (absSize >= 0.5) {
    return("M")
  } else if (absSize >= 0.2) {
    return("S")
  } else{
    return("N")
  }
}


performWilcoxTest <- function(a, b, shouldBeHigherThanBaseline) { # a hiegher than b ?
  if (shouldBeHigherThanBaseline) { # for AUC, BACC, Recall
    wilcoxResult = wilcoxsign_test(a ~ b, distribution="exact", alternative="greater")
    effectSize = getEffSize(attr(wilcoxResult@statistic, "teststatistic"), length(a))
    pvalue = coin::pvalue(wilcoxResult)[1]
  } else { # for FAR, D2H
    wilcoxResult = wilcoxsign_test(a ~ b, distribution="exact", alternative="less")
    effectSize = getEffSize(attr(wilcoxResult@statistic, "teststatistic"), length(b))
    pvalue = coin::pvalue(wilcoxResult)[1]
  }
  return(c(effectSize, pvalue))
}



getComparingText <- function(baseline, lrm, svm, rf, cart, shouldBeHigherThanBaseline) {
  # Convert effect size into 3 classes according to Rule of Thumb (getEffSizeText())
  # Convert p-value as well, i.e., *** p < 0.001, ** p < 0.01, * p < 0.05, ◦ p ≥ 0.05
  # Noted that, this will convert into (my)Latex command, i.e., \three = ***, \two = **, \one = *, \zero = ◦
  
  lrmWilcox = performWilcoxTest(lrm, baseline, shouldBeHigherThanBaseline)
  lrmEff = lrmWilcox[1]
  lrmP = lrmWilcox[2]
  lrmText = paste0(getEffSizeText(lrmEff), getPStar(lrmP))
  
  svmWilcox = performWilcoxTest(svm, baseline, shouldBeHigherThanBaseline)
  svmEff = svmWilcox[1]
  svmP = svmWilcox[2]
  svmText = paste0(getEffSizeText(svmEff), getPStar(svmP))
  
  rfWilcox = performWilcoxTest(rf, baseline, shouldBeHigherThanBaseline)
  rfEff = rfWilcox[1]
  rfP = rfWilcox[2]
  rfText = paste0(getEffSizeText(rfEff), getPStar(rfP))
  
  cartWilcox = performWilcoxTest(cart, baseline, shouldBeHigherThanBaseline)
  cartEff = cartWilcox[1]
  cartP = cartWilcox[2]
  cartText = paste0(getEffSizeText(cartEff), getPStar(cartP))
  return(c(lrmText, svmText, rfText, cartText))
}

getValue = function(projectKey, testName, random, oner, lrm, svm, rf, cart, ldaRf, ldaCart, ldaNb, majorityVote, textNb) {
  # use wilcoxsign test coz it return z statistic
  shouldBeHigherThanBaseline = !(testName == "D2H" || testName == "FAR")
  if ((testName == "D2H" || testName == "FAR") && shouldBeHigherThanBaseline == TRUE) {
    print("ASSERT WRONG")
  }
  
  compareWithOneR = getComparingText(oner, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmOneR = compareWithOneR[1]
  svmOneR = compareWithOneR[2]
  rfOneR = compareWithOneR[3]
  cartOneR = compareWithOneR[4]
  
  compareWithRandom = getComparingText(random, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmRandom = compareWithRandom[1]
  svmRandom = compareWithRandom[2]
  rfRandom = compareWithRandom[3]
  cartRandom = compareWithRandom[4]
  
  compareWithMajorityVote = getComparingText(majorityVote, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmMajorityVote = compareWithMajorityVote[1]
  svmMajorityVote = compareWithMajorityVote[2]
  rfMajorityVote = compareWithMajorityVote[3]
  cartMajorityVote = compareWithMajorityVote[4]
  
  compareWithLdaRf = getComparingText(ldaRf, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmLdaRf = compareWithLdaRf[1]
  svmLdaRf = compareWithLdaRf[2]
  rfLdaRf = compareWithLdaRf[3]
  cartLdaRf = compareWithLdaRf[4]
  
  compareWithLdaCart = getComparingText(ldaCart, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmLdaCart = compareWithLdaCart[1]
  svmLdaCart = compareWithLdaCart[2]
  rfLdaCart = compareWithLdaCart[3]
  cartLdaCart = compareWithLdaCart[4]
  
  compareWithLdaNb = getComparingText(ldaNb, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmLdaNb = compareWithLdaNb[1]
  svmLdaNb = compareWithLdaNb[2]
  rfLdaNb = compareWithLdaNb[3]
  cartLdaNb = compareWithLdaNb[4]
  
  compareWithTextNb = getComparingText(textNb, lrm, svm, rf, cart, shouldBeHigherThanBaseline)
  lrmTextNb = compareWithTextNb[1]
  svmTextNb = compareWithTextNb[2]
  rfTextNb = compareWithTextNb[3]
  cartTextNb = compareWithTextNb[4]
  
  
  df = data.frame(
    Proj = projectKey,
    Measure = testName,
    lrmLdaRf = lrmLdaRf,
    lrmMajorityVote = lrmMajorityVote,
    lrmOneR = lrmOneR,
    lrmTextNb = lrmTextNb,
    svmOneR = svmOneR,
    svmLdaRf = svmLdaRf,
    svmMajorityVote = svmMajorityVote,
    svmTextNb = svmTextNb,
    rfLdaRf = rfLdaRf,
    rfMajorityVote = rfMajorityVote,
    rfOneR = rfOneR,
    rfTextNb = rfTextNb,
    cartLdaRf = cartLdaRf,
    cartMajorityVote = cartMajorityVote,
    cartOneR = cartOneR,
    cartTextNb = cartTextNb
  )
  
  return(df)
}

projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
resultDf = data.frame()
for (projectKey in projectList) {
  performance <- read.csv(paste0("../data/RQ4_performance/", projectKey, "_LIST.csv"))
  #performance[is.na(performance)] <- 0
  auc = getValue(projectKey, "AUC", 
                 performance$random_auc, 
                 performance$oner_auc, 
                 performance$lrm_auc, 
                 performance$svm_auc, 
                 performance$rf_auc,
                 performance$cart_auc,
                 performance$ldaRf_auc,
                 performance$ldaCart_auc,
                 performance$ldaNb_auc,
                 performance$majorityVote_auc, 
                 performance$bernoulliNb_auc)
  resultDf = rbind(resultDf, auc)
  
  acc = getValue(projectKey, "BACC", 
                 performance$random_acc,
                 performance$oner_acc, 
                 performance$lrm_acc, 
                 performance$svm_acc, 
                 performance$rf_acc,
                 performance$cart_acc,
                 performance$ldaRf_acc,
                 performance$ldaCart_acc,
                 performance$ldaNb_acc,
                 performance$majorityVote_acc,
                 performance$bernoulliNb_acc)
  resultDf = rbind(resultDf, acc)
  
  d2h = getValue(projectKey, "D2H",
                 performance$random_d2h,
                 performance$oner_d2h, 
                 performance$lrm_d2h,
                 performance$svm_d2h, 
                 performance$rf_d2h,
                 performance$cart_d2h,
                 performance$ldaRf_d2h,
                 performance$ldaCart_d2h,
                 performance$ldaNb_d2h,
                 performance$majorityVote_d2h,
                 performance$bernoulliNb_d2h
  )
  resultDf = rbind(resultDf, d2h)
  
  
  far = getValue(projectKey, "FAR",
                 performance$random_far,
                 performance$oner_far, 
                 performance$lrm_far,
                 performance$svm_far, 
                 performance$rf_far,
                 performance$cart_far,
                 performance$ldaRf_far,
                 performance$ldaCart_far,
                 performance$ldaNb_far,
                 performance$majorityVote_far,
                 performance$bernoulliNb_far
  )
  resultDf = rbind(resultDf, far)
  
  
  recall = getValue(projectKey, "RECALL",
                    performance$random_recall,
                    performance$oner_recall, 
                    performance$lrm_recall,
                    performance$svm_recall, 
                    performance$rf_recall,
                    performance$cart_recall,
                    performance$ldaRf_recall,
                    performance$ldaCart_recall,
                    performance$ldaNb_recall,
                    performance$majorityVote_recall,
                    performance$bernoulliNb_recall
  )
  resultDf = rbind(resultDf, recall)
}


writeCsv(dataFrame = resultDf, "../data/RQ4_performance/stattest_paper.csv")
