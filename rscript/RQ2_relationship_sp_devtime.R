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

########################################################################################
# Table 3: Statistics summary of our models in RQ2. The larger the χ2 value is, 
# the better the ﬁt of our models (and the SP values). ”LRχ 2” indicates LR χ2 
# of the model and ”Waldχ 2” indicates the normalized Wald χ2 of SP.
########################################################################################
library(dplyr)
library(plyr)
library(car)
library(multilevel) 
library(lmtest)
library(nlme)
getLRatio = function(fullModel) {
  lrResult = lmtest::lrtest(fullModel) # LR test compare our (full) model with a null model
  lRatioChisq = lrResult$Chisq[2] # calculate Log-likelihood ratio
  lRatioPvalue = lrResult$`Pr(>Chisq)`[2] # extract p-value
  
  aResult = Anova(fullModel, type = 3)  # perform ANOVA type 3
  spChisq = aResult$Chisq[2] # extract Chi-square value of SP
  normalizeSpChisq = spChisq / lRatioChisq # normalized the SP chisq by the LRatio chisq
  spPvalue = aResult$`Pr(>Chisq)`[2]
  
  lRatioPLevel = getPStar(lRatioPvalue)
  spPLevel = getPStar(spPvalue)
  
  return(list(round(lRatioChisq, 2), lRatioPLevel, round(normalizeSpChisq, 2), spPLevel))
}


getMMRE = function(model, data) {
  predicted = predict(model, newdata=data, level=0)
  mre = abs(data$time_in_status - predicted) / data$time_in_status
  mmre = mean(mre)
  return(mmre)
}

getSa = function(model, data) {
  predicted = predict(model, newdata=data, level=0)
  
  max = max(data$time_in_status)
  min = min(data$time_in_status)
  set.seed(0)
  random <- sample(min:max, nrow(data))
  maeRandom = sum(abs(data$time_in_status - random))/length(random)
  maePredicted = sum(abs(data$time_in_status - predicted))/length(predicted)
  
  sa = (1 - maePredicted/maeRandom) * 100
  return(sa)
}

run = function(project) {
  timeDf = read.csv(paste0("../data/rq2_inprogresstime/", project, ".csv"))
  revertedDf = read.csv(paste0("../data/reverted/", project, ".csv"))
  mergeDf = merge(revertedDf, dplyr::select(timeDf, issue_key, time_in_status), by = "issue_key")
  mergeDf = filterIssueValidSprint(mergeDf)
  mergeDf$priority = as.factor(mergeDf$priority)
  
  spChgDf <<- mergeDf[mergeDf$sp_change == "True",]
  spChgDf <<- na.exclude(object = spChgDf)
  noSpChgDf <<- mergeDf[mergeDf$sp_change == "False",]
  
  # build Linear Mixed-Effects Models between time that a work item's status is set to 'in progress' and Story Points
  
  # for unchangedSP
  unchangedSpModel = lme(time_in_status~storypoints, random=list(T_sp_year=~1, issuetype=~1, priority=~1), data=noSpChgDf)
  temp = getLRatio(unchangedSpModel)
  unchangedSp_L = temp[1] # Log Likelihood Ratio of the model
  unchangedSp_L_P = temp[2] # P-value of Log Likelihood Ratio
  unchangedSp_Sp = temp[3] # Normalize Wald-Chisqareof Story Points
  unchangedSp_P = temp[4] # P-value of Normalize Wald-Chisqareof Story Points
  
  unchangedSp_mmre = getMMRE(unchangedSpModel, noSpChgDf)
  unchangedSp_sa = round(getSa(unchangedSpModel, noSpChgDf), digits=2)
  
  # for changedSP sprint
  changedSpSprintModel = lme(time_in_status~sp_sprint, random=list(T_sp_year=~1, issuetype=~1, priority=~1), data=spChgDf)
  temp = getLRatio(changedSpSprintModel)
  changedSpSprint_L = temp[1]
  changedSpSprint_L_P = temp[2]
  changedSpSprint_Sp = temp[3]
  changedSpSprint_P = temp[4]
  
  changedSpSprint_mmre = getMMRE(changedSpSprintModel, spChgDf)
  changedSpSprint_sa = round(getSa(changedSpSprintModel, spChgDf), digits = 2)
  
  # for changedSP last
  changedSpLastModel = lme(time_in_status~storypoints, random=list(T_sp_year=~1, issuetype=~1, priority=~1), data=spChgDf)
  temp = getLRatio(changedSpLastModel)
  changedSpLast_L = temp[1]
  changedSpLast_L_P = temp[2]
  changedSpLast_Sp = temp[3]
  changedSpLast_Sp_P = temp[4]
  
  changedSpLast_mmre = getMMRE(changedSpLastModel, spChgDf)
  changedSpLast_sa = round(getSa(changedSpLastModel, spChgDf), digits = 2)
  
  
  ##############        unchgSP LRχ2 | Pvalue             chgSP-sprint LRχ2 | Pvalue                chgSP-last LR χ2 | Pvalue                         unchgSP χ2 | Pvalue                chgSP-sprint Wald χ2 | Pvalue                chgSP-last Wald χ2 | Pvalue 
  cat(sprintf(paste0(unchangedSp_L, " \\", unchangedSp_L_P, " & ", changedSpSprint_L, " \\", changedSpSprint_L_P, " & ", changedSpLast_L, " \\", changedSpLast_L_P, " & ", unchangedSp_Sp,  " \\", unchangedSp_P, " & ", changedSpSprint_Sp,  " \\", changedSpSprint_P, " & ", changedSpLast_Sp,  " \\", changedSpLast_Sp_P)))
  cat("\n")
  cat(sprintf(paste0(unchangedSp_sa, " & ", changedSpSprint_sa, " & ", changedSpLast_sa)))
}

run("TISTUD")
run("DM")
run("MESOS")
run("MULE")
run("TIMOB")
run("USERGRID")
run("XD")
