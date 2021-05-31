# Developer experience is the second most inﬂuential metric for our datasets.
# This file observe the effect of the change of reporter-stable-sp-rate to the odds (likelihood) that the SP will be changed.
# To do so, we use summary() function to observe 'Effect' of the metrics.
# We do this for each of 100 classifiers and find average odds.

# Value X in 'Effect' column means that, if value of a metric increasing from column 'Low' (i.e., 1st Q.) to column 'High' (i.e., 3rd Q.),
# the value of 'y' (i.e., the likelihood that SP will be changed) will be increased by X.


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

library(data.table)
library(ScottKnottESD)


projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "XD")

allOddsDf = data.frame()
for (projectKey in projectList) {
  oddsDf = data.frame()
  
  for (i in 0:99) {
    # load Logistic Regression classifiers built from RQ4.r
    print(paste0(projectKey, "-", i))
    load(paste0("../data/rq4_model/", projectKey, "_lrm_", i, ".rda"))
    summaryDf = data.frame(summary(lrm))
    setDT(summaryDf, keep.rownames = TRUE)[]
    
    # extract "Effect" value of each metrics (i.e., Odds / Likelihood)
    tempVarName = ""
    for (j in 1:nrow(summaryDf)) {
      if (j %% 2 == 1) {
        tempVarName = summaryDf$rn[j]
      } else {
        odds = summaryDf[j, "Effect"] - 1
        tempDf = data.frame(var = tempVarName, odds = odds)
        oddsDf = rbind(oddsDf, tempDf)
      }
    }
  }
  
  # calculate mean odds 
  for (v in unique(oddsDf$var)) {
    oddsMean = mean(oddsDf[oddsDf$var == v,]$Effect)
    oddsMean = format(oddsMean, digits = 4)
    toAdd = data.frame(project = projectKey, var = v, oddsMean = oddsMean)
    allOddsDf = rbind(allOddsDf, toAdd)
  }
}
# We observe only reporter_stable_sp_rate
allOddsDf = allOddsDf[allOddsDf$var == "reporter_stable_sp_rate",]
writeCsv(allOddsDf, paste0('../data/rq4_metrics_influence/odds_reporter_stable_sp_rate.csv'))

