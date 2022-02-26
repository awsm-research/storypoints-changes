# This file use to Examining the Inﬂuence of Metrics (Section 5.4.2 in the paper)
# First, we load classifiers that we built using RQ4.r.
# We extract the Wald Chi-square statistics of each metrics using anova test
# Then, we normalize the Wald χ2 of each metrics by the overall χ2 of our classifiers
# Finally, we use ScottKnott-ESD to rank the variables across 100 classifiers of the 100 bootstrap datasets

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

library(ScottKnottESD)
library(data.table)

projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
run <- function() {
  for (projectKey in projectList) {
    allAnovaDf = data.frame()
    for (i in 0:99) {
      print(paste0(projectKey, "-", i))
      
      # load prebuilt Logistic Regression classifiers (RQ4.r)
      load(paste0("../data/rq4_model/", projectKey, "_lrm_", i, ".rda"))
      # use anova to calculate Wald Chi-square statistics of each metrics
      anovaDf = data.frame(anova(lrm))
      setDT(anovaDf, keep.rownames = TRUE)[]
      totalChiSq = anovaDf[anovaDf$rn == "TOTAL",]$Chi.Square
      # normalize the Wald χ2 of each metrics by the overall χ2 of our classifiers
      anovaDf$chisqPercent = anovaDf$Chi.Square / totalChiSq * 100 
      anovaDf = anovaDf[anovaDf$rn != "TOTAL",]
      
      allAnovaDf = rbind(allAnovaDf, anovaDf)
    }
    skEsdTable = data.frame(matrix(0, ncol = 0, nrow = 100))
    for (var in unique(allAnovaDf$rn)) {
      skEsdTable[var] = allAnovaDf[allAnovaDf$rn == var,]$chisqPercent
    }
    # use ScottKnott-ESD to rank the variables across 100 classifiers of the 100 bootstrap datasets
    sk = sk_esd(skEsdTable)
    
    resultDf = data.frame(sk$groups)
    resultDf['project'] = shorternProjectKey(projectKey)
    setDT(resultDf, keep.rownames = TRUE)[]
    colnames(resultDf) = c("var", "rank", "project")
    writeCsv(resultDf, paste0('../data/rq4_metrics_influence/varrank_', projectKey, '.csv'))
  }
}

run()

