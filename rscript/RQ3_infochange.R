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
# RQ3: The contingency table of the work items with changed story points and the 
# work items with changed information
########################################################################################
# contingency tble
projectList = c("TISTUD", "DM", "MESOS", "MULE", "TIMOB", "USERGRID", "XD")
library(scales)
library(effsize)

for (project in projectList) {
  print(paste0("==== ", project, "===="))
  df = read.csv(paste0("../data/reverted/", project, ".csv"))
  df = filterIssueValidSprint(df)
  
  print("desc/summ changed")
  infochg = df$info_change == "True"
  infochg = factor(infochg, levels = c(TRUE, FALSE))
  spchg = df$sp_change == "True"
  spchg = factor(spchg, levels = c(TRUE, FALSE))
  tab = table(infochg, spchg)
  
  print(tab)
  fisher = fisher.test(tab, alternative = "greater")
  print(paste0("fisher odds ratio : ", format(fisher$estimate, digits=4), "  p :", format(fisher$p.value, digits =4)))
}

