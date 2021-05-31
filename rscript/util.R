

writeCsv <- function(dataFrame, filePath) {
  print(paste0("write to ", filePath))
  write.csv(dataFrame, file = filePath ,row.names=FALSE, na="")
}

strSplit <- function(str, splitterChar) {
  return(unlist(strsplit(selectedCluster,split='-',fixed = TRUE)))
}

#Log Likelihood ratio test, if L.Ratio > 0, model2 is better than mode11
LRT <- function (model1, model2) { 
  L0 <- logLik(model1) #Get LogLikelihood of model1
  L1 <- logLik(model2) #Get LogLikelihood of model2
  LR <- as.vector(- 2 * (L0 - L1)) #Get -2LogLik ratio
  df <- attr(L1, "df") - attr(L0, "df") #Difference of D.F.
  p_value <- pchisq(LR, df, lower.tail = FALSE) # pchisq() compute area under the χ² curve for a specific number of degrees of freedom
  return(data.frame(L.Ratio = LR, Diff_D.F. = df, 
                    "p-value" =  p_value) )
} 

shorternProjectKey <- function(str) {
  r = ""
  if (str == "DM") {
    r = "DM"
  } else if (str == "MESOS") {
    r = "ME"
  }else if (str == "MULE") {
    r = "MU"
  }else if (str == "TDQ") {
    r = "TD"
  }else if (str == "TIMOB") {
    r = "TI"
  }else if (str == "TISTUD") {
    r = "AS"
  }else if (str == "XD") {
    r = "XD"
  }
  return(r)
}

filterIssueValidSprint <- function(dataFrame) {
  # Just to ensure that every work item have SP assigned not later than one hour after being assigned to a sprint 
  return(dataFrame[dataFrame['isValid_sprint'] == "True",])
}

getPStar = function(p) {
  # Noted that, this will convert into (my)Latex command, i.e., \three = *** (p<0.001), \two = ** (p<0.01), \one = * (p<0.05), \zero = ◦ (p>=0.05)
  if (p < 0.001) {
    return("\\three{}") 
  } else if (p < 0.01) {
    return("\\two{}") 
  } else if (p < 0.05) {
    return("\\one{}")
  } else{
    return("\\zero{}")
  }
}


na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
