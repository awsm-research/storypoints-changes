library(cutpointr)
library(stringr)
findMajorityFactorRule <- function(df, varname) {
  # ถ้าเป็น true หมด จะได้ acc เท่าไหร่
  # ถ้าเป็น false หมด จะได้ acc เท่าไหร่
  trueVarPredicted = df[df[varname] == TRUE,]$y
  trueVarTruePositivePercent = length(trueVarPredicted[trueVarPredicted==TRUE])/length(trueVarPredicted) # sensitivity
  if (is.na(trueVarTruePositivePercent)) {
    trueVarTruePositivePercent = 0
  }
  falseVarPredicted = df[df[varname] == FALSE,]$y
  falseVarTruePositivePercent = length(falseVarPredicted[falseVarPredicted==TRUE])/length(falseVarPredicted) # sensitivity
  if (is.na(falseVarTruePositivePercent)) {
    falseVarTruePositivePercent = 0
  }
  
  if (trueVarTruePositivePercent >= falseVarTruePositivePercent) {
    trueValueClass = TRUE
    falseValueClass = FALSE
  } else {
    trueValueClass = FALSE
    falseValueClass = TRUE
  }
  return(paste0(TRUE,":",trueValueClass,",", FALSE, ":", falseValueClass))
}

findMajorityNumericalRule <- function(df, varname) {
  df$temp = df[[varname]]
  cp <- cutpointr(df, 'temp', 'y', method = maximize_metric, metric = sens_constrain, pos_class=TRUE, neg_class=FALSE)
  #e.g., >=4 is TRUE
  if (cp$direction == "<=") {
    direction = "lowerOrEqual"
  } else {
    direction = "greaterOrEqual"
  }
    
  return(paste0(direction, ',', cp$optimal_cutpoint))
}


createMajorityVoteRule <- function(df, ind_vars) {
  columnDf = data.frame(sapply(df[,ind_vars], class))
  columnDf <- cbind(newColName = rownames(columnDf), columnDf)
  rownames(columnDf) <- NULL
  colnames(columnDf) <- c('varname', 'cls')
  
  columnDf$rule = apply(columnDf, 1, function(x) { 
    if (x['cls'][[1]] == 'numeric' || x['cls'][[1]] == 'integer') { return(findMajorityNumericalRule(df, x['varname'][[1]])) }
    else {findMajorityFactorRule(df, x['varname'][[1]])}
  })
  return(columnDf)
}

predictMajorityVote <- function(ruleDf, predictingDf, ind_vars) {
  #ruleDf = read.csv(rulePath, stringsAsFactors = TRUE)
  resultDf = apply(predictingDf, 1, function(x) { 
    #print("__________________________________________________")
    #print(x['issue_key'])
    votes = c()
  
    for (var in ind_vars) {
      ruleRow = ruleDf[ruleDf$varname == var,]
      if (ruleRow$cls == 'numeric' || ruleRow$cls == 'integer') {
        rowVarValue = as.numeric(trimws(x[var]))
        ruleDirection = str_split(ruleRow$rule, ",")[[1]][1]
        ruleNumValue = as.numeric(trimws(str_split(ruleRow$rule, ",")[[1]][2]))
        if (ruleDirection == "lowerOrEqual") {
          vote = rowVarValue <= ruleNumValue
        } else {
          vote = rowVarValue >= ruleNumValue
        }
        votes = c(votes, vote)
        #print(paste0("NUM, x['", var, "'] ", ruleRow$rule, " (", rowVarValue, "=", vote, ")" ))
      } else {
        rowVarValue = x[var][[1]] == "TRUE"
        stopifnot(rowVarValue == FALSE || rowVarValue == TRUE)
        rules = str_split(ruleRow$rule, ",")[[1]]
        ruleTrue = str_split(rules[1], ":")[[1]][2] == "TRUE"
        ruleFalse = str_split(rules[2], ":")[[1]][2] == "TRUE"
        if (rowVarValue == TRUE) {
          vote = ruleTrue
        } else {
          vote = ruleFalse
        }
        votes = c(votes, vote)
        #print(paste0("FAC, x['", var, "'] ", ruleRow$rule, " (", rowVarValue, "=", vote, ")" ))
      }
    }
    #print(length(votes[votes == TRUE])/length(votes))
    return(length(votes[votes == TRUE])/length(votes))
  })
  return(resultDf)
}




