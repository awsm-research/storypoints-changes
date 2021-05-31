# sp-change-replication

A replication package for a journal paper entitled "Understanding and Predicting the change of Story Points".

This package contains two main parts:

1) 'data'
- issueLists: List of work items used in this study.
- reverted: Information extracted from JIRA for each projects. The information of each work items is reverted to the time when they are assigned to a sprint (+1 hour).
- rq2_inprogresstime: the duration that a work item's status was set to 'in progress'
- rq3_coding_result: Outcome of our open coding and manual classifications in RQ3-Qualitative 
- rq4_feature: Feature extracted from work items to be used to train the classifiers in RQ4
- rq4_metrics_influence: A list of metrics ranked by how much each of them influencing the model's prediction.
- rq4_model: Trained models in RQ4.
- rq4_performance: Performance of the model (for each model and mean value). The performance is calculated for each of the 100 models (built based on each bootstrap dataset).

2) 'rscript'
- RQ2_relationship_sp_devtime.R: Statistics summary of our models in RQ2.
- RQ3_infochange.R: The contingency table of the work items with changed story points and the work items with changed information (with Fisher's exact test).
- RQ4.R: Builds OneR (baseline), Logistic Regression, Support Vector Machine, and Random Forest classifiers
- RQ4_performance_measure.R: Measures the model's performance of our classifiers.
- RQ4_performance_comparison.R: Evaluate the model's performance of our classifiers with baseline approach.
- RQ4_influence_metric_odds.R: Observe the effect of the change of reporter-stable-sp-rate to the odds (likelihood) that the SP will be changed. 
- RQ4_influence_metric_varrank.R: Examining the influence of metrics in the models


REQUIREMENTS

* We recommend to run the code on Anaconda.
* Requires R version 3.6.1 (2019-07-05)

On Conda command line, run:
- conda install -c r rstudio

Main package required: 
- dplyr 0.8.0.1
- plyr 1.8.4
- car 3.0-10
- multilevel 2.6
- lmtest 0.9-38
- nlme 3.1-139
- scales 1.0.0
- effsize 0.8.1
- tidyr 0.8.3
- fps 2.2-9
- rms 6.2-0
- randomForest 4.6-14
- e1071 1.7-7
- ROCR 1.0-11
- caret 6.0-88
- OneR 2.2
- coin 1.4-1
- rlang 0.4.11 
- ScottKnottESD 2.0.3

After that, install Rnalytica 0.1.1 using these following instructions:
1) Install R package "devtools" 2.4.1
2) Install DMwR from the file by ...
- RUN R command: install.packages('[PATH_TO_THIS_PACKPAGE]/DMwR_0.4.1.tar.gz', repos = NULL, type="source")
3) Install Rnalytica using devtools by ..
- RUN R command: devtools::install_github('software-analytics/Rnalytica')


Full list of installed packages can be found in R-requirements.csv

