# RSVMLLSS
Here you can find the R scripts, RDS files, csv datasets, used in Machine Learning–Based Limited Sampling Strategy to Predict Rosuvastatin Exposure

Combinations is the script used for testing different ML algorithms, combination sizes, concentrations at different times, inclusion of covariates and derived predictors from concentrations, developed for predicting RAL AUC from 3 samples. In this script the algorithms available are xgboost, glm, rf, and svm.

RSVscript contains the POPPK models implementation, PK profiles and patients simulation, datasets management, plots, metrics and initial feature selection

The RDS file named "model_rfrsv2021.rds" is the RDS file of the xgboost model developed in this study and the "explainer_externalrfrsv2021.rds" was used for the shiny webapplication development. The RDS can be loaded directly to get predictions from other datasets run the script and get the RSV AUC predictions.

The CSV file named "rsvsim2021.csv" is the simulated patients dataset used for training the ML models. 
The CSV files named "" and "" are the datasets, simulated and real world patients, repectivelly, used for the external validation of the ML models.
