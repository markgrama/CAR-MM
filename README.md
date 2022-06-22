# Structure induced by a multiple membership transformation on the conditional autoregressive model

This repo contains the code to reproduce the simulation study and the data analysis in the paper.
*In each script it is necessary to set the working directory to the repo*

## Simulation study

* **DATA** the folder `Data` contains the script to generate the data, namely the $10^4$ datasets
* **MODELS** the RStan models for both parameterisations are in the folder `Models`
* The script `R_CALL.R` contains calls the data and MCMC samplers for one simulation under one of the parameterisations used for generating the data. If you wish to launch the MCMC as a single batch, please contact me and I will share bash scripts to do it. 
* `RESULTS_SBC_PLOTS.R` creates the calibration plots
* `RESULTS_BIAS.R` retrieves results for Bias, Absolute Bias and RMSE

## Data Analysis

All the code for the Data Analysis is contained in the `DataAnalysis` folder, included the RStan models with Negative Binomial likelihood. The script `SEL_analysis.R` calls all the MCMC samplers and computes the posterior checks used in the paper.


