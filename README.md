Test Version<br>

#the example only use 30 perturbation for quick run, set ptbTime=400 or larger for actual run;
# Load necessary libraries;
lapply(c("survival", "Hmisc", "rms", "cmprsk", "readr"), require, character.only = TRUE)

library(cmpCurve)

data(example_data)
cmpCurve(Surv(Y,delta) ~ Z1+Z2, data = example_data, cvtimes=5, tau=4,ptbTime=30,option="both")
