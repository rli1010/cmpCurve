Test Version<br>

The example only use 30 perturbation for quick run, set ptbTime=400 or larger for actual run;
## Install
library(devtools)<br>
install_github("rli1010/cmpCurve")<br><br>


## Load necessary libraries;
lapply(c("survival", "Hmisc", "rms", "cmprsk", "readr"), require, character.only = TRUE) <br>


library(cmpCurve)<br>

data(example_data)<br><br>
cmpCurve(Surv(Y,delta) ~ Z1+Z2, data = example_data, cvtimes=5, tau=4,ptbTime=30,option="both")
