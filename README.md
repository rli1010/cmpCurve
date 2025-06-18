## Install
library(devtools)<br>
install_github("rli1010/cmpCurve")<br><br>

## Example
#### The example use 30 perturbation for quick run, set ptbTime=400 or larger for actual run;

lapply(c("survival", "Hmisc", "rms", "cmprsk", "readr"), require, character.only = TRUE) <br>


library(cmpCurve)<br>

data(example_data)<br><br>
cmpCurve(Surv(Y,delta) ~ Z1+Z2, data = example_data, cvtimes=5, tau=4,ptbTime=30,option="both")<br>

## Ref
Tao, W., Ning, J., Li, W., Chan, W., Luo, X., Li, R., Evaluating the Predictiveness Curve for Risk Prediction Models with Competing Risks Data. Manuscript.  
