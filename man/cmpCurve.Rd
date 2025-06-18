\name{cmpCurve}
\alias{cmpCurve}
\title{Predictiveness Curve for Competing Risks Regression Model}
\description{
Predictiveness curve for cumulative incidence based on the Fine & Gray regression model.
}
\usage{
cmpCurve(formula, data,  tau, option="rcs", cvtimes=10, df = 3, ptbTime=400,seed0=1234)
}
\arguments{
  \item{formula}{Survival formula in the form of Surv(Yvar,delta)~Z1+Z2, where delta takes value in 0,1,2, with 0=random censoring, 1=event of interest, 2=all other competing risks events. Used for fitting the model.}
  \item{data}{a dataset with the variables.}
  \item{tau}{predictive horizon.}
  \item{option}{"rcs",or "glm", or "both". Used for plotting.}
  \item{cvtimes}{number of repetitions for cross-validation.}
  \item{df}{df for rcs, usually 3-5.}
  \item{ptbTime}{number of perturbations}
  \item{seed0}{random seed for cross-validation data split and subsequent perturbation.}
}
\details{
Returns the estimated predictiveness curve.
}
\value{
  \item{Result}{Result with estimated curve and 95\% CI for a range of v.}}
}
\references{

}
\author{}
\examples{
#the example only use 30 perturbation for quick run, set ptbTime=400 or larger for actual run;
# Load necessary libraries;
lapply(c("survival", "Hmisc", "rms", "cmprsk", "readr"), require, character.only = TRUE)


data(example_data)
cmpCurve(Surv(Y,delta) ~ Z1+Z2, data = example_data, cvtimes=5, tau=4,ptbTime=30,option="both")

}
