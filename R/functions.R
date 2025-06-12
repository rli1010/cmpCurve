# Load necessary libraries
#library(survival)
#library(Hmisc)
#library(rms)
#library(cmprsk)
#library(readr)

expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

# Main function
cmpCurve <- function(formula, data, option="rcs", tau, cvtimes=10, df = 3, ptbTime=400,seed0=1234) {
  ##For event type, code the cause of interest as 1 and censor as 0, other causes as 2
  ##df is for RCS splines, commonly 3-5
  ##cvtimes is the number of repetitions for cross-validation
  ##ptbTime is the number of perturbations
  ##option can be rcs of glm for parameterization

  call <- match.call()
  formula <- as.formula(formula)
  Time=as.character(formula[[2]][[2]]) ###extract Y from formula
  status=as.character(formula[[2]][[3]])

  # Define option ranges
  options_list <- list(rcs = 1:99, glm = 100:198) #the first 1-99 is rcs with grid seq(0.01,0.99,0.01)

  # Validate and map the option
  if (is.character(option) && option %in% names(options_list)) {
    option <- options_list[[option]]
  } else if (!is.numeric(option) || !(option %in% unlist(options_list))) {
    stop("Invalid option. Choose 'rcs' or 'glm'.")}

  # point estimate, seed0 used for cross-validation and subsequent perturbation
  set.seed(seed0)
  predict <- rpt.risk(data, ptb.weight = rep(1, nrow(data)), cvtimes, formula, df, tau, Time, status)[option]

  # perturbation ptbTime=400 times and create CI
  data_rep <- replicate(ptbTime, rpt.risk(data, ptb.weight = rexp(nrow(data), rate = 1),
                                          cvtimes, formula, df, tau, Time, status)[option])
  df <- data.frame(matrix(unlist(data_rep), nrow = ptbTime, byrow = TRUE))
  trans <- logit(df)
  se.trans <- apply(trans, 2, sd)
  lower_bound_trans <- logit(predict) - 1.96 * se.trans
  upper_bound_trans <- logit(predict) + 1.96 * se.trans
  lower_bound <- expit(lower_bound_trans)
  upper_bound <- expit(upper_bound_trans)
  #combine result
  result <- data.frame(v = seq(0.01, 0.99, by = 0.01), Risk = predict, L_bound = lower_bound, U_bound = upper_bound,
                       SE=apply(df,2,sd))

  # Plotting
  plot(result$v, result$Risk, type = "l", col = "blue", lty = 1, ylim = c(0, max(upper_bound)),
       xlab = "v", ylab = "Risk", main = "Predictiveness Curve")
  lines(result$v, result$L_bound, col = "black", lty = 2)
  lines(result$v, result$U_bound, col = "black", lty = 2)
  set.seed(NULL) #free the seed
  return(result)
}


# Generate test data, using train to generate test cause1.Xi
test.data <- function(train, test, formula, df, Time, status) {
  event <- factor(train[[status]], levels = 0:2, labels = c("censor", "death", "caus2"))

  pdata <- finegray(Surv(train[[Time]], event) ~ ., data = train)  ###use the formula to fit FG
  pdata$pdata.weight <- pdata$fgwt * pdata$ptbweight

  formula_rhs <- as.formula(paste("~", paste(attr(terms(formula), "term.labels"), collapse = " + ")))
  fg <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                               paste(attr(terms(formula_rhs), "term.labels"), collapse = " + "))),
              weights = pdata$pdata.weight, data = pdata)

  formula0=paste(Time,as.character(formula_rhs),collapse="")
  test_matrix <- model.matrix(stats::formula(formula)[-2], data = test)[, -1]  # remove intercept
  # Calculate caus1.xi using matrix multiplication
  caus1.xi <- test_matrix %*% fg$coefficients
  return(cbind(test, caus1.xi))
}

# Generate numerator w, mintime of Y and tau for denominator, inverse delta for censoring, and ptb weight
weight <- function(data, ptb.weight, tau, Time, status) {
  mintime <- pmin(data[[Time]],tau)
  w1.1 <- (data[[Time]]<=tau) & (data[[status]]>0)
  w1.2 <- (data[[Time]]>tau)
  w1 <- (w1.1+w1.2)
  w.data <- cbind(data, w1, mintime)

  w.data$ivsdlt <- (w.data[[status]] == 0) + 0
  w.data$ptbweight <- ptb.weight
  fit <- survfit(Surv(w.data[[Time]], w.data$ivsdlt) ~ 1, weights = ptbweight, data = w.data)
  g.hat <- stepfun(x = fit$time, y = c(1, fit$surv))
  w.data$g.hat <- pmax(g.hat(w.data$mintime),0.025)

  w.data$ivs.weight <- (w.data$w1) / (w.data$g.hat)
  w.data$Event1 <- (w.data[[Time]] <= tau & w.data[[status]] == 1) + 0
  w.data$weight <- w.data$ivs.weight * w.data$ptbweight
  wt.data <- w.data[, !names(w.data) %in% c("w1", "mintime", "ivsdlt", "Tt")]
  return(wt.data)}

# RCS model, for testing data, get R(v) for each quantile point
RCS.weight <- function(train, test, formula, df, Time, status) {
  ptb.data <- test.data(train, test, formula, df, Time, status) # Test data
  nwdata <- data.frame(caus1.xi = wtd.quantile(ptb.data$caus1.xi, weights = ptb.data$ptbweight,
                                               probs = c(seq(0.01, 0.99, by = 0.01))))
  d <- datadist(ptb.data$caus1.xi)
  options(datadist = d)
  m3 <- suppressWarnings(lrm((Event1 == 1) ~ rcs(caus1.xi, df), weights = weight, data = ptb.data))
  ptb.LRM <- predict(m3, newdata = nwdata, type = "fitted")
  predc=cbind(seq(0.01, 0.99, by = 0.01), ptb.LRM)[, 2]
  return(predc)
}

# GLM parameterization, for testing data, get R(v) for each quantile point
GLM.weight <- function(train, test, formula, df, Time, status) {
  ptb.data <- test.data(train, test, formula, df, Time, status) # Test data
  nwdata <- data.frame(caus1.xi = wtd.quantile(ptb.data$caus1.xi, weights = ptb.data$ptbweight,
                                               probs = c(seq(0.01, 0.99, by = 0.01))))
  d <- datadist(ptb.data$caus1.xi)
  options(datadist = d)
  m1 <- suppressWarnings(lrm((Event1 == 1) ~ caus1.xi, weights = weight, data = ptb.data)) #only linear term
  ptb.GLM <- predict(m1, newdata = nwdata, type = "fitted")
  predc=cbind(seq(0.01, 0.99, by = 0.01), ptb.GLM)[, 2]
  return(predc)
}

# Split data, generate ave risk for each data
sig.risk <- function(data, ptb.weight, formula, df, tau, Time, status) {
  sim.data <- weight(data, ptb.weight, tau, Time, status)
  sample <- sample.int(n = nrow(sim.data), size = floor(0.5 * nrow(sim.data)), replace = FALSE) #two fold CV
  sim1.data <- sim.data[sample, ]
  sim2.data <- sim.data[-sample, ]

  risk1 <- RCS.weight(train = sim1.data, test = sim2.data, formula, df, Time, status)
  risk2 <- RCS.weight(train = sim2.data, test = sim1.data, formula, df, Time, status)
  RCS.risk <- rowMeans(cbind(risk1, risk2))

  risk3 <- GLM.weight(train = sim1.data, test = sim2.data, formula, df, Time, status)
  risk4 <- GLM.weight(train = sim2.data, test = sim1.data, formula, df, Time, status)
  glm.risk <- rowMeans(cbind(risk3, risk4))
  cbind(rcs = RCS.risk, glm = glm.risk)
}


# Repeat CV on same data, calculate multiple times avg
rpt.risk <- function(data, ptb.weight, cvtimes, formula, df, tau, Time, status) {
  rep.ptb <- replicate(cvtimes, sig.risk(data, ptb.weight, formula, df, tau, Time, status))
  avg.ptb <- data.frame(matrix(unlist(rep.ptb), nrow = cvtimes, byrow = TRUE))
  rpt.perturb <- colMeans(avg.ptb)
  return(rpt.perturb)
}
















