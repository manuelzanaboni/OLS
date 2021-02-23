# Initial package installation
install.packages(c("tibble", "rstudioapi", "knitr"))

# Load needed packages
library("rstudioapi") # automatic working directory setup
library("tibble") # tibble function output
library("knitr") # table function output

# Path, starting from project root, in which data is stored
dataPathFromRoot <- "/wages_in_Belgium/bwages.dat"

# Automatic working directory setup and data loading
setwd(dirname(getActiveDocumentContext()$path))
dataPath <- paste0(getwd(), dataPathFromRoot)
data <- read.delim(dataPath, header = TRUE, sep = "\t")

# OLS - main function
# performs a complete ordinary least squares computation
# INPUT:  dependent - (matrix-like object) dependent variable, ground truth to be fitted
#         independent - (matrix-like object) independent variables, regressors
#         useQRdecomposition - (logical) whether to use QR-decomposition or not, in order to compute coefficients
#         returnTable - (logical) whether to return a formatted table as an output or tibble-like data
#         skipDataValidityTest - (logical) whether to skip initial data validity check
#         robustErrors - (logical) whether to calculate heteroskedasticity-robust errors
# OUTPUT: table-like or tibble-like regression results (coefficients, errors, statistics, fitted and residuals' values)

OLS <- function(dependent, independent, 
                useQRdecomposition = TRUE,
                returnTable = FALSE, 
                skipDataValidityTest = FALSE,
                robustErrors = FALSE){
  
  # regressors (independent variables)
  X <- as.data.frame(independent)
  
  if(!skipDataValidityTest)
    TestDataValidity(X)
  
  X <- as.matrix(data.frame(cbind(INTERCEPT = 1, X)))
  
  # dependent variable
  y <- as.matrix(dependent)
  
  if(useQRdecomposition) {
    # solve using formula: b = R^-1 * Q' * y
    qrResult <- QRdecomposition(X)
    b <- as.vector(backsolve(qrResult$R, t(qrResult$Q) %*% y))
  }
  else {
    # solve using formula: b = (X'X)^−1 * X'y
    Xtransposed <- t(X)
    XtransposedX <- Xtransposed %*% X
    
    # calculate the inverse of X'X
    XtransposedXinverse <- solve(XtransposedX)
    
    # calculate X′y
    Xtransposedy <- Xtransposed %*% y
    
    # calculate regression coefficients
    b <- XtransposedXinverse %*% Xtransposedy
  }
  
  # fitted values
  fitted <- X %*% b
  
  # residuals
  residuals <- y - fitted
  
  # calculate statistics
  if(useQRdecomposition)
    stats <- CalculateStatistics(X, y, fitted, residuals, b, R = qrResult$R)
  else
    stats <- CalculateStatistics(X, y, fitted, residuals, b, robustErrors, 
                                 XtransposedXinverse = XtransposedXinverse)
  
  # build return values
  # vars
  vars <- tibble(
    term = colnames(X),
    coefficients = as.numeric(b),
    standardError = as.numeric(stats$standardError),
    tStatistic = as.numeric(stats$tStatistic),
    pValue = as.numeric(stats$pValue)
  )
  # statistics
  statistics <- tibble(
    ResidualsStandardError = as.numeric(stats$RSE),
    R2 = as.numeric(stats$R2),
    R2adj = as.numeric(stats$R2adj),
    fStatistic = as.numeric(stats$fStatistic)
  )
  # values (fitted and residuals)
  values <- list("fitted" = fitted, 
                 "residuals" = residuals)
  if(returnTable) {
    vars$pValue <- format(vars$pValue, digits = 3)
    return(kable(list(vars, statistics), "simple"))
  } else
      return(list("vars" = vars,
                  "statistics" = statistics, 
                  "values" = values))
}

# TestDataValidity - auxiliary function
# checks the validity of the data and gives errors and warnings
# INPUT:  data - (matrix-like object) regressors
# OUTPUT: NULL

TestDataValidity <- function(data) {
  # Firstly test observation num. vs variables num.
  if(nrow(data) < ncol(data))
    stop("Number of observations is less than number of variables.")
  
  # Test for special values (NA, Inf)
  special <- apply(data, 2, function(x) any(is.na(x) | is.infinite(x)))
  if(any(special == TRUE))
    warning("Data contains special values.")
    
  # Bivariate correlation
  if(ncol(data)>1)
    TestCorrelation(data)
}

# TestCorrelation - auxiliary function
# tests bivariate correlation
# INPUT:  data - (matrix-like object) regressors
# OUTPUT: NULL

TestCorrelation <- function(data) {
  names <- colnames(data)
  observationsNumber <- nrow(data)
  correlations <- NULL
  
  for(var1 in names) {
    for(var2 in names) {
      if(var1 != var2) {
        val1 <- as.numeric(unlist(data[var1]))
        val2 <- as.numeric(unlist(data[var2]))
        numerator <- sum((val1 - mean(val1, na.rm = TRUE)) * (val2 - mean(val2, na.rm = TRUE)))
        r <- (1 / (observationsNumber - 1)) * (numerator / (sd(val1, na.rm = TRUE) * sd(val2, na.rm = TRUE)))
        correlations <- rbind(correlations, r)
      }
    }
  }
  
  for(val in correlations) {
    if(val > 0.8 || val < -0.8)
      stop("Data is highly correlated.")
  }
}

# QRdecomposition - auxiliary function
# performs QR decomposition over regressors data
# INPUT:  X - (matrix-like object) regressors
# OUTPUT: an orthogonal matrix Q and an upper triangular matrix R
# code found at the following link: https://rpubs.com/akshatdwivedi/qr

QRdecomposition <- function(X) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  q1 <- matrix(0, nrow = n, ncol = p)
  r1 <- matrix(0, nrow = p, ncol = p) 
  u <- matrix(0, nrow = n, ncol = p)
  
  u[, 1] <- X[, 1]
  
  for (k in 2:ncol(X)) {
    
    u[, k] <- X[, k]
    
    # successive orthogonalization
    for (ctr in seq(1, (k - 1), 1)) {
      # classical GS
      u[, k] <- u[, k] - ((sum(u[, ctr] * X[, k]) / sum(u[, ctr] * u[, ctr])) * (u[, ctr]))
    }
  }
  
  # dividing each column by respective vector (l2) norm (rescaling each vector to have length 1)
  q1 <- apply(u, 2, function(x) { x / sqrt(sum (x %*% x)) })
  r1 = t(q1) %*% X

  return(list(Q = q1, R = r1))
}

# CalculateStatistics - auxiliary function
# calculates regression statistics
# INPUT:  X - (matrix-like object) independent variables, regressors
#         y - (matrix-like object) dependent variable
#         fitted - (vector-like object) fitted values
#         residuals - (vector-like object) residuals values
#         b - (vector-like object) regression coefficients
#         robustErrors - (logical) whether to calculate heteroskedasticity-robust errors
#         XtransposedXinverse - (matrix) result of (X'X)^−1, using classical formula
#         R - (matrix) upper triangular matrix, using QR decomposition
# OUTPUT: list containing: R2, R2 adjusted, Residuals Standard Error, coefficients standard errors,
#                          t-statistics, F-statistics and p-values

CalculateStatistics <- function(X, y, 
                                fitted,
                                residuals,
                                b,
                                robustErrors = FALSE,
                                XtransposedXinverse = NULL, 
                                R = NULL) {
  observationsNumber <- nrow(X)
  variablesNumber <- ncol(X)
  degreesOfFreedom <- observationsNumber - variablesNumber
  
  # R2
  SSR <- sum(residuals^2)
  SST <- sum((y - mean(y, na.rm = TRUE))^2)
  R2 <- 1 - (SSR / SST)
  
  # Adjusted R2
  R2adj <- 1 - (((1 - R2) * (observationsNumber - 1)) / (degreesOfFreedom - 1))
  
  SE2 <- SSR / degreesOfFreedom
  
  # Residual standard error
  RSE <- sd(residuals)
  
  # Standard Errors
  if(robustErrors) {
    resDiag <- diag(c(matrix(residuals^2, 1)))
    
    # heteroscedasticity-corrected covariance matrix
    V <- XtransposedXinverse %*% t(X) %*% resDiag %*% X %*% XtransposedXinverse
  }
  else {
    if(is.null(XtransposedXinverse) && !is.null(R)) {
      # covariance matrix
      V <- SE2 * chol2inv(R)
    } else {
      # covariance matrix
      V <- XtransposedXinverse * SE2
    }
  }
  stdErrors <- sqrt(diag(V))

  # t-statistic
  tStatistic <- as.numeric(b) / as.numeric(stdErrors)
  
  # F-statistic
  fStatistic <- (R2 / (variablesNumber - 1)) / ((1 - R2) / degreesOfFreedom)
  
  # P-value 2 tailed
  pValue <- 2 * pt(abs(tStatistic), degreesOfFreedom, lower.tail = FALSE)
  
  return(list("R2" = R2,
              "R2adj" = R2adj,
              "RSE" = RSE,
              "standardError" = stdErrors,
              "tStatistic" = tStatistic,
              "fStatistic" = fStatistic,
              "pValue" = pValue))
}

# CustomPlot - auxiliary function
# wrapper function to customize fitted-residuals plot
# INPUT:  fitted - (vector-like object) fitted values
#         residuals - (vector-like object) residuals values
# OUTPUT: NULL

CustomPlot <- function(fitted, residuals) {
  plot(fitted, residuals, pch = 16, xlab = "Fitted", ylab = "Residuals")
}

##### Initial functionality test on simple data ##### 

# data from https://machinelearningmastery.com/solve-linear-regression-using-linear-algebra/
X <- matrix(c(0.05, 0.18, 0.31, 0.42, 0.5))
y <- matrix(c(0.12, 0.22, 0.35, 0.38, 0.49))
res <- OLS(y, X, useQRdecomposition = FALSE)
res$vars
plot(X, y)
abline(coef = res$vars$coefficients, col = "red")

# Table output
OLS(y, X, useQRdecomposition = FALSE, returnTable = TRUE)

# Less observations than variables test
X1 <- matrix(c(0.05))
X2 <- 2 * X1
X <- matrix(c(X1,X2), nrow = 1)

tryCatch({
  OLS(y, X)
}, error = function(e) {
  message("Error catched:")
  message(e)
})

# Collinearity test
X1 <- matrix(c(0.05, 0.18, 0.31, 0.42, 0.5))
X2 <- 2 * X1
X <- matrix(c(X1,X2), nrow = 5)

tryCatch({
  OLS(y, X)
}, error = function(e) {
  message("Error catched:")
  message(e)
})

##### Wages in Belgium - Data summary #####
str(data)
summary(data)

hist(data$WAGE, xlab = "WAGE", main = "Histogram of WAGE", col = "green")

hist(data$EXPER, xlab = "EXPER", main = "Histogram of EXPER")

table(data["MALE"] == 1) # 893 males, 579 females

plot(data$EDUC, ylab = "EDUC")

cor(data)

## First regression of WAGE on MALE
y <- subset(data, select = WAGE)
X <- subset(data, select =  MALE)
OLS(y, X, returnTable = TRUE)


# New FEMALE variable
data$FEMALE <- as.integer(!as.logical(data$MALE))

## Regress of WAGE on MALE and FEMALE, variables are correlated (singular)
X <- subset(data, select = c(MALE, FEMALE))

tryCatch({
  OLS(y, X)
}, error = function(e) {
  message("Error catched:")
  message(e)
})

# check with R built-in lm() function
# FEMALE variable gets ignored
lm(WAGE ~ MALE + FEMALE, data)


## Regress WAGE on MALE, EDUC, EXPER, EXPER^2
data$EXPER2 <- data$EXPER ^ 2
X <- subset(data, select = c(MALE, EDUC, EXPER, EXPER2))

tryCatch({
  OLS(y, X)
}, error = function(e) {
  message("Error catched:")
  message(e)
})
# EXPER and EXPER^2 are obviously correlated, so let's skip data validity test
out <- OLS(y, X, skipDataValidityTest = TRUE)
out$vars
out$statistics

# Plot the residuals versus the fitted values
# increased variation in the residuals for higher fitted values 
# --> heteroskedasticity
CustomPlot(out$values$fitted, out$values$residuals)


## Regress LNWAGE on MALE, LNEDUC, LNEXPER, LNEXPER^2
data$LNEXPER2 <- data$LNEXPER ^ 2
y <- subset(data, select = LNWAGE)
X <- subset(data, select = c(MALE, LNEDUC, LNEXPER, LNEXPER2))

tryCatch({
  OLS(y, X)
}, error = function(e) {
  message("Error catched:")
  message(e)
})
# LNEXPER and LNEXPER^2 are obviously correlated, so let's skip data validity test
out <- OLS(y, X, skipDataValidityTest = TRUE)
out$vars
out$statistics

# Plot the residuals versus the fitted values
# heteroskedasticity still present, but residual variation is much less pronounced
CustomPlot(out$values$fitted, out$values$residuals)


## Heteroskedasticity robust errors
y <- subset(data, select = WAGE)
X <- subset(data, select = c(MALE, EDUC, EXPER, EXPER2))
out <- OLS(y, X, skipDataValidityTest = TRUE, 
           useQRdecomposition = FALSE, robustErrors = TRUE)
# Check corrected Std. errors
out$vars
out$statistics

# OneHotEncode
# performs one hot encoding over a categorical (factors) variable
# INPUT:  variable - (vector-like object) categorical variable
#         oldName - (vector-like object) input variable's name
# OUTPUT: one hot encoded variable (dummy) representation with headers vector
#         (data frame of n columns, n = number of unique values)

OneHotEncode <- function(variable, oldName) {
  uniqueValues <- unique(variable)
  encoded <- NULL
  names <- NULL
  for(val in uniqueValues) {
    newVar <- as.integer(variable == val)
    encoded <- cbind(encoded, newVar)
    names <- rbind(names, paste0(oldName, "_", val))
  }
  encoded <- as.data.frame(encoded)
  colnames(encoded) <- names
  return(list("encoded" = encoded, "names" = as.vector(names)))
}

## Experiment with EDUC
encodedData <- OneHotEncode(data$EDUC, "EDUC")
data <- cbind.data.frame(data, encodedData$encoded)

y <- subset(data, select = WAGE)
# consider only 4 of the 5 new variables
dummyVariables <- encodedData$names[2:length(encodedData$names)]

# compute cross terms EDUC - MALE
crossTerms <- list()
for(var in dummyVariables) {
  header <- paste0(var, "_", "MALE")
  crossTerms[header] <- data[var] * data$MALE
}
crossTerms <- as.data.frame(crossTerms)
crossTermsHeaders <- colnames(crossTerms)
data <- cbind.data.frame(data, crossTerms)

# final linear regression 
finalVariables <- c(dummyVariables, "MALE", "EXPER", "EXPER2", crossTermsHeaders)

X <- subset(data, select = finalVariables)
out <- OLS(y, X, skipDataValidityTest = TRUE,
           useQRdecomposition = FALSE, robustErrors = TRUE)
out$vars
out$statistics