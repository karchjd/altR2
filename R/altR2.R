hypergeoApprox <- function(Rsquared,N,p,klimit){
  sum <- 1
  for (k in 1:klimit){
    factor1 <- gamma(1+k)^2*gamma((N-p+1)/2)/gamma((N-p+1)/2+k)
    factor2 <- (1-Rsquared)^k/factorial(k)
    sum <- sum+factor1*factor2
  }
  return(sum)
}

OPKEstimator <- function(Rsquared,N,p,k){
  factor1 <- (N-3)/(N-p-1)*(1-Rsquared)
  factor2 <- hypergeoApprox(Rsquared,N,p,k)
  return(1-factor1*factor2)
}

OP2Estimator <- purrr::partial(OPKEstimator,k=2)
OP5Estimator <- purrr::partial(OPKEstimator,k=5)

PEstimator <- function(Rsquared,N,p){
  factor1 <- (N-3)*(1-Rsquared)/(N-p-1)
  factor2 <- 1+2*(1-Rsquared)/(N-p-2.3)
  return(1-factor1*factor2)
}

OPExactEstimator <- function(Rsquared,N,p){
  factor1 <- (N-3)/(N-p-1)*(1-Rsquared)
  factor2 <- gsl::hyperg_2F1(1,1,(N-p+1)/2,1-Rsquared)
  res <- 1-factor1*factor2
  return(res)
}

checkInput <- function(lmOut,N,p){
  if (class(lmOut)!="lm"){
    stop("lmOut must be an output of the lm function")
  }

  if (!attr(lmOut$terms, "intercept")){
    stop("The linear model must contain an intercept")
  }

  if (length(coef(lmOut))!=lmOut$rank){
    stop("At least one indepedent variable is a linear combination of the remaining independent variables")
  }

}



#' @importFrom purrr partial
#' @importFrom gsl hyperg_2F1
#' @export
altR2 <- function(lmOut) {
  checkInput(lmOut,N,p)

  #get needed quantities
  lmSum <-summary(lmOut)
  Rsquared <- lmSum$r.squared
  N <- nrow(lmOut$model)
  p <- lmOut$rank-1


  if (!(N-p)>=2){
    stop(sprintf("Sample size %d is not at least 2 bigger than number of predictors %d",N,p))
  }

  #create results by calling the respective shrinkage functions
  result <- numeric(10)
  esNames <- c("Ezekiel","Olkin & Pratt, K=2", "Olkin & Pratt, K=5", "Pratt", "Olkin & Pratt, Exact")
  esNames <- c(esNames,paste0('Positive, ',esNames))
  names(result) <- esNames
  result[1] <- lmSum$adj.r.squared
  result[2] <- OP2Estimator(Rsquared,N,p)
  result[3] <- OP5Estimator(Rsquared,N,p)
  result[4] <- PEstimator(Rsquared,N,p)
  result[5] <- OPExactEstimator(Rsquared,N,p)
  for (i in 6:10){
    result[i] <- ifelse(result[i-5]>=0,result[i-5],0)
  }
  return(result)
}