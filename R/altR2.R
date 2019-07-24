hypergeoApprox <- function(Rsquared,N,p,klimit){
  sum <- 1
  for (k in 1:klimit){
    factor1 <- gamma(1+k)^2
    factor2 <- 1/prod(seq(from=(N-p+1)/2,by=1,to=(N-p+1)/2+k-1))
    factor3 <- (1-Rsquared)^k/factorial(k)
    sum <- sum+factor1*factor2*factor3
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

#implements effective marginal likelihood function as defined on page 174 in https://doi.org/10.1111/j.1467-842X.1985.tb00559.x
effectiveLH <- function(Rsquared,N,p,rhoSquared){
  -1*(1-rhoSquared)^(N/2)*gsl::hyperg_2F1(0.5*N,0.5*N,0.5*p,rhoSquared*Rsquared)
}


mlEstimator <- function(Rsquared,N,p){
  if(Rsquared<=p/N){
    return(0)
  }else{
    likelihoodFunction <- purrr::partial(effectiveLH,Rsquared=Rsquared,N=N,p=p)
    res <- optim(Rsquared,likelihoodFunction,method="Brent",lower=0,upper=1)
    return(res$par)
  }
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
  result <- numeric(11)
  esNames <- c("Ezekiel","Olkin_Pratt_K_2", "Olkin_Pratt_K_5", "Pratt", "Olkin_Pratt_Exact","Maximum_Likelihood")
  esNames <- c(esNames,paste0(setdiff(esNames,"Maximum_Likelihood"),'_Positive'))
  names(result) <- esNames
  result[1] <- lmSum$adj.r.squared
  result[2] <- OP2Estimator(Rsquared,N,p)
  result[3] <- OP5Estimator(Rsquared,N,p)
  result[4] <- PEstimator(Rsquared,N,p)
  result[5] <- OPExactEstimator(Rsquared,N,p)
  result[6] <- mlEstimator(Rsquared,N,p)
  for (i in 7:11){
    result[i] <- ifelse(result[i-6]>=0,result[i-6],0)
  }
  return(result)
}
