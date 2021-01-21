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

OP1Estimator <- purrr::partial(OPKEstimator,k=1)
OP2Estimator <- purrr::partial(OPKEstimator,k=2)
OP5Estimator <- purrr::partial(OPKEstimator,k=5)

PEstimator <- function(Rsquared,N,p){
  factor1 <- (N-3)*(1-Rsquared)/(N-p-1)
  factor2 <- 1+2*(1-Rsquared)/(N-p-2.3)
  return(1-factor1*factor2)
}

WEstimator <-function(Rsquared,N,p){
  return(1-(N-1)/(N-p)*(1-Rsquared))
}

SEstimator <-function(Rsquared,N,p){
  return(1-N/(N-p)*(1-Rsquared))
}

CEstimator <- function(Rsquared,N,p){
  factor1 <- (N-4)*(1-Rsquared)/(N-p-1)
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
    res <- stats::optim(Rsquared,likelihoodFunction,method="Brent",lower=0,upper=1)
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

  if (length(stats::coef(lmOut))!=lmOut$rank){
    stop("At least one indepedent variable is a linear combination of the remaining independent variables")
  }

}


#' Obtain estimates of the multiple squared correlation
#'
#' Returns different estimates of the multiple squared correlation.
#'
#' @param lmOut object of class "lm" as returned by the function \code{\link{lm}}
#
#' @return A named vector with the different estimates
#' @examples
#'## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
#'## Page 9: Plant Weight Data.
#'ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#'trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#'group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#'weight <- c(ctl, trt)
#'lm.D9 <- lm(weight ~ group)
#'estimates <- altR2(lm.D9)
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
  result <- estimate_adj_R2(Rsquared, N, p)
  return(result)
}

estimate_adj_R2 <- function(Rsquared, N, p){
  result <- numeric(20)
  esNames <- c("Rsquared","Smith","Ezekiel","Wherry","Olkin_Pratt_K_1","Olkin_Pratt_K_2", "Olkin_Pratt_K_5", "Pratt", "Claudy", "Olkin_Pratt_Exact","Maximum_Likelihood")
  esNames <- c(esNames,paste0(setdiff(esNames,c("Maximum_Likelihood","Rsquared")),'_Positive'))
  names(result) <- esNames
  result["Rsquared"] <- Rsquared
  result["Smith"] <- SEstimator(Rsquared,N,p)
  result["Ezekiel"] <- 1-(N-1)/(N-p-1)* (1-Rsquared)
  result["Wherry"] <- WEstimator(Rsquared,N,p)
  result["Olkin_Pratt_K_1"] <- OP1Estimator(Rsquared,N,p)
  result["Olkin_Pratt_K_2"] <- OP2Estimator(Rsquared,N,p)
  result["Olkin_Pratt_K_5"] <- OP5Estimator(Rsquared,N,p)
  result["Pratt"] <- PEstimator(Rsquared,N,p)
  result["Claudy"] <- CEstimator(Rsquared,N,p)
  result["Olkin_Pratt_Exact"] <- OPExactEstimator(Rsquared,N,p)
  result["Maximum_Likelihood"] <- mlEstimator(Rsquared,N,p)
  positiveEstimators <- names(result)[grepl("*_Positive",names(result))]
  for (posEst in positiveEstimators){
    normalEst <- gsub("_Positive","",posEst)
    result[posEst] <- ifelse(result[normalEst]>=0,result[normalEst],0)
  }
  return(result)
}
