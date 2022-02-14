library(MASS)
OP2EstimatorManual <- function(Rsquared,N,p,k){
  factor1 <- (N-3)/(N-p-1)*(1-Rsquared)
  factor2 <- 1+2*(1-Rsquared)/(N-p+1)+8*(1-Rsquared)^2/((N-p+1)*(N-p+3))
  return(1-factor1*factor2)
}

MyDataGenerationFull <- function(n, rho, p) {
  # constants
  varResidual <- 10
  theCov <- 0.3
  theVar <- 1

  # distribution predictors
  popSigma <- matrix(theCov, nrow = p, ncol = p)
  diag(popSigma) <- theVar
  X <- mvrnorm(n, rep(0, p), popSigma)

  # get residual variance and variance of f(x) from rho and total variance (varianceY)
  if (rho == 0) {
    varY <- varResidual
  } else {
    varY <- varResidual / (1 - rho)
  }
  varPreds <- varY - varResidual

  # calculate coefficients
  oneCoef <- sqrt(varPreds / (p * theVar + p*(p-1) * theCov))
  coefs <- as.matrix(rep(oneCoef, p))
  y <- 100 + X %*% coefs + rnorm(n, sd = sqrt(varResidual))

  XY <- cbind(X, y)
  return(XY)
}

test_that("adj.R.Squared sanity", {
  N <- nrow(x$model)
  p <- x$rank-1
  expect_equivalent(normalRes["Ezekiel"],1-(N-1)/(N-1-p)*(1-summary(x)$r.squared))
 })

test_that("OPK+New Sanity", {
  OP50Estimator <- purrr::partial(altR2:::OPKEstimator,k=50)
  N <- nrow(x$model)
  p <- x$rank-1
  op50res <- OP50Estimator(summary(x)$r.squared,N,p)
  expect_equivalent(op50res,normalRes["Olkin_Pratt_Exact"])
})

test_that("All Similar", {
  set.seed(21399)
  testData <- MyDataGeneration(150,0.5,0,2,TRUE)
  testData <- as.data.frame(cbind(testData$y,testData$X))
  x <- lm(V1 ~ .,data=testData)
  normalRes <- altR2(x)
  expect(all(abs(normalRes-mean(normalRes))< 0.005),'Not all similiar')
})


test_that("Maximum likelihood sanity",{
  rhoRange <- seq(from=0,by=0.01,to=1)
  lh <- -1*altR2:::effectiveLH(.25,40,10,rhoRange)
  # plot(rhoRange,ln) visually checked against figure 2 of https://doi.org/10.3102%2F10769986027003223
  expect_equal(which.max(lh),1)
})

test_that("Positive and nonpositive are the same",{
  normalEstimators <- setdiff(names(normalRes)[!grepl("*_Positive",names(normalRes))],c("Maximum_Likelihood","Rsquared"))
  for (estimator in normalEstimators){
    print(estimator)
    expect_equivalent(normalRes[estimator],normalRes[paste0(estimator,"_Positive")])
  }
})

test_that("Positive versions are positive",{
  set.seed(9041249)
  testData <- MyDataGeneration(12,0,0,10,TRUE)
  testData <- as.data.frame(cbind(testData$y,testData$X))
  x <- lm(V1 ~ .,data=testData)
  normalRes <- altR2(x)
  posEstimators <- grepl("*_Positive",names(normalRes))
  expect(all(posEstimators)>=0,"not all positive estimators >0")
})

test_that("All different",{
  normalEstimators <- names(normalRes)[!grepl("*_Positive",names(normalRes))]
  for (estimator in normalEstimators){
    expect(all(normalRes[estimator] != normalRes[setdiff(normalEstimators,estimator)]),sprintf('%s is equal to at least one other estimator',estimator))
  }
})

test_that("new Estimator more elaborate",{
  julsSeq <- function(from,to,by){
    if (to<from){
      return(integer())
    }else{
      return(seq(from=from,to=to,by=by))
    }
  }

  julsHypergeo <- function(a,b,c,z){
    stopifnot(a==1)
    stopifnot(b==1)
    stopifnot(c %% 0.5 == 0)
    #see http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/06/07/02/0005/
    #works for c=3. if not not. with numerical problems
    if (c %% 1==0){
      prefactor <- (c-1)*z*(z-1)^{-2}
      theSum<-0
      shared <- ((z-1)/z)
      for (k in julsSeq(from=2,to=(c-1),by=1)){
        term1 <- shared^k/(c-k)
        # print(sprintf('theSum %f, term1 %f',theSum,term1))
        theSum <- theSum+term1
      }
      result <- prefactor*(theSum-shared^c*log(1-z))
    }else{
      # see https://mathoverflow.net/questions/335424/evaluate-gaussian-hypergeometric-function-2f-111cz?noredirect=1#comment837983_335424
      # works with numerical problems
      curRes <- 1/(1-z)*(1+(sqrt(z)*asin(sqrt(z)))/sqrt(1-z))
      for (i in julsSeq(from=1.5,to=c,by=1)){
        curC <- i-1
        # print(sprintf('My %f, Gnu %f',curRes,hyperg_2F1(1,1,curC,as.double(z))))
        curRes <- (curC-curC*(1-z)*curRes)/(z*(curC-1))
      }
      result <- curRes
    }
    return(result)
  }

  OPExactEstimatorJul <- function(Rsquared,N,p){
    factor1 <- (N-3)/(N-p-1)*(1-Rsquared)
    factor2 <- julsHypergeo(1,1,(N-p+1)/2,1-Rsquared)
    res <- 1-factor1*factor2
    return(res)
  }

  N <- 10
  p <- 8
  rsquaredRange <- seq(from=0.01,by=0.01,to=0.99)
  for (rsquared in rsquaredRange){
    expect_equivalent(OPExactEstimatorJul(rsquared,N,p),OPExactEstimator(rsquared,N,p),tolerance=10^-5)
  }
})

test_that("ML estimator corner case", {
  XY <- MyDataGenerationFull(500, .9, 2)
  X <- XY[,-ncol(XY)]
  Y <- XY[,ncol(XY)]
  model <- lm(Y ~ X)
  r2res <- altR2(model)
})


