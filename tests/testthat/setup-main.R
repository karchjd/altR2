library(MASS)
MyDataGeneration <- function(n,rho,mc,p,allCor){
  varResidual <- 10
  theCov <- mc
  theVar <- 1
  popSigma <- matrix(theCov,nrow=p,ncol=p)
  diag(popSigma) <-  theVar
  X <- mvrnorm(n,rep(0,p),popSigma)
  if (rho==0){
    varY <- varResidual
  }else{
    varY <- varResidual/(1-rho)
  }
  varPreds <- varY-varResidual
  if (allCor){
    oneCoef <- sqrt(varPreds/(p*theVar+(p^2-p)*theCov))
    coefs <- as.matrix(rep(oneCoef,p))
  }else{
    firstCoef <- sqrt(varPreds/theVar)
    coefs <- as.matrix(c(firstCoef,rep(0,p-1)))
  }
  y <- 100+X%*%coefs+rnorm(n,sd=sqrt(varResidual))
  return(list(y=y,X=X,p=p))
}
set.seed(21399)
testData <- MyDataGeneration(20,0.01,0,18,TRUE)
testData <- as.data.frame(cbind(testData$y,testData$X))
x <- lm(V1 ~ .,data=testData)
normalRes <- altR2(x)
