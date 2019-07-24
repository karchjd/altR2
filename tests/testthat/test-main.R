OP2EstimatorManual <- function(Rsquared,N,p,k){
  factor1 <- (N-3)/(N-p-1)*(1-Rsquared)
  factor2 <- 1+2*(1-Rsquared)/(N-p+1)+8*(1-Rsquared)^2/((N-p+1)*(N-p+3))
  return(1-factor1*factor2)
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
