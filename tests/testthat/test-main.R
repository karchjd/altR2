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
  OP50Estimator <- purrr::partial(OPKEstimator,k=50)
  N <- nrow(x$model)
  p <- x$rank-1
  op50res <- OP50Estimator(summary(x)$r.squared,N,p)
  expect_equivalent(op50res,normalRes["Olkin & Pratt, Exact"])
})
