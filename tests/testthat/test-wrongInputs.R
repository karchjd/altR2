test_that("wrong object", {
  x <- 4
  expect_error(altR2(x),"lmOut must be an output of the lm function")
})

test_that("no intercept", {
  x <- lm(V1 ~ . -1,data=testData)
  expect_error(altR2(x),"The linear model must contain an intercept")
})

test_that("rank deficient", {
  tmpTestData <- testData
  tmpTestData[,3] <- tmpTestData[,4]
  x <- lm(V1 ~ .,data=tmpTestData)
  expect_error(altR2(x),"At least one indepedent variable is a linear combination of the remaining independent variables")
})

test_that("invalid estimator inputs", {
  expect_error(estimate_adj_R2(NA_real_, 10, 1), "Rsquared must be")
  expect_error(estimate_adj_R2(1.1, 10, 1), "Rsquared must be")
  expect_error(estimate_adj_R2(0.5, 10.5, 1), "N must be")
  expect_error(estimate_adj_R2(0.5, 10, -1), "p must be")
  expect_error(estimate_adj_R2(0.5, 10, 9), "Sample size 10")
})
