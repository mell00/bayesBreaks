set.seed(123)

test_that("test data generators return well-formed data frames", {
  d1 <- test_data_1()
  expect_s3_class(d1, "data.frame")
  expect_equal(names(d1), c("time", "data_1"))
  expect_equal(nrow(d1), 90)
  expect_true(all(diff(d1$time) == 1))

  d44 <- test_data_44()
  expect_s3_class(d44, "data.frame")
  expect_equal(names(d44), c("time", "data_44"))
  expect_equal(nrow(d44), 300)
  expect_true(all(diff(d44$time) == 1))
})
