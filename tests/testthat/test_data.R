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

test_that("test data generators are reproducible and preserve timelines", {
  set.seed(123)
  first <- test_data_1()
  set.seed(123)
  second <- test_data_1()

  expect_equal(first, second)

  set.seed(321)
  shifted <- test_data_1()
  expect_equal(shifted$time, first$time)
  expect_false(identical(shifted$data_1, first$data_1))

  segment_means <- with(test_data_44(), tapply(data_44, rep(1:3, each = 100), mean))
  expect_equal(length(segment_means), 3)
  expect_true(all(segment_means > 7 & segment_means < 33))
})
