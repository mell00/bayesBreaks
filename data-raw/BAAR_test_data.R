## code to prepare `BAAR_test_data` dataset

## Fix seed for reproducibility
set.seed(600)

data_0_a   <- test_data_0_a()
data_0_b   <- test_data_0_b()
data_1     <- test_data_1()
data_2     <- test_data_2()
data_3     <- test_data_3()
data_4     <- test_data_4()
data_5     <- test_data_5()
data_6     <- test_data_6()
data_7     <- test_data_7()
data_8     <- test_data_8()
data_9     <- test_data_9()
data_10    <- test_data_10()
data_11    <- test_data_11()
data_44    <- test_data_44()
data_100   <- test_data_100()
data_200   <- test_data_200()
data_300   <- test_data_300()

## Save as .rda files in data/
usethis::use_data(
  data_0_a, data_0_b,
  data_1, data_2, data_3, data_4,
  data_5, data_6, data_7, data_8, data_9,
  data_10, data_11, data_44,
  data_100, data_200, data_300,
  overwrite = TRUE
)
