test_that("runShiny launches the app correctly", {
  appDir <- system.file("shiny", "app", package = "MuPETFlow")
  # Make sure the app directory exists
  expect_true(dir.exists(appDir))
  # Test that the function runs without errors
  expect_no_error(runMuPETFlow())
})
