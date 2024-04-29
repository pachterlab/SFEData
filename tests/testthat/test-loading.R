test_that("McKellar muscle data loads properly", {
  sfe <- McKellarMuscleData("small")
  expect_s4_class(sfe, "SpatialFeatureExperiment")
  expect_true(all(dim(sfe) > 0))
})

test_that("When it returns a file path", {
    expect_output(fp <- CosMXOutput(file_path = "foo"),
                  "The downloaded files are in .+/foo/cosmx")
    expect_true(grepl("/foo/cosmx$", fp))
    unlink("foo", recursive = TRUE)
})
