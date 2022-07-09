test_that("McKellar muscle data loads properly", {
  sfe <- McKellarMuscleData("small")
  expect_s4_class(sfe, "SpatialFeatureExperiment")
  expect_true(all(dim(sfe) > 0))
})
