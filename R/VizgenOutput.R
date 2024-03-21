#' Vizgen MERFISH output from human brain cancer
#'
#' This dataset is unpublished, and is a small subset of the original dataset.
#' When downloaded, the files are in the output format of Vizgen, not an SFE
#' object. The purpose of this dataset is to demonstrate and test
#' \code{\link{readVizgen}}.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use, can be "hdf5" where cell segmentations
#' are stored in hdf5 files or "cellpose" where cell segmentations are stored in
#' a parquet file.
#' @return Path to the tarball containing the output directory.
#' @export
VizgenOutput <- .make_data_fun(datasets = c("hdf5", "cellpose"), ids = 0)
