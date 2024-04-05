#' CosMX output from mouse brain
#'
#' This is a small subset of the CosMX output from the mouse quarter brain
#' downloaded from the
#' \href{https://nanostring.com/resources/coronal-hippocampus-and-cortex-basic-data-files/}{Nanostring
#' website}. When downloaded, the files are in the output format of Vizgen, not
#' an SFE object. The purpose of this dataset is to demonstrate and test
#' \code{\link{readVizgen}}.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Can only be "subset1".
#' @return Path to the tarball containing the output directory.
#' @export
CosMXOutput <- .make_data_fun(datasets = "subset1", ids = 9480)
