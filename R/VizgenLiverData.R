#' Vizgen MERFISH mouse liver data
#'
#' This is one of the example datasets from Vizgen's website, downloaded from
#' \href{https://console.cloud.google.com/storage/browser/vz-liver-showcase/Liver1Slice1;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false&pli=1}{here}.
#' The gene count matrix, cell metadata provided by Vizgen, QC metrics, and cell
#' segmentation in one z-plane are included in the SFE object. While it appears
#' that 7 z-planes are in the cell boundaries in the hdf5 files on the website,
#' all 7 z-planes are the same, so only one was used here.
#'
#' This is the largest SFE example dataset thus far. Whlie the SFE object can
#' fit into memory due to the relatively small number of genes, we are
#' considering making a \code{HDF5Array} version of this example dataset.
#' Furthermore, the geometries of the large number of cells can also consume a
#' lot of memory. We are considering using
#' \href{https://github.com/apache/incubator-sedona}{Apache Sedona} and possibly
#' \href{https://bioconductor.org/packages/release/bioc/html/SQLDataFrame.html}{SQLDataFrame}
#' for on disk geometries and geometric operations in a future version of
#' \code{SpatialFeatureExperiment} and \code{Voyager}.
#'
#' While Vizgen provides transcript spot locations, we don't yet know what to do
#' with the huge dataset, so this is not included in the SFE object.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use, for now can only be "Liver1Slice1".
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
VizgenLiverData <- .make_data_fun(datasets = "Liver1Slice1", ids = 7746)
