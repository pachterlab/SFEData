#' Xenium FFPE human breast cancer data
#'
#' This dataset was downloaded from the
#' \href{https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast}{10X
#' website}, and described in the paper High resolution mapping of the breast
#' cancer tumor microenvironment using integrated single cell, spatial and in
#' situ analysis of FFPE tissue,
#' \href{https://doi.org/10.1101/2022.10.06.510405}{Janesick et al}. The dataset
#' might not be representative of later Xenium data. There are two samples,
#' which can both be downloaded with this package. For each sample, the raw gene
#' counts, QC metrics, cell and nuclei segmentation polygons in one z-plane, and
#' cell centroids are included in the SFE object. The two samples are in
#' separate SFE objects. A small number of nuclei polygons are invalid due to
#' self-intersection; these cases were resolved by making a buffer of distance 0
#' and then removing the holes. Additional cell metadata provided by 10X, such
#' as cell area, are also included.
#'
#' As the SFE and Voyager packages are in the experimental stage and they were
#' originally developed and tested on relatively small Visium datasets, they are
#' not yet very scalable to larger smFISH datasets. While 10X provided
#' transcript spot locations, these are not included in the SFE objects for now
#' as we are not sure if \code{spatstat} can work with such a large dataset for
#' spatial point process analyses, nor does SFE integrate with \code{spatstat}.
#' In a future version of this package, the transcript locations might be added
#' as a separate dataset, but this is not guaranteed.
#'
#' @inheritParams McKellarMuscleData
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @param dataset Which dataset to use, must be one of "rep1" and "rep2".
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
JanesickBreastData <- function(dataset = c("rep1", "rep2"), force = FALSE,
                               verbose = TRUE) {
    dataset <- match.arg(dataset)
    bfc <- BiocFileCache()
    url <- if (dataset == "rep1") 
        "https://caltech.box.com/public/static/9d9s1ixs379q16m6zunpxpeuw1u2o5j5"
    else
        "https://caltech.box.com/public/static/y9f16wiaz70ojvteavctcu3k9zs6ctfl"
    fn <- bfcrpath(bfc, url)
    readRDS(fn)
}
