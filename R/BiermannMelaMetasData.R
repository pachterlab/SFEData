.make_data_fun <- function(datasets, ids) {
    function(dataset = datasets,
             force = FALSE, verbose = TRUE) {
        dataset <- match.arg(dataset)
        ids <- ids
        names(ids) <- datasets
        id <- paste0("EH", ids[dataset])
        eh <- ExperimentHub()
        ds <- query(eh, "SFEData")
        ds[[id, force = force, verbose = verbose]]
    }
}

#' Melanoma brain metastasis slide-seq2 data
#' 
#' This is one of the human melanoma brain metastasis samples from the paper
#' Dissecting the treatment-naive ecosystem of human melanoma brain metastasis,
#' \href{https://doi.org/10.1016/j.cell.2022.06.007}{Biermann et al}. The
#' dataset is \code{GSM6025935_MBM05_rep1}. The raw counts and cell metadata
#' were downloaded from GEO. The raw counts, QC metrics such as number of UMIs
#' and genes detected per barcode, and centroid coordinates as \code{sf} 
#' POINT geometry, are included in the SFE object.
#' 
#' @inheritParams McKellarMuscleData
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
#' @examples
#' sfe <- BiermannMelaMetasData()
BiermannMelaMetasData <- .make_data_fun(datasets = "MBM05_rep1", ids = 0)
# placehold for id
