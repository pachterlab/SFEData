.make_data_fun <- function(datasets, ids) {
    function(dataset = datasets, file_path = ".",
             force = FALSE, verbose = TRUE) {
        dataset <- match.arg(dataset)
        ids <- ids
        names(ids) <- datasets
        id <- paste0("EH", ids[dataset])
        eh <- ExperimentHub()
        ds <- query(eh, "SFEData")
        out <- ds[[id, force = force, verbose = verbose]]
        if (is.character(out)) {
            untar(out, exdir = file_path)
            # To temporarily deal with my stupid mistake before Lori corrects it
            if (list.dirs(file_path, recursive = FALSE) == file.path(file_path, "inst")) {
                browser()
                path_use <- list.dirs(file.path(file_path, "inst", "extdata"), recursive = FALSE)
                out <- file.path(file_path, basename(path_use))
                file.rename(path_use, out)
                unlink(file.path(file_path, "inst"), recursive = TRUE)
            }
            cat("The downloaded files are in", out, "\n")
        }
        return(out)
    }
}

#' Melanoma metastasis slide-seq2 data
#'
#' This function can download one of the human melanoma brain metastasis (MBM)
#' samples and one of the melanoma extracranial metastasis (ECM) samples from
#' the paper Dissecting the treatment-naive ecosystem of human melanoma brain
#' metastasis, \href{https://doi.org/10.1016/j.cell.2022.06.007}{Biermann et
#' al}. The datasets are \code{GSM6025935_MBM05_rep1} and
#' \code{GSM6025946_ECM01_rep1}. The raw counts and cell metadata were
#' downloaded from GEO. The raw counts, QC metrics such as number of UMIs and
#' genes detected per barcode, and centroid coordinates as \code{sf} POINT
#' geometry, are included in the SFE object.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use, must be one of "MBM05_rep1" and
#' "ECM01_rep1".
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
#' @examples
#' sfe <- BiermannMelaMetasData()
BiermannMelaMetasData <- .make_data_fun(datasets = c("MBM05_rep1", "ECM01_rep1"),
                                        ids = 7741:7742)
