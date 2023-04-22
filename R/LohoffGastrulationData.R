#' seqFISH mouse gastrulation dataset
#'
#' This dataset was downloaded from the
#' \href{https://marionilab.cruk.cam.ac.uk/SpatialMouseAtlas/}{companion website}
#' titled, Integration of spatial and single-cell transcriptomic data elucidates
#' mouse organogenesis \href{https://doi.org/10.1038/s41587-021-01006-2}{Lohoff et al}.

#' There are three biological replicates available in this dataset, each
#' representing a different embryo. For each dataset, the raw gene counts,
#' metadata, and cell segmentation in one z-plane are provided in the SFE object.
#' Segmentation data were not provided for provided for all cells in the count
#' matrix for embryos 1 and 2. In these cases, the segmentation data are
#' represented by empty polygons.
#'
#' While the authors provided the spot location for each mRNA molecule, these
#' are not included in the SFE objects.
#'
#' @inheritParams McKellarMuscleData
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @param dataset Which dataset to use, must be one of "rep1", "rep2", and "rep3".
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
LohoffGastrulationData <- .make_data_fun(datasets = c("rep1", "rep2", "rep3"),
                                         ids = 7800:7802)

