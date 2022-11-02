#' Nanostring FFPE CosMX human NSCLC data
#'
#' One of the CosMX example datasets for human non small cell lung cancer
#' (NSCLC, Lung5_Rep1) from the
#' \href{https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/}{Nanostring
#' website} was downloaded and formatted into an SFE object. The dataset is
#' described in the paper High-plex Multiomic Analysis in FFPE at Subcellular
#' Level by Spatial Molecular Imaging,
#' \href{https://doi.org/10.1101/2021.11.03.467020}{He et al}. Since there's no
#' easy way to get the cell segmentation polygon coordinates from the Nanostring
#' website, the polygon coordinates were downloaded from Seurat's vignette. The
#' raw count matrix, QC metrics, cell segmentation in one z-plane, and other
#' cell attributes such as area, aspect ratio, mean DAPI level, mean
#' immunofluorescence signal, and etc. are included in the SFE object.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use, for now can only be "Lung5_Rep1".
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
HeNSCLCData <- .make_data_fun(datasets = "Lung5_Rep1", ids = 7743)
