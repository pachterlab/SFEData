#' Download McKellar et al. mouse skeletal muscle data
#'
#' In the first version of this package, only the first time point, 2 days after
#' notexin injury, is available. We may add the later time points in later
#' versions of this package.
#'
#' All datasets are \code{SpatialFeatureExperiment} (SFE) objects, with a
#' \code{counts} assay for the raw gene counts. Column metadata includes total
#' UMI counts (\code{nCounts}) and number of genes (\code{nGenes}) detected per
#' spot. Row metadata includes means, variances, and CV2 of each gene in the
#' full dataset. Column geometry includes Visium spot polygons
#' (\code{spotPoly}). Annotation geometry includes tissue boundary
#' (\code{tissueBoundary}), myofiber segmentation (full resolution
#' \code{myofiber_full} and simplified \code{myofiber_simplified}), nuclei
#' segmentation (\code{nuclei}), and nuclei centroids (\code{nuclei_centroid}).
#'
#' Myofibers were segmented manually with the LabKit ImageJ plugin on a 4x
#' downsized H&E image, downsized so the image can be loaded into LabKit, and
#' the \code{terra} R package was used to convert the TIFF segmentation masks
#' into polygons. Coordinates in the SFE objects are in pixels in the full
#' resolution H&E image. Hence the coordinates of the myofiber segmentations
#' were scaled up to match the other coordinates. The full resolution myofiber
#' segmentation looks pixelated; the \code{mapshaper} R package was used to
#' simplify polygons while conserving contiguity. Morphological (area,
#' perimeter, eccentricity, angle) and Haralick (see
#' \code{EBImage::computeFeatures.haralick}) metrics were computed for the
#' myofibers with the \code{EBImage} R package.
#'
#' Nuclei were segmentated with StarDist. About 3000 nuclei from randomly
#' selected regions in the H&E image from this and later time points were
#' manually annotated with LabKit to train the StarDist model, which was then
#' used to segment all nuclei. OpenCV was used to convert segmentation masks
#' into polygons and compute morphological metrics.
#'
#' Tissue boundary was obtained by first thresholding the H&E image by grascale
#' intensity and then converting the mask into polygons with OpenCV. Small
#' pieces which are debris were removed.
#'
#' @param dataset Which dataset to use. Whether the full dataset, the first
#' small subset, or the second small subset. The second small subset has a
#' different \code{sample_id}.
#' @param force Logical, whether to force redownload if the files are already
#' present. Defaults to \code{FALSE}.
#' @param verbose Whether to display progress of download.
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
#' @examples
#' sfe <- McKellarMuscleData("small")
McKellarMuscleData <- .make_data_fun(datasets = c("full", "small", "small2"),
                                     ids = 7560:7562)
