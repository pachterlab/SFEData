#' Seurat objects to test coercion function
#'
#' Some of the other test datasets in this package are made as Seurat object to
#' unit test function to convert Seurat objects to SFE.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use. Description of the datasets:
#'\describe{
#' \item{Visium}{From `SeuratData` (ie `stxBrain.SeuratData`), subsetted to keep first 50 genes}
#' \item{VisiumHD}{Visium HD mouse brain data from the \href{https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he}{10X website},
#' with the first 50 genes in 8 um and 16 um bins}
#' \item{Vizgen}{Same dataset in \code{VizgenOutput} with \code{dataset = "hdf5"}}
#' \item{VizgenMulti}{Same as in \code{Vizgen} but with a subset of the data used
#' as if it's another sample to test the coercion function for multiple samples}
#' \item{Xenium}{Same as in \code{XeniumOutput}}
#' \item{XeniumMulti}{Same as in \code{Xenium} but with a subset of the data used
#' as if it's another sample to test the coercion function for multiple samples}
#' }
#' @return A Seurat object
#' @export
SeuratTestData <- .make_data_fun(datasets = c("Visium", "VisiumHD", "Vizgen", "VizgenMulti",
                                              "Xenium", "XeniumMulti"),
                                 ids = 0) # Placeholder
