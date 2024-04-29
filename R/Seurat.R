#' Seurat objects to test coercion function
#'
#' Some of the other test datasets in this package are made as Seurat object to
#' unit test function to convert Seurat objects to SFE.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Which dataset to use. Description of the datasets:
#'\describe{
#' \item{Visium}{From `SeuratData` (ie `stxBrain.SeuratData`), subsetted to keep first 50 genes}
#' \item{VisiumHD8}{Visium HD mouse brain data from the \href{https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he}{10X website},
#' with the first 50 genes in 8 um bins}
#' \item{VisiumHDmulti}{Visium HD mouse brain data from the \href{https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he}{10X website},
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
SeuratTestData <- function(dataset = c("Visium", "VisiumHD8", "VisiumHDmulti",
                                       "Vizgen", "VizgenMulti",
                                       "Xenium", "XeniumMulti"),
                           file_path = ".", force = FALSE, verbose = TRUE) {
    dataset <- match.arg(dataset)
    bfc <- BiocFileCache()
    url <- switch(dataset,
                  Visium = "https://caltech.box.com/shared/static/zni0tn9zttbe5rf88u8im6syl49rdp1w",
                  VisiumHD8 = "https://caltech.box.com/shared/static/xlgsmskh7iqih6d8m83rze4281fr78cq",
                  VisiumHDmulti = "https://caltech.box.com/shared/static/y4vcocisqh55psbbrkcv2pbrvyttyj7s",
                  Vizgen = "https://caltech.box.com/shared/static/fg9v2awggcjkjbw2hxmhm3awa91z1n1p",
                  VizgenMulti = "https://caltech.box.com/shared/static/fn0ip859bbz8xtkjzpl8lz7zecsxe7tr",
                  Xenium = "https://caltech.box.com/shared/static/aphhxmw4ycz5r9sjv64grvm6xw08msub",
                  XeniumMulti = "https://caltech.box.com/shared/static/b0tds2olets0vegix0e72gotes9yk35l"
                  )
    fn <- bfcrpath(bfc, url)
    readRDS(fn)
}
#SeuratTestData <- .make_data_fun(datasets = c("Visium", "VisiumHD8", "VisiumHDmulti",
#                                              "Vizgen", "VizgenMulti",
#                                              "Xenium", "XeniumMulti"),
#                                 ids = 0) # Placeholder
