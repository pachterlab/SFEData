.XeniumOutput <- .make_data_fun(datasets = c("v1", "v2"), ids = c(9481,0),
                                fn = c("xenium_lr", "xenium2"))
# 0 for placeholder before I upload the data
.Xenium_get_file <- function(file_path, url) {
    bfc <- BiocFileCache()
    fn <- bfcrpath(bfc, url)
    untar(fn, exdir = file_path)
    out <- file.path(file_path, "xenium2") |> normalizePath()
    cat("The downloaded files are in", out, "\n")
    out
}

#' Xenium output
#'
#' These are small subsets of Xenium output from data downloaded from the 10X
#' website used to test \code{readXenium()} in SFE. The first subset comes from
#' the
#' \href{https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip}{mouse
#' brain}, generated with Xenium Onboarding Analysis (XOA) v1. The second subset
#' comes from the
#' \href{https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Pancreas_FFPE/Xenium_V1_human_Pancreas_FFPE_outs.zip}{human
#' pancreas}, generated with XOA v2. XOA v1 and v2 have different output
#' formats, hence the separate subsets used for testing. To reduce download
#' time, the zarr files in the original output were removed since they are not
#' used by the SFE package. Also, only lower resolution images are kept.
#'
#' For the v1 output, the cell and nucleus boundary parquet files have two
#' versions, one with the arrow raw bytes and the other without. There are no
#' longer raw bytes since XOA v1.4.
#'
#' To make the test data smaller, the transcript spots have been down sampled
#' in the v1 data but not in v2.
#'
#' @inheritParams McKellarMuscleData
#' @param dataset Either "v1" or "v2", as described above.
#' @return Path to the tarball containing the output directory.
#' @export
XeniumOutput <- function(dataset = c("v1", "v2"), file_path = ".", force = FALSE,
                         verbose = TRUE) {
    dataset <- match.arg(dataset)
    switch(dataset,
           v1 = .XeniumOutput("v1", file_path = file_path),
           v2 = .Xenium_get_file(file_path,
                                 "https://caltech.box.com/shared/static/efqxgiqc89tum4ms6jbpulwvqje782uc"))

}