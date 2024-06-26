% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XeniumOutput.R
\name{XeniumOutput}
\alias{XeniumOutput}
\title{Xenium output}
\usage{
XeniumOutput(
  dataset = datasets,
  file_path = ".",
  force = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{dataset}{Either "v1" or "v2", as described above.}

\item{file_path}{Path to save downloaded files for the \code{*Output}
functions which don't return an SFE object.}

\item{force}{Logical, whether to force redownload if the files are already
present. Defaults to \code{FALSE}.}

\item{verbose}{Whether to display progress of download.}
}
\value{
Path to the tarball containing the output directory.
}
\description{
These are small subsets of Xenium output from data downloaded from the 10X
website used to test \code{readXenium()} in SFE. The first subset comes from
the
\href{https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip}{mouse
brain}, generated with Xenium Onboarding Analysis (XOA) v1. The second subset
comes from the
\href{https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Pancreas_FFPE/Xenium_V1_human_Pancreas_FFPE_outs.zip}{human
pancreas}, generated with XOA v2. XOA v1 and v2 have different output
formats, hence the separate subsets used for testing. To reduce download
time, the zarr files in the original output were removed since they are not
used by the SFE package. Also, only lower resolution images are kept.
}
\details{
For the v1 output, the cell and nucleus boundary parquet files have two
versions, one with the arrow raw bytes and the other without. There are no
longer raw bytes since XOA v1.4.

To make the test data smaller, the transcript spots have been down sampled
in the v1 data but not in v2.
}
