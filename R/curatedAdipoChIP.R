#' \code{curatedAdipoChIP} package
#'
#' A Curated ChIP-Seq Dataset of MDI-induced Differentiated Adipocytes (3T3-L1)
#'
#' A curated dataset of publicly available ChIP-sequencing of transcription
#' factors, chromatin remodelers and histone modifications in the 3T3-L1
#' pre-adipocyte cell line. The package document the data collection,
#' pre-processing and processing of the data. In addition to the documentation,
#' the package contains the scripts that was used to generated the data.
#'
#' @docType package
#' @name curatedAdipoChIP
#'
#' @details The dataset can be accessed through the
#' \code{\link{ExperimentHubData}} as a \code{RangedSummarizedExperiment}
#' object contains:
#' \describe{
#' \item{assay}{The read counts \code{matrix}.}
#' \item{colData}{The phenotype data and quality control data of the samples.}
#' \item{rowRanges}{The feature data and annotation of the peaks.}
#' \item{metadata}{The study level metadata which contains one object called
#' \code{studies}. This is a \code{data.frame} of bibliography information of
#' the studies from which the samples were collected.}
#' }
#'
#' @examples
#' \dontrun{
#' # load the data object
#' library(ExperimentHub)
#' hub <- ExperimentHub()
#' x <- query(hub, c("curatedAdipoChIP", "peak_counts"))
#' x
#'
#' ## download resource
#' peak_counts = x[[1]]
#' }
NULL
