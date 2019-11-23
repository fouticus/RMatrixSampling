#' RMatrixSampling
#'
#' Provides Various Matrix Sampling Routines
#'
#' @docType package
#' @author Alex Fout
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib RMatrixSampling
#' @name RMatrixSampling
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("RMatrixSampling", libpath)
}
