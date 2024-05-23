#' @keywords internal
"_PACKAGE"

#' mehenDi: Tree-based differrential testing for RNA-Seq
#'
#' This package provides a tree based differential testing method for RNA-Seq
#' data, where the tree is obtained from \strong{TreeTerminus}. The
#' candidate nodes can consist of both leaves and inner nodes, with
#' inner nodes comprising a subset of transcripts that exhibit high uncertainty.
#'
#' The main function is:
#' \itemize{
#' \item \code{\link{mehenDi}} - performs tree-based differential testing
#' }
#' @import fishpond
#' @import TreeSummarizedExperiment
#' @import beaveR
#' @importFrom methods is as slot
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assayNames<- assay assay<- assays assays<-
#' colData colData<- mcols mcols<-
#' @importFrom matrixStats rowRanks rowMedians rowVars colVars rowQuantiles
#' @importFrom Matrix rowMeans
#' @importFrom svMisc progress
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom phangorn Descendants Ancestors
#' @importFrom parallel mclapply
#' @importFrom stats p.adjust
#' @importFrom abind abind
#' @name mehenDi-package
#' @aliases trenDi-package
NULL
