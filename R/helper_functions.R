estimatePThresh <- function(pvalues, adjPval = 0.05) {
  adPval <- which.min(abs(p.adjust(pvalues, method = "BH") - adjPval))
  pvalues[adPval]

}

# The code for this function has been taken from the package fishpond, where infRepError is not an exported function
infRepError <- function(infRepIdx) {
  if (length(infRepIdx) == 0) {
    stop("there are no inferential replicates in the assays of 'y';
see Quick Start in the swish vignette for details")
  }
}

# The code for this function has been taken from the package fishpond, where getInfReps is not an exported function
getInfReps <- function(ys) {
  infRepIdx <- grep("infRep",assayNames(ys))
  infRepError(infRepIdx)
  infReps <- assays(ys)[infRepIdx]
  abind(as.list(infReps), along=3)
}

#' Assigns a sign change between directions to a transcript/node
#'
#' @param y SummarizedExperiment containing the scaled inferential replicates
#' @param x character which is the column in colData
#' that contains the condition information
#' @param minP float value b/w 0 and 1, denoting how much proportion
#' of the inferential replicates should hold the same sign change
#' @param pc numeric, pseudocount, default is 5
#'
#' @return numeric vector, of length equal to the number of rows in y, with
#' 1 indicating positive sign change, -1 negative sign change and 0 unsure
#'

#' @export
computeSign <- function(y, x, minP = 0.70, pc = 5) {
  stopifnot(is(y, "SummarizedExperiment"))
  stopifnot("infRepsScaled" %in% names(metadata(y)))
  stopifnot(metadata(y)[["infRepsScaled"]])
  stopifnot(is(x, "character"))
  stopifnot(x %in% colnames(colData(y)))

  infRepsArray <- getInfReps(y)
  condition <- colData(y)[[x]]
  stopifnot(is.factor(condition))
  stopifnot(nlevels(condition) == 2)
  stopifnot(!anyNA(condition))

  dims <- dim(infRepsArray)
  cond1 <- condition == levels(condition)[1]
  cond2 <- condition == levels(condition)[2]
  diffs <- matrix(nrow=dims[1],ncol=dims[3])
  for (k in seq_len(dims[3])) {
    diffs[,k] <- log2(rowMeans(infRepsArray[,cond2,k]) + pc) -
      log2(rowMeans(infRepsArray[,cond1,k]) + pc)
  }
  signs <- rep(0, dims[1])
  pos <- diffs>0
  #print(head(pos))
  signs[rowMeans(pos) >= minP] = 1
  signs[rowMeans(!pos) >= minP] = -1
  return(signs)
}

defunct = function(msg = "This function is depreciated") function(...) return(stop(msg))
#' Function is deprecated in favor of mehenDi
#' @inheritParams mehenDi
#' @export
trenDi <- function(tse, x, pvalues, minP=0.70, mIRVThresh=0.4, alpha = 0.01, pThresh = NULL, cores=1) {
  defunct("Name of trenDi has been changed to mehenDi, call mehenDi")
}
