estimatePThresh <- function(pvalues, adjPval = 0.05) {
  adPval <- which.min(abs(p.adjust(pvalues, method = "BH") - adjPval))
  pvalues[adPval]
    
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
computeSign <- function(y, x, minP = 0.8, pc = 5) {
  stopifnot(is(y, "SummarizedExperiment"))
  stopifnot("infRepsScaled" %in% names(metadata(y)))
  stopifnot(metadata(y)[["infRepsScaled"]])
  stopifnot(is(x, "character"))
  stopifnot(x %in% colnames(colData(y)))

  infRepsArray <- fishpond:::getInfReps(y)
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
