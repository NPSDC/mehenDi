test_that("test_mehenDi", {
  dir <- system.file("extdata", package="brainSimSmall")
  clustFile <- file.path(dir, "cluster_nwk.txt")
  quantDir <- file.path(dir, "salmon_quants")
  samples <- as.vector(outer(seq(6), c(1, 2), function(x, y) paste(x, y, sep = "_")))
  quantFiles <- file.path(quantDir, samples, "quant.sf")
  coldata <- data.frame(files = quantFiles, names=samples, condition=factor(rep(c(1,2), each=6)))

  clustFile <- file.path(dir, "cluster_nwk.txt")
  tse <- suppressPackageStartupMessages(beaveR::buildTSE(clustFile, coldata))
  tree <- suppressPackageStartupMessages(TreeSummarizedExperiment::rowTree(tse))

  ts <- beaveR::computeSizeFactors(tse)
  ts <- beaveR::scInfReps(ts)
  x <- "condition"
  pvals <- runif(nrow(tse))
  expect_error(mehenDi(pvals, "condition", pvals, 0.70, 0.40, 0.01))
  expect_error(mehenDi(tse, "codition", pvals, 0.70, 0.40, 0.01))
  expect_error(mehenDi(tse, "condition", pvals+100, 0.40, 0.40, 0.01))
  expect_error(mehenDi(tse, "condition", pvals, 0.70, 0.40, 1.2))
  expect_error(mehenDi(tse, "condition", pvals[1:100], 0.70, 0.40, 0.01))
  expect_error(mehenDi(tse, "condition", pvals, 1.2, 0.40, 0.01))
  expect_error(mehenDi(tse, "condition", pvals, "0.2", 0.40, 0.01))

  ts <- fishpond::labelKeep(tse)
  metadata(ts)[["infRepsScaled"]] <- TRUE
  tC <- mehenDi(ts, "condition", pvals, 0.7, 0.40, 0.01)
  expect_equal(length(tC), 2)
})
