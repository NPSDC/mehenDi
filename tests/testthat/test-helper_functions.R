test_that("data", {
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
  expect_error(computeSign(ts, "condition"))

  ts <- beaveR::scInfReps(ts)
  expect_error(computeSign(ts, 12))
  expect_error(computeSign(ts, "candition"))
  signs = computeSign(ts, "condition")
  expect_equal(sum(signs==0), 3933)
  expect_equal(sum(signs==-1), 3192)
  expect_equal(sum(signs==1), 1817)
})