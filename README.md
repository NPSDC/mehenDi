# trenDi

## trenDi: Tree-based differential testing using inferential uncertainty for RNA-Seq
The **trenDi** method performs differential testing use the transcripts-tree obtained from **TreeTerminus**.
The transcripts are arranged in a way such that the uncertainty decreases on ascending the tree, with the
inner nodes comprising transcripts groups. **trenDi** aims to leverage the tree structure to find differential
signal in the tree which might get lost due to uncertainty. The candidate nodes can consist of both transcripts
and inner nodes, with inner nodes representing a set of transcripts that exhibit high uncertainty.

### Installation
```
devtools::install_github("NPSDC/trenDi")
```

### Quick start

#### Loading Data
```
dir <- system.file("extdata", package="brainSimSmall")
clustFile <- file.path(dir, "cluster_nwk.txt")
quantDir <- file.path(dir, "salmon_quants")
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep='_')))
quantFiles <- file.path(quantDir, samples, 'quant.sf')
coldata <- data.frame(files=quantFiles, names=samples, condition=factor(rep(1:2, each=6)))

tse <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata) ## could be the path to any Salmon quantified RNA-Seq experiment and the corresponding forest file produced by TreeTerminus
print(tse)
tree <- TreeSummarizedExperiment::rowTree(tse)
print(tree)
l <- length(tree$tip)
```

#### Running trenDi
```
suppressPackageStartupMessages(library(trenDi))
set.seed(1)
yTxps <- tse[1:l,]
yTxps <- fishpond::swish(yTxps, x="condition")

yAll <- beaveR::computeSizeFactors(tse)
yAll <- beaveR::scInfReps(yAll)
yAll <- fishpond::labelKeep(yAll)
metadata(yAll)$infRepsScaled <- TRUE
yInn <- yAll[(l+1):nrow(yAll),]
set.seed(1)
yInn <- fishpond::swish(yInn, x="condition")

pvals <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])
tD <- trenDi::trenDi(yAll, x="condition", pvalues = pvals,
                    minP=0.7, mIRVThresh=0.4, alpha=0.01)
print(length(tD$candNodes)) ## Total number of candidate nodes
print(sum(tD$candNodes > l)) ## Number of inner nodes that are candidates
```
The input to **trenDi** is the *TreeSummarizedExperiment* (tse) object. This object can be created using the package [beaveR](https://github.com/NPSDC/beaver), 
which requires the same input as [tximeta](https://github.com/thelovelab/tximeta) and the forest file produced by [TreeTerminus](https://github.com/COMBINE-lab/TreeTerminus).
Among the other inputs are the p-values corresponding to each node in tree belonging to the `tse` object, which here in this example have been generated using `swish` in the
package [fishpond](https://github.com/thelovelab/fishpond).

#### Exploration of the nodes output by trenDi
```{r}
fishpond::plotInfReps(yAll, tD[["candNodes]][1], x="condition")
nodeInf <- beaveR::findNodeInformation(yAll, node = tD[["candNodes]][1], type = "tips")
print(nodeInf)
```
