## Checking the three base criteria at each node for it to be called selected
meetCriteria <- function(tse, nInd, signs, descAll, childNodes, mIRVThresh, pThresh) {
    tree <- rowTree(tse)
    pvalues <- mcols(tse)[["pvalue"]]
    mIRV <- mcols(tse)[["meanInfRV"]]

    if(pvalues[nInd] > pThresh | is.na(pvalues[nInd])) ## Check pvalue
        return(F)
    if(nInd <= length(tree$tip)) {  ## leaf node that is significant
       if(is.na(pvalues[nInd]))
           return(F)
       if(pvalues[nInd] <= pThresh)
           return(T)
    }
    desc <- descAll[[nInd]]
    children <- childNodes[[nInd]]

    if((all(signs[desc] >= 0) | all(signs[desc] <= 0)) & all(mIRV[children] > mIRVThresh)) { ##same sign and mIRV
        if(sum(is.na(pvalues[desc])) > 0 | !all(pvalues[desc] <= pThresh, na.rm=T)) {
            if(sum(is.na(pvalues[children])) > 0)
                return(T)
            if(!all(pvalues[children] <= pThresh))  ##all children can't be significant
                return(T)
        }
    }
    return(F)
}

## Return the underlying selected node/s (could be the node itself) for a node
findSelNodes <- function(tse, nInd, signs, descAll, childNodes, mIRVThresh, pThresh, cores=1) {
    tree <- rowTree(tse)
    pvalues <- mcols(tse)[["pvalue"]]
    mIRV <- mcols(tse)[["meanInfRV"]]

    if(meetCriteria(tse, nInd, signs, descAll, childNodes, mIRVThresh, pThresh)) {
        return(nInd)
    }
    else {
        allDescNodes <- descAll[[nInd]]
        if(all(pvalues[allDescNodes] > pThresh, na.rm=T)) {
            return(c())
        }
        cNodes <- c()
        children <- childNodes[[nInd]]
        cNodes <- unlist(mclapply(children, function(child) {
            findSelNodes(tse, child, signs, descAll, childNodes, mIRVThresh, pThresh, cores = cores)
        }, mc.cores=cores))
        return(cNodes)
    }
}

#' mehenDi method: Finding differential nodes in the tree with no two nodes having an ancestor/descendant relationship
#' @param tse TreeSummarizedExperiment which contains the scaled inferential replicates
#' @param x character, the name of the condition variable. A factor for two group analysis
#' @param pvalues numeric, pvalues for all the nodes in the tree
#' @param minP numeric, value betweem 0 and 1, the proportion of the total inferential replicates
#' which should have the same sign change between conditions (default 0.7)
#' @param mIRVThresh numeric, minimum meanInfRV that a node should have for it to be considered for aggregation (default 0.4)
#' @param alpha numeric, the rate for the BH correction on the leaves which will be the pvalue threshold for deeming a node significant, wont be used if pThresh is set to some value
#' @param pThresh numeric, pvalue threshold for climbing, default is NULL since pThresh is computed automatically based on alpha, if not set to NULL then this pvalue threshold will used for deeming a node significant
#' @param cores numeric, the number of cores that will be used during parallelization
#'
#' @return list that contains the selected nodes and the pvalue threshold used for deeeming
#' a node significant
#'
#' @examples
#'
#' dir <- system.file("extdata", package="brainSimSmall")
#' clustFile <- file.path(dir, "cluster_nwk.txt")
#' quantDir <- file.path(dir, "salmon_quants")
#' samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep='_')))
#' quantFiles <- file.path(quantDir, samples, 'quant.sf')
#' coldata <- data.frame(files=quantFiles, names=samples, condition=factor(rep(1:2, each=6)))
#' tse <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata)
#' tree <- TreeSummarizedExperiment::rowTree(tse)
#' l <- length(tree$tip)
#'
#' yAll <- beaveR::computeSizeFactors(tse)
#' yAll <- beaveR::scInfReps(yAll)
#' yAll <- fishpond::labelKeep(yAll)

#' set.seed(1)
#' yTxps <- fishpond::swish(yAll[1:l,], x="condition")
#' yInn <- fishpond::swish(yAll[(l+1):nrow(yAll),], x="condition")
#' pvals <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])
#' tD <- mehenDi(yAll, x="condition", pvalues = pvals,
#'                    minP=0.7, mIRVThresh=0.4, alpha=0.01)
#' @export
mehenDi <- function(tse, x, pvalues, minP=0.70, mIRVThresh=0.4, alpha = 0.01, pThresh = NULL, cores=1) {
    stopifnot(is(tse, "TreeSummarizedExperiment"))
    stopifnot("meanInfRV" %in% colnames(mcols(tse)))
    stopifnot(is(pvalues, "numeric"))
    stopifnot(all(pvalues <= 1, na.rm=T) & all(pvalues >= 0, na.rm=T))
    stopifnot(is(minP,"numeric"))
    stopifnot(minP>=0 & minP <= 1)
    stopifnot(is(mIRVThresh,"numeric"))
    stopifnot(is(alpha,"numeric"))
    stopifnot(alpha>=0 & alpha <= 1)
    stopifnot(is(x, "character"))
    stopifnot(x %in% colnames(colData(tse)))
    stopifnot(length(pvalues)==nrow(tse))
    if(!is.null(pThresh)){
        stopifnot(pThresh >= 0 & pThresh <= 1)
    }
    mcols(tse)[["pvalue"]] <- pvalues

    tree <- rowTree(tse)
    l <- length(tree$tip.label)
    pThresh <- ifelse(is.null(pThresh), estimatePThresh(pvalues[1:l], alpha), pThresh)

    cNodes <- Descendants(tree, l+1, "child") ## Child nodes of root (l+1)
    leaves <- cNodes[cNodes <= l] ## children of root that are leaf
    sigL <- leaves[which(pvalues[leaves] <= pThresh)] ## significant nodes with pvalue less than alpha
    innNodes <- setdiff(cNodes, leaves)

    descAll <- Descendants(tree, 1:nrow(tse), "all")
    childNodes <- Descendants(tree, 1:nrow(tse), "child")
    innSigN <- innNodes[which(sapply(descAll[innNodes], function(nodes) { ## branches that contain atleast 1 sig node
        sum(pvalues[nodes] <= pThresh, na.rm=T)>0
    }))]

    signs <- computeSign(tse, x, minP=minP)
    selNodes <- unlist(mclapply(innSigN, function(node) {
        findSelNodes(tse, node, signs, descAll, childNodes, mIRVThresh, pThresh, cores=1)
    }, mc.cores=cores))

    selNodes <- c(sigL, selNodes)
    return(list("selNodes" = selNodes, "pThresh"=pThresh))
}
