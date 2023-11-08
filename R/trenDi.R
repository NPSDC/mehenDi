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

findCandNodes <- function(tse, nInd, signs, descAll, childNodes, mIRVThresh, pThresh, cores=1) {
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
            findCandNodes(tse, child, signs, descAll, childNodes, mIRVThresh, pThresh, cores = cores)
        }, mc.cores=cores))
        return(cNodes)
    }
}

#' @export
trenDi <- function(tse, x, pvalues, minP=0.70, mIRVThresh=0.4, alpha = 0.01, cores=1) {
    stopifnot(is(tse, "TreeSummarizedExperiment"))
    stopifnot("meanInfRV" %in% colnames(mcols(tse)))
    stopifnot(is(pvalues, "numeric"))
    stopifnot(all(pvalues <= 1, na.rm=T) & all(pvalues >= 0, na.rm=T))
    stopifnot(is(minP,"numeric"))
    stopifnot(is(mIRVThresh,"numeric"))
    stopifnot(is(alpha,"numeric"))
    stopifnot(alpha>=0 & alpha <= 1)
    stopifnot(is(x, "character"))
    stopifnot(x %in% colnames(colData(tse)))
    stopifnot(length(pvalues)==nrow(tse))
    mcols(tse)[["pvalue"]] <- pvalues

    tree <- rowTree(tse)
    pThresh <- estimatePThresh(pvalues[1:l], alpha)
    print(pThresh)

    cNodes <- Descendants(tree, l+1, "child") ## Child nodes of root (l+1)
    leaves <- cNodes[cNodes <= l] ## children of root that are leaf
    sigNodes <- leaves[which(pvalues[leaves] <= pThresh)] ## significant nodes with pvalue less than alpha
    innNodes <- setdiff(cNodes, leaves)

    descAll <- Descendants(tree, 1:nrow(tse), "all")
    childNodes <- Descendants(tree, 1:nrow(tse), "child")
    innSigN <- innNodes[which(sapply(descAll[innNodes], function(nodes) { ## branches that contain atleast 1 sig node
        sum(pvalues[nodes] <= pThresh, na.rm=T)>0 
    }))]


    l <- length(tree$tip.label)
    
    signs <- computeSign(tse, x, minP=minP)
    sigs <- unlist(mclapply(innSigN, function(node) {
        findCandNodes(tse, node, signs, descAll, childNodes, mIRVThresh, pThresh, cores)
    }, mc.cores=cores))
        
    sigNodes <- c(sigNodes, sigs)
    sigNodes
}