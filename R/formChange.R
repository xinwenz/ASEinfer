formChange <- function(hyList,dHy,coList,dCo) { #hyList=list of dataframe,dHy=vector d
    nRepHy <<- length(hyList) # add to global evn
    nRepCo <<- length(coList)
    nRep <<- nRepHy+nRepCo     # add to global evn

    if (nRepHy != length(dHy)) {stop("hybrid replicates is not the same number as DNA proportion measurement replicates")}
    if (nRepCo != length(dCo)) {stop("coculture replicates is not the same number as DNA proportion measurement replicates")}


    #combine all hybrid replicate with DNA proportion
    dfHy <- hyList[[1]]
    nGenes <<- nrow(dfHy) # add to global evn
    dfHy <- cbind(dfHy,dHy[1])
    colnames(dfHy) <- c("hyA.1","hyB.1","hyd.1")

    for(i in 2:nRepHy) {
        if (nGenes != nrow(hyList[[i]])) {stop("gene numbers is different in Hybird replicates")}
        tmp <- cbind(hyList[[i]],dHy[i])
        colnames(tmp) <- c(paste0("hyA.",i),paste0("hyB.",i),paste0("hyd.",i))
        dfHy <- transform(merge(dfHy,tmp,by=0),row.names=Row.names,Row.names=NULL)
    }

    #combine all cocultrue replicate with DNA proportion
    dfCo <- coList[[1]]
    if(nGenes != nrow(dfCo)) {stop("gene number in co-culture replicates is different from hybrid replicates")}
    dfCo <- cbind(dfCo,dCo[1])
    colnames(dfCo) <- c("coA.1","coB.1","cod.1")

    for(i in 2:nRepCo) {
        if (nGenes != nrow(coList[[i]])) {stop("gene numbers is different in CO-culture replicates")}
        tmp <- cbind(coList[[i]],dCo[i])
        colnames(tmp) <- c(paste0("coA.",i),paste0("coB.",i),paste0("cod.",i))
        dfCo <- transform(merge(dfCo,tmp,by=0),row.names=Row.names,Row.names=NULL)
    }

    return(transform(merge(dfHy,dfCo,by=0),row.names=Row.names,Row.names=NULL))
}


