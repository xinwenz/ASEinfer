 organize <- function(expDat) {
    ### expDat should be a list of data frames,each one represent one sample
    spc1 <- names(expDat[[1]])[2]
    spc2 <- names(expDat[[1]])[3]
    
    ### choose hybrid,get replicate number,merge to df,get d
    hybridDat <- expDat[sapply(expDat,function(x) attributes(x)$src == "hybrid" )]
    hybridNum <- length(hybridDat)
    hybridDf <- Reduce(function(x) merge(x,by="name",all=T),hybridDat) ## merged dataframe
    hybridRatio <- sapply(hybridDat,function(x) attributes(x)$d)
    
    ###choose coculture,get replicate number,merge to df,get d
    cocultureDat <- expDat[sapply(expDat,function(x) attributes(x)$src == "coculture" )]
    cocultureNum <- length(cocultureDat)
    cocultureDf <- Reduce(function(x) merge(x,by="name",all=T),cocultureDat) ## merged dataframe
    cocultureRatio <- sapply(cocultureDat,function(x) attributes(x)$d)
    
    ### merge hybrid and Cocu
    expDatDf <- merge(hybridDf,cocultureDf,by="name",all=TRUE)
    
    #put exp data for one gene in one object, and put all my object into a list
    allGene <- list()
    for (i in 1:nrow(expDatDf)) {
        ###
        tmpHybrid <- expDatDf[i,2:(2*hybridNum+1)]
        rawHybrid <- as.data.frame(matrix(tmpHybrid,ncol=2,byrow=TRUE))
        names(rawHybrid) <- c(spc1,spc2)
        rawHybrid <- cbind(rawHybrid,hybridRatio)
        rawHybrid <- na.omit(rawHybrid)
        attributes(rawHybrid)$replicate <- nrow(rawHybrid)
        ###
        tmpCocu <- expDatDf[i,-(1:(2*hybridNum+1))]
        rawCocu <- as.data.frame(matrix(tmpCocu,ncol=2,byrow=TRUE))
        names(rawCocu) <- c(spc1,spc2)
        rawCocu <- cbind(rawCocu,cocultureRatio)
        rawCocu <- na.omit(rawCocu)
        attributes(rawCocu)$replicate <- nrow(rawCocu)
        ###
        tmp <- list("rawHybrid"=rawHybrid,"rawCocu"=rawCocu)
        class(tmp) <- "oneGeneRaw"
        assign(as.character(expDatDf[i,1]),tmp)
        allGene[[as.character(expDatDf[i,1])]] <- tmp
        print(i) # see the speed
    }
    return(allGene)
}


### put in a list of data frames with attributes---src <- (hybrid || coculture), and ratio d <- (any number ), 
###give out a list of objects, each objects belongs to class "oneGeneRaw", one object is constitute by two data.frame: one is hybrid raw data and the other is coculture raw data, each dataframe has a special attribute : replicate: show how many samples for hybrid and Coculture. 