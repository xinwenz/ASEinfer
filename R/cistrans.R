args <- commandArgs(TRUE)
cistrans <- function(expDat="asesim.txt",Model="both",highPower=FALSE) {
    fileDf <- read.table(expDat,header=T)
    geneLstmp <- split(fileDf,fileDf$name)
    geneLs <- lapply(geneLstmp,function(x) split(x,x$expt))


    source("readin.R",local=TRUE)
    geneLsReady <- lapply(geneLs,formChange)
    geneLsOrg <<- lapply(geneLs,getinfo)

    source("raw_data_report.R",local=TRUE)
    allGeneModel <- check_model(geneLsReady)

    source("mle_functions.R",local=TRUE)
    allGenePara <- lapply(allGeneModel,paraGet)

    source("final_result.R",local=TRUE)
    #allGeneFinal <- lapply(allGenePara,function(x) paraRead(x))
    
    allGeneFinal <- vector(mode = "list", length = length(allGenePara))
    names(allGeneFinal) <- names(allGenePara)
    #for (i in 1:3) {
    for (i in 1:length(allGenePara)) {
        out <- tryCatch(
        {
        allGeneFinal[[i]] <- paraRead(allGenePara[[i]])
        },
        error=function(cond){
    	   message(cond)
            message(paste("the problem is in this gene:",i))
        },
        warnings=function(cond) {
    	   message(cond)
    	   message(i)
    	}
    	)
	}
    return(allGeneFinal)
}

ansgenes <- cistrans(expDat=args[1])
posit <- paste(args[1],".RData",sep="")
save(ansgenes,geneLsOrg,file=posit)