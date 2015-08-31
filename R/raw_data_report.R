check_model <- function(allGene) {
    numGene <- length(allGene)
    reportLs <- lapply(allGene,function(x) c(attributes(x$rawHybrid)$replicate,attributes(x$rawCocu)$replicate))
    reportDf <- do.call(rbind.data.frame,reportLs)
    names(reportDf) <- c("hybrid_Rp","Cocu_Rp")
    reportDf <- data.frame(reportDf,Model=as.character(Model),highPower,stringsAsFactors=FALSE)
     
    ### conditions unavoid correlated error 
    lh <- which(reportDf$hybrid_Rp < 2 & reportDf$highPower==FALSE)
    numLackHybrid <- length(lh)
    if (numLackHybrid>0) {
        warning(numLackHybrid," gene(s) don't have enough Hybird data input. Correlated error cannot be avoided for the ",numLackHybrid, " gene(s).")
    reportDf[lh,"highPower"] <- TRUE
    }
    
    ### conditions can't use betabinomial model
    lhc <- which((reportDf$hybrid_Rp < 4 | reportDf$Cocu_Rp < 2) & (Model == "both" | Model=="betaBinomial"))
    numLackHC <- length(lhc)
    if (numLackHC > 0) {
        warning(numLackHC," gene(s) don't have enough Hybrid and Coculture input. betaBinomial model can't be used. binomial model will be used for the ",numLackHC ," gene(s) instead.")
    reportDf[lhc,"Model"] <- "binomial" 
    }
        
    ### add model into each gene's Raw data
    ModelLs <- split(reportDf[c("Model","highPower")],rownames(reportDf))
    allGeneModel <- mapply(c,allGene,ModelLs,SIMPLIFY = FALSE)
    return(allGeneModel)
}