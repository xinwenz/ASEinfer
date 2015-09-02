#' allele specific expression inference
#'
#' Based on allel specif expression results, give out allle expression difference due to cis effect or trans effect.
#'
#' @param hyList a list of dataframes which contain hybrid replicates of experiment, each dataframe contains two collums, each represents the expression of one allel. The rownames should be gene name.
#' @param dHy a vector of numbers which the ratio of the two allels' DAN content in Hybrid replicates data. dHy=(allel A genomic amount)/(allel B genomic amount)
#' @param coList same as the hyList, but it's for the co-culture replicates.
#' @param dCo similar to dHy, it's the DNA content ratio for co-culture replicates.
#' @param pe parallel environment, how many threads you want to use. defalt is 3.
#' @param mdl which model you want to use, the default is:"BetaB", which is a beta-binomial model. Another choise is "B", which is a binomial model.
#'
#' @return Return a dataframe of numbers, each row is the result for each gene, and the collumns are the log of effect of cis/trans and their confidence intervals,log of the over-disperse parameter and the corresponding confidence intervals.
#'
#' @examples
#' aseinfer(hyList=list(aspHy01,aspHy02,aspHy03,aspHy04,aspHy05,
#'          aspHy06,aspHy07,aspHy08,aspHy09,aspHy10),
#'          dHy = c(1,1.1,1.5,1,1.23,1,1.42,1.27,1.3,1.23),
#'          coList=list(aspCo17,aspCo18,aspCo19,aspCo20,aspCo21),
#'          dCo=c(1,1.1,1.5,1,1.2),
#'          pe=3,model="BetaB")
#'
#' aseinfer(hyList=list(aspHy01,aspHy02,aspHy03,aspHy04,aspHy05,
#'          aspHy06,aspHy07,aspHy08),
#'          dHy = c(1,1.1,1.5,1,1.23,1,1.42,1.27),
#'          coList=list(aspCo17,aspCo18,aspCo19,aspCo20),
#'          dCo=c(1,1.1,1.5,1),
#'          pe=3,model="B")
#' @export
aseinfer <- function(hyList,dHy,coList,dCo,pe=3,mdl="BetaB") {
    require("foreach")
    require("doParallel")
    require("iterators")

    nRepHy <<- length(hyList) # add to global evn
    nRepCo <<- length(coList)
    nRep <<- nRepHy+nRepCo     # add to global evn
    dvd <<- floor(nRepHy/2)
    nGenes <<- nrow(hyList[[1]])

    if (nRepHy != length(dHy)) {stop("hybrid replicates is not the same number as DNA proportion measurement replicates")}
    if (nRepCo != length(dCo)) {stop("coculture replicates is not the same number as DNA proportion measurement replicates")}
    if (nRepHy < 6 ){stop("Replicates of Hybrid experiment is not enough to avoid correlated error")}
    if (nRepCo <3 ) {stop("Replicates of CO-culture less than 3, not enough to get parameters")}

    df <- formChange(hyList,dHy,coList,dCo)
    # check form of the ReadData data.frame, any one of the following conditions should be stictly true.
    if (ncol(dfdf) %% 3 != 0 ) {stop("incorrect form of ReadData")}
    if (sum(grepl("hy",names(dfdf) != 3*nRepHy))) {stop("hybrid part of ReadData has incorrect number of columns")}
    if (sum(grepl("co",names(dfdf) != 3*nRepCo))) {stop("co-culture part of ReadData has incorrect number columns")}

    readDataMt <- as.matrix(dfdf)

    # beta binomial model
    if (mdl=="BetaB") {
        dbetabinom <<- function(k, n, u, rho, log = FALSE) {
            a <- u/rho
            b <- (1-u)/rho
            if (log == TRUE) {
                return(lchoose(n = n, k = k) + lbeta(a = k+a, b = n-k+b) - lbeta(a = a, b = b))
            } else {
                return(choose(n = n, k = k) * (beta(a = k+a, b = n-k+b) /  beta(a = a, b = b)))
            }
        }

        neglhbetaBinomial <<- function(log.ec, log.et, log.rCo, log.rHy, dCo,dHy,xCo,nCo,xHy,nHy) {
            ec  <- 2^log.ec
            et  <- 2^log.et
            rHy <- 2^log.rHy
            rCo <- 2^log.rCo
            uHy <- (dHy*ec)/(dHy*ec + 1)
            uCo <- (dCo*ec*et)/(dCo*ec*et+1)
            result <-
                -sum(dbetabinom(k = xHy, n = nHy, u = uHy, rho = rHy, log = TRUE)) -
                sum(dbetabinom(k = xCo, n = nCo, u = uCo, rho = rCo, log = TRUE))
            return(result)
        }

        neglhbetaBinomial_cisonly <<- function(log.ec,log.rHy,dHy,xHy,nHy) {
            ec  <- 2^log.ec
            rHy <- 2^log.rHy
            uHy <- (dHy*ec)/(dHy*ec + 1)
            result <- -sum(dbetabinom(k = xHy, n = nHy, u = uHy, rho = rHy, log = TRUE))
            return(result)
        }


        # the number of cores
        cl <- makeCluster(pe)
        registerDoParallel(cl)# initiate cluster # not for windows users
        clusterExport(cl,c("nRepCo","nRepHy","nRep","dvd","dbetabinom","neglhbetaBinomial","neglhbetaBinomial_cisonly"))
        ans <- foreach(d=iter(readDataMt,by="row"),
                       .combine=rbind,
                       .packages='bbmle') %dopar%
            paraGetBB(d)
        stopCluster(cl)
        row.names(ans) <- row.names(dfdf)
        colnames(ans) <- c("log.ecis","ecis-","ecis+","log.etrans","etrans-","etans+","log.rHy","rHy-","rhy+","log.rCo","rCo-","rCo+")
        rm(list=c("nRepCo","nRepHy","nRep","dvd","dbetabinom","neglhbetaBinomial","neglhbetaBinomial_cisonly","nGenes"),envir=.GlobalEnv)
        return(ans)
    }



#############################..Binomial Model..###################

    # binomial model
    if (mdl=="B") {
        neglhBinomial <<- function(log.ec, log.et, dCo,dHy,xCo,nCo,xHy,nHy) {
            ec <- 2^log.ec
            et <- 2^log.et
            pHy <- (dHy*ec)/(dHy*ec + 1)
            pCo <- (dCo*ec*et)/(dCo*ec*et+1)
            result <-
                -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE) +
                         dbinom(x = xCo, size = nCo, prob = pCo, log = TRUE))
            return(result)
        }

        neglhBinomial_cisonly <<- function(log.ec,dHy,xHy,nHy){
            ec <- 2^log.ec
            pHy <- (dHy*ec)/(dHy*ec + 1)
            result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
            return(result)
        }

        # the number of cores
        cl <- makeCluster(pe)
        registerDoParallel(cl)# initiate cluster # not for windows users
        clusterExport(cl,c("nRepHy","nRepCo","nRep","dvd","neglhBinomial","neglhBinomial_cisonly")) # export from gloabl Env
        ans <- foreach(d=iter(readDataMt,by="row"),
                .combine=rbind,
                .verbose=T,
                .packages='bbmle') %dopar%
            paraGetB(d)
        stopCluster(cl)
        row.names(ans) <- row.names(dfdf)
        colnames(ans) <- c("log.ecis","ecis-","ecis+","log.etrans","etrans-","etans+")
        rm(list=c("nRepHy","nRepCo","nRep","dvd","neglhBinomial","neglhBinomial_cisonly","nGenes"),envir=.GlobalEnv)
        return(ans)
    }
}
