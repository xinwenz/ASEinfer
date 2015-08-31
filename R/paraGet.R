require("foreach")
require("doParallel")
require("iterators")

paraGet <- function(readData,pe=2,model="BetaB") {
    # check form of the ReadData data.frame, any one of the following conditions should be stictly true.
    if (nRepHy < 6 ){stop("Replicates of Hybrid experiment is not enough to avoid correlated error")}
    dvd <- floor(nRepHy/2)
    if (nRepCo <3 ) {stop("Replicates of CO-culture less than 3, not enough to get parameters")}
    if (ncol(readData) %% 3 != 0 ) {stop("incorrect form of ReadData")}
    if (sum(grepl("hy",names(readData) != 3*nRepHy))) {stop("hybrid part of ReadData has incorrect columns")}
    if (sum(grepl("co",names(readData) != 3*nRepCo))) {stop("co-culture part of ReadData has incorrect columns")}

    readDataMt <- as.matrix(readData)
    # beta binomial model
    if (model=="BetaB") {
    dbetabinom <- function(k, n, u, rho, log = FALSE) {
        a <- u/rho
        b <- (1-u)/rho
        if (log == TRUE) {
            return(lchoose(n = n, k = k) + lbeta(a = k+a, b = n-k+b) - lbeta(a = a, b = b))
        } else {
            return(choose(n = n, k = k) * (beta(a = k+a, b = n-k+b) /  beta(a = a, b = b)))
        }
    }

    neglhbetaBinomial <- function(log.ec, log.et, log.rCo, log.rHy, dCo,dHy,xCo,nCo,xHy,nHy) {
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

    neglhbetaBinomial_cisonly <- function(log.ec,log.rHy,dHy,xHy,nHy) {
        ec  <- 2^log.ec
        rHy <- 2^log.rHy
        uHy <- (dHy*ec)/(dHy*ec + 1)
        result <- -sum(dbetabinom(k = xHy, n = nHy, u = uHy, rho = rHy, log = TRUE))
        return(result)
    }


    no_cores <- pe  # the number of cores
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)# initiate cluster # not for windows users
    ans <- foreach(d=iter(readDataMt,by="row"),
            .combine=rbind,
            .export=c("nRepCo","nRepHy","nRep","dvd","paraGetBB","dbetabinom","neglhbetaBinomial","neglhbetaBinomial_cisonly"),
            .packages='bbmle') %dopar%
        paraGetBB(d)
    stopCluster(cl)
    row.names(ans) <- row.names(readData)
    colnames(ans) <- c("log.ecis","ecis-","ecis+","log.etrans","etrans-","etans+","log.rHy","rHy-","rhy+","log.rCo","rCo-","rCo+")
    return(ans)
    }





    # binomial model
    if (model=="B") {
    neglhBinomial <- function(log.ec, log.et, dCo,dHy,xCo,nCo,xHy,nHy) {
        ec <- 2^log.ec
        et <- 2^log.et
        pHy <- (dHy*ec)/(dHy*ec + 1)
        pCo <- (dCo*ec*et)/(dCo*ec*et+1)
        result <-
            -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE) +
                     dbinom(x = xCo, size = nCo, prob = pCo, log = TRUE))
        return(result)
    }

    neglhBinomial_cisonly <- function(log.ec,dHy,xHy,nHy){
        ec <- 2^log.ec
        pHy <- (dHy*ec)/(dHy*ec + 1)
        result <- -sum(dbinom(x = xHy, size = nHy, prob = pHy, log = TRUE))
        return(result)
    }

    no_cores <- pe  # the number of cores
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)# initiate cluster # not for windows users
    ans <- foreach(d=iter(readDataMt,by="row"),
            .combine=rbind,
            .export=c("nRepCo","nRepHy","nRep","dvd","paraGetB","neglhBinomial","neglhBinomial_cisonly"),
            .packages='bbmle') %dopar%
        paraGetB(d)
    stopCluster(cl)
    row.names(ans) <- row.names(readData)
    colnames(ans) <- c("log.ecis","ecis-","ecis+","log.etrans","etrans-","etans+")
    return(ans)
    }
}
