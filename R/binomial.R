library(bbmle)
paraGet <- function(oneGeneModel) {
    # binomial model
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

    # beta binomial model
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


    #the first easy model
    if (oneGeneModel$highPower == TRUE){
        experiment <- list(dCo=unlist(oneGeneModel$rawCocu$cocultureRatio),
                           dHy=unlist(oneGeneModel$rawHybrid$hybridRatio),
                           xCo=unlist(oneGeneModel$rawCocu[[1]]),
                           nCo=unlist(oneGeneModel$rawCocu[[1]])+unlist(oneGeneModel$rawCocu[[2]]),
                           xHy=unlist(oneGeneModel$rawHybrid[[1]]),
                           nHy=unlist(oneGeneModel$rawHybrid[[1]])+unlist(oneGeneModel$rawHybrid[[2]])
        )

        fit.binom <- mle2(
            minuslogl = neglhBinomial,
            optimizer = "nlminb",
            lower     = list(log.ec = -20, log.et = -20),
            start     = list(log.ec =   0, log.et =   0),
            upper     = list(log.ec =  20, log.et =  20),
            data      = experiment)

        if (oneGeneModel$Model == "binomial"){
            oneGeneResult <- list(Parameter=fit.binom)
            attributes(oneGeneResult)$hyRep <- nrow(oneGeneModel$rawHybrid)
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCocu)
            class(oneGeneResult) <- c("Binomial","HighPower")
            return(oneGeneResult)
        }


        fit.beta <- mle2(
            minuslogl = neglhbetaBinomial,
            optimizer = "nlminb",
            lower     = list(log.ec = -20, log.et = -20, log.rHy = -20, log.rCo = -20),
            start     = list(log.ec =   0, log.et =   0, log.rHy =   0, log.rCo =   0),
            upper     = list(log.ec =  20, log.et =  20, log.rHy =  20, log.rCo =  20),
            data      = experiment)

        if (oneGeneModel$Model == "betaBinomial"){
            oneGeneResult <- list(Parameter=fit.beta)
            class(oneGeneResult) <- c("BetaBinomial","HighPower")
            attributes(oneGeneResult)$hyRep <- nrow(oneGeneModel$rawHybrid)
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCocu)
            return(oneGeneResult)
        }

        if(oneGeneModel$Model=="both"){
            oneGeneResult <- list(Parameter1=fit.binom,Parameter2=fit.beta)
            class(oneGeneResult) <- c("Both","HighPower")
            attributes(oneGeneResult)$hyRep <- nrow(oneGeneModel$rawHybrid)
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCocu)
            return(oneGeneResult)
        }
    }


    ### the The hard one to avoid correlated error
    if (oneGeneModel$highPower == FALSE){
        hybridNum <- nrow(oneGeneModel$rawHybrid)
        split <- floor(hybridNum/2)

        experiment1 <- list(dHy=unlist(oneGeneModel$rawHybrid$hybridRatio)[1:split],
                            xHy=unlist(oneGeneModel$rawHybrid[[1]])[1:split],
                            nHy=unlist(oneGeneModel$rawHybrid[[1]])[1:split]+unlist(oneGeneModel$rawHybrid[[2]])[1:split])

        experiment2 <- list(dCo=unlist(oneGeneModel$rawCocu$cocultureRatio),
                            dHy=unlist(oneGeneModel$rawHybrid$hybridRatio)[(split+1):hybridNum],
                            xCo=unlist(oneGeneModel$rawCocu[[1]]),
                            nCo=unlist(oneGeneModel$rawCocu[[1]])+unlist(oneGeneModel$rawCocu[[2]]),
                            xHy=unlist(oneGeneModel$rawHybrid[[1]])[(split+1):hybridNum],
                            nHy=unlist(oneGeneModel$rawHybrid[[1]])[(split+1):hybridNum]+unlist(oneGeneModel$rawHybrid[[2]])[(split+1):hybridNum] )

        #         print("secod exp")
        #         print(experiment2)

        fit.binom.cis <-  mle2(
            minuslogl = neglhBinomial_cisonly,
            optimizer = "nlminb",
            lower     = list(log.ec = -20),
            start     = list(log.ec =  0),
            upper     = list(log.ec =  20),
            data      = experiment1)

        fit.binom.all <- mle2(
            minuslogl = neglhBinomial,
            optimizer = "nlminb",
            lower     = list(log.ec = -20, log.et = -20),
            start     = list(log.ec =   0, log.et =   0),
            upper     = list(log.ec =  20, log.et =  20),
            data      = experiment2)

        if (oneGeneModel$Model == "binomial"){
            oneGeneResult <- list(Parameter.cis = fit.binom.cis, Parameter.all=fit.binom.all)
            class(oneGeneResult) <- c("Binomial","LowPower")
            attributes(oneGeneResult)$hyRep.cis <- split
            attributes(oneGeneResult)$hyRep.all <- hybridNum - split
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCocu)
            return(oneGeneResult)
        }

        fit.beta.cis <- mle2(
            minuslogl = neglhbetaBinomial_cisonly,
            optimizer = "nlminb",
            lower     = list(log.ec = -20,log.rHy = -20),
            start     = list(log.ec =   0,log.rHy =   0),
            upper     = list(log.ec =  20,log.rHy =  20),
            data      = experiment1)

        fit.beta.all <- mle2(
            minuslogl = neglhbetaBinomial,
            optimizer = "nlminb",
            lower     = list(log.ec = -20, log.et = -20, log.rHy = -20, log.rCo = -20),
            start     = list(log.ec =   0, log.et =   0, log.rHy =   0, log.rCo =   0),
            upper     = list(log.ec =  20, log.et =  20, log.rHy =  20, log.rCo =  20),
            data      = experiment2)

        if (oneGeneModel$Model == "betaBinomial"){
            oneGeneResult <- list(Parameter.cis =fit.beta.cis,Parameter.all = fit.beta.all)
            class(oneGeneResult) <- c("BetaBinomial","LowPower")
            attributes(oneGeneResult)$hyRep.cis <- split
            attributes(oneGeneResult)$hyRep.all <- hybridNum - split
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCo)
            return(oneGeneResult)
        }

        if (oneGeneModel$Model == "both"){
            oneGeneResult <- list(Parameter1.cis =fit.binom.cis,Parameter1.all = fit.binom.all,Parameter2.cis=fit.beta.cis,Parameter2.all=fit.beta.all)
            class(oneGeneResult) <- c("Both","LowPower")
            attributes(oneGeneResult)$hyRep.cis <- split
            attributes(oneGeneResult)$hyRep.all <- hybridNum - split
            attributes(oneGeneResult)$coRep <- nrow(oneGeneModel$rawCocu)
            return(oneGeneResult)
        }
    }
}
