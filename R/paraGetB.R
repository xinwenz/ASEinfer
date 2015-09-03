paraGetB <- function(x) {
    expAll <- list(dCo= x[(nRepHy+1):nRep * 3],
                   dHy=x[(dvd+1):nRepHy *3 ],
                   xCo=x[nRepHy:(nRep-1) * 3 +1 ],
                   nCo=x[nRepHy:(nRep-1) * 3 +1 ] + x[nRepHy:(nRep-1) * 3 +2 ],
                   xHy=x[dvd:(nRepHy-1)*3 +1],
                   nHy=x[dvd:(nRepHy-1)*3 +1] + x[dvd:(nRepHy-1)*3 +2])


    expCis <- list(dHy=x[ (1:dvd)*3  ],
                   xHy=x[ (0:(dvd-1))*3 +1 ],
                   nHy=x[ (0:(dvd-1))*3 +1 ] + x[ (0:(dvd-1))*3 +2 ])

    fit.binom.cis <-  mle2(
        minuslogl = neglhBinomial_cisonly,
        optimizer = "nlminb",
        lower     = list(log.ec = -20),
        start     = list(log.ec =  0),
        upper     = list(log.ec =  20),
        data      = expCis)

    fit.binom.all <- mle2(
        minuslogl = neglhBinomial,
        optimizer = "nlminb",
        lower     = list(log.ec = -20, log.et = -20),
        start     = list(log.ec =   0, log.et =   0),
        upper     = list(log.ec =  20, log.et =  20),
        data      = expAll)


    tmp1 <- coef(fit.binom.cis)

    tmp2 <- tryCatch({
        confint(fit.bimom.cis)
    },error=function(e){
        c(0,0)
    })

    tmp3 <- coef(fit.binom.all)

    tmp4 <- tryCatch({
        confint(fit.binom.all)
    },error=function(e){
        rbind(log.ec=c(0,0),log.et=c(0,0))
    })

    c(tmp1["log.ec"],tmp2[1],tmp2[2],
      tmp3["log.et"],tmp4["log.et",1],tmp4["log.et",2])
}
