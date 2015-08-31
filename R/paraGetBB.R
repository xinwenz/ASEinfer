ParaGetBB <- function(x) {
    expCis <- list(dHy=x[ (1:dvd)*3  ],
                    xHy=x[ (0:(dvd-1))*3 +1 ],
                    nHy=x[ (0:(dvd-1))*3 +1 ] + x[ (0:(dvd-1))*3 +2 ])

    expAll <- list(dCo= x[(nRepHy+1):nRep * 3],
                    dHy=x[(dvd+1):nRepHy *3 ],
                    xCo=x[nRepHy:(nRep-1) * 3 +1 ],
                    nCo=x[nRepHy:(nRep-1) * 3 +1 ] + x[nRepHy:(nRep-1) * 3 +2 ],
                    xHy=x[dvd:(nRepHy-1)*3 +1],
                    nHy=x[dvd:(nRepHy-1)*3 +1] + x[dvd:(nRepHy-1)*3 +2])

    fit.beta.cis <- mle2(
        minuslogl = neglhbetaBinomial_cisonly,
        optimizer = "nlminb",
        lower     = list(log.ec = -20,log.rHy = -20),
        start     = list(log.ec =   0,log.rHy =   0),
        upper     = list(log.ec =  20,log.rHy =  20),
        data      = expCis)

    fit.beta.all <- mle2(
        minuslogl = neglhbetaBinomial,
        optimizer = "nlminb",
        lower     = list(log.ec = -20, log.et = -20, log.rHy = -20, log.rCo = -20),
        start     = list(log.ec =   0, log.et =   0, log.rHy =   0, log.rCo =   0),
        upper     = list(log.ec =  20, log.et =  20, log.rHy =  20, log.rCo =  20),
        data      = expAll)

    tmp1 <- coef(fit.beta.cis)
    tmp2 <- confint(fit.beta.cis)
    tmp3 <- coef(fit.beta.all)
    tmp4 <- confint(fit.beta.all)

    c(tmp1["log.ec"],tmp2["log.ec",1],tmp2["log.ec",2],
      tmp3["log.et"],tmp4["log.et",1],tmp4["log.et",2],
      tmp1["log.rHy"],tmp2["log.rHy",1],tmp2["log.rHy",2]
      tmp3["log.rCo"],tmp4["log.rCo",1],tmp4["log.rCo",2])
}
