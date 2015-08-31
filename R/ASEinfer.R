ASEinfer <- function(hyList,dHy,coList,dCo,pe=3,mdl="BetaB") {
    require("foreach")
    require("doParallel")
    require("iterators")
    df <- formChange(hyList,dHy,coList,dCo)
    print(dfdf)
    ans <- paraGet(dfdf,paraEnv=pe,model=mdl)
    return(ans)
}
