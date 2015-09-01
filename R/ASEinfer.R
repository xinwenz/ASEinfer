#' Allele Specific Expression data inference
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
#'  @return a dataframe of numbers, each row is the result for each gene, and the collumns are the log of effect of cis/trans and their confidence intervals,log of the over-disperse parameter and the corresponding confidence intervals
#'
#' @examples
#' ASEinfer(hyList=list(aspHy01,aspHy02,aspHy03,aspHy04,aspHy05,aspHy06,aspHy07,aspHy08,aspHy09,aspHy10),
#' dHy = c(1,1.1,1.5,1,1.23,1,1.42,1.27,1.3,1.23),
#' coList=list(aspCo17,aspCo18,aspCo19,aspCo20,aspCo21),
#' dCo=c(1,1.1,1.5,1,1.2),pe=3,model="BetaB")
#'
#' @export
ASEinfer <- function(hyList,dHy,coList,dCo,pe=3,mdl="BetaB") {
    require("foreach")
    require("doParallel")
    require("iterators")
    df <- formChange(hyList,dHy,coList,dCo)
    dfdf <- df[1:10,]
    ans <- paraGet(dfdf,paraEnv=pe,model=mdl)
    return(ans)
}
