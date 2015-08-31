paraRead <- function(oneGene,Method,Parameters,Confint,level) UseMethod(generic = "paraRead", object = oneGene)

### correlated error result
paraRead.HighPower <- function(oneGene,Method=TRUE,Parameters=TRUE,Confint=TRUE,level=0.95) {
    method_summary <- NULL
    para_summary <- NULL
    confint_summary <- NULL
    
    ### Method
    if (Method==TRUE) {
        method_summary <- data.frame(Model=class(oneGene)[1],Co_Error="Exist",Repeat_Hybrid=attributes(oneGene)$hyRep,Repeat_Coculture=attributes(oneGene)$coRep)
    }

    ### Parameters
    if (Parameters==TRUE) {
        Binomial <- data.frame()
        BetaBinomial <- data.frame()
        if (class(oneGene)[1] == "Binomial") { 
            Binomial <- as.data.frame(c(coef(oneGene$Parameter),MaxLogLik=logLik(oneGene$Parameter)))
            colnames(Binomial) <- "Binomial"
            para_summary <- Binomial
        }
    
        if (class(oneGene)[1] == "BetaBinomial") {
            BetaBinomial <- as.data.frame(c(coef(oneGene$Parameter),MaxLogLik=logLik(oneGene$Parameter)))
            colnames(BetaBinomial) <- "BetaBinomial"
            para_summary <- BetaBinomial
        }
    
    	if (class(oneGene)[1] == "Both" ) { 
            Binomial <- as.data.frame(c(coef(oneGene$Parameter1),MaxLogLik=logLik(oneGene$Parameter1)))
			colnames(Binomial) <- "Binomial"
            BetaBinomial <- as.data.frame(c(coef(oneGene$Parameter2),MaxLogLik=logLik(oneGene$Parameter2)))
            colnames(BetaBinomial) <- "BetaBinomial"
            para_summary <- transform(merge(Binomial,BetaBinomial,by=0,all=TRUE,sort=FALSE),row.names=Row.names,Row.names=NULL)
        }
    }
    
    ### confidence_interval
    if (Confint == TRUE) {
        Binomial <- data.frame()
        BetaBinomial <- data.frame()
        if (class(oneGene)[1] == "Binomial") { 
            Binomial <- as.data.frame(confint(oneGene$Parameter,level=level))   
            colnames(Binomial) <- paste(sep="","Binomial_",colnames(Binomial))
            confint_summary <- Binomial
        }
    
        if (class(oneGene)[1] == "BetaBinomial") {
            BetaBinomial <-as.data.frame(confint(oneGene$Parameter,level=level)) 
            colnames(BetaBinomial) <- paste(sep="","BetaBino_",colnames(BetaBinomial))
            confint_summary <- BetaBinomial
        }

        if (class(oneGene)[1] == "Both" ) { 
            Binomial <- as.data.frame(confint(oneGene$Parameter1,level=level))   
            colnames(Binomial) <- paste(sep="","Binomial_",colnames(Binomial))
            BetaBinomial <-as.data.frame(confint(oneGene$Parameter2,level=level)) 
            colnames(BetaBinomial) <- paste(sep="","BetaBino_",colnames(BetaBinomial))
            confint_summary <- transform(merge(Binomial,BetaBinomial,by=0,all=TRUE,sort=FALSE),row.names=Row.names,Row.names=NULL)        }
    }
    
    ### return something
    return(list(Method=method_summary,Parameters=para_summary,Confint=confint_summary))
}    


### low power part but without correlated error
paraRead.LowPower <- function(oneGene,Method=TRUE,Parameters=TRUE,Confint=TRUE,level=0.95) {
    method_summary <- NULL
    para_summary <- NULL
    confint_summary <- NULL
    
    ### Method
    if (Method==TRUE) {
        method_summary <- data.frame(Model=class(oneGene)[1],Co_Error="No",Repeat_Hy_Cis=attributes(oneGene)$hyRep.cis,Repeat_Hy_Trans=attributes(oneGene)$hyRep.all,Repeat_Co_Trans = attributes(oneGene)$coRep)
    }
    
    ### Parameters
    if (Parameters==TRUE) {
        Binomial <- data.frame()
        BetaBinomial <- data.frame()
        if (class(oneGene)[1] == "Binomial") { 
            Binomial <- as.data.frame(c(coef(oneGene$Parameter.cis),MaxLogLik.cis=logLik(oneGene$Parameter.cis),coef(oneGene$Parameter.all)["log.et"],MaxLogLik=logLik(oneGene$Parameter.all)),optional=TRUE)
            colnames(Binomial) <- "Binomial"
            para_summary <- Binomial

        }
        
        if (class(oneGene)[1] == "BetaBinomial") {
            BetaBinomial <- as.data.frame(c(coef(oneGene$Parameter.cis),MaxLogLik.cis=logLik(oneGene$Parameter.cis),coef(oneGene$Parameter.all)[c("log.et","log.rCo")],MaxLogLik=logLik(oneGene$Parameter.all)),optional=TRUE)
            colnames(BetaBinomial) <- "BetaBinomial"
            para_summary <- BetaBinomial    
        }



        if (class(oneGene)[1] == "Both" ) { 
            Binomial <- as.data.frame(c(coef(oneGene$Parameter1.cis),MaxLogLik.cis=logLik(oneGene$Parameter1.cis),coef(oneGene$Parameter1.all)["log.et"],MaxLogLik=logLik(oneGene$Parameter1.all)),optional=TRUE)
            colnames(Binomial) <- "Binomial"

            BetaBinomial <- as.data.frame(c(coef(oneGene$Parameter2.cis),MaxLogLik.cis=logLik(oneGene$Parameter2.cis),coef(oneGene$Parameter2.all)[c("log.et","log.rCo")],MaxLogLik=logLik(oneGene$Parameter2.all)),optional=TRUE)
            colnames(BetaBinomial) <- "BetaBinomial"  
            para_summary <- transform(merge(Binomial,BetaBinomial,by=0,all=TRUE,sort=FALSE),row.names=Row.names,Row.names=NULL)
        }
    }
    
    ### confidence_interval
    if (Confint == TRUE) {
        Binomial <- data.frame()
        BetaBinomial <- data.frame()
        if (class(oneGene)[1] == "Binomial") { 
            Binomial <- as.data.frame(rbind(confint(oneGene$Parameter.cis,level=level),confint(oneGene$Parameter.all,level=level)["log.et",]))  
            colnames(Binomial) <- paste(sep="","Binomial_",colnames(Binomial))
            rownames(Binomial) <- c("log.ec","log.et")
            confint_summary <- Binomial
        }
        
        if (class(oneGene)[1] == "BetaBinomial") {
            BetaBinomial <-as.data.frame(rbind(confint(oneGene$Parameter.cis,level=level),confint(oneGene$Parameter.all,level=level)[c("log.et","log.rCo"),]))   
            names(BetaBinomial) <- paste(sep="","BetaBino_",colnames(BetaBinomial))
            confint_summary <- BetaBinomial
        }


        if (class(oneGene)[1] == "Both" ) { 
        	Binomial <- as.data.frame(rbind(confint(oneGene$Parameter1.cis,level=level),confint(oneGene$Parameter1.all,level=level)["log.et",]))  
            colnames(Binomial) <- paste(sep="","Binomial_",colnames(Binomial))
            rownames(Binomial) <- c("log.ec","log.et")

            BetaBinomial <-as.data.frame(rbind(confint(oneGene$Parameter2.cis,level=level),confint(oneGene$Parameter2.all,level=level)[c("log.et","log.rCo"),]))   
            names(BetaBinomial) <- paste(sep="","BetaBino_",colnames(BetaBinomial))

            confint_summary <- transform(merge(Binomial,BetaBinomial,by=0,all=TRUE,sort=FALSE),row.names=Row.names,Row.names=NULL)        
        }       
    }
    
    ### return something
    return(list(Method=method_summary,Parameters=para_summary,Confint=confint_summary))
} 