confintplot <- function(sim, para)  { # sim is the list of the orignal data in global evir, para the list of events
	nGenes <- length(sim)
	par(pty="s",mfrow=c(2,2))
	res <- list()

	#for log.ec

	logec.sim <- unlist(lapply(sim,function(x) x$logec))	
	logec.para <- unlist(lapply(para,function(x) subset(x$Parameters,grepl("log.ec",rownames(x$Parameters)),grepl("Beta",colnames(x$Parameters)))))
	logec.para.up <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.ec",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[2]))
	logec.para.dow <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.ec",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[1]))

	res$wrongEs_ec <- names(logec.para)[(logec.sim > logec.para.up)| (logec.sim < logec.para.dow)]
	res$correct_ec <- 1-length(res$wrongEs_ec)/nGenes
	res$lm_ec <- lm(logec.para ~ logec.sim)
	
	plot(logec.para ~ logec.sim,col='red',pch=18,xlim=range(-3,3),ylim=range(-3,3),asp=1)
	segments(logec.sim,logec.para.dow,logec.sim,logec.para.up,col='red')
	logec.sim1 <- logec.sim
	points(logec.sim1 ~ logec.sim, col='green',pch=20)
	
	
	# for log.et
	loget.sim <- unlist(lapply(sim,function(x) x$loget))	
	loget.para <- unlist(lapply(para,function(x) subset(x$Parameters,grepl("log.et",rownames(x$Parameters)),grepl("Beta",colnames(x$Parameters)))))
	loget.para.up <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.et",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[2]))
	loget.para.dow <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.et",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[1]))

	res$wrongEs_et <- names(loget.para)[(loget.sim > loget.para.up) | (loget.sim < loget.para.dow)]
	res$correct_et <- 1-length(res$wrongEs_et)/nGenes
	res$lm_et <- lm(loget.para ~ loget.sim)

	plot(loget.para ~ loget.sim,col='red',xlim=range(-3,3),ylim=range(-3,3),pch=18,asp=1)
	segments(loget.sim,loget.para.dow,loget.sim,loget.para.up,col='red')
	loget.sim1 <- loget.sim
	points(loget.sim1 ~ loget.sim, col='green',pch=20)


	# for rhoHy

	rhoHy.sim <- unlist(lapply(sim,function(x) x$rhoHy))	
	
	rhoHy.para <- 2^unlist(lapply(para,function(x) subset(x$Parameters,grepl("log.rHy",rownames(x$Parameters)),grepl("Beta",colnames(x$Parameters)))))
	
	rhoHy.para.up.tmp <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.rHy",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[2]))
	rhoHy.para.up <- ifelse(!is.na(rhoHy.para.up.tmp), 2^rhoHy.para.up.tmp, 1)

	
	rhoHy.para.dow.tmp <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.rHy",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[1]))
	rhoHy.para.dow <- ifelse(!is.na(rhoHy.para.dow.tmp), 2^rhoHy.para.dow.tmp,0)

	res$wrongEs_rhoHy <- names(rhoHy.para)[(rhoHy.sim > rhoHy.para.up) | (rhoHy.sim < rhoHy.para.dow)]
	res$correct_rhoHy <- 1-length(res$wrongEs_rhoHy)/nGenes
	res$lm_rhoHy <- lm(rhoHy.para ~ rhoHy.sim)

	plot(rhoHy.para ~ rhoHy.sim,col='red',xlim=range(0,0.6),ylim=range(0,0.6),pch=18,asp=1)
	segments(rhoHy.sim,rhoHy.para.dow,rhoHy.sim,rhoHy.para.up,col='red')
	rhoHy.sim1 <- rhoHy.sim
	points(rhoHy.sim1 ~ rhoHy.sim, col='green',pch=20)

	#for rhoCo
	rhoCo.sim <- unlist(lapply(sim,function(x) x$rhoCo))	
	
	rhoCo.para <- 2^unlist(lapply(para,function(x) subset(x$Parameters,grepl("log.rCo",rownames(x$Parameters)),grepl("Beta",colnames(x$Parameters)))))
	
	rhoCo.para.up.tmp <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.rCo",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[2]))
	rhoCo.para.up <- ifelse(!is.na(rhoCo.para.up.tmp), 2^rhoCo.para.up.tmp, 1)

	
	rhoCo.para.dow.tmp <- unlist(lapply(para,function(x) subset(x$Confint,grepl("log.rCo",rownames(x$Confint)),grepl("Beta",colnames(x$Confint)))[1]))
	rhoCo.para.dow <- ifelse(!is.na(rhoCo.para.dow.tmp), 2^rhoCo.para.dow.tmp,0)

	res$wrongEs_rhoCo <- names(rhoCo.para)[(rhoCo.sim > rhoCo.para.up) | (rhoCo.sim < rhoCo.para.dow)]
	res$correct_rhoCo <- 1-length(res$wrongEs_rhoCo)/nGenes
	res$lm_rhoCo <- lm(rhoCo.para ~ rhoCo.sim)

	plot(rhoCo.para ~ rhoCo.sim,col='red',xlim=range(0,0.6),ylim=range(0,0.6),pch=18,asp=1)
	segments(rhoCo.sim,rhoCo.para.dow,rhoCo.sim,rhoCo.para.up,col='red')
	rhoCo.sim1 <- rhoCo.sim
	points(rhoCo.sim1 ~ rhoCo.sim, col='green',pch=20)

	result <<- res
	return(NULL)
}