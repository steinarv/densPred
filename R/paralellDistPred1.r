# Performes rolling window estimateion and n step ahead predicitions on multiple series simultaniously

# In these function ... is to controll the solver options (solver.method & solver.control)

#In simpleDoRoll: Some functions only return quantiles, some only variance. Needs to be fixed

simpleDoRoll <- function(x, f=stdGARCH, pval=c(0.01, 0.025), n.ahead=10, win.size=500, dates=NULL, 
						aggf=aggprod, ... ){
	n <- length(x)
	if(!is.null(dates))if(length(dates)!=n){
							warning("Length of dates array do not match length of data array")
							dates=NULL
						}
	
	num <- floor((n-win.size)/n.ahead)
  
	dfOut <- data.frame(Convergence = logical(num), Emp_ret = numeric(num), Exp_var = numeric(num),
						Date=as.Date(NA), matrix(NA, ncol=length(pval), nrow=num))
	colnames(dfOut)[-c(1,2,3,4)] <- paste0(pval)
  
	a <- 1
	b <- win.size
	i <- 1
	while(b <= n-n.ahead){
    
		# f() standard optmizing function for package densPred
		fit <- f(x[a:b], n.ahead=n.ahead, ... )
		if(fit$convergence!=0){
			warning(paste0("Convergence is: ", fit$convergence, 
							" for window: ", a, " to ", b, ". Using paramters from previous optimization."))
			if(a==1)stop("Convergance not possible for first window")
			else fit <- f(x[a:b], n.ahead=n.ahead, param=prevPar, doOptim=FALSE)
			bConv = FALSE
		}else{
			prevPar <- fit$par
			bConv = TRUE
		}
    
		vHt <- sum(fit$HtOut)
		
		if(!is.null(dates))dfOut[i, "Date"] <- dates[b]		#The date for the last obs in sample
		dfOut[i, c("Convergence", "Emp_ret", "Exp_var", paste0(pval))] <- 
								c(bConv, aggf(x[(b+1):(b+n.ahead)]), vHt, qnorm(pval, sd=sqrt(vHt)))
                      
		i <- i+1
		a <- a+n.ahead
		b <- b+n.ahead
	}
  
	dfOut
}


# Takes a series and do rolling window optmiziation and out of sample predictions and evaluation with (possible) several models
parallelTest.byFunc <- function(ser, serName="series1", vf=c(stdGARCH, eGARCH), vfNames=paste0("model", 1:length(vf)), dates=NULL,
								pval=c(0.01, 0.025), n.ahead=10, win.size=500, ret.type="arithmetic", freq="freq", 
								... ){
	n <- length(vf)
	if(length(vfNames)!=n){
		warning("Number of function names do not match the number of functions to be evaluated")
		vfNames=rep("", n)
	}
	
	if(!is.null(dates))if(length(dates)!=length(ser)){
							warning("Length of dates array do not match length of data array")
							dates=NULL
						}else{
							daterange <- format(range(dates), "%Y%m%d")
						}
	else{
		daterange <- rep("00000000" ,2)
	}
	if(any(is.na(daterange)))daterange <- rep("00000000" ,2)
	
	
	aggf = switch(ret.type, "arithmetic"=aggprod, "logarithmic"=aggsum)
	
	ptm <- proc.time()
	cl <- makeCluster(rep("localhost", n), type = "SOCK") #n should not exceed number of cores
	registerDoSNOW(cl)

	lOut <- foreach(i=1:n, .packages="densPred", .inorder=TRUE) %dopar% {  
		
		f <- vf[[i]]
		
		dfX <- simpleDoRoll(ser, f=f, pval=pval, n.ahead=n.ahead, win.size=win.size, aggf=aggf, dates, ...)
		res <- varTest(pval, dfX[, "Emp_ret"], dfX[, paste0(pval)])
		
		# Collection of evaluation scores for distributional predictions
		mPass <- mean(res[ ,-c(1,2)]>0.05)
		medAPE <- median(abs((dfX$Emp_ret[-1]^2-dfX$Exp_var[-1])/dfX$Emp_ret[-1]^2))
		AMAPE <- mean(abs((dfX$Emp_ret[-1]^2-dfX$Exp_var[-1])/(dfX$Emp_ret[-1]^2+dfX$Exp_var[-1])))
		
		# Identification string for this analysis. 
		rObjName <- paste0( serName, "_", vfNames[i], "_",
							paste(daterange, collapse=""), "_",
							n.ahead, "_", win.size, "_", freq, "_",
							ret.type)
  
		list(series=serName, model=vfNames[i], rObj=rObjName,
			num.pred= nrow(dfX), n.non.conv=sum(!dfX[ ,"Convergence"]), 
			qScore=res, mPass=mPass, medAPE=medAPE, AMAPE=AMAPE, 
			daterange=daterange, win.size=win.size, ret.type=ret.type, n.ahead=n.ahead, 
			tot.obs=length(ser), freq=freq, allData=dfX) 
	}#foreach

	stopCluster(cl)
	print(proc.time() - ptm)
	
	names(lOut) <- vfNames
	lOut
}

# Takes (possible) several series and do rolling window optmiziation and out of sample predictions and evaluation
parallelTest.bySeries <- function(ser, serNames=paste0("series", 1:length(ser)), f=stdGARCH, fName="model", dates=NULL,
									pval=c(0.01, 0.025), n.ahead=10, win.size=500, ret.type="arithmetic", 
									freq="freq", ... ){
	n <- length(ser)
	if(length(serNames)!=n){
		warning("Number of series names do not match the number of series")
		serNames=rep("", n)
	}
	
	aggf = switch(ret.type, "arithmetic"=aggprod, "logarithmic"=aggsum)
	
	ptm <- proc.time()
	cl <- makeCluster(rep("localhost", n), type = "SOCK") #n should not exceed number of cores
	registerDoSNOW(cl)

	lOut <- foreach(i=1:n, .packages="densPred", .inorder=TRUE) %dopar% {  

		dfX <- simpleDoRoll(ser[[i]], f=f, pval=pval, n.ahead=n.ahead, win.size=win.size, aggf=aggf, dates[[i]], ...)
		res <- varTest(pval, dfX[, "Emp_ret"], dfX[, paste0(pval)])
		
		#Dates must be handled within loop since the date
		if(!is.null(dates[[i]]))if(length(dates[[i]])!=length(ser[[i]])){
							warning("Length of dates array do not match length of data array")
							dates[[i]]=NULL
						}else{
							daterange <- format(range(dates[[i]]), "%Y%m%d")
						}
		else{
			daterange <- rep("00000000" ,2)
		}
		if(any(is.na(daterange)))daterange <- rep("00000000" ,2)		
		#------------------------------------------------
		
		
		# Collection of evaluation scores for distributional predictions
		mPass <- mean(res[ ,-c(1,2)]>0.05)
		medAPE <- median(abs((dfX$Emp_ret[-1]^2-dfX$Exp_var[-1])/dfX$Emp_ret[-1]^2))
		AMAPE <- mean(abs((dfX$Emp_ret[-1]^2-dfX$Exp_var[-1])/(dfX$Emp_ret[-1]^2+dfX$Exp_var[-1])))
		
		# Identification string for this analysis. 
		rObjName <- paste0( serNames[i], "_", fName, "_",
							paste(daterange, collapse=""), "_",
							n.ahead, "_", win.size, "_", freq, "_",
							ret.type)
  
		list(series=serNames[i], model=fName, rObj=rObjName,
			num.pred= nrow(dfX), n.non.conv=sum(!dfX[ ,"Convergence"]), 
			qScore=res, mPass=mPass, medAPE=medAPE, AMAPE=AMAPE, 
			daterange=daterange, win.size=win.size, ret.type=ret.type, n.ahead=n.ahead, 
			tot.obs=length(ser[[i]]), freq=freq, allData=dfX) 
	}#foreach

	stopCluster(cl)
	print(proc.time() - ptm)
	
	names(lOut) <- paste0("l", serNames)
	lOut
}