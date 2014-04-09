# -------------------------------- Basic variance models -------------------------------------------
stdNormVar <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n), HtOut=rep(h, n.ahead), convergence=0)
}

#optTdf(x, dfrange=c(1,00))sum(log(dt(x, df=df)))
#stdTVar <- function(x, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
#	n <- length(x)
#	h <- var(x)
#	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
#}

stdHS <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}

stdEWMA <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}

filterHS <- function(x, n.ahead=1, ... ){
	n <- length(x)
	h <- var(x)
	list(HtIn=rep(h, n) , HtOut=rep(h, n.ahead) , convergence=0)
}
# --------------------------   Estimates simple GARCHish models -------------------------------------

OPTstdGARCH <- function(x, h0, param){
	.Call("LLstdGARCH", x, h0, param, PACKAGE = "densPred" )
}

OPTeGARCH <- function(x, h0, param){
	.Call("LLeGARCH", x, h0, param, PACKAGE = "densPred" )
}

testEgarch <- function(x){
	.Call("HTeGARCH", x, var(x), c( log(var(x))*0.25-0.016, 0.02, 0, 0.75 ), 1, PACKAGE = "densPred" )
}

stdGARCH <- function(x, h0=NULL, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	
	if(is.null(h0))h0 <- var(x);
	if(is.null(param)){param <-  log(c(var(x)*0.05, 0.2, 0.75))}else{param <- log(param)}
	if(doOptim){
		opt <- optim(param, OPTstdGARCH, x=x, h0=h0, method=solver.method, control=solver.control)
		param <- opt$par
		opt$par <- exp(opt$par)
		opt$value <- -opt$value
	}
	
	Ht <- .Call("HTstdGARCH", x, h0, param, n.ahead, PACKAGE = "densPred" )
	
	
	lOut <- list(HtIn=Ht[1:n], HtOut=Ht[(n+1):(n+n.ahead)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
	
}

eGARCH <- function(x, h0=NULL, param=NULL, doOptim=TRUE, n.ahead=1, solver.method="Nelder-Mead", solver.control=list()){
	n <- length(x)
	
	if(is.null(h0))h0 <- var(x);
	if(is.null(param)) param <- c( log(var(x))*0.25-0.016, 0.02, 0, 0.75 ) #a0=log(var(x))*(1-b0)-a1*sqrt(2/pi)
	if(doOptim){
		opt <- optim(param, OPTeGARCH, x=x, h0=h0, method=solver.method, control=solver.control)
		param <- opt$par
		opt$value <- -opt$value
	}

	Ht <- exp(.Call("HTeGARCH", x, h0, param, n.ahead, PACKAGE = "densPred" ))
	
	
	lOut <- list(HtIn=Ht[1:n], HtOut=Ht[(n+1):(n+n.ahead)])
	if(doOptim)lOut <- c(lOut, opt)
	
	lOut
	
}