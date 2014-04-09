# Simple function for aggregating returns accross time, and making step ahead predictions.
# For n days ahead predictions based on daily returns, the returns are aggreagetd to make n
# days return, and one step a head forcasting functins is called. 

simpleAggMain <- function(x, pval=0.025, n.ahead=1, win.size=500, f=qreg, Date=NULL,
                 agfunc = aggprod, ... ){
  
  n <- length(x)
  
  rest <- win.size%%n.ahead
  if(rest>0){
    warning("n.ahead not a multiple of win.size")
    win.size <- (win.size-rest)/n.ahead
    
  }
  
  rest <- length(x)%%n.ahead
  if(rest>0){
    x <- x[-((n-rest+1):n)]
    if(!is.null(Date))Date <- Date[-((n-rest+1):n)]
  }
  
  x <- aggregate(x, list(rep(1:(length(x)/n.ahead), each=n.ahead)), 
                 FUN=agfunc)$x
  n <- length(x)
  
  if(!is.null(Date)){
    Date <- aggregate(Date, by=list(rep(1:(length(Date)/n.ahead), each=n.ahead)),
                      FUN=max)
    Date <- as.Date(Date$x, origin="1970-01-01")
  }
  
  
  n.out <- n-win.size #Number of out of sample obs
  dfOut <- data.frame(Date=as.Date(NA), matrix(NA, ncol=1+length(pval), 
                                               nrow=n.out)) 
                     #Date, y, var1, var2...
  if(!is.null(Date))dfOut$Date <- Date[(win.size+1):n]
  
  for(i in win.size:(n-1)){
    j <- i-win.size+1
    dfOut[j, 2] <- x[i+1] #Predicted observation (not handed to f())
    
    for(pv in 1:length(pval)){
    dfOut[j, 2+pv] <- f(x=x[j:i], pval=pval[pv], ... )
    }
                   
  }
  
  dfOut
}