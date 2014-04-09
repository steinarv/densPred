#include "varianceModels.h"

using namespace Rcpp ;

// ---------------------------- Returns log likelihood for standard GARCH model --------------------------------------
SEXP LLstdGARCH(SEXP X, SEXP H0, SEXP Param) {
	double dH0 = as<double>(H0);
	NumericVector nvX(X); NumericVector nvParam(Param);
	  
	int n = nvX.size();
	double a0 = ::exp(nvParam(0)); double a1 = ::exp(nvParam(1));
	double b0 = ::exp(nvParam(2));
	
	if((a1+b0)>1)return(wrap(100000));
	  
	NumericVector H(n);
	H(0) = dH0;
	  
	double logL=0;
	  
	for(int i=1;i<n;i++){
		H(i)=a0+a1*pow(nvX(i-1), 2)+b0*H(i-1);
		logL+=::log(H(i))+pow(nvX(i),2)/H(i); //The correct way is to subtract rather than add 
	}
	  
	return wrap(logL);
}
// -------------------------------------------------------------------------------------------------------------------

// ---------------------------- Returns conditional variance for standard GARCH model --------------------------------------
SEXP HTstdGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) {
	double dH0 = as<double>(H0); double dNout = as<double>(Nout);
	NumericVector nvX(X); NumericVector nvParam(Param);
	
	int n = nvX.size();
	double a0 = ::exp(nvParam(0)); double a1 = ::exp(nvParam(1));
	double b0 = ::exp(nvParam(2));
  
	NumericVector H(n+dNout);
	H(0) = dH0;
  
	for(int i=1;i<(n+dNout);i++){
		//std::cout << "a0: " << a0 << " a1: " << a1 << " b0: " << b0 << std::endl;
		//std::cout << "pow(nvX(i-1),2)   " << pow(nvX(i-1),2) << std::endl;
		//std::cout << "a0+a1*pow(nvX(i-1),2)+b0*H(i-1)   " << a0+a1*pow(nvX(i-1),2)+b0*H(i-1) << std::endl;
		if(i<=n)H(i)=a0+a1*pow(nvX(i-1),2)+b0*H(i-1);
		else H(i)=a0+(a1+b0)*H(i-1); //Out of sample predictions
	}
  
	return wrap(H);
}
// -------------------------------------------------------------------------------------------------------------------

// --------------------------------- Returns log likelihood for EGARCH model -----------------------------------------
SEXP LLeGARCH(SEXP X, SEXP H0, SEXP Param) {
	double dH0 = as<double>(H0);
	NumericVector nvX(X); NumericVector nvParam(Param);
	
	int n = nvX.size();
	double a0 = nvParam(0); double a1 = nvParam(1);
	double gamma = nvParam(2); double b0 = nvParam(3);
  
	NumericVector logH(n);
	logH(0) = ::log(dH0);
	
	double logL=0;
  
	for(int i=1;i<n;i++){
		logH(i)=a0+a1*(std::abs(nvX(i-1))+gamma*nvX(i-1))/sqrt(::exp(logH(i-1)))+b0*logH(i-1);
		logL+=logH(i)+pow(nvX(i),2)/::exp(logH(i));		//The correct way is to subtract rather than add 
	}
  
	return wrap(logL);
}

// -------------------------------- Returns times series of conditional EGARCH variance ---------------------------------
SEXP HTeGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) {
	double dH0 = as<double>(H0); double dNout = as<double>(Nout);
	NumericVector nvX(X); NumericVector nvParam(Param);	
	
	int n = nvX.size();
	double a0 = nvParam(0); double a1 = nvParam(1);
	double gamma = nvParam(2); double b0 = nvParam(3);
  
	NumericVector logH(n+dNout);
	logH(0) = ::log(dH0);
	
	//std::cout << "dH0    " << dH0 << std::endl;
	//std::cout << "::log(dH0)   " << ::log(dH0)  << std::endl;

	for(int i=1;i<(n+dNout);i++){
		if(i<=n) logH(i)=a0+a1*(std::abs(nvX(i-1))+gamma*nvX(i-1))/sqrt(::exp(logH(i-1)))+b0*logH(i-1);
		else logH(i)=a0+a1*sqrt(2/M_PI)+b0*logH(i-1);
	}
  
	return wrap(logH);
}