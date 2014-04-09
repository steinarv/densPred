#ifndef _densPred_var_H
#define _densPred_var_H

#include <RcppArmadillo.h>

RcppExport SEXP LLstdGARCH(SEXP X, SEXP H0, SEXP Param) ;
RcppExport SEXP HTstdGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) ;

RcppExport SEXP LLeGARCH(SEXP X, SEXP H0, SEXP Param) ;
RcppExport SEXP HTeGARCH(SEXP X, SEXP H0, SEXP Param, SEXP Nout) ;


#endif