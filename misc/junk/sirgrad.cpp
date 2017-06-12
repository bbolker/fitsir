#include <Rcpp.h>     // header files, definitions etc.
using namespace Rcpp; // can use Rcpp functions without doing Rcpp::

// we can't use with() any more, we have to define
//  the positions of the parameters; indexing from 0
// "compiler"/macro definitions, very fast
#define beta pars[0]
#define gamma pars[1]
#define N pars[2]
#define S y[0]
#define logI y[1]
#define Sdot res[0]
#define logIdot res[1]

// some stuff I looked at while I was figuring this out
// http://gallery.rcpp.org/articles/sugar-for-high-level-vector-operations/
// http://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf

// telling compiler that we want to use this function in R
// [[Rcpp::export]]  // magic 

List sirgrad(double t, NumericVector y, NumericVector pars) {
    int n = y.size();
    NumericVector res( n );
    Sdot = -beta/N*exp(logI)*S;
    logIdot = beta/N*S-gamma;
    List ret;   // make a new list
    ret["grad"] = res;  // put the result into it
    return(ret);

}

// initializer (verbatim from sirgrad.c)
// how do I make this look like an Rcpp function?

static double parms[3];
void initmod ( void (*odeparms)(int *, double *) )
{
    int num=3;
    odeparms (&num, parms) ;
}
