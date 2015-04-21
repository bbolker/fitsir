/* SIR_grad_logI .c */
/* Author : Irena Papst */
/* SIR gradient in C for use with the fitsir package in R */
# include <R.h> // include R library
static double parms[3];
# define beta parms[0]
# define gamma parms[1]
# define N parms[2]

/* initializer */
void initmod ( void (*odeparms)(int *, double *) )
{
    int num=3;
    odeparms (&num, parms) ;
}

/* Derivatives and 1 output variable */
void derivs (int *neq , double *t , double *vars , double *varsdot ,
	     double *varsout , int *ip )
{
    if (ip[0] <1) error ( "nout should be at least 1" ) ;
    varsdot[0] = - beta * vars[0]* exp (vars[1]) / N ; // equation for S
    varsdot[1] = beta * vars[0]/ N - gamma ; // equation for log (I)
}

# define vbeta vparms[0]
# define vgamma vparms[1]
# define vN vparms[2]

void derivs0 (double *t , double *vars , double *varsdot, double *vparms ) 
{
    varsdot[0] = - vbeta * vars[0]* exp (vars[1]) / vN ; // equation for S
    varsdot[1] = vbeta * vars[0]/ vN - vgamma ; // equation for log (I)
}
