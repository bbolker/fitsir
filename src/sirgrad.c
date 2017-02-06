/* SIR_grad_logI .c */
/* Author : Irena Papst */
/* SIR gradient in C for use with the fitsir package in R */
# include <R.h> // include R library
# define NPARMS 3
static double parms[NPARMS];
# define beta parms[0]
# define gamma parms[1]
# define N parms[2]

# define Sgrad varsdot[0]
# define Igrad varsdot[1]
# define dnu_S_b varsdot[2]
# define dnu_S_g varsdot[3]
# define dnu_S_N varsdot[4]
# define dnu_S_i varsdot[5]
# define dnu_I_b varsdot[6]
# define dnu_I_g varsdot[7]
# define dnu_I_N varsdot[8]
# define dnu_I_i varsdot[9]

# define S vars[0]
# define logI vars[1]
/* additional state variables for sensitivity equations */
# define nu_S_b vars[2]
# define nu_S_g vars[3]
# define nu_S_N vars[4]
# define nu_S_i vars[5]
# define nu_I_b vars[6]
# define nu_I_g vars[7]
# define nu_I_N vars[8]
# define nu_I_i vars[9]

/* initializer */
void initmod ( void (*odeparms)(int *, double *) )
{

    int num=NPARMS; /* number of parameters, hard-coded */
    odeparms (&num, parms) ;
}

/* Derivatives and 1 output variable */
void derivs (int *neq , double *t , double *vars , double *varsdot ,
	     double *varsout , int *ip )
{
    double scaledincid = beta*S/N;

    Sgrad =  - scaledincid*exp(logI) ; // equation for S
    Igrad =    scaledincid - gamma ;  // equation for log (I)
}

/* variant gradient function for sensitivity equations */
void derivs_sens (int *neq , double *t , double *vars , double *varsdot ,
	     double *varsout , int *ip )
{
    double scaledincid = beta*S/N;
    double I = exp(logI);
    
    Sgrad = -scaledincid*I;
    Igrad = scaledincid - gamma;
        
    double grad_SS = - beta * I/N;
    double grad_SI = - beta * S/N;
    double grad_IS = beta*I/N;
    double grad_II = beta*S/N-gamma;
    
    dnu_S_b = grad_SS * nu_S_b + grad_SI * nu_I_b - S*I/N;
    dnu_S_N = grad_SS * nu_S_N + grad_SI * nu_I_N + beta*S*I/(N * N);
    dnu_S_g = grad_SS * nu_S_g + grad_SI * nu_I_g;
    dnu_S_i = grad_SS * nu_S_i + grad_SI * nu_I_i;
    dnu_I_b = grad_IS * nu_S_b + grad_II * nu_I_b + S*I/N;
    dnu_I_N = grad_IS * nu_S_N + grad_II * nu_I_N - beta*S*I/(N * N);
    dnu_I_g = grad_IS * nu_S_g +  grad_II * nu_I_g - I;
    dnu_I_i = grad_IS * nu_S_i + grad_II * nu_I_i;
}