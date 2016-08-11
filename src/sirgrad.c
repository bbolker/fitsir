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
# define dnu_beta_S varsdot[2]
# define dnu_gamma_S varsdot[3]
# define dnu_N_S varsdot[4]
# define dnu_I0_S varsdot[5]
# define dnu_beta_I varsdot[6]
# define dnu_gamma_I varsdot[7]
# define dnu_N_I varsdot[8]
# define dnu_I0_I varsdot[9]

# define S vars[0]
# define logI vars[1]
/* additional state variables for sensitivity equations */
# define nu_beta_S vars[2]
# define nu_gamma_S vars[3]
# define nu_N_S vars[4]
# define nu_I0_S vars[5]
# define nu_beta_I vars[6]
# define nu_gamma_I vars[7]
# define nu_N_I vars[8]
# define nu_I0_I vars[9]

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

    if (ip[0] <1) error ( "nout should be at least 1" ) ;
    Sgrad =  - scaledincid*exp(logI) ; // equation for S
    Igrad =    scaledincid - gamma ;  // equation for log (I)
}

/* variant gradient function for sensitivity equations */
void derivs_sens (int *neq , double *t , double *vars , double *varsdot ,
	     double *varsout , int *ip )
{
    double scaledincid = beta*S/N;
    double I = exp(logI);

    if (ip[0] <1) error ( "nout should be at least 1" ) ;

    Sgrad = -scaledincid*I;
    Igrad = scaledincid - gamma;
        
    double grad_SS = - beta * I/N;
    double grad_SI = - beta * S/N;
    double grad_IS = beta*I/N;
    double grad_II = beta*S/N-gamma;
    
    dnu_beta_S = grad_SS * nu_beta_S + grad_SI * nu_beta_I - S*I/N;
    dnu_N_S = grad_SS * nu_N_S + grad_SI * nu_N_I + beta*S*I/(N * N);
    dnu_gamma_S = grad_SS * nu_gamma_S + grad_SI * nu_gamma_I;
    dnu_I0_S = grad_SS * nu_I0_S + grad_SI * nu_I0_I;
    dnu_beta_I = grad_IS * nu_beta_S + grad_II * nu_beta_I + S*I/N;
    dnu_N_I = grad_IS * nu_N_S + grad_II * nu_N_I - beta*S*I/(N * N);
    dnu_gamma_I = grad_IS * nu_gamma_S +  grad_II * nu_gamma_I - I;
    dnu_I0_I = grad_IS * nu_I0_S + grad_II * nu_I0_I;
}