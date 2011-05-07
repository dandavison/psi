#define MATHLIB_STANDALONE TRUE
#include <Rmath.h>
// #include "digamma.h"

#define DIGAMMA my_digamma // digamma_lookup

double my_digamma(double x) ;
double vbam1(const unsigned char *x, double *px, int n, int L, int K, const double *lambda0, double *lambda1, const double *alpha0, double *alpha1) ;
void vbnam0(const unsigned char *x, int *_n, int *_L, int *_K, bool *_haploid,
	    double const *lambda_prior, double *lambda1,
	    double const *alpha_prior, double *alpha1,
	    double *pz,
	    double *_free_energy, 
	    double *tol, 
	    int *niters) ;
void vbam0(unsigned char *x, int *_n, int *_L, int *_K, bool *_haploid,
	   const double *lambda_prior, double *lambda,
	   const double *alpha_prior, double *alpha, 
	   double *_free_energy) ;
double vbnam1(const unsigned char *x, double *px, int n, int L, int K, const double *lambda0, double *lambda1, const double *alpha0, double *alpha1, double *pz) ;
int get_allele(const unsigned char **x, int a) ;
void write_current_state_vbam(double *lambda, double *alpha, double free_energy, int iter, int K, int n, int L) ;
void write_current_state_vbnam(double *pz, double *lambda, double *alpha, double free_energy, int iter, int K, int n, int L) ;
double KL_Dirichlets(double *w, const double *v, int K) ; // second one happens to be prior in VB application; thus const. Don't really understand how to use const.
double free_energy(const unsigned char *x, double *pz,  const double *alpha0, double *alpha1, const double *lambda0, double *lambda1, 
		   int n, int K, int L, int iter) ;
double *initialise_lambda_prior_vbam(int K, int n) ;
double *initialise_lambda_prior_vbnam(int K, int n) ;
double *initialise_alpha_prior(int K, int L) ;
