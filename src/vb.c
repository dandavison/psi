#include "psi.h"
#include "vb.h"

double my_digamma(double x) {
    return digamma(x) ;
}

double *initialise_alpha_prior(int K, int L) {
    int k, l, a ;
    double *alpha0 = ALLOC(L * K * NALS, double) ;
    const double DELTA = 1e-1 ;
    
    /* set hyperparameters of priors */
    for(l = 0 ; l < L ; l++)
	for(k = 0 ; k < K ; k++)
	    for(a = 0 ; a < NALS ; a++)
		alpha0[l*K*NALS + k*NALS + a] = runif(1 - DELTA, 1 + DELTA) ;
    
    return alpha0 ;
}

double KL_Dirichlets(double *w, const double *v, int K) { // second one happens to be prior in VB application; thus const. Don't really understand how to use const.
    /* Kullback-Leibler divergence between two Dirichlets with parameters w_1,...,w_K and v_1,...,v_K
       Rezek & Roberts et al. variational Bayes HMM book chapter
       http://www.robots.ox.ac.uk/~irezek/Outgoing/Papers/varhmm.ps.gz */
    
    double d=0.0, sumw=0.0, sumv=0.0 ;
    int k ;

    for(k = 0 ; k < K ; ++k, ++w, ++v) { 
	sumw += *w ; 
	sumv += *v ;
	d += lgammafn(*v) - lgammafn(*w)  +  (*w - *v) * DIGAMMA(*w) ;
    }

    d += lgammafn(sumw) - lgammafn(sumv)  -  (sumw - sumv) * DIGAMMA(sumw) ;

    return d ;
}
