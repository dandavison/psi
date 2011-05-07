#define MATHLIB_STANDALONE TRUE
#include <Rmath.h>

void mcmcam(const unsigned char *x0, int n, int L, int K, double *qmean) ;
void write_current_state_mcmc(double *q, int iter) ;
void rdirichlet_update(double *observed, double *prior, int K, double *p) ;
void update_skip_indicators() ;
