#ifndef _HAVE_STRUCTURE_H
#define _HAVE_STRUCTURE_H

#include <dan.h>
#include <rdan.h>
#include <unistd.h>
#include <assert.h>
#include <sys/time.h>
#include <ran.h>
#include <io.h>
#include <iogeno.h>

#define PARALLEL FALSE

#if PARALLEL
#include <mpi.h>
#endif

#define BIT_ENCODING FALSE // don't use BIT_ENCODING with STANDALONE!
#define NALS 2 // diallelism is hard-coded at the moment
#define OUTPUT_PERIOD 1000 // write current parameter values every OUTPUT_PERIOD iterations
#define VERBOSE TRUE

#define ROOT 0
#define RUNIF rnd
#define RDIRICHLET RDirichlet
#define FIXED_INIT FALSE
enum MSG_TAGS { SNP_RANGE_MSG_TAG } ;
enum MODEL { EMAM, MCMCAM, VBAM0, VBAM1, VBNAM0, VBNAM1} ; // EM Admixture Model, VB No Admixture Model

struct options {
    int n ;
    int L ;
    int K ;
    double q_prior_strength ; // VBAM: 0 < prior_strength < 1: 0 => prior counts = 0; 1 => prior counts = 2L/K
    double stop_tol ;
    double *fit ;
    double *p ;
    int plateau_size ;
    bool obey_stop_tol ;
    bool uncertain_data ;
    enum MODEL model ;
    bool diploid ;
    char *outdir ;
    char fnamebuff[1000] ;
    bool show_output ;
    int output_period ;
    int print_period ;
    double *snploglike ; // temporarily here
    
    int niters ;
    int burnin ;
    int thin ;
    int ndraws ;
    
    char *popfile ;
    bool *updatez ; // flags indicating whether individuals' populations are fixed
    int *fixpop ;
    double *w ; // population mean admixture proportions
    bool popmean ;
    bool spatial ;
    double *skernel ;
    FILE *skernelfile ;

    /* skip loci approximation*/
    bool *skip_locus ;
    int nskip ;
    double *skipfreq ;
    int skip_update_interval ;
    double skip_prior ;
    double *freqs ;
    double *snploglike_nostructure ;
    double skip_tol ;

    FILE *timefile, *statefile, *wstatefile ;
    
} opt ;


#define SMALL 1e-50
#define NUDGE 1e-4

/*
  #if BIT_ENCODING
  #define MASK_234 192
  #define MASK_134 48
  #define MASK_124 12
  #define MASK_123 3
  const unsigned char MASK[4] = {MASK_234,MASK_134,MASK_124,MASK_123} ;
  #endif
*/

int get_allele(const unsigned char **x, int a) ;
void set_first_snps(int L, int num_processes, int *first_snps) ;
bool reached_plateau(int iter, double *fit) ;
double *rprior_mu(int L, int K) ;
double *rprior_q(int n, int K) ;
void compute_freqs_and_snploglike_nostructure(unsigned char *x, double *p, double *loglike) ;
char *model_name() ;

extern int L_tot, thisproc ;

#endif // _HAVE_STRUCTURE_H
