#include "psi.h"
#include "em.h"
#include "mcmc.h"
#include "vb.h"


/*
  Entry routine for fitting the structure Admixture Model, by EM, MCMC or VB.
  
  Including parallel algorithm via MPI:

  0. Process 0 allocates SNPs to the other processes, and informs them of their allocation.
  1. Each process reads in its subset of the data, and sets starting values of parameters.
  2. Each process updates (beliefs about) the portion of the missing data z corresponding to its own SNPs,
  and also updates (beliefs about) the cluster allele frequencies at those SNPs.
  3. At the end of each iteration, process 0 gathers information about the admixture proportions q from each 
  other process, averages them, and sends the averaged information back out to the processes.
*/



int L_tot, thisproc ;
// double memory=0.0, gigs_per_byte = 9.313226e-10 ; //pow(2.0, -30.0)

int main(int argc, char *argv[]) {
    int n, twon, K ;
    int i, l ;
    int *first_snps=NULL, snp_range[2] ;
    unsigned char *x=NULL ;
    double *px=NULL ;
    char *genotypes_file=NULL, *chiamo_file=NULL, *q_infile=NULL, *mu_infile=NULL, fnamebuff[1000] ;
    int nprocs, c ;
    bool *snp_include ;
    FILE *f ;

    srand(time(NULL) + getpid()) ;
    set_seed(time(NULL) + getpid(), time(NULL) + getpid() % 17) ;

    // srand(1) ;
    // set_seed(1,1) ;

#if PARALLEL
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &thisproc) ;
#else
    nprocs = 1 ;
    thisproc = ROOT ;
#endif
    
    opt.popmean = FALSE ;
    opt.spatial = FALSE ;

    opt.stop_tol = 1e-7 ;
    opt.obey_stop_tol = TRUE ;
    opt.diploid = TRUE ;
    opt.show_output = TRUE ;
    opt.output_period = 100 ;
    opt.uncertain_data = FALSE ;
    opt.plateau_size = 5 ;
    opt.burnin = 500 ;
    opt.thin = 1 ;
    opt.skip_prior = 0 ;
    opt.skip_update_interval = 1000 ;
    opt.skip_tol = 0.01 ;
    opt.print_period = 100 ;
    opt.popfile = NULL ;

    while((c = getopt(argc, argv, "a:b:f:g:h:K:L:m:n:o:p:q:r:s:t:v:z:")) != -1) {
	switch(c) {
	case 'a':
	    opt.print_period = atoi(optarg) ; break ;
	case 'b':
	    opt.burnin = atoi(optarg) ; break ;
	case 'g':
	    assert(chiamo_file == NULL) ;
	    genotypes_file = optarg ;
	    opt.uncertain_data = FALSE ;
	    break ;
	case 'h':
	    opt.thin = atoi(optarg) ; break ;
	case 'K':
	    K = opt.K = atoi(optarg) ; break ;
	case 'L':
	    L_tot = atoi(optarg) ; break ;
	case 'm':
	    mu_infile = optarg ; break ;
	case 'n':
	    n = opt.n = atoi(optarg) ; break ;
	case 'o':
	    opt.outdir = optarg ; break ;
	case 'p':
	    assert(genotypes_file == NULL) ;
	    chiamo_file = optarg ;
	    opt.uncertain_data = TRUE ;
	    break ;
	case 'q':
	    q_infile = optarg ; break ;
	case 'r':
	    opt.skip_prior = atof(optarg) ; break ;
	case 's':
	    opt.spatial = TRUE ;
	    opt.skernelfile = optarg ;
	    break ;
	case 't':
	    opt.niters = atoi(optarg) ; break ;
	case 'v':
	    opt.model = atoi(optarg) ; break ;
	case 'z':
	    opt.popfile = optarg ; break ;
	case '?':
	    ERROR("Unrecognised option") ;
	}
    }
    twon = 2*n ;
    opt.fit = ALLOC(opt.niters, double) ;
    opt.q_prior_strength = 0.001 ;

    if(thisproc == ROOT) {
	/* Record arguments */
	sprintf(fnamebuff, "%s/argv", opt.outdir) ;
	f = fopen(fnamebuff, "w") ;
	for(i = 0 ; i < argc ; i++) fprintf(f, "%s ", argv[i]) ;
	fputc('\n', f) ;
	fclose(f) ;
	
	if(VERBOSE) {
	    PRINT("argv: ") ; for(i = 0 ; i < argc ; i++) PRINT("%s ", argv[i]) ; PRINT("\n") ;
	    PRINT("n = %d\nL_tot = %d\nK = %d\nnprocs = %d\nuncertain_data = %s\n",
		  n, L_tot, K, nprocs, opt.uncertain_data ? "TRUE":"FALSE") ;
	    PRINT("data_file = %s\nq_infile = %s\nmu_infile = %s\noutdir = %s\n", 
		  opt.uncertain_data ? chiamo_file : genotypes_file, 
		  q_infile != NULL ? q_infile : "NULL", mu_infile != NULL ? mu_infile : "NULL", opt.outdir) ;
	    PRINT("thin = %d\n", opt.thin) ;
	    PRINT("fixpopfile = %s\n", opt.popfile) ;
	    PRINT("popmean = %s\n", opt.popmean ? "TRUE": "FALSE") ;
	}
	    
	/* Process 0 allocates SNPs to the other processes, and informs them of their allocation. */
	first_snps = ALLOC((nprocs+1), int) ;
	set_first_snps(L_tot, nprocs, first_snps) ;
	first_snps[nprocs] = L_tot ;
    }


#if PARALLEL
    MPI_Scatter(first_snps,     1, MPI_INT, snp_range,     1, MPI_INT, ROOT, MPI_COMM_WORLD) ;
    MPI_Scatter(first_snps + 1, 1, MPI_INT, snp_range + 1, 1, MPI_INT, ROOT, MPI_COMM_WORLD) ;
#else
    snp_range[0] = first_snps[0] ;
    snp_range[1] = first_snps[1] ;
#endif
    
    /*  Each process reads in its subset of the data, and sets starting values of parameters */
    opt.L = snp_range[1] - snp_range[0] ;

    snp_include = ALLOC(L_tot, int) ;
        
    for(l = 0 ; l < L_tot ; l++) snp_include[l] = (l >= snp_range[0] && l < snp_range[1]) ;
    
    if(opt.uncertain_data) {
	px = ALLOC(3 * n * opt.L, double) ; 
	f = fopen(chiamo_file, "r") ;
	read_chiamo(f, n, L_tot, snp_include, px) ;
	fclose(f) ;
    }
    else x = read_genotypes_snp_range_all_indivs(genotypes_file, n, snp_range[0], snp_range[1]) ;
    
    if(VERBOSE) PRINT("%s model fit: (n,L,K) = (%d, %d, %d)\n", model_name(), n, opt.L, K) ;

    opt.snploglike = ALLOC(opt.L, double) ;
    opt.skip_locus = CALLOC(opt.L, bool) ;
    opt.skipfreq = CALLOC(opt.L, double) ;
    opt.snploglike_nostructure = ALLOC(opt.L, double) ;
    opt.freqs = ALLOC(opt.L, double) ;
    compute_freqs_and_snploglike_nostructure(x, opt.freqs, opt.snploglike_nostructure) ;
    sprintf(opt.fnamebuff, "%s/times", opt.outdir) ;
    opt.timefile = fopen(opt.fnamebuff, "w") ;
    sprintf(opt.fnamebuff, "%s/states", opt.outdir) ;
    opt.statefile = fopen(opt.fnamebuff, "w") ;
    sprintf(opt.fnamebuff, "%s/wstates", opt.outdir) ;
    opt.wstatefile = fopen(opt.fnamebuff, "w") ;
    
    if(opt.model == EMAM) {

	/* EM admixture model fit */

	double *mu, *mu_anc=NULL ;
	double logpost ;
	double *F = initialise_F(K) ;
	double *q ;
	
	if(mu_infile != NULL) mu = read_mu(mu_infile, K, opt.L, L_tot, snp_include) ;
	else mu = rprior_mu(opt.L, K) ;
	
	if(q_infile != NULL) q = read_matrix_double(q_infile, n, K, "%lf") ;
	else q = rprior_q(n, K) ;
    
	// mu_anc = ALLOC(opt.L, double) ;
	
	PRINT("starting EMAM model fit\n") ;
	logpost = emam(x, px, &n, &opt.L, &K, q, mu, mu_anc, F) ;
	write_current_state_emam(q, mu, logpost, -1, K, n, opt.L) ;
	free(mu) ;
	free(q) ;
	free(F) ;
	// free(mu_anc) ;
    }
    else if(opt.model == MCMCAM) {
	double *qmean = CALLOC(n * K, double) ;
	opt.niters = opt.niters * opt.thin ;

	opt.updatez = ALLOC(n, bool) ;
	if(opt.popfile != NULL) {
	    opt.fixpop = read_matrix_int(opt.popfile, n, 1, "%d") ;
	    for(i = 0 ; i < n ; i++) opt.updatez[i] = (opt.fixpop[i] < 0 || opt.fixpop[i] >= K) ;
	}
	else for(i = 0 ; i < n ; i++) opt.updatez[i] = TRUE ;

	if(opt.spatial) opt.skernel = read_matrix_double(opt.skernelfile, n, n, "%lf") ;

	mcmcam(x, n, opt.L, K, qmean) ;
	write_current_state_mcmc(qmean, -1) ;
	free(qmean) ;
    }
    else if(opt.model >= VBAM0 && opt.model <= VBNAM1){
	/* VB model fit */
	
	if(opt.show_output) {
	    PRINT("%5s %15s %15s %15s %15s %15s %15s\n", 
		  "iter", "d_KL" , "E log q(z)" , "E log like" , "entropy", "free energy", "% increase") ;
	    PRINT("%5s %15s %15s %15s %15s %15s %15s\n", 
		  "----", "----" , "----------" , "----------" , "-------", "-----------", "-------------") ;
	}
	const double *alpha0 = initialise_alpha_prior(K, opt.L) ;
	double *alpha1 = ALLOC(opt.L * K * NALS, double) ;                         // [1:J, 1:K, 1:L]
	double free_energy ;

	// make_digamma_table() ;
	
	if(opt.model == VBAM0 || opt.model == VBAM1) {
	    
	    /* admixture model */
	    
	    const double *lambda0 = initialise_lambda_prior_vbam(K, n) ;
	    double *lambda1 = ALLOC(n * K, double) ;
	
	    if(opt.model == VBAM0) {
		bool haploid = FALSE ;
		vbam0(x, &n, &opt.L, &K, &haploid, lambda0, lambda1, alpha0, alpha1, &free_energy) ;
	    }
	    else if(opt.model == VBAM1)
		free_energy = vbam1(x, px, n, opt.L, K, lambda0, lambda1, alpha0, alpha1) ;
	    else ERROR("can't get here") ;

	    write_current_state_vbam(lambda1, alpha1, free_energy, -1, K, n, opt.L) ;
	
	    free(lambda0) ;
	    free(lambda1) ;
	    free(alpha0) ;
	    free(alpha1) ;
	}
	else if(opt.model == VBNAM0 || opt.model == VBNAM1) {
	    
	    /* no admixture model */
	    
	    const double *lambda0 = initialise_lambda_prior_vbnam(K, n) ;
	    double *lambda1 = ALLOC(K, double) ;
	    double *pz = ALLOC(n * K, double) ; // pz[i,1:K] is a prob. dist. over the missing indicator for individual i
	    
	    if(opt.model == VBNAM0) {
		bool haploid = FALSE ;
		vbnam0(x, &n, &opt.L, &K, &haploid, lambda0, lambda1, alpha0, alpha1, pz, 
		       &free_energy, &opt.stop_tol, &opt.niters) ;
	    }
	    else
		free_energy = vbnam1(x, px, n, opt.L, K, lambda0, lambda1, alpha0, alpha1, pz) ;
	    
	    write_current_state_vbnam(pz, lambda1, alpha1, free_energy, -1, K, n, opt.L) ;
	
	    fclose(opt.timefile) ;
	    fclose(opt.statefile) ;
	    fclose(opt.wstatefile) ;
	    free(lambda0) ;
	    free(lambda1) ;
	    free(alpha0) ;
	    free(alpha1) ;
	    free(pz) ;
	}
    }
    else ERROR("can't get here") ;
    
#if PARALLEL
    MPI_Finalize() ; 
#endif
    if(thisproc == ROOT) free(first_snps) ;
    if(opt.uncertain_data) free(px) ;
    else free(x) ;
    free(snp_include) ;
    free(opt.snploglike) ;
    free(opt.fit) ;
    free(opt.updatez) ;
    if(opt.popfile != NULL) free(opt.fixpop) ;
    
    return 0 ;
}

void set_first_snps(int L, int nprocs, int *first_snps) {
    int i, proc, chunk = L / nprocs ;
    
    for(i = 0, proc = 0 ; i < L ; i++)
	if( i % chunk == 0 ) first_snps[proc++] = i ;
}

int get_allele(const unsigned char **x, int a) {
    int geno = **x - '0' ;
    if(a == 1) (*x)++ ;
    
    assert(geno == 0 || geno == 1 || geno == 2) ;
    assert(geno != MISSING) ;
    return (geno == 2) || (geno == 1 && a == 1) ; // have allele 1 if homozygote or second allele at heterozygote
}

bool reached_plateau(int iter, double *fit) {
    double rel_inc, begin, end ;
    const double SMALLVAL = 1e-9 ;
    
    begin = fit[iter - opt.plateau_size] ;
    end = fit[iter] ;

    rel_inc = (end - begin) / -begin ;

    if(!(rel_inc > 0 || fabs(rel_inc) < SMALLVAL || iter == 0)) {
	write_matrix_double(fit + iter - opt.plateau_size, stdout, 1, opt.plateau_size, "%lf") ;
	ERROR("reached_plateau(): begin = %.12lf, end = %.12lf, rel_inc = %.12lf:", begin, end, rel_inc) ;
    }
    
    return rel_inc < opt.stop_tol ;
}

double *rprior_mu(int L, int K) {
    double *mu, *alpha, tmp[NALS] ;
    int a, l, k, lk ;

    mu = ALLOC(K * L, double) ;
    
    alpha = ALLOC(NALS, double) ;
    for(a = 0 ; a < NALS ; a++)
	alpha[a] = 1.0 / NALS ;
    
    for(l = lk = 0 ; l < L ; l++)
	for(k = 0 ; k < K ; k++, lk++) {
	    rdirichlet(alpha, NALS, tmp) ;
	    mu[lk] = tmp[1] ; // diallelism is hard-coded at the moment
	    if(FIXED_INIT) mu[lk] = 0.5 ; // fixed initialisation
	}
    
    free(alpha) ;
    return mu ;
}

double *rprior_q(int n, int K) {
    double *q, *lambda ;
    int i, k ;

    q = ALLOC(n * K, double) ;
    
    lambda = ALLOC(K, double) ;
    for(k = 0 ; k < K ; k++) {
	lambda[k] = 1 ; // 2008-06-02
	// lambda[k] = 1.0 / K ;
    }

    for(i = 0 ; i < n ; i++) {
	rdirichlet(lambda, K, q + (i * K)) ;
	if(FIXED_INIT) for(k = 0 ; k < K ; k++)
	    q[k + (i*K)] = 1.0 / K ; // fixed initialisation
    }
    
    free(lambda) ;
    return q ;
}

void compute_freqs_and_snploglike_nostructure(unsigned char *x, double *p, double *loglike) {
    int l, i, li, x_li, nA_l, twon_l ;

    for(l = 0, li = 0 ; l < opt.L ; l++) {
	nA_l = twon_l = 0 ;
	for(i = 0 ; i < opt.n ; i++, li++) {
	    x_li = x[li] - '0' ;
	    if(x_li != MISSING) {
		nA_l += x_li ;
		twon_l += 2 ;
	    }
	}
	p[l] = nA_l / (double) twon_l ;

	if(nA_l == 0 || nA_l == twon_l) loglike[l] = 0 ;
	else loglike[l] = 
		 nA_l            * log(    p[l]) + 
		 (twon_l - nA_l) * log(1 - p[l]) ;
    }
}

char *model_name() {
    switch(opt.model) {
    case EMAM: return "emam" ;
    case MCMCAM: return "mcmcam" ;
    case VBAM0: return "vbam0" ;
    case VBAM1: return "vbam1" ;
    case VBNAM0: return "vbnam0" ;
    case VBNAM1: return "vbnam1" ;
    }
    ERROR("can't get here") ;
    return NULL ;
}
