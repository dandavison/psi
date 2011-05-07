#include "psi.h"
#include "vb.h"


double *initialise_lambda_prior_vbnam(int K, int n) {
    int k ;
    double *lambda0 = ALLOC(K, double) ;

    for(k = 0 ; k < K ; k++)
	// lambda0[k] = n / (double) K ; // runif(.9,1.1) ;
	lambda0[k] = 1 / (double) K ;

    return lambda0 ;
}

double vbnam1(const unsigned char *x, double *px, int n, int L, int K, const double *lambda0, 
	      double *lambda1, const double *alpha0, double *alpha1, double *pz) {
    
    FILE *f = fopen("alphavals", "w") ;

    const unsigned char *x0 = x ;
    int i, j, a, k, l, iter, x_lj, geno ;
    double *pz_i, *logpz, *logpz_i, *zeros, *enla, *alpha1_l ;
    const double *alpha0_l ;
    double entropy, d_post_prior, E_log_pz, E_log_like, sum_lambda1, tot, tmp, fac ;
    
    logpz = pz ;
    zeros = CALLOC(n * K * NALS, double) ;
    enla = ALLOC(K, double) ;
    
    /* Initialise posteriors equal to priors */
    memcpy(alpha1, alpha0, L * K * NALS * sizeof(double)) ;
    memcpy(lambda1, lambda0, K * sizeof(double)) ;
    
    sum_lambda1 = 0 ;
    for(k = 0 ; k < K ; k++) sum_lambda1 += lambda1[k] ;

    iter = 0 ; 
    do {
	memcpy(enla, zeros, K * sizeof(double)) ;
	memcpy(logpz, zeros, n * K * sizeof(double)) ;
	d_post_prior = entropy = E_log_pz = E_log_like = 0.0 ; 
	
	d_post_prior += KL_Dirichlets(lambda1, lambda0, K) ;

	/* E STEP */
	x = x0 ;
	for(l = 0 ; l < L ; l++) {
	    alpha1_l = alpha1 + l*K*NALS ;
	    alpha0_l = alpha0 + l*K*NALS ;
	    
	    for(k = 0 ; k < K ; k++)
		// 2008-05-25 *NALS was ommitted
		d_post_prior += KL_Dirichlets(alpha1_l + k*NALS, alpha0_l + k*NALS, NALS) ;
	    
	    for(j = 0 ; j < 2*n ; j++) {
		i = j / 2 ;
		a = j % 2 ;
		
		// x_lj = get_allele(&x, a) ;
		geno = *x - '0' ;
		x_lj = (geno == 2) || (geno == 1 && a == 1) ;
		if(a == 1) x++ ;

		logpz_i = logpz + i*K ;
		
		for(k = 0 ; k < K ; k++) {
		    
		    // fprintf(f, "%lf\n", alpha1_l[k*NALS + x_lj]) ;
		    tmp = DIGAMMA(alpha1_l[k*NALS + x_lj]) - DIGAMMA( alpha1_l[k*NALS + 0] + alpha1_l[k*NALS + 1] ) ;
		    
		    // if( isnan(tmp) || tmp > 0 ) ERROR("tmp = %lf: l = %d, i = %d, k = %d\n", tmp, l, i, k) ;
		    
		    logpz_i[k] += tmp ;
		    /* if( isnan(logpz_i[k]) || logpz_i[k] > 0 ) 
		       ERROR("logpz[i=%d,k=%d] = %lf, tmp = %lf", i, k, logpz_i[k], tmp) ; */
		}
	    }
	}

	for(i = 0 ; i < n ; i++) {
	    fac = logpz[i*K + 0] + DIGAMMA(lambda1[0]) - DIGAMMA(sum_lambda1) ;
	    for(k = 0 ; k < K ; k++) {
		logpz[i*K + k] += DIGAMMA(lambda1[k]) - DIGAMMA(sum_lambda1) ;
		fac = MAX(fac, logpz[i*K + k]) ;
	    }
	    tot = 0 ;
	    for(k = 0 ; k < K ; k++) {
		pz[i*K + k] = exp( logpz[i*K + k] - fac ) ;
		tot += pz[i*K + k] ;
	    }
	    for(k = 0 ; k < K ; k++) {
		pz[i*K + k] /= tot ;
		// if(isnan(pz[i*K + k])) ERROR("pz[i=%d,k=%d] = %lf, tot = %lf", i, k, pz[i*K + k], tot) ;
		enla[k] += pz[i*K + k] ;
		entropy -= pz[i*K + k] * log( pz[i*K + k] ) ;
	    }
	}
	x = x0 ;
	opt.fit[iter] = free_energy(x, pz, alpha0, alpha1, lambda0, lambda1, n, K, L, iter) ;
	if(opt.obey_stop_tol && iter >= opt.plateau_size && reached_plateau(iter, opt.fit)) break ;
	
	/* M STEP */
	x = x0 ;
	memcpy(alpha1, alpha0, L*K*NALS*sizeof(double)) ; // 2008-05-22
	for(l = 0 ; l < L ; l++) {
	    alpha1_l = alpha1 + l*K*NALS ;
	    // memcpy(alpha1_l, alpha0, K*NALS*sizeof(double)) ; 2008-05-22
	    for(j = 0 ; j < 2*n ; j++) {
		i = j / 2 ;
		a = j % 2 ;
		
		x_lj = get_allele(&x, a) ;
		pz_i = pz + i * K ;

		for(k = 0 ; k < K ; k++)
		    alpha1_l[k*NALS + x_lj] += pz_i[k] ;
	    }
	}
	
	/* insert free energy computation here if you want to do it after M step */

	sum_lambda1 = 0 ;
	memcpy(lambda1, lambda0, K * sizeof(double)) ;
	for(k = 0 ; k < K ; k++) {
	    // lambda1[k] += enla[k] ; 2008-05-26: not updating lambda, as in Pritchard et al. 2000
	    sum_lambda1 += lambda1[k] ;
	}
    }
    while(++iter < opt.niters) ;
    
    free(zeros) ;
    free(enla) ;

    fclose(f) ;

    return opt.fit[iter - 1] ;
}

/**************************************************************************************************/

double free_energy(const unsigned char *x, double *pz, const double *alpha0, double *alpha1,
		   const double *lambda0, double *lambda1, int n, int K, int L, int iter) {

    int l, i, j, k, c, x_lj, geno ;
    double F, rel_inc, entropy, d_post_prior, E_log_pz, E_log_like, sum_lambda, digamma_sum_alpha ;
    double *enlka_l, *enk, *alpha1_l, *pz_i ;
    const double *alpha0_l ;

    enk = CALLOC(K, double) ;
    enlka_l = CALLOC(K*NALS, double) ;
    
    d_post_prior = E_log_like = 0.0 ; 
    
    d_post_prior += KL_Dirichlets(lambda1, lambda0, K) ;

    for(l = 0 ; l < L ; l++) {
	alpha1_l = alpha1 + l*K*NALS ;
	alpha0_l = alpha0 + l*K*NALS ;

	for(j = 0 ; j < 2*n ; j++) {
	    i = j / 2 ;
	    c = j % 2 ;
	    
	    // x_lj = get_allele(&x, c) ;
	    geno = *x - '0' ;
	    x_lj = (geno == 2) || (geno == 1 && c == 1) ;
	    if(c == 1) x++ ;

	    pz_i = pz + i * K ;

	    for(k = 0 ; k < K ; k++)
		enlka_l[k*NALS + x_lj] += pz_i[k] ;
	}
	for(k = 0 ; k < K ; k++) {
	    digamma_sum_alpha = DIGAMMA( alpha1_l[k*NALS + 0] + alpha1_l[k*NALS + 1] ) ;
	    
	    E_log_like += enlka_l[k*NALS + 0] * ( DIGAMMA(alpha1_l[k*NALS + 0]) - digamma_sum_alpha ) ;
	    E_log_like += enlka_l[k*NALS + 1] * ( DIGAMMA(alpha1_l[k*NALS + 1]) - digamma_sum_alpha ) ;
	    
	    enlka_l[k*NALS + 0] = enlka_l[k*NALS + 1] = 0 ;
	    
	    d_post_prior += KL_Dirichlets(alpha1_l + k*NALS, alpha0_l + k*NALS, NALS) ;
   	}
    }
    entropy = 0 ;
    for(i = 0 ; i < n ; i++)
	for(k = 0 ; k < K ; k++) {
	    enk[k] += pz[i*K + k] ;
	    if(pz[i*K + k] > 0) // prob=0 doesn't contribute to entropy
		entropy -= log(pz[i*K + k]) * pz[i*K + k] ;
	}
    
    E_log_pz = sum_lambda = 0 ;
    for(k = 0 ; k < K ; k++) {
	E_log_pz += enk[k] * DIGAMMA(lambda1[k]) ;
	sum_lambda += lambda1[k] ;
    }
    E_log_pz -= n * DIGAMMA(sum_lambda) ;
    
    F = E_log_pz + E_log_like + entropy - d_post_prior ;

    rel_inc = iter > 0 ? (F - opt.fit[iter - 1]) / (-opt.fit[iter - 1]) : log(-1) ;

    if(opt.show_output && !(iter % opt.print_period) && thisproc == ROOT)
	PRINT("%5d %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf\n", 
	      iter, d_post_prior , E_log_pz , E_log_like , entropy, F, 100 * rel_inc) ;
    
    if(isnan(F)) {
	PRINT("pz:\n") ;
	write_matrix_double(pz, stdout, K, n, "%8.3lf") ;
	ERROR("nan generated") ;
    }

    return F ;
}


void write_current_state_vbnam(double *pz, double *lambda, double *alpha, double F, int iter, int K, int n, int L) {
    FILE *f ;
    bool final = (iter < 0) ;
    

    if(thisproc == ROOT) {
	if(final) sprintf(opt.fnamebuff, "%s/lambda", opt.outdir) ;
	else sprintf(opt.fnamebuff, "%s/lambda-%05d", opt.outdir, iter) ;
	f = fopen(opt.fnamebuff, "w") ;
	write_matrix_double(lambda, f, K, 1, "%-10.5lf") ;
	fclose(f) ;
	
	if(final) sprintf(opt.fnamebuff, "%s/F", opt.outdir) ;
	else sprintf(opt.fnamebuff, "%s/F-%05d", opt.outdir, iter) ;
	f = fopen(opt.fnamebuff, "w") ;
	fprintf(f, "%lf\n", F) ;
	fclose(f) ;
    
	if(final) sprintf(opt.fnamebuff, "%s/pz", opt.outdir) ;
	else sprintf(opt.fnamebuff, "%s/pz-%05d", opt.outdir, iter) ;
	f = fopen(opt.fnamebuff, "w") ;
	write_matrix_double(pz, f, K, n, "%-7.4lf") ;
	fclose(f) ;
    }
    
#if PARALLEL
    /* Every process writes out allele frequencies at its own set of SNPs */
    if(final) sprintf(opt.fnamebuff, "%s/alpha-%02d", opt.outdir, thisproc) ;
    else sprintf(opt.fnamebuff, "%s/alpha-%05d-%02d", opt.outdir, iter, thisproc) ;
#else
    if(final) sprintf(opt.fnamebuff, "%s/alpha", opt.outdir) ;
    else sprintf(opt.fnamebuff, "%s/alpha-%05d", opt.outdir, iter) ;
#endif
    f = fopen(opt.fnamebuff, "w") ;
    write_matrix_double(alpha, f, K*NALS, L, "%-8.4lf ") ;
    fclose(f) ;
}

void vbnam0(const unsigned char *x, int *_n, int *_L, int *_K, bool *_haploid,
	    const double *lambda_prior, double *lambda,
	    const double *alpha_prior, double *alpha,
	    double *pz,
	    double *_free_energy, 
	    double *stop_tol, 
	    int *niters) {

    const int 
	n = *_n,                             // number of individuals
	twon = 2*n,                          // number of allele copies at a locus
	L = *_L,                             // number of loci
	K = *_K ;                             // number of components in mixture
    const bool haploid = *_haploid ;
    double
	*m = CALLOC(K, double),           // m_ik     = E #{individuals that derive from population k}
	*mm_l = CALLOC(NALS * K, double),    // mm_l[1:NALS,1:K] = E #{alleles of type j at locus l that derive from population k}
	*zeros = CALLOC(n*K, double) ;
    int l, j, i, pair, a, k, x_lj, geno, iter=0 ;
    double 
	*alpha_l, *alpha_lk,
	*mm_lk,
	*logpz = pz,
	digamma_sum_lambda,
	entropy, E_log_like, E_log_p_z, d_KL,
	tot, fac ;
    const double *alpha_prior_l, *alpha_prior_lk ;
    const unsigned char *x_li ;
    
    if(VERBOSE) PRINT("structure VB no-admixture version 0.0: n=%d, L=%d, K=%d, %s, NALS=%d, tol=%lf, niters=%d\n", 
		      n, L, K, haploid ? "haploid" : "diploid", NALS, opt.stop_tol, *niters) ;    
    
    assert(alpha_prior != NULL) ;
    memcpy(alpha, alpha_prior, L * K * NALS * sizeof(double)) ;
    
    assert(lambda_prior != NULL) ;
    memcpy(lambda, lambda_prior, K * sizeof(double)) ; 
    
    tot = 0.0 ;
    for(k = 0 ; k < K ; k++) tot += lambda[k] ;
    digamma_sum_lambda = DIGAMMA(tot) ;
    
    do {
	x_li = x ; pair = 0 ; 
	memcpy(logpz, zeros, n * K * sizeof(double)) ;  // 2008-04-24
	entropy = d_KL = E_log_like = E_log_p_z = 0 ;

	for(l = 0 ; l < L ; ++l) {
	    
	    alpha_l = alpha + (l * NALS * K) ; alpha_prior_l = alpha_prior + (l * NALS * K) ;
	    
	    for(j = 0 ; j < twon ; ++j) {
		
		i = j / 2 ;                            // j indexes chromosomes; i indexes individuals
		// m_i = m + i * K ;
		// lambda_i = lambda + i * K ;
		// digamma_sum_lambda_i = digamma_sum_lambda + i ;
		
		a = j % 2 ;                            // a == 0 if first allele, a == 1 if second allele
		
		/**************************************************************** 
		 WORK OUT WHAT ALLELE IS PRESENT ON THIS CHROMOSOME AT THIS LOCUS 
		****************************************************************/

		//x_lj = get_allele(&x_li, a, pair, haploid) ; assert(x_lj != MISSING) ;
		// x_lj = get_allele(&x_li, a) ; assert(x_lj != MISSING) ;
		geno = *x_li - '0' ;
		x_lj = (geno == 2) || (geno == 1 && a == 1) ;
		if(a == 1) x_li++ ;

		/***************************************
		              E-STEP
                 INCREMENT THE APPROXIMATE LOG-POSTERIOR OVER THE 
                 MISSING INDICATORS FOR THIS INDIVIDUAL
		****************************************/
		
		if(x_lj != MISSING)
		    for(k = 0, alpha_lk = alpha_l ; k < K ; ++k, alpha_lk += NALS)
			logpz[i*K + k] += DIGAMMA(alpha_lk[x_lj]) - DIGAMMA(alpha_lk[0] + alpha_lk[1]); // 2008-05-22 E{likelihood given missing indicators}
		//logpz[i + k*n] += DIGAMMA(alpha_lk[x_lj]) - DIGAMMA(alpha_lk[0] + alpha_lk[1]); // E{likelihood given missing indicators}
	    }
	}

	/***************************************
		        E-STEP
	  NORMALISE PROBABILITY DISTRIBUTIONS 
          OVER MISSING INDICATORS FOR EACH 
                      INDIVIDUAL
	****************************************/

	for(k = 0 ; k < K ; ++k) m[k] = 0.0 ;
	for(i = 0 ; i < n ; ++i) {
	    fac = logpz[i*K + 0] + DIGAMMA(lambda[0]) - digamma_sum_lambda ;
	    for(k = 0 ; k < K ; ++k) {
		logpz[i*K + k] += DIGAMMA(lambda[k]) - digamma_sum_lambda ; // 2008-05-25 2008-05-22
		fac = MAX(fac, logpz[i*K + k]) ;
	    }
	    tot = 0.0 ;
	    for(k = 0 ; k < K ; ++k) {
		pz[i*K + k] = exp(logpz[i*K + k] - fac) ; // 2008-05-25
		tot += pz[i*K + k] ;
	    }
	    for(k = 0 ; k < K ; ++k) {
		pz[i*K + k] /= tot ; // 2008-05-25
		m[k] += pz[i*K + k] ;
		if(pz[i*K + k] > 0) entropy -= pz[i*K + k] * log(pz[i*K + k]) ;
	    }
	}
	
	
	/****************************************************************************************
		                                M-STEP
         USE APPROXIMATE POSTERIOR ON MISSING INDICATORS TO INCREMENT UPDATES TO HYPERPARAMETERS
	****************************************************************************************/
	

	//PRINT("%5d ", iter) ;
	//x_tmp = x ;
	//free_energy(x_tmp, pz, alpha_prior, alpha, lambda_prior, lambda, n, K, L, iter) ;

	
	/* Update lambda */

	//PRINT("enla:\n") ;
	//write_matrix_double(m, stdout, K, 1, "%.3lf ") ;

	d_KL += KL_Dirichlets(lambda, lambda_prior, K) ;
	E_log_p_z -= n * digamma_sum_lambda ; // 2008-04-23
	
	tot = 0.0 ;
	for(k = 0 ; k < K ; ++k) {
	    E_log_p_z += m[k] * DIGAMMA(lambda[k]) ;
	    tot += lambda[k] = lambda_prior[k] ; // + m[k] ; not updating lambda
	}
	digamma_sum_lambda = DIGAMMA(tot) ;

	/* Update mu */
	x_li = x ; pair = 0 ; 
	for(l = 0 ; l < L ; ++l) {
	    
	    alpha_l = alpha + (l * NALS * K) ; alpha_prior_l = alpha_prior + (l * NALS * K) ;
	    
	    for(j = 0 ; j < twon ; ++j) {
		
		i = j / 2 ;                            // j indexes chromosomes; i indexes individuals
		a = j % 2 ;                            // a == 0 if first allele, a == 1 if second allele
		
		/**************************************************************** 
		 WORK OUT WHAT ALLELE IS PRESENT ON THIS CHROMOSOME AT THIS LOCUS 
		****************************************************************/
		
		//x_lj = get_allele(&x_li, a, pair, haploid) ; assert(x_lj != MISSING) ;
		// x_lj = get_allele(&x_li, a) ; assert(x_lj != MISSING) ;
		geno = *x_li - '0' ;
		x_lj = (geno == 2) || (geno == 1 && a == 1) ;
		if(a == 1) x_li++ ;


		for(k = 0 ; k < K ; ++k) mm_l[x_lj + k*NALS] += pz[i*K + k] ; // 2008-05-25
		//for(k = 0 ; k < K ; ++k) mm_l[x_lj + k*NALS] += pz[i + k*n] ;
	    }

	    for(alpha_lk = alpha_l, alpha_prior_lk = alpha_prior_l, mm_lk = mm_l, 
		    k = 0 ; k < K ; ++k, 
		    // alpha_lk += NALS, mm_lk += NALS) { 
		    alpha_lk += NALS, alpha_prior_lk += NALS, mm_lk += NALS) { // 2008-04-24
				
		d_KL += KL_Dirichlets(alpha_lk, alpha_prior_lk, NALS) ;
		E_log_like += 
		    mm_lk[0] * DIGAMMA(alpha_lk[0]) + mm_lk[1] * DIGAMMA(alpha_lk[1]) -
		    (mm_lk[0] + mm_lk[1]) * DIGAMMA(alpha_lk[0] + alpha_lk[1]) ;
		
		alpha_lk[0] = mm_lk[0] + alpha_prior_lk[0] ; // assumes
		alpha_lk[1] = mm_lk[1] + alpha_prior_lk[1] ; // biallelic
		
		mm_lk[0] = mm_lk[1] = 0.0 ;
	    }
	}

	opt.fit[iter] = -d_KL + E_log_p_z + E_log_like + entropy ;

	if(opt.show_output && !(iter % opt.print_period) && thisproc == ROOT)
	    PRINT("%5d %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf\n", 
		  iter, d_KL , E_log_p_z , E_log_like , entropy, opt.fit[iter], 
		  iter > 0 ? (100 * (opt.fit[iter] - opt.fit[iter - 1]) / -opt.fit[iter - 1]) : log(-1)) ;

	if(opt.obey_stop_tol && iter >= opt.plateau_size && reached_plateau(iter, opt.fit)) break ;
	// free_energy(x_li, pz, alpha_prior, alpha, lambda_prior, lambda, n, K, L, iter) ;
    }
    while(++iter < *niters) ; 
    
    *_free_energy = opt.fit[iter - 1] ;
}
