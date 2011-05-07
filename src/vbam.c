#include "psi.h"
#include "vb.h"

double *initialise_lambda_prior_vbam(int K, int n) {
    int i, k ;
    double *lambda0 = ALLOC(n * K, double) ;
	
    for(i = 0 ; i < n ; i++)
	for(k = 0 ; k < K ; k++)
	    lambda0[i*K + k] = 1 ; //opt.q_prior_strength * 2 * opt.L / (double) K ;
    
    return lambda0 ;
}

double vbam1(const unsigned char *x, double *px, int n, int L, int K, 
	     const double *lambda0, double *lambda1, const double *alpha0, double *alpha1) {
    
    const unsigned char *x0 = x ;
    int i, j, a, k, l, ik, k0, k1, iter, x_lj, geno ;
    double *pz_lj, *logpz_lj, *enla, *enia_l, *alpha1_l, *sum_lambda1, *digamma_lambda1_minus_sum ;
    const double *alpha0_l, *zeros ;
    double entropy, d_post_prior, E_log_pz, E_log_like, tot, tmp, fac ;
    
    if(VERBOSE) PRINT("vbam1: n=%d, L=%d, K=%d\n", n, L, K) ;

    pz_lj = ALLOC(K, double) ;
    logpz_lj = pz_lj ;
    enla = ALLOC(n * K, double) ; // E #{allele copies in individual i that derive from population k}
    enia_l = ALLOC(NALS * K, double) ;// E #{alleles of type j at locus l that derive from population k}
    sum_lambda1 = ALLOC(n, double) ;
    zeros = CALLOC(n * K * NALS, double) ;
    digamma_lambda1_minus_sum = ALLOC(n*K, double) ;
    
    /* Initialise posteriors equal to priors */
    memcpy(alpha1, alpha0, L * K * NALS * sizeof(double)) ;
    memcpy(lambda1, lambda0, n * K * sizeof(double)) ;
    
    iter = 0 ; 
    do {
	
	d_post_prior = entropy = E_log_pz = E_log_like = 0.0 ; 
	memcpy(enla, zeros, n * K * sizeof(double)) ;
	x = x0 ;

	for(i = 0 ; i < n ; i++)
	    d_post_prior += KL_Dirichlets(lambda1 + i*K, lambda0 + i*K, K) ;

	for(i = 0 ; i < n ; i++) {
	    sum_lambda1[i] = 0 ;
	    for(k = 0 ; k < K ; k++) sum_lambda1[i] += lambda1[i*K + k] ;
	}

	/*
	  efficiency:
	  digamma(lambda) stuff doesn't need to be recomputed at every locus
	*/
	for(i = 0 ; i < n ; i++)
	    for(k = 0 ; k < K ; k++)
		digamma_lambda1_minus_sum[i*K + k] = DIGAMMA(lambda1[i*K + k]) - DIGAMMA(sum_lambda1[i]) ;
	
	for(l = 0 ; l < L ; l++) {
	    alpha1_l = alpha1 + l*K*NALS ;
	    alpha0_l = alpha0 + l*K*NALS ;
	    memcpy(enia_l, zeros, NALS * K * sizeof(double)) ;
	    
	    for(k = 0 ; k < K ; k++)
		d_post_prior += KL_Dirichlets(alpha1_l + k*NALS, alpha0_l + k*NALS, NALS) ;

	    /* E STEP */

	    for(j = 0 ; j < 2*n ; j++) {
		i = j / 2 ;
		a = j % 2 ;
		
		// x_lj = get_allele(&x, a) ;
		geno = *x - '0' ;
		x_lj = (geno == 2) || (geno == 1 && a == 1) ;
		if(a == 1) x++ ;
		
		tot = 0 ;
		fac = log(0) ;
		for(k = 0 ; k < K ; k++) {
		  
		    // logpz_lj[k] = DIGAMMA(lambda1[i*K + k]) - DIGAMMA(sum_lambda1[i]) ;
		    logpz_lj[k] = digamma_lambda1_minus_sum[i*K + k] ;

		    tmp = DIGAMMA(alpha1_l[NALS * k + x_lj]) - DIGAMMA( alpha1_l[NALS * k + 0] + alpha1_l[NALS * k + 1] ) ;
		    // if( isnan(tmp) || tmp > 0 ) ERROR("tmp = %lf: l = %d, i = %d, k = %d\n", tmp, l, i, k) ;
		    logpz_lj[k] += tmp ;
		    // if( isnan(logpz_lj[k]) || logpz_lj[k] > 0 )
		    // 	ERROR("logpz[l=%d,j=%d,k=%d] = %lf, tmp = %lf", l, j, k, logpz_lj[k], tmp) ;
		    
		    fac = MAX(fac, logpz_lj[k]) ;
		}
		
		/* Normalise pz_lj & increment expected counts in preparation for M STEP */
		
		for(k = 0 ; k < K ; k++)
		    tot += pz_lj[k] = exp(logpz_lj[k] - fac) ;

		for(k = 0 ; k < K ; k++) {
		    pz_lj[k] /= tot ;
		    if(pz_lj[k]) {
			enla[i*K + k] += pz_lj[k] ;
			enia_l[k*NALS + x_lj] += pz_lj[k] ;
			entropy -= log(pz_lj[k]) * pz_lj[k] ;
		    }
		}
	    }
	    
	    /* M STEP */
	    
	    /* Done all alleles at this locus; now use enia and alpha1 to increment E_log_like,
	       before updating alpha at this locus */

	    for(k = 0 ; k < K ; k++) {
		k0 = k*NALS + 0 ;
		k1 = k*NALS + 1 ;
		
		E_log_like += 
		    enia_l[k0] * DIGAMMA(alpha1_l[k0])   + 
		    enia_l[k1] * DIGAMMA(alpha1_l[k1])   -
		    (enia_l[k0] + enia_l[k1]) * DIGAMMA(alpha1_l[k0] + alpha1_l[k1]) ;
		
		alpha1_l[k0] = alpha0_l[k0] + enia_l[k0] ;
		alpha1_l[k1] = alpha0_l[k1] + enia_l[k1] ;
	    }
	}
	
	/* Done all loci for this iteration; now update lambda and
	   compute remaining components of free energy term */
	
	for(i = 0, ik = 0 ; i < n ; i++) {
	    E_log_pz -= 2 * L * DIGAMMA(sum_lambda1[i]) ;
	    
	    for(k = 0 ; k < K ; k++, ik++) {
		E_log_pz += enla[ik] * DIGAMMA(lambda1[ik]) ;
		lambda1[ik] = lambda0[ik] + enla[ik] ;
	    }
	}
	
	opt.fit[iter] = E_log_pz + E_log_like + entropy - d_post_prior ;

	if(opt.show_output && !(iter % opt.print_period) && thisproc == ROOT)
	    PRINT("%5d %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf\n", 
		  iter, d_post_prior , E_log_pz , E_log_like , entropy, opt.fit[iter], 
		  iter > 0 ? (100 * (opt.fit[iter] - opt.fit[iter - 1]) / -opt.fit[iter - 1]) : log(-1)) ;

	write_current_state_vbam(lambda1, NULL, -1, iter, K, n, L) ;
	if(opt.obey_stop_tol && iter >= opt.plateau_size && reached_plateau(iter, opt.fit)) break ;
    }
    while(++iter < opt.niters) ;
    
    PRINT("%d iters\n", iter+1) ;
    free(pz_lj) ;
    free(enla) ;
    free(enia_l) ;
    free(sum_lambda1) ;
    free(digamma_lambda1_minus_sum) ;
    free(zeros) ;
    
    return opt.fit[iter - 1] ;
}

void write_current_state_vbam(double *lambda, double *alpha, double F, int iter, int K, int n, int L) {
    FILE *f ;
    bool final = (iter < 0) ;
    struct timeval tv ;
    double now ;


    if(thisproc == ROOT) {
	if(final) {
	    sprintf(opt.fnamebuff, "%s/lambda", opt.outdir) ;
	    // else sprintf(opt.fnamebuff, "%s/lambda-%05d", opt.outdir, iter) ;
	    f = fopen(opt.fnamebuff, "w") ;
	    write_matrix_double(lambda, f, K, n, "%-10.5lf ") ;
	    fclose(f) ;
	    
	    sprintf(opt.fnamebuff, "%s/F", opt.outdir) ;
	    // else sprintf(opt.fnamebuff, "%s/F-%05d", opt.outdir, iter) ;
	    f = fopen(opt.fnamebuff, "w") ;
	    fprintf(f, "%lf\n", F) ;
	    fclose(f) ;
	}
	else {
	    gettimeofday(&tv, NULL) ;
	    now = (double) tv.tv_sec + 1e-6 * (double) tv.tv_usec;
	    fprintf(opt.timefile, "%lf\n", now) ;
	    write_matrix_double(lambda, opt.statefile, K*n, 1, "%-10.5lf ") ;
	}
    }

    if(final) {
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
}


/* 
       Fit structure admixture model via variational Bayes. 
       version 0.0

       Original version, written to be called via .C() in R

       Dan Davison
        
       The code can be compiled differently according to the way in which the genotype data is stored:

       #if BIT_ENCODING:
       Genotypes \in {NA,AA,Aa,aa} are encoded using 2 bits, so 1 byte (e.g. a char) contains 4 genotypes.
       The encoding is (NA,AA,Aa,aa) <--> (00,01,10,11)

       #else:
       Genotypes \in {NA,AA,Aa,aa} are encoded using single bytes.
       1 byte = 8 bits = 2^8 possibilities = 16^2 possibilities, so they can be writen as two hexadecimal numbers, e.g. 3a, 03
       The encoding is (NA,AA,Aa,aa) <--> (00,01,02,03).

       Pointers with no underscore subscript (e.g. double *lambda) are constant placeholders for the beginning of the array; 
       the versions with subscripts are used to move through the array.

       This is intended to extend easily to multiallelic data (J > 2), but a couple of things assume biallelic, including:
           - There is no digamma_sum_alpha, instead this is formed by digamma(alpha[0] + alpha[1])
*/
    
void vbam0(unsigned char *x, int *_n, int *_L, int *_K, bool *_haploid,
	   const double *lambda_prior, double *lambda,
	   const double *alpha_prior, double *alpha, 
	   double *_free_energy) {

    const int 
	n = *_n,                             // number of individuals
	twon = 2*n,                          // number of allele copies at a locus
	L = *_L,                             // number of loci
	J = 2,                               // number of alleles at a locus (this code for biallelic SNP data)
	K = *_K ;                             // number of components in mixture
    bool haploid = *_haploid ;
    int l, j, i, pair, a, k, x_lj, iter=0, ploidy = haploid ? 1 : 2 ;
    double 
	*alpha_l, *alpha_lk,
	*lambda_i, *lambda_ik, 
	*m, *m_i, *m_ik,                     // m_ik     = E #{allele copies in individual i that derive from population k}
	*mm_l, *mm_lk,                       // mm_lk[j] = E #{alleles of type j at locus l that derive from population k}
	*digamma_sum_lambda, *digamma_sum_lambda_i, 
	*p_lj, *p_ljk,                          // p_lj[1:K] is a probability distribution over the missing indicator at a single allele copy
	entropy, E_log_like, E_log_p_z, d_KL,
	tmp, tot ;
    const double *alpha_prior_l, *alpha_prior_lk, *lambda_prior_i ;
    unsigned char *x_li, geno ;
    
    if(VERBOSE) PRINT("structure VB version 0.0: n=%d, L=%d, K=%d, %s, J=%d, tol=%lf\n", n, L, K, haploid ? "haploid" : "diploid", J, opt.stop_tol) ;
    
#if BIT_ENCODING
    if(VERBOSE) PRINT("2 bits per genotype\n") ;
#else
    if(VERBOSE) PRINT("1 byte per genotype\n") ;
#endif
    
    // alpha = ALLOC(L * K * J, double) ;                         // [1:J, 1:K, 1:L]
    memcpy(alpha, alpha_prior, L * K * J * sizeof(double)) ;
    
    // lambda = ALLOC(n * K, double) ;                            // [1:K, 1:n]
    memcpy(lambda, lambda_prior, n * K * sizeof(double)) ; 

    m = CALLOC(n * K, double) ;
    mm_l = CALLOC(J * K, double) ;
    
    p_lj = ALLOC(K, double) ;
    digamma_sum_lambda = CALLOC(n, double) ; 

    // make_digamma_table() ;
    
    // if(FALSE && version > 0) compute_allele_frequencies(mu_anc, x, n, L) ;
    

    // bug fix start

/*     for(k = 0, F_k = F ; k < K ; ++k, ++F_k) { */
/* 	*F_k = (1.0 - *F_k) / *F_k ; */
/* 	for(i = 0 ; i < n ; ++i) digamma_sum_lambda[i] += lambda[k + i*K] ; */
/*     } */
    
    // for(k = 0, F_k = F ; k < K ; ++k, ++F_k) *F_k = (1.0 - *F_k) / *F_k ;
    for(i = 0 ; i < n ; ++i) {
	for(k = 0 ; k < K ; ++k) digamma_sum_lambda[i] += lambda[k + i*K] ;
	digamma_sum_lambda[i] = DIGAMMA(digamma_sum_lambda[i]) ;
    }
    // end bug fix

    do {
	entropy = 0.0 ; E_log_like = 0.0 ; d_KL = 0.0 ; x_li = x ; pair = 0 ; 
	
	for(l = 0 ; l < L ; ++l) {
	    
	    alpha_l = alpha + (l * J * K) ; alpha_prior_l = alpha_prior + (l * J * K) ;
	    
	    for(j = 0 ; j < twon ; ++j) {
		
		i = j / 2 ;                            // j indexes chromosomes; i indexes individuals
		m_i = m + i * K ;
		lambda_i = lambda + i * K ;
		digamma_sum_lambda_i = digamma_sum_lambda + i ;
		
		a = j % 2 ;                            // a == 0 if first allele, a == 1 if second allele
		
		/**************************************************************** 
		 WORK OUT WHAT ALLELE IS PRESENT ON THIS CHROMOSOME AT THIS LOCUS 
		 ****************************************************************/

#if BIT_ENCODING
		geno = *x_li & MASK[pair] ;
		if(geno ==  4 || geno == 16 || geno ==  64) geno = 1 ;
		else if(geno ==  8 || geno == 32 || geno == 128) geno = 2 ;
		else if(geno == 12 || geno == 48 || geno == 192) geno = 3 ;
		if(a && ++pair == 4) { pair = 0 ; ++x_li ; }
#else
		geno = *x_li ;
		if(a == 1) ++x_li ;
#endif
		if(a == 0) {                           
		    if(geno == '9') x_lj = MISSING ;      // missing
		    else if(geno == '0') x_lj = 0 ;       // AA
		    else if(geno == '1') x_lj = 0 ;       // Aa -- heterozygote is always (0,1)
		    else if(geno == '2') x_lj = 1 ;       // aa
		}
		else if(haploid) {
		    assert(geno != '1') ;         // haploids are coded as one or other homozygote
		    x_lj = MISSING ;            // haploidy implemented by always treating second allele as missing
		}
		else if(geno == '1') x_lj = 1 ;   // only need to change allele if heterozygote
		
		assert(x_lj != MISSING) ;
		

		/***************************************
		              E-STEP
                 FORM THE APPROXIMATE POSTERIOR OVER THE 
                 MISSING INDICATORS FOR THIS ALLELE COPY
		****************************************/
		
		for(alpha_lk = alpha_l, lambda_ik = lambda_i, p_ljk = p_lj, tot = 0.0,
			k = 0 ;  k < K ; ++k, 
			alpha_lk += J, ++lambda_ik, ++p_ljk) {
		    
		    tmp = DIGAMMA(*lambda_ik) - *digamma_sum_lambda_i ; // E{missing indicators prior} 
		    if(x_lj != MISSING)
			tmp += DIGAMMA(alpha_lk[x_lj]) - DIGAMMA(alpha_lk[0] + alpha_lk[1]) ;// E{ like | missing indicators }
		    
		    tot += *p_ljk = exp(tmp) ;
		}

		
		/***************************************************
		                     M-STEP
                 USE APPROXIMATE POSTERIOR ON MISSING INDICATOR AT
                 THIS ALLELE TO INCREMENT UPDATES TO HYPERPARAMETERS
		****************************************************/

		for(p_ljk = p_lj, m_ik = m_i, mm_lk = mm_l,
			k = 0 ; k < K ; ++k,
			++m_ik, ++p_ljk, mm_lk += J) {
		    
		    *p_ljk /= tot ;
		    *m_ik += *p_ljk ;                   // used to update lambda when we've gone through all the loci
		    if(x_lj != MISSING) mm_lk[x_lj] += *p_ljk ; // used to update alpha when we've finished with this locus
		    if(*p_ljk) entropy -= *p_ljk * log(*p_ljk) ;
		}
	    }
	    
	    /* Done all alleles at this locus; now update alpha at this locus and record the KL divergence */
	    
	    // start with 'prior counts', expected 'real counts' added later
	    // memcpy(alpha_l, alpha_prior_l, J * K * sizeof(double)) ; 
	    
	    for(alpha_lk = alpha_l, alpha_prior_lk = alpha_prior_l, mm_lk = mm_l,
		    k = 0 ; k < K ; ++k, 
		    alpha_lk += J, alpha_prior_lk += J, mm_lk += J) {
		
		// alpha_lk[0] += mm_lk[0] ;  // assumes biallelic
		// alpha_lk[1] += mm_lk[1] ;  // assumes biallelic
			
		// 2008-05-26 NB in vbam1 I record KL divergence and E_log_like at outset of E step
		// so that they correspond to prior on iteration 0
		// 2008-05-26 OK I've effected that change now by switching the order of things in this loop
		E_log_like += 
		    mm_lk[0] * DIGAMMA(alpha_lk[0])   + 
		    mm_lk[1] * DIGAMMA(alpha_lk[1])   -
		    (mm_lk[0] + mm_lk[1]) * DIGAMMA(alpha_lk[0] + alpha_lk[1]) ;
	
		d_KL += KL_Dirichlets(alpha_lk, alpha_prior_lk, J) ;

		alpha_lk[0] = alpha_prior_lk[0] + mm_lk[0] ;  // assumes biallelic
		alpha_lk[1] = alpha_prior_lk[1] + mm_lk[1] ;  // assumes biallelic
		mm_lk[0] = mm_lk[1] = 0.0 ;
	    }
	}
	
	/* Done all loci for this iteration; now update lambda and compute remaining components of free energy term */
	
	// start with 'prior counts', expected 'real counts' are added during loop
	// memcpy(lambda, lambda_prior, n * K * sizeof(double)) ;
	
	for(lambda_i = lambda, lambda_prior_i = lambda_prior, 
		m_ik = m, digamma_sum_lambda_i = digamma_sum_lambda, E_log_p_z = 0.0,
		i = 0 ; i < n ; ++i, 
		++digamma_sum_lambda_i, lambda_i += K, lambda_prior_i += K) {
	    
	    d_KL += KL_Dirichlets(lambda_i, lambda_prior_i, K) ;
	    E_log_p_z -= ploidy * L * *digamma_sum_lambda_i ; // 2008-05-30 added factor of K

	    for(lambda_ik = lambda_i, tot = 0.0,
		    k = 0 ; k < K ; ++k, 
		    ++lambda_ik, ++m_ik) { 
		
		// *lambda_ik += *m_ik ;    // add expected 'real data counts' to 'prior counts'
		E_log_p_z += *m_ik * DIGAMMA(*lambda_ik) ; 
		*lambda_ik = lambda_prior_i[k] + *m_ik ;    // add expected 'real data counts' to 'prior counts'
		tot += *lambda_ik ;
		*m_ik = 0.0 ;
	    }

	    *digamma_sum_lambda_i = DIGAMMA(tot) ;
	    //E_log_p_z -= ploidy * L * K * *digamma_sum_lambda_i ;
	    // E_log_p_z -= ploidy * L * *digamma_sum_lambda_i ;
	    // d_KL += KL_Dirichlets(lambda_i, lambda_prior_i, K) ;
	}
	
	opt.fit[iter] = -d_KL + E_log_p_z + E_log_like + entropy ;
	
	if(opt.show_output && !(iter % opt.print_period) && thisproc == ROOT)
	    PRINT("%5d %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf %15.4lf\n", 
		  iter, d_KL , E_log_p_z , E_log_like , entropy, opt.fit[iter], 
		  iter > 0 ? (100 * (opt.fit[iter] - opt.fit[iter - 1]) / -opt.fit[iter - 1]) : log(-1)) ;

	if(opt.obey_stop_tol && iter >= opt.plateau_size && reached_plateau(iter, opt.fit)) break ;
    }
    while(++iter < opt.niters) ; 

    *_free_energy = opt.fit[iter - 1] ;
}
