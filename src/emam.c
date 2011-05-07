#include "psi.h"
#include "em.h"

double emam(unsigned char *x0, double *px0, int *_n, int *_L, int *_K, double *qold0, double *mu0, 
		    double *alpha, double *F) {
    
    const int n = *_n, L = *_L, K = *_K, twon = 2*n ;
#if PARALLEL
    extern int L_tot, thisproc ;
    const double one_over_2L_tot = 1.0 / (2.0 * L_tot) ;
#else 
    const double one_over_2L_tot = 1.0 / (2.0 * L) ;
#endif
    int l, j, i, a, k, ik,
	geno, iter ;
    double allele_like, like, logpost, logpost_prev, logpost_prev_snp,
	*mu_l, *munum_l, *mudenom_l, 
	*f, *qold_i, *qnew_i, *qnew0, *qtmp, 
	*pz_lj, // probability distribution over {1,...,K} at allele (l,j), i.e. the quantity computed during the E-step
	pz_ljk, // temp storage for pz_lj[k] 
	*px, // pointer to probabilistic data (the beginning of which is always pointed to by px0)
	priorx_lj, // prior probability of reference allele (1) at allele copy (l,j)
	px_lj, // Pr(x_lj = 1), i.e. prob that current allele is a 1, with pnull component contributing priorx_lj 
	pnull // 
	;
    unsigned char *x ;

    /* 
       Fit structure admixture model via EM. 
       As parameters are updated, the likelihood is computed for the old set of parameters.
       To do so, two copies of q are required, but only one copy of mu.
       At the end of each iteration, the likelihood is that of qold and the (now erased) 'previous' value of mu.
    */

    munum_l = ALLOC(K, double) ;
    mudenom_l = ALLOC(K, double) ;
    qnew0 = ALLOC(n * K, double) ;
    pz_lj = ALLOC(K, double) ;
    f = ALLOC(K, double) ;
    
    for(k = 0 ; k < K ; ++k) {
	f[k] = (1.0 - F[k]) / F[k] ;
	munum_l[k] = mudenom_l[k] = 0.0 ;
    }
    logpost = log(0.0) ;
    
    iter = 0 ;
    do {
	for(ik = 0 ; ik < n*K ; ++ik) qnew0[ik] = 0.0 ;
	qnew_i = qnew0 ;
	
	logpost_prev = logpost ; logpost = 0.0 ; logpost_prev_snp = 0.0 ; like = 1.0 ;
	
	x = x0 ; px = px0 ;
	for(l = 0 ; l < L ; ++l) {
	    
	    if(FALSE) if(opt.skip_locus[l] || 
			 opt.snploglike[l] - opt.snploglike_nostructure[l] < opt.skip_tol * opt.snploglike_nostructure[l]) {
		   
		    opt.skip_locus[l] = TRUE ;
		    logpost += opt.snploglike_nostructure[l] ;
		    logpost_prev_snp = logpost ;
		    opt.snploglike[l] = opt.snploglike_nostructure[l] ;
		    if(opt.uncertain_data) px += 3 * n ;
		    else x += n ;
		    continue ;
		}
	    
	    mu_l = mu0 + l * K ;
		
	    for(j = 0 ; j < twon ; ++j) {
		
		/**************************************************************** 
		 WORK OUT WHAT ALLELE IS PRESENT ON THIS CHROMOSOME AT THIS LOCUS 
		 ****************************************************************/

		i = j / 2 ;                            // j indexes chromosomes; i indexes individuals
		a = j % 2 ;                            // a == 0 if first allele, a == 1 if second allele
		qold_i = qold0 + i * K ;
		
		priorx_lj = 0.0 ;
		for(k = 0 ; k < K ; k++)
		    priorx_lj += qold_i[k] * mu_l[k] ;

		if(opt.uncertain_data) {
		    pnull = 1 - px[0] - px[1] - px[2] ;
		    px_lj = px[2]  +  px[1] * (a == 1)  +  pnull * priorx_lj ;
		    if(a == 1) px += 3 ;
		}
		else {
		    geno = *x - '0' ;
		    if(a == 1) ++x ;
		
		    if(geno != MISSING)
			px_lj = (geno == 2) || (geno == 1 && a == 1) ; // have allele 1 with prob 1 if homozygote 
		                                                       // or if second allele at heterozygote
		    else
			px_lj = priorx_lj ;
		}
		
		/***********************************************************************************************
		 E-STEP: FORM THE POSTERIOR OVER THE MISSING DATA FOR THIS ALLELE USING CURRENT PARAMETER VALUES
		************************************************************************************************/
		
		allele_like = 0.0 ;
		for(k = 0 ; k < K ; k++)
		    allele_like += pz_lj[k] = ( (1-px_lj) * (1-mu_l[k])  +  px_lj * mu_l[k] ) * qold_i[k] ;
		like *= allele_like ;
		
		if(like < SMALL) {
		    logpost += log(like) ;
		    like = 1.0 ;
		}
		
		/***************************************************************************************************** 
		 M-STEP: USE CURRENT POSTERIOR ON MISSING DATA AT THIS ALLELE TO INCREMENT NEW ESTIMATES OF PARAMETERS
		******************************************************************************************************/
		
		qnew_i = qnew0 + i * K ;
		for(k = 0 ; k < K ; ++k) {
		    pz_ljk = pz_lj[k] /= allele_like ;
		    qnew_i[k] += one_over_2L_tot * pz_ljk ;
		    munum_l[k] += pz_ljk * px_lj ;
		    mudenom_l[k] += pz_ljk ;
		}
	    }
	    
	    if(TRUE || opt.model == EMAM) {
		for(k = 0 ; k < K ; ++k) {
		    mu_l[k] = munum_l[k] / mudenom_l[k] ;
		    munum_l[k] = mudenom_l[k] = 0.0 ;
		}
	    }
	    else {
		ERROR("Not doing correlated model\n") ;
		for(k = 0 ; k < K ; ++k) {
		    /*
		      NB code in this section has been altered in keeping with the move to standard indexing rather
		      than pointer increments. However the alterations have not been tested.
		    */
		    // Following two lines for mu prior without pseudocounts of 1 of each allele in each cluster
		    // logpost += (alpha[l] * f[k] - 1) * log(mu[k]) + ((1 - alpha[l]) * f[k] - 1) * log(1 - mu[k]) ;
		    // mu[k] = (munum_l[k] + alpha[l] * f[k] - 1) / (mudenom_l[k] + f[k] - 2) ;
		    // Following two lines for mu prior with pseudocounts of 1 of each allele in each cluster
		    logpost += (alpha[l] * f[k]) * log(mu_l[k]) + ((1 - alpha[l]) * f[k]) * log(1 - mu_l[k]) ;
		    mu_l[k] = (munum_l[k] + alpha[l] * f[k]) / (mudenom_l[k] + f[k]) ;
		    munum_l[k] = mudenom_l[k] = 0.0 ;
		}
	    }
	    logpost += log(like) ;
	    like = 1.0 ;
	    opt.snploglike[l] = logpost - logpost_prev_snp ;
	    logpost_prev_snp = logpost ;
	}
#if PARALLEL
	MPI_Allreduce(&logpost, &logpost_tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ; // sum log posteriors and 
	logpost = logpost_tmp ;
	MPI_Allreduce(qnew0, qold0, n * K, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ; // ancestry proportions across processes
#else
	qtmp = qold0 ; qold0 = qnew0 ; qnew0 = qtmp ;
#endif
	
	if(!(iter % opt.print_period) && thisproc == ROOT)
	    PRINT("%s\tlog%s[%d] = %.12f\n", timestring(), opt.model == EMAM ? "like" : "post", iter, logpost) ;
	
	if(iter > 0 && iter + 1 < opt.niters && !(iter % opt.output_period))
	    write_current_state_emam(qold0, mu0, logpost, iter, K, n, L) ;
    }
    while(++iter < opt.niters) ;
    
#if !PARALLEL // returned value of q is associated with returned likelihood; not so for mu
    memcpy(qold0, qnew0, n * K * sizeof(double)) ; 
#endif
    
    free(munum_l) ;
    free(mudenom_l) ;
    free(opt.niters % 2 ? qold0 : qnew0) ; // they get swapped each iteration...
    free(pz_lj) ;
    free(f) ;

    return logpost ;
}

void write_current_state_emam(double *q, double *mu, double logpost, int iter, int K, int n, int L) {
    FILE *f ;
    bool final = (iter < 0) ;
    
    if(thisproc == ROOT) {
	if(final) sprintf(opt.fnamebuff, "%s/q", opt.outdir) ;
	else sprintf(opt.fnamebuff, "%s/q-%05d", opt.outdir, iter) ;
	f = fopen(opt.fnamebuff, "w") ;
	write_matrix_double(q, f, K, n, "%-8.5lf") ;
	fclose(f) ;
	
	if(final) sprintf(opt.fnamebuff, "%s/log%s", opt.outdir, opt.model == EMAM ? "like" : "post") ;
	else sprintf(opt.fnamebuff, "%s/log%s-%05d", opt.outdir, opt.model == EMAM ? "like" : "post", iter) ;
	f = fopen(opt.fnamebuff, "w") ;
	fprintf(f, "%lf\n", logpost) ;
	fclose(f) ;
    }
    
    /* Every process writes out allele frequencies at its own set of SNPs */
#if PARALLEL
    if(final) sprintf(opt.fnamebuff, "%s/mu-%02d", opt.outdir, thisproc) ;
    else sprintf(opt.fnamebuff, "%s/mu-%05d-%02d", opt.outdir, iter, thisproc) ;
#else
    if(final) sprintf(opt.fnamebuff, "%s/mu", opt.outdir) ;
    else sprintf(opt.fnamebuff, "%s/mu-%05d", opt.outdir, iter) ;
#endif
    f = fopen(opt.fnamebuff, "w") ;
    write_matrix_double(mu, f, K, L, "%-7.4lf") ;
    fclose(f) ;

#if PARALLEL
    /* Every process writes out log likes at its own set of SNPs */
    if(final) sprintf(opt.fnamebuff, "%s/snploglike-%02d", opt.outdir, thisproc) ;
    else sprintf(opt.fnamebuff, "%s/snploglike-%05d-%02d", opt.outdir, iter, thisproc) ;
#else
    if(final) sprintf(opt.fnamebuff, "%s/snploglike", opt.outdir) ;
    else sprintf(opt.fnamebuff, "%s/snploglike-%05d", opt.outdir, iter) ;
#endif
    f = fopen(opt.fnamebuff, "w") ;
    write_matrix_double(opt.snploglike, f, 1, L, "%-7.4lf") ;
    fclose(f) ;
}

double *read_mu(char *fname, int K, int L_sub, int L_tot, bool *snp_include) {
    // cat mu* > mu; mu_infile is mus at all SNPs
    FILE *f ;
    int k ;
    double *mu = ALLOC(L_sub * K, double) ;
    bool *trues = ALLOC(K, bool) ;
    if(mu == NULL) ERROR("proc %d: mu array malloc error:", thisproc) ;
    for(k = 0 ; k < K ; k++) trues[k] = TRUE ;

    f = fopen(fname, "r") ;
    read_submatrix_double(f, K, trues, L_tot, snp_include, "%lf", mu) ;
    fclose(f) ;
    
    return mu ;
}


double *initialise_F(int K) {
    int k ;
    double *F ;

    F = ALLOC(K, double) ;
    for(k = 0 ; k < K ; k++)
	F[k] = 1.0 ;

    return F ;
}

