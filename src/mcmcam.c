#include "psi.h"
#include "mcmc.h"

void increment_qmean(double *incr, double weight, double *qmean) {
    int ik, nK = opt.n * opt.K ;
    
    for(ik = 0 ; ik < nK ; ik++)
	qmean[ik] += incr[ik] * weight ;
}

void update_skip_indicators() {
    int l ;
    double p0, p1, U ;

    opt.nskip = 0 ; 
    for(l = 0 ; l < opt.L ; l++) {
	
	p1 = opt.skip_prior * exp(opt.snploglike_nostructure[l] - opt.snploglike[l]) ;
	p0 = 1 - opt.skip_prior ;
	
	U = runif(0, 1) ;

	opt.skip_locus[l] = U > p0 / (p0 +p1) ;

	if(FALSE && opt.p != NULL)
	    PRINT("l=%d: (%.4lf, %.4lf) %lf vs. %lf\tp0 = %.4lf\tU = %lf\t%s\n", 
		  l, opt.p[l*opt.K + 0], opt.p[l*opt.K + 1], 
		  opt.snploglike[l], opt.snploglike_nostructure[l], p0 / (p0 + p1), U, opt.skip_locus[l] ? "YES":"NO") ;

	if(opt.skip_locus[l]) opt.nskip++ ;
    }
}

bool update_skip_indicator(double loglike0, double loglike1, double prior1) {
    double p0, p1 ;

    p1 = prior1 * exp(loglike1 - loglike0) ;
    p0 = 1 - prior1 ;
    
    return runif(0,1)  >  p0/(p0 + p1) ;
}


void mcmcam(const unsigned char *x0, int n, int L, int K, double *qmean) {
    
    const int twon = 2*n ;
#if PARALLEL
    extern int L_tot, thisproc ;
#endif
    int l, j, i, ii, a, k, al, x_lj, z_lj, iter, geno ;
    double *nla, *mla, *nila, *nia_l ;
    double
	*prxz_lj, // joint probability of observed x_lj and z_lj over {1,...,K} at allele (l,j)
	*p, *p_l,// p_l[k] = cluster k allele frequencies at locus l
	*q,      // q[i,k] = proportion of individual i ancestry from cluster k
	*priorz_i,
	*w ;
    const unsigned char *x ;
    double p_prior_counts[NALS], *q_prior_counts, *w_prior_counts, prx_lj, loglike, like ;
    const int *zeros ;
    // FILE *q1 = fopen("/tmp/psi/out/q1", "w"), *qn = fopen("/tmp/psi/out/qn", "w") ;
    
    prxz_lj = ALLOC(K, double) ;
    nla = ALLOC(n * K, double) ;      // nla[i,k]    = #{(l,a)  : z_{ila} = k}
    nia_l = ALLOC(K * NALS, double) ; // nia_l[k,al] = #{(i,a)  : z_{ila} = k, x_{ila} = al}
    nila = ALLOC(K, double) ;         // nila[k]     = #{(i,l,a): z_{ila} = k}

    opt.w = w = ALLOC(K, double) ;
    w_prior_counts = ALLOC(K, double) ;
    mla = ALLOC(n * K, double) ;     /* mla[i,k] = sum_{i'la} I(z_{i'la} = k) skernel(i,i')
					proportion of alleles near i that are assigned to cluster k */
    
    p = ALLOC(K * L, double) ;
    q = ALLOC(n * K, double) ;
    zeros = CALLOC(n * K * NALS, double) ;
    q_prior_counts = ALLOC(K, double) ;

    opt.p = p ;

    /* sample initial values of p and q from prior */
    for(k = 0 ; k < K ; k++) {
	q_prior_counts[k] = 1 ;
	w_prior_counts[k] = 1 ;
    }
    for(al = 0 ; al < NALS ; al++) p_prior_counts[al] = 1 ;
    
    for(l = 0 ; l < L ; l++)
	for(k = 0 ; k < K ; k++)
	    p[l*K + k] = rbeta(p_prior_counts[1], p_prior_counts[0]) ; // NB p is frequency of second of two alleles!

    for(i = 0 ; i < n ; i++)
	rdirichlet(q_prior_counts, K, q + i*K) ;
    rdirichlet(w_prior_counts, K, w) ;
    
    /* sample skip indicators from prior */
    memcpy(opt.snploglike, opt.snploglike_nostructure, opt.L * sizeof(double)) ;
    update_skip_indicators() ;
    for(iter = -opt.burnin ; iter < opt.niters ; iter++) {
	
	opt.nskip = 0 ;
	x = x0 ;
	memcpy(nla, zeros, n * K * sizeof(double)) ;
	memcpy(mla, zeros, n * K * sizeof(double)) ;
	memcpy(nila, zeros, K * sizeof(double)) ;

	loglike = 0 ;
	
	for(l = 0 ; l < L ; l++) {
	    
	    if(FALSE && opt.skip_locus[l] && (iter >= 0) && (iter % opt.skip_update_interval)) {
		loglike += opt.snploglike[l] ;
		x += n ;
		opt.nskip++ ;
		opt.skipfreq[l] += 1.0 / opt.niters ;
		continue ;
	    }

	    memcpy(nia_l, zeros, K * NALS * sizeof(double)) ;
	    p_l = p + l*K ;
	    opt.snploglike[l] = 0 ;
	    like = 1 ;
	    
	    for(j = 0 ; j < twon ; j++) {
		
		i = j / 2 ;                            // j indexes chromosomes; i indexes individuals
		a = j % 2 ;                            // a == 0 if first allele, a == 1 if second allele

		// x_lj = get_allele(&x, a) ;
		geno = *x - '0' ;
		if(geno != MISSING) x_lj = (geno == 2) || (geno == 1 && a == 1) ;
		if(a == 1) x++ ;

		// MISSING DATA STUFF INSERTED WHILE TIRED Nov 2008

		priorz_i = opt.popmean ? w : q + i*K ;
		
		/* Update ancestry indicator z_lj and increment latent counts */
		prx_lj = 0 ;
		for(k = 0 ; k < K ; k++) {
		    prxz_lj[k] = priorz_i[k] ;
		    if(geno != MISSING) prxz_lj[k] *= (x_lj ? p_l[k] : 1 - p_l[k]) ;
		    prx_lj += prxz_lj[k] ;
		}
		like *= prx_lj ;
		if(like < SMALL) {
		    opt.snploglike[l] += log(like) ;
		    like = 1 ;
		}
		
		if(opt.updatez[i]) {
		    z_lj = sample(prxz_lj, K, prx_lj) ;
		    nila[z_lj]++ ;
		}
		else z_lj = opt.fixpop[i] ;
		
		if(opt.spatial)
		    for(ii = 0 ; ii < n ; ii++)
			mla[z_lj + ii*K] += opt.skernel[ii + i*n] ;
		
		if(geno != MISSING) nia_l[z_lj*NALS + x_lj]++ ;
		nla[i*K + z_lj]++ ;
	    }
	    opt.snploglike[l] += log(like) ;
	    loglike += opt.snploglike[l] ;

	    assert(!opt.skip_locus[l] || iter < 0 || !(iter % opt.skip_update_interval)) ;

	    if(FALSE && !(iter % opt.skip_update_interval)) {
		opt.skip_locus[l] = update_skip_indicator(opt.snploglike[l], opt.snploglike_nostructure[l], opt.skip_prior) ;
		if(opt.skip_locus[l]) {
		    opt.nskip++ ;
		    opt.skipfreq[l] += 1.0 / opt.niters ;
		}
	    }

	    /* Update allele frequencies at this locus */
	    for(k = 0 ; k < K ; k++)
		p_l[k] = rbeta(p_prior_counts[1] + nia_l[k*NALS + 1],
			       p_prior_counts[0] + nia_l[k*NALS + 0]) ;
	}
	
	/* Update ancestry proportions q */
	if(opt.spatial)
	    for(i = 0 ; i < n ; i++)
		rdirichlet_update(mla + i*K, q_prior_counts, K, q + i*K) ;
	else
	    for(i = 0 ; i < n ; i++)
		rdirichlet_update(nla + i*K, q_prior_counts, K, q + i*K) ;

	/* Update mean ancestry proportions */
	rdirichlet_update(nila, w_prior_counts, K, w) ;

	if(iter >= 0 && !(iter % opt.thin)) {
	    if(!(iter % opt.print_period)) PRINT("%d ", iter / opt.thin) ;
	    write_current_state_mcmc(q, iter) ; 
	}
	if(iter >= 0) increment_qmean(q, 1.0 / opt.niters, qmean) ;
	// if(iter >= 0) increment_qmean(nla, 1.0 / (2 * (L - opt.nskip) * opt.niters), qmean) ;
	
	// if(!(iter % opt.skip_update_interval)) PRINT("skipfreq[%d] = %.3lf\n", iter, opt.nskip / (double) opt.L) ;
    }
    PRINT("\n") ;
    // fclose(q1) ;
    // fclose(qn) ;
    
    free(zeros) ;
    free(q) ;
    free(prxz_lj) ;
    free(p) ;
    free(nla) ;
    free(nila) ;
    free(mla) ;
    free(w) ;
    free(nia_l) ;
}

void rdirichlet_update(double *observed, double *prior, int K, double *p) {
    int k ;
    double sum = 0 ;

    for(k = 0 ; k < K ; k++) {
	p[k] = rgamma(observed[k] + prior[k], 1) ;
	sum += p[k] ;
    }
    for(k = 0 ; k < K ; k++)
	p[k] /= sum ;
}

void write_current_state_mcmc(double *q, int iter) {
    FILE *f ;
    bool final = (iter < 0) ;
    int i ;
    struct timeval tv ;
    double now ;

    if(thisproc == ROOT) {
	if(final) {
	    sprintf(opt.fnamebuff, "%s/q", opt.outdir) ;
	    //else sprintf(opt.fnamebuff, "%s/q-%05d", opt.outdir, iter) ;
	    f = fopen(opt.fnamebuff, "w") ;
	    write_matrix_double(q, f, opt.K, opt.n, "%-8.5lf") ;
	    fclose(f) ;
	}
	else {
	    gettimeofday(&tv, NULL) ;
	    now = (double) tv.tv_sec + 1e-6 * (double) tv.tv_usec;
	    fprintf(opt.timefile, "%lf\n", now) ;
	    
	    write_matrix_double(q, opt.statefile, opt.K * opt.n, 1, "%-10.5lf ") ;
	    write_matrix_double(opt.w, opt.wstatefile, opt.K, 1, "%-10.5lf ") ;
	    
	    if(TRUE) for(i = 0 ; i < opt.n ; i++) {
		sprintf(opt.fnamebuff, "%s/qsam-%05d", opt.outdir, i + 1) ;
		f = fopen(opt.fnamebuff, "a") ;
		write_matrix_double(q + i*opt.K, f, opt.K, 1, "%-8.5lf") ;
		fclose(f) ;
	    }
	}
    }
    if(final) {
	sprintf(opt.fnamebuff, "%s/skipfreq", opt.outdir) ;
	f = fopen(opt.fnamebuff, "w") ;
	write_matrix_double(opt.skipfreq, f, 1, opt.L, "%.5lf") ;
	fclose(f) ;
    }
}
