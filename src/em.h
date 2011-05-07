double emam(unsigned char *x0, double *px, int *_n, int *_L, int *_K, 
		    double *qold0, double *mu0, double *alpha0, double *F0) ;
void write_current_state_emam(double *q, double *mu, double logpost, int iter, int K, int n, int L) ;
double *read_mu(char *fname, int K, int L_sub, int L_tot, bool *snp_include) ;
double *initialise_F(int K) ;
