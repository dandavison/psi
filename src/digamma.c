#include <dan.h>
#include "digamma.h"

double DIGAMMA_SMALL_TABLE_INCR = (DIGAMMA_TABLE_MID - DIGAMMA_TABLE_MIN) / (DIGAMMA_TABLE_SIZE - 1.0) ;
double DIGAMMA_BIG_TABLE_INCR = (DIGAMMA_TABLE_MAX - DIGAMMA_TABLE_MID) / (DIGAMMA_TABLE_SIZE - 1.0) ;
double *digamma_small_table, *digamma_big_table, DIGAMMA_SMALL_TABLE_FAC, DIGAMMA_BIG_TABLE_FAC ;


void *make_digamma_table() {
    int i ;
    double *small_tab, *big_tab, x ;
    
    digamma_small_table = ALLOC(DIGAMMA_TABLE_SIZE, double) ;
    digamma_big_table = ALLOC(DIGAMMA_TABLE_SIZE, double) ;

    for(i = 0, x = DIGAMMA_TABLE_MIN ; i < DIGAMMA_TABLE_SIZE ; ++i, x += DIGAMMA_SMALL_TABLE_INCR)
	digamma_small_table[i] = digamma(x) ;

    for(i = 0 ; i < DIGAMMA_TABLE_SIZE ; ++i, x += DIGAMMA_BIG_TABLE_INCR)
	digamma_big_table[i] = digamma(x) ;
    
    DIGAMMA_SMALL_TABLE_FAC = 1.0 / DIGAMMA_SMALL_TABLE_INCR ;
    DIGAMMA_BIG_TABLE_FAC = 1.0 / DIGAMMA_BIG_TABLE_INCR ; 
}

double digamma_lookup(double x) {
    int el ;

    if(x <= DIGAMMA_TABLE_MID) {
	//PRINT("%lf -> %d\n", x, (int) ROUND(x * DIGAMMA_SMALL_TABLE_FAC) - 1) ;
	// return digamma_small_table[(int) ROUND(x * DIGAMMA_SMALL_TABLE_FAC) - 1] ;
	// return digamma(x) ;
	el = (int) ROUND(x * DIGAMMA_SMALL_TABLE_FAC) ;
	return 0.5 * (digamma_small_table[el - 1] + digamma_small_table[el]) ;
    }
    else if(FALSE && x < DIGAMMA_TABLE_MAX) {
	el = (int) ROUND((x - DIGAMMA_TABLE_MID) * DIGAMMA_BIG_TABLE_FAC) ;
	return 0.5 * (digamma_big_table[el - 1] + digamma_big_table[el]) ;
    }
    else return digamma(x) ;
}


void digamma_test(double *x, int *n, double *ans) {
    int i ;
    make_digamma_table() ;
    
    // PRINT("%lf\t%lf\n", DIGAMMA_SMALL_TABLE_INCR, DIGAMMA_BIG_TABLE_INCR) ;
    // PRINT("%lf\t%lf\n", DIGAMMA_SMALL_TABLE_FAC, DIGAMMA_BIG_TABLE_FAC) ;
    // print_dbl(digamma_small_table, 10) ;
    // print_dbl(digamma_big_table, 10) ;
    
    for(i = 0 ; i < *n ; ++i)
	ans[i] = digamma_lookup(x[i]) ;
}

