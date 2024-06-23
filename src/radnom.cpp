#include "random.h"

gsl_rng * radnomGenerator(){
	gsl_rng * r;
    const gsl_rng_type * T;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
    return r;
}

double randomUniform(gsl_rng * r){
    return gsl_rng_uniform(r);
}

int randomDir(gsl_rng * r){
    // UP = -2, down = 2, left = -1, right = 1
	int D = gsl_rng_uniform_int(r, 4);
    if (D >= 2){
        D = D - 1;
    }
    else{
        D = D -2;
    }
    return D;
}

int radnomInt(gsl_rng * r, int size){
    return gsl_rng_uniform_int(r, size);
}
