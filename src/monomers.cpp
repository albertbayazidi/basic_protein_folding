#include "monomers.h"
#include "utils.h"
#include "random.h"
#include "visual.h"


void monomerValues(gsl_vector *monomers, gsl_rng *r, int nrMonomers){
    int monomerValue;
    for (int i = 0; i < nrMonomers; i++){
       monomerValue =  radnomInt(r, diffMonomers) + 1;
       gsl_vector_set(monomers,i,monomerValue);
    }
}



void monomerEnergis(gsl_matrix *monomerE, gsl_rng *r){
    // Energy might have to be multiplyed with boltzmann constant
    long unsigned int i;
    long unsigned int j;
    for (i = 0; i < monomerE->size1; i++){
        for (j = 0; j <= i; j++){
            int interactionEnergi = -(radnomInt(r, 3) + 2);
            gsl_matrix_set(monomerE,i,j,interactionEnergi);
            gsl_matrix_set(monomerE,j,i,interactionEnergi);
            
        }
    }
}



