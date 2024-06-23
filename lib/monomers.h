#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#define diffMonomers 20
//#define diffMonomers 10

void monomerValues(gsl_vector *monomers, gsl_rng *r, int nrMonomers);


void monomerEnergis(gsl_matrix *monomerE, gsl_rng *r);

