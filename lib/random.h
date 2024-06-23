#pragma once
#include <gsl/gsl_rng.h>

gsl_rng * radnomGenerator();

double randomUniform(gsl_rng * r);

int randomDir(gsl_rng * r);

int radnomInt(gsl_rng * r, int size);