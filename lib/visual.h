#pragma once
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <gsl/gsl_linalg.h>
#include "utils.h"



void printVec(gsl_vector *v);

void printMat(gsl_matrix *A);

void saveMat(gsl_matrix *A, char *filename);

char* fileloc(char *filename, char *condition);