#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "visual.h"
#include "random.h"
#include "monomers.h"

void checkFilePointer(FILE *fptr);

void computeDirection(int *updatedKoords, int *currKoords, int Dir);

int findDirBasedOnKoords(int xSum,int ySum);

void checkOccupied(gsl_vector *koords, int* currkords, int* updatedkords, int index);

int interpRandomDir(int dir, int randomDir);

bool ileagalBoundaryPlacment(int *domainSize, int *nextKoords);

bool checkOppositeDir(int privDir, int nextDir);

void accecptUpdate(gsl_vector *koords, int *updatedKoords, int currMonomer);

bool covalentBond(int currMonomer,int nextMonomer, int *V1 ,int *V2); // not needed (I THINK)

int occupiedBound(gsl_matrix *neighbours, int v1Monomer);

void checkUp(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer);

void checkDown(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer);

void checkLeft(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer);

void checkRight(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer);

bool isPlaneConstructed2D(int privDir, int nextDir);

void possibleKoord2D(gsl_vector *koords, int* currkords,int* updatedkords, int privDir, int nextDir, int index);

bool isUnique(gsl_vector * last_dirs, int nextDir);

int findAbsSum(gsl_vector *V);
