// classes
#include "proteins2D.h"

// functions
#include "visual.h"
#include "random.h"
#include "utils.h"
#include "monomers.h"
#include "sim.h"

int main (void){
	printf ("Chosen seed = %lu\n", gsl_rng_default_seed);

	int sizeXY[2] = {100,100}; 	
	int nrMonomers = 100; 		
	int folded = 1;				// 1 = folded, 0 = unfolded
	int T = 10;					// Temperature
	int X = 40;				// nr of MC steps

	//char MainSim[150] = "_misc/data/monte_carlo_sim/task_2.5/main_sim/Protein2D.txt";

	char folderPath[100] = "_misc/data/monte_carlo_sim/task_2.5/";

	gsl_vector *koords = gsl_vector_alloc(2*nrMonomers);
	gsl_matrix *EnergyInteraction = gsl_matrix_alloc(diffMonomers, diffMonomers); 

	gsl_rng * r = radnomGenerator();
	
	
	// init obj
	monomerEnergis(EnergyInteraction, r);		// Makes interaction matrix

	Protein2D protein = Protein2D(r,koords,EnergyInteraction,sizeXY, folded); 

	//Run2DSim( r, EnergyInteraction, firstBoy,  sizeXY, folderPath, folded,  T, X );


	// Save inital protein structure and NN
    char P_Inital[100] =  "inital_and_last/inital";
    char* ProteinInital = fileloc(folderPath, P_Inital);
    char NN_Inital[100] =  "inital_and_last/NN_inital";
	char* NNInital = fileloc(folderPath, NN_Inital);


	FILE *firstptr;
	firstptr = fopen(ProteinInital, "w");
	protein.saveProtein2D(firstptr);
	saveMat(protein.neighbours, NNInital);
	fclose(firstptr);


	// Run MC
    char mainS[100] =  "main_sim/Protein2D";
    char* mainSimPath = fileloc(folderPath, mainS);

	printf("Energy = %d at t = 0\n", protein.energy);
	protein.MCStep(r, EnergyInteraction, T, X, mainSimPath); // implement RoG also
	printf("Energy = %d at t = %d \n", protein.energy, X);


	
	// saving last protein structure and NN
    char P_Last[100] =  "inital_and_last/last";
    char* ProteinLast = fileloc(folderPath, P_Last);
    char NN_Last[100] =  "inital_and_last/NN_last";
	char* NNLast= fileloc(folderPath, NN_Last);

	FILE *fptr;
	fptr = fopen(ProteinLast, "w");
	protein.saveProtein2D(fptr);
    saveMat(protein.neighbours, NNLast);
	fclose(fptr);

	// Free
    free(ProteinInital);
    free(NNInital);
    free(mainSimPath);
    free(ProteinLast);
    free(NNLast);
	
	gsl_vector_free(koords);
	gsl_matrix_free(EnergyInteraction);


	return 0;
}
