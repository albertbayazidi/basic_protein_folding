/*

#include "sim.h"

void Run2DSim(gsl_rng * r,gsl_matrix *EnergyInteraction, Protein2D protein, int * sizeXY, char * folderPath, int folded, int T, int X ){



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


    free(ProteinInital);
    free(NNInital);
    free(mainSimPath);
    free(ProteinLast);
    free(NNLast);

}


*/
