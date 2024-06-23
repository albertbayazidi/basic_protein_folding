#pragma once
#include <gsl/gsl_matrix.h>
#include <iostream>

#include "utils.h"
#include "random.h"
#include "monomers.h"
#include "visual.h"

//MADE WITH 2D IN MIND
class Protein2D {
    public:
        int energy = 0;             
        gsl_matrix * neighbours;    //MADE WITH 2D IN MIND
        gsl_vector * koords;        //MADE WITH 2D IN MIND
        gsl_vector * monomers;      //MADE WITH 2D IN MIND
        gsl_vector * dirArray;      //MADE WITH 2D IN MIND
        Protein2D(gsl_rng * r,gsl_vector *koords, gsl_matrix *EnergyInteraction, int * sizeXY, int folded);                     //constructor 
        
        // ~Protein2D();               //destructor should free the allocated memory

        // methods
        void initMonomer(gsl_rng *r, int *sizeXY);
        void initUnfoldMonomer(gsl_rng *r, int *sizeXY);
        bool checkSelfIntersection(int *currKoords, long unsigned int currMonomer);
        void initDirArrayFolded();
        void nearestNeighbour();
        void copmuteEnergy(gsl_matrix *EnergyInteraction);

        // MONTE CARLO METHODS
        void MonteCarloSweep(gsl_rng *r, gsl_matrix *EnergyInteraction,FILE *fptr, int T);
        void MCStep(gsl_rng *r, gsl_matrix *EnergyInteraction, int T, int X, char * folderPath);

        // HELPER FUNCTIONS
        void __UpdatedMonomerKoordsNeighbour__(int *updatedkords, int randomAminoAcidIndex, int nextDir, int privDir);
        int __HandleEndsMonomers__(gsl_matrix *EnergyInteraction, gsl_rng *r, int *updatedkords, int T,int dir,int dirDir,int index);
        void __UpdatedEndsMonomerKoordsNeighbour__(int *updatedkords, int index, int dir);

        // Logger related
        void logger(FILE *fptr);    
        void saveProtein2D(FILE *fptr);
        void saveNeighbours2D(FILE *fptr);
        void endToEndDistance(FILE *fptr);
        void saveEnergy(FILE *fptr);

        // NOT USED
        //int checkEnergy(gsl_matrix *EnergyInteraction, int *updatedkoords, int privDir, int nextDir,int index);    // BRUKER IKKE DENNE MED MINDRE TING GÅR SAKTE
        //int __HandleLastMonomers__(gsl_matrix *EnergyInteraction, gsl_rng *r, int *updatedkords, int T);
};

void Protein2D::saveProtein2D(FILE *fptr){
    fprintf(fptr,"%ld\n", this->koords->size);
    for (unsigned int i = 0; i < this->koords->size; i+=2)
    {
        fprintf(fptr,"%f,",gsl_vector_get(this->koords,i));
        fprintf(fptr,"%f,",gsl_vector_get(this->koords,i+1));
        fprintf(fptr,"%f\n",gsl_vector_get(this->monomers,int(i/2)));
    }
}

void Protein2D::saveNeighbours2D(FILE *fptr){
    fprintf(fptr,"%ld\n", this->neighbours->size1);
    fprintf(fptr,"%ld\n", this->neighbours->size2);
    // KAN HENDE DETTE BØR BYGGES BORT OG IKKE RETT NED
    for (unsigned int i = 0; i < this->neighbours->size1; i++)
    {   
        for (unsigned int j = 0; j < this->neighbours->size2; j++)
        {
            fprintf(fptr,"%f,", gsl_matrix_get(this->neighbours,i,j));
        }
        fprintf(fptr, "\n");

    }
}

void Protein2D::endToEndDistance(FILE *fptr){
    int x1 = gsl_vector_get(this->koords,0);
    int y1 = gsl_vector_get(this->koords,1);
    int x2 = gsl_vector_get(this->koords,this->koords->size-2);
    int y2 = gsl_vector_get(this->koords,this->koords->size-1);
    double distance = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
    fprintf(fptr,"%f\n", distance);
}

void Protein2D::saveEnergy(FILE *fptr){
    fprintf(fptr,"%d\n", this->energy);
}

void Protein2D::logger(FILE *fptr){
    saveProtein2D(fptr);
    saveNeighbours2D(fptr);
    endToEndDistance(fptr);
    // SHOULD ADD ROG
    saveEnergy(fptr);
}

Protein2D::Protein2D(gsl_rng * r, gsl_vector *koords, gsl_matrix *EnergyInteraction, int * sizeXY, int folded){
    int nr_monomers = (int)((koords->size)/2);
    this->koords = koords;
    this->monomers = gsl_vector_calloc(nr_monomers);
    this->dirArray = gsl_vector_calloc(nr_monomers-1);
    this->neighbours = gsl_matrix_calloc((int)(this->koords->size/2), 2);
    gsl_matrix_set_all(this->neighbours,-1); 

    monomerValues(this->monomers, r, nr_monomers);         

    if (folded == 0){
        initUnfoldMonomer(r, sizeXY);
        gsl_vector_set_all(this->dirArray,1);
        this->energy = 0;
    }
    else{
        initMonomer(r, sizeXY);
        //initDirArrayFolded();
        nearestNeighbour();
        copmuteEnergy(EnergyInteraction);   
    }
}

bool Protein2D::checkSelfIntersection(int *currKoords, long unsigned int currMonomer){ //should also take care of opposite directions
    int index = 0;
    for (long unsigned int i = 0; i < (currMonomer); i++){  
        int koordX = gsl_vector_get(this->koords,index);
        int koordY = gsl_vector_get(this->koords,index+1); 
        if ((currKoords[0] == koordX) && (currKoords[1] == koordY)){
            return true;
        }
        index += 2;
    }
    return false;
}

void Protein2D::initUnfoldMonomer(gsl_rng *r,int *sizeXY){
    int currMonomerIndex = 0;
    for (long unsigned int  i =0; i < this->monomers->size; i++){
        gsl_vector_set(this->koords, currMonomerIndex, (int)(sizeXY[0]/2)+i);
        gsl_vector_set(this->koords, currMonomerIndex+1, (int)(sizeXY[0]/2));
        currMonomerIndex +=2;
    }
}

void Protein2D::initMonomer(gsl_rng *r, int *sizeXY){
    gsl_vector *laste_dirs = gsl_vector_calloc(4);
    int choosenDir = 0; 
    int xSum = 0 , ySum = 0;
    int koordX = radnomInt(r, sizeXY[0]);
    int koordY = radnomInt(r, sizeXY[1]);
    int currKoords[2]= {0,0};
    currKoords[0] = koordX;
    currKoords[1] = koordY;

    // set the first monomer
    accecptUpdate(this->koords, currKoords, 0);

    // pick a random direction
    int privDir = randomDir(r);

    // keep track of the current monomer
    long unsigned int currMonomer = 1;
    int currMonomerIndex = 2; 

    // seting up the next monomer
    int updatedKoords[2] = {0,0};
    currKoords[0] = koordX;
    currKoords[1] = koordY;

    // check out the move
    computeDirection(updatedKoords, currKoords, privDir); 
    while (ileagalBoundaryPlacment(sizeXY, updatedKoords)){
        privDir = randomDir(r);
        computeDirection(updatedKoords, currKoords, privDir); 
    }

    // accept the move
    currKoords[0] = updatedKoords[0];
    currKoords[1] = updatedKoords[1];
    accecptUpdate(this->koords, currKoords, currMonomerIndex);

    // finds the first direction based on displacement
    xSum = currKoords[0] - gsl_vector_get(this->koords,currMonomerIndex-2);
    ySum = currKoords[1]-gsl_vector_get(this->koords,currMonomerIndex-1);
    choosenDir = findDirBasedOnKoords(xSum,ySum);
    
    // accept the first direction
    gsl_vector_set(this->dirArray, 0, choosenDir);

    // keep track of the current monomer
    while (currMonomer < (this->monomers->size-1)){
    
        // pick a random direction and moves node to the new koords
        int nextDir = randomDir(r); 
        computeDirection(updatedKoords, currKoords, nextDir);

        gsl_vector_set_all(laste_dirs,0);
        
        // check if the move is valid
        while (Protein2D::checkSelfIntersection(updatedKoords,currMonomer) 
        || ileagalBoundaryPlacment(sizeXY, updatedKoords)){
            if (Protein2D::checkSelfIntersection(updatedKoords,currMonomer)){

                gsl_vector_set(laste_dirs,0,nextDir);
                long unsigned int indexDir = 1;      

                int last_dir_abs_sum = 0;
                last_dir_abs_sum = findAbsSum(laste_dirs);
                
                while (last_dir_abs_sum != 6 && Protein2D::checkSelfIntersection(updatedKoords,currMonomer)){
                    
                    if (isUnique(laste_dirs, nextDir)){
                        gsl_vector_set(laste_dirs,indexDir,nextDir);
                        indexDir++;
                    }

                    last_dir_abs_sum = 0;
                    last_dir_abs_sum = findAbsSum(laste_dirs);

                    
                    nextDir = randomDir(r);
                    computeDirection(updatedKoords, currKoords, nextDir);

                }
                if (last_dir_abs_sum == 6){
                    // undo move and pick a dir that is not the same as the one we just came from and take the move and give a new dir
                    int stuckVal = (int)gsl_vector_get(this->dirArray,currMonomer-1);

                    // revert the move
                    updatedKoords[0] = 0;
                    updatedKoords[1] = 0;
                    
                    currKoords[0] = (int)gsl_vector_get(this->koords,currMonomerIndex-2);
                    currKoords[1] = (int)gsl_vector_get(this->koords,currMonomerIndex-1);
                    
                    gsl_vector_set(this->koords,currMonomerIndex,0);
                    gsl_vector_set(this->koords,currMonomerIndex,0);

                    nextDir = randomDir(r);
                    computeDirection(updatedKoords, currKoords, nextDir);

                    while (nextDir == stuckVal 
                    || Protein2D::checkSelfIntersection(updatedKoords,currMonomer) 
                    || ileagalBoundaryPlacment(sizeXY, updatedKoords)){

                        nextDir = randomDir(r);
                        computeDirection(updatedKoords, currKoords, nextDir);
                    }
                    // accept the move
                    currKoords[0] = updatedKoords[0];
                    currKoords[1] = updatedKoords[1];

                    gsl_vector_set(this->dirArray, currMonomer-1, nextDir);
                    accecptUpdate(this->koords, currKoords, currMonomerIndex); 

                    // give next dir
                    nextDir = randomDir(r);
                    computeDirection(updatedKoords, currKoords, nextDir);
                   
                    // reset the laste_dirs
                    gsl_vector_set_all(laste_dirs,0);
                }
                
            }
            else{ // ileagalBoundaryPlacment
                nextDir = randomDir(r);
                computeDirection(updatedKoords, currKoords, nextDir);
            }
        }

        // update the monomer to node location
        currMonomer++;
        currMonomerIndex += 2;  
        currKoords[0] = updatedKoords[0];
        currKoords[1] = updatedKoords[1];
        

        // accept the move and store the direction
        gsl_vector_set(this->dirArray, currMonomer-1, nextDir);
        accecptUpdate(this->koords, currKoords, currMonomerIndex);  
    }
    printf("Finished constructing the Protein2D\n");
    gsl_vector_free(laste_dirs);
}

void Protein2D::initDirArrayFolded(){
    for (long unsigned int i = 0; i < this->koords->size-2; i+=2){
        int xSum = 0;
        int ySum = 0;
        int xCurr = gsl_vector_get(this->koords,i);
        int yCurr = gsl_vector_get(this->koords,i+1);
        int xNext = gsl_vector_get(this->koords,i+2);
        int yNext = gsl_vector_get(this->koords,i+3);
        xSum = xNext - xCurr;
        ySum = yNext - yCurr;
        int Dir = findDirBasedOnKoords(xSum, ySum);
        gsl_vector_set(this->dirArray, (int)(i/2), Dir);
    }

}

void Protein2D::nearestNeighbour(){
    // function stores the nearest neighbour for each monomer.
    // the output is a matrix where the nr of each monomer is stored not the specific koords, nor the specific amino acid.
    // MABYE STORE THE NEAREST NEIGHBOUR COORDS 

    this->neighbours = gsl_matrix_calloc((int)(this->koords->size/2), 2);
    gsl_matrix_set_all(this->neighbours,-1); 
    int size1 = this->koords->size;

    //find out if a koords has a neighbour
    for (int i = 0; i < size1; i+=2){        // vector one
    int v1Monomer = (int)(i/2);
    int v1X = gsl_vector_get(this->koords, i);
    int v1Y = gsl_vector_get(this->koords, i+1);
    int V1[2] = {0,0};
    V1[0]  = v1X;
    V1[1]  = v1Y;

        for (int j = i + 4; j < size1; j+=2){    // vector two
            int v2Monomer = (int)(j/2);
            int v2X = gsl_vector_get(this->koords, j);
            int v2Y = gsl_vector_get(this->koords, j+1);
            int V2[2] = {0,0};
            V2[0]  = v2X;
            V2[1]  = v2Y;
            // MABYE STORE THE NEAREST NEIGHBOUR COORDS 
            checkUp(this->neighbours, V1, V2, v1Monomer, v2Monomer);
            checkDown(this->neighbours, V1, V2, v1Monomer, v2Monomer);
            checkLeft(this->neighbours, V1, V2, v1Monomer, v2Monomer);
            checkRight(this->neighbours, V1, V2, v1Monomer, v2Monomer);

        }
    }
}

void Protein2D::copmuteEnergy(gsl_matrix *EnergyInteraction){
    // made for 2D 
    // setting the energy to zero so that the next time the function is called the energy is not added to the previous energystate 
    int EnergyJ = 0;
    int EnergyK = 0;
    this->energy = 0;
    for (long unsigned int i = 0; i < this->monomers->size; i++){
        int neighbourMonomerFirst = gsl_matrix_get(this->neighbours,i,0);
        int neighbourMonomerSecond = gsl_matrix_get(this->neighbours,i,1);
        int aminoAcid = gsl_vector_get(this->monomers,i);

        if (neighbourMonomerFirst != -1){           
            int neighbourAminoAcidFirst = gsl_vector_get(this->monomers,neighbourMonomerFirst);        
            EnergyJ = gsl_matrix_get(EnergyInteraction, aminoAcid-1, neighbourAminoAcidFirst-1);
        }
        else{
            EnergyJ = 0;
        }
        if (neighbourMonomerSecond != -1){
            int neighbourAminoAcidSecond = gsl_vector_get(this->monomers,neighbourMonomerSecond);
            EnergyK = gsl_matrix_get(EnergyInteraction, aminoAcid-1, neighbourAminoAcidSecond-1);
        }
        else{
            EnergyK = 0;
        }
        
        this->energy += EnergyJ + EnergyK;
    }
}

/* UNUSED METHOD FOR CHECKING ENERGY UPDATED
int Protein2D::checkEnergy(gsl_matrix *EnergyInteraction, int *updatedkoords, int privDir, int nextDir,int index){
    // KANSKJE BARE SEND RETUR NY VERDI OG IKKE OPPDATER DEN HER, GJØR DET HELLER NÅR VI SER AT DEN ER AKKSEPTERT
    // FOR AT DENNE METODEN SKAL FUNGERE TRENGER VI ET STEG SOM FJERNER ENERGIEN OM VI BRYTER ET BOND ( FORSETTER DERFOR MEG MER KOMPUTER INTESIV METODE)
    int updateDir = - nextDir; 
    int TempDir[2] = {privDir,updateDir};
    //int bond_index[2] = {0,0};
    //int nr_new_bonds = 0;
    for (int i = 0; i < 1; i++){

        int checkKords[2] = {0,0};
        checkKords[0] = updatedkoords[0];       // not actually needed (I THINK)
        checkKords[1] = updatedkoords[1];       // not actually needed (I THINK)
        int dir = TempDir[i];
        computeDirection(checkKords, updatedkoords, dir);

        int currMonomer = 0;
        for (int j = 0; j < (int)this->koords->size/2; j+=2){
            if (currMonomer == index){
                continue;
            }
            int currMonomerX = gsl_vector_get(this->koords,j);
            int currMonomerY = gsl_vector_get(this->koords,j+1);

            if (currMonomerX == checkKords[0] && currMonomerY == checkKords[1]){
                //bond_index[i] = currMonomer;
                //nr_new_bonds++;
                break;
            }
            else{
                //do nothing
            }
            currMonomer ++;
        }
        // found a new bond, now we need to check the energy 
        int aminoAcid = gsl_vector_get(this->monomers,index);
        int neighbourAminoAcid = gsl_vector_get(this->monomers,currMonomer);

        int EnergyJ = gsl_matrix_get(EnergyInteraction, aminoAcid-1, neighbourAminoAcid-1);
        this->energy += EnergyJ;

    }    
}
*/

void Protein2D::MonteCarloSweep(gsl_rng *r, gsl_matrix *EnergyInteraction, FILE *fptr, int T){
    int failedIndex = -1;
    int updatedkords[2] = {0,0};
    for (int n = 0; n < (int)this->koords->size/2; n++){ // change back to
        int randomAminoAcidIndex = radnomInt(r, (int)this->koords->size/2);
        while (failedIndex == randomAminoAcidIndex){
            randomAminoAcidIndex = radnomInt(r, (int)this->koords->size/2); 
        }
        if (randomAminoAcidIndex == 0 ){  // handle the first monomer
            int index = 0;
            int nextDir = gsl_vector_get(dirArray,0);
            int nextNextDir = gsl_vector_get(dirArray,1); 
            failedIndex = __HandleEndsMonomers__(EnergyInteraction, r, updatedkords, T, nextDir, nextNextDir, index);
        }
        else if (randomAminoAcidIndex == (int)koords->size/2-1){ // handle the last monomer
            int index = (int)koords->size/2-1;
            int privDir = gsl_vector_get(dirArray,index-1); 
            int privPrivDir = gsl_vector_get(dirArray,index-2);
            failedIndex = __HandleEndsMonomers__(EnergyInteraction, r, updatedkords, T, privDir, privPrivDir, index);

        }
        else{// handle internal monomers
            int privDir = gsl_vector_get(dirArray,randomAminoAcidIndex-1);
            int nextDir = gsl_vector_get(dirArray,randomAminoAcidIndex); 
            if (isPlaneConstructed2D(privDir,nextDir)){ 
                int privEnergy = this->energy;
                int currkords[2] = {0,0};
                currkords[0] = gsl_vector_get(this->koords,2*randomAminoAcidIndex);
                currkords[1] = gsl_vector_get(this->koords,2*randomAminoAcidIndex+1);

                // finn mulig koords for forflytting og sjekk validitet
                possibleKoord2D(this->koords, currkords, updatedkords, privDir,  nextDir, randomAminoAcidIndex);

                if (updatedkords[0] == currkords[0] && updatedkords[1] == currkords[1]){
                    //printf("no update\n");
                    failedIndex = randomAminoAcidIndex;
                }
                else{
                    // update koords og dirArray
                    __UpdatedMonomerKoordsNeighbour__(updatedkords, randomAminoAcidIndex, nextDir, privDir);
                    copmuteEnergy(EnergyInteraction);

                    int newEnergy = this->energy;
                    double u = randomUniform(r);
                    double a = std::min(1.0,exp((double)(privEnergy-newEnergy)/T));
                    
                    if (a > u){
                        // accept the move
                        failedIndex = -1;
                    }
                    else{// revert the changes
                        __UpdatedMonomerKoordsNeighbour__(currkords, randomAminoAcidIndex, privDir, nextDir);
                        this->energy = privEnergy;
                        failedIndex = -1;
                    }
                }
            }
            else{
                //printf("not on a plane\n");
                failedIndex = randomAminoAcidIndex;
            }
        }

        fprintf(fptr,"%d\n", n);
        logger(fptr);
        
    }
}

void Protein2D::MCStep(gsl_rng *r, gsl_matrix *EnergyInteraction, int T, int X, char * folderPath){
    FILE *fptr;
    fptr = fopen(folderPath, "w");
    checkFilePointer(fptr);
    fprintf(fptr,"%d\n", X);

    for (int i = 0; i < X; i++){
        printf("\rMCStep nr: %d / %d ", i+1 , X);
        fflush(stdout);
        MonteCarloSweep(r, EnergyInteraction, fptr, T);
    }
    printf("\n");

    fclose(fptr);
}

void Protein2D::__UpdatedMonomerKoordsNeighbour__(int *updatedkords, int randomAminoAcidIndex, int nextDir, int privDir){
    accecptUpdate(this->koords, updatedkords, 2*randomAminoAcidIndex); 
    gsl_vector_set(this->dirArray,randomAminoAcidIndex-1,nextDir);
    gsl_vector_set(this->dirArray,randomAminoAcidIndex,privDir);
    nearestNeighbour();
}

void Protein2D::__UpdatedEndsMonomerKoordsNeighbour__(int *updatedkords, int index, int dir){
    if (index==0){
        accecptUpdate(this->koords, updatedkords, 0); 
        gsl_vector_set(this->dirArray,0,dir);
        nearestNeighbour();
    }
    else{
        accecptUpdate(this->koords, updatedkords, 2*index); 
        gsl_vector_set(this->dirArray,index-1,dir); 
        nearestNeighbour();
    } 
}

int Protein2D::__HandleEndsMonomers__(gsl_matrix *EnergyInteraction, gsl_rng *r, int *updatedkords, int T, int dir, int dirDir, int index){
    // BE CAREFULL THAT THE CHANGES ARE REVERTED (IF nessacery) BEFORE RETURNING
    int currkords[2] = {0,0};
    int privEnergy = this->energy;
    if (abs(dir) + abs(dirDir) != 3){ 
        // two directions
        int randomDir = radnomInt(r,2);
        int newDir = interpRandomDir(dir, randomDir);

        // check if occupied
        currkords[0] = gsl_vector_get(this->koords,2*index);
        currkords[1] = gsl_vector_get(this->koords,2*index+1);

        if (index == 0){
            possibleKoord2D(this->koords, currkords, updatedkords, newDir, dir, index); 
        }
        else{
            possibleKoord2D(this->koords, currkords, updatedkords, dirDir, newDir, index);
        } 

        // MUST check if the new koords are valid and compute energy and such
        if (updatedkords[0] == currkords[0] && updatedkords[1] == currkords[1]){
            int failedIndex = index;
            return failedIndex;
        }
        else{
            __UpdatedEndsMonomerKoordsNeighbour__(updatedkords, index, newDir);
            copmuteEnergy(EnergyInteraction);

            int newEnergy = this->energy;
            double u = randomUniform(r);
            double a = std::min(1.0,exp((double)(privEnergy-newEnergy)/T));
            
            if (a > u){
                // accept the move
                //printf("accepted change at ends\n");
                return -2;
            }
            else{// revert the changes
                //printf("reverting changes at ends\n");
                __UpdatedEndsMonomerKoordsNeighbour__(currkords, index, dir);
                this->energy = privEnergy;
                return -1;          // since this move is based on probability we need to reset the failedIndex as the move has a chance of being accepted next time 
            }
        }
    }
    else{ //FEIL INNI HER ET STED!!!!!!!!!1
        // one direction
        // check if occupied and sett failedIndex = randomAminoAcidIndex;
        currkords[0] = gsl_vector_get(this->koords,2*index);
        currkords[1] = gsl_vector_get(this->koords,2*index+1);

        //printf("curr koords: %d, %d\n", currkords[0], currkords[1]);

        // finn mulig koords for forflytting og sjekk validitet
        if (index == 0){
            possibleKoord2D(this->koords, currkords, updatedkords, dirDir, dir, index); 
        }
        else{
            possibleKoord2D(this->koords, currkords, updatedkords, dir, dirDir, index);
        }   

        if (updatedkords[0] == currkords[0] && updatedkords[1] == currkords[1]){
            int failedIndex = 0;
            return failedIndex;
        }
        else{
            // update koords and dirArray
            
            __UpdatedEndsMonomerKoordsNeighbour__(updatedkords, index, dirDir); 
            copmuteEnergy(EnergyInteraction);

            int newEnergy = this->energy;
            double u = randomUniform(r);
            double a = std::min(1.0,exp((double)(privEnergy-newEnergy)/T));
            
            if (a > u){
                // accept the move
                return -1;
            }
            else{// revert the changes
                //printf("reverting changes at ends\n");
                __UpdatedEndsMonomerKoordsNeighbour__(currkords, index, dir);
                this->energy = privEnergy;
                return -1;          // since this move is based on probability we need to reset the failedIndex as the move has a chance of being accepted next time 
            }
        }
    }
}




