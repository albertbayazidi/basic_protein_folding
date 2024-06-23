#include "utils.h"

void checkFilePointer(FILE *fptr){
    if (fptr == NULL)
    {
        printf("Error: no such file/filepath\n ");
        exit(1);
    }
}

void computeDirection(int *updatedKoords, int *currKoords, int Dir){
    // update one of the koords based on the direction
    // UP = -2, down = 2, left = -1, right = 1
    if (Dir == -2){
        //printf("Dir = up\n");
        updatedKoords[1] = currKoords[1] + 1;
        updatedKoords[0] = currKoords[0];

    }
    else if (Dir == 2){
        //printf("Dir = down\n");
        updatedKoords[1] = currKoords[1] - 1;
        updatedKoords[0] = currKoords[0];
    }
    else if (Dir == -1){
        //printf("Dir = left\n");
        updatedKoords[0] = currKoords[0] - 1;
        updatedKoords[1] = currKoords[1];
    }
    else if (Dir == 1){
        //printf("Dir = right\n");
        updatedKoords[0] = currKoords[0] + 1;
        updatedKoords[1] = currKoords[1];
    }
    else{
        printf("Error: computeDirection \n");
    }
}

int findDirBasedOnKoords(int xSum,int ySum){
    // UP = -2, down = 2, left = -1, right = 1
    if (xSum == 0 && ySum == 1){
        return -2;
    }
    else if (xSum == 0 && ySum == -1){
        return 2;
    }
    else if (xSum == -1 && ySum == 0){
        return -1;
    }
    else if (xSum == 1 && ySum == 0){
        return 1;
    }
    else{
        printf("Error: findDirBasedOnKoords\n");
        return 0;
    }
}

int interpRandomDir(int dir, int randomDir){
    int newDir = 0;
    if (randomDir == 0 && abs(dir) == 2){
        newDir = 1;
    }
    else if (randomDir == 1 && abs(dir) == 2){
        newDir = -1;
    }
    else if (randomDir == 0 && abs(dir) == 1){
        newDir = 2;
    }
    else if (randomDir == 1 && abs(dir) == 1){
        newDir = -2;
    }
    else{
        printf("error in randomDir\n");
    } 
    return newDir;
}

bool checkOppositeDir(int privDir, int nextDir){
    if (-nextDir == privDir){
        return true;
    }
    else{
        return false;
    }
}

bool ileagalBoundaryPlacment(int *domainSize, int *nextKoords){
    // nextKoords[0] = x
    // nextKoords[1] = y
    if (nextKoords[0] == domainSize[0]){        //checking above boundry         // KAN HENDE DENNE MÅ VÆRE +1
        return true;
    }
    else if(nextKoords[0] == -1){               //checking bellow boundry
        return true;
    }
    else if(nextKoords[1] == domainSize[1]){    //checking right boundry        // KAN HENDE DENNE MÅ VÆRE +1
        return true;
    }
    else if(nextKoords[0] == -1){               //checking left boundry
        return true;
    }
    else{
        return false;
    }
}

void accecptUpdate(gsl_vector *koords, int *updatedKoords, int currMonomerIndex){
    gsl_vector_set(koords,currMonomerIndex,updatedKoords[0]);
    gsl_vector_set(koords,currMonomerIndex+1,updatedKoords[1]);
}

bool covalentBond(int currMonomer,int nextMonomer){
    if ((currMonomer +1 == nextMonomer) || (currMonomer -1 == nextMonomer)){
        return false;
    }
    else{
        return true;
        }
}

int occupiedBound(gsl_matrix *neighbours, int v1Monomer){
    // should check more then first slot, but should be ok in 2d case
    int bound = gsl_matrix_get(neighbours,v1Monomer,0);
    if (bound > 0){
        return 1;
    }
    else{
        return 0;
    }
}

void checkUp(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer){
    int boundPlacement = occupiedBound(neighbours, v1Monomer);
    if ((V1[0] == V2[0]) && (V1[1] == V2[1]+1)){
        gsl_matrix_set(neighbours, v1Monomer, boundPlacement, v2Monomer);
    }
}

void checkDown(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer){
    int boundPlacement = occupiedBound(neighbours, v1Monomer);
    if ((V1[0] == V2[0]) && (V1[1] == V2[1]-1)){
        gsl_matrix_set(neighbours, v1Monomer, boundPlacement, v2Monomer);
    }
}

void checkLeft(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer){
    int boundPlacement = occupiedBound(neighbours, v1Monomer);
    if ((V1[0] == V2[0]-1) && (V1[1] == V2[1])){
        gsl_matrix_set(neighbours, v1Monomer, boundPlacement, v2Monomer);
    }
}

void checkRight(gsl_matrix *neighbours, int *V1, int *V2, int v1Monomer, int v2Monomer){
    int boundPlacement = occupiedBound(neighbours, v1Monomer);
    if ((V1[0] == V2[0]+1) && (V1[1] == V2[1])){
        gsl_matrix_set(neighbours, v1Monomer, boundPlacement, v2Monomer);
    }
}

bool isPlaneConstructed2D(int privDir, int nextDir){
    if (abs(privDir) - abs(nextDir) == 0){
        return false;
    } 
    else{
        return true;
    }
}

void possibleKoord2D(gsl_vector *koords, int* currkords, int* updatedkords, int privDir, int nextDir, int index){         
    int updateDir = - privDir;
    updatedkords[0] = currkords[0];
    updatedkords[1] = currkords[1];

    computeDirection(updatedkords, currkords, updateDir);
    computeDirection(updatedkords, updatedkords, nextDir);
    checkOccupied(koords, currkords, updatedkords, index);
}

void checkOccupied(gsl_vector *koords, int* currkords, int* updatedkords, int index){
    // checks if the update koords is occupied
    // takes in the koords, the current koords, the updated koords and the index of the monomer that is being updated
    int nr_checks = -1;
    for (long unsigned int i = 0; i < koords->size; i+=2){
        if ((int)i/2 == index){
            nr_checks = nr_checks + 1;
            continue;
        }
        int ithKoord = gsl_vector_get(koords,i);
        int ip1thKoord = gsl_vector_get(koords,i+1);
 
        // increment nr_checks if the update koords is not occupied
        if (updatedkords[0] == ithKoord && updatedkords[1] == ip1thKoord){
            nr_checks = nr_checks + 0;
        }
        else{
            nr_checks = nr_checks + 1;
        }

        // checks if the update has failed, if so the update is is not valid. mabye save the monomer index so we can avoid it imidiatly
        if (nr_checks != (int)(i/2)){
            // don't update koords it is occupied, undo the update
            //printf("ilegal move \n");
            updatedkords[0] = currkords[0];
            updatedkords[1] = currkords[1]; 
            break;
        }
        else if(nr_checks == (int)(koords->size/2)-1){
            // keep update koords
            //printf("legal move \n");
            
        }
        else{
            // do nothing
        }
    }
}

bool isUnique(gsl_vector * last_dirs, int nextDir){
    for (long unsigned int i = 0; i < last_dirs->size; i++){
        if (gsl_vector_get(last_dirs,i) == nextDir){
            return false;
        }
    }
    return true;
}


int findAbsSum(gsl_vector *V){
    int last_dir_abs_sum = 0;
    for (long unsigned int i = 0; i < V->size; i++){
    last_dir_abs_sum += abs(gsl_vector_get(V,i));
    }
    return last_dir_abs_sum;
}






