#include <stdio.h>
#include "visual.h"


void printVec(gsl_vector *v){
    // print a vector
    for (unsigned int i = 0; i < v->size; i++)
    {
        printf("[%d], %f \n",i,gsl_vector_get(v,i));
    }
    
}

void printMat(gsl_matrix *A){
    // print a matrix
    for (unsigned int i = 0; i < A->size1; i++)
    {   
        printf("[%d], ",i);
        for (unsigned int j = 0; j < A->size2; j++)
        {
            printf("%f ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
}


void saveMat(gsl_matrix *A, char *filename){

    FILE *fptr;

    // Open file
    fptr = fopen(filename, "w");
    checkFilePointer(fptr);

    // Write 

    fprintf(fptr,"%ld\n", A->size1);
    fprintf(fptr,"%ld\n", A->size2);
    for (unsigned int i = 0; i < A->size1; i++)
    {   
        for (unsigned int j = 0; j < A->size2; j++)
        {
            fprintf(fptr,"%f,", gsl_matrix_get(A,i,j));
        }
        fprintf(fptr, "\n");

    }
    // close file
    fclose(fptr);
}


char* fileloc(char *filename, char *condition){
    // input: filename, folder
    // output: file location REMEMBERED TO FREE MEMORY OF RESULT

    char *result = static_cast<char*>(malloc(200));
    strcpy(result, filename);  // Copy str1 into result
    strcat(result, "/");   // Add a space to result
    strcat(result, condition);  // Concatenate str2 to result
    strcat(result, ".txt");  // Concatenate str2 to result

    return result;
}

