#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>
#include<math.h>





//mpicc -lm transpose.c -o transpose.exe
//mpirun -n 4 ./transpose.exe


// block size
#define BS 5
/**
 * Perform the finalize function
 * @return 0 if done correctly else return -1
 */
int finalize(){
	if(MPI_Finalize() != MPI_SUCCESS){
		printf("Error in the finalize function\n");
		MPI_Abort(MPI_COMM_WORLD,1); //todo check it up if there is need to use abort
		return -1;
	}
	return 0;
}

/* generate a random floating point number from min to max */
double randfrom(double min, double max)
{
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

// generates a random sizeXsize matrix of doubles
void initMatrix(double** m, int size, int rank)
{
	// initiate a generator for rand function
	srand((unsigned int)(time(NULL) ^ rank));
	for (int i = 0; i < size; ++i){
		for (int j = 0; j < size; ++j) {
			m[i][j] = randfrom(0.,100.);
		}
	}
}


void readMatrix(double** m, char* relPath, int size){
    FILE * file;
    int i,j;
	double val;
    char fileName[8] = "/block_";
    strcat(relPath,fileName);
    for(int p=0; p<size; p++){
        char rank_str[3];
        sprintf(rank_str, "%d", p);
        char blockPath[64] = "";
        strcat(blockPath,relPath);
        strcat(blockPath,rank_str);
        //Open the file
        file= fopen(blockPath,"r+");
        if (file == NULL){
            printf("error opening file");
            finalize();
            exit(1);
        }
        for (int b=0 ; b<(BS*BS); b++){
            fread(&i, sizeof(int),1,file);
            fread(&j, sizeof(int),1,file);
            fread(&val, sizeof(double),1,file);
            m[i][j] = val;

        }
        int errC = fclose(file);

        if (errC != 0){
            finalize();
            exit(1);
        }
    }
}





/**
 * save the matrix to a file
 * @param matrix A 2D array of size 5*5
 * @param rank The process rank
 */
void saveInput(double** matrix, int rank, int size){
	FILE * file;
	//Create the name of the file
	char rank_str[3];
	sprintf(rank_str, "%d", rank);
	char name[32]= "transpose/input/block_";
	strcat(name,rank_str);

	//Open the file
	file= fopen(name,"w+");
	if (file == NULL){
		printf("error opening file");
		finalize();
		exit(1);
	}
	//Write the matrix to the file and close it

    int p = (int)sqrt(size);
    int rank_x = rank%p;
    int rank_y = (int)floor(rank/p);
    int start_x = rank_x * BS;
    int start_y = rank_y * BS;
    for(int i=0; i<BS; i++){
        int globe_i = start_y + i;
        for(int j=0; j<BS;j++){
            int globe_j = start_x + j;
            fwrite(&globe_i, sizeof(int),1,file);
            fwrite(&globe_j, sizeof(int),1,file);
            fwrite(&matrix[i][j], sizeof(double),1,file);
        }
    }

	int errC = fclose(file);

	if (errC != 0){
		finalize();
		exit(1);
	}
}

void transposeSubMatrix(double **A, int size, int rank){
    FILE * file;
    //Create the name of the file
    char rank_str[3];
    sprintf(rank_str, "%d", rank);
    char name[32]= "transpose/output/block_";
    strcat(name,rank_str);

    //Open the file
    file= fopen(name,"w+");
    if (file == NULL){
        printf("error opening file");
        finalize();
        exit(1);
    }

	int globe_i, globe_j;
	int p = (int)sqrt(size);
	int rank_x = rank%p;
	int rank_y = (int)floor(rank/p);
	int start_x = rank_x * BS;
	int start_y = rank_y * BS;
	for(int i=0; i<BS; i++){
		for(int j=0; j<BS;j++){
			globe_i = start_y + i;
			globe_j = start_x + j;
            fwrite(&globe_j, sizeof(int),1,file);
            fwrite(&globe_i, sizeof(int),1,file);
            fwrite(&A[i][j], sizeof(double),1,file);
		}
	}
    int errC = fclose(file);

    if (errC != 0){
        finalize();
        exit(1);
    }
}

int validateTranspose(double **A, double **B, int n){
	for (int i=0; i<n ; i++){
		for (int j=0; j<n ; j++){
			if ((A[i][j] - B[j][i] > 0.000000001) || (A[i][j] - B[j][i] < -0.000000001) ){
				printf("transpose failed! :( ");
				return -1;
			}
		}
	}
	printf("transpose succeeded! :) \n");
	return 0;
}

void printMatrix(double** m, int size){
	for (int i = 0; i < size; ++i){
		for (int j = 0; j < size; ++j) {
			printf("%f\t",m[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, char **argv){

	int rank, size;
	double **M;
	if(MPI_Init(&argc, &argv ) != MPI_SUCCESS){
		printf("Error in init function\n");
		return -1;
	}
	//Get the total size
	if(MPI_Comm_size(MPI_COMM_WORLD,&size) != MPI_SUCCESS){
		printf("Error in the size function\n");
		finalize();
	}
	//get the rank
	if(MPI_Comm_rank(MPI_COMM_WORLD,&rank) != MPI_SUCCESS){
		printf("Error in the rank function\n");
		finalize();
	}

    if (rank == 0) {
        if (sqrt(size) - (int) sqrt(size) > 0.000000001) {
            printf("\n\n\nunvalid size\n\n\n\n");
            exit(1);
        }
    }
    if(MPI_Barrier(MPI_COMM_WORLD) !=  MPI_SUCCESS){
        printf("Error in the barrier function\n");
        finalize();
    }

	//Allocating the matrix
	M = (double **)malloc(BS * sizeof(double *));
	for (int i=0; i<BS; i++)
		M[i] = (double *)malloc(BS * sizeof(double));
	initMatrix(M, BS, rank);
//	printMatrix(A, n);
    saveInput(M, rank, size);

	transposeSubMatrix(M, size, rank);

    for (int i=0; i<BS; i++){
        free(M[i]);
    }
    free(M);

	if(MPI_Barrier(MPI_COMM_WORLD) !=  MPI_SUCCESS){
		printf("Error in the barrier function\n");
		finalize();
	}


	if (rank == 0){
        int n = (int)(BS*sqrt(size));
		double **A = (double **)malloc(n * sizeof(double *));
        for (int i=0; i<n; i++)
            A[i] = (double *)malloc(n * sizeof(double));
		double **B = (double **)malloc(n * sizeof(double *));
		for (int i=0; i<n; i++)
			B[i] = (double *)malloc(n * sizeof(double));

		char inputPath[50] = "transpose/input";
        readMatrix(A, inputPath, size);
		//printMatrix(A, n);
		//printf("\n");
        char outputPath[50] = "transpose/output";
        readMatrix(B, outputPath, size);
		//printMatrix(B, n);
		validateTranspose(A, B, n);

        for (int i=0; i<BS; i++) {
            free(A[i]);
            free(B[i]);
        }
        free(A);
        free(B);
	}

	finalize();

	return 0;
}



