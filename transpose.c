#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>

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

void initMatrix(int** m, int size, int rank)
{
	// initiate a generator for rand function
	srand((unsigned int)(time(NULL))^(unsigned int)rank);
	for (int i = 0; i < size; ++i){
		for (int j = 0; j < size; ++j) {
			m[i][j] = rand()%100;
		}
	}char str[10];

	sprintf(str, "%d", 12);
}

void saveToAFile(int** matrix,int rank, int size){
	FILE * file;
	char str[3];
	sprintf(str, "%d", rank);
	char name[9]= "matrix_";
	strcat(name,str);

	file= fopen("name","w+");
	if (file == NULL){
		printf("error opening file");//todo deal with it
	}
	int matrix_size = size*size;
	if(fwrite(matrix, sizeof(int),(size_t)size*size,file)!= matrix_size){
		printf("error writing to the file");//todo deal with it
	}
	fclose(file); //todo incase it doesn't work
}

int main(int argc, char **argv){

	int rank,size;

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
	//Allocating the matrix
	int **matrix = (int **)malloc(size * sizeof(int *));
	for (int i=0; i<size; i++)
		matrix[i] = (int *)malloc(size * sizeof(int));

	//Initiate the matrix
	initMatrix(matrix, size, rank);
	//save the matrix to a file
	saveToAFile(matrix,rank, size);
	return 0;
}

