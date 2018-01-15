#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>
#include<math.h>


#define VECTOR_SIZE 5
/**
 * Perform the finalize function
 * @return 0 if done correctly else return -1
 */
int finalize(){
	if(MPI_Finalize() != MPI_SUCCESS){
		printf("Error in the finalize function\n");
		MPI_Abort(MPI_COMM_WORLD,1);
		return -1;
	}
	return 0;
}

/**
 * generates a random vector of 5 integers
 * @param v The vector
 * @param rank The rank of the processor
 */
void initVector(int* v, int rank)
{
	// initiate a generator for rand function
	srand((unsigned int)(time(NULL) ^ rank));
	for (int i = 0; i < VECTOR_SIZE; ++i){
		v[i] = rand()%100;
	}
}

/**
 * Create a periodic grid topology
 */
void createGridTopology(){
	MPI_Comm old_comm, new_comm;
	int ndims, reorder, periods[2], dim_size[2];
	old_comm = MPI_COMM_WORLD;
	ndims = 2; /* 2D matrix/grid */
	dim_size[0] = 3; /* rows */
	dim_size[1] = 2; /* columns */
	periods[0] = 1; /* row periodic (each column forms a ring) */
	periods[1] = 1; /* column periodic (each row forms a ring) */
	reorder = 1; /* allows processes reordered for efficiency */
	if(MPI_Cart_create(old_comm, ndims, dim_size, periods, reorder, &new_comm)!= MPI_SUCCESS){
		finalize();
	}
}

int main(int argc, char **argv) {

	int rank, size, n;
	int* v;
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
	// Check the if the size of the processors is proper
	if (rank == 0) {
		if (sqrt(size) - (int) sqrt(size) > 0.000000001) {
			printf("\n\n\nunvalid size\n\n\n\n");
			exit(1);
		}
	}
	// Synchronize all of the processors so they would'nt initiate the function if the size isn't proper
	if(MPI_Barrier(MPI_COMM_WORLD) !=  MPI_SUCCESS){
		printf("Error in the barrier function\n");
		finalize();
	}

	//allocate the vector for each process
	v = (int *)malloc(VECTOR_SIZE *sizeof(int));
	//initiate the vector for every process
	initVector(v, rank);

	//Create a grid
	createGridTopology();




}
