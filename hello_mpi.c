
#include <stdio.h>
#include "mpi.h"

// mpicc hello_mpi.c –o hello.exe
// mpirun -n 5 ./hello.exe


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

int main(int argc, char **argv){
	int rank, size;
	if(MPI_Init(&argc, &argv ) != MPI_SUCCESS){
		printf("Error in init function\n");
		return -1;
	}
	//get the rank
	if(MPI_Comm_rank(MPI_COMM_WORLD,&rank) != MPI_SUCCESS){
		printf("Error in the rank function\n");
		finalize();
	}
	//Get the total size
	if(MPI_Comm_size(MPI_COMM_WORLD,&size) != MPI_SUCCESS){
		printf("Error in the size function\n");
		finalize();
	}
	printf("Rank %d of %d total \n", rank, size);

	// Synchronize all of the processor up to this point
	if(MPI_Barrier(MPI_COMM_WORLD) !=  MPI_SUCCESS){
		printf("Error in the barrier function\n");
		finalize();
	}

	if (rank == 0) {
		printf("Hello World\n");
	}

	finalize();
	return 0;
}
