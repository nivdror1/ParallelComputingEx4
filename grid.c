#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>
#include<math.h>

#define VS 5
#define NDIM 2

//mpicc -lm grid.c -o grid.exe
//mpirun -n 4 ./grid.exe

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


/**
 * save the vector to a file
 * @param vector array of size 5
 * @param rank The process rank
 */
void saveInput(int* vector, int rank_x, int rank_y, int id){
    FILE * file;
    //Create the name of the file
    char rank_x_str[3],rank_y_str[3];
    sprintf(rank_x_str, "%d", rank_x);
    sprintf(rank_y_str, "%d", rank_y);
    char name[32]= "kernel/";
    if (id==0){
        strcat(name,"input");
    }else{
        strcat(name,"output");
    }
    strcat(name,"/vector_");
    strcat(name,rank_x_str);
    strcat(name,"_");
    strcat(name,rank_y_str);

    //Open the file
    file= fopen(name,"w+");
    if (file == NULL){
        printf("error opening file");
        finalize();
        exit(1);
    }
    //Write the matrix to the file and close it
    size_t errW = 0;

    for(int i=0; i<VS; i++){
        errW = fwrite(&vector[i], sizeof(int),1,file);
    }

    int errC = fclose(file);

    if (errW != 1 || errC != 0){
        finalize();
        exit(1);
    }
}


/**
 * read the vector to a file
 * @param vector array of size 5
 * @param rank The process rank
 */
void readInput(int*** mat, int p, int id){
    FILE * file;
    int sqrt_p = (int)sqrt(p);
    //Create the name of the file

    for(int rank_x=0;rank_x<sqrt_p;rank_x++){
        for(int rank_y=0;rank_y<sqrt_p;rank_y++){
            char rank_x_str[3],rank_y_str[3];
            sprintf(rank_x_str, "%d", rank_x);
            sprintf(rank_y_str, "%d", rank_y);
            char name[32]= "kernel/";
            if (id==0){
                strcat(name,"input");
            }else{
                strcat(name,"output");
            }
            strcat(name,"/vector_");
            strcat(name,rank_x_str);
            strcat(name,"_");
            strcat(name,rank_y_str);
            //Open the file
            file= fopen(name,"r+");
            if (file == NULL){
                printf("error opening file");
                finalize();
                exit(1);
            }
            //read the vector from the file and close it
	        size_t errW = 0;
            for(int i=0; i<VS; i++){
                errW = fread(&mat[rank_x][rank_y][i], sizeof(int),1,file);
            }
            int errC = fclose(file);

            if (errW != 1 || errC != 0){
                finalize();
                exit(1);
            }
        }
    }

}

// multiply the vector v with scalar a
void scalarMul(int* v, int a){
    for (int i = 0; i < VS; ++i){
        v[i] = a * v[i];
    }
}

// sums 5 given vectors into sumVec
void vecSum(int* v1, int* v2, int* v3, int* v4, int* v5, int* sumVec){
    for (int i = 0; i < VS; ++i){
        sumVec[i] = v1[i] + v2[i] + v3[i] + v4[i] + v5[i];
    }
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
    for (int i = 0; i < VS; ++i){
        v[i] = rand()%5;
    }
}

/**
 * Create a periodic grid topology
 */
void createGridTopology(MPI_Comm *new_comm, int p){
    int ndims, reorder, periods[2], dim_size[2];
    ndims = 2; /* 2D matrix/grid */
    dim_size[0] = (int)(sqrt(p)); /* rows */
    dim_size[1] = (int)(sqrt(p)); /* columns */
    periods[0] = 1; /* row periodic (each column forms a ring) */
    periods[1] = 1; /* column periodic (each row forms a ring) */
    reorder = 1; /* allows processes reordered for efficiency */
    if(MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, new_comm)!= MPI_SUCCESS){
        finalize();
    }
}

void computeConv(int ***inM, int ***outM, int p){
    int a = (int)(sqrt(p));
    for (int i=0 ; i<a ; i++){
        for (int j=0 ; j<a ; j++){
            for (int k=0 ; k<VS ; k++){
                outM[i][j][k] = inM[i][j][k] + 1*inM[(i+a-1)%a][j][k] + 2*inM[(i+1)%a][j][k]+ 3*inM[i][(j+a-1)%a][k] + 4*inM[i][(j+1)%a][k];
            }
        }
    }
}

int compMat(int ***A, int ***B, int p){
    int a = (int)(sqrt(p));
    for (int i=0 ; i<a ; i++){
        for (int j=0 ; j<a ; j++){
            for (int k=0 ; k<VS ; k++){
                if (A[i][j][k] != B[i][j][k]){
                    printf("error! \n");
                    return -1;
                }
            }
        }
    }
    printf("well done :) \n");
    return 0;
}

void malloc3dint(int ****array, int k) {
    *array = (int ***)malloc(k*sizeof(int**));
    for (int i=0 ; i<k ; i++){
        (*array)[i] = (int **)malloc(k*sizeof(int*));
        for (int j=0 ; j<k ; j++){
            (*array)[i][j] = (int *)malloc(VS*sizeof(int));
        }
    }
}

void free3dint(int ****array, int k) {
    for (int i=0 ; i<k ; i++){
        for (int j=0 ; j<k ; j++){
            free((*array)[i][j]);
        }
        free((*array)[i]);
    }
    free((*array));
}

void print3DMatrix(int*** m, int size){
    for (int i = 0; i < size; ++i){
        for (int j = 0; j < size; ++j) {
            printf("%d\t",m[i][j][0]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    int rank, p;
    int* v;
    if(MPI_Init(&argc, &argv ) != MPI_SUCCESS){
        printf("Error in init function\n");
        return -1;
    }
    //Get the total size
    if(MPI_Comm_size(MPI_COMM_WORLD,&p) != MPI_SUCCESS){
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
        if (sqrt(p) - (int) sqrt(p) > 0.000000001) {
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
    v = (int *)malloc(VS *sizeof(int));
    //initiate the vector for every process
    initVector(v, rank);



    //Create a grid
    MPI_Comm grid_comm;
    createGridTopology(&grid_comm, p);

    //get the new rank
    if(MPI_Comm_rank(grid_comm,&rank) != MPI_SUCCESS){
        printf("Error in the rank function\n");
        finalize();
    }
    int coords[2];
    // get cartesian coordinates
    MPI_Cart_coords(grid_comm, rank, NDIM, coords);
    //save vector to disk
    saveInput(v, coords[0],coords[1],0);
//    printf("coords: %d %d  in: %d \n",coords[0], coords[1], v[0]);

    int members[5];

    MPI_Cart_shift(grid_comm, 0, 1, &members[0], &members[1]);
    MPI_Cart_shift(grid_comm, 1, 1, &members[2], &members[3]);
    members[4] = rank;
//    printf("my rank %d my neighbors: %d %d %d %d\n",members[4], members[0],members[1],members[2],members[3]);

    MPI_Request r1, r2, r3, r4, s1, s2, s3, s4;
    MPI_Isend(v, VS,MPI_INT, members[0],0 , grid_comm, &r1);
    MPI_Isend(v, VS,MPI_INT, members[1],0 , grid_comm, &r2);
    MPI_Isend(v, VS,MPI_INT, members[2],0 , grid_comm, &r3);
    MPI_Isend(v, VS,MPI_INT, members[3],0 , grid_comm, &r4);

    int rv1[VS], rv2[VS], rv3[VS], rv4[VS];
    MPI_Irecv(rv1, VS,MPI_INT, members[0],0 , grid_comm, &s1);
    MPI_Irecv(rv2, VS,MPI_INT, members[1],0 , grid_comm, &s2);
    MPI_Irecv(rv3, VS,MPI_INT, members[2],0 , grid_comm, &s3);
    MPI_Irecv(rv4, VS,MPI_INT, members[3],0 , grid_comm, &s4);

    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    MPI_Wait(&r2, MPI_STATUS_IGNORE);
    MPI_Wait(&r3, MPI_STATUS_IGNORE);
    MPI_Wait(&r4, MPI_STATUS_IGNORE);
    MPI_Wait(&s1, MPI_STATUS_IGNORE);
    MPI_Wait(&s2, MPI_STATUS_IGNORE);
    MPI_Wait(&s3, MPI_STATUS_IGNORE);
    MPI_Wait(&s4, MPI_STATUS_IGNORE);

    int res_vec[VS];
    scalarMul(rv2, 2);
    scalarMul(rv3, 3);
    scalarMul(rv4, 4);
    vecSum(v, rv1, rv2, rv3, rv4, res_vec);

    saveInput(res_vec, coords[0],coords[1],1);
//    printf("coords: %d %d  out: %d \n",coords[0], coords[1], res_vec[0]);

    MPI_Barrier(grid_comm);
    if(rank==0){

        printf("validating... \n");
        int ***inM, ***outM,***valM;
        int a = (int)(sqrt(p));
        malloc3dint(&inM, a);
        malloc3dint(&outM, a);
        malloc3dint(&valM, a);
        readInput(inM,p,0);
        readInput(outM,p,1);
        computeConv(inM, valM, p);
        compMat(outM, valM, p);
        free3dint(&inM, a);
        free3dint(&outM, a);
        free3dint(&valM, a);
    }
//    printf("rank %d my vector (%d,%d), result (%d,%d)\n",rank,v[0],v[1],res_vec[0], res_vec[1]);

    free(v);
    finalize();
}

