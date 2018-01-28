#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>
#include<math.h>

#define N 4
#define NDIM 3

//mpicc -lm 3d-mult.c -o 3d-mult.exe
//mpirun -n 8 ./3d-mult.exe



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

void readMatrix(double** m, int id, int a, int k){
    FILE * file;
    int i,j;
	double val;
    char relPath[32] =  "matMul/ /block_";
    if (id==0){
        relPath[7] =  'A';
    }else{
        if (id==1) {
            relPath[7] = 'B';
        }else{
            relPath[7] = 'C';
        }
    }
    for(int b=0; b<a*a; b++){
        char rank_str[3];
        sprintf(rank_str, "%d", b);
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
        for (int c=0 ; c<(k*k); c++){
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


void saveMatrix(double** matrix, int a, int id, int k, int rank_x, int rank_y){
    FILE * file;
    //Create the name of the file
    char rank_str[3];
    sprintf(rank_str, "%d", a*rank_x+rank_y);
    char relPath[32] =  "matMul/ /block_";
    if (id==0){
        relPath[7] =  'A';
    }else{
        if (id==1) {
            relPath[7] = 'B';
        }else{
            relPath[7] = 'C';
        }
    }
    strcat(relPath,rank_str);

    //Open the file
    file= fopen(relPath,"w+");
    if (file == NULL){
        printf("error opening file");
        finalize();
        exit(1);
    }
    //Write the matrix to the file and close it
    int start_x = rank_x * k;
    int start_y = rank_y * k;
    for(int i=0; i<k; i++){
        int globe_i = start_y + i;
        for(int j=0; j<k;j++){
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



int malloc2Ddouble(double ***array, int k) {

    /* allocate the n*m contiguous items */
	double *p = (double *)malloc(k*k*sizeof(double));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (double **)malloc(k*sizeof(double*));
    if (!(*array)) {
        free(p);
        return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<k; i++)
        (*array)[i] = &(p[i*k]);

    return 0;
}

int free2Ddouble(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

// generates a random sizeXsize matrix of doubles
void initBlock(double** m, int size)
{
    // initiate a generator for rand function
    for (int i = 0; i < size; ++i){
        for (int j = 0; j < size; ++j) {
            m[i][j] = randfrom(0.,100.);
        }
    }
}


void matMul(double **A, double **B, double **C, int n) {
    for (int i = 0; i != n; i++) {
        for (int j = 0; j != n; j++) {
            for (int k = 0; k != n; k++) {
                C[i][j] = 0;
            }
        }
    }
    for (int i = 0; i != n; i++) {
        for (int j = 0; j != n; j++) {
            for (int k = 0; k != n; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

/**
 * Create a periodic grid topology
 */
void create3DTopology( MPI_Comm *new_comm, int a){
    int ndims, reorder, periods[3], dim_size[3];
    ndims = NDIM; /* 3D matrix/grid */
    dim_size[0] = a; /* x axis */
    dim_size[1] = a; /* y axis */
    dim_size[2] = a; /* z axis */
    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;
    reorder = 1; /* allows processes reordered for efficiency */

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, new_comm);
}

int initializeMPI(int argc, char **argv, int *p, int *a, double *t1){
    int rank;
    MPI_Init(&argc, &argv );
	*t1 =MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD,p);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // Synchronize all of the processors so they would'nt initiate the function if the size isn't proper
    *a = (int)(pow(*p, (1./NDIM)));
    if (rank == 0) {
        if (pow(*p, 1./NDIM) - *a > 0.000000001 || N/(*a) - (int)(N/(*a)) > 0.000000001) {
            printf("\n\n\ninvalid size\n\n\n\n");
            exit(1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void initGrid( MPI_Comm *grid_comm, int a, int *rank, int *coords){
    //Create a grid
    create3DTopology(grid_comm, a);
    //get the new rank
    if(MPI_Comm_rank(*grid_comm, rank) != MPI_SUCCESS){
        printf("Error in the rank function\n");
        finalize();
    }
    // get cartesian coordinates
    MPI_Cart_coords(*grid_comm, *rank, NDIM, coords);
}

void printMatrix(double** m, int size){
    for (int i = 0; i < size; ++i){
        for (int j = 0; j < size; ++j) {
            printf("%f\t",m[i][j]);
        }
        printf("\n");
    }
}


void transpose(double **A, double **TA, int n){
    for (int i=0 ; i<n ; i++){
        for (int j=0 ; j<n ; j++) {
            TA[i][j] = A[j][i];
        }
    }
}

int validate(double **A, double **B, int n){
    for (int i=0 ; i<n ; i++){
        for (int j=0 ; j<n ; j++) {
	        if ((A[i][j] - B[i][j] > 0.000000001) || (A[i][j] - B[i][j] < -0.000000001) ){
                printf("error!\n");
                return -1;
            }
        }
    }
    printf("well done :)\n");
}

int main(int argc, char **argv) {
    // rank is the process id
    // p is the total number of processes
    // a in the number of processes on each axis
    int rank, arank, brank, crank, p, a;
	double t1,t2;
    int coords[NDIM];
    MPI_Comm grid_comm;

    initializeMPI(argc, argv, &p, &a,&t1);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    initGrid( &grid_comm, a, &rank, coords);


    int k = (int) N/a;
    double **A, **B, **C, **Cout, **TB;
    malloc2Ddouble(&A, k);
    malloc2Ddouble(&TB, k);
    malloc2Ddouble(&B, k);
    malloc2Ddouble(&C, k);
    malloc2Ddouble(&Cout, k);

    srand((unsigned int)(time(NULL) ^ rank));

    if (coords[1] == 0){
        initBlock(A, k);
        saveMatrix(A, a, 0, k, coords[0], coords[2]);
//        printf("coords (%d,%d,%d) Ain %d\n", coords[0],coords[1], coords[2], A[0][0]);

    }
    if (coords[0] == 0){
        initBlock(B, k);
        saveMatrix(B, a, 1, k, coords[1], coords[2]);
//        printf("coords (%d,%d,%d) Bin %d\n", coords[0],coords[1], coords[2], B[0][0]);
    }

    MPI_Comm xComm, yComm, zComm, aComm, bComm, cComm;
    MPI_Comm_split(grid_comm, coords[0], rank, &xComm);
    MPI_Comm_split(xComm, coords[2], coords[1], &aComm);

    MPI_Comm_split(grid_comm, coords[1], rank, &yComm);
    MPI_Comm_split(yComm, coords[2], coords[0], &bComm);

    MPI_Comm_split(xComm, coords[1], coords[2], &cComm);


    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(A[0], k*k, MPI_DOUBLE, 0, aComm);
    MPI_Bcast(B[0], k*k, MPI_DOUBLE, 0, bComm);

    transpose(B, TB, k);
    matMul(TB, A, C, k);

//    printf("coords (%d,%d,%d) Aval %d\n", coords[0],coords[1], coords[2], A[0][0]);
//    printf("coords (%d,%d,%d) Bval %d\n", coords[0],coords[1], coords[2], B[0][0]);
//    printf("coords (%d,%d,%d) Cval %d\n", coords[0],coords[1], coords[2], C[0][0]);
    MPI_Barrier(MPI_COMM_WORLD);
    free2Ddouble(&TB);
    free2Ddouble(&A);
    free2Ddouble(&B);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(C[0], Cout[0], k*k, MPI_DOUBLE, MPI_SUM, 0, cComm);


    if (coords[2]==0){
        saveMatrix(Cout, a, 2, k, coords[0], coords[1]);
//        printf("coords (%d,%d,%d) Cout %d\n", coords[0],coords[1], coords[2], Cout[0][0]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free2Ddouble(&C);
    free2Ddouble(&Cout);


    if (rank==0){
	    double **TB;
        malloc2Ddouble(&A, N);
        malloc2Ddouble(&B, N);
        malloc2Ddouble(&C, N);
        malloc2Ddouble(&Cout, N);
        malloc2Ddouble(&TB, N);

        readMatrix(A, 0, a, k);
        readMatrix(B, 1, a, k);
        readMatrix(C, 2, a, k);


        transpose(B, TB, N);
        matMul(TB, A, Cout, N);


//        transpose(A, TA, N);
//        inverse(A, IA, N);
//        matMul(A, B, Cout, N);
//        matMul(B, A, Cout, N);
//        matMul(TA, TB, Cout, N);
//        matMul(TB, TA, Cout, N);
//        matMul(TA, B, Cout, N);
//        matMul(B, TA, Cout, N);
//        matMul(A, TB, Cout, N);


        validate(C, Cout, N);

//        printMatrix(A, N);
//        printf("\n");
//        printMatrix(B, N);
//        printf("\n");
//        printMatrix(C, N);
//        printf("\n");
//        printMatrix(Cout, N);

        free2Ddouble(&A);
        free2Ddouble(&B);
        free2Ddouble(&C);
        free2Ddouble(&Cout);
        free2Ddouble(&TB);
    }

    MPI_Barrier(MPI_COMM_WORLD);
	t2= MPI_Wtime();
	printf("MPI_WTIME measured is %f\n",t2-t1);
    finalize();
    return 0;


}



















