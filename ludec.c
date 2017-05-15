/*
 LU decomposition is mostly used algorithm to find the determinant of a square matrix. It divides the
 matrix into two matrices i.e. lower triangular matrix and upper triangular matrix.
 The lower triangular matrix Aij, has zeros in all indices where i<j and 1 at i=j.
 The Upper triangular matrix Aij, has zeros in all indices where i>j
 The determinant of the total matrix A is equal to the determinant of both upper and lower triangular
 matrices. And the determinant of the upper triangular matrix is the product of all the diagonal elements
 of the matrix and the determinant of the lower triangular matrix is one. So the overall determinant of
 the final matrix A, is the product of diagonal elements of the upper triangular matrix.
 So we find the lower and upper triangular matrices, and product the diagonal elements to get the
 determinant.
 The algorithm is based on 1-D partitioning where the rows are divided among the process to get the
 final output. The computed L and U matrices are stored in the matrix A itself, thereby saving space.
 The high level steps involved in the algorithm are,
    1. Process 0 generates a random nxn matrix.
    2. Then Process 0 takes its share of rows to compute and sends the remaining rows to the other
    processors.
    3. Other processors receive the rows to compute from process 0.
    4. After each processor receives their work, they perform LU decomposition.
    5. In LU decomposition, the steps are as follows:
        In Computation phase:
        A[k,j]= A[k,j]/A[k,k]
        A[i , j]=A[i , j] – A[i , k] x A[k , j] for k < i < n and k < j < n
    In communication phase:
        One to all broadcast of rows in A
    6. After the LU decomposition, then the computed values are send to process 0
    7. At the root process, the determinant of the matrix is calculated.
 */

#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
int main(int argc, char **argv) {
  int rank, world_size, l, a, b, v, w,i,j,k,matrix_size=4,size_per_proc;
  double total_time,total_send_time,total_recv_time,total_comp_time;
  float pivot,*recv_comp_Matrix,det=1;
  total_time=total_send_time=total_recv_time=total_comp_time=0;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;
  total_time -= MPI_Wtime();
  float matrix[matrix_size][matrix_size];
//generate a random matrix at process 0 
 if (rank == 0) {
    printf("The generated random matrix is=\n");
    for (i = 0; i < matrix_size; i++) { 
	for (j = 0; j < matrix_size; j++) {
        matrix[i][j] = rand()%10;
	printf("%f ",matrix[i][j]);
      }
	printf("\n");
    }
  }
  size_per_proc = matrix_size / world_size;
  int extra_size_per_proc= matrix_size % world_size;
  float Matrix_per_proc[size_per_proc][matrix_size];

//Process 0 divides the matrix and sends it to the mesh of processors for computation based on size_per_proc(a)  and tag(b)
printf("");
  if(rank == 0) {
      for (i = 0; i < size_per_proc + extra_size_per_proc; i++) {
        for (j = 0; j < matrix_size; j++) {
          Matrix_per_proc[i][j] = matrix[i][j];
        }
      }
      if (matrix_size % world_size == 0) { //if matrix size in multiples of number of processors
        for (i = 0; i < matrix_size; i++) {
	  a=i/size_per_proc;
          if (a != 0) {
            b = i % size_per_proc ;
	    total_send_time -= MPI_Wtime();
            MPI_Send(&matrix[i][0], matrix_size, MPI_FLOAT, a, b, MPI_COMM_WORLD);
	    total_send_time += MPI_Wtime();
          }
        }
      }
      else
	{
	for(i=0;i<matrix_size;i++)
	{
	a= i / (size_per_proc + extra_size_per_proc);
          if (a != 0) {
            b = i % size_per_proc ;
	    total_send_time -= MPI_Wtime();
            MPI_Send(&matrix[i][0], matrix_size, MPI_FLOAT, a, b, MPI_COMM_WORLD);
	    total_send_time += MPI_Wtime();
	}	
	}	 
	}      
     }

//All other processes  receieving data from process 0 based on tags	
    else { 
      if ( matrix_size % world_size == 0) {
        for (i = 0; i < size_per_proc; i++) {
	  total_recv_time -= MPI_Wtime();
          MPI_Recv(&Matrix_per_proc[i][0], matrix_size, MPI_FLOAT, 0, i , MPI_COMM_WORLD,&status);
	  total_recv_time += MPI_Wtime();
        }
      }
    }

//LU decomposition
  recv_comp_Matrix = (float * ) malloc(sizeof(float) * matrix_size);
   for (i = 0; i < world_size; i++) {
    for (j = 0; j < size_per_proc; j++) {
  
      if (rank == i && rank>=0) {
        for (k = 0; k < matrix_size; k++) {
          recv_comp_Matrix[k] = Matrix_per_proc[j][k];
        }
	total_send_time -= MPI_Wtime();
        MPI_Bcast(recv_comp_Matrix, matrix_size, MPI_FLOAT, rank, MPI_COMM_WORLD);
	total_send_time += MPI_Wtime();
        for (k = j; k < size_per_proc - 1; k++) {
           pivot = Matrix_per_proc[k + 1][i * size_per_proc + j] / recv_comp_Matrix[i * size_per_proc + j];
          for (l = i * size_per_proc + j; l < matrix_size; l++) {
            Matrix_per_proc[k + 1][l] -= recv_comp_Matrix[l] * pivot;
          }
          Matrix_per_proc[k + 1][i * size_per_proc + j] = pivot;
        }
      }

      else {
	total_send_time -= MPI_Wtime();
        MPI_Bcast(recv_comp_Matrix, matrix_size, MPI_FLOAT, i, MPI_COMM_WORLD);
	total_send_time += MPI_Wtime();
      }
  
      if (rank > i) {
        for (k = 0; k < size_per_proc; k++) {
          pivot = Matrix_per_proc[k][i * size_per_proc + j] / recv_comp_Matrix[i * size_per_proc + j];
          for (l = i * size_per_proc + j; l < matrix_size; l++) {
            Matrix_per_proc[k][l] -= recv_comp_Matrix[l] * pivot;
          }
          Matrix_per_proc[k][i * size_per_proc + j] = pivot;
        }
      }
    }
  }

//If not process 0, send the computed result to process 0
  if (rank != 0) {
    for (i = 0; i < size_per_proc; i++) {
      total_send_time -= MPI_Wtime();	
      MPI_Send(&Matrix_per_proc[i][0], matrix_size, MPI_FLOAT, 0, i, MPI_COMM_WORLD);
      total_send_time += MPI_Wtime();	
    }
  }

//If process 0, assign the computed work of process 0 to the final matrix
  if (rank == 0) {
    for (i = 0; i < size_per_proc; i++) {
      for (j = 0; j < matrix_size; j++) {
        matrix[i][j] = Matrix_per_proc[i][j];
      }
    }
//Receive all the computed work from other processes and assign it to final matrix
    for (i = 1; i < world_size; i++) {
      for (j = 0; j < size_per_proc; j++) {
       total_recv_time -= MPI_Wtime(); 
        MPI_Recv(recv_comp_Matrix, matrix_size, MPI_FLOAT, i, j, MPI_COMM_WORLD, & status);
	total_recv_time += MPI_Wtime();
        l = j +( i * size_per_proc);
        for (k = 0; k < matrix_size; k++) {
          matrix[l][k] = recv_comp_Matrix[k];
        }
      }
    }
   }

  if(rank==0){
//Calculation of determinant from the diagonal elements of the upper  triangular matrix
    for (i = 0; i < matrix_size; i++) {
      det = det * matrix[i][i];
    }
	printf("\n");
 
    printf("\nThe Lower and Upper matrix embbedded in Original matrix is=\n");
    for (i = 0; i < matrix_size; i++){
      for (j = 0; j < matrix_size; j++){
	printf("%f ",matrix[i][j]);
	}
      printf("\n");
	}
    printf("\n");

    printf("Determinant of the matrix is = %f\n", det);
    total_time += MPI_Wtime();
    printf("\nThe total running time = %f seconds\n",total_time);
    printf("\nThe total send time = %f seconds\n",total_send_time);
    printf("\nThe total receive time = %f seconds\n",total_recv_time);
    printf("\nThe total computation time = %f seconds\n",total_time-total_send_time-total_recv_time);
   
  }
  MPI_Finalize();
}
