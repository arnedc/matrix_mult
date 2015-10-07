#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <time.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include "matrix.h"
#include <stddef.h>

int main(int argc, char **argv) {
  int n;
  clock_t c0, c1, c2,c3;
  seed=(unsigned) time(NULL);
  double max;
  int blocksize, new_size;
  dense_matrix matA, matB, matprod, matprod_blas, matprod_blocked;
  if (argc < 4){
    printf("Invalid number of input parameters. \n");
    printf("Correct use: ./matrix_mult <matrix size> <block size> <max value> \n");
    return;
  }
  n=atoi(*++argv);
  printf("Matrix size: %d \n ",n);
  blocksize=atoi(*++argv);
  printf("Block size: %d \n",blocksize);
  max=atof(*++argv);
  printf("Maximum value: %g \n ",max);
  if (n%blocksize==0)
    new_size=n;
  else
    new_size= (n/blocksize +1) * blocksize;
  printf("new size of matrices: %d \n", new_size);
  matA=(double *) calloc(new_size*new_size,sizeof(double));
  create_matrix(matA,seed,n,max,new_size);
  printdense(new_size,matA,"matrixA");
  seed++;
  matB=(double *) calloc(new_size*new_size,sizeof(double));
  create_matrix(matB,seed,n,max,new_size);
  printdense(new_size,matB,"matrixB");
  printf("matrices created \n");
  matprod=(double *) calloc(new_size*new_size,sizeof(double));
  matprod_blocked=(double *) calloc(new_size*new_size,sizeof(double));
  matprod_blas=(double *) calloc(new_size*new_size,sizeof(double));
  c0=clock();
  multiply_matrix(n,matA,matB,new_size,matprod);
  c1=clock();
  multiply_blas_matrix(new_size,matA,matB,matprod_blas);
  c2=clock();
  multiply_matrix_blocked(new_size,matA,matB,blocksize,new_size,matprod_blocked);
  c3=clock();
  printf ("\t elapsed CPU time normal:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  printf ("\t elapsed CPU time blocked:       %f\n", (float) (c3 - c2)/CLOCKS_PER_SEC);
  printf ("\t elapsed CPU time blas:          %f\n", (float) (c2 - c1)/CLOCKS_PER_SEC);
  
  //  printdense(n,matprod_blocked,"matprod_blocked");
  //printdense(n,matprod,"matprod");
  //printdense(n,matprod_blas,"matprod_blas");
  return 0;
}

/*int main(int argc, char **argv) {
  long n=10;
  clock_t c0, c1, c2;
  seed=(unsigned) time(NULL);
  double max=1;
/*  printf("Geef grootte matrix: ");
  scanf("%d",&n);
  sparse2_matrix sparse2A, sparse2B, sparse2prod,sparse2Bprod;
  struct sparse_CRS *matA, *matB, *matp;
/*  printf("Geef maximum getal in matrix: ");
  scanf("%lf",&max);
  printf("gelezen: %g \n",max);
  sparseA=sparse2sparse(n,sparse2A);
  sparse2A=create_sparse2_matrix(seed,n,max);
  ++seed;
  sparse2B=create_sparse2_matrix(seed,n,max);
  matA=create_sparse_CRS(seed,n,max);
  ++seed;
  matB=create_sparse_CRS(seed,n,max);
/*  sparseB=sparse2sparse(n,sparse2B);
/*  printf("matrices created \n");
  c0=clock();
  matp=multiply_CRS_dense(n, sparse2A, sparse2B);
  c1=clock();
  sparse2Bprod=multiply_sparse2B_matrix(n, sparse2A, sparse2B);
  c2=clock();
  printf ("\t elapsed CPU time sparse2:       %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  printf ("\t elapsed CPU time sparse:        %f\n", (float) (c2 - c1)/CLOCKS_PER_SEC);
  printsparse_CRS(n,matA,"matrix");
  /*mem1=n * (sizeof (struct cellrow)+sizeof(struct cell*)+sizeof(struct cell));
  mem2=sizeof(struct cellrow *);
  printf ("\t used memory sparse 2: %d \n", mem1);
  printf ("\t used memory sparse:   %d \n", mem2);

  printsparse2(n,sparse2A,"sparse2A");
  printsparse2(n,sparse2B,"sparse2B");
  printsparse2(n,sparse2prod,"sparse2prod");
  printsparse(n,sparseA,"sparseA");
  printsparse(n,sparseB,"sparseB");
  printsparse(n,sparseprod,"sparseprod");

  destroy_sparse2(n,sparse2A);
  destroy_sparse2(n,sparse2B);
  destroy_sparse2(n,sparse2Bprod);
  destroy_sparse2(n,sparse2prod);
  destroy_sparse_CRS(matA);
  return 0;
}*/
